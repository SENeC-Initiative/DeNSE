#include "gc_critical.hpp"
#include "config_impl.hpp"
#include "kernel_manager.hpp"
#include <math.h>

// Integrate the differential equation
// da1 = -k*a1*dt
namespace growth
{
GrowthCone_Critical::GrowthCone_Critical()
    : GrowthCone("competitive")
    , stochastic_tmp_(0)
    , received_(0)
    , sqrt_corr_(0)
    , noise_(0)
    , stored_(CRITICAL_ELONGATION_TH)
    , leakage_(CRITICAL_LEAKAGE)
    , use_ratio_(CRITICAL_USE_RATIO)
    , correlation_(CRITICAL_CORRELATION)
    , variance_(CRITICAL_VARIANCE)
    , elongation_factor_(CRITICAL_ELONGATION_FACTOR)
    , elongation_th_(CRITICAL_ELONGATION_TH)
    , retraction_factor_(CRITICAL_RETRACTION_FACTOR)
    , retraction_th_(CRITICAL_RETRACTION_TH)
    , weight_diameter_(CRITICAL_WEIGHT_DIAMETER)
    , weight_centrifugal_(CRITICAL_WEIGHT_CENTRIFUGAL)
{
    observables_.push_back("resource");
    consumption_rate_ = use_ratio_ + 1. / leakage_;
}


GrowthCone_Critical::GrowthCone_Critical(const GrowthCone_Critical &copy)
    : GrowthCone(copy)
    , stochastic_tmp_(0)
    , received_(copy.received_)
    , sqrt_corr_(copy.sqrt_corr_)
    , noise_(copy.noise_)
    , stored_(copy.stored_)
    , leakage_(copy.leakage_)
    , use_ratio_(copy.use_ratio_)
    , consumption_rate_(copy.consumption_rate_)
    , correlation_(copy.correlation_)
    , variance_(copy.variance_)
    , elongation_factor_(copy.elongation_factor_)
    , elongation_th_(copy.elongation_th_)
    , retraction_factor_(copy.retraction_factor_)
    , retraction_th_(copy.retraction_th_)
    , weight_diameter_(copy.weight_diameter_)
    , weight_centrifugal_(copy.weight_centrifugal_)
{
    normal_ = std::normal_distribution<double>(0, 1);
    observables_.insert(observables_.end(), copy.observables_.begin(),
                        copy.observables_.end());
}


GCPtr GrowthCone_Critical::clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                                 double distanceToParent, std::string binaryID,
                                 const Point &position, double angle)
{
#ifndef NDEBUG
    printf(" It's calling Critical->clone! with direction %f\n", angle);
#endif
    auto newCone = std::make_shared<GrowthCone_Critical>(*this);
    int omp_id   = kernel().parallelism_manager.get_thread_local_id();
    newCone->update_topology(parent, neurite, distanceToParent, binaryID,
                             position, angle);

    // update containing area
    newCone->current_area_ =
        using_environment_
            ? kernel().space_manager.get_containing_area(position, omp_id)
            : "";

    if (using_environment_)
    {
        newCone->update_growth_properties(current_area_);
    }

    return newCone;
}

void GrowthCone_Critical::initialize_CR()
{
    sqrt_corr_ = sqrt(1 - correlation_ * correlation_);
    stored_    = (stored_ == 0) ? 1.5 * elongation_th_ : stored_;
}


/**
 * @brief General purpose function run by the
 * Neurite before the split
 *
 * This function is overwritten by each mode.
 */
void GrowthCone_Critical::prepare_for_split()
{
    stored_ *= 0.5;
    received_ *= 0.5;
}


/**
 * @brief General purpose function run by the
 * Neurite after the split
 *
 * This function is overwritten by each model.
 */
void GrowthCone_Critical::after_split() {}


void GrowthCone_Critical::compute_CR_received(double substep)
{
    received_ = biology_.own_neurite->get_available_cr() *
                // local demand (weighted)
                get_CR_demand() *
                // normalization factor for the neurite
                biology_.own_neurite->get_quotient_cr();
}


/**
 * @brief Compute the demand of CR at the actual step
 *
 * @param rnd_engine
 *
 * @return CR_demand
 */
double GrowthCone_Critical::get_CR_demand()
{
    double current_demand = stored_ * consumption_rate_;
    // weight by centrifugal order and diameter if required
    //~ +(1+weight_centrifugal_*get_centrifugal_order())
    //~ +(1+weight_diameter_*get_diameter());

    return current_demand;
}


/**
 * @brief Compute the speed of next step
 *
 * Translate the critical resource amount in the actual speed.
 * Verify the amount of CR in respect to the thresholds
 *
 * @param rnd_engine
 */
double GrowthCone_Critical::compute_cr_speed(mtPtr rnd_engine, double substep)
{
    double speed = 0.;

    // if it's over the elongation threshold the neurite will extend
    if (stored_ < retraction_th_)
    {
        speed =
            retraction_factor_ * (stored_ - retraction_th_) / retraction_th_;
    }
    else if (stored_ >= elongation_th_)
    {
        speed = elongation_factor_ * (stored_ - elongation_th_) /
                (stored_ + elongation_th_);
    }

    return speed;
}


/**
 * @brief Compute the speed of next step
 *
 * NOTE: for the critical model, the speed is computed by the neurite as the
 * average speed over the substep, i.e. it has be set previously using the
 * results of compute_cr_speed and does not need to be recomputed here.
 */
void GrowthCone_Critical::compute_speed(mtPtr rnd_engine, double substep)
{
    move_.speed = 0;

    // if it's over the elongation threshold the neurite will extend
    if (stored_ < retraction_th_)
    {
        move_.speed =
            retraction_factor_ * (stored_ - retraction_th_) / retraction_th_;
    }
    else if (stored_ >= elongation_th_)
    {
        move_.speed = elongation_factor_ * (stored_ - elongation_th_) /
                      (stored_ + elongation_th_);
    }
}


double GrowthCone_Critical::compute_CR(mtPtr rnd_engine, double substep)
{
    if (not stuck_)
    {
        // compute received CR from the soma, with respect to other GC.
        compute_CR_received(substep);

        // correlated gaussian (unit standard deviation)
        noise_ =
            noise_ * correlation_ + sqrt_corr_ * normal_(*(rnd_engine).get());

        // then use Euler for deterministic term and corrected Wiener term
        stored_ = stored_ +
                  substep * (received_ - stored_ * consumption_rate_) +
                  sqrt(substep) * variance_ * noise_;

        // amount of molecule cannot be negative
        if (stored_ < 0.)
        {
            stored_ = 0.;
        }
    }

    return stored_;
}

void GrowthCone_Critical::reset_CR_demand()
{
    //.demand = critical_.initial_demand;
}

double GrowthCone_Critical::get_speed() const { return elongation_factor_; }


double GrowthCone_Critical::get_CR_received() const { return received_; }


double GrowthCone_Critical::get_CR_speed_factor() const
{
    return elongation_factor_;
}


void GrowthCone_Critical::printinfo() const
{
    printf("################ \n");
    printf("CR stored %f \n", stored_);
    printf("CR received %f \n", received_);
    printf("CR elongation_th:  %f \n", elongation_th_);
    printf("CR elongation_th:  %f \n", retraction_th_);
    printf("CR elongation_factor:  %f \n", elongation_factor_);
    printf("CR retraction_factor:  %f \n", retraction_factor_);
    printf("CR use_ratio:  %f \n", use_ratio_);
    //~ printf("CR available:  %f \n",
    //biology_.own_neurite->get_available_cr(substep));

    printf("################ \n");
}


void GrowthCone_Critical::set_status(const statusMap &status)
{
    // state parameters
    get_param(status, names::resource, stored_);

    // speed-related stuff
    get_param(status, names::CR_elongation_factor, elongation_factor_);
    get_param(status, names::CR_retraction_factor, retraction_factor_);
    get_param(status, names::CR_elongation_th, elongation_th_);
    get_param(status, names::CR_retraction_th, retraction_th_);

    // use and leakage
    get_param(status, names::CR_use_ratio, use_ratio_);
    get_param(status, names::CR_leakage, leakage_);
    get_param(status, names::CR_correlation, correlation_);
    get_param(status, names::CR_variance, variance_);
    get_param(status, names::CR_weight_diameter, weight_diameter_);
    get_param(status, names::CR_weight_centrifugal, weight_centrifugal_);

    consumption_rate_ = use_ratio_ + 1. / leakage_;

    initialize_CR();

#ifndef NDEBUG
    printinfo();
#endif
}


void GrowthCone_Critical::get_status(statusMap &status) const
{
    GrowthCone::get_status(status);

    // state parameters
    set_param(status, names::resource, stored_);

    // speed-related
    set_param(status, names::CR_elongation_factor, elongation_factor_);
    set_param(status, names::CR_retraction_factor, retraction_factor_);
    set_param(status, names::CR_elongation_th, elongation_th_);
    set_param(status, names::CR_retraction_th, retraction_th_);

    // use and leakage
    set_param(status, names::CR_use_ratio, use_ratio_);
    set_param(status, names::CR_leakage, leakage_);
    set_param(status, names::CR_correlation, correlation_);
    set_param(status, names::CR_variance, variance_);
    set_param(status, names::CR_weight_diameter, weight_diameter_);
    set_param(status, names::CR_weight_centrifugal, weight_centrifugal_);
}


/**
 * @brief Get the current value of one of the observables
 */
double GrowthCone_Critical::get_state(const char *observable) const
{
    double value = 0.;

    value = GrowthCone::get_state(observable);

    TRIE(observable)
    CASE("resource")
    value = stored_;
    ENDTRIE;

    return value;
}

} // namespace growth
