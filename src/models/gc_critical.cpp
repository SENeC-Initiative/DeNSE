#include "gc_critical.hpp"
#include "kernel_manager.hpp"
#include "config_impl.hpp"
#include <math.h>

// Integrate the differential equation
// da1 = -k*a1*dt
namespace growth
{
GrowthCone_Critical::GrowthCone_Critical()
    : GrowthCone("competitive")
    , critical_{0, 0, 0, 0,
                CRITICAL_ELONGATION_TH,
                CRITICAL_LEAKAGE,
                CRITICAL_USE_RATIO,
                CRITICAL_CORRELATION,
                CRITICAL_VARIANCE}
    , cr_speed_{CRITICAL_ELONGATION_FACTOR,
                CRITICAL_ELONGATION_TH,
                CRITICAL_RETRACTION_FACTOR,
                CRITICAL_RETRACTION_TH}
    , demand_{CRITICAL_WEIGHT_DIAMETER,
              CRITICAL_WEIGHT_CENTRIFUGAL,0,0}

{
    observables_.push_back("resource");
}


GrowthCone_Critical::GrowthCone_Critical(const GrowthCone_Critical &copy)
    : GrowthCone(copy)
    , critical_(copy.critical_)
    , demand_(copy.demand_)
    , cr_speed_(copy.cr_speed_)


{
    neurite_dyn = biology_.own_neurite->get_branching_model();
    normal_     = std::normal_distribution<double>(0, 1);
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
    critical_.sqrt_corr = sqrt(1 - critical_.correlation * critical_.correlation);
    critical_.correlation = sqrt(critical_.correlation);
    // critical_.used = critical_.elongation_th;
    critical_.stored = cr_speed_.elongation_th/ critical_.use_ratio;
    neurite_dyn    = biology_.own_neurite->get_branching_model();
}

/**
 * @brief General purpose function run by the
 * Neurite before the split
 *
 * This function is overwritten by each mode.
 */
void GrowthCone_Critical::prepare_for_split()
{
    //critical_.stored /= 2.;
}


/**
 * @brief General purpose function run by the
 * Neurite after the split
 *
 * This function is overwritten by each model.
 */
void GrowthCone_Critical::after_split() {}

void GrowthCone_Critical::compute_CR_received()
{

    critical_.received =
        neurite_dyn->get_CR_available()*
        // local demand
        demand_.demand *
        // normalization factor for the neurite
        neurite_dyn->get_CR_quotient();
}
/**
 * @brief Compute the demand of CR at the actual step
 *
 * @param rnd_engine
 *
 * @return CR_demand
 */
void GrowthCone_Critical::compute_CR_demand(mtPtr rnd_engine)
{
    demand_.demand = critical_.stored
                   +(1+demand_.weight_centrifugal*get_centrifugal_order())
                   +(1+demand_.weight_diameter*get_diameter());
    // printf("demand is %f \n", critical_.demand);
}


/**
 * @brief Compute the speed of next step
 *
 * Translate the critical resource amount in the actual speed.
 * Verify the amount of CR in respect to the thresholds
 *
 * @param rnd_engine
 */
void GrowthCone_Critical::compute_speed(mtPtr rnd_engine, double substep)
{
    compute_CR(rnd_engine);
    move_.speed = 0;

    // if it's over the elongation threshold the neurite will outgrowth
    if (critical_.stored * critical_.use_ratio  < cr_speed_.retraction_th)
    {
        GrowthCone::move_.speed =cr_speed_.retraction_factor*
                                (critical_.stored * critical_.use_ratio - cr_speed_.retraction_th)
                                /cr_speed_.retraction_th ;
        //* (critical_.used - critical_.elongation_th);
    }
    if (critical_.stored * critical_.use_ratio > cr_speed_.elongation_th)
    {
        GrowthCone::move_.speed = cr_speed_.elongation_factor*
                                (critical_.stored * critical_.use_ratio - cr_speed_.elongation_th)/
                                (neurite_dyn->get_CR_available()*critical_.use_ratio - cr_speed_.elongation_th);
    }
    //printf(" speed computed is %f, elongation_factor is %f\n", move_.speed, cr_speed_.elongation_factor);
}


double GrowthCone_Critical::CR_differential(double a, double k,  double noise, double received, double dt)
{
    return - a*k*(dt + noise) + received;
}

void GrowthCone_Critical::compute_CR(mtPtr rnd_engine)
{

    if (not stuck_)
    {
    timestep_= kernel().simulation_manager.get_resolution();
    //compute received CR from the soma, in respect to other GC.
    compute_CR_received();
    // correlated gaussian
    critical_.stochastic_tmp = critical_.stochastic_tmp * critical_.correlation +
                           critical_.sqrt_corr * normal_(*(rnd_engine).get());
    critical_.noise = (critical_.variance * critical_.stochastic_tmp)*sqrt(timestep_);

    //consumption_rate is the stochastic correlated variable.
    //integrate by the Heun method:
    //first step: integrate with Euler
    euler_step_= critical_.stored + CR_differential(critical_.stored,
                                                    demand_.consumption_rate,
                                                    critical_.noise,
                                                    critical_.received,
                                                    timestep_);

    //second step: correct the prediction with trapezoidal rule
    critical_.stored = critical_.stored +
                    timestep_*0.5 *(
                            CR_differential(critical_.stored,
                                            demand_.consumption_rate,
                                            critical_.noise,
                                            critical_.received,
                                            timestep_)
                            +
                            //@TODO here the euler step is computed with the received
                            //      cr computed with a_i at the beginning of the time interval,
                            //      the whole process should be updated and recomputed at the end
                            //      of the interval, this means to communicate with the neurite
                            //      and recompute the total demand
                            CR_differential(euler_step_,
                                            demand_.consumption_rate,
                                            critical_.noise,
                                            critical_.received,
                                            timestep_)
                            );
    }
    //printinfo();

}

void GrowthCone_Critical::reset_CR_demand()
{
    //.demand = critical_.initial_demand;
}

double GrowthCone_Critical::get_speed() const {
return cr_speed_.elongation_factor;
}


double GrowthCone_Critical::get_CR_received() const
{
    return critical_.received;
}


double GrowthCone_Critical::get_CR_speed_factor() const
{
    return cr_speed_.elongation_factor;
}


double GrowthCone_Critical::get_CR_demand() const { return demand_.demand; }

void GrowthCone_Critical::printinfo() const
{
    printf("################ \n");
    printf("CR stored %f \n", critical_.stored);
    printf("CR demand %f \n", demand_.demand);
    printf("CR received %f \n", critical_.received);
    printf("CR elongation_th:  %f \n", cr_speed_.elongation_th);
    printf("CR elongation_th:  %f \n", cr_speed_.retraction_th);
    printf("CR elongation_factor:  %f \n", cr_speed_.elongation_factor);
    printf("CR retraction_factor:  %f \n", cr_speed_.retraction_factor);
    printf("CR use_ratio:  %f \n", critical_.use_ratio);
    printf("CR available:  %f \n", neurite_dyn->get_CR_available());

    printf("################ \n");
}


void GrowthCone_Critical::set_status(const statusMap &status)
{
    get_param(status, names::CR_elongation_factor, cr_speed_.elongation_factor);
    get_param(status, names::CR_retraction_factor, cr_speed_.retraction_factor);
    get_param(status, names::CR_elongation_th,     cr_speed_.elongation_th);
    get_param(status, names::CR_retraction_th,     cr_speed_.retraction_th);

    //use and leakage
    get_param(status, names::CR_use_ratio, critical_.use_ratio);
    get_param(status, names::CR_leakage, critical_.leakage);
    get_param(status, names::CR_correlation, critical_.correlation);
    get_param(status, names::CR_variance, critical_.variance);
    get_param(status, names::CR_weight_diameter, demand_.weight_diameter);
    get_param(status, names::CR_weight_centrifugal, demand_.weight_centrifugal);

    demand_.consumption_rate = critical_.use_ratio + 1./critical_.leakage;

    initialize_CR();

#ifndef NDEBUG
    printinfo();
#endif
}
void GrowthCone_Critical::get_status(statusMap &status) const
{
    GrowthCone::get_status(status);
    set_param(status, names::CR_elongation_factor, cr_speed_.elongation_factor);
    set_param(status, names::CR_retraction_factor, cr_speed_.retraction_factor);
    set_param(status, names::CR_elongation_th,     cr_speed_.elongation_th);
    set_param(status, names::CR_retraction_th,     cr_speed_.retraction_th);

    //use and leakage
    set_param(status, names::CR_use_ratio, critical_.use_ratio);
    set_param(status, names::CR_leakage, critical_.leakage);
    set_param(status, names::CR_correlation, critical_.correlation);
    set_param(status, names::CR_variance, critical_.variance);
    set_param(status, names::CR_weight_diameter, demand_.weight_diameter);
    set_param(status, names::CR_weight_centrifugal, demand_.weight_centrifugal);

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
    value = critical_.stored;
    ENDTRIE;

    return value;
}

} // namespace

