#include "gc_critical.hpp"
#include "config_impl.hpp"
#include <math.h>

namespace growth
{

GrowthCone_Critical::GrowthCone_Critical()
    : GrowthCone()
    , critical_(0, CRITICAL_INITIAL_DEMAND, CRITICAL_INITIAL_DEMAND,
                CRITICAL_LEAKAGE, CRITICAL_RETRACTION_TH,
                CRITICAL_ELONGATION_TH, CRITICAL_SPEED_FACTOR,
                CRITICAL_USE_RATIO)
    , demand_(CR_DEMAND_MEAN, CR_DEMAND_CORRELATION, 0, CR_DEMAND_STDDEV,
              CR_DEMAND_MEAN)
{
    observables_.push_back("resource");
}


GrowthCone_Critical::GrowthCone_Critical(const GrowthCone_Critical &copy)
    : GrowthCone(copy)
    , critical_(copy.critical_)
    , demand_(copy.demand_)


{
    neurite_dyn = biology_.own_neurite->get_branching_model();
    normal_     = std::normal_distribution<double>(0, 1);
    observables_.insert(observables_.cend(), copy.observables_.cbegin(),
                        copy.observables_.cend());
}


GCPtr GrowthCone_Critical::clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                                 double distanceToParent, std::string binaryID,
                                 const Point &position, double angle)
{
#ifndef NDEBUG
    printf(" It's calling Critical->clone! with direction %f\n", angle);
#endif
    auto newCone = std::make_shared<GrowthCone_Critical>(*this);
    newCone->update_topology(parent, neurite, distanceToParent, binaryID,
                             position, angle);
    return newCone;
}

void GrowthCone_Critical::initialize_CR()
{
    demand_.sqrt_corr = sqrt(1 - demand_.correlation * demand_.correlation);
    // critical_.used = critical_.elongation_th;
    critical_.left = critical_.elongation_th / critical_.use_ratio;
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
    critical_.left /= 2.;
    demand_.local_demand /= 2.;
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
        // local demand
        critical_.demand *
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
double GrowthCone_Critical::compute_CR_demand(mtPtr rnd_engine)
{

    // implement a correlated Gaussian whit fixed -user defined- correlation
    // coefficient
    demand_.local_demand = demand_.local_demand * demand_.correlation +
                           demand_.sqrt_corr * normal_(*(rnd_engine).get());

    critical_.demand = (demand_.mean + demand_.std_dev * demand_.local_demand);

    // attenuation of growth cone demanding
    // strength with topological distance
    // printf("name is %s \n", get_treeID().c_str());
    // printf("local demand is %f, sqrt %f\n", demand_.local_demand,
    // demand_.sqrt_corr);
    // printf("demand is %f \n", critical_.demand);

    return critical_.demand;
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
    compute_CR();
    move_.speed = 0;

    // if it's over the elongation threshold the neurite will outgrowth
    if (critical_.used > critical_.elongation_th)
    {
        GrowthCone::move_.speed = critical_.used * critical_.speed_factor;
        //* (critical_.used - critical_.elongation_th);
    }
    if (critical_.used < critical_.retraction_th)
    {
        GrowthCone::move_.speed = -critical_.speed_factor *
                                  (critical_.retraction_th - critical_.used);
    }
    // printf(" speed computed is %f \n", move_.speed);
}

void GrowthCone_Critical::compute_CR() {}

void GrowthCone_Critical::reset_CR_demand()
{
    //.demand = critical_.initial_demand;
}


double GrowthCone_Critical::get_CR_received() const
{
    return critical_.received;
}

double GrowthCone_Critical::get_CR_left() const { return critical_.left; }

double GrowthCone_Critical::get_CR_used() const { return critical_.used; }


double GrowthCone_Critical::get_CR_speed_factor() const
{
    return critical_.speed_factor;
}


double GrowthCone_Critical::get_CR_demand() const { return critical_.demand; }

void GrowthCone_Critical::printinfo() const
{
    printf("CR demand %f \n", critical_.demand);
    printf("CR received %f \n", critical_.received);
    printf("CR elongation_th:  %f \n", critical_.elongation_th);
    printf("CR retraction_th:  %f \n", critical_.retraction_th);
    printf("CR use_ratio:  %f \n", critical_.use_ratio);
    printf("after CR used:  %f \n", critical_.used);
    printf("after CR left:  %f \n", critical_.left);
    printf("################ \n");
}

void GrowthCone_Critical::set_status(const statusMap &status)
{

    get_param(status, names::CR_speed_factor, critical_.speed_factor);
    get_param(status, names::CR_elongation_th, critical_.elongation_th);
    get_param(status, names::CR_retraction_th, critical_.retraction_th);

    get_param(status, "left", critical_.left);
    get_param(status, names::CR_use_ratio, critical_.use_ratio);
    get_param(status, names::CR_leakage, critical_.leakage);

    get_param(status, names::CR_demand_correlation, demand_.correlation);
    if (get_param(status, names::CR_initial_demand, critical_.initial_demand))
    {
        critical_.demand = critical_.initial_demand;
    }
    get_param(status, names::CR_demand_mean, demand_.mean);
    get_param(status, names::CR_demand_stddev, demand_.std_dev);

    initialize_CR();

#ifndef NDEBUG
    printf("\n"
           "CRITICAL RESOURCE MODEL \n"
           "CR_leakage %f \n"
           "CR_elongation_th %f \n"
           "CR_retraction_th %f \n"
           "CR_initial_demand %f \n"
           "CR_speed_factor  %f \n"
           "Demand model \n"
           "CR_demand:  %f +- %f with correlation %f \n ",
           critical_.leakage, critical_.elongation_th, critical_.retraction_th,
           critical_.initial_demand, critical_.speed_factor, demand_.mean,
           demand_.std_dev, demand_.correlation);
#endif
}
void GrowthCone_Critical::get_status(statusMap &status) const
{
    GrowthCone::get_status(status);
    set_param(status, names::CR_speed_factor, critical_.speed_factor);
    set_param(status, names::CR_elongation_th, critical_.elongation_th);
    set_param(status, names::CR_retraction_th, critical_.retraction_th);

    set_param(status, "left", critical_.left);
    set_param(status, names::CR_use_ratio, critical_.use_ratio);
    set_param(status, names::CR_leakage, critical_.leakage);

    set_param(status, names::CR_demand_correlation, demand_.correlation);
    set_param(status, names::CR_initial_demand, critical_.demand);

    set_param(status, names::CR_demand_mean, demand_.mean);
    set_param(status, names::CR_demand_stddev, demand_.std_dev);
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
    value = critical_.left;
    ENDTRIE;

    return value;
}


//============================
// Detailed
//============================


GrowthCone_Critical_Langevin::GrowthCone_Critical_Langevin()
    : GrowthCone_Critical()
{
}


GrowthCone_Critical_Langevin::GrowthCone_Critical_Langevin(
    const GrowthCone_Critical &copy)
    : GrowthCone_Critical(copy)
{
}


/**
 * @brief Compute the critical_resource received from the neurite
 *
 * Compare the demand of the growth cone to the demand of the other growth cones
 * in the neurite
 * The summation of the demand is done in Branching.cpp.
 * The competition behaviour is here, and  the attenuation is computed too.
 * the sum of the critical_resource_received from all the growth cones in the
 * neurite is
 * equal to the critical_resource
 * amount of the neurite.
 *
 */
void GrowthCone_Critical_Langevin::compute_CR()
{

    // compute received normalization on the neurite.
    compute_CR_received();

    critical_.left +=
        (critical_.received) - (critical_.left) * (critical_.leakage);
    critical_.used = (critical_.left) * critical_.use_ratio;
    critical_.left *= (1 - critical_.use_ratio);


    /*    printf("centrifugal_order: %f x %f x %f x %f \n = %f"*/
    //,powf(2, -centrifugalOrder_ * critical_resource_topo_coeff_ )
    //, neurite_dyn->get_critical_resource_quotient()
    //, critical_resource_demand_
    //, neurite_dyn->get_critical_resource_amount(),
    /*critical_resource_received_);*/
}

GrowthCone_Critical_Lurd::GrowthCone_Critical_Lurd()
    : GrowthCone_Critical()
{
}


GrowthCone_Critical_Lurd::GrowthCone_Critical_Lurd(
    const GrowthCone_Critical &copy)
    : GrowthCone_Critical(copy)
{
}
void GrowthCone_Critical_Lurd::compute_CR()
{

    compute_CR_received();


    critical_.used =
        (critical_.left + critical_.received) * critical_.use_ratio;
    critical_.left =
        (critical_.received + critical_.left) * (1 - critical_.use_ratio);

    // printinfo();

    /*    printf("centrifugal_order: %f x %f x %f x %f \n = %f"*/
    //,powf(2, -centrifugalOrder_ * critical_resource_topo_coeff_ )
    //, neurite_dyn->get_critical_resource_quotient()
    //, critical_resource_demand_
    //, neurite_dyn->get_critical_resource_amount(),
    /*critical_resource_received_);*/
}

GrowthCone_Critical_Gaussian::GrowthCone_Critical_Gaussian()
    : GrowthCone_Critical()
{
}

GrowthCone_Critical_Gaussian::GrowthCone_Critical_Gaussian(
    const GrowthCone_Critical &copy)
    : GrowthCone_Critical(copy)
{
}

void GrowthCone_Critical_Gaussian::compute_CR()
{


    compute_CR_received();

    critical_.left = 0;
    critical_.used = critical_.received * critical_.use_ratio;

    // printinfo();
}

} // namespace
