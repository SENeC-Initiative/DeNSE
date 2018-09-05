#include "gc_run_tumble.hpp"

// C++ includes
#include <assert.h>
#include <cmath>
#include <memory>

// kernel include
#include "kernel_manager.hpp"

// elements includes
#include "Neurite.hpp"

// libgrowth includes
#include "config_impl.hpp"
#include "growth_names.hpp"
#include "growth_time.hpp"


namespace growth
{

/**
 * This class provides a Random walk model with:
 *
 * - an history weighted direction, which is parametrized into memory_
 * - an exponential decaying correlation between gaussian, which is
 *   parametrized into corr_rw_
 * - a tunable sensing angle: the variance of the RW which determines the
 *   environment interaction.
 *
 * If any more sophisticated model for elongation isn't offered in
 * growth::set_status then a gaussian distribution for the speed is adopted.
 */
GrowthCone_RunTumble::GrowthCone_RunTumble()
    : GrowthCone("run_tumble")
    /* Model variables:: will be passed to statusMap
     defaults are initialized:
    */
    , tumbling_(false)
    , num_tumbles_(0)
{
    observables_.push_back("num_tumbles");
    initialize_RT();
}


GrowthCone_RunTumble::GrowthCone_RunTumble(const GrowthCone_RunTumble &copy)
    : GrowthCone(copy)
    , tau_(copy.tau_)
    , tumbling_(false)
    , num_tumbles_(0)
{
    exponential_rt_ = std::exponential_distribution<double>(tau_);

    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    mtPtr rng  = kernel().rng_manager.get_rng(omp_id);

    next_tumble_ = exponential_rt_(*(rng).get());
}


// Clone function

GCPtr GrowthCone_RunTumble::clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                                  double distanceToParent, std::string binaryID,
                                  const Point &position, double angle)
{
    auto newCone = std::make_shared<GrowthCone_RunTumble>(*this);
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

void GrowthCone_RunTumble::initialize_RT()
{
    // this renormalization of the "tumbling rate" is necessary to obtain
    // the correct persistence length
    tau_ =
        24. / (max_sensing_angle_ * max_sensing_angle_ * persistence_length_);

    exponential_rt_ = std::exponential_distribution<double>(tau_);

    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    mtPtr rng  = kernel().rng_manager.get_rng(omp_id);

    next_tumble_ = exponential_rt_(*(rng).get());
}


/**
 * @brief Get the current value of one of the observables
 */
double GrowthCone_RunTumble::get_state(const char *observable) const
{
    double value = 0.;

    value = GrowthCone::get_state(observable);

    TRIE(observable)
    CASE("num_tumbles")
    value = num_tumbles_;
    ENDTRIE;

    return value;
}


/**
 * @brief Intrinsic direction is not necessary for run and tumble
 *
 * Just compute the total_proba_
 */
void GrowthCone_RunTumble::compute_intrinsic_direction(
    std::vector<double> &directions_weights, double substep)
{
    total_proba_ = 0.;
    stuck_       = true;

    for (double weight : directions_weights)
    {
        if (not std::isnan(weight))
        {
            total_proba_ += weight;
            stuck_ = false;
        }
    }
}


/**
 * @brief Compute the new angle of the run and tumble
 */
Point GrowthCone_RunTumble::compute_target_position(
    const std::vector<double> &directions_weights, mtPtr rnd_engine,
    double &substep, double &new_angle)
{
    // test whether tumbling happens
    if (tumbling_)
    {
        // make a tumble (leave delta_angle_ unchanged)
        num_tumbles_++;
        // compute next tumble and reset tumbling
        next_tumble_ = exponential_rt_(*(rnd_engine).get());
        tumbling_    = false;
    }
    else
    {
        // if not tumbling check for interactions
        if (interacting_)
        {
            // take the pull into account: just scale delta_angle by the substep
            //~ delta_angle_ /=  sqrt(substep);
            //~ delta_angle_ = 0.;
        }
        else
        {
            delta_angle_ = 0.;
        }
    }

    // test whether a new tumble will happen before the end of the substep
    if (next_tumble_ - move_.module <= 0)
    {
        // we will tumble before the full module, change distance done
        double distance_done = next_tumble_;
        // compute how much time remains until tumble and update the current
        // substep accordingly.
        //~ printf("distance done: %f, module %f\n", distance_done, move_.module);
        substep = substep * (distance_done / move_.module);
        // we are still running, so set delta_angle_ to zero

        // IMPORTANT: update move_.module and set tumbling for next substep
        move_.module = distance_done;
        tumbling_    = true;
    }
    else
    {
        // no spontaneous tumble during substep, keep straight
        next_tumble_ -= move_.module;
    }

    // compute target position
    new_angle = move_.angle + delta_angle_;
    Point target_pos =
        Point(geometry_.position.at(0) + cos(new_angle) * move_.module,
              geometry_.position.at(1) + sin(new_angle) * move_.module);

    return target_pos;
}


void GrowthCone_RunTumble::prepare_for_split() {}


void GrowthCone_RunTumble::after_split() {}


//#############################################
//          Get Set status
//#############################################
void GrowthCone_RunTumble::set_status(const statusMap &status)
{
    GrowthCone::set_status(status);

    initialize_RT();
}

} // namespace growth
