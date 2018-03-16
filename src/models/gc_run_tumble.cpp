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
    , deterministic_angle_(0)
    , persistence_length_(RT_PERSISTENCE_LENGTH)
{
    initialize_RT();
}


GrowthCone_RunTumble::GrowthCone_RunTumble(const GrowthCone_RunTumble &copy)
    : GrowthCone(copy)
    , deterministic_angle_(copy.deterministic_angle_)
    , persistence_length_(copy.persistence_length_)
    , tau_(copy.tau_)
{
    initialize_RT();
}


// Clone function

GCPtr GrowthCone_RunTumble::clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                                   double distanceToParent,
                                   std::string binaryID, const Point &position,
                                   double angle)
{
#ifndef NDEBUG
    printf(" It's calling RunTumble->clone! with direction %f\n", angle);
#endif
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
    tau_ = 1./persistence_length_;
    exponential_ = std::exponential_distribution<double>(tau_);

    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    mtPtr rng = kernel().rng_manager.get_rng(omp_id);

    next_tumble_ = exponential_(*(rng).get());
}


/**
 * @brief Compute the new angle of the run and tumble
 */
Point GrowthCone_RunTumble::compute_target_position(
        const std::vector<double> &directions_weights, mtPtr rnd_engine,
        double& substep, double &new_angle)
{
    // the tumble occurs when the distance left before next_tumble_ reaches 0
    if (next_tumble_ - move_.module <= 0)
    {
        // we tumbled before the full module, change distance done
        double distance_done = next_tumble_;
        next_tumble_         = exponential_(*(rnd_engine).get());
        // tumble occured, compute how much time it took using the average speed
        // and update the current substep accordingly.
        substep = distance_done * substep / move_.module;
        //! IMPORTANT: update move_.module
        move_.module = distance_done;
    }
    else
    {
        // no spontaneous tumble occured, check for interactions
        next_tumble_ -= move_.module;
        if (interacting_)
        {
            // take the pull into account: just scale delta_angle by the substep
            //~ delta_angle_ /=  sqrt(substep);
        }
        else
        {
            delta_angle_ = 0.;
        }
    }

    // compute target position
    new_angle = move_.angle + delta_angle_;
    Point target_pos =
        Point(geometry_.position.at(0) + cos(new_angle) * move_.module,
              geometry_.position.at(1) + sin(new_angle) * move_.module);

    return target_pos;
}


void GrowthCone_RunTumble::prepare_for_split()
{
}


void GrowthCone_RunTumble::after_split()
{
}


//#############################################
//          Get Set status
//#############################################
void GrowthCone_RunTumble::set_status(const statusMap &status)
{
    GrowthCone::set_status(status);

    get_param(status, names::rt_persistence_length, persistence_length_);

    initialize_RT();
}


void GrowthCone_RunTumble::get_status(statusMap &status) const
{
    GrowthCone::get_status(status);

    set_param(status, names::rt_persistence_length, persistence_length_);
}
}
