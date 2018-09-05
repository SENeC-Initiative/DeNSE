#include "gc_self_referential_forces.hpp"

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
#include "kernel_manager.hpp"


namespace growth
{

double generateGaussianNoise(double mean, double stdDev)
{
    double u, v, s;

    do
    {
        u = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        v = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
        s = u * u + v * v;
    } while ((s >= 1.0) or (s == 0.0));

    s = sqrt(-2.0 * log(s) / s);

    return mean + stdDev * u * s;
}


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
GrowthCone_SelfReferentialForces::GrowthCone_SelfReferentialForces()
    : GrowthCone("self_referential_forces")
    /* Model variables:: will be passed to statusMap
     defaults are initialized:
    */
    , inertial_{SRF_INERTIAL_FORCE, SRF_INERTIAL_DECAY}
    , somatropic_{SRF_SOMATROPIC_FORCE, SRF_SOMATROPIC_DECAY}
    // self avoidance or self affinity can be introduced together with
    // fasciculation model
    , selfavoidance_{SRF_AVOIDANCE_FORCE, SRF_AVOIDANCE_DECAY}
{
}


GrowthCone_SelfReferentialForces::GrowthCone_SelfReferentialForces(
    const GrowthCone_SelfReferentialForces &copy)
    : GrowthCone(copy)
    , inertial_(copy.inertial_)
    , somatropic_(copy.somatropic_)
    , selfavoidance_(copy.selfavoidance_)
{
}


// Clone function

GCPtr GrowthCone_SelfReferentialForces::clone(
    BaseWeakNodePtr parent, NeuritePtr neurite, double distanceToParent,
    std::string binaryID, const Point &position, double angle)
{
#ifndef NDEBUG
    printf(" It's calling RandomWalk->clone! with direction %f\n", angle);
#endif
    auto newCone = std::make_shared<GrowthCone_SelfReferentialForces>(*this);
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


/**
 * @brief Compute the new step with self_referential_forces model
 *  This model relay on choosing next step inside a sphere (or circle),
 *  the sphere center is positioned by inertial, environmental, etc... drifts.
 *  the sphere radius is proportional to the variance of the process.
 *  The sphere size is rescaled with module length:
 *  forces and variance are defined for the unit module and then expanded.
 *
 *
 */
Point GrowthCone_SelfReferentialForces::compute_target_position(
    const std::vector<double> &directions_weights, mtPtr rnd_engine,
    double &substep, double &new_angle)
{
    // center the position of the "sphere" the growth cone actual position
    // double weight =1;
    Point new_center    = Point(0, 0);
    Point actual_center = biology_.branch->get_last_xy();
    size_t idx          = 1;

    // run over previous direction and weight over the whole path
    // while (weight > 0.01)
    //{
    // PointArray previous = biology_.branch->at(idx);
    // Point previous_point = Point(previous[0], previous[1]);

    // double distance = previous[2];
    // Point vector = actual_center - previous_point;

    // idx += 1;

    //}

    // printf("force is %f %f\n",inertial_.force, inertial_.decay );
    // get neurite information to compute self referential forces
    double distanceToSoma = biology_.branch->final_distance_to_soma();
    double soma_angle     = biology_.own_neurite->get_soma_angle();

    // apply drift from previous points
    new_center.shift(inertial_.force, move_.angle);
    // apply drift from soma
    new_center.shift(somatropic_.force *
                         pow(distanceToSoma, -somatropic_.decay),
                     soma_angle);
    // aply drift from environment

    // rescale step with step_length:
    double step         = sqrt(pow(new_center[0], 2) + pow(new_center[1], 2));
    double scale_factor = (move_.speed / step);
    new_center =
        Point(new_center[0] * scale_factor, new_center[1] * scale_factor);

    double sigma_angle = sqrt(substep) * move_.sigma_angle;

    double x = generateGaussianNoise(new_center[0], scale_factor * sigma_angle);
    double y = generateGaussianNoise(new_center[1],
                                     scale_factor * scale_factor * sigma_angle);

    // compute the angle and fallback to polar coordinate system already
    // implemented
    new_angle = atan2(y, x);

    // check that this step is allowed, otherwise move towards default_angle
    Point target_pos =
        Point(geometry_.position.at(0) + cos(new_angle) * move_.module,
              geometry_.position.at(1) + sin(new_angle) * move_.module);

    return target_pos;
}


void GrowthCone_SelfReferentialForces::prepare_for_split() {}


void GrowthCone_SelfReferentialForces::after_split() {}


//#############################################
//          Get Set status
//#############################################
void GrowthCone_SelfReferentialForces::set_status(const statusMap &status)
{
    GrowthCone::set_status(status);
    get_param(status, names::srf_avoidance_decay, selfavoidance_.decay);
    get_param(status, names::srf_avoidance_force, selfavoidance_.force);
    get_param(status, names::srf_inertial_force, inertial_.force);
    get_param(status, names::srf_inertial_decay, inertial_.decay);
    get_param(status, names::srf_somatropic_force, somatropic_.force);
    get_param(status, names::srf_somatropic_decay, somatropic_.decay);
}


void GrowthCone_SelfReferentialForces::get_status(statusMap &status) const
{

    GrowthCone::get_status(status);
    set_param(status, names::srf_avoidance_decay, selfavoidance_.decay, "");
    set_param(status, names::srf_avoidance_force, selfavoidance_.force, "");

    set_param(status, names::srf_inertial_force, inertial_.force, "");
    set_param(status, names::srf_inertial_decay, inertial_.decay, "");
    set_param(status, names::srf_somatropic_force, somatropic_.force, "");
    set_param(status, names::srf_somatropic_decay, somatropic_.decay, "");
}

} // namespace growth
