#include "GrowthCone.hpp"

// C++ includes
#define _USE_MATH_DEFINES
#include <assert.h>
#include <cmath>
#include <memory>
#include <set>

#include <fstream>
#include <iostream>

// kernel include
#include "kernel_manager.hpp"

// elements includes
#include "Branch.hpp"
#include "growth_names.hpp"

#define SQRT_FRAC_1_2PI 0.3989422804014327


namespace growth
{

/*
 * Test function to see whether some directions are accessible for the
 * growth cone
 */
bool allnan(const std::vector<double> &weights)
{
    bool all_nan = true;
    for (double w : weights)
    {
        all_nan *= std::isnan(w);
    }
    return all_nan;
}


GrowthCone::GrowthCone()
    : TopologicalNode()
    , observables_({"length", "speed", "angle"})
    , delta_angle_(0)
    , stuck_(false)
    , filopodia_{{},
                 {},
                 FILOPODIA_ANGULAR_RES,
                 FILOPODIA_FINGER_LENGTH,
                 FILOPODIA_WALL_AFFINITY}
    , move_()
    , speed_growth_cone_(RW_SPEED_GROWTH_CONE)
    , speed_variance_(0)
    , rw_sensing_angle_(RW_SENSING_ANGLE)
    , average_speed_(0)
{
    // initialize move variables
    move_.sigma_angle = rw_sensing_angle_;
    move_.speed       = speed_growth_cone_;
    // random distributions
    normal_  = std::normal_distribution<double>(0, 1);
    uniform_ = std::uniform_real_distribution<double>(0., 1.);
    update_kernel_variables();
}


GrowthCone::GrowthCone(const GrowthCone &copy)
    : TopologicalNode(copy)
    , observables_(copy.observables_)
    , delta_angle_(0)
    , stuck_(false)
    , filopodia_(copy.filopodia_)
    , move_(copy.move_)
    , speed_growth_cone_(copy.speed_growth_cone_)
    , speed_variance_(copy.speed_variance_)
    , rw_sensing_angle_(copy.rw_sensing_angle_)
    , average_speed_(copy.average_speed_)
{
    normal_  = std::normal_distribution<double>(0, 1);
    uniform_ = std::uniform_real_distribution<double>(0., 1.);
    update_kernel_variables();
}


GrowthCone::~GrowthCone()
{
    assert(biology_.branch.use_count() == 1);
#ifndef NDEBUG
    printf("#### cone %s deletion #####\n", topology_.binaryID.c_str());
#endif
}


/**
 * Update the geometrical and topological informations of the
 * GrowthCone.
 * This is called to update a GrowthCone after it has been cloned from
 * when branching occurs.
 */
void GrowthCone::update_topology(BaseWeakNodePtr parent, NeuritePtr own_neurite,
                                 float distance_to_parent,
                                 const std::string &binaryID, Point position,
                                 double angle)
{
    topology_ = NodeTopology(parent, parent.lock()->get_centrifugal_order() + 1,
                             false, 0, binaryID);
    geometry_ = NodeGeometry(
        position, parent.lock()->get_distance_to_soma() + distance_to_parent,
        distance_to_parent);
    biology_ = NodeBiology(
        false, std::make_shared<Branch>(position, geometry_.dis_to_soma),
        own_neurite, 1);
    move_.angle = angle;
}


// ###########################################################
//      Grow, Retraction and Pruning
// ###########################################################


/**
 * @brief let's compute next step
 * Core function for the simulation, it shouldn't be overwritten
 *
 * Explanation:
 * 1.   Compute the module of the next step
 *      It will be done from the elongation model, critical_resourcee or
 * random_walk
 *
 * 2.   If module is positive, compute the angle and elongate:
 *      To grow it has to sense the environment and then convolve with the
 *      particular intrinsic distribution (e.g. random walk model)
 *      as set in the grwth cone model.
 *      'sense_environment' and 'accesible_environment'
 *      are defined in the abstract class GrowthCone (here),
 *      but it is required they are called inside the function
 *      'compute_directions' overriden by the GrowthCone model.
 *      This requirement is due to:
 *          + the functions need the sigma (persistence_length) defined in the
 *          models
 *          + some models required to perform a rotation before sense the
 *          environment
 *
 * 3.   If module is negative retract:
 *      The retraction of a 0 length branch will imply the pruning of the
 *      GrowthCone
 *
 * @param rnd_engine
 */
void GrowthCone::grow(mtPtr rnd_engine, size_t cone_n, double substep)
{
    compute_speed(rnd_engine, substep);
    compute_module(substep);
    if (move_.module > 0)
    {
        // check environment and mechanical interactions
        std::vector<double> directions_weights;
        compute_pull_and_accessibility(
            directions_weights, rnd_engine, substep);

        if (not stuck_)
        {
            // weight the directions by the intrinsic preference of the growth
            // cone
            compute_intrinsic_direction(directions_weights);
            // set delta_angle_ as main pulling direction
            choose_pull_direction(directions_weights, rnd_engine);
            // compute new direction based on delta_angle_
            compute_new_direction(rnd_engine, substep);
            // store new position
            geometry_.position = Point(
                geometry_.position.at(0) + move_.module * cos(move_.angle),
                geometry_.position.at(1) + move_.module * sin(move_.angle));
            biology_.branch->add_point(geometry_.position, move_.module);
        }
        // reset move_.sigma_angle to its default value
        move_.sigma_angle = rw_sensing_angle_;
    }
    if (move_.module < 0)
    {
        printf("retracting for cone %lu\n", cone_n);
        retraction(move_.module);
        if (biology_.branch->size() == 0)
        {
            prune(cone_n);
        }
        geometry_.position = biology_.branch->get_last_xy();
    }
}


void GrowthCone::compute_module(double substep) {
    move_.module = move_.speed * substep;
}


void GrowthCone::retraction(double module)
{
    double subtracted = 0;

    while (module + subtracted < 0)
    {
        if (biology_.branch->size() > 0)
        {
            subtracted += biology_.branch->points[2].back();
            biology_.branch->retract();
            subtracted -= biology_.branch->points[2].back();
        }
        else
        {
            break;
        }
    }

    move_.module = -subtracted;
}


void GrowthCone::prune(size_t cone_n)
{
    biology_.own_neurite->delete_cone(cone_n);
}


// ###########################################################
//              Space management
// ###########################################################

/**
 * @brief compute likelihood of stepping in each available direction.
 *
 * Scan surrounding environment and assess whether a step in each of the
 * possible directions is possible, and, if so, how likely it is.
 */
void GrowthCone::compute_pull_and_accessibility(
    std::vector<double> &directions_weights, mtPtr rnd_engine, double substep)
{
    if (using_environment_)
    {
        // unstuck neuron
        stuck_ = false;
        // initialize direction test variables
        bool all_nan = true;
        for (int i = 0; i < filopodia_.size; i++)
        {
            directions_weights.push_back(1.);
        }
        while (all_nan and not stuck_)
        {
            // initialize directions_weights
            // test the possibility of the step (set to NaN if position is not
            // accessible)
            kernel().space_manager.sense(directions_weights, filopodia_,
                                         geometry_.position, move_,
                                         move_.module, 1., nan(""));
            all_nan = allnan(directions_weights);
            if (all_nan)
            {
                if (abs(move_.sigma_angle - M_PI) < 1e-6)
                {
                // completely stuck
#ifndef NDEBUG
                    printf("\nNeurite completely stucked\n");
#endif
                    stuck_ = true;
                }
                else
                {
                    // increase sensing angle and reset weights
                    move_.sigma_angle = std::min(2 * move_.sigma_angle, M_PI);
                    for (int i = 0; i < directions_weights.size(); i++)
                    {
                        directions_weights[i] = 1.;
                    }
                }
            }
            else
            {
                // test the presence of walls to which the growth cone could
                // be attracted
                kernel().space_manager.sense(
                    directions_weights, filopodia_, geometry_.position, move_,
                    filopodia_.finger_length, 1., filopodia_.wall_affinity);

                // check stronger interaction if filopodia is long enough
                if (0.5 * filopodia_.finger_length > move_.module)
                {
                    kernel().space_manager.sense(
                        directions_weights, filopodia_, geometry_.position,
                        move_, filopodia_.finger_length * 0.5, 1.,
                        filopodia_.wall_affinity * 2);
                }
            }
        }
    }
    else
    {
        // NB this is a temporary fix while neurite-neurite interactions are
        // not implemented, because of the way directions_weights is initialized
        // MUST be removed afterwards, investigate also change in the way
        // directions_weights is initialized then remove
        directions_weights =
            std::vector<double>(filopodia_.normal_weights.size(), 1.);
    }
}


/**
 * @brief Include the intrinsic probability from the growth cone mechanics
 *
 * Multiply the stepping probability obtained from the environment by the
 * intrinsic probability of the growth cone, based on its rigidity (it
 * prefers going straight in its current direction because of the elastic
 * properties of the cytoskeleton).
 */
void GrowthCone::compute_intrinsic_direction(
    std::vector<double> &directions_weights)
{
#ifndef NDEBUG
    double sum = 0.;
#endif
    // if we're using the environment, this will pick the new direction,
    // otherwise filopodia is empty so this will be skipped
    for (int n = 0; n < filopodia_.size; n++)
    {
        directions_weights[n] *= filopodia_.normal_weights[n];
#ifndef NDEBUG
        if (not std::isnan(directions_weights[n]))
        {
            sum += directions_weights[n];
        }
#endif
    }
    assert(sum > 0 && not std::isnan(sum));
}


/**
 * @brief Choose the direction that will exerts strongest pull at current step
 *
 * From the likelihood distribution, pick one direction at random that will
 * represent the direction exerting the strongest pull during the current step.
 * Note that though it is likely to be one of the direction where the affinity
 * is highest, this is not necessarily the case.
 */
void GrowthCone::choose_pull_direction(std::vector<double> &directions_weights,
                                       mtPtr rnd_engine)
{
    double sum = 0;
    assert(directions_weights.size() == filopodia_.size);
    assert(filopodia_.directions.size() == filopodia_.size);

    for (auto weight : directions_weights)
    {
        if (not std::isnan(weight))
        {
            sum += weight;
        }
    }

    double x = uniform_(*(rnd_engine).get()) * sum;

    assert(sum > 0 && not std::isnan(sum));
    int n         = 0;
    double tot    = 0.;
    double weight = 0.;

    for (n = 0; n < filopodia_.size; n++)
    {
        weight = directions_weights[n];
        if (not std::isnan(weight))
        {
            tot += weight;
            if (x < tot)
            {
                delta_angle_ = filopodia_.directions[n];
                break;
            }
        }
    }
}


/**
 * @brief Compute the new step direction from strongest pull direction
 *
 * From the direction exerting the strongest pull (delta_angle_), we update
 * the angle describing the growth cone direction.
 */
void GrowthCone::compute_new_direction(mtPtr rnd_engine, double substep)
{
    move_.angle += move_.sigma_angle * delta_angle_;
}


double GrowthCone::init_filopodia()
{
    // here we are creating an array of angles in the range [-4,4]
    // this array is the base to compute the gaussian distribution,
    // it's computed as each point is equidistant from the other
    // and span from [-4,4]
    double bin = 4. / (filopodia_.size / 2. - 0.5);
    filopodia_.directions.clear();

#ifndef NDEBUG
    printf("spanning angle: sigma * [-4,4] \n"
           "with resolution: sigma * %f \n",
           bin);
#endif

    for (int n_angle = 0; n_angle < filopodia_.size; n_angle++)
    {
        double angle = bin * (n_angle + (1 - filopodia_.size) / 2.);
        filopodia_.directions.push_back(angle);
    }

    // printf("\n end angle array \n");
    // this function will create a normal distribution N(0,1)
    // uniformly sampled over filopodia_.size number
    // in the range  [-4,4] which is 4th-sigma range and contains the
    // 0.99993665751 of probability
    // the sampling is computed as the integral of the polygonal chain
    // approximating the curve
    double half_bin = bin / 2.;
    filopodia_.normal_weights.clear();
    for (int n_angle = 0; n_angle < filopodia_.size; n_angle++)
    {
        double x = filopodia_.directions[n_angle];
        double P = SQRT_FRAC_1_2PI * bin *
                   (0.5 * exp(-0.5 * (x - half_bin) * (x - half_bin)) +
                    0.5 * exp(-0.5 * (x + half_bin) * (x + half_bin)));
        filopodia_.normal_weights.push_back(P);
    }
}


// ###########################################################
//              Interface functions
// ###########################################################


double GrowthCone::compute_CR_demand(mtPtr rnd_engine)
{
    return 0;
}


void GrowthCone::compute_speed(mtPtr rnd_engine, double substep)
{
    if (speed_variance_ > 0)
    {
        move_.speed = speed_growth_cone_
                    + speed_variance_ * sqrt(substep)
                                      * normal_(*(rnd_engine).get());
    }
}


/**
 * @brief  It's a standard constructor which allows to inherit the model from
 * another 'GrowtCone'.
 *
 *
 * @param parent
 * @param neurite
 * @param distanceToParent
 * @param binaryID
 * @param position
 * @param angle
 *
 * @return
 */
GCPtr GrowthCone::clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle)
{
    auto newCone = std::make_shared<GrowthCone>(*this);
    newCone->update_topology(parent, neurite, distanceToParent, binaryID,
                             position, angle);
    return newCone;
}


void GrowthCone::prepare_for_split() {}


void GrowthCone::after_split() {}


double GrowthCone::get_growth_cone_speed() const { return 0; }


double GrowthCone::get_CR_received() const { return -1; }


double GrowthCone::get_CR_left() const { return -1; }


double GrowthCone::get_CR_used() const { return -1; }


void GrowthCone::reset_CR_demand() {}


double GrowthCone::get_CR_speed_factor() const { return -1; }


double GrowthCone::get_CR_topo_coeff() const { return 0; }


double GrowthCone::get_module() const { return move_.module; }


void GrowthCone::set_angle(double angle) { move_.angle = angle; }


void GrowthCone::set_cone_ID()
{
    gc_ID_ = biology_.own_neurite->get_and_increment_gc_ID();
}


size_t GrowthCone::get_cone_ID() const { return gc_ID_; }


// ###########################################################
//              Set Get status
// ###########################################################

void GrowthCone::set_status(const statusMap &status)
{

    get_param(status, names::filopodia_wall_affinity, filopodia_.wall_affinity);
    get_param(status, names::filopodia_finger_length, filopodia_.finger_length);
    assert(filopodia_.finger_length > 0);
    get_param(status, names::filopodia_angular_resolution, filopodia_.size);
    assert(filopodia_.size > 0);

#ifndef NDEBUG
    printf("\n"
           "ENVIRONMENT INTERACTION \n");
#endif
    init_filopodia();

    // other models can set the average speed
    get_param(status, names::speed_growth_cone, speed_growth_cone_);
    average_speed_ = speed_growth_cone_;
    move_.speed    = speed_growth_cone_;

    assert(average_speed_ > 0);

    get_param(status, names::speed_variance, speed_variance_);

    if (speed_variance_ < 0)
    {
        throw std::runtime_error("`speed_variance` must be positive.");
    }

    get_param(status, names::rw_sensing_angle, rw_sensing_angle_);
    move_.sigma_angle = rw_sensing_angle_;

#ifndef NDEBUG
    printf("angular resolution: %i \n"
           "finger length:  %f\n"
           "wall affinity:  %f\n",
           filopodia_.size, filopodia_.finger_length, filopodia_.wall_affinity);
#endif
}


void GrowthCone::get_status(statusMap &status) const
{
    set_param(status, names::filopodia_wall_affinity, filopodia_.wall_affinity);
    set_param(status, names::filopodia_finger_length, filopodia_.finger_length);
    set_param(status, names::filopodia_angular_resolution, filopodia_.size);
    set_param(status, names::speed_growth_cone, speed_growth_cone_);
    set_param(status, names::rw_sensing_angle, rw_sensing_angle_);

    // update observables
    std::vector<std::string> tmp;
    get_param(status, names::observables, tmp);
    // use set to keep only unique values
    std::set<std::string> obs(tmp.cbegin(), tmp.cend());
    obs.insert(observables_.cbegin(), observables_.cend());
    // set the full vector and update status
    tmp = std::vector<std::string>(obs.begin(), obs.end());
    set_param(status, names::observables, tmp);
}


/**
 * @brief Get the current value of one of the observables
 */
double GrowthCone::get_state(const char* observable) const
{
    double value = 0.;

    TRIE(observable)
    CASE("length")
        value = biology_.branch->get_distance_to_soma();
    CASE("speed")
        value = move_.speed;
    CASE("angle")
        value = move_.angle;
    ENDTRIE;

    return value;
}


void GrowthCone::update_kernel_variables()
{
    using_environment_ = kernel().using_environment();
    timestep_          = kernel().simulation_manager.get_resolution();
}

} // namespace
