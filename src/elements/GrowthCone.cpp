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

// lib include
#include "growth_names.hpp"

// spatial include
#include "Area.hpp"


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
                 FILOPODIA_SUBSTRATE_AFINITY,
                 FILOPODIA_WALL_AFFINITY}
    , move_()
    , avg_speed_(SPEED_GROWTH_CONE)
    , speed_variance_(0)
    , sensing_angle_(SENSING_ANGLE)
    , current_area_("")
    , duration_retraction_(DURATION_RETRACTION)
    , proba_retraction_(PROBA_RETRACTION)
    , retracting_todo_(0)
    , speed_ratio_retraction_(SPEED_RATIO_RETRACTION)
    , proba_down_move_(PROBA_DOWN_MOVE)
    , max_sensing_angle_(MAX_SENSING_ANGLE)
    , scale_up_move_(SCALE_UP_MOVE)
{
    // initialize move variables
    move_.sigma_angle = sensing_angle_;
    move_.speed       = avg_speed_;

    // random distributions
    normal_  = std::normal_distribution<double>(0, 1);
    uniform_ = std::uniform_real_distribution<double>(0., 1.);

    update_kernel_variables();

    // only 'initial models' are directly created, other are cloned, so
    // we don't need to get the area.
    update_growth_properties(current_area_);

    init_filopodia();
}


GrowthCone::GrowthCone(const GrowthCone &copy)
    : TopologicalNode(copy)
    , observables_(copy.observables_)
    , delta_angle_(0)
    , current_area_(copy.current_area_)
    , stuck_(false)
    , filopodia_(copy.filopodia_)
    , move_(copy.move_)
    , avg_speed_(copy.avg_speed_)
    , speed_variance_(copy.speed_variance_)
    , sensing_angle_(copy.sensing_angle_)
    , scale_up_move_(copy.scale_up_move_)
    , duration_retraction_(copy.duration_retraction_)
    , proba_retraction_(copy.proba_retraction_)
    , retracting_todo_(0)
    , speed_ratio_retraction_(copy.speed_ratio_retraction_)
    , proba_down_move_(copy.proba_down_move_)
    , max_sensing_angle_(copy.max_sensing_angle_)
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
                                 const std::string &binaryID,
                                 const Point &position, double angle)
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
    int omp_id = kernel().parallelism_manager.get_thread_local_id();

    compute_speed(rnd_engine, substep);
    compute_module(substep);

    // if the neurite is stuck, there is a given probability that it starts
    // retracting
    if (stuck_)
    {
        if (uniform_(*rnd_engine.get()) < proba_retraction_ * substep)
        {
            retracting_todo_ = duration_retraction_;
        }
    }

    if (retracting_todo_ > 0)
    {
        retracting_todo_ = std::max(0., retracting_todo_ - substep);
        if (move_.module > 0)
        {
            move_.module = -move_.module * speed_ratio_retraction_;
        }
    }

    if (move_.module > 0)
    {
        // check environment and mechanical interactions
        std::vector<double> directions_weights;
        std::vector<std::string> new_pos_area;

        compute_pull_and_accessibility(directions_weights, new_pos_area,
                                       rnd_engine, substep);

        if (not stuck_)
        {
            // weight the directions by the intrinsic preference of the growth
            // cone
            compute_intrinsic_direction(directions_weights);
            // set delta_angle_ as main pulling direction
            size_t n_direction = choose_pull_direction(
                directions_weights, new_pos_area, rnd_engine);

            if (not stuck_)
            {
                // compute new direction based on delta_angle_
                compute_new_direction(rnd_engine, substep);
                // store new position
                geometry_.position = Point(
                    geometry_.position.at(0) + move_.module * cos(move_.angle),
                    geometry_.position.at(1) + move_.module * sin(move_.angle));
                biology_.branch->add_point(geometry_.position, move_.module);
                // check if we switched to a new area
                if (new_pos_area.at(n_direction) != current_area_)
                {
                    current_area_ = kernel().space_manager.get_containing_area(
                        geometry_.position, omp_id);
                    update_growth_properties(current_area_);
                }
                else
                {
                    // reset move_.sigma_angle to its default value
                    AreaPtr area =
                        kernel().space_manager.get_area(current_area_);
                    move_.sigma_angle =
                        sensing_angle_ *
                        area->get_property(names::sensing_angle);
                }
            }
        }
    }
    if (move_.module < 0)
    {
        retraction(move_.module);
        if (biology_.branch->size() == 0)
        {
            prune(cone_n);
        }
        geometry_.position = biology_.branch->get_last_xy();
    }
}


void GrowthCone::compute_module(double substep)
{
    move_.module = move_.speed * substep;
}


void GrowthCone::retraction(double module)
{
    double subtracted = 0;
    double x0, y0, x1, y1;

    while (module < 0)
    {
        if (biology_.branch->size() > 0)
        {
            module += biology_.branch->points[2].back();
            biology_.branch->retract();
        }
        else
        {
            break;
        }
    }

    // set the new growth cone angle
    auto points = biology_.branch->points;
    size_t last = points[0].size();

    if (last > 0)
    {
        x1 = points[0][last - 1];
        y1 = points[1][last - 1];

        if (last > 1)
        {
            x0 = points[0][last - 2];
            y0 = points[1][last - 2];
        }
        else
        {
            Point p = get_parent().lock()->get_position();
            x0      = p[0];
            y0      = p[1];
        }

        move_.angle = atan2(y1 - y0, x1 - x0);
    }
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
    std::vector<double> &directions_weights,
    std::vector<std::string> &new_pos_area, mtPtr rnd_engine, double substep)
{
    if (using_environment_)
    {
        // unstuck neuron
        stuck_ = false;
        //~ if (stuck_)
        //~ {
        //~ if (abs(move_.sigma_angle - 0.75*M_PI) > 1e-6)
        //~ {
        //~ // increase sensing angle and reset weights
        //~ move_.sigma_angle = std::min(
        //~ 1.5 * move_.sigma_angle, 0.75*M_PI);
        //~ }
        //~ stuck_ = false;
        //~ }

        // initialize direction test variables
        bool all_nan = true;

        directions_weights = std::vector<double>(filopodia_.size, 1.);
        new_pos_area = std::vector<std::string>(filopodia_.size, current_area_);

        while (all_nan and not stuck_)
        {
            // test the possibility of the step (set to NaN if position is not
            // accessible)
            kernel().space_manager.sense(
                directions_weights, new_pos_area, filopodia_,
                geometry_.position, move_, move_.module, substep,
                filopodia_.substrate_affinity, nan(""), current_area_,
                proba_down_move_, scale_up_move_);

            all_nan = allnan(directions_weights);

            if (all_nan)
            {
                if (nonmax_sensing_angle())
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
                    widen_sensing_angle();
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
                kernel().space_manager.sense_walls(
                    directions_weights, filopodia_, geometry_.position, move_,
                    filopodia_.finger_length, substep,
                    filopodia_.wall_affinity * substep, current_area_);

                // check stronger interaction if filopodia is long enough
                if (0.5 * filopodia_.finger_length > move_.module)
                {
                    kernel().space_manager.sense_walls(
                        directions_weights, filopodia_, geometry_.position,
                        move_, filopodia_.finger_length * 0.5, substep,
                        filopodia_.wall_affinity * 2 * substep, current_area_);
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
        directions_weights = std::vector<double>(filopodia_.size, 1.);
        new_pos_area = std::vector<std::string>(filopodia_.size, current_area_);
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
    //~ printf("entering GrowthCone::intrinsic directions\n");
    double sum = 0.;
    int count  = 0;
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
            count += 1;
        }
#endif
    }
#ifndef NDEBUG
    if (sum <= 0 || std::isnan(sum))
    {
        printf("invalid sum %f; stuck %u; count %i;\n", sum, stuck_, count);
    }
#endif
    //~ assert(sum > 0 && not std::isnan(sum));
}


/**
 * @brief Choose the direction that will exerts strongest pull at current step
 *
 * From the likelihood distribution, pick one direction at random that will
 * represent the direction exerting the strongest pull during the current step.
 * Note that though it is likely to be one of the direction where the affinity
 * is highest, this is not necessarily the case.
 */
int GrowthCone::choose_pull_direction(
    const std::vector<double> &directions_weights,
    const std::vector<std::string> &new_pos_area, mtPtr rnd_engine)
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
#ifndef NDEBUG
        if (weight < 0)
        {
            printf("negative weight: %f\n", weight);
        }
#endif
    }

    //~ #ifndef NDEBUG
    //~ printf("sum: %f\n", sum);
    //~ #endif

    double x = uniform_(*(rnd_engine).get()) * sum;

    int n         = 0;
    double tot    = 0.;
    double weight = 0.;

    // we test whether we should make the next move or stop moving if moving
    // anywhere is too unlikely
    if (uniform_(*rnd_engine.get()) < sum)
    {
        for (n = 0; n < filopodia_.size; n++)
        {
            weight = directions_weights.at(n);
            if (not std::isnan(weight))
            {
                tot += weight;
                if (x < tot)
                {
                    delta_angle_ = filopodia_.directions.at(n);
                    return n;
                }
            }
        }
    }
    else
    {
        stuck_ = true;
        if (nonmax_sensing_angle())
        {
            widen_sensing_angle();
        }
        return -1;
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
    move_.angle += move_.sigma_angle * sqrt(substep) * delta_angle_;
}


double GrowthCone::init_filopodia()
{
    // here we are creating an array of angles in the range [-4,4]
    // this array is the base to compute the gaussian distribution,
    // it's computed as each point is equidistant from the other
    // and span from [-4,4]
    double bin = 4. / (filopodia_.size / 2. - 0.5);
    filopodia_.directions.clear();

    //~ #ifndef NDEBUG
    //~ printf("spanning angle: sigma * [-4,4] \n"
    //~ "with resolution: sigma * %f \n",
    //~ bin);
    //~ #endif

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
    // total probability is 2 so that a growth cone following a wall still has
    // a 100% probablity of continuing.
    double half_bin = bin / 2.;
    filopodia_.normal_weights.clear();
    for (int n_angle = 0; n_angle < filopodia_.size; n_angle++)
    {
        double x = filopodia_.directions[n_angle];
        double P = SQRT_FRAC_1_2PI * bin *
                   (exp(-0.5 * (x - half_bin) * (x - half_bin)) +
                    exp(-0.5 * (x + half_bin) * (x + half_bin)));
        filopodia_.normal_weights.push_back(P);
    }
}


// ###########################################################
//              Interface functions
// ###########################################################


double GrowthCone::compute_CR_demand(mtPtr rnd_engine) { return 0; }


void GrowthCone::compute_speed(mtPtr rnd_engine, double substep)
{
    if (speed_variance_ > 0)
    {
        move_.speed = avg_speed_ + speed_variance_ * sqrt(substep) *
                                       normal_(*(rnd_engine).get());
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
    int omp_id   = kernel().parallelism_manager.get_thread_local_id();

    // update topological properties
    newCone->update_topology(parent, neurite, distanceToParent, binaryID,
                             position, angle);

    // update containing area
    newCone->current_area_ =
        using_environment_
            ? kernel().space_manager.get_containing_area(position, omp_id)
            : "";

    newCone->update_growth_properties(current_area_);

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
    get_param(status, names::scale_up_move, scale_up_move_);
    get_param(status, names::filopodia_angular_resolution, filopodia_.size);
    assert(filopodia_.size > 0);
    init_filopodia();

    // other models can set the average speed
    get_param(status, names::speed_growth_cone, avg_speed_);

    assert(avg_speed_ > 0);

    get_param(status, names::speed_variance, speed_variance_);

    get_param(status, names::proba_retraction, proba_retraction_);
    get_param(status, names::duration_retraction, duration_retraction_);
    get_param(status, names::speed_ratio_retraction, speed_ratio_retraction_);
    get_param(status, names::proba_down_move, proba_down_move_);

    if (speed_variance_ < 0)
    {
        throw std::runtime_error("`speed_variance` must be positive.");
    }

    get_param(status, names::sensing_angle, sensing_angle_);
    move_.sigma_angle = sensing_angle_; // is it necessary to keep this?

    //~ #ifndef NDEBUG
    //~ printf("angular resolution: %i \n"
    //~ "finger length:  %f\n"
    //~ "wall affinity:  %f\n",
    //~ filopodia_.size, filopodia_.finger_length, filopodia_.wall_affinity);
    //~ #endif

    // set growth properties
    update_growth_properties(current_area_);
}


void GrowthCone::get_status(statusMap &status) const
{
    set_param(status, names::filopodia_wall_affinity, filopodia_.wall_affinity);
    set_param(status, names::filopodia_finger_length, filopodia_.finger_length);
    set_param(status, names::filopodia_angular_resolution, filopodia_.size);
    set_param(status, names::speed_growth_cone, avg_speed_);
    set_param(status, names::sensing_angle, sensing_angle_);

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
double GrowthCone::get_state(const char *observable) const
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


void GrowthCone::update_growth_properties(const std::string &area_name)
{
    // update the growth properties depending on the area
    AreaPtr area = kernel().space_manager.get_area(area_name);

    // test because first "model" growth cones are not in any Area
    if (area != nullptr)
    {
        // set area name
        current_area_ = area_name;
        // speed and sensing angle may vary
        move_.sigma_angle =
            sensing_angle_ * area->get_property(names::sensing_angle);
        move_.speed = avg_speed_ * area->get_property(names::speed_growth_cone);

        // substrate affinity depends on the area
        filopodia_.substrate_affinity =
            area->get_property(names::substrate_affinity);
    }
}


void GrowthCone::update_kernel_variables()
{
    using_environment_ = kernel().using_environment();
}

} // namespace
