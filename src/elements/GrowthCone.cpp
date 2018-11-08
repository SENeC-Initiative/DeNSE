#include "GrowthCone.hpp"

// C++ includes
#define _USE_MATH_DEFINES
#include <cassert>
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


GrowthCone::GrowthCone(const std::string &model)
    : TopologicalNode()
    , model_(model)
    , observables_({"length", "speed", "angle", "retraction_time", "stopped"})
    , total_proba_(0.)
    , delta_angle_(0)
    , stuck_(false)
    , stopped_(false)
    , interacting_(false)
    , turning_(0)
    , turned_(0.)
    , update_filopodia_(false)
    , max_stop_(10)
    , current_stop_(0)
    , filopodia_{{},
                 {},
                 FILOPODIA_MIN_NUM,
                 FILOPODIA_FINGER_LENGTH,
                 FILOPODIA_SUBSTRATE_AFINITY,
                 FILOPODIA_WALL_AFFINITY}
    , num_filopodia_(FILOPODIA_MIN_NUM)
    , min_filopodia_(FILOPODIA_MIN_NUM)
    , move_()
    , avg_speed_(SPEED_GROWTH_CONE)
    , local_avg_speed_(SPEED_GROWTH_CONE)
    , speed_variance_(0)
    , local_speed_variance_(0)
    , sensing_angle_(SENSING_ANGLE)
    , current_area_("")
    , duration_retraction_(DURATION_RETRACTION)
    , proba_retraction_(PROBA_RETRACTION)
    , retracting_todo_(-1)
    , speed_ratio_retraction_(SPEED_RATIO_RETRACTION)
    , proba_down_move_(PROBA_DOWN_MOVE)
    , max_sensing_angle_(MAX_SENSING_ANGLE)
    , scale_up_move_(SCALE_UP_MOVE)
{
    // random distributions
    normal_      = std::normal_distribution<double>(0, 1);
    uniform_     = std::uniform_real_distribution<double>(0., 1.);
    exponential_ = std::exponential_distribution<double>(proba_retraction_);

    update_kernel_variables();

    // initialize speed (sensing is initialized by filopodia)
    move_.speed       = local_avg_speed_;
    current_diameter_ = get_diameter();

    // only 'initial models' are directly created, other are cloned, so
    // we don't need to get the area.
    update_growth_properties(current_area_);

    init_filopodia();
}


GrowthCone::GrowthCone(const GrowthCone &copy)
    : TopologicalNode(copy)
    , model_(copy.model_)
    , observables_(copy.observables_)
    , total_proba_(0.)
    , delta_angle_(0)
    , current_area_(copy.current_area_)
    , stuck_(false)
    , stopped_(false)
    , turning_(0)
    , turned_(0.)
    , max_stop_(copy.max_stop_)
    , current_stop_(0)
    , interacting_(false)
    , update_filopodia_(false)
    , filopodia_(copy.filopodia_)
    , move_(copy.move_)
    , avg_speed_(copy.avg_speed_)
    , local_avg_speed_(copy.local_avg_speed_)
    , speed_variance_(copy.speed_variance_)
    , local_speed_variance_(copy.local_speed_variance_)
    , sensing_angle_(copy.sensing_angle_)
    , scale_up_move_(copy.scale_up_move_)
    , duration_retraction_(copy.duration_retraction_)
    , proba_retraction_(copy.proba_retraction_)
    , retracting_todo_(-1)
    , speed_ratio_retraction_(copy.speed_ratio_retraction_)
    , proba_down_move_(copy.proba_down_move_)
    , max_sensing_angle_(copy.max_sensing_angle_)
    , min_filopodia_(copy.min_filopodia_)
    , num_filopodia_(copy.num_filopodia_)
    , current_diameter_(copy.current_diameter_)
    , persistence_length_(copy.persistence_length_)
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

    if (using_environment_)
    {
        newCone->update_growth_properties(current_area_);
    }

    return newCone;
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
        own_neurite, 0);
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
    double total_distance = 0.;
    // useful values
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    double current_time(0.), local_substep(substep), tmp, old_substep;
    double angle_widening, unstuck_angle, waiting_time;
    size_t n_direction, i;

    std::vector<double> directions_weights;
    std::vector<std::string> new_pos_area;
    std::vector<bool> wall_presence;

    // we make local substeps until current time is equal to substep
    int loop_index = 0;

    while (current_time < substep)
    {
        loop_index++;
        if (loop_index > 100)
        {
            printf("long loop for gc of neuron %lu on %i - substep %f - ctime %f - next %f\n", biology_.own_neurite->get_parent_neuron().lock()->get_gid(), omp_id, substep, current_time, local_substep);
            throw std::runtime_error("Neuron stuck in an infinite loop");
        }

        local_substep = substep - current_time;

        // reset stopped status
        stopped_ = false;

        // initialize direction test variables
        directions_weights = std::vector<double>(filopodia_.size, 1.);
        wall_presence      = std::vector<bool>(filopodia_.size, false);
        new_pos_area = std::vector<std::string>(filopodia_.size, current_area_);

        // =================================================== //
        // Are we retracting because we were previously stuck? //
        // =================================================== //

        // test stuck/stopped-induced retraction (compute retraction time if
        // necessary)
        if (retracting_todo_ > 0)
        {
            // check for end of retraction during local step
            if (retracting_todo_ <= local_substep)
            {
                local_substep    = retracting_todo_;
                retracting_todo_ = -1;
            }
            else
            {
                retracting_todo_ -= local_substep;
            }

            // compute speed and module
            compute_speed(rnd_engine, local_substep);
            compute_module(local_substep);

            // retract
            retraction(std::abs(move_.module), cone_n, omp_id);
        }
        else
        {
            // ============================= //
            // Extending, compute local pull //
            // ============================= //

            // check environment and mechanical interactions
            try
            {
                interacting_ = compute_pull(directions_weights, wall_presence,
                                            substep, rnd_engine);
            }
            catch (const std::exception &except)
            {
                std::throw_with_nested(std::runtime_error(
                    "Passed from `GrowthCone::compute_pull`."));
            }

            // if interacting with obstacles, and adaptive timestep, switch
            // local_substep down
            if (interacting_ and adaptive_timestep_ > 1)
            {
                double substep_tmp = 1.;

                if (local_substep > adaptive_timestep_)
                {
                    substep_tmp =
                        std::max(timestep_divider_ * local_substep, 1.);
                }

                // always check remaining time >= 1.
                if (substep - current_time - substep_tmp >= 1.)
                {
                    local_substep = substep_tmp;
                }
            }

            // compute speed and module
            compute_speed(rnd_engine, local_substep);
            compute_module(local_substep);

            // Assess current situation and decide on following move base on
            // that. Check and change the sensing_angle if necessary
            // update current_time and set next local_substep
            if (move_.module > 0)
            {
                // ========================= //
                // Forward! Can we go there? //
                // ========================= //

                if (using_environment_)
                {
                    try
                    {
                        compute_accessibility(directions_weights, new_pos_area,
                                              local_substep);
                    }
                    catch (const std::exception &except)
                    {
                        std::throw_with_nested(std::runtime_error(
                            "Passed from "
                            "`GrowthCone::compute_accessibility`."));
                    }
                }

                // check for stuck/total_proba_
                // take optional GC rigidity into account
                compute_intrinsic_direction(directions_weights, local_substep);

                if (not stuck_)
                {
                    // ========================= //
                    // Yes, making forward move! //
                    // ========================= //

                    // set delta_angle_ as main pulling direction
                    try
                    {
                        make_move(directions_weights, new_pos_area,
                                  local_substep, rnd_engine, omp_id);
                    }
                    catch (const std::exception &except)
                    {
                        std::throw_with_nested(std::runtime_error(
                            "Passed from `GrowthCone::make_move`."));
                    }

                    // assess stopped state (computed in make_move)
                    if (not stopped_)
                    {
                        // reset retraction
                        retraction_time_ = -1;
                        // reset turning
                        turning_ = 0;
                        turned_  = 0.;
                    }
                }

                if ((stuck_ or stopped_) and turning_ == 0)
                {
                    // if we are stuck/stopped, we will start widening the
                    // sensing angle and turning in a random direction so
                    // we choose it
                    turning_ = uniform_(*rnd_engine.get()) > 0.5 ? 1 : -1;
                }
            }
            else
            {
                // =========== //
                // Back we go! //
                // =========== //

                // retracting distance is the opposite of the (negative module)
                retraction(-move_.module, cone_n, omp_id);
            }
        }

        // ============================= //
        // Check for stalled growth cone //
        // ============================= //

        if (using_environment_)
        {
            if (stuck_)
            {
                // ============================== //
                // We're stuck. Widen or retract? //
                // ============================== //

                // compute the angle widening necessary to unstuck
                unstuck_angle = kernel().space_manager.unstuck_angle(
                    geometry_.position, move_.angle, filopodia_.finger_length,
                    current_area_, omp_id);

                // compute the time necessary to reach that angle
                local_substep = std::max(
                    (unstuck_angle / 4. - move_.sigma_angle) / ONE_DEGREE, 1.);

                if (unstuck_angle > 0.5 * max_sensing_angle_)
                {
                    local_substep = substep - current_time;
                }
                else if (local_substep > substep - current_time)
                {
                    local_substep = substep - current_time;
                }

                // check whether we would retract in that time
                local_substep = check_retraction(local_substep, rnd_engine);

                // we do not retract, widen angle to unstuck
                change_sensing_angle(ONE_DEGREE * local_substep);
            }
            else if (stopped_)
            {
                // check whether we would retract in that time
                local_substep = check_retraction(local_substep, rnd_engine);
            }
            else
            {
                // reset retraction_time_
                retraction_time_ = -1.;
            }

            if (total_proba_ < 1. and not stuck_)
            {
                // widen angle depending on total_proba_
                angle_widening = 1. - total_proba_;
                change_sensing_angle(ONE_DEGREE * angle_widening *
                                     local_substep);
            }
            else
            {
                // bring move_.sigma_angle towards its default value based on
                // what was done during previous step
                change_sensing_angle(-local_substep * ONE_DEGREE *
                                     local_substep);
            }
        }

        // update current time
        current_time += local_substep;

        total_distance += move_.module;

        // update retraction time if necessary
        if (retraction_time_ > 0.)
        {
            retraction_time_ -= local_substep;
        }
    }
}


void GrowthCone::compute_module(double substep)
{
    move_.module = move_.speed * substep;
}


double GrowthCone::check_retraction(double substep, mtPtr rnd_engine)
{
    if (retraction_time_ < 0)
    {
        retraction_time_ = 1. + exponential_(*rnd_engine.get());
    }
    //~ printf("checking retraction: time %f vs substep %f\n", retraction_time_, substep);

    if (retraction_time_ <= substep)
    {
        // we will retract, subtract the time we waited
        substep = substep - retraction_time_;
        // set retraction duration
        retracting_todo_ = duration_retraction_;
        // reset retraction time
        retraction_time_ = -1;
    }

    assert(substep >= 0);

    return substep;
}


void GrowthCone::retraction(double distance, size_t cone_n, int omp_id)
{
    assert(distance >= 0.);

    // cannot be stuck_ or on low proba mode when retracting, so reset all
    stuck_       = false;
    total_proba_ = filopodia_.size;
    // also reset turning
    turning_ = 0;
    turned_  = 0.;

    // remove the points
    double distance_done;
    //~ printf("distance to retract %f vs length %f\n", distance, get_branch()->get_length());

    while (distance > 0)
    {
        if (biology_.branch->size() > 1)
        {
            distance_done = biology_.branch->points[2].back() -
                       biology_.branch->points[2][biology_.branch->size() - 2];

            if (distance_done < distance)
            {
                distance -= distance_done;
                biology_.branch->retract();
            }
            else
            {
                double remaining = distance_done - distance;
                Point p1 = biology_.branch->xy_at(biology_.branch->size() - 2);
                Point p2 = biology_.branch->get_last_xy();

                double new_x =
                    (p2[0] * remaining + p1[0] * (distance_done - remaining)) / distance_done;
                double new_y =
                    (p2[1] * remaining + p1[1] * (distance_done - remaining)) / distance_done;
                Point new_p = Point(new_x, new_y);

                biology_.branch->retract();
                biology_.branch->add_point(new_p, remaining);

                distance = 0.;
            }
        }
        else
        {
            break;
        }
    }

    // set the new growth cone angle
    auto points = biology_.branch->points;
    size_t last = points[0].size();
    double x0, y0, x1, y1;

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
            Point p = TopologicalNode::get_position();
            x0      = p[0];
            y0      = p[1];
        }

        move_.angle = atan2(y1 - y0, x1 - x0);
    }

    geometry_.position = biology_.branch->get_last_xy();

    // prune growth cone if necessary
    if (biology_.branch->size() == 1)
    {
        prune(cone_n);
    }

    // check if we changed area
    if (using_environment_)
    {
        current_area_ = kernel().space_manager.get_containing_area(
            geometry_.position, omp_id);
        update_growth_properties(current_area_);
        // reset move_.sigma_angle to its default value
        AreaPtr area = kernel().space_manager.get_area(current_area_);
        move_.sigma_angle =
            sensing_angle_ * area->get_property(names::sensing_angle);
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
 * Scan surrounding environment and assess how likely to be chosen each
 * direction is, depending on how strongly it pulls on the GC.
 *
 * @return distance between GC and closest obstacle
 */
bool GrowthCone::compute_pull(std::vector<double> &directions_weights,
                              std::vector<bool> &wall_presence, double substep,
                              mtPtr rnd_engine)
{
    if (using_environment_)
    {
        double up_move = scale_up_move_ == 0 ? std::nan("") : scale_up_move_;

        return kernel().space_manager.sense(
            directions_weights, wall_presence, filopodia_, geometry_.position,
            move_, current_area_, proba_down_move_, up_move, substep,
            sqrt_resol_, num_filopodia_ - min_filopodia_);
    }

    return false;
}


/**
 * @brief compute the possibility of the next move.
 *
 * Scan surrounding environment and assess whether a step in each of the
 * possible directions is possible given the desired step length.
 */
void GrowthCone::compute_accessibility(std::vector<double> &directions_weights,
                                       std::vector<std::string> &new_pos_area,
                                       double substep)
{
    // unstuck neuron
    stuck_ = false;

    // test the possibility of the step (set to NaN if position is not
    // accessible)
    kernel().space_manager.move_possibility(
        directions_weights, new_pos_area, filopodia_, geometry_.position, move_,
        substep, sqrt_resol_, num_filopodia_ - min_filopodia_);
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
    std::vector<double> &directions_weights, double substep)
{
    total_proba_ = 0.;
    stuck_       = true;

    double weight;

    //~ for (unsigned int n = 0; n < filopodia_.size; n++)
    for (unsigned int n = 0; n < directions_weights.size(); n++)
    {
        weight = directions_weights[n] * filopodia_.normal_weights[n];
        directions_weights[n] = weight;

        if (not std::isnan(weight))
        {
            total_proba_ += weight;
            stuck_ = false;
        }
    }
}


/**
 * @brief Choose the direction that will exerts strongest pull at current step
 *
 * From the likelihood distribution, pick one direction at random that will
 * represent the direction exerting the strongest pull during the current step.
 * Note that though it is likely to be one of the direction where the affinity
 * is highest, this is not necessarily the case.
 */
void GrowthCone::make_move(const std::vector<double> &directions_weights,
                           const std::vector<std::string> &new_pos_area,
                           double &substep, mtPtr rnd_engine, int omp_id)
{
    assert(directions_weights.size() == filopodia_.size);
    assert(filopodia_.directions.size() == filopodia_.size);

    double x = uniform_(*(rnd_engine).get()) * total_proba_;

    double bin = max_sensing_angle_ / filopodia_.size;

    int n         = -1;                       // index of chosen filopodia
    double tot    = 0.;                       // cumulant of the probability
    double frac   = 0.;                       // fraction of current
    double weight = directions_weights.at(0); // initial weight
    double new_angle, default_angle, mean;
    GeomPtr line;
    Point p;

    // check whether we're stopped
    stopped_ = (uniform_(*(rnd_engine.get())) > total_proba_);

    // we test whether we should make the next move or stop moving if moving
    // anywhere is too unlikely
    if (not stopped_)
    {
        for (n = 0; n < filopodia_.size; n++)
        {
            weight = directions_weights.at(n);

            if (not std::isnan(weight))
            {
                tot += weight;

                if (x < tot)
                {
                    frac = tot - x;

                    delta_angle_ =
                        filopodia_.directions[0] + bin * (n + frac / weight);

                    // compute the target position base on the model
                    p = compute_target_position(directions_weights, rnd_engine,
                                                substep, new_angle);

                    line = kernel().space_manager.geosline_from_points(
                        geometry_.position, p);

                    // set the new angle (chose a valid one if target position
                    // is outside the environment)
                    if (using_environment_)
                    {
                        default_angle = move_.angle + filopodia_.directions[n];

                        // move towards default_angle
                        while (
                            kernel().space_manager.env_intersects(line, omp_id))
                        {
                            new_angle = 0.5 * (new_angle + default_angle);
                            p         = Point(geometry_.position.at(0) +
                                          cos(new_angle) * move_.module,
                                      geometry_.position.at(1) +
                                          sin(new_angle) * move_.module);
                            line = kernel().space_manager.geosline_from_points(
                                geometry_.position, p);
                        }

                        // check if we switched to a new area
                        std::string new_area =
                            kernel().space_manager.get_containing_area(
                                geometry_.position, omp_id);

                        if (new_area != current_area_)
                        {
                            update_growth_properties(new_area);
                        }
                    }

                    // update angle
                    delta_angle_ = new_angle - move_.angle;
                    move_.angle  = new_angle;

                    // store new position
                    geometry_.position = p;
                    biology_.branch->add_point(geometry_.position,
                                               move_.module);

                    break;
                }
            }
        }
    }

    // check for unsucessful move
    if (tot == total_proba_)
    {
        stopped_ = true;
    }
}


/**
 * @brief Compute the new step direction from strongest pull direction
 *
 * From the direction exerting the strongest pull (delta_angle_), we update
 * the angle describing the growth cone direction and set the new position
 */
Point GrowthCone::compute_target_position(
    const std::vector<double> &directions_weights, mtPtr rnd_engine,
    double &substep, double &new_angle)
{
    new_angle = move_.angle + delta_angle_;
    Point target_pos =
        Point(geometry_.position.at(0) + cos(new_angle) * move_.module,
              geometry_.position.at(1) + sin(new_angle) * move_.module);

    return target_pos;
}


/**
 * Init the filopodia's kernel
 *
 * @warning must always be preceeded by a call to update_growth_properties
 */
void GrowthCone::init_filopodia()
{
    double dtheta, std_norm, proba_norm, bin, angle, P;

    double resol = kernel().simulation_manager.get_resolution();

    // special case for the simple random_walk model (angle varies with the
    // resolution, hence so does the number of filopodia)
    if (model_ == "simple_random_walk")
    {
        // update max sensing angle depending on resolution
        //! IMPORTANT: MUST COME BEFORE dtheta
        double old_max = max_sensing_angle_;
        max_sensing_angle_ =
            std::min(max_sensing_angle_ * sqrt_resol_, 2 * M_PI);

        // we add three quarters of what's needed, this is the lowest increase
        // giving correct results for the random walk behavior
        bin = max_sensing_angle_ / (min_filopodia_ - 1);
        unsigned int add_filo =
            static_cast<int>((max_sensing_angle_ - old_max) / bin);

        //! IMPORTANT: MUST COME BEFORE dtheta
        num_filopodia_     = min_filopodia_ + add_filo;
        max_sensing_angle_ = (num_filopodia_ - 1) * bin;

        // update filopodia number
        dtheta     = max_sensing_angle_ / num_filopodia_;
        std_norm   = 0.5 / (move_.sigma_angle * move_.sigma_angle);
        proba_norm = SQRT_FRAC_1_2PI * dtheta / (move_.sigma_angle);

        filopodia_.directions     = std::vector<double>(num_filopodia_);
        filopodia_.normal_weights = std::vector<double>(num_filopodia_);

        // set the angles and the weights
        for (int n_angle = 0; n_angle < num_filopodia_; n_angle++)
        {
            // fill directions
            angle = bin * n_angle - 0.5 * max_sensing_angle_;
            filopodia_.directions[n_angle] = angle;
            // set normal weights
            angle = filopodia_.directions[0] + n_angle * dtheta;
            P     = proba_norm *
                (exp(-std_norm * angle * angle) +
                 exp(-std_norm * (angle + dtheta) * (angle + dtheta)));
            filopodia_.normal_weights[n_angle] = P;
        }

        filopodia_.size = num_filopodia_;
    }
    else
    {
        // move sensing angle must have been updated by previous call to
        // update_growth_properties

        // compute angles
        num_filopodia_ = min_filopodia_;
        dtheta         = max_sensing_angle_ / num_filopodia_;
        std_norm       = 0.5 / (move_.sigma_angle * move_.sigma_angle);
        proba_norm     = SQRT_FRAC_1_2PI * dtheta / move_.sigma_angle;
        bin            = max_sensing_angle_ / (min_filopodia_ - 1);

        filopodia_.directions     = std::vector<double>(num_filopodia_);
        filopodia_.normal_weights = std::vector<double>(num_filopodia_);

        // set the angles
        for (int n_angle = 0; n_angle < num_filopodia_; n_angle++)
        {
            angle = bin * n_angle - 0.5 * max_sensing_angle_;
            filopodia_.directions[n_angle] = angle;
        }

        filopodia_.size = num_filopodia_;

        // update the filopodia normal weights
        update_filopodia(resol);
    }
}


// ###########################################################
//              Interface functions
// ###########################################################

void GrowthCone::compute_speed(mtPtr rnd_engine, double substep)
{
    if (speed_variance_ > 0)
    {
        move_.speed = local_avg_speed_ + local_speed_variance_ * sqrt(substep) *
                                             normal_(*(rnd_engine).get());
    }
    else
    {
        move_.speed = local_avg_speed_;
    }
}


void GrowthCone::prepare_for_split() {}


void GrowthCone::after_split() {}


double GrowthCone::get_growth_cone_speed() const { return move_.speed; }


double GrowthCone::get_module() const { return move_.module; }


void GrowthCone::set_angle(double angle) { move_.angle = angle; }


void GrowthCone::set_diameter(double diameter)
{
    current_diameter_ = diameter;
}


// ###########################################################
//              Set Get status
// ###########################################################

void GrowthCone::set_status(const statusMap &status)
{

    get_param(status, names::filopodia_wall_affinity, filopodia_.wall_affinity);
    get_param(status, names::filopodia_finger_length, filopodia_.finger_length);
    assert(filopodia_.finger_length > 0);
    get_param(status, names::scale_up_move, scale_up_move_);
    bool update_filo =
        get_param(status, names::filopodia_min_number, min_filopodia_);

    if (min_filopodia_ < 10)
    {
        throw std::invalid_argument(
            "`filopodia_min_number` should be at least 10 to ensure physically "
            "relevant behaviors.");
    }

    // other models can set the average speed
    bool speed= get_param(status, names::speed_growth_cone, avg_speed_);

    get_param(status, names::speed_variance, speed_variance_);
    if (speed_variance_ < 0)
    {
        throw std::invalid_argument("`speed_variance` must be positive.");
    }

    get_param(status, names::proba_retraction, proba_retraction_);

    // update exponential distribution
    exponential_ = std::exponential_distribution<double>(proba_retraction_);

    get_param(status, names::duration_retraction, duration_retraction_);
    get_param(status, names::speed_ratio_retraction, speed_ratio_retraction_);
    get_param(status, names::proba_down_move, proba_down_move_);

    // sensing angle and persistence length
    double sa_tmp, lp, old_lp(persistence_length_);
    bool set_lp = get_param(status, names::persistence_length, lp);
    bool set_sa = get_param(status, names::sensing_angle, sa_tmp);

    if (model_ == "simple_random_walk")
    {
        if (set_sa and set_lp)
        {
            throw std::runtime_error("Cannot set both `sensing_angle` and "
                                     "`persistence_length` in the "
                                     "`simple_random_walk` model.");
        }
        else if (set_sa)
        {
            sensing_angle_     = sa_tmp;
            sensing_angle_set_ = true;
        }
        else
        {
            persistence_length_ = lp;
            sensing_angle_set_  = false;
        }
    }
    else
    {
        sensing_angle_      = sa_tmp;
    }

    update_filo += set_sa;

    double msa = max_sensing_angle_;
    update_filo += get_param(status, names::max_sensing_angle, msa);

    if (msa < 0.25 * M_PI)
    {
        throw std::invalid_argument("`max_sensing_angle` must be greater "
                                    "than 90°.");
    }
    max_sensing_angle_ = std::min(msa, 2 * M_PI);

    // set growth properties
    if (using_environment_)
    {
        update_growth_properties(current_area_);
    }
    else
    {
        local_avg_speed_      = avg_speed_;
        local_speed_variance_ = speed_variance_;
    }

    // persistence length
    // IMPORTANT: MUST COME AFTER local_avg_speed_ HAS BEEN SET!
    if (set_lp)
    {
        if (lp < resol_ * local_avg_speed_)
        {
            persistence_length_ = old_lp;

            throw InvalidParameter(
                names::persistence_length, std::to_string(lp),
                ">= " + std::to_string(resol_ * local_avg_speed_) +
                    " (the minimal step)",
                __FUNCTION__, __FILE__, __LINE__);
        }

        if (model_ == "simple_random_walk")
        {
            sensing_angle_ =
                std::sqrt(2 * local_avg_speed_ * resol_ / persistence_length_);
            move_.sigma_angle = sensing_angle_;
        }

        persistence_length_ = lp;
    }

    // reset filopodia parameters
    if (update_filo)
    {
        init_filopodia();
    }
}


void GrowthCone::get_status(statusMap &status) const
{
    set_param(status, names::filopodia_wall_affinity, filopodia_.wall_affinity,"");
    set_param(status, names::filopodia_finger_length, filopodia_.finger_length, "micrometer");
    set_param(status, names::filopodia_min_number, min_filopodia_, "");
    set_param(status, names::speed_growth_cone, avg_speed_, "micrometer / minute");

    // @todo change behavior of sensing_angle and max_sensing angle:
    // sensing angle set the min/max of the filopodia angle and is limited by
    // max_sensing_angle
    set_param(status, names::sensing_angle, sensing_angle_, "rad");
    set_param(status, names::max_sensing_angle, max_sensing_angle_, "rad");

    set_param(status, names::scale_up_move, scale_up_move_, "micrometer");
    set_param(status, names::proba_down_move, proba_down_move_, "");

    set_param(status, names::persistence_length, persistence_length_, "micrometer");

    set_param(status, names::proba_retraction, proba_retraction_, "");
    set_param(status, names::duration_retraction, duration_retraction_, "minute");

    // update observables
    std::vector<std::string> tmp;
    get_param(status, names::observables, tmp);
    // use set to keep only unique values
    std::set<std::string> obs(tmp.cbegin(), tmp.cend());
    obs.insert(observables_.cbegin(), observables_.cend());
    // set the full vector and update status
    tmp = std::vector<std::string>(obs.begin(), obs.end());
    set_param(status, names::observables, tmp, "");
}


/**
 * @brief Get the current value of one of the observables
 */
double GrowthCone::get_state(const char *observable) const
{
    double value = 0.;

    TRIE(observable)
    CASE("length")
    value = biology_.branch->get_length();
    CASE("speed")
    value = move_.speed;
    CASE("angle")
    value = move_.angle;
    CASE("stopped")
    value = 2 * stuck_ + stopped_; // 0: moving, 1: stopped, 2: stuck
    CASE("retraction_time")
    value = retraction_time_;
    CASE("diameter")
    value = current_diameter_;
    ENDTRIE;

    return value;
}


double GrowthCone::get_diameter() const
{
    return current_diameter_;
}


void GrowthCone::update_growth_properties(const std::string &area_name)
{
    // update the growth properties depending on the area
    AreaPtr area = kernel().space_manager.get_area(area_name);

    // test because first "model" growth cones are not in any Area
    if (area != nullptr)
    {
        current_area_ = area_name;
        // speed and sensing angle may vary
        if (model_ == "simple_random_walk" and not sensing_angle_set_)
        {
            sensing_angle_ = std::sqrt(2 * local_avg_speed_ * resol_ /
                                       persistence_length_);
        }
        move_.sigma_angle =
            sensing_angle_ * area->get_property(names::sensing_angle);
        local_avg_speed_ =
            avg_speed_ * area->get_property(names::speed_growth_cone);
        local_speed_variance_ =
            speed_variance_ * area->get_property(names::speed_growth_cone);

        // substrate affinity depends on the area
        filopodia_.substrate_affinity =
            area->get_property(names::substrate_affinity);
    }
    else
    {
        move_.sigma_angle     = sensing_angle_;
        local_avg_speed_      = avg_speed_;
        local_speed_variance_ = speed_variance_;
    }
}


void GrowthCone::update_kernel_variables()
{
    using_environment_ = kernel().using_environment();

    // check change in resolution
    double old_resol = resol_;
    resol_           = kernel().simulation_manager.get_resolution();

    // check adaptive timestep
    adaptive_timestep_ = kernel().get_adaptive_timestep();
    timestep_divider_  = 1. / adaptive_timestep_;

    // check adaptive timestep
    adaptive_timestep_ = kernel().get_adaptive_timestep();
    timestep_divider_  = 1. / adaptive_timestep_;

    // check adaptive timestep
    adaptive_timestep_ = kernel().get_adaptive_timestep();
    timestep_divider_  = 1. / adaptive_timestep_;

    // check if filopodia should be updated
    if (old_resol != resol_)
    {
        sqrt_resol_ = sqrt(resol_);
        int omp_id  = kernel().parallelism_manager.get_thread_local_id();

        update_growth_properties(
            kernel().space_manager.get_containing_area(get_position(), omp_id));

        init_filopodia();
    }
}


void GrowthCone::change_sensing_angle(double angle)
{
    AreaPtr area     = kernel().space_manager.get_area(current_area_);
    double old_sigma = move_.sigma_angle;

    if (angle > 0.)
    {
        move_.sigma_angle =
            std::min(move_.sigma_angle + angle, max_sensing_angle_);

        // when not moving, start turning
        if (turning_ != 0 and turned_ < 0.39269908169872414)
        {
            move_.angle += turning_ * angle;
            turned_ += angle;
        }
    }
    else
    {
        move_.sigma_angle =
            std::max(move_.sigma_angle + angle,
                     sensing_angle_ * area->get_property(names::sensing_angle));
    }

    // change filopodia normal weights if necessary
    if (old_sigma != move_.sigma_angle)
    {
        update_filopodia_ = true;
    }
}


void GrowthCone::update_filopodia(double substep)
{
    assert(model_ != "simple_random_walk");

    double dtheta     = max_sensing_angle_ / filopodia_.size;
    double std_norm   = 0.5 / (move_.sigma_angle * move_.sigma_angle);
    double proba_norm = SQRT_FRAC_1_2PI * dtheta / move_.sigma_angle;

    double angle, P;

    for (unsigned int n = 0; n < num_filopodia_; n++)
    {
        // fill normal weights
        angle = filopodia_.directions[0] + n * dtheta;
        P     = proba_norm * (exp(-std_norm * angle * angle) +
                          exp(-std_norm * (angle + dtheta) * (angle + dtheta)));
        filopodia_.normal_weights[n] = P;
    }
}

} // namespace growth
