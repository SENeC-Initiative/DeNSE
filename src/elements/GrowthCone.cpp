/*
 * GrowthCone.cpp
 *
 * This file is part of DeNSE.
 *
 * Copyright (C) 2019 SeNEC Initiative
 *
 * DeNSE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * DeNSE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DeNSE. If not, see <http://www.gnu.org/licenses/>.
 */

#include "GrowthCone.hpp"

// C++ includes
#define _USE_MATH_DEFINES
#include <algorithm>
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
#include "tools.hpp"

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
    , neuron_id_(0)
    , neurite_name_("")
    , active_(true)
    , model_(model)
    , observables_({"length", "speed", "angle", "retraction_time", "status",
                    "stepping_probability"})
    , total_proba_(0.)
    , delta_angle_(0)
    , stuck_(false)
    , stopped_(false)
    , interacting_(false)
    , just_retracted_(false)
    , turning_(0)
    , turned_(0.)
    , update_filopodia_(false)
    , aff_()
    , max_stop_(50)
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
    , sensing_angle_(SENSING_ANGLE)
    , current_area_("")
    , duration_retraction_(DURATION_RETRACTION)
    , proba_retraction_(PROBA_RETRACTION)
    , retracting_todo_(-1)
    , retraction_time_(-1.)
    , speed_ratio_retraction_(SPEED_RATIO_RETRACTION)
    , proba_down_move_(PROBA_DOWN_MOVE)
    , max_sensing_angle_(MAX_SENSING_ANGLE)
    , scale_up_move_(SCALE_UP_MOVE)
    , old_angle_(0.)
{
    // random distributions
    normal_      = std::normal_distribution<double>(0, 1);
    uniform_     = std::uniform_real_distribution<double>(0., 1.);
    exponential_ = std::exponential_distribution<double>(proba_retraction_);

    update_kernel_variables();

    init_filopodia();
}


GrowthCone::GrowthCone(const GrowthCone &copy)
    : TopologicalNode(copy)
    , neuron_id_(copy.neuron_id_)
    , neurite_name_(copy.neurite_name_)
    , active_(true)
    , model_(copy.model_)
    , observables_(copy.observables_)
    , total_proba_(0.)
    , delta_angle_(0)
    , current_area_(copy.current_area_)
    , stuck_(false)
    , stopped_(false)
    , just_retracted_(false)
    , turning_(0)
    , turned_(0.)
    , max_stop_(copy.max_stop_)
    , current_stop_(0)
    , interacting_(false)
    , update_filopodia_(false)
    , aff_(copy.aff_)
    , filopodia_(copy.filopodia_)
    , move_(copy.move_)
    , sensing_angle_(copy.sensing_angle_)
    , scale_up_move_(copy.scale_up_move_)
    , duration_retraction_(copy.duration_retraction_)
    , proba_retraction_(copy.proba_retraction_)
    , retracting_todo_(-1)
    , retraction_time_(-1.)
    , speed_ratio_retraction_(copy.speed_ratio_retraction_)
    , proba_down_move_(copy.proba_down_move_)
    , max_sensing_angle_(copy.max_sensing_angle_)
    , min_filopodia_(copy.min_filopodia_)
    , num_filopodia_(copy.num_filopodia_)
    , old_angle_(copy.old_angle_)
{
    normal_  = std::normal_distribution<double>(0, 1);
    uniform_ = std::uniform_real_distribution<double>(0., 1.);
    update_kernel_variables();
}


GrowthCone::~GrowthCone() { assert(branch_.use_count() == 1); }


/**
 * Update the geometrical and topological informations of the
 * GrowthCone.
 * This is called to update a GrowthCone after it has been cloned from
 * when branching occurs.
 */
void GrowthCone::update_topology(BaseWeakNodePtr parent, NeuritePtr own_neurite,
                                 double distance_to_parent,
                                 const BPoint &position, double angle)
{
    parent_            = parent;
    centrifugal_order_ = parent.lock()->get_centrifugal_order() + 1;

    position_ = position;

    dist_to_soma_ = parent.lock()->get_distance_to_soma() + distance_to_parent;
    dist_to_parent_ = distance_to_parent;

    has_child_   = false;
    branch_      = std::make_shared<Branch>(position, dist_to_soma_);
    own_neurite_ = own_neurite;

    neuron_id_    = own_neurite->get_parent_neuron().lock()->get_gid();
    neurite_name_ = own_neurite->get_name();

    move_.angle = angle;
    old_angle_  = angle;
}


// ###########################################################
//      Grow, Retraction and Pruning
// ###########################################################


/**
 * @brief let's compute next step
 * Core function for the simulation, it shouldn't be overwritten
 *
 * Explanation:
 * 1. Compute the module of the next step
 *    It will be done from the extension model, critical_resourcee or
 *    random_walk
 *
 * 2. If module is positive, compute the angle and elongate:
 *    To grow it has to sense the environment and then convolve with the
 *    particular intrinsic distribution (e.g. random walk model)
 *    as set in the grwth cone model.
 *    'sense_environment' and 'accesible_environment'
 *    are defined in the abstract class GrowthCone (here),
 *    but it is required they are called inside the function
 *    'compute_directions' overriden by the GrowthCone model.
 *    This requirement is due to:
 *    - the functions need the sigma (persistence_length) defined in the
 *      models
 *    - some models required to perform a rotation before sense the
 *      environment
 *
 * 3. If module is negative retract:
 *    The retraction of a 0 length branch will imply the pruning of the
 *    GrowthCone
 *
 * @param rnd_engine
 */
void GrowthCone::grow(mtPtr rnd_engine, stype cone_n, double substep)
{
    // useful values
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    double current_time(0.), local_substep(substep), tmp, old_substep;
    double angle_widening, unstuck_angle, waiting_time;
    stype n_direction, i;

    std::vector<double> directions_weights;
    std::vector<std::string> new_pos_area;
    std::vector<bool> wall_presence;

    // remember distance done
    double final_distance_done = 0.;

    // we make local substeps until current time is equal to substep
    unsigned int loop_index = 0;

    while (current_time < substep)
    {
        loop_index++;
        if (loop_index > 1000)
        {
            printf("long loop for gc of neuron %lu on %i - substep %f "
                   "- current time %f - next substep %f\n",
                   own_neurite_->get_parent_neuron().lock()->get_gid(),
                   omp_id, substep, current_time, local_substep);
            throw std::runtime_error("Neuron stuck in an infinite loop");
        }

        local_substep = substep - current_time;

        // reset stopped status
        stopped_ = false;

        // initialize direction test variables
        directions_weights = std::vector<double>(filopodia_.size, 1.);
        wall_presence      = std::vector<bool>(filopodia_.size, false);

        new_pos_area = std::vector<std::string>(filopodia_.size,
                                                current_area_);

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
                // signal that we just finished retracting
                just_retracted_ = true;
            }
            else
            {
                retracting_todo_ -= local_substep;
            }

            // compute speed and module
            compute_speed(rnd_engine, local_substep);
            compute_module(local_substep);

            // retraction => module has to be negative
            move_.module = -std::abs(move_.module);

            // retract (expects positive distance)
            retraction(-move_.module, cone_n, omp_id);
        }
        else
        {
            // ============================= //
            // Extending, compute local pull //
            // ============================= //

            // check environment and mechanical interactions
            try
            {
                interacting_ = sense_surroundings(
                    directions_weights, wall_presence, substep, rnd_engine);
            }
            catch (const std::exception &except)
            {
                std::throw_with_nested(std::runtime_error(
                    "Passed from `GrowthCone::sense_surroundings`."));
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
                // ======== //
                // Forward! //
                // ======== //

                // check accessibility
                kernel().space_manager.check_accessibility(
                    directions_weights, filopodia_, get_position(), move_,
                    get_branch()->get_last_segment());

                // check for stuck/total_proba_
                // take optional GC rigidity into account
                compute_direction_probabilities(directions_weights,
                                                local_substep);

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
                    catch (...)
                    {
                        std::throw_with_nested(std::runtime_error(
                            "Passed from `GrowthCone::grow`."));
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

                // we elongated so we did not just get out of a retraction
                just_retracted_ = false;
            }
            else if (move_.module < 0)
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

        if (stuck_ or stopped_)
        {
            // ============================== //
            // We're stuck. Widen or retract? //
            // ============================== //

            // compute the angle widening necessary to unstuck
            unstuck_angle = kernel().space_manager.unstuck_angle(
                position_, move_.angle, filopodia_.finger_length, current_area_,
                omp_id);

            // compute the time necessary to reach that angle
            local_substep = std::min(
                std::max((unstuck_angle / 4. - move_.sigma_angle) / ONE_DEGREE,
                         1.),
                substep - current_time);

            // check whether we would retract in that time
            local_substep = check_retraction(local_substep, rnd_engine);

            // we do not retract, widen angle to unstuck
            change_sensing_angle(ONE_DEGREE * local_substep);
        }
        else
        {
            // reset retraction_time_
            retraction_time_ = -1.;
            final_distance_done += move_.module;
        }

        if (total_proba_ < 1. and not stuck_)
        {
            // widen angle depending on total_proba_
            angle_widening = 1. - total_proba_;
            change_sensing_angle(ONE_DEGREE * angle_widening * local_substep);
        }
        else if (move_.sigma_angle != sensing_angle_)
        {
            // bring move_.sigma_angle towards its default value based on
            // what was done during previous step
            change_sensing_angle(sgn(sensing_angle_ - move_.sigma_angle) *
                                 local_substep * ONE_DEGREE);
        }

        // update current time
        current_time += local_substep;

        // update retraction time if necessary
        if (retraction_time_ > 0.)
        {
            retraction_time_ -= local_substep;
        }

        // clear neighbors
        current_neighbors_.clear();

        // retraction step stop once the whole retraction is done (growth cone
        // pauses until the step is finished)
        if (just_retracted_)
        {
            break;
        }
    }

    // set the value of move_.module to the distance the GC did
    move_.module = final_distance_done;
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


void GrowthCone::retraction(double distance, stype cone_n, int omp_id)
{
    assert(distance >= 0.);

    // remove the points
    double distance_done;

    while (distance > 0 and branch_->size() > 1)
    {
        distance_done    = branch_->get_last_segment_length();
        BPolygonPtr poly = branch_->get_last_segment();
        BBox box;
        ObjectInfo info;

        if (poly != nullptr)
        {
            box = bg::return_envelope<BBox>(*(poly.get()));
            // there is one less segment than point, so size - 2 for segment
            info = std::make_tuple(neuron_id_, neurite_name_, get_node_id(),
                                   branch_->size() - 2);
        }

        double remaining = distance_done - distance;

        if (remaining < 1e-6) // distance is greater than what we just did
        {
            distance -= distance_done;
            branch_->retract();

            // we also remove the tree box
            if (poly != nullptr)
            {
                kernel().space_manager.remove_object(box, info, omp_id);
            }
        }
        else
        {
            BPoint p1 = branch_->xy_at(branch_->size() - 2);
            BPoint p2 = branch_->get_last_xy();

            double new_x =
                (p2.x() * remaining + p1.x() * (distance_done - remaining)) /
                distance_done;
            double new_y =
                (p2.y() * remaining + p1.y() * (distance_done - remaining)) /
                distance_done;

            BPoint new_p = BPoint(new_x, new_y);

            branch_->retract();

            // we remove the previous object and add the new, shorter one
            if (poly != nullptr)
            {
                try
                {
                    kernel().space_manager.remove_object(box, info, omp_id);
                }
                catch (...)
                {
                    std::throw_with_nested(std::runtime_error(
                        "Passed from `GrowthCone::retract`, coming from "
                        "space_manager::remove_object."));
                }

                try
                {
                    kernel().space_manager.add_object(
                        p1, new_p, get_diameter(), remaining,
                        own_neurite_->get_taper_rate(), info, branch_, omp_id);
                }
                catch (...)
                {
                    std::throw_with_nested(std::runtime_error(
                        "Passed from `GrowthCone::retract` coming from "
                        "space_manager::add_object (after retraction)."));
                }
            }

            distance = 0.;
        }
    }

    // set the new growth cone angle
    stype last = branch_->size();
    double x0, y0, x1, y1;
    BPoint p;

    if (last > 0)
    {
        p  = branch_->xy_at(last - 1);
        x1 = p.x();
        y1 = p.y();

        if (last > 1)
        {
            p  = branch_->xy_at(last - 2);
            x0 = p.x();
            y0 = p.y();
        }
        else
        {
            p  = TopologicalNode::get_position();
            x0 = p.x();
            y0 = p.y();
        }

        move_.angle = atan2(y1 - y0, x1 - x0);
        old_angle_  = move_.angle;
    }

    set_position(branch_->get_last_xy());

    // prune growth cone if necessary
    if (branch_->size() == 1)
    {
        prune(cone_n);
    }

    // check if we changed area
    if (using_environment_)
    {
        current_area_ = kernel().space_manager.get_containing_area(position_);
        update_growth_properties(current_area_);
        // reset move_.sigma_angle to its default value
        AreaPtr area = kernel().space_manager.get_area(current_area_);

        double old_sigma = move_.sigma_angle;
        move_.sigma_angle =
            sensing_angle_ * area->get_property(names::sensing_angle);

        if (old_sigma != move_.sigma_angle)
        {
            update_filopodia();
        }
    }

    // cannot be stuck_ or on low proba mode when retracting, so reset all
    stuck_       = false;
    total_proba_ = 1.;
    // also reset turning
    turning_ = 0;
    turned_  = 0.;
}


void GrowthCone::prune(stype cone_n)
{
#ifndef NDEBUG
    printf("pruning %lu %s %lu because (%i, %i, %f)\n", neuron_id_,
           neurite_name_.c_str(), get_node_id(), stuck_, stopped_,
           total_proba_);
    printf("prop: retraction %f, %i\n", retraction_time_, just_retracted());
#endif
    own_neurite_->delete_cone(cone_n);
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
bool GrowthCone::sense_surroundings(std::vector<double> &directions_weights,
                                    std::vector<bool> &wall_presence,
                                    double substep, mtPtr rnd_engine)
{
    if (sensing_required_)
    {
        double up_move = scale_up_move_ == 0 ? std::nan("") : scale_up_move_;

        return kernel().space_manager.sense(
            directions_weights, wall_presence, filopodia_, position_, move_,
            current_area_, proba_down_move_, up_move, aff_, substep,
            0.5 * get_diameter(), shared_from_this(), current_neighbors_);
    }

    return false;
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

    // check whether we're stopped
    stopped_ = (uniform_(*(rnd_engine.get())) > total_proba_);

    // we test whether we should make the next move or stop moving if moving
    // anywhere is too unlikely
    if (not stopped_)
    {
        double new_angle = move_.angle;
        stype default_direction;

        // select direction updates the angle by delta_angle without
        // renormalizing it
        select_direction(directions_weights, rnd_engine, substep, new_angle,
                         default_direction);

        BPolygonPtr last_segment = get_branch()->get_last_segment();

        // default angle
        double default_angle =
            filopodia_.directions[default_direction] + move_.angle;

        BPoint p(position_.x() + cos(new_angle) * move_.module,
                 position_.y() + sin(new_angle) * move_.module);

        BLineString line =
            kernel().space_manager.line_from_points(position_, p);

        // check that delta_angle is less than PI/2
        if (std::abs(new_angle - old_angle_) > 0.5 * M_PI)
        {
            stopped_ = true;
        }
        else
        {
            bool update         = false;
            double min_distance = std::numeric_limits<double>::max();
            double radius       = 0.5 * get_diameter();
            double distance;

            // check forbidden self overlap

            if (std::isnan(aff_.affinity_self) and last_segment != nullptr)
            {
                while (bg::covered_by(p, *(last_segment.get())) and
                       new_angle != default_angle)
                {
                    new_angle = 0.5 * (new_angle + default_angle);

                    p = BPoint(position_.x() + cos(new_angle) * move_.module,
                               position_.y() + sin(new_angle) * move_.module);
                }
            }

            // set the new angle (chose a valid position if target position
            // is outside the environment)
            if (using_environment_)
            {
                if (kernel().space_manager.intersects("environment", line))
                {
                    // stop at `radius` from the environment
                    bool intsct = kernel().space_manager.get_point_at_distance(
                        line, "environment", radius, p, distance);

                    if (intsct)
                    {
                        // check if we moved a bit or if we were stopped
                        if (distance > 0 and distance < min_distance)
                        {
                            // we moved: modify move_.module and substep
                            // accordingly
                            substep *= distance / move_.module;
                            move_.module = distance;
                            min_distance = distance;
                            update       = true;
                        }
                        else
                        {
                            stopped_ = true;
                        }
                    }
                }
            }

            // check forbidden overlap with other neurites
            if (not stopped_)
            {
                for (auto obstacle : current_neighbors_)
                {
                    // Note: for unknown reasons, this was previously
                    // "intersects", switched to covered_by but does not prevent
                    // cases where "intersection" in "get_point_at_distance" is
                    // empty.
                    if (bg::covered_by(p, *(obstacle.second.get())))
                    {
                        // stop at "radius" from the obstacle
                        bool intsct =
                            kernel().space_manager.get_point_at_distance(
                                line, obstacle.second, radius, p, distance);

                        // check if the intersection was indeed obtained
                        if (intsct)
                        {
                            // check if we moved a bit or if we were stopped
                            if (distance > 0 and distance < min_distance)
                            {
                                // we moved: modify move_.module and substep
                                // accordingly
                                substep *= distance / move_.module;
                                move_.module = distance;
                                min_distance = distance;
                                update       = true;
                            }
                            else
                            {
                                stopped_ = true;
                            }
                        }
                    }
                }

                if (update)
                {
                    if (move_.module < 1e-2)
                    {
                        stopped_ = true;
                    }
                    else
                    {
                        p = BPoint(
                            position_.x() + cos(new_angle) * move_.module,
                            position_.y() + sin(new_angle) * move_.module);
                    }
                }

                if (not stopped_)
                {
                    // update angle
                    delta_angle_ = new_angle - move_.angle;
                    move_.angle  = new_angle;

                    // send the new segment to the space manager
                    // note the size - 1 in the tuple because there is always
                    // one less segment than the number of points
                    try
                    {
                        kernel().space_manager.add_object(
                            position_, p, get_diameter(), move_.module,
                            own_neurite_->get_taper_rate(),
                            std::make_tuple(neuron_id_, neurite_name_,
                                            get_node_id(), branch_->size() - 1),
                            branch_, omp_id);
                    }
                    catch (...)
                    {
                        printf("module %f - delta angle: %f - old angle %f - "
                               "new angle %f on OMP %i\n",
                               move_.module, delta_angle_, old_angle_,
                               new_angle, omp_id);

                        std::cout << "p " << bg::wkt(p) << std::endl;

                        BPoint p2 = branch_->get_last_xy();
                        std::cout << "p2 " << bg::wkt(p2) << std::endl;

                        if (branch_->size() > 1)
                        {
                            BPoint p1 = branch_->xy_at(branch_->size() - 2);
                            std::cout << "p1 " << bg::wkt(p1) << std::endl;

                            double old_angle =
                                std::atan2(p2.y() - p1.y(), p2.x() - p1.x());
                            printf("old angle 2: %f on OMP %i\n", old_angle, omp_id);
                        }

                        std::throw_with_nested(std::runtime_error(
                            "Passed from `GrowthCone::make_move` on "
                            "OMP " + std::to_string(omp_id) + "."));
                    }

                    // store new position and angle
                    set_position(p);
                    old_angle_ = move_.angle;

                    // check if we switched to a new area
                    std::string new_area =
                        kernel().space_manager.get_containing_area(p);

                    if (new_area != current_area_)
                    {
                        update_growth_properties(new_area);
                        update_filopodia();
                    }
                }
            }
        }
    }
}


/**
 * Set GrowthCone position and update distances to parent and soma
 */
void GrowthCone::set_position(const BPoint &pos)
{
    position_ = pos;

    dist_to_parent_ = branch_->get_length();
    dist_to_soma_   = branch_->final_distance_to_soma();
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

    // move sensing angle must have been updated by previous call to
    // update_growth_properties
    num_filopodia_ = min_filopodia_;

    filopodia_.directions     = std::vector<double>(num_filopodia_);
    filopodia_.normal_weights = std::vector<double>(num_filopodia_);

    filopodia_.size = num_filopodia_;

    // update the filopodia normal weights
    update_filopodia();
}


// ###########################################################
//              Interface functions
// ###########################################################

double GrowthCone::get_growth_cone_speed() const { return move_.speed; }


double GrowthCone::get_module() const { return move_.module; }


void GrowthCone::set_angle(double angle) { move_.angle = angle; }


// ###########################################################
//              Set Get status
// ###########################################################

void GrowthCone::set_status(const statusMap &status)
{

    get_param(status, names::filopodia_wall_affinity, filopodia_.wall_affinity);
    get_param(status, names::scale_up_move, scale_up_move_);
    get_param(status, names::filopodia_min_number, min_filopodia_);

    double finger_length(filopodia_.finger_length);
    get_param(status, names::filopodia_finger_length, finger_length);

    if (finger_length < MIN_FILOPODIA_FINGER_LENGTH)
    {
        throw std::invalid_argument(
            names::filopodia_finger_length + " must " +
            "be greater or equal to " +
            std::to_string(MIN_FILOPODIA_FINGER_LENGTH) + ".");
    }

    filopodia_.finger_length = finger_length;

    if (min_filopodia_ < 10)
    {
        throw std::invalid_argument(
            "`filopodia_min_number` should be at least 10 to ensure physically "
            "relevant behaviors.");
    }

    get_param(status, names::proba_retraction, proba_retraction_);

    // update exponential distribution
    exponential_ = std::exponential_distribution<double>(proba_retraction_);

    get_param(status, names::duration_retraction, duration_retraction_);
    get_param(status, names::speed_ratio_retraction, speed_ratio_retraction_);
    get_param(status, names::proba_down_move, proba_down_move_);

    // sensing angle and persistence length
    double sa_tmp;
    bool b_sa = get_param(status, names::sensing_angle, sa_tmp);

    if (b_sa)
    {
        if (sa_tmp < 0)
        {
            throw std::invalid_argument("`sensing_angle` must be positive.");
        }

        sensing_angle_ = sa_tmp;
    }

    double msa = max_sensing_angle_;
    get_param(status, names::max_sensing_angle, msa);

    if (msa < sensing_angle_)
    {
        throw std::invalid_argument(
            "`sensing_angle` is greater than `max_sensing_angle`, please "
            "either increase the latter, or reduce `sensing_angle`.");
    }

    if (msa > M_PI)
    {
        throw std::invalid_argument(
            "`max_sensing_angle` must be less than 180°.");
    }

    max_sensing_angle_ = msa;

    // set affinities
    if (own_neurite_->get_type() == "axon")
    {
        get_param(status, names::affinity_axon_self, aff_.affinity_self);
        get_param(status, names::affinity_axon_axon_other_neuron,
                  aff_.affinity_axon_other_neuron);
        get_param(status, names::affinity_axon_dendrite_same_neuron,
                  aff_.affinity_dendrite_same_neuron);
        get_param(status, names::affinity_axon_dendrite_other_neuron,
                  aff_.affinity_dendrite_other_neuron);
        get_param(status, names::affinity_axon_soma_same_neuron,
                  aff_.affinity_soma_same_neuron);
        get_param(status, names::affinity_axon_soma_other_neuron,
                  aff_.affinity_soma_other_neuron);
    }
    else
    {
        get_param(status, names::affinity_dendrite_self, aff_.affinity_self);
        get_param(status, names::affinity_dendrite_axon_same_neuron,
                  aff_.affinity_axon_same_neuron);
        get_param(status, names::affinity_dendrite_axon_other_neuron,
                  aff_.affinity_axon_other_neuron);
        get_param(status, names::affinity_dendrite_dendrite_same_neuron,
                  aff_.affinity_dendrite_same_neuron);
        get_param(status, names::affinity_dendrite_dendrite_other_neuron,
                  aff_.affinity_dendrite_other_neuron);
        get_param(status, names::affinity_dendrite_soma_same_neuron,
                  aff_.affinity_soma_same_neuron);
        get_param(status, names::affinity_dendrite_soma_other_neuron,
                  aff_.affinity_soma_other_neuron);
    }

    // reset filopodia parameters
    init_filopodia();
}


void GrowthCone::get_status(statusMap &status) const
{
    set_param(status, names::filopodia_wall_affinity, filopodia_.wall_affinity,
              "");
    set_param(status, names::filopodia_finger_length, filopodia_.finger_length,
              "micrometer");
    set_param(status, names::filopodia_min_number, min_filopodia_, "");

    // @todo change behavior of sensing_angle and max_sensing angle:
    // sensing angle set the min/max of the filopodia angle and is limited by
    // max_sensing_angle
    set_param(status, names::sensing_angle, sensing_angle_, "rad");
    set_param(status, names::max_sensing_angle, max_sensing_angle_, "rad");

    set_param(status, names::scale_up_move, scale_up_move_, "micrometer");
    set_param(status, names::proba_down_move, proba_down_move_, "");

    set_param(status, names::proba_retraction, proba_retraction_, "");
    set_param(status, names::duration_retraction, duration_retraction_,
              "minute");

    // set affinities
    if (own_neurite_ != nullptr)
    {
        if (own_neurite_->get_type() == "axon")
        {
            set_param(status, names::affinity_axon_self, aff_.affinity_self,
                      "");
            set_param(status, names::affinity_axon_axon_other_neuron,
                      aff_.affinity_axon_other_neuron, "");
            set_param(status, names::affinity_axon_dendrite_same_neuron,
                      aff_.affinity_dendrite_same_neuron, "");
            set_param(status, names::affinity_axon_dendrite_other_neuron,
                      aff_.affinity_dendrite_other_neuron, "");
            set_param(status, names::affinity_axon_soma_same_neuron,
                      aff_.affinity_soma_same_neuron, "");
            set_param(status, names::affinity_axon_soma_other_neuron,
                      aff_.affinity_soma_other_neuron, "");
        }
        else
        {
            set_param(status, names::affinity_dendrite_self, aff_.affinity_self,
                      "");
            set_param(status, names::affinity_dendrite_axon_same_neuron,
                      aff_.affinity_axon_same_neuron, "");
            set_param(status, names::affinity_dendrite_axon_other_neuron,
                      aff_.affinity_axon_other_neuron, "");
            set_param(status, names::affinity_dendrite_dendrite_same_neuron,
                      aff_.affinity_dendrite_same_neuron, "");
            set_param(status, names::affinity_dendrite_dendrite_other_neuron,
                      aff_.affinity_dendrite_other_neuron, "");
            set_param(status, names::affinity_dendrite_soma_same_neuron,
                      aff_.affinity_soma_same_neuron, "");
            set_param(status, names::affinity_dendrite_soma_other_neuron,
                      aff_.affinity_soma_other_neuron, "");
        }
    }
}


/**
 * @brief Get the current value of one of the observables
 */
double GrowthCone::get_state(const std::string &observable) const
{
    double value = std::nan("");

    if (observable == "length")
    {
        value = branch_->get_length();
    }
    else if (observable == "speed")
    {
        value = move_.speed;
    }
    else if (observable == "angle")
    {
        value = move_.angle;
    }
    else if (observable == "status")
    {
        value = 2 * stuck_ + stopped_; // 0: moving, 1: stopped, 2: stuck
    }
    else if (observable == "retraction_time")
    {
        value = retraction_time_;
    }
    else if (observable == "diameter")
    {
        value = get_diameter();
    }
    else if (observable == "stepping_probability")
    {
        value = total_proba_;
    }

    return value;
}


/**
 * @brief Get the current value of one of the observables
 */
double GrowthCone::get_state(const std::string &observable,
                             std::string &unit) const
{
    double value = std::nan("");

    if (observable == "length")
    {
        value = branch_->get_length();
        unit  = "micrometer";
    }
    else if (observable == "speed")
    {
        value = move_.speed;
        unit  = "um/minute";
    }
    else if (observable == "angle")
    {
        value = move_.angle;
        unit  = "radian";
    }
    else if (observable == "status")
    {
        value = 2 * stuck_ + stopped_; // 0: moving, 1: stopped, 2: stuck
        unit  = "";
    }
    else if (observable == "retraction_time")
    {
        value = retraction_time_;
        unit  = "minute";
    }
    else if (observable == "diameter")
    {
        value = get_diameter();
        unit  = "micrometer";
    }
    else if (observable == "stepping_probability")
    {
        value = total_proba_;
        unit  = "";
    }

    return value;
}


double GrowthCone::get_self_affinity() const { return aff_.affinity_self; }


bool GrowthCone::just_retracted() const { return just_retracted_; }


bool GrowthCone::is_active() const { return active_; }


void GrowthCone::update_kernel_variables()
{
    using_environment_ = kernel().using_environment();
    sensing_required_ =
        using_environment_ or kernel().space_manager.interactions_on();

    // check change in resolution
    double old_resol = resol_;
    resol_           = kernel().simulation_manager.get_resolution();

    // check adaptive timestep
    adaptive_timestep_ = kernel().get_adaptive_timestep();
    timestep_divider_  = 1. / adaptive_timestep_;
}


void GrowthCone::change_sensing_angle(double angle)
{
    double old_sigma = move_.sigma_angle;
    // set the local modifier (1. if not using env)
    double area_modifier = 1.;

    if (using_environment_)
    {
        AreaPtr area  = kernel().space_manager.get_area(current_area_);
        area_modifier = area->get_property(names::sensing_angle);
    }

    if (angle > 0.)
    {
        move_.sigma_angle =
            std::min(move_.sigma_angle + angle, max_sensing_angle_);

        double dangle = std::abs(move_.angle - old_angle_);

        // when not moving, start turning
        if (turning_ != 0 and dangle < 0.39269908169872414)
        {
            move_.angle += turning_ * angle;
            turned_ += turning_ * angle;
        }
    }
    else
    {
        move_.sigma_angle = std::min(
            std::max(move_.sigma_angle + angle, sensing_angle_ * area_modifier),
            max_sensing_angle_);
    }

    // change filopodia normal weights if necessary
    if (old_sigma != move_.sigma_angle)
    {
        update_filopodia_ = true;
    }
}


void GrowthCone::update_filopodia()
{
    double dtheta = move_.sigma_angle / (filopodia_.size - 1);
    double angle;

    // set the angles
    for (int n_angle = 0; n_angle < num_filopodia_; n_angle++)
    {
        angle = dtheta * n_angle - 0.5 * move_.sigma_angle;
        filopodia_.directions[n_angle] = angle;
    }
}


stype GrowthCone::get_neuron_id() const { return neuron_id_; }


const std::string &GrowthCone::get_neurite_name() const
{
    return neurite_name_;
}


const std::string &GrowthCone::get_model_name() const { return model_; }


const BPolygonPtr GrowthCone::get_last_segment() const
{
    return branch_->get_last_segment();
}

} // namespace growth
