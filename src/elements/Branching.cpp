/*
 * Branching.cpp
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

#include "Branching.hpp"

// C++ includes
#include <cmath>
#include <functional>
#include <memory>

// kernel includes
#include "kernel_manager.hpp"

// elements include
#include "GrowthCone.hpp"
#include "Neurite.hpp"
#include "Node.hpp"

// spatial include
#include "search.hpp"


namespace growth
{

const Event invalid_ev(std::make_tuple(Time(), 0, std::string(""), -1, -1));


Branching::Branching()
    : Branching(nullptr)
{
}

Branching::Branching(NeuritePtr neurite)
    : neurite_(neurite)
    // parameters for van Pelt model
    , use_van_pelt_(false)
    , next_vanpelt_event_(invalid_ev)
    , B_(VP_B)
    , E_(VP_S)
    , S_(VP_E)
    , T_(VP_T)
    // parameters for uniform branching
    , use_uniform_branching_(false)
    , uniform_branching_rate_(UNIFORM_BRANCHING_RATE)
    , next_uniform_event_(invalid_ev)
    , latbranch_dist_(20.) // todo: make user-defined param
    // parameters for front lateral branching
    , use_flpl_branching_(false)
    , flpl_branching_rate_(UNIFORM_BRANCHING_RATE)
    , next_flpl_event_(invalid_ev)
    // parameters for uniform split
    , use_uniform_split_(false)
    , uniform_split_rate_(UNIFORM_BRANCHING_RATE)
    , next_usplit_event_(invalid_ev)
{
    exponential_uniform_ =
        std::exponential_distribution<double>(uniform_branching_rate_);
    exponential_flpl_ =
        std::exponential_distribution<double>(flpl_branching_rate_);
    exponential_usplit_ =
        std::exponential_distribution<double>(uniform_split_rate_);
}


/**
 * @brief Initialize the scheduler mechanism for event generation
 * @details The event generation needs to be initialized for each branching
 * model. The initialization take care of reconnecting the simulation with the
 * previous one
 * even if the resolution of the simulation has changed
 *
 * @param rnd_engine
 *
 * @param new_resolution_ratio ratio between previous resolution and new
 * resolution
 *
 * @param last_simulation_n_steps steps of simulation in the previous session
 */
void Branching::initialize_next_event(mtPtr rnd_engine)
{
    if (use_uniform_branching_ and next_uniform_event_ == invalid_ev)
    {
        compute_uniform_event(rnd_engine);
    }
    if (use_flpl_branching_ and next_flpl_event_ == invalid_ev)
    {
        compute_flpl_event(rnd_engine);
    }
    if (use_uniform_split_ and next_usplit_event_ == invalid_ev)
    {
        compute_usplit_event(rnd_engine);
    }
    if (use_van_pelt_ and next_vanpelt_event_ == invalid_ev)
    {
        compute_vanpelt_event(rnd_engine);
    }
}


/**
 * Set the event content based on the duration between the current time and
 * the next event
 */
void Branching::set_branching_event(Event &ev, signed char ev_type,
                                    double duration)
{
    Time ev_time = kernel().simulation_manager.get_time();

    // separate duration into days, hours, minutes, seconds
    double total_hours = std::floor(duration / 60.);

    double seconds = (duration - std::floor(duration)) * 60;
    ev_time.set_sec(ev_time.get_sec() + seconds);

    unsigned char minutes = duration - 60 * total_hours;
    ev_time.set_min(ev_time.get_min() + minutes);

    stype days = total_hours / 24.;
    ev_time.set_day(ev_time.get_day() + days);

    unsigned char hours = total_hours - 24 * days;
    ev_time.set_hour(ev_time.get_hour() + hours);

    // set the informations of the event
    auto n                   = neurite_->get_parent_neuron().lock();
    stype neuron_gid         = n->get_gid();
    std::string neurite_name = neurite_->get_name();

    ev = std::make_tuple(ev_time, neuron_gid, neurite_name, -1, ev_type);
}


/**
 * @brief Verify and actuate if a branching event is scheduled for the present
 * step
 *
 * @param step present step of the simulator
 * @param rnd_engine
 *
 * @return whether the branching event was successful or not.
 */
bool Branching::branching_event(mtPtr rnd_engine, const Event &ev)
{
    // uniform_branching_event
    bool uniform_occurence =
        std::get<edata::TIME>(next_uniform_event_) == std::get<edata::TIME>(ev);
    bool flpl_occurence =
        std::get<edata::TIME>(next_flpl_event_) == std::get<edata::TIME>(ev);
    bool usplit_occurence =
        std::get<edata::TIME>(next_usplit_event_) == std::get<edata::TIME>(ev);
    bool van_pelt_occurence =
        std::get<edata::TIME>(next_vanpelt_event_) == std::get<edata::TIME>(ev);
    bool res_occurence = std::get<edata::GC>(ev) != -1;

    // prepare pointers to nodes and branching_point
    TNodePtr branching_node(nullptr);
    GCPtr second_cone(nullptr);
    NodePtr new_node = nullptr;
    stype branching_point;

    bool success(false);

    // check resource-based split
    if (res_occurence)
    {
        success = res_new_branch(branching_node, new_node, branching_point,
                                 rnd_engine, second_cone, ev);
    }
    // check uniform event
    else if (use_uniform_branching_ and uniform_occurence)
    {
        success = uniform_new_branch(branching_node, new_node, branching_point,
                                     rnd_engine);
    }
    // check flpl event
    else if (use_flpl_branching_ and flpl_occurence)
    {
        success = flpl_new_branch(branching_node, new_node, branching_point,
                                  rnd_engine);
    }
    // check usplit event
    else if (use_uniform_split_ and usplit_occurence)
    {
        success = usplit_new_branch(branching_node, new_node, branching_point,
                                    rnd_engine, second_cone);
    }
    // verify vanpelt event
    else if (use_van_pelt_ and van_pelt_occurence)
    {
        success = vanpelt_new_branch(branching_node, new_node, branching_point,
                                     rnd_engine, second_cone);
    }

    if (success)
    {
        // update the R-tree: all points from branch_point on need to be
        // removed and reassigned to the new point.
        int omp_id = kernel().parallelism_manager.get_thread_local_id();

        kernel().space_manager.update_objects_branching(
            branching_node, new_node, branching_point,
            neurite_->get_parent_neuron().lock()->get_gid(),
            neurite_->get_name(), omp_id);

        // for split events
        if (van_pelt_occurence or res_occurence or usplit_occurence)
        {
            update_splitting_cones(branching_node, second_cone, new_node);
        }
    }

    // note that critical resource splitting is instantaneous and unplanned
    // thus it is never communicated through branching events
    // @todo change that

    return success;
}


void Branching::update_splitting_cones(TNodePtr branching_cone,
                                       GCPtr second_cone, NodePtr new_node)
{
    GCPtr old_cone   = std::dynamic_pointer_cast<GrowthCone>(branching_cone);
    BranchPtr branch = new_node->get_branch();

    // to prevent overlap between the two second_cones, we put them on the two
    // corners (last points) of the branch
    std::pair<BPoint, BPoint> lps = branch->get_last_points();
    double new_angle              = second_cone->get_state("angle");
    double module                 = branch->get_last_segment_length();

    BPoint lp1(lps.first), lp2(lps.second);
    BPoint pos = second_cone->get_position();
    BPoint tmp = BPoint(pos.x() + module * cos(new_angle),
                        pos.y() + module * sin(new_angle));

    double d1(bg::distance(tmp, lp1)), d2(bg::distance(tmp, lp2));
    BPoint pos1 = BPoint(0.5 * (pos.x() + lp1.x()), 0.5 * (pos.y() + lp1.y()));
    BPoint pos2 = BPoint(0.5 * (pos.x() + lp2.x()), 0.5 * (pos.y() + lp2.y()));

    if (d2 > d1)
    {
        std::swap(pos1, pos2);
    }

    // delete last segment from parent branch and make the two cones
    // start from the previous positions and end on lp1/lp2

    int omp_id = kernel().parallelism_manager.get_thread_local_id();

    // remove last point (previous position) from old branch
    BPolygonPtr last_seg = branch->get_last_segment();
    branch->retract();

    if (last_seg != nullptr)
    {
        BBox box;
        bg::envelope(*(last_seg.get()), box);

        // the last segment is initial size - 2 i.e. new size - 1
        ObjectInfo info = std::make_tuple(
            neurite_->get_parent_neuron().lock()->get_gid(),
            neurite_->get_name(), new_node->get_node_id(), branch->size() - 1);

        kernel().space_manager.remove_object(box, info, omp_id);
    }

    tmp          = branch->get_last_xy();
    double d_som = branch->final_distance_to_soma();

    branching_cone->set_first_point(tmp, d_som);
    second_cone->set_first_point(tmp, d_som);

#ifndef NDEBUG
    printf("retracted\n");
    std::cout << bg::wkt(*(new_node->get_branch()->get_last_segment().get()))
              << std::endl;
    std::cout << bg::wkt(pos2) << "\n"
              << bg::wkt(lp1) << "\n"
              << bg::wkt(lp2) << "\n"
              << bg::wkt(pos) << std::endl;
    // std::cout << bg::wkt(*(second_cone->get_branch()->get_last_segment())) <<
    // std::endl;
#endif

    // update the second cone and its branch if possible
    if (kernel().space_manager.env_contains(pos2))
    {
        try
        {
            second_cone->geometry_.position = pos2;

            double module = bg::distance(tmp, pos2);

            kernel().space_manager.add_object(
                tmp, pos2, second_cone->get_diameter(), module,
                neurite_->get_taper_rate(),
                std::make_tuple(neurite_->get_parent_neuron().lock()->get_gid(),
                                neurite_->get_name(),
                                second_cone->get_node_id(), 0),
                second_cone->get_branch(), omp_id);
        }
        catch (...)
        {
            std::throw_with_nested(std::runtime_error(
                "Passed from `Neurite::growth_cone_split`."));
        }
    }

    // move the branching cone if possible and update its branch
    if (kernel().space_manager.env_contains(pos1))
    {
        try
        {
            old_cone->geometry_.position = pos1;

            double module = bg::distance(tmp, pos1);

            kernel().space_manager.add_object(
                tmp, pos1, branching_cone->get_diameter(), module,
                neurite_->get_taper_rate(),
                std::make_tuple(neurite_->get_parent_neuron().lock()->get_gid(),
                                neurite_->get_name(),
                                branching_cone->get_node_id(), 0),
                branching_cone->get_branch(), omp_id);
        }
        catch (...)
        {
            std::throw_with_nested(std::runtime_error(
                "Passed from `Neurite::growth_cone_split`."));
        }
    }

#ifndef NDEBUG
    printf("after branching corrections\n");
    std::cout << bg::wkt(
                     *(branching_cone->get_branch()->get_last_segment().get()))
              << std::endl;
    std::cout << bg::wkt(old_cone->get_position()) << std::endl;
    std::cout << bg::wkt(*(second_cone->get_branch()->get_last_segment().get()))
              << std::endl;
    std::cout << bg::wkt(second_cone->get_position()) << std::endl;
#endif
}


//#######################################################
//              Uniform Branching Model
//#######################################################

/**
 * @brief Compute the next uniform branching event
 * @details
 * The segment where next branching event will happen is computed
 * with an uniform distribution with respect to the length of segments
 * considered. The time (in step unit) of next lateral branching event is
 * computed with an exponential distribution whose rate is set with
 * Branching::set_status(...) Values are stored in the Branching instance
 * itself.
 *
 * @param rnd_engine
 */
void Branching::compute_uniform_event(mtPtr rnd_engine)
{
    // here we compute next event with exponential distribution
    // we add one to save the case in which the distance between
    // two branching event is so short to appear zero.
    if (not neurite_->growth_cones_.empty())
    {
        double duration = exponential_uniform_(*(rnd_engine).get());

        set_branching_event(next_uniform_event_, names::lateral_branching,
                            duration);

        // send it to the simulation and recorder managers
        kernel().simulation_manager.new_branching_event(next_uniform_event_);
    }
}


/**
 * @brief Branch the neurite in a uniformly randoim chosen node
 *
 * This function is enable with the flag 'use_lateral_branching'
 * This function implement the branchin throught
 * @function Neurite::lateral_branching
 *
 * @param rnd_engine
 */
bool Branching::uniform_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                                   stype &branching_point, mtPtr rnd_engine)
{
#ifndef NDEBUG
    printf("@@@@@@@ Lateral branching @@@@@@@@\n");
#endif
    branching_node = nullptr;
    new_node       = nullptr;
    GCPtr branching_cone;

    // Compute the total length of dendritic branch for nodes and gcs
    //###################################################################
    if (not neurite_->growth_cones_.empty())
    {
        bool success        = false;
        double total_length = 0;

        for (auto &cone : neurite_->gc_range())
        {
            if (not cone.second->is_dead() and
                cone.second->get_branch_size() > 2 * latbranch_dist_)
            {
                total_length += cone.second->get_branch()->get_length();
            }
        }
        for (auto &node : neurite_->nodes_)
        {
            if (node.second->get_branch()->size() > 2 * latbranch_dist_)
            {
                total_length += node.second->get_branch()->get_length();
            }
        }
        //###################################################################


        // Pick up a random number and check which interval it belongs: length_i
        // < random < length_i+1 It is equivalent to make a weight choice over
        // the branches
        //###################################################################
        double random_length  = total_length * uniform_(*(rnd_engine).get());
        double current_length = 0.;

        for (auto &cone : neurite_->gc_range())
        {
            while (branching_node == nullptr)
            {
                if (not cone.second->is_dead() and
                    cone.second->get_branch_size() > 2 * latbranch_dist_)
                {
                    current_length += cone.second->get_branch()->get_length();
                    if (current_length < random_length)
                    {
                        branching_cone = cone.second;
                        branching_node = branching_cone;
                    }
                }
            }
        }
        for (auto &node : neurite_->nodes_)
        {
            while (branching_node == nullptr)
            {
                if (not node.second->is_dead() and
                    node.second->get_branch_size() > 2 * latbranch_dist_)
                {
                    current_length += node.second->get_branch()->get_length();
                    if (current_length < random_length)
                    {
                        branching_node = node.second;
                    }
                }
            }
        }

        //###################################################################

        // if no node was suited for lateral branching skip the branching.
        if (branching_node == nullptr)
        {
            // choose the point uniformly on the branch, except for first 2 and
            // last 2 points.
            branching_point = uniform_(*(rnd_engine).get()) *
                                  (branching_node->get_branch()->size() - 4) +
                              2;

            // actuate lateral branching on the elected node through the
            // NEURITE.
            success = neurite_->lateral_branching(
                branching_node, branching_point, new_node, rnd_engine);
            next_uniform_event_ = invalid_ev;
        }

        compute_uniform_event(rnd_engine);

        return success;
    }


    next_uniform_event_ = invalid_ev;

    return false;
}


/**
 * @brief Compute the next non-uniform flpl branching event
 * @details
 * The segment where next branching event will happen is computed
 * with an power law distribution with respect to the distance from the branch
 * tip.
 * The time (in step unit) of next flpl branching event is computed with an
 * exponential distribution
 * whose rate is set with Branching::set_status(...)
 * Values are stored in the Branching instance itself.
 *
 * @param rnd_engine
 */
void Branching::compute_flpl_event(mtPtr rnd_engine)
{
    if (not neurite_->growth_cones_.empty())
    {
        // here we compute next event with exponential distribution
        double duration = exponential_flpl_(*(rnd_engine).get());

        set_branching_event(next_flpl_event_, names::lateral_branching,
                            duration);

        // send it to the simulation and recorder managers
        kernel().simulation_manager.new_branching_event(next_flpl_event_);
    }
}


/**
 * @brief Branch on a random node following power law distribution from the tips
 *
 * This function is enabled through the flag 'use_flpl_branching'
 * This function implements the branchin throught
 * @function Neurite::lateral_branching
 *
 * @param rnd_engine
 */
bool Branching::flpl_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                                stype &branching_point, mtPtr rnd_engine)
{
#ifndef NDEBUG
    printf("@@@@@@@ Lateral branching (FLPL) @@@@@@@@\n");
#endif
    branching_node = nullptr;
    GCPtr branching_cone;

    // Compute the total length of dendritic branch for nodes and gcs
    //###################################################################
    if (not neurite_->growth_cones_.empty())
    {
        bool success        = false;
        double total_length = 0;

        for (auto &cone : neurite_->gc_range())
        {
            if (not cone.second->is_dead() and
                cone.second->get_branch_size() > 2 * latbranch_dist_)
            {
                total_length += cone.second->get_branch()->get_length();
            }
        }
        for (auto &node : neurite_->nodes_)
        {
            if (node.second->get_branch()->size() > 2 * latbranch_dist_)
            {
                total_length += node.second->get_branch()->get_length();
            }
        }
        //###################################################################


        // Pick up a random number and check which interval it belongs: length_i
        // < random < length_i+1 It is equivalent to make a weight choice over
        // the branches
        //###################################################################
        double random_length  = total_length * uniform_(*(rnd_engine).get());
        double current_length = 0.;

        for (auto &cone : neurite_->gc_range())
        {
            while (branching_node == nullptr)
            {
                if (not cone.second->is_dead() and
                    cone.second->get_branch_size() > 2 * latbranch_dist_)
                {
                    current_length += cone.second->get_branch()->get_length();
                    if (current_length < random_length)
                    {
                        branching_cone = cone.second;
                        branching_node = branching_cone;
                    }
                }
            }
        }
        for (auto &node : neurite_->nodes_)
        {
            while (branching_node == nullptr)
            {
                if (not node.second->is_dead() and
                    node.second->get_branch_size() > 2 * latbranch_dist_)
                {
                    current_length += node.second->get_branch()->get_length();
                    if (current_length < random_length)
                    {
                        branching_node = node.second;
                    }
                }
            }
        }
        //###################################################################

        // if no node was suited for lateral branching skip the branching.
        if (branching_node == nullptr)
        {
            // choose the point with a power law distribution over the branch
            // length,
            // where y is a uniform variate, n is the distribution power,
            // x0 and x1 define the range of the distribution, and x is your
            // power-law distributed variate.
            double y = uniform_(*(rnd_engine).get());
            int x_0  = latbranch_dist_;
            int x_1 =
                branching_node->get_branch()->get_length() - latbranch_dist_;
            // TODO check if we should not make n user-defined.
            int n    = 2;
            double a = (powf(x_1, (n + 1)) - powf(x_0, (n + 1))) * y +
                       powf(x_0, (n + 1));
            double branching_dist = powf(a, 1. / (n + 1));

            branching_point = get_closest_point(branching_node, branching_dist);

            // lateral branching on the elected node through the NEURITE.
            success = neurite_->lateral_branching(
                branching_node, branching_point, new_node, rnd_engine);

            next_flpl_event_ = invalid_ev;
        }

        compute_flpl_event(rnd_engine);

        return success;
    }

    next_flpl_event_ = invalid_ev;

    return false;
}


/**
 * @brief Compute the next uniform split event
 * @details
 * The time (in step unit) of next uniform split event is computed with an
 * exponential distribution
 * whose rate is set with Branching::set_status(...)
 * Values are stored in the Branching instance itself.
 *
 * @param rnd_engine
 */
bool Branching::usplit_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                                  stype &branching_point, mtPtr rnd_engine,
                                  GCPtr &second_cone)
{
    branching_node = nullptr;
    new_node       = nullptr;

    if (not neurite_->growth_cones_.empty())
    {
        // simple implementation of a reservoir sampling algorithm for weigthed
        // choice
        GCPtr next_usplit_cone;
        double total_weight = 0.;
        std::unordered_map<stype, double> weights;

        for (auto cone : neurite_->growth_cones_)
        {
            weights[cone.first] = 1.;
            total_weight += 1.;
        }
        double extracted = uniform_(*(rnd_engine).get()) * total_weight;

        total_weight = 0;
        for (auto cone : neurite_->growth_cones_)
        {
            total_weight += weights[cone.first];
            next_usplit_cone = cone.second;

            if (total_weight >= extracted)
            {
                break;
            }
        }

        // set branching node and point
        branching_node  = next_usplit_cone;
        branching_point = branching_node->get_branch()->size() - 1;

        //@TODO define new legth for van pelt
        double new_length = 0.;
        double new_angle, old_angle;
        double old_diameter = next_usplit_cone->get_diameter();
        double new_diameter = old_diameter;
        neurite_->gc_split_angles_diameter(rnd_engine, new_angle, old_angle,
                                           new_diameter, old_diameter);
        bool success = neurite_->growth_cone_split(
            next_usplit_cone, new_length, new_angle, old_angle, new_diameter,
            old_diameter, new_node, second_cone);

        next_usplit_event_ = invalid_ev;
        compute_usplit_event(rnd_engine);
        return success;
    }

    next_usplit_event_ = invalid_ev;
    return false;
}


/**
 * @brief Compute the next uniform split event
 * @details
 * The time (in step unit) of next uniform split event is computed with an
 * exponential distribution
 * whose rate is set with Branching::set_status(...)
 * Values are stored in the Branching instance itself.
 * @brief Branch the neurite in a uniformly randoim chosen node
 *
 * @param rnd_engine
 */
void Branching::compute_usplit_event(mtPtr rnd_engine)
{
    if (not neurite_->growth_cones_.empty())
    {
        // here we compute next event with exponential distribution
        double duration = exponential_usplit_(*(rnd_engine).get());

        set_branching_event(next_usplit_event_, names::gc_splitting, duration);

        // send it to the simulation and recorder managers
        kernel().simulation_manager.new_branching_event(next_usplit_event_);
    }
}


//###################################################
//                 Van Pelt Model
//###################################################

/**
 * @brief Compute the next vanpelt branching event
 * @details
 * The time of branching is computed with the statistical branching algorithm
 * from Van Pelt
 * (see article for details) the relevant parameters are B, E, T
 * The growth cone who is going to branch is computed with a reservoir
 * sample algorithm; the weight are defined in the Van Pelt model and
 * the model parameter S is set in the Branching::set_status
 *
 * This distribution is verified to reproduce the Van Pelt results.
 *
 * @param rnd_engine
 */
void Branching::compute_vanpelt_event(mtPtr rnd_engine)
{
    if (not neurite_->growth_cones_.empty())
    {
        // get the current second to compute the time-dependent
        // exponential decreasing probability of having a branch.
        double t_0    = kernel().simulation_manager.get_current_minutes();
        stype num_gcs = neurite_->growth_cones_.size() +
                        neurite_->growth_cones_inactive_.size();

        double delta = exp((t_0 + 1) / T_) * T_ / B_ * powf(num_gcs, E_);

        //-t_0 - log(exp(-T_ * t_0) -
        // powf(neurite_->growth_cones_.size(), E_ - 1) / B_) /T_;
        exponential_ = std::exponential_distribution<double>(1. / delta);

        double duration = exponential_(*(rnd_engine).get());

        if (duration > std::numeric_limits<stype>::max())
        {
            next_vanpelt_event_ = invalid_ev;
        }
        else
        {
            set_branching_event(next_vanpelt_event_, names::gc_splitting,
                                duration);

            // send it to the simulation and recorder managers
            kernel().simulation_manager.new_branching_event(
                next_vanpelt_event_);
        }
    }
}


/**
 * @brief Compute the Van Pelt normalization factor and update GrowthCones
 * This function will be run every time a branching happen if the
 * 'use_van_pelt' flag is set True.
 * This function implement the branchin throught
 * @function growth_cone_split(...)
 *
 */
bool Branching::vanpelt_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                                   stype &branching_point, mtPtr rnd_engine,
                                   GCPtr &second_cone)
{
    branching_node = nullptr;
    new_node       = nullptr;

    if (not neurite_->growth_cones_.empty())
    {
        // simple implementation of a reservoir sampling algorithm for weigthed
        // choice
        GCPtr nex_vanpelt_cone;
        double weight, total_weight(0);
        std::unordered_map<stype, double> weights;

        for (auto cone : neurite_->growth_cones_)
        {
            weight = powf(2, -cone.second->get_centrifugal_order() * S_);
            weights[cone.first] = weight;
            total_weight += weight;
        }
        double extracted = uniform_(*(rnd_engine).get()) * total_weight;

        total_weight = 0;
        for (auto cone : neurite_->growth_cones_)
        {
            total_weight += weights[cone.first];
            nex_vanpelt_cone = cone.second;
            if (total_weight >= extracted)
            {
                break;
            }
        }

        // set branching node and point
        branching_node  = nex_vanpelt_cone;
        branching_point = branching_node->get_branch()->size() - 1;

        //@TODO define new legth for van pelt
        double new_length = 0.;
        double new_angle, old_angle;
        double old_diameter = nex_vanpelt_cone->get_diameter();
        double new_diameter = old_diameter;
        neurite_->gc_split_angles_diameter(rnd_engine, new_angle, old_angle,
                                           new_diameter, old_diameter);
        bool success = neurite_->growth_cone_split(
            nex_vanpelt_cone, new_length, new_angle, old_angle, new_diameter,
            old_diameter, new_node, second_cone);

        next_vanpelt_event_ = invalid_ev;
        compute_vanpelt_event(rnd_engine);
#ifndef NDEBUG
        printf("VP branching on %lu, %s, %lu at %lu\n",
               neurite_->get_parent_neuron().lock()->get_gid(),
               neurite_->get_name().c_str(), branching_node->get_node_id(),
               branching_point);
#endif
        return success;
    }

    next_vanpelt_event_ = invalid_ev;

    return false;
}


//###################################################
//                  Tubuline model
//###################################################
/**
 * @brief Growth Cone split when critical_resource value exceed the treshold
 * @details
 * The GC is selected in Branching::branching_event with the chosen condition
 * This function implement the branchin throught
 * @function growth_cone_split(...)
 * The splitting angles are computed with a experimentally derived distribution
 * whose
 * parameters are set in Branching::set_status(...)
 *
 * @param rnd_engine
 * @param splitting_cone the cone who is going to split
 */
bool Branching::res_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                               stype &branching_point, mtPtr rnd_engine,
                               GCPtr &second_cone, const Event &ev)
{
    int gc_id = std::get<edata::GC>(ev);
    GCPtr gc  = neurite_->growth_cones_[gc_id];

    // set branching node and point
    branching_node  = gc;
    branching_point = gc->get_branch()->size() - 1;

    double new_length = 0.;
    double new_angle, old_angle;
    double old_diameter = gc->get_diameter();
    double new_diameter = old_diameter;

    neurite_->gc_split_angles_diameter(rnd_engine, new_angle, old_angle,
                                       new_diameter, old_diameter);

    bool success = neurite_->growth_cone_split(
        gc, new_length, new_angle, old_angle, new_diameter, old_diameter,
        new_node, second_cone);

    return success;
}


void Branching::set_status(const statusMap &status)
{

    //                 Van Pelt Params
    //###################################################
    get_param(status, names::use_van_pelt, use_van_pelt_);

    get_param(status, names::B, B_);
    get_param(status, names::E, E_);
    get_param(status, names::S, S_);
    get_param(status, names::T, T_);

    if (not use_van_pelt_)
    {
        next_vanpelt_event_ = invalid_ev;
    }

    //                 Uniform_branching Params
    //###################################################
    get_param(status, names::use_uniform_branching, use_uniform_branching_);

    get_param(status, names::uniform_branching_rate, uniform_branching_rate_);

    exponential_uniform_ =
        std::exponential_distribution<double>(uniform_branching_rate_);

    if (not use_uniform_branching_)
    {
        next_uniform_event_ = invalid_ev;
    }

    //                 FLPL branching Params
    //###################################################

    get_param(status, names::use_flpl_branching, use_flpl_branching_);

    get_param(status, names::flpl_branching_rate, flpl_branching_rate_);

    exponential_flpl_ =
        std::exponential_distribution<double>(flpl_branching_rate_);

    //                 Uniform split Params
    //###################################################

    get_param(status, names::use_uniform_split, use_uniform_split_);

    get_param(status, names::uniform_split_rate, uniform_split_rate_);

    exponential_usplit_ =
        std::exponential_distribution<double>(uniform_split_rate_);

    if (not use_uniform_split_)
    {
        next_usplit_event_ = invalid_ev;
    }
}


void Branching::get_status(statusMap &status) const
{
    set_param(status, names::use_van_pelt, use_van_pelt_, "");
    if (use_van_pelt_)
    {
        set_param(status, names::B, B_, "count / minute");
        set_param(status, names::E, E_, "");
        set_param(status, names::S, S_, "");
        set_param(status, names::T, T_, "minute");
    }

    set_param(status, names::use_uniform_branching, use_uniform_branching_, "");
    if (use_uniform_branching_)
    {
        set_param(status, names::uniform_branching_rate,
                  uniform_branching_rate_, "count / minute");
    }

    set_param(status, names::use_flpl_branching, use_flpl_branching_, "");
    if (use_flpl_branching_)
    {
        set_param(status, names::flpl_branching_rate, flpl_branching_rate_,
                  "count / minute");
    }

    set_param(status, names::use_uniform_split, use_uniform_split_, "");
    if (use_uniform_split_)
    {
        set_param(status, names::uniform_split_rate, uniform_split_rate_,
                  "count / minute");
    }
}

} // namespace growth
