/*
 * Neurite.cpp
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

#define _USE_MATH_DEFINES

#include "Neurite.hpp"

// c++ includes
#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>

// elements includes
#include "GrowthCone.hpp"
#include "Neuron.hpp"

// kernel includes
#include "growth_time.hpp"
#include "kernel_manager.hpp"
#include "search.hpp"

// models includes
#include "ActinWave.hpp"
#include "Node.hpp"
//~ #include "gc_critical.hpp"
#include "growth_names.hpp"

// debug
#include <cassert>
#include <typeinfo>


namespace growth
{

/**
 * Neurite constructor:
 * \param neurite_type th etype of the neurite (either "axon" or "dendrite").
 * \param p pointer to parent class (Neuron).
 *
 * The neurite has no points but branches, so the constructor will call a
 * growthConebuild which will instanciate an branch with a random angle.
 */
Neurite::Neurite(const std::string &name, const std::string &neurite_type,
                 const std::string &gc_model, NeuronWeakPtr p)
    : parent_(p)
    , branching_model_(std::make_shared<Branching>())
    , name_(name)
    , observables_(
          {"length", "speed", "num_growth_cones", "retraction_time", "status"})
    , num_created_nodes_(0)
    , growth_cone_model_(gc_model)
    , neurite_type_(neurite_type)
    , active_(true)
    , max_gc_num_(MAX_GC_NUM)
    , max_arbor_len_(MAX_ARBOR_LENGTH)
    , fixed_arbor_len_(0.)
    // parameters for van Pelt branching
    , lateral_branching_angle_mean_(LATERAL_BRANCHING_ANGLE_MEAN)
    , lateral_branching_angle_std_(LATERAL_BRANCHING_ANGLE_STD)
    , gc_split_angle_mean_(GC_SPLIT_ANGLE_MEAN)
    , gc_split_angle_std_(GC_SPLIT_ANGLE_STD)
    , diameter_eta_exp_(DIAMETER_ETA_EXP)
    , taper_rate_(THINNING_RATIO)
    , min_diameter_(MIN_DIAMETER)
    , diameter_ratio_std_(DIAM_RATIO_STD)
    , diameter_ratio_avg_(DIAMETER_RATIO_AVG)
    , diam_frac_lb_(DIAM_FRAC_LB)
    // parameters for critical_resource-driven growth
    , cr_neurite_{CRITICAL_GENERATED,
                  CRITICAL_GENERATED,
                  CRITICAL_GEN_TAU,
                  CRITICAL_DEL_TAU,
                  CRITICAL_GEN_VAR,
                  CRITICAL_GEN_CORR,
                  CRITICAL_GC_SUPPORT,
                  CRITICAL_SLOPE,
                  0,
                  0,
                  0,
                  0}
{
    // check if the growth cone model requires resource
    if (growth_cone_model_.find("resource-based") == 0)
    {
        use_critical_resource_ = true;
    }
    else
    {
        use_critical_resource_ = false;
    }

    uniform_   = std::uniform_real_distribution<double>(0., 1.);
    poisson_   = std::poisson_distribution<>(0);
    normal_    = std::normal_distribution<double>(0, 1);
    cr_normal_ = std::normal_distribution<double>(0, CRITICAL_GEN_VAR);
}


//! destructor
Neurite::~Neurite()
{
    // it's important to clear actin waves and growth cones
    // to use properly the SmartPointers
    growth_cones_.clear();
    growth_cones_tmp_.clear();
    growth_cones_inactive_.clear();
    growth_cones_inactive_tmp_.clear();
    dead_cones_.clear();
    actinDeck_.clear();
}


//#######################################################
//                  Initialization
//#######################################################

/**
 * @brief Initialize the first node of the Neurite
 *
 * It's called from neuron during neurite initialization,
 * it creates a new node with same position and attribute of the soma,
 * Now the branching algorithm is standardized for all the nodes
 *
 * @param soma
 * @param pos position of the soma
 * @param neurite_name
 */
void Neurite::init_first_node(BaseWeakNodePtr soma, const BPoint &pos,
                              const std::string &name, double soma_radius,
                              double neurite_diameter)
{
    auto firstNode = std::make_shared<Node>(
        soma, 0., pos, neurite_diameter, shared_from_this());

    firstNode->centrifugal_order_  = 0;
    firstNode->dist_to_soma_       = soma_radius;

    add_node(firstNode);

    // also initialize branching model
    if (branching_model_->neurite_ == nullptr)
    {
        branching_model_ = std::make_shared<Branching>(shared_from_this());
    }
}


void Neurite::finalize()
{
    // include the growth cones that were created during the previous substep
    if (growth_cones_tmp_.size() > 0)
    {
        growth_cones_.insert(growth_cones_tmp_.begin(),
                             growth_cones_tmp_.end());
    }

    growth_cones_tmp_.clear();
}


void Neurite::set_soma_angle(const double angle) { soma_angle_ = angle; }


double Neurite::get_soma_angle() const { return soma_angle_; }


/**
 * @brief store kernel parameters locally to avoid repeated calls
 *
 * This functions stores the time resolution from the kernel, it is called
 * every time the user changes that variable in the kernel.
 * It also updates these variables in all active growth cones.
 */
void Neurite::update_kernel_variables()
{
    for (auto &gc : growth_cones_)
    {
        gc.second->update_kernel_variables();
    }
}


//#######################################################
//                  Growth
//#######################################################

/**
 * @brief Neurite grow up with branching and elongation
 *
 * All the dynamical event of the neurite are computed here,
 * at this stage it manages for:
 * - growth cone splitting
 * - lateral branching event
 * - actin wave propagation
 * This function is called once for each simulation step
 *
 * @param rnd_engine
 */
void Neurite::grow(mtPtr rnd_engine, stype current_step, double substep)
{
    // include the growth cones that were created during the previous substep
    if (not growth_cones_tmp_.empty())
    {
        growth_cones_.insert(growth_cones_tmp_.begin(),
                             growth_cones_tmp_.end());
    }

    growth_cones_tmp_.clear();

    // call the branching model specific update
    update_growth_cones(rnd_engine, substep);

    // grow all the growth cones
    double diameter, b_length, total_b_length(0.);

    for (auto &gc : growth_cones_)
    {
        assert(gc.second.use_count() == 2);

        diameter = gc.second->get_diameter();

        // compute and check growth cones' diameters for stop
        if (diameter > min_diameter_)
        {
            try
            {
                gc.second->grow(rnd_engine, gc.first, substep);
            }
            catch (...)
            {
                std::throw_with_nested(
                    std::runtime_error("Passed from `Neurite::grow`."));
            }

            if (gc.second->stopped_ or gc.second->stuck_)
            {
                gc.second->current_stop_ += 1;
            }
            else
            {
                gc.second->current_stop_ = 0;
            }

            b_length = gc.second->get_branch()->get_length();
            total_b_length += b_length;

            diameter -= taper_rate_*gc.second->get_module();

            // negative diameter can be reach at the end of the growth if step
            // was too long
            if (diameter < 0)
            {
                BranchPtr bp = gc.second->get_branch();
                double old_length = bp->get_segment_length_at(bp->size() - 2) - bp->initial_distance_to_soma();
                printf("grew %f, to length %f, new diam %f; Delta diam %f - old length %f - computed old diam %f\n", gc.second->move_.module, b_length, diameter, taper_rate_*gc.second->move_.module, old_length, gc.second->get_diameter() + taper_rate_*(b_length - old_length));
                printf("gc speed was %f (avg %f) and substep %f\n", gc.second->move_.speed, gc.second->local_avg_speed_, substep);

                // compute where min_diameter was reached and retract up to
                // that position
                double retract = (min_diameter_ - diameter) / taper_rate_;
                printf("%s: length %f, retract %f\n", get_name().c_str(), b_length, retract);

                int omp_id     = kernel().parallelism_manager.get_thread_local_id();
                gc.second->retraction(retract, gc.first, omp_id);
                diameter = min_diameter_;
            }

            gc.second->set_diameter(diameter);
        }

        if (diameter <= min_diameter_ or
            gc.second->current_stop_ >= gc.second->max_stop_)
        {
            gc.second->active_ = false;
            growth_cones_inactive_tmp_.insert(gc);
        }
    }

    // move gc that are too thin to inactive
    if (not growth_cones_inactive_tmp_.empty())
    {
        for (auto &gc : growth_cones_inactive_tmp_)
        {
            growth_cones_inactive_.insert(gc);
            growth_cones_.erase(gc.first);

            // tell recorders that this growth cone became inactive
            // (this is equivalent to dying as far as recorders are concerned)
            kernel().record_manager.gc_died(parent_.lock()->get_gid(), name_,
                                            gc.first);
        }
    }

    growth_cones_inactive_tmp_.clear();

    // We remove the growth cones afterwards because death occurs while looping
    // on growth_cones_, hence we cannot modify it.
    // if dead_cones_ is not ordered, then this will fail!
    std::sort(dead_cones_.begin(), dead_cones_.end(), reverse_sorting);

    for (auto &cone_n : dead_cones_)
    {
        // check that no additional pointer exists (one in dead_cones_, one in
        // growth_cones_)
        assert(growth_cones_[cone_n].use_count() == 2);
        growth_cones_.erase(cone_n);
        // tell recorders that this growth cone died
        kernel().record_manager.gc_died(parent_.lock()->get_gid(), name_,
                                        cone_n);
    }

    dead_cones_.clear();

    // check total length
    if (fixed_arbor_len_ + total_b_length >= max_arbor_len_)
    {
        active_ = false;
    }
}


/**
 * @brief Update the growth cones depending on their model
 * @details Each model of neurite, like critical_resource has it's own
 * parameters to update, this function will be overriden by neurite's models
 *
 * @param rnd_engine
 */
void Neurite::update_growth_cones(mtPtr rnd_engine, double substep)
{
    // if using critical_resource model it's necessary to recompute the amount
    // of critical_resource required from each growth cone.
    //~ if (use_critical_resource_)
    //~ {
    //~ std::shared_ptr<GrowthCone_Critical> gcc;
    //~ cr_neurite_.tot_demand = 0;

    //~ stype i(0), j(0);
    //~ double ministep = 0.1;
    //~ double elapsed  = 0.;
    //~ double norm     = 1. / substep;
    //~ std::vector<double> avg_speed(growth_cones_.size(), 0);

    //~ while (elapsed < substep)
    //~ {
    //~ i += 1;
    //~ j = 0;

    //~ if (elapsed + ministep > substep)
    //~ {
    //~ ministep = substep - elapsed;
    //~ elapsed  = substep;
    //~ }
    //~ else
    //~ {
    //~ elapsed += ministep;
    //~ }

    //~ // compute the total demand
    //~ for (auto &gc : growth_cones_)
    //~ {
    //~ gcc = std::dynamic_pointer_cast<GrowthCone_Critical>(gc.second);
    //~ cr_neurite_.tot_demand += gcc->get_res_demand();
    //~ }

    //~ assert(cr_neurite_.tot_demand >= 0.);

    //~ if (cr_neurite_.tot_demand != 0)
    //~ {
    //~ cr_neurite_.tot_demand = 1. / cr_neurite_.tot_demand;
    //~ }

    //~ // compute the evolution of the resource
    //~ for (auto &gc : growth_cones_)
    //~ {
    //~ gcc = std::dynamic_pointer_cast<GrowthCone_Critical>(gc.second);
    //~ gcc->compute_CR(rnd_engine, ministep);
    //~ avg_speed[j] +=
    //~ gcc->compute_cr_speed(rnd_engine, ministep) * ministep;
    //~ j++;
    //~ }

    //~ // update the total resource
    //~ cr_neurite_.available =
    //~ cr_neurite_.available +
    //~ ministep * (cr_neurite_.eq_cr - cr_neurite_.available) /
    //~ cr_neurite_.tau +
    //~ sqrt(ministep) * cr_normal_(*(rnd_engine).get());
    //~ }

    //~ // set the gcs' speed to the average value
    //~ j = 0;

    //~ for (auto &gc : growth_cones_)
    //~ {
    //~ gc.second->move_.speed = avg_speed[j] * norm;
    //~ j++;
    //~ }
    //~ }
}


/**
 * Return the inverse of the total demand
 * (precomputed and inversed despite the misleading name)
 */
double Neurite::get_quotient_cr() const { return cr_neurite_.tot_demand; }


/**
 * Returns maximum amount of CR that can be delivered to a GC.
 */
double Neurite::get_available_cr() const
{
    return cr_neurite_.available / cr_neurite_.tau_delivery;
}


/**
 * Remove the parent node of two growth cone, one of which just died.
 * Transfer the remaining child to the grand-parent.
 *
 * @param parent    - the parent of the dead growth cone
 * @param child_id  - the id of the remaining, living growth cone child
 */
void Neurite::delete_parent_node(NodePtr parent, int living_child_id)
{
    auto child             = parent->children_[living_child_id];
    stype grand_parent_ID  = parent->get_parent().lock()->get_node_id();
    NodePtr grand_parent   = nodes_[grand_parent_ID];

    // reconcile the branch (remove parent length from fixed length)
    fixed_arbor_len_ -= parent->get_branch()->get_length();
    parent->branch_->append_branch(child->get_branch());

    // change parental relations
    for (stype i = 0; i < grand_parent->children_.size(); i++)
    {
        if (grand_parent->get_child(i)->get_node_id() == parent->get_node_id())
        {
            grand_parent->children_[i] = child;
            child->update_branch_and_parent(grand_parent,
                                            parent->branch_);
        }
    }

    // update tree structure
    update_tree_structure(grand_parent);

    nodes_.erase(parent->get_node_id());
    // assert(nodes_.size() == old_nodes_num - 1);
}


/**
 * Remove a growth cone that was absorbed back inside the parent branch,
 * then call :cpp:func:`delete_parent_node` to remove the lone parent
 * which has only one child.
 */
void Neurite::delete_cone(stype cone_n)
{
    // check if not already dead (can call this function several times in a
    // single timestep)
    GCPtr dead_cone = growth_cones_[cone_n];
    bool alive = not dead_cone->dead_;

    // delete only if not last growth cone (neurites cannot die)
    if (growth_cones_.size() - dead_cones_.size() > 1 and alive)
    {
        dead_cone->dead_ = true;

        assert(dead_cone->get_branch()->size() == 1);
#ifndef NDEBUG
        if (dead_cone->parent_.expired())
        {
            printf("invalid pointer coming\n");
        }
#endif

        stype parent_ID  = dead_cone->get_parent().lock()->get_node_id();
        auto parent_node = nodes_[parent_ID];
        dead_cones_.push_back(cone_n);

        // @todo: error in delete_parent_node (incorrect distances or diameters)

        //~ for (stype i = 0; i < parent_node->children_.size(); i++)
        //~ {
        //~ if (parent_node->get_child(i)->is_dead() == false)
        //~ {
        //~ delete_parent_node(parent_node, i);
        //~ }
        //~ }

        //~ assert(parent_node.use_count() == 1);
    }
}


//#######################################################
//                  Branching
//#######################################################

/**
 * @brief Compute the angle and the diameter of branching neurite
 * the output variables (old|new)_(diameter|angle) represent the branching cone
 * and the newly created cone
 * both are reference& to the growthcone and are actively used in the
 * computation too.
 * YOU NEED TO PASS THE ACTUAL ANGLE AND DIAMETER TO COMPUTE THE NEXT
 *
 * This function would require a clear model for the branching which relates the
 * angle to the diameter of the neurite.
 * Such a reference it's missing at this moment, excepet for this article: DOI
 * 10.1002/neu.20108
 *
 * @param rnd_engine random enginer
 * @param old_angle
 * @param new_angle
 * @param old_diameter
 * @param new_diameter
 */
void Neurite::gc_split_angles_diameter(mtPtr rnd_engine, double &old_angle,
                                       double &new_angle, double &old_diameter,
                                       double &new_diameter)
{
    // draw the branching angle from a gaussian distribution as hypothized in
    // reference article
    double branching_angle = gc_split_angle_mean_ +
                             gc_split_angle_std_ * normal_(*(rnd_engine).get());

    // ratio between the diameters of the two neurites,
    // it's a gaussian distributed value arround 1.
    double ratio;

    do
    {
        ratio = diameter_ratio_avg_ +
                diameter_ratio_std_ * normal_(*(rnd_engine).get());
    } while (ratio <= 0);

    // The diameters are computed on the base of: d^eta = d_1^eta + d_2^eta
    // looks weird but it's a simple optimization, next two lines compute the
    // new diameter
    // from the old diameter of the neurite and the ratio between the rising
    // cones.
    new_diameter = old_diameter * powf(1. + powf(ratio, diameter_eta_exp_),
                                       -1. / diameter_eta_exp_);
    old_diameter = powf(powf(old_diameter, diameter_eta_exp_) -
                            powf(new_diameter, diameter_eta_exp_),
                        1. / diameter_eta_exp_);

    // the diameter affect the branching angle: the largest cone goes straighter
    // than the other.
    double eps = (old_diameter - new_diameter) / (old_diameter + new_diameter);
    double half_tang = eps * tan(branching_angle / 2);
    new_angle        = old_angle;

    new_angle = -branching_angle / 2. - half_tang;
    old_angle = +branching_angle / 2. - half_tang;
}


/**
 * Set the parental relationships and add the new node to the neurite
 *
 * Adding the node sets its node_id_
 */
void Neurite::update_parent_nodes(NodePtr new_node, TNodePtr branching)
{
    // update parent node
    assert(new_node->get_parent().lock() == branching->get_parent().lock());

    NodePtr parent_node = nodes_[new_node->get_parent().lock()->get_node_id()];

    assert(parent_node->has_child() == true);

    for (stype i = 0; i < parent_node->children_.size(); i++)
    {
        if (parent_node->get_child(i) == branching)
        {
            parent_node->children_[i] = new_node;
            break;
        }
    }

    branching->parent_   = new_node;
    new_node->has_child_ = true;

    add_node(new_node);
}


/**
 * Create a new growth cone from a new parent node `new_node` created after a
 * branching event.
 * The new growth cone is a at a distance `new_length` of the parent `new_node`.
 */
GCPtr Neurite::create_branching_cone(const TNodePtr branching_node,
                                     NodePtr new_node, double dist_to_parent,
                                     double new_diameter, const BPoint &xy,
                                     double new_cone_angle, bool split)
{
    // check that the new position is inside the environment and does not
    // overlap with another branch from the same element.
    BPoint p(cos(new_cone_angle) * dist_to_parent,
             sin(new_cone_angle) * dist_to_parent);
    bg::add_point(p, xy);

    if (new_diameter - taper_rate_ * dist_to_parent <= 0)
    {
        return nullptr;
    }

    if (not kernel().space_manager.env_contains(p))
    {
        return nullptr;
    }

    if (kernel().space_manager.interactions_on())
    {
        if (std::isnan(growth_cones_.begin()->second->get_self_affinity()))
        {
            // check that the target point is not inside its own neurite
            BPolygon poly;

            if (kernel().space_manager.is_inside(
                    p, parent_.lock()->get_gid(), name_,
                    0.5 * nodes_[0]->get_diameter(), poly))
            {
                return nullptr;
            }
        }
    }

    // create a new growth cone
    GCPtr sibling = growth_cones_.begin()->second->clone(
        new_node, shared_from_this(), dist_to_parent, p, new_cone_angle);

    // Here we copy model and status from a random growth cone
    // in the neurite since all the growth cones have same status and
    // model. It's not possible to copy from the splitting cone since this
    // function is common with lateral branching
    // growth_cones_.back()->get_status(status);
    // sibling->set_status(status);
    // insert new elements in the tree
    new_node->children_.push_back(branching_node);
    new_node->children_.push_back(sibling);

    // make the branch (start from parent branch then add new point just
    // outside the old branch radius)
    double parent_to_soma = new_node->get_distance_to_soma();
    BranchPtr b           = std::make_shared<Branch>(xy, parent_to_soma);

    int omp_id      = kernel().parallelism_manager.get_thread_local_id();
    ObjectInfo info = std::make_tuple(parent_.lock()->get_gid(), name_,
                                      num_created_nodes_, 0);

    // update diameter properties
    new_diameter -= dist_to_parent*taper_rate_;
    sibling->set_diameter(new_diameter);

    // update branch and assign it to the new growth cone
    // this is only done for lateral branching; for splitting, the procedure
    // is more complex and has to be performed after the R-tree update, so
    // in Branching::update_splitting_cones, called from
    // Branching::branching_event after the call to this function.
    // @todo see whether both things cannot be done at the same place to
    // simplify things
    if (not split)
    {
        // update the parent and the branch
        sibling->update_branch_and_parent(new_node, b);

        try
        {
            kernel().space_manager.add_object(xy, p, new_diameter,
                                              dist_to_parent, taper_rate_, info,
                                              b, omp_id);
        }
        catch (...)
        {
            std::throw_with_nested(std::runtime_error(
                "Passed from `Neurite::create_branching_cone`."));
        }
    }

    add_cone(sibling);

    return sibling;
}


/**
 * @brief Branch from a node of the neuritic tree
 *
 * The branching event can happen wherever along the branch of `branching_node`.
 * Both the internal nodes and the leaves (the GrowthCones) can have a
 * lateral branching event.
 * This function will create a new node at the branching point `branch_point`
 * along the branch of `branching_point`.
 * A new Node will be created there and a new growth cone will start from this
 * node, cloned from the Neurite's default model in 'models_manager'.
 *
 * @param branching_node the node which is going to branch
 * @param branchpoint index of the point in the branch where the
 *    branching occurs.
 * @param new_node (passed as nullptr) will be created if branching can
 *     occur.
 * @param new_length the length of the newborn GrowthCone
 * @param rnd_engine
 */
bool Neurite::lateral_branching(TNodePtr branching_node, stype branch_point,
                                NodePtr &new_node, mtPtr rnd_engine)
{
    if (not branching_node->is_dead() and active_)
    {
        // Locate the event and the event parameters
        double branch_direction(0), distance_to_soma(0);
        BranchPtr branch = branching_node->get_branch();
        BPoint xy;

        locate_from_idx(xy, branch_direction, distance_to_soma, branch,
                        branch_point);

        char branching_side = 1;
        if (uniform_(*(rnd_engine).get()) < 0.5)
        {
            branching_side = -1;
        }

        double angle =
            branching_side * (lateral_branching_angle_mean_) +
            lateral_branching_angle_std_ * normal_(*(rnd_engine).get());

        // compute the distance from branching_node (towards the soma)
        double distance_from_branching =
            branch->final_distance_to_soma() - distance_to_soma;

        double distance_to_parent =
            distance_to_soma - branch->initial_distance_to_soma();

        // compute the local diameter on the branch
        double new_diam = branching_node->get_diameter()
                          + taper_rate_*distance_from_branching;

        // create the new node at position xy
        new_node = std::make_shared<Node>(
            branching_node->parent_, distance_to_parent, xy, new_diam,
            shared_from_this());

        // create the new growth cone just outside the local diameter
        // to do so, first get the polygons around the branching point
        // get the properties of the neighboring geometries
        BPoint line_end = BPoint(2 * new_diam * cos(branch_direction + angle),
                                 2 * new_diam * sin(branch_direction + angle));
        bg::add_point(line_end, xy);
        BLineString ls({xy, line_end});

        std::vector<ObjectInfo> neighbors_info;
        // kernel().space_manager.get_objects_in_range(line_end, 1.5*new_diam,
        //                                             neighbors_info);
        kernel().space_manager.get_intersected_objects(xy, line_end,
                                                       neighbors_info);

        stype gid = parent_.lock()->get_gid();
        stype nid = branching_node->get_node_id();

        BMultiPoint mp;

        for (auto info : neighbors_info)
        {
            if (std::get<ndata::NEURON>(info) == gid and
                std::get<ndata::NEURITE>(info) == name_ and
                std::get<ndata::NODE>(info) == nid)
            {
                stype seg_id = std::get<ndata::SEGMENT>(info);

                // check if the segment still exists (the gc may be retracting
                // just next to the branching point)
                if (seg_id + 1 < branch->size())
                {
                    bg::intersection(*(branch->get_segment_at(seg_id).get()),
                                     ls, mp);
                }
            }
        }

        double max_distance = std::numeric_limits<double>::min();
        double dist;

        for (BPoint p : mp)
        {
            dist = bg::distance(xy, p);

            if (dist > max_distance)
            {
                max_distance = dist;
            }
        }

        // get distance (a little more to avoid being inside)
        double dist_to_bp = 1.05 * max_distance;
        // @todo set diameter fraction
        double new_cone_diam = diam_frac_lb_ * new_diam;

        auto sibling = create_branching_cone(branching_node, new_node,
                                             dist_to_bp, new_cone_diam, xy,
                                             branch_direction + angle, false);

        // check if sibling could indeed be created
        if (sibling != nullptr)
        {
            // this will add the nodes to the neurite and set their node_id_
            update_parent_nodes(new_node, branching_node);

            // update the existing node
            new_node->branch_ = std::make_shared<Branch>(
                *(branching_node->branch_.get()));

            branching_node->branch_ =
                new_node->branch_->resize_head(branch_point);

            branching_node->set_first_point(xy, distance_to_soma);

            // update the new node's branch
            // NB2: the + 1 is necessary for the segment between the branches
            new_node->get_branch()->resize_tail(branch_point + 1);

            // modify subsequent nodes
            update_tree_structure(new_node);

            assert(sibling->get_centrifugal_order() ==
                   branching_node->get_centrifugal_order());
            assert(new_node->get_child(0) == branching_node);
            assert(new_node->get_child(1) == sibling);

            return true;
        }
    }

    return false;
}


/**
 * Manage the "growth cone splitting" branching event, computing the new
 * directions of the two GrowthCone objects.
 *
 * GrowthCone splitting means that an existing GrowthCone will divide into
 * two new GrowthCone objects.
 * This results in a new TopologicalNode being created at the position
 * where the split occured. This new node becomes the new parent of the
 * existing GrowthCone and of the new GrowthCone that is created.
 */
bool Neurite::growth_cone_split(GCPtr branching_cone, double new_length,
                                double new_angle, double old_angle,
                                double new_diameter, double old_diameter,
                                NodePtr &new_node, GCPtr &sibling)
{
    if (not branching_cone->is_dead() and active_)
    {
        double direction = branching_cone->move_.angle;

        // prepare growth cone variables for split
        branching_cone->prepare_for_split();

        // create new node as branching point
        new_node = std::make_shared<Node>(
            branching_cone->parent_,
            branching_cone->get_branch_length(),
            branching_cone->get_position(),
            branching_cone->get_diameter(),
            shared_from_this());

        // create the sibling
        sibling = create_branching_cone(
            branching_cone, new_node, new_length, new_diameter,
            branching_cone->get_position(), direction + new_angle, true);

        if (sibling != nullptr)
        {
            // this will add the nodes to the neurite and set their node_id_
            update_parent_nodes(new_node, branching_cone);

            // move the old growth cone --> growth cone split specific
            branching_cone->move_.angle = direction + old_angle;
            branching_cone->set_diameter(old_diameter);
            branching_cone->topological_advance();

            // copy the existing node branch to new_node
            new_node->branch_ = std::make_shared<Branch>(
                *(branching_cone->branch_.get()));

            // make a new empty branch for the existing branching_cone
            branching_cone->branch_ = std::make_shared<Branch>(
                branching_cone->get_position(),
                new_node->get_branch()->final_distance_to_soma());

            // update growth cones' variables after split
            branching_cone->after_split();
            sibling->after_split();

            assert(sibling->get_centrifugal_order() ==
                   branching_cone->get_centrifugal_order());
            assert(new_node->get_child(0) == branching_cone);
            assert(nodes_[new_node->get_node_id()]->has_child() == true);

            return true;
        }
    }

    return false;
}


void Neurite::update_tree_structure(TNodePtr root)
{
    std::deque<TNodePtr> nodes{root};
    root->centrifugal_order_ =
        root->parent_.lock()->get_centrifugal_order() + 1;

    while (not nodes.empty())
    {
        TNodePtr node = nodes.front();
        nodes.pop_front();
        if (node->has_child())
        {
            NodePtr mynode = std::dynamic_pointer_cast<Node>(node);
            for (stype i = 0; i < mynode->children_.size(); i++)
            {
                mynode->children_[i]->centrifugal_order_ =
                    mynode->get_centrifugal_order() + 1;
                nodes.push_back(mynode->children_[i]);
            }
        }
    }
}


//#######################################################
//              Actin Model
//#######################################################

void Neurite::add_actin(ActinPtr actin) { actinDeck_.push_front(actin); }


void Neurite::update_actin_waves(mtPtr rnd_engine, double substep)
{
    // count the number of actin elements because we only want to loop over
    // them once (they put themselves or new AW back in the deque in step.
    unsigned long aw_count = 0;
    auto size              = actinDeck_.size();
    while (aw_count < size)
    {
        aw_count++;
        auto aw = actinDeck_.front();
        aw->step(rnd_engine, substep);
        actinDeck_.pop_front();
    }
}


void Neurite::start_actin_wave(double actin_content)
{
    actinDeck_.push_back(std::make_shared<ActinWave>(
        get_first_node(), actin_content, shared_from_this()));
}


//#######################################################
//              Utilities Functions
//#######################################################

unsigned int Neurite::num_growth_cones() const
{
    return growth_cones_.size() + growth_cones_tmp_.size();
}


/**
 * @brief Add a GrowthCone to the neurite
 *
 * @param GCPtr pointer to the GrowthCone
 */
void Neurite::add_cone(GCPtr cone)
{
    cone->node_id_ = num_created_nodes_;

    // first gc (only soma node exists) is added directly
    if (num_created_nodes_ == 1)
    {
        growth_cones_[num_created_nodes_] = cone;
    }
    else
    {
        growth_cones_tmp_[num_created_nodes_] = cone;
    }

    num_created_nodes_++;

    // if using critical resource model, update total production
    if (use_critical_resource_)
    {
        double ratio  = num_created_nodes_ / cr_neurite_.typical_gc_support;
        double factor = cr_neurite_.increase_slope *
                        cr_neurite_.typical_gc_support * cr_neurite_.target_cr;

        cr_neurite_.eq_cr = (cr_neurite_.tau / cr_neurite_.tau_generation) *
                            factor * std::tanh(ratio);
    }

    // check if we do not reach the growth cone number limit
    if (growth_cones_.size() >= max_gc_num_)
    {
        active_ = false;
    }
}


/**
 * @brief Add a node to the neurite
 *
 * @param NodePtr pointer to the Node
 */
void Neurite::add_node(NodePtr node)
{
    node->node_id_ = num_created_nodes_;

    nodes_[num_created_nodes_] = node;

    num_created_nodes_++;

    // update the fixed arbor length (add only if all children are GCs)
    bool all_gcs = true;

    for (auto child : node->children_)
    {
        if (child->has_child())
        {
            all_gcs = false;
        }
    }

    if (all_gcs)
    {
        fixed_arbor_len_ += node->get_branch()->get_length();
    }
}


bool Neurite::walk_tree(NodeProp &np) const
{
    static bool reset   = false;
    static auto gcrange = gc_range();
    static auto gc_it   = gcrange.begin();
    static auto n_it    = nodes_.cbegin();

    if (reset)
    {
        gcrange = gc_range();
        gc_it   = gc_range().begin();
        n_it    = nodes_.cbegin();

        reset = false;
    }

    stype nid, pid;
    double diam, dtp;

    if (n_it != nodes_.cend())
    {
        // get node id
        nid = n_it->second->get_node_id();

        // get diameter (average between parent and current node)
        // get_diameter returns the value at the root of the branch.
        diam = n_it->second->get_diameter();

        if (n_it->first != 0)
        {
            diam -=
                0.5 * taper_rate_ * n_it->second->get_branch()->get_length();
            // get parent id
            pid = n_it->second->get_parent().lock()->get_node_id();
        }
        else
        {
            // no real parent for soma
            pid = 0;
        }

        // get distance to parent
        dtp = n_it->second->get_distance_parent();

        if (dtp != n_it->second->get_branch_length())
        {
            printf("node %lu dtp %f != branch length %f\n", n_it->first, dtp, n_it->second->get_branch_length());
            std::cout << bg::wkt(n_it->second->get_position()) << std::endl;
        }

        // get position (root of the branch)
        BPoint p = n_it->second->get_position();
        std::vector<double> coords({p.x(), p.y()});

        np = NodeProp(nid, pid, diam, dtp, coords);

        n_it++;

        return true;
    }
    else if (gc_it != gcrange.end())
    {
        // get node id
        nid = gc_it->second->get_node_id();

        // get parent id
        // @todo check for nullptr error
        if (gc_it->second->get_parent().lock() != nullptr)
        {
            pid = gc_it->second->get_parent().lock()->get_node_id();
        }
        else
        {
            pid = nid;
        }

        // get distance to parent
        dtp = gc_it->second->get_branch_length();

        if (dtp != gc_it->second->get_branch_length())
            printf("gc dtp %f != branch length %f\n", dtp, gc_it->second->get_distance_parent());

        // get diameter (average between root and tip)
        diam = gc_it->second->get_diameter() + 0.5*dtp*taper_rate_;
        // get position
        BPoint p = gc_it->second->get_position();
        std::vector<double> coords({p.x(), p.y()});

        np = NodeProp(nid, pid, diam, dtp, coords);

        gc_it++;

        return true;
    }

    // reset iterators
    reset = true;

    return false;
}


simple_gc_range Neurite::active_gc_range() const
{
    return boost::make_iterator_range(growth_cones_);
}


joined_gc_range Neurite::gc_range() const
{
    return boost::range::join(growth_cones_, growth_cones_inactive_);
}


std::unordered_map<stype, NodePtr>::const_iterator Neurite::nodes_cbegin() const
{
    return nodes_.cbegin();
}


std::unordered_map<stype, NodePtr>::const_iterator Neurite::nodes_cend() const
{
    return nodes_.cend();
}


NodePtr Neurite::get_first_node() const { return nodes_.at(0); }


NeuronWeakPtr Neurite::get_parent_neuron() const { return parent_; }


const std::string &Neurite::get_name() const { return name_; }


const std::string &Neurite::get_type() const { return neurite_type_; }


void Neurite::get_distances(stype node, stype segment, double &dist_to_parent,
                            double &dist_to_soma) const
{
    TNodePtr tnode;

    auto n = nodes_.find(node);

    if (n == nodes_.end())
    {
        auto gc = growth_cones_.find(node);
        tnode   = gc->second;
    }
    else
    {
        tnode = n->second;
    }

    BranchPtr branch         = tnode->get_branch();
    double init_dist_to_soma = branch->initial_distance_to_soma();

    dist_to_soma   = branch->at(segment)[2];
    dist_to_parent = dist_to_soma - init_dist_to_soma;
}


//###################################################
//                  get/set status
//###################################################
//

void Neurite::set_status(const statusMap &status)
{
    get_param(status, names::max_gc_number, max_gc_num_);
    get_param(status, names::max_arbor_length, max_arbor_len_);

    get_param(status, names::diameter_fraction_lb, diam_frac_lb_);

    get_param(status, names::diameter_eta_exp, diameter_eta_exp_);
    get_param(status, names::gc_split_angle_mean, gc_split_angle_mean_);
    get_param(status, names::gc_split_angle_std, gc_split_angle_std_);

    double tr;
    bool tr_set = get_param(status, names::taper_rate, tr);

    if (tr_set and tr != taper_rate_)
    {
        if (kernel().simulation_manager.get_current_minutes() > 0.)
        {
            throw std::invalid_argument("Cannot change `taper_rate` after "
                                        "simulation start.");
        }

        taper_rate_ = tr;
    }

    double sd_std(diameter_ratio_std_), sd_avg(diameter_ratio_avg_);
    get_param(status, names::diameter_ratio_std, sd_std);
    if (sd_std < 0)
    {
        throw InvalidArg("`diameter_ratio_std` must be positive.", __FUNCTION__,
                         __FILE__, __LINE__);
    }
    diameter_ratio_std_ = sd_std;
    get_param(status, names::diameter_ratio_avg, sd_avg);
    if (sd_avg < 0)
    {
        throw InvalidArg("`diameter_ratio_avg` must be positive.", __FUNCTION__,
                         __FILE__, __LINE__);
    }
    diameter_ratio_avg_ = sd_avg;

    get_param(status, names::lateral_branching_angle_mean,
              lateral_branching_angle_mean_);

    get_param(status, names::lateral_branching_angle_std,
              lateral_branching_angle_std_);

    if (branching_model_->neurite_ == nullptr)
    {
        branching_model_ = std::make_shared<Branching>(shared_from_this());
    }

    //                 Critical Resource Params
    //###################################################

    if (use_critical_resource_)
    {
        get_param(status, names::res_neurite_generated, cr_neurite_.target_cr);
        get_param(status, names::res_neurite_variance, cr_neurite_.var);
        get_param(status, names::res_neurite_generated_tau,
                  cr_neurite_.tau_generation);
        get_param(status, names::res_neurite_delivery_tau,
                  cr_neurite_.tau_delivery);
        get_param(status, names::res_typical_gc_support,
                  cr_neurite_.typical_gc_support);
        get_param(status, names::res_increase_slope,
                  cr_neurite_.increase_slope);
        // get_param(status, names::res_neurite_correlation, cr_neurite_.tau);
        //
        // optimize variable for less computation:
        cr_neurite_.tau = 1. / (1. / cr_neurite_.tau_delivery +
                                1. / cr_neurite_.tau_generation);
        double ratio    = num_created_nodes_ / cr_neurite_.typical_gc_support;
        double factor   = cr_neurite_.increase_slope *
                        cr_neurite_.typical_gc_support * cr_neurite_.target_cr;

        cr_neurite_.eq_cr = (cr_neurite_.tau / cr_neurite_.tau_generation) *
                            factor * std::tanh(ratio);
        // @TODO this is not obvious!!
        cr_neurite_.available = cr_neurite_.eq_cr;

        cr_normal_ = std::normal_distribution<>(0, cr_neurite_.var);
#ifndef NDEBUG
        printf("\n"
               " CRITICAL RESOURCE BRANCHING \n"
               "%s : %f \n"
               "%s : %f \n"
               "%s : %f \n",
               names::res_neurite_available.c_str(), cr_neurite_.available,
               names::gc_split_angle_mean.c_str(),
               gc_split_angle_mean_ * 180 / M_PI,
               names::gc_split_angle_std.c_str(),
               gc_split_angle_std_ * 180 / M_PI);
#endif
    }

    if (use_critical_resource_)
    {
        observables_.push_back("A"); // add A as observable
    }
    else if (observables_.back() == "A")
    {
        observables_.pop_back(); // delete A as observable (last)
    }

    branching_model_->set_status(status);

    for (auto &gc : growth_cones_)
    {
        gc.second->set_status(status);
    }
    for (auto &gc : growth_cones_inactive_)
    {
        gc.second->set_status(status);
    }
}


//@TODO
void Neurite::get_status(statusMap &status, const std::string &level) const
{
    if (level == "neurite")
    {
        set_param(status, names::neurite_type, neurite_type_, "");
        set_param(status, names::diameter_fraction_lb, diam_frac_lb_, "");
        set_param(status, names::gc_split_angle_mean, gc_split_angle_mean_,
                  "rad");
        set_param(status, names::gc_split_angle_std, gc_split_angle_std_,
                  "rad");
        set_param(status, names::lateral_branching_angle_std,
                  lateral_branching_angle_std_, "rad");
        set_param(status, names::lateral_branching_angle_mean,
                  lateral_branching_angle_mean_, "rad");
        set_param(status, names::observables, observables_, "");

        set_param(status, names::taper_rate, taper_rate_, "1 / micrometer");
        set_param(status, names::diameter_ratio_std, diameter_ratio_std_, "");
        set_param(status, names::diameter_ratio_avg, diameter_ratio_avg_, "");
        set_param(status, names::diameter_eta_exp, diameter_eta_exp_, "");

        set_param(status, names::max_gc_number, max_gc_num_, "");
        set_param(status, names::max_arbor_length, max_arbor_len_,
                  "micrometer");
        set_param(status, names::active, active_, "");

        // branching properties
        branching_model_->get_status(status);

        // critical resource properties
        if (use_critical_resource_)
        {
            set_param(status, names::res_neurite_generated,
                      cr_neurite_.target_cr, "micromole / liter");
            set_param(status, names::res_neurite_variance, cr_neurite_.var,
                      "micromole / liter / minute**0.5");
            set_param(status, names::res_neurite_generated_tau,
                      cr_neurite_.tau_generation, "minute");
            set_param(status, names::res_neurite_delivery_tau,
                      cr_neurite_.tau_delivery, "minute");
            set_param(status, names::res_typical_gc_support,
                      cr_neurite_.typical_gc_support, "");
            set_param(status, names::res_increase_slope,
                      cr_neurite_.increase_slope, "");
        }
    }

    // growth_cone properties
    set_param(status, names::growth_cone_model, growth_cone_model_, "");
    set_param(status, names::taper_rate, taper_rate_, "");

    if (not growth_cones_.empty())
    {
        growth_cones_.begin()->second->get_status(status);
    }
    else
    {
        growth_cones_inactive_.begin()->second->get_status(status);
    }
}


/**
 * @brief Get the current value of one of the observables
 */
double Neurite::get_state(const std::string &observable) const
{
    double value = 0.;

    if (observable == "A")
    {
        value = cr_neurite_.available;
    }
    else
    {
        // default case, just sum up
        for (const auto &gc : growth_cones_)
        {
            value += gc.second->get_state(observable);
        }
        for (const auto &gc : growth_cones_inactive_)
        {
            value += gc.second->get_state(observable);
        }
    }

    return value;
}


double Neurite::get_taper_rate() const { return taper_rate_; }


double Neurite::get_max_resol() const
{
    // check speed limit
    if (not growth_cones_.empty())
    {
        double max_speed = growth_cones_.begin()->second->avg_speed_ +
                           5 * growth_cones_.begin()->second->speed_variance_;
        double length = growth_cones_.begin()->second->filopodia_.finger_length;

        return length / max_speed;
    }

    return std::numeric_limits<double>::infinity();
}

} // namespace growth
