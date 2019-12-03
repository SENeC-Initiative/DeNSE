/*
 * Node.cpp
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

#include "Node.hpp"

#include <memory>

#include "exceptions.hpp"

#include "Branch.hpp"
#include "GrowthCone.hpp"


namespace growth
{

BaseNode::BaseNode()
  : position_(BPoint())
  , dist_to_soma_(-1)
  , dist_to_parent_(-1)
{}

BaseNode::BaseNode(const BPoint &position, double distance_to_parent,
         double distance_to_soma)
  : position_(position)
  , dist_to_soma_(distance_to_soma)
  , dist_to_parent_(distance_to_parent)
{}

BaseNode::BaseNode(const BaseNode &copy)
  : position_(copy.position_)
  , dist_to_soma_(copy.dist_to_soma_)
  , dist_to_parent_(copy.dist_to_parent_)
{}

// BaseNode functions

void BaseNode::set_position(const BPoint &pos) { position_ = pos; }


double BaseNode::get_distance_to_soma() const { return dist_to_soma_; }


BPoint BaseNode::get_position() const { return position_; }


int BaseNode::get_centrifugal_order() const { return -1; }


stype BaseNode::get_node_id() const { return 0; }


double BaseNode::get_distance_parent() const
{
    return dist_to_parent_;
};


// TopologicalNode functions

TopologicalNode::TopologicalNode()
    : BaseNode()
    // topological properties
    , parent_(std::make_shared<BaseNode>())
    , centrifugal_order_(-1)
    , has_child_(false)
    , node_id_(0)
    // biological properties
    , dead_(false)
    , branch_(std::make_shared<Branch>())
    , own_neurite_(nullptr)
    , diameter_(0)
{}


TopologicalNode::TopologicalNode(const TopologicalNode &tnode)
    : BaseNode(tnode)
    , parent_(tnode.parent_)
    , centrifugal_order_(tnode.centrifugal_order_)
    , has_child_(false)
    , node_id_(0)
    , dead_(false)
    , branch_(std::make_shared<Branch>())
    , own_neurite_(tnode.own_neurite_)
    , diameter_(tnode.diameter_)
{}


TopologicalNode::TopologicalNode(BaseWeakNodePtr parent,
                                 double distance_to_parent,
                                 const BPoint &position,
                                 double diameter,
                                 NeuritePtr neurite)
    : BaseNode(position, distance_to_parent,
               distance_to_parent +
               parent.lock()->get_distance_to_soma())
    , parent_(parent)
    , centrifugal_order_(parent.lock()->get_centrifugal_order() + 1)
    , has_child_(false)
    , node_id_(0)
    , dead_(false)
    , branch_(std::make_shared<Branch>())
    , own_neurite_(neurite)
    , diameter_(diameter)
{}


void TopologicalNode::topological_advance()
{
    centrifugal_order_++;
}


void TopologicalNode::topological_advance() { topology_.centrifugal_order++; }


/**
 * @brief Overwrite the first element of the owned branch
 * This function update the dist_to_parent of the TNode and
 * call the homonym function in the Branch
 * @param xy
 * @param dist_to_soma
 */
void TopologicalNode::set_first_point(const BPoint &p,
                                      double dist_to_soma)
{
    dist_to_parent_ = dist_to_soma_ - dist_to_soma;

    branch_->set_first_point(p, dist_to_soma);
}


void TopologicalNode::set_position(const BPoint &pos)
{
    LogicError("Cannot set only position for a TopologicalNode",
               __FUNCTION__, __FILE__, __LINE__);
}


void TopologicalNode::update_branch_and_parent(BaseNodePtr parent,
                                               BranchPtr b)
{
    // update parent
    parent_ = parent;

    double old_dts = dist_to_parent_;

    dist_to_parent_ =
        dist_to_soma_ - parent->get_distance_to_soma();

    if (dist_to_parent_ < 0)
        printf("Updated a node with negative dts %f from %f; source were %f and %f\n", dist_to_parent_, old_dts, dist_to_soma_, parent->get_distance_to_soma());

    branch_ = b;
}


void TopologicalNode::set_diameter(double diameter)
{
    diameter_ = diameter;
}


stype TopologicalNode::get_branch_size() const
{
    return branch_->size();
}


double TopologicalNode::get_branch_length() const
{
    return branch_->get_length();
}


seg_range TopologicalNode::segment_range() const
{
    return branch_->segment_range();
}


bool TopologicalNode::support_AW() const { return false; }

//##########################################################
//                   Node functions
//##########################################################
/**
This parameters
Growth cone constructor:
\param parent_ is the parent node in the dendritic tree.
\param distance_to_parent is the distance in micrometers from the parent node.
\param position is the real space point of the node
\param BranchParentID is the relative name respect to the parent,
necessary for build the Id
*/
Node::Node(BaseWeakNodePtr parent, double distance_to_parent,
           const BPoint &pos, double diameter, NeuritePtr neurite)
  : TopologicalNode(parent, distance_to_parent, pos, diameter, neurite)
{
    has_child_ = true;
}


TNodePtr Node::get_child(int index) const { return children_.at(index); }

} // namespace growth
