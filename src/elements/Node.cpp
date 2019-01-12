#include "Node.hpp"

#include <memory>

#include "exceptions.hpp"

#include "Branch.hpp"
#include "GrowthCone.hpp"


namespace growth
{

BaseNode::BaseNode() : geometry_{}
{}

BaseNode::BaseNode(const BPoint &position, double distance_to_parent,
         double distance_to_soma)
  : geometry_{position, distance_to_soma, distance_to_parent}
{}

BaseNode::BaseNode(const BaseNode &copy)
  : geometry_(copy.geometry_)
{}

// BaseNode functions

void BaseNode::set_position(const BPoint &pos) { geometry_.position = pos; }


double BaseNode::get_distance_to_soma() const { return geometry_.dis_to_soma; }


BPoint BaseNode::get_position() const { return geometry_.position; }


int BaseNode::get_centrifugal_order() const { return -1; }


size_t BaseNode::get_nodeID() const { return 0; }


double BaseNode::get_distance_parent() const
{
    return geometry_.dis_to_parent;
};


// TopologicalNode functions
//

TopologicalNode::TopologicalNode()
    : BaseNode()
    , topology_{}
    , biology_{}
{}


TopologicalNode::TopologicalNode(const TopologicalNode &tnode)
    : BaseNode(tnode)
    , topology_(tnode.topology_)
    , biology_(tnode.biology_)
{
    topology_.has_child = false;
}


TopologicalNode::TopologicalNode(BaseWeakNodePtr parent, float distanceToParent,
                                 const BPoint &position)
    : BaseNode(position, distanceToParent,
               distanceToParent + parent.lock()->get_distance_to_soma())
    , topology_{parent, parent.lock()->get_centrifugal_order() + 1, false, 0}
    , biology_{false, std::make_shared<Branch>(), nullptr, 1}
{}


void TopologicalNode::topological_advance()
{
    topology_.centrifugal_order++;
}


/**
 * @brief Overwrite the first element of the owned branch
 * This function update the distanceToSoma of the TNode and
 * call the omonoom function in the Branch
 * @param xy
 * @param distanceToSoma
 */
void TopologicalNode::set_first_point(const BPoint &p, double length)
{
    biology_.branch->set_first_point(p, length);
}


void TopologicalNode::set_position(const BPoint &pos)
{
    LogicError("Cannot set only position for a TopologicalNode", __FUNCTION__,
               __FILE__, __LINE__);
}


void TopologicalNode::set_position(const BPoint &pos, double dist_to_soma,
                                   BranchPtr b)
{
    assert(pos.x() == b->get_last_xy().x() and
           pos.y() == b->get_last_xy().y());

    geometry_.position    = pos;
    geometry_.dis_to_soma = dist_to_soma;
    geometry_.dis_to_parent =
        dist_to_soma - topology_.parent.lock()->get_distance_to_soma();

    biology_.branch = b;
}


void TopologicalNode::set_diameter(double diameter)
{
    biology_.diameter = diameter;
}


size_t TopologicalNode::get_branch_size() const
{
    return biology_.branch->size();
}


seg_range TopologicalNode::segment_range() const
{
    return biology_.branch->segment_range();
}


bool TopologicalNode::support_AW() const { return false; }

//##########################################################
//                   Node functions
//##########################################################
/**
This parameters
Growth cone constructor:
\param parent_ is the parent node in the dendritic tree.
\param distanceToParent is the distance in micrometers from the parent node.
\param OwnNeurite is the neurite which the GC belongs to
\param position is the real space point of the node
\param BranchParentID is the relative name respect to the parent,
necessary for build the Id
*/
Node::Node(BaseWeakNodePtr parent, float distance, const BPoint &pos)
  : TopologicalNode(parent, distance, pos)
{
    topology_.has_child = true;
}


/**lateral_branching_angle_mean
 * Modified Node copy constructor,
 * used to copy a GrowthCone into a standard Node after a branching
 * event.
 */
Node::Node(const TopologicalNode &copyTopoNode)
  : TopologicalNode(copyTopoNode)
{
}


Node::Node(const Node &copyNode)
  : TopologicalNode(copyNode)
{
    topology_.has_child = true;
}


TNodePtr Node::get_child(int index) const { return children_.at(index); }

} // namespace growth
