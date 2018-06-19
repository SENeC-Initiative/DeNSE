#include "Node.hpp"

#include <memory>

#include "exceptions.hpp"

#include "Branch.hpp"
#include "GrowthCone.hpp"


namespace growth
{

// BaseNode functions

void BaseNode::set_position(const Point &pos) { geometry_.position = pos; }


double BaseNode::get_distance_to_soma() const { return geometry_.dis_to_soma; }


Point BaseNode::get_position() const { return geometry_.position; }


int BaseNode::get_centrifugal_order() const { return -1; }


std::string BaseNode::get_treeID() const { return somaID_; }


size_t BaseNode::get_nodeID() const { return 0; }


double BaseNode::get_distance_parent() const
{
    return geometry_.dis_to_parent;
};


// TopologicalNode functions
//

TopologicalNode::TopologicalNode()
    : topology_{}
    , geometry_{}
    , biology_{}
{
}


TopologicalNode::TopologicalNode(const TopologicalNode &tnode)
    : topology_(tnode.topology_)
    , geometry_(tnode.geometry_)
    , biology_(tnode.biology_)
{
    topology_.has_child = false;
}


TopologicalNode::TopologicalNode(BaseWeakNodePtr parent, float distanceToParent,
                                 const Point &position,
                                 const std::string &binaryID)
    : topology_{parent, parent.lock()->get_centrifugal_order() + 1, false, 0,
                binaryID}
    , geometry_{position, distanceToParent,
                parent.lock()->get_distance_parent() +
                    parent.lock()->get_distance_to_soma()}
    , biology_{false, std::make_shared<Branch>(), nullptr, 1}

{
}


void TopologicalNode::topological_advance()
{
    topology_.binaryID = topology_.binaryID + "0";
    topology_.centrifugal_order++;
}


/**
 * @brief Overwrite the first element of the owned branch
 * This function update the distanceToSoma of the TNode and
 * call the omonoom function in the Branch
 * @param xy
 * @param distanceToSoma
 */
void TopologicalNode::set_first_point(const Point p, double length)
{
    biology_.branch->set_first_point(p, length);
}


void TopologicalNode::set_position(const Point &pos)
{
    LogicError("Cannot set only position for a TopologicalNode", __FUNCTION__,
               __FILE__, __LINE__);
}


void TopologicalNode::set_position(const Point &pos, double dist_to_soma,
                                   BranchPtr b)
{
    assert(pos == b->get_last_xy());

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
Node::Node(BaseWeakNodePtr parent, float distance, Point position,
           std::string binaryID)
    : TopologicalNode(parent, distance, position, binaryID)
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
