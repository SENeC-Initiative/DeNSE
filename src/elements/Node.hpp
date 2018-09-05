#ifndef NODE_H
#define NODE_H

#include <random>

#include "elements_types.hpp"
#include "spatial_types.hpp"

namespace growth
{

class Neurite;
class Neuron;


typedef struct NodeTopology
{
    BaseWeakNodePtr parent;
    int centrifugal_order;
    bool has_child;
    size_t nodeID;
    std::string binaryID;
    NodeTopology()
        : parent(std::make_shared<BaseNode>())
        , centrifugal_order(-1)
        , has_child(false)
        , nodeID(0)
        , binaryID("default")
    {
    }
    NodeTopology(BaseWeakNodePtr parent, int centrifugal_order, bool has_child,
                 int nodeID, std::string binaryID)
        : parent(parent)
        , centrifugal_order(centrifugal_order)
        , has_child(has_child)
        , nodeID(nodeID)
        , binaryID(binaryID)
    {
    }
} NodeTopology;


typedef struct NodeGeometry
{
    Point position;
    double dis_to_soma;
    double dis_to_parent;
    NodeGeometry()
        : position(Point())
        , dis_to_soma(-1)
        , dis_to_parent(-1)
    {
    }
    NodeGeometry(Point position, double dis_to_soma, double dis_to_parent)
        : position(position)
        , dis_to_soma(dis_to_soma)
        , dis_to_parent(dis_to_parent)
    {
    }
} NodeGeometry;


typedef struct NodeBiology
{
    bool dead;
    BranchPtr branch;
    NeuritePtr own_neurite;
    double diameter;

    NodeBiology()
        : dead(true)
        , branch(std::make_shared<Branch>())
        , own_neurite(nullptr)
        , diameter(0)
    {
    }

    NodeBiology(bool dead, BranchPtr branch, NeuritePtr own_neurite,
                double diameter)
        : dead(dead)
        , branch(branch)
        , own_neurite(own_neurite)
        , diameter(diameter)
    {
    }
} NodeBiology;


/**
 * @brief BaseNode is the root of neurite topological tree
 *
 * This class is needed to initialize the neurite tree.
 * Soma is an istance of BaseNode
 */
class BaseNode
{
    friend class Neuron;

  protected:
    NodeGeometry geometry_;

  private:
    std::string somaID_;

  public:
    /**
     * @brief Set the position of the node
     *
     * @param Point xy: x and y of the node.
     */
    virtual void set_position(const Point &);

    //! getters function
    virtual int get_centrifugal_order() const;
    virtual Point get_position() const;
    virtual double get_distance_to_soma() const;
    virtual double get_distance_parent() const;
    virtual std::string get_treeID() const;
    virtual size_t get_nodeID() const;
};


class TopologicalNode : public BaseNode
{
    friend class Neurite;

  protected:
    NodeTopology topology_;
    NodeGeometry geometry_;
    NodeBiology biology_;

  public:
    TopologicalNode();
    TopologicalNode(const TopologicalNode &tnode);
    TopologicalNode(BaseWeakNodePtr parent, float distanceToParent,
                    const Point &position, const std::string &binaryID);

    // typedef of method pointer, used in Branching.cpp
    /**
     * @brief Update the binaryId and the centrifugal order
     */
    void topological_advance();
    void set_first_point(const Point p, double length);
    virtual void set_diameter(double diameter);
    void set_position(const Point &) override;
    void set_position(const Point &pos, double dist_to_soma, BranchPtr b);

    // geometry getter functions
    inline Point get_position() const override { return geometry_.position; }
    inline double get_distance_to_soma() const override
    {
        return geometry_.dis_to_soma;
    }
    inline double get_distance_parent() const override
    {
        return geometry_.dis_to_parent;
    }

    // topology getter functions
    inline BaseWeakNodePtr get_parent() const { return topology_.parent; }
    inline int get_centrifugal_order() const override
    {
        return topology_.centrifugal_order;
    }
    inline bool has_child() const { return topology_.has_child; }
    inline size_t get_nodeID() const override { return topology_.nodeID; }
    inline std::string get_treeID() const override
    {
        return topology_.binaryID;
    }

    // biology getter functions
    inline bool is_dead() const { return biology_.dead; }
    inline BranchPtr get_branch() const { return biology_.branch; }
    size_t get_branch_size() const;
    inline virtual double get_diameter() const { return biology_.diameter; }

    bool support_AW() const;
};


class Node : public TopologicalNode
{
    friend class Neurite;
    friend class Neuron;

  protected:
    std::vector<TNodePtr> children_;

  public:
    typedef size_t (Node::*Get_method)(void) const;
    Node(BaseWeakNodePtr parent, float distanceToParent, Point position,
         std::string binaryID);

    Node(const Node &);

    Node(const TopologicalNode &);

    TNodePtr get_child(int) const;

}; // Node
} // namespace growth
#endif
