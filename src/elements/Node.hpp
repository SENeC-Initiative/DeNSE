/*
 * Node.hpp
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

    NodeTopology()
        : parent(std::make_shared<BaseNode>())
        , centrifugal_order(-1)
        , has_child(false)
        , nodeID(0)
    {
    }
    NodeTopology(BaseWeakNodePtr parent, int centrifugal_order, bool has_child,
                 int nodeID)
        : parent(parent)
        , centrifugal_order(centrifugal_order)
        , has_child(has_child)
        , nodeID(nodeID)
    {
    }
} NodeTopology;


typedef struct NodeGeometry
{
    BPoint position;
    double dis_to_soma;
    double dis_to_parent;
    NodeGeometry()
      : position(BPoint())
      , dis_to_soma(-1)
      , dis_to_parent(-1)
    {
    }
    NodeGeometry(const BPoint &position, double dis_to_soma, double dis_to_parent)
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
        : dead(false)
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

  public:
    BaseNode();
    BaseNode(const BPoint &position, double distance_to_parent,
             double distance_to_soma);
    BaseNode(const BaseNode &copy);

    /**
     * @brief Set the position of the node
     *
     * @param Point xy: x and y of the node.
     */
    virtual void set_position(const BPoint &);

    //! getters function
    virtual int get_centrifugal_order() const;
    virtual BPoint get_position() const;
    virtual double get_distance_to_soma() const;
    virtual double get_distance_parent() const;
    virtual size_t get_node_id() const;
};


class TopologicalNode : public BaseNode
{
    friend class Neurite;

  protected:
    NodeTopology topology_;
    NodeBiology biology_;

  public:
    TopologicalNode();
    TopologicalNode(const TopologicalNode &tnode);
    TopologicalNode(BaseWeakNodePtr parent, float distanceToParent,
                    const BPoint &position);

    /**
     * @brief Update the centrifugal order
     */
    void topological_advance();
    void set_first_point(const BPoint &p, double length);
    virtual void set_diameter(double diameter);
    void set_position(const BPoint &) override;
    void set_position(const BPoint &pos, double dist_to_soma, BranchPtr b);

    // geometry getter functions
    inline BPoint get_position() const override { return geometry_.position; }
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
    inline size_t get_node_id() const override { return topology_.nodeID; }

    seg_range segment_range() const;

    // biology getter functions
    inline bool is_dead() const { return biology_.dead; }
    inline BranchPtr get_branch() const { return biology_.branch; }
    size_t get_branch_size() const;
    double get_branch_length() const;
    inline virtual double get_diameter() const { return biology_.diameter; }

    bool support_AW() const;
};


class Node : public TopologicalNode
{
    friend class Neurite;
    friend class Neuron;

  protected:
    std::vector<TNodePtr> children_;
    double diameter_;  // diameter at the node's position

  public:
    Node(BaseWeakNodePtr parent, float distanceToParent, const BPoint &pos);

    Node(const Node &copy);

    Node(const TopologicalNode &copy);

    TNodePtr get_child(int) const;

}; // Node
} // namespace growth
#endif
