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


/**
 * @brief BaseNode is the root of neurite topological tree
 *
 * This class is needed to initialize the neurite tree.
 * Soma is an instance of BaseNode
 */
class BaseNode
{
    friend class Neuron;

  protected:
    BPoint position_;
    double dist_to_soma_;
    double dist_to_parent_;

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
    BPoint get_position() const;
    virtual double get_distance_to_soma() const;
    virtual double get_distance_parent() const;
    virtual stype get_node_id() const;
};


class TopologicalNode : public BaseNode
{
    friend class Neurite;

  protected:
    BaseWeakNodePtr parent_;
    int centrifugal_order_;
    bool has_child_;
    stype node_id_;

    bool dead_;
    BranchPtr branch_;
    NeuritePtr own_neurite_;
    double diameter_;

  public:
    TopologicalNode();
    TopologicalNode(const TopologicalNode &tnode);
    TopologicalNode(BaseWeakNodePtr parent, double distance_to_parent,
                    const BPoint &position, double diameter,
                    NeuritePtr neurite);

    /**
     * @brief Update the centrifugal order
     */
    void topological_advance();
    void set_first_point(const BPoint &p, double length);
    void set_diameter(double diameter);
    virtual void set_position(const BPoint &) override;
    void update_branch_and_parent(BaseNodePtr parent, BranchPtr b);

    // geometry getter functions
    inline double get_distance_to_soma() const override
    {
        return dist_to_soma_;
    }
    inline double get_distance_parent() const override
    {
        return dist_to_parent_;
    }

    // topology getter functions
    inline BaseWeakNodePtr get_parent() const { return parent_; }
    inline int get_centrifugal_order() const override
    {
        return centrifugal_order_;
    }
    inline bool has_child() const { return has_child_; }
    inline stype get_node_id() const override { return node_id_; }

    seg_range segment_range() const;

    // biology getter functions
    inline bool is_dead() const { return dead_; }
    inline BranchPtr get_branch() const { return branch_; }
    stype get_branch_size() const;
    double get_branch_length() const;
    inline virtual double get_diameter() const { return diameter_; }

    bool support_AW() const;
};


class Node : public TopologicalNode
{
    friend class Neurite;
    friend class Neuron;

  protected:
    std::vector<TNodePtr> children_;

  public:
    Node(BaseWeakNodePtr parent, double distanceToParent, const BPoint &pos,
         double diameter, NeuritePtr neurite);

    virtual void set_position(const BPoint &) override final;

    TNodePtr get_child(int) const;

}; // Node
} // namespace growth
#endif
