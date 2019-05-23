/*
 * Neurite.hpp
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

#ifndef NEURITE_H
#define NEURITE_H

// c++ includes
#include <deque>
#include <random>
#include <string>
#include <vector>

#include <boost/range/join.hpp>
#include <boost/range/iterator_range.hpp>

// elements includes
#include "Branch.hpp"
#include "Branching.hpp"
#include "growth_names.hpp"

// libgrowth includes
#include "elements_types.hpp"
#include "spatial_types.hpp"


namespace growth
{

typedef std::unordered_map<size_t, GCPtr> gc_map;
typedef boost::iterator_range< gc_map::const_iterator > simple_gc_range;
typedef boost::range::joined_range<const gc_map, const gc_map> joined_gc_range;


// Neuron forward declaration
class GrowthConeContinuousRecorder;
class Branching;
class Neuron;


typedef struct res_Neurite
{
    // res_params
    double target_cr;
    double eq_cr;
    double tau_generation;
    double tau_delivery;
    double var;
    double correlation;

    // evolution with number of growth cones
    double typical_gc_support;
    double increase_slope;

    // res_run time variable
    double available;
    double tot_demand;
    double tau; // t^-1 = t_A ^-1 t_d^-1
    double stochastic_tmp;
} res_Neurite;


//! Neurite class,
class Neurite : public std::enable_shared_from_this<Neurite>
{

    friend class GrowthConeContinuousRecorder;
    friend class Branching;
    friend class Neuron;

  public:
    Neurite() = delete;
    Neurite(const std::string &name, const std::string &neurite_type_,
            const std::string &gc_model, NeuronWeakPtr neuron);
    ~Neurite();

    // Init and finalize functions
    void init_first_node(BaseWeakNodePtr soma, const BPoint &pos,
                         const std::string &neurite_name, double soma_radius,
                         double neurite_diameter);
    void set_soma_angle(const double angle);
    double get_soma_angle() const;
    void finalize();

    // Growth functions
    void update_growth_cones(mtPtr rnd_engine, double substep);
    void update_resource(mtPtr rnd_engine, double substep);
    void grow(mtPtr rnd_engine, size_t current_step, double substep);
    void delete_cone(size_t cone_n);

    // critical resource
    double get_quotient_cr() const;
    double get_available_cr() const;

    // Branching functions
    bool lateral_branching(TNodePtr branching_node, size_t branch_point,
                           NodePtr &new_node, mtPtr rnd_engine);
    bool growth_cone_split(GCPtr branching_cone, double new_length,
                           double new_angle, double old_angle,
                           double new_diameter, double old_diameter,
                           NodePtr &new_node, GCPtr& sibling);
    GCPtr create_branching_cone(const TNodePtr branching_node, NodePtr new_node,
                                double new_length, double new_diameter,
                                const BPoint &xy, double new_cone_angle,
                                bool split);
    void update_parent_nodes(NodePtr new_node, TNodePtr oldnode);
    void update_tree_structure(TNodePtr root);
    void delete_parent_node(NodePtr parent, int child_id);

    void gc_split_angles_diameter(mtPtr rnd_engine, double &angle_1,
                                  double &angle_2, double &old_diameter,
                                  double &new_diameter);

    const Branching *get_branching_model() const;

    void get_distances(size_t node, size_t segment, double &dist_to_parent,
                       double &dist_to_soma) const;

    //@TODO
    void update_actin_waves(mtPtr rnd_engine, double substep);
    void start_actin_wave(double actin_content);
    void add_actin(ActinPtr);

    // Get/set functions
    // void init_status(const statusMap &status);
    const std::string & get_type() const;
    void set_status(const statusMap &);
    void get_status(statusMap &, const std::string &level) const;
    double get_state(const char *observable) const;
    unsigned int num_growth_cones() const;
    NodePtr get_first_node() const;
    NeuronWeakPtr get_parent_neuron() const;
    const std::string& get_name() const;
    double get_taper_rate() const;
    double get_max_resol() const;

    void update_kernel_variables();
    void add_cone(GCPtr);
    void add_node(NodePtr);

    bool walk_tree(NodeProp& np) const;
    simple_gc_range active_gc_range() const;
    joined_gc_range gc_range() const;
    std::unordered_map<size_t, NodePtr>::const_iterator nodes_cbegin() const;
    std::unordered_map<size_t, NodePtr>::const_iterator nodes_cend() const;

  private:
    //! Neuron parent
    NeuronWeakPtr parent_;
    BranchingPtr branching_model_;
    std::string name_;
    bool active_;
    double max_arbor_len_;
    double fixed_arbor_len_;
    // keep track of how many nodes were created to set the ids
    size_t num_created_nodes_;
    // observables for recorders
    std::vector<std::string> observables_;

    //! timestep and distributions for random number generation
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    std::poisson_distribution<> poisson_;
    std::exponential_distribution<> exponential_;
    std::normal_distribution<> cr_normal_;

    //! growth cones
    std::string growth_cone_model_;
    gc_map growth_cones_;
    gc_map growth_cones_tmp_;
    gc_map growth_cones_inactive_;
    gc_map growth_cones_inactive_tmp_;

    std::vector<size_t> dead_cones_;
    std::deque<ActinPtr> actinDeck_;
    std::unordered_map<size_t, NodePtr> nodes_;
    std::vector<size_t> dead_nodes_;
    size_t max_gc_num_;

    //! declare the type of neurite (dendrite or axon)
    std::string neurite_type_;
    double taper_rate_;     // diameter thinning with distance
    double min_diameter_;

    // competition
    bool use_critical_resource_;
    res_Neurite cr_neurite_;

    //! branch direction parameters
    double diameter_ratio_avg_;
    double diameter_eta_exp_;
    double diameter_ratio_std_;
    double soma_angle_;
    double lateral_branching_angle_mean_;
    double lateral_branching_angle_std_;
    double gc_split_angle_mean_;
    double gc_split_angle_std_;
    double diam_frac_lb_;

}; // Neurite


inline bool reverse_sorting(int i, int j) { return (i > j); }

} // namespace growth

#endif // NEURITE_H
