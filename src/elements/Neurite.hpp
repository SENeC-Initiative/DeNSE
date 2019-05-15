#ifndef NEURITE_H
#define NEURITE_H

// c++ includes
#include <deque>
#include <random>
#include <string>
#include <vector>

// elements includes
#include "Branch.hpp"
#include "Branching.hpp"
#include "growth_names.hpp"

// libgrowth includes
#include "elements_types.hpp"
#include "spatial_types.hpp"


namespace growth
{


// Neuron forward declaration
class Neuron;
class Branching;

//! Neurite class,
class Neurite : public std::enable_shared_from_this<Neurite>
{

    friend class Neuron;
    friend class Branching;

  private:
    //! Neuron parent
    NeuronWeakPtr parent_;
    Branching branching_model_;
    // keep track of how many nodes were created to set the ids
    size_t num_created_nodes_;

    //! timestep and distributions for random number generation
    double timestep_;
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    std::poisson_distribution<> poisson_;
    std::exponential_distribution<> exponential_;

    std::vector<GCPtr> growth_cones_tmp_;
    //! string list of currently GrowthCones
    std::vector<GCPtr> growth_cones_;
    std::vector<size_t> dead_cones_;
    std::deque<ActinPtr> actinDeck_;
    std::unordered_map<size_t, NodePtr> nodes_;
    std::vector<size_t> dead_nodes_;

    //! declare the type of neurite (dendrite or axon)
    std::string neurite_type_;
    //! branch direction parameters
    double diameter_eta_exp_;
    double diameter_variance_;
    double lateral_branching_angle_mean_;
    double lateral_branching_angle_std_;
    double gc_split_angle_mean_;
    double gc_split_angle_std_;

  public:
    size_t num_created_cones_;
    Neurite();
    Neurite(const std::string &neurite_type_, NeuronWeakPtr neuron);
    ~Neurite();

    // Init functions
    void init_first_node(BaseWeakNodePtr soma, Point pos,
                         std::string neurite_name, double soma_radius,
                         double neurite_diameter);

    // Growth functions
    void grow(mtPtr rnd_engine);
    void update_growth_cones(mtPtr rnd_engine);
    void delete_cone(size_t cone_n);

    // Branching functions
    void lateral_branching(TNodePtr branching_node, size_t branch_point,
                           double new_length, mtPtr rnd_engine);
    void growth_cone_split(GCPtr branching_cone, double new_length,
                           double new_angle, double old_angle,
                           double new_diameter, double old_diameter);
    GCPtr create_branching_cone(const TNodePtr branching_node, NodePtr new_node,
                                double new_length, Point xy,
                                double new_cone_angle);
    void update_parent_nodes(NodePtr new_node, TNodePtr oldnode);
    void update_tree_structure(TNodePtr root);
    void delete_parent_node(NodePtr parent, int child_id);

    void gc_split_angles_diameter(mtPtr rnd_engine, double &angle_1,
                                  double &angle_2, double &old_diameter,
                                  double &new_diameter);

    const Branching *get_branching_model() const;

    //@TODO
    void update_actin_waves(mtPtr rnd_engine);
    void start_actin_wave(double actin_content);
    void add_actin(ActinPtr);

    // Get/set functions
    // void init_status(const statusMap &status);
    void set_status(const statusMap &);
    void get_status(statusMap &) const;
    int num_growth_cones() const;
    NodePtr get_first_node() const;
    void update_kernel_variables();
    size_t get_and_increment_gc_ID();
    void add_cone(GCPtr);
    std::vector<GCPtr>::const_iterator gc_cbegin() const;
    std::vector<GCPtr>::const_iterator gc_cend() const;

}; // Neurite

inline bool reverse_sorting(int i, int j) { return (i > j); }
} // namespace

#endif // NEURITE_H
