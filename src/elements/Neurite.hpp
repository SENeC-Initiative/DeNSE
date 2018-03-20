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
class GrowthConeContinuousRecorder;
class Branching;
class Neuron;


//! Neurite class,
class Neurite : public std::enable_shared_from_this<Neurite>
{

    friend class GrowthConeContinuousRecorder;
    friend class Branching;
    friend class Neuron;

  public:
    Neurite() = delete;
    Neurite(std::string name, const std::string &neurite_type_,
            NeuronWeakPtr neuron);
    ~Neurite();

    // Init and finalize functions
    void init_first_node(BaseWeakNodePtr soma, Point pos,
                         std::string neurite_name, double soma_radius,
                         double neurite_diameter);
    void set_soma_angle(const double angle);
    double get_soma_angle() const;
    void finalize();

    // Growth functions
    void grow(mtPtr rnd_engine, size_t current_step, double substep);
    void delete_cone(size_t cone_n);


    // Branching functions
    void lateral_branching(TNodePtr branching_node, size_t branch_point,
                           double new_length, mtPtr rnd_engine);
    bool growth_cone_split(GCPtr branching_cone, double new_length,
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
    void update_actin_waves(mtPtr rnd_engine, double substep);
    void start_actin_wave(double actin_content);
    void add_actin(ActinPtr);

    // Get/set functions
    // void init_status(const statusMap &status);
    void set_status(const statusMap &);
    void get_status(statusMap &, const std::string& level) const;
    double get_state(const char *observable) const;
    unsigned int num_growth_cones() const;
    NodePtr get_first_node() const;
    NeuronWeakPtr get_parent_neuron() const;
    std::string get_name() const;
    double get_max_resol() const;
    void update_kernel_variables();
    size_t get_and_increment_gc_ID();
    void add_cone(GCPtr);
    std::unordered_map<size_t, GCPtr>::const_iterator gc_cbegin() const;
    std::unordered_map<size_t, GCPtr>::const_iterator gc_cend() const;

  private:
    //! Neuron parent
    NeuronWeakPtr parent_;
    Branching branching_model_;
    std::string name_;
    // keep track of how many nodes were created to set the ids
    size_t num_created_nodes_;
    size_t num_created_cones_;
    // observables for recorders
    std::vector<std::string> observables_;

    //! timestep and distributions for random number generation
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    std::poisson_distribution<> poisson_;
    std::exponential_distribution<> exponential_;

    //! growth cones
    std::string growth_cone_model_;
    std::unordered_map<size_t, GCPtr> growth_cones_;
    std::unordered_map<size_t, GCPtr> growth_cones_tmp_;
    std::vector<size_t> dead_cones_;
    std::deque<ActinPtr> actinDeck_;
    std::unordered_map<size_t, NodePtr> nodes_;
    std::vector<size_t> dead_nodes_;

    //! declare the type of neurite (dendrite or axon)
    std::string neurite_type_;
    //! branch direction parameters
    double diameter_eta_exp_;
    double diameter_variance_;
    double soma_angle_;
    double lateral_branching_angle_mean_;
    double lateral_branching_angle_std_;
    double gc_split_angle_mean_;
    double gc_split_angle_std_;

}; // Neurite


inline bool reverse_sorting(int i, int j) { return (i > j); }

} // namespace

#endif // NEURITE_H
