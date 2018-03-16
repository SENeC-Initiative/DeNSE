/**
 * \file Neuron.hpp
 * This file is part of the Growth project.
 */

#ifndef NEURON_H
#define NEURON_H

// c++ includes
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

// elements includes
#include "Branch.hpp"
#include "Neurite.hpp"

// libgrowth includes
#include "elements_types.hpp"
#include "spatial_types.hpp"


namespace growth
{

typedef struct NeuronDetails
{
    double soma_radius;
    double dendrite_diameter;
    double axon_diameter;
    NeuronDetails()
        : soma_radius(SOMA_RADIUS)
        , dendrite_diameter(DENDRITE_DIAMETER)
        , axon_diameter(AXON_DIAMETER)
    {
    }
} NeuronDetails;


/*
 *  Friend classes forward declaration
 */

class GrowthConeContinuousRecorder;
class NeuriteContinuousRecorder;
class NeuriteDiscreteRecorder;
class Skeleton;
class Swc;


/**
 * Implementation of the main container ``Neuron`` for the creation of
 * neuronal networks.
 */
class Neuron : public std::enable_shared_from_this<Neuron>
{

    friend class GrowthConeContinuousRecorder;
    friend class NeuriteContinuousRecorder;
    friend class NeuriteDiscreteRecorder;
    friend class Skeleton;
    friend class Swc;

  public:
    Neuron() = delete;
    Neuron(size_t gid);

    // Init functions
    void init_status(const statusMap &status, const statusMap &astatus,
                     const statusMap &dstatus, mtPtr rnd_engine);
    void initialize_next_event(mtPtr rnd_engine, double new_resolution,
                               size_t previous_step);
    void finalize();

    // Growth functions
    void grow(mtPtr rnd_engine, size_t current_step, double substep);
    bool branch(mtPtr rnd_engine, const Event &ev);
    void next_actin_event(mtPtr rnd_engine);

    // New neurite function
    std::string new_neurite(const std::string &name,
                            const std::string &neurite_type,
                            const GCPtr gc_model, mtPtr rnd_engine);


    // Getter/setter functions
    BaseNodePtr get_soma() const;
    Point get_position() const;
    size_t get_gid() const;
    std::string get_gc_model() const;
    double get_state(const char *observable) const;
    void get_status(statusMap &status) const;
    int get_num_neurites() const;
    double get_soma_radius() const;

    void set_status(const statusMap &status);
    void set_neurite_status(const std::string &neurite_type,
                            const statusMap &status);
    void get_neurite_status(statusMap &status, std::string neurite_type,
                            const std::string& level);
    void update_kernel_variables();

    // constant iterator to neurites map
    inline std::unordered_map<std::string, NeuritePtr>::const_iterator
    neurite_cbegin() const
    {
        return neurites_.cbegin();
    }

    inline std::unordered_map<std::string, NeuritePtr>::const_iterator
    neurite_cend() const
    {
        return neurites_.cend();
    }

  private:
    size_t gid_;
    std::string description_;
    //! Container for the ``NeuritePtr`` objects
    NeuriteMap neurites_;
    BaseNodePtr soma_;
    //! Center of mass of the neuron's soma
    NeuronDetails details;
    // obserables for recorders
    std::vector<std::string> observables_;
    // Actin waves
    bool use_actin_waves_;
    size_t aw_generation_step_;
    double actin_content_;
    size_t next_actin_event_;
    double axon_angle_;
    bool axon_angle_set_;

    //! Simulator space variables, passed with statusMap
    //!\param actin_waves_triggered is bool value to activate or not the actin
    //! wave
    //! defalut is ''False''
    //! \param growth_cone_model is the model that will be used for the
    //! ''Dendrite'
    //  objects of this ''Neuron''
    std::string growth_cone_model_;
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
};

}

#endif // NEURON_H
