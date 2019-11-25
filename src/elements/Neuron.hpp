/*
 * Neuron.hpp
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


#ifndef NEURON_H
#define NEURON_H

#define _USE_MATH_DEFINES

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
class NeuronManager;
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
    friend class NeuronManager;
    friend class Skeleton;
    friend class Swc;

  public:
    Neuron() = delete;
    Neuron(stype gid);
    ~Neuron();

    // Init functions
    void init_status(const statusMap &status, const statusMap &astatus,
                     const statusMap &dstatus, mtPtr rnd_engine);
    void initialize_next_event(mtPtr rnd_engine);
    void finalize();

    // Growth functions
    void grow(mtPtr rnd_engine, stype current_step, double substep);
    bool branch(mtPtr rnd_engine, const Event &ev);
    void next_actin_event(mtPtr rnd_engine);

    // New neurite function
    std::string new_neurite(const std::string &name,
                            const std::string &neurite_type,
                            const GCPtr gc_model, mtPtr rnd_engine);

    void delete_neurites(const std::vector<std::string> &names);

    // Getter/setter functions
    BaseNodePtr get_soma() const;
    BPoint get_position() const;
    stype get_gid() const;
    std::string get_gc_model() const;
    NeuriteWeakPtr get_neurite(const std::string &name) const;
    double get_state(const std::string &observable) const;
    void get_status(statusMap &status) const;
    int get_num_neurites() const;
    double get_soma_radius() const;
    bool is_neurite(const std::string &neurite);

    void set_status(const statusMap &status);
    void set_neurite_status(const std::string &neurite,
                            const statusMap &status);
    void get_neurite_status(statusMap &status, std::string neurite_type,
                            const std::string &level);
    void update_kernel_variables();

    bool has_axon() const;

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
    stype gid_;
    std::string description_;
    //! Container for the ``NeuritePtr`` objects
    NeuriteMap neurites_;
    BaseNodePtr soma_;
    bool has_axon_;
    //! Center of mass of the neuron's soma
    NeuronDetails details;
    // obserables for recorders
    std::vector<std::string> observables_;
    // Actin waves
    bool use_actin_waves_;
    stype aw_generation_step_;
    double actin_content_;
    stype next_actin_event_;
    double axon_angle_;
    std::unordered_map<std::string, double> neurite_angles_;
    double axon_polarization_weight_;
    double polarization_strength_;
    bool axon_angle_set_;
    bool random_rotation_angles_;
    double rnd_angle_;

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

} // namespace growth

#endif // NEURON_H
