/*
 * neuron_manager.hpp
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

#ifndef NEURON_M_H
#define NEURON_M_H

#include <string>
#include <unordered_map>

#include "config.hpp"
#include "elements_types.hpp"
#include "manager_interface.hpp"

namespace growth
{

// typedefs
typedef std::unordered_map<stype, NeuronPtr> gidNeuronMap;
typedef std::unordered_map<stype, int> gidThreadMap;
typedef std::unordered_map<std::string, GCPtr> modelMap;
typedef std::vector<std::vector<NeuronPtr>> threadNeurons;


// forward declarations
class Neuron;
class GrowthCone;


class NeuronManager : public ManagerInterface
{
  public:
    NeuronManager();

    virtual void initialize();
    virtual void finalize();

    /**
     * Create neurons.
     */
    stype create_neurons(const std::vector<statusMap> &neuron_params,
                         const std::vector<statusMap> &axon_params,
                         const std::vector<statusMap> &dendrites_params);

    void delete_neurons(const std::vector<stype> &gids);

    NeuronPtr get_neuron(stype gid);
    void get_all_neurons(std::vector<NeuronPtr> &);
    std::vector<stype> get_gids() const;
    gidNeuronMap get_local_neurons(int local_thread_id);
    int get_neuron_thread(stype gid) const;

    void init_neurons_on_thread(unsigned int num_local_threads);
    void update_kernel_variables();

    void get_defaults(statusMap &status, const std::string &object,
                      bool detailed) const;
    const statusMap get_neuron_status(stype gid) const;
    const statusMap get_neurite_status(stype gid, const std::string &type,
                                       const std::string &level) const;

    bool is_neuron(stype gid) const;

    stype num_neurons() const;

    gidNeuronMap::const_iterator iter_neurons();

    void set_max_resol(stype neuron, double max_resol);
    double get_max_resol() const;

    //! here we call the models file and register each model in the map
    void register_model(std::string, GCPtr);
    modelMap model_map_;

  private:
    stype num_created_neurons_;
    NeuronPtr model_neuron_;          // unused model neuron for get_defaults
    gidNeuronMap neurons_;            // get neuron from gid
    threadNeurons neurons_on_thread_; // group neurons by thread
    gidThreadMap thread_of_neuron_;   // get thread from gid
    std::vector<std::unordered_map<stype, double>>
        max_resolutions_; // max allowed resol
};


void _fill_skel(const SkelNeurite &source_container,
                SkelNeurite &target_container, bool add_nan);
} // namespace growth

#endif // NEURON_M_H
