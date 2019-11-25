/*
 * module.hpp
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

#ifndef MODULE_H
#define MODULE_H


// C++ include
#include <array>
#include <string>
#include <vector>

// include GEOS C-API
#include <geos_c.h>

// kernel include
#include "recorders.hpp"

// libgrowth include
#include "config.hpp"
#include "elements_types.hpp"
#include "growth_time.hpp"

// spatial include
#include "space_manager.hpp"
#include "spatial_types.hpp"


namespace growth
{

/* Init and finalize */

void finalize_growth_();


void init_growth_(int *argc, char **argv[]);


void reset_kernel_();


/* Creation and deletion */

stype create_objects_(const std::string &object_name,
                      const std::vector<statusMap> &obj_params);


stype create_neurons_(const std::vector<statusMap> &neuron_params,
                      const std::vector<statusMap> &axon_params,
                      const std::vector<statusMap> &dendrites_params);


void create_neurites_(const std::vector<stype> &neurons, stype num_neurites,
                      const std::vector<statusMap> &params,
                      const std::vector<std::string> &neurite_types,
                      const std::vector<double> &angles,
                      const std::vector<std::string> &names);


void delete_neurons_(const std::vector<stype> &gids);


void delete_neurites_(const std::vector<stype> &gids,
                      const std::vector<std::string> &names);


/* Setters */

void set_kernel_status_(const statusMap &status_dict,
                        std::string simulation_ID);


void set_environment_(
    GEOSGeometry *environment, const std::vector<GEOSGeometry *> &areas,
    std::vector<double> heights, const std::vector<std::string> &names,
    const std::vector<std::unordered_map<std::string, double>> &properties);


void set_status_(stype gid, statusMap status, statusMap axon_status,
                 statusMap dendrites_status);


void set_neurite_status_(stype neuron, std::string neurite, statusMap status);


/* Simulation */

void simulate_(const Time &simtime);


void test_random_generator_(Random_vecs &values, stype size);


/* Getters functions */

void get_environment_(
    GEOSGeometry *&environment, std::vector<GEOSGeometry *> &areas,
    std::vector<double> &heights, std::vector<std::string> &names,
    std::vector<std::unordered_map<std::string, double>> &properties);


const Time get_current_time_();


std::string get_simulation_id_();


void get_skeleton_(SkelNeurite &axon, SkelNeurite &dendrites,
                   SkelNeurite &nodes, SkelNeurite &growth_cones,
                   SkelSomas &somas, std::vector<stype> gids,
                   unsigned int resolution);


void get_swc_(std::string output_file, std::vector<stype> gids,
              unsigned int resolution);


statusMap get_kernel_status_();


Property get_kernel_status_(const std::string &key);


stype get_num_objects_();


stype get_num_created_objects_();


double get_state_(stype gid, const std::string &level,
                  const std::string &variable);


statusMap get_status_(stype gid);


statusMap get_neurite_status_(stype gid, const std::string &neurite_type,
                              const std::string &level);


bool is_neuron_(stype gid);

bool is_neurite_(stype gid, const std::string &neurite);


// neuron- and structure-related

std::vector<stype> get_neurons_();


std::vector<std::string> get_neurites_(stype gid);


bool neuron_has_axon_(stype gid);


void get_branches_data_(stype neuron, const std::string &neurite,
                        std::vector<std::vector<std::vector<double>>> &points,
                        std::vector<double> &diameters,
                        std::vector<int> &parents, std::vector<stype> &nodes,
                        stype start_point);


void get_geom_skeleton_(std::vector<stype> gids,
                        std::vector<GEOSGeometry *> &axons,
                        std::vector<GEOSGeometry *> &dendrites,
                        std::vector<stype> &dendrite_gids,
                        std::vector<std::vector<double>> &somas);


void generate_synapses_(
    bool crossings_only, double density, bool only_new_syn,
    bool autapse_allowed, const std::set<stype> &presyn_pop,
    const std::set<stype> &postsyn_pop, std::vector<stype> &presyn_neurons,
    std::vector<stype> &postsyn_neurons,
    std::vector<std::string> &presyn_neurites,
    std::vector<std::string> &postsyn_neurites,
    std::vector<stype> &presyn_nodes, std::vector<stype> &postsyn_nodes,
    std::vector<stype> &presyn_segments, std::vector<stype> &postsyn_segments,
    std::vector<double> &pre_syn_x, std::vector<double> &pre_syn_y,
    std::vector<double> &post_syn_x, std::vector<double> &post_syn_y);


void get_distances_(stype gid, const std::string &neurite_name, stype node,
                    stype segment, double &dist_to_parent,
                    double &dist_to_soma);


// parameters and recordings

void get_defaults_(const std::string &object_name,
                   const std::string &object_type, const std::string &gc_model,
                   bool detailed, statusMap &status);


void get_models_(std::unordered_map<std::string, std::string> &models,
                 bool abbrev);


std::vector<std::string> get_elongation_types_();


std::vector<std::string> get_steering_methods_();


std::vector<std::string> get_direction_selection_methods_();


void get_recorder_type_(stype gid, std::string &level, std::string &event_type);


bool get_next_recording_(stype gid, std::vector<Property> &ids,
                         std::vector<double> &values);


bool get_next_time_(stype gid, std::vector<Property> &ids,
                    std::vector<double> &values, const std::string &time_units);


std::string get_default_model_();


/* tools */

std::string object_type_(stype gid);


void _fill_skel(const SkelNeurite &source_container,
                SkelNeurite &target_container, bool add_nan);


bool walk_neurite_tree_(stype neuron, std::string neurite, NodeProp &np);


void is_valid_timestep_(double timestep);


void get_backtrace_(std::string &msg, int depth);

} // namespace growth

#endif // MODULE_H
