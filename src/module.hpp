#ifndef MODULE_H
#define MODULE_H


// C++ include
#include <array>
#include <string>
#include <vector>

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


/* Creation */

size_t create_objects_(const std::string &object_name,
                       const std::vector<statusMap> &obj_params);


size_t create_neurons_(const std::vector<statusMap> &neuron_params,
                       const std::vector<statusMap> &axon_params,
                       const std::vector<statusMap> &dendrites_params);


/* Setters */

void set_kernel_status_(const statusMap &status_dict, std::string simulation_ID);


void set_environment_(
    GEOSGeometry *environment, const std::vector<GEOSGeometry *> &areas,
    std::vector<double> heights, const std::vector<std::string> &names,
    const std::vector<std::unordered_map<std::string, double>> &properties);


void set_status_(size_t gid, statusMap status, statusMap axon_status,
                statusMap dendrites_status);


/* Simulation */

void simulate_(const Time &simtime);


void test_random_generator_(Random_vecs &values, size_t size);


/* Getters functions */

void get_environment_(
    GEOSGeometry * &environment, std::vector<GEOSGeometry *> &areas,
    std::vector<double> &heights, std::vector<std::string> &names,
    std::vector<std::unordered_map<std::string, double>> &properties);


const Time get_current_time_();


std::string get_simulation_id_();


void get_skeleton_(SkelNeurite &axon, SkelNeurite &dendrites,
                   SkelNeurite &nodes, SkelNeurite &growth_cones,
                   SkelSomas &somas,  std::vector<size_t> gids,
                   unsigned int resolution);


void get_swc_(std::string output_file, std::vector<size_t> gids,
              unsigned int resolution);


statusMap get_kernel_status_();


Property get_kernel_status_(const std::string &key);


size_t get_num_objects_();


double get_state_(size_t gid, const std::string& level,
                  const std::string& variable);


statusMap get_status_(size_t gid);


statusMap get_neurite_status_(size_t gid, const std::string &neurite_type,
                              const std::string &level);


bool is_neuron_(size_t gid);

bool is_neurite_(size_t gid, const std::string& neurite);


// neuron- and structure-related

std::vector<size_t> get_neurons_();


std::vector<std::string> get_neurites_(size_t gid);


void get_branches_data_(size_t neuron, const std::string &neurite,
                        std::vector<std::vector<std::vector<double>>> &points,
                        std::vector<double> &diameters,
                        std::vector<int> &parents, std::vector<size_t> &nodes,
                        size_t start_point);


void get_geom_skeleton_(std::vector<size_t> gids,
                        std::vector<GEOSGeometry*>& axons,
                        std::vector<GEOSGeometry*>& dendrites,
                        std::vector< std::vector<double> >& somas);


// parameters and recordings

void get_defaults_(const std::string &object_name,
                   const std::string &object_type, const std::string &gc_model,
                   bool detailed, statusMap &status);


void get_models_(std::unordered_map<std::string, std::string> &models,
                 bool abbrev);


std::vector<std::string> get_elongation_types_();


std::vector<std::string> get_steering_methods_();


std::vector<std::string> get_direction_selection_methods_();


void get_recorder_type_(size_t gid, std::string &level, std::string &event_type);


bool get_next_recording_(size_t gid, std::vector<Property> &ids,
                         std::vector<double> &values);


bool get_next_time_(size_t gid, std::vector<Property> &ids,
                    std::vector<double> &values, const std::string &time_units);


std::string get_default_model_();


/* tools */

std::string object_type_(size_t gid);


void _fill_skel(const SkelNeurite &source_container,
                SkelNeurite &target_container, bool add_nan);


bool walk_neurite_tree_(size_t neuron, std::string neurite, NodeProp& np);


void is_valid_timestep_(double timestep);


void get_backtrace_(std::string &msg, int depth);

} // namespace growth

#endif // MODULE_H
