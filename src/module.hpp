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

void finalize_growth();


void init_growth(int *argc, char **argv[]);


void reset_kernel();


/* Creation */

size_t create_objects(const std::string &object_name,
                      const std::vector<statusMap> &obj_params);


size_t create_neurons(const std::vector<statusMap> &neuron_params,
                      const std::vector<statusMap> &axon_params,
                      const std::vector<statusMap> &dendrites_params);


/* Setters */

void set_kernel_status(const statusMap &status_dict, std::string simulation_ID);


void set_environment(
    GEOSGeom environment, const std::vector<GEOSGeom> &areas,
    std::vector<double> heights, const std::vector<std::string> &names,
    const std::vector<std::unordered_map<std::string, double>> &properties);


void set_status(size_t gid, statusMap status, statusMap axon_status,
                statusMap dendrites_status);


/* Simulation */

void simulate(const Time &simtime);


void test_random_generator(Random_vecs &values, size_t size);


/* Getters functions */

void get_environment(
    GEOSGeom &environment, std::vector<GEOSGeom> &areas,
    std::vector<double> &heights, std::vector<std::string> &names,
    std::vector<std::unordered_map<std::string, double>> &properties);


const Time get_current_time();


std::string get_simulation_ID();


void get_skeleton(SkelNeurite &axon, SkelNeurite &dendrites, SkelNeurite &nodes,
                  SkelNeurite &growth_cones, SkelSomas &somas,
                  std::vector<size_t> gids, unsigned int resolution);


void get_swc(std::string output_file, std::vector<size_t> gids,
             unsigned int resolution);


statusMap get_kernel_status();


Property get_kernel_status(const std::string &key);


size_t get_num_objects();


double get_state(size_t gid, const std::string& level,
                 const std::string& variable);


statusMap get_status(size_t gid);


statusMap get_neurite_status(size_t gid, const std::string &neurite_type,
                             const std::string &level);


bool is_neurite(size_t gid, const std::string& neurite);


// neuron- and structure-related

std::vector<size_t> get_neurons();


std::vector<std::string> get_neurites(size_t gid);


void get_branches_data(size_t neuron, const std::string &neurite,
                       std::vector<std::vector<std::vector<double>>> &points,
                       std::vector<double> &diameters,
                       std::vector<int> &parents, std::vector<size_t> &nodes,
                       size_t start_point);


// parameters and recordings

void get_defaults(const std::string &object_name,
                  const std::string &object_type, statusMap &status);


void get_models(std::vector<std::string> &models,
                const std::string &object_type);


void get_recorder_type(size_t gid, std::string &level, std::string &event_type);


bool get_next_recording(size_t gid, std::vector<Property> &ids,
                        std::vector<double> &values);


bool get_next_time(size_t gid, std::vector<Property> &ids,
                   std::vector<double> &values, const std::string &time_units);


/* tools */

std::string object_type(size_t gid);


void _fill_skel(const SkelNeurite &source_container,
                SkelNeurite &target_container, bool add_nan);


bool walk_neurite_tree(size_t neuron, std::string neurite, NodeProp& np);


void is_valid_timestep(double timestep);

} // namespace growth

#endif // MODULE_H
