#ifndef MODULE_H
#define MODULE_H

#include <array>
#include <string>
#include <vector>

#include "config.hpp"
#include "elements_types.hpp"
#include "growth_time.hpp"
#include "space_manager.hpp"
#include "spatial_types.hpp"

namespace growth
{

size_t create_objects(const std::string &object_name,
                      const std::vector<statusMap> &obj_params);

/**
 * @brief Create neurons and set status
 *
 * The number of created neurons depend on the size of the vector
 *
 * @param neuron_params set the parameters for all neurites
 * @param axon_params overwrite the params for the axon, empty vector is default
 * @param dendrites_params overwrite the params for the densrite, empty vecor is
 * default
 *
 * @return
 */
size_t create_neurons(const std::vector<statusMap> &neuron_params,
                      const std::vector<statusMap> &axon_params,
                      const std::vector<statusMap> &dendrites_params);

void finalize_growth();


void init_growth(int *argc, char **argv[]);


void reset_kernel();

void set_kernel_status(const statusMap &status_dict, std::string simulation_ID);

void set_environment(GEOSGeom environment);


void set_status(size_t gid, statusMap status, statusMap axon_status,
                statusMap dendrites_status);

void simulate(const Time &simtime);

void test_random_generator(Random_vecs &values, size_t size);

/* Getters functions */

void get_environment(GEOSGeom &environment);

const Time get_current_time();

std::string get_simulation_ID();

void get_skeleton(SkelNeurite &axon, SkelNeurite &dendrites, SkelNeurite &nodes,
                  SkelNeurite &growth_cones, SkelSomas &somas,
                  std::vector<size_t> gids);


void get_swc(std::string output_file, std::vector<size_t> gids,
             unsigned int resolution);

statusMap get_kernel_status();

Property get_kernel_status(const std::string &key);

size_t get_num_objects();

statusMap get_status(size_t gid);

statusMap get_neurite_status(size_t gid, const std::string &neurite_type);

std::vector<size_t> get_neurons();

void get_defaults(const std::string &object_name,
                  const std::string &object_type, statusMap &status);

void get_models(std::vector<std::string> &models,
                const std::string &object_type);

// tool

std::string object_type(size_t gid);

void _fill_skel(const SkelNeurite &source_container,
                SkelNeurite &target_container, bool add_nan);
}

#endif // MODULE_H
