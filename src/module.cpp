#include "module.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

// IO includes
#include "GrowthCone.hpp"
#include "Neuron.hpp"
#include "Skeleton.hpp"
#include "Swc.hpp"

// kernel include
#include "kernel_manager.hpp"
#include "neuron_manager.hpp"
#include "space_manager.hpp"

// models include
#include "models_manager.hpp"


namespace growth
{

size_t create_objects(const std::string &object_name,
                      const std::vector<statusMap> &obj_params)
{
    throw InvalidParameter("KernelManager::create_objects",
                           "Creating something other than a neuron (`" +
                               object_name + "`) is not implemented yet.");

    return 0;
}

size_t create_neurons(const std::vector<statusMap> &neuron_params,
                      const std::vector<statusMap> &axon_params,
                      const std::vector<statusMap> &dendrites_params)
{
    return kernel().neuron_manager.create_neurons(neuron_params, axon_params,
                                                  dendrites_params);
}

/*
 * Getter functions
 */

const Time get_current_time() { return kernel().simulation_manager.get_time(); }


statusMap get_kernel_status() { return kernel().get_status(); }


statusMap get_status(size_t gid)
{
    if (kernel().neuron_manager.is_neuron(gid))
        return kernel().neuron_manager.get_neuron_status(gid);
    else
        throw std::runtime_error("Only neurons are supprted so far.");
}


size_t get_num_objects() { return kernel().get_num_objects(); }


statusMap get_neurite_status(size_t gid, const std::string &neurite_type)
{
    return kernel().neuron_manager.get_neurite_status(gid, neurite_type);
}

void get_defaults(const std::string &object_name,
                  const std::string &object_type, statusMap &status)
{
    GCPtr gc = nullptr;

    if (object_type == "growth_cones")
    {
        gc = kernel().neuron_manager.get_model(object_name);
        gc->get_status(status);
    }
}


void get_models(std::vector<std::string> &models,
                const std::string &object_type)
{
    if (object_type == "all" || object_type == "growth_cones")
    {
        kernel().neuron_manager.get_models(models);
    }
}


/*
 * Init functions
 */

void init_growth(int *argc, char **argv[])
{
    KernelManager::create_kernel_manager();
    kernel().parallelism_manager.init_mpi(argc, argv);
    kernel().initialize();
    init_models();
}

void finalize_growth()
{
    kernel().parallelism_manager.mpi_finalize();
    KernelManager::destroy_kernel_manager();
}

std::string object_type(size_t gid)
{
    if (kernel().neuron_manager.is_neuron(gid))
        return "neuron";
    else
        throw std::runtime_error("Object does not exist.");
}

void reset_kernel() { kernel().reset(); }

void set_kernel_status(const statusMap &status_dict, std::string simulation_ID)
{
    kernel().set_simulation_ID(simulation_ID);
    kernel().set_status(status_dict);
}

std::string get_simulation_ID() { return kernel().get_simulation_ID(); }

void set_environment(GEOSGeom environment)
{
    kernel().space_manager.set_environment(environment);
}

void get_environment(GEOSGeom &environment)
{
    kernel().space_manager.get_environment(environment);
}

void set_status(size_t gid, statusMap neuron_status, statusMap axon_status,
                statusMap dendrites_status)
{
    auto neuron = kernel().neuron_manager.get_neuron(gid);
    neuron->set_status(neuron_status);

    auto local_params = neuron_status;
    for (auto &param : axon_status)
    {
        local_params[param.first] = param.second;
    }
    neuron->set_neurite_status("axon", local_params);

    local_params = neuron_status;
    for (auto &param : dendrites_status)
    {
        local_params[param.first] = param.second;
    }
    neuron->set_neurite_status("dendrites", local_params);
}

void test_random_generator(Random_vecs &values, size_t size)
{
    kernel().simulation_manager.test_random_generator(values, size);
    printf("%lu number generated from rng\n", values[0].size());
}

void get_skeleton(SkelNeurite &axon, SkelNeurite &dendrites, SkelNeurite &nodes,
                  SkelNeurite &growth_cones, SkelSomas &somas,
                  std::vector<size_t> gids)


{
#ifndef NDEBUG
    printf("Getting skeleton\n");
#endif
    std::vector<NeuronPtr> neurons_vector;
    for (const auto &neuron_gid : gids)
    {
        neurons_vector.push_back(
            kernel().neuron_manager.get_neuron(neuron_gid));
    }

    for (auto const &neuron : neurons_vector)
    {
        Skeleton neuron_skel = Skeleton(neuron.get());

        // fill the neurites
        _fill_skel(neuron_skel.axon, axon, true);
        _fill_skel(neuron_skel.dendrites, dendrites, true);
        _fill_skel(neuron_skel.branching_points, nodes, false);
        _fill_skel(neuron_skel.growth_cones, growth_cones, false);

        somas[0].push_back(neuron_skel.soma_position.at(0));
        somas[1].push_back(neuron_skel.soma_position.at(1));
        somas[2].push_back(neuron_skel.soma_radius);
    }
#ifndef NDEBUG
    printf("the soma is in %f, %f", somas[0], somas[1]);
    printf(" %lu neurons has been imported for visualization \n"
           " the size of neurites vector is: %lu \n"
           " the size of axon vector is :    %lu \n"
           " the size of soma vector is :    %lu \n"
           " the size of growth_cones is:    %lu \n",
           neurons_vector.size(), dendrites.first.size(), axon.first.size(),
           somas[0].size(), growth_cones.first.size());
#endif
}

void get_swc(std::string output_file, std::vector<size_t> gids,
             unsigned int resolution)
{
#ifndef NDEBUG
    printf("Getting skeleton\n");
#endif
    std::sort(gids.begin(), gids.end());
    Swc swc(output_file, resolution);
    for (const auto &neuron_gid : gids)
    {
        const NeuronPtr neuron = kernel().neuron_manager.get_neuron(neuron_gid);
        // @todo: pass non default arguments
        swc.to_swc(neuron.get(), neuron_gid);
    }
    swc.close_file();
#ifndef NDEBUG
    printf("swc to file %s", output_file.c_str());
#endif
}


void simulate(const Time &simtime)
{
    kernel().simulation_manager.simulate(simtime);
}

void _fill_skel(const SkelNeurite &source_container,
                SkelNeurite &target_container, bool add_nan)
{
    auto it_x = target_container.first.end();
    auto it_y = target_container.second.end();
    target_container.first.insert(it_x, source_container.first.begin(),
                                  source_container.first.end());
    target_container.second.insert(it_y, source_container.second.begin(),
                                   source_container.second.end());
    if (add_nan)
    {
        target_container.first.push_back(std::nan(""));
        target_container.second.push_back(std::nan(""));
    }
}
void _fill_swc(const SkelNeurite &source_container,
               SkelNeurite &target_container, bool add_nan)
{
    auto it_x = target_container.first.end();
    auto it_y = target_container.second.end();
    target_container.first.insert(it_x, source_container.first.begin(),
                                  source_container.first.end());
    target_container.second.insert(it_y, source_container.second.begin(),
                                   source_container.second.end());
    if (add_nan)
    {
        target_container.first.push_back(std::nan(""));
        target_container.second.push_back(std::nan(""));
    }
}
std::vector<size_t> get_neurons() { return kernel().neuron_manager.get_gids(); }

// tool function
}
