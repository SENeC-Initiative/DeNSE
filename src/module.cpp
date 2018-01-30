#include "module.hpp"


// C++ includes
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

// elements includes
#include "GrowthCone.hpp"
#include "Neuron.hpp"
#include "Skeleton.hpp"
#include "Swc.hpp"

// kernel include
#include "kernel_manager.hpp"
#include "neuron_manager.hpp"
#include "space_manager.hpp"


namespace growth
{

size_t create_objects(const std::string &object_name,
                      const std::vector<statusMap> &obj_params)
{
    size_t gid = 0;

    if (object_name == "recorder")
    {
        gid = kernel().record_manager.create_recorder(obj_params);
    }
    else
    {
        throw InvalidParameter("Creating something other than a neuron or a "
                               " recorder (`" +
                                   object_name +
                                   "`) is not "
                                   "currently supported.",
                               __FUNCTION__, __FILE__, __LINE__);
    }

    return gid;
}


/**
 * @brief Create neurons and set status
 *
 * The number of created neurons depends on the size of the vector.
 *
 * @param neuron_params set the parameters for neuron and all neurites.
 * @param axon_params overwrite the main parameters from `neuron_params` for
 *                    the axon, empty vector is default.
 * @param dendrites_params overwrite the main parameters from `neuron_params`
 *                         for dendrites, empty vector is default.
 *
 * @return the gid of the neuron created.
 */
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
    {
        return kernel().neuron_manager.get_neuron_status(gid);
    }
    else if (kernel().record_manager.is_recorder(gid))
    {
        return kernel().record_manager.get_recorder_status(gid);
    }
    else
    {
        throw std::runtime_error("Only neurons an recorders are supprted so "
                                 "far.");
    }
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

    if (object_type == "growth_cone")
    {
        gc = kernel().neuron_manager.get_model(object_name);
        gc->get_status(status);
        gc->get_status(status);
    }
    else if (object_type == "neuron" || object_type == "neurite")
    {
        kernel().neuron_manager.get_defaults(status, object_name);
    }
    else if (object_type == "recorder")
    {
        // return default NeuronContinuousRecorder
        kernel().record_manager.get_defaults(status);
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
    kernel().parallelism_manager.mpi_init(argc, argv);
    kernel().initialize();
    printf("models initialized\n");
}


void finalize_growth()
{
    kernel().parallelism_manager.mpi_finalize();
    KernelManager::destroy_kernel_manager();
}


std::string object_type(size_t gid)
{
    if (kernel().neuron_manager.is_neuron(gid))
    {
        return "neuron";
    }
    else if (kernel().record_manager.is_recorder(gid))
    {
        return "recorder";
    }
    else
    {
        throw std::runtime_error("Object does not exist.");
    }
}


void reset_kernel() { kernel().reset(); }


void set_kernel_status(const statusMap &status_dict, std::string simulation_ID)
{
    kernel().set_simulation_ID(simulation_ID);
    kernel().set_status(status_dict);
}


std::string get_simulation_ID() { return kernel().get_simulation_ID(); }


void set_environment(
    GEOSGeom environment, const std::vector<GEOSGeom> &walls,
    const std::vector<GEOSGeom> &areas, std::vector<double> heights,
    const std::vector<std::string> &names,
    const std::vector<std::unordered_map<std::string, double>> &properties)
{
    kernel().space_manager.set_environment(environment, walls, areas, heights,
                                           names, properties);
}


void get_environment(
    GEOSGeom &environment, std::vector<GEOSGeom> &areas,
    std::vector<double> &heights, std::vector<std::string> &names,
    std::vector<std::unordered_map<std::string, double>> &properties)
{
    kernel().space_manager.get_environment(environment, areas, heights, names,
                                           properties);
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


std::vector<size_t> get_neurons() { return kernel().neuron_manager.get_gids(); }


bool get_next_recording(size_t gid, std::vector<Property> &ids,
                        std::vector<double> &values)
{
    return kernel().record_manager.get_next_recording(gid, ids, values);
}


bool get_next_time(size_t gid, std::vector<Property> &ids,
                   std::vector<double> &values, const std::string &time_units)
{
    return kernel().record_manager.get_next_time(gid, ids, values, time_units);
}


void get_recorder_type(size_t gid, std::string &level, std::string &event_type)
{
    kernel().record_manager.get_recorder_type(gid, level, event_type);
}


/* tool functions */

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

} /* namespace */
