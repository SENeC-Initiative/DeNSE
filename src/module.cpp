/*
 * module.cpp
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

#include "module.hpp"


// C++ include
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>

#include <boost/range/adaptor/strided.hpp>

// elements include
#include "Branch.hpp"
#include "GrowthCone.hpp"
#include "Neurite.hpp"
#include "Neuron.hpp"
#include "Skeleton.hpp"
#include "Swc.hpp"

// kernel include
#include "kernel_manager.hpp"
#include "neuron_manager.hpp"
#include "models_manager.hpp"
#include "space_manager.hpp"


namespace ba = boost::adaptors;

const int INVALID_AXON_NAME(1);


namespace growth
{

size_t create_objects_(const std::string &object_name,
                      const std::vector<statusMap> &obj_params)
{
    size_t gid = 0;

    if (object_name == "recorder")
    {
        gid = kernel().record_manager.create_recorder(obj_params);
    }
    else
    {
        throw InvalidParameter("Creating something other than a recorder (`" +
                               object_name + "`) is not currently supported.",
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
size_t create_neurons_(const std::vector<statusMap> &neuron_params,
                       const std::vector<statusMap> &axon_params,
                       const std::vector<statusMap> &dendrites_params)
{
    size_t num_created = kernel().neuron_manager.create_neurons(
        neuron_params, axon_params, dendrites_params);

    // update max_resolution
    kernel().simulation_manager.set_max_resolution();

    return num_created;
}


void create_neurites_(const std::vector<size_t> &neurons, size_t num_neurites,
                      const std::vector<statusMap> &params,
                      const std::vector<std::string> &neurite_types,
                      const std::vector<double> &angles,
                      const std::vector<std::string> &names)
{
    int num_omp = kernel().parallelism_manager.get_num_local_threads();
    std::vector< std::vector<size_t> > omp_neuron_vec(num_omp);

    for (size_t i=0; i < neurons.size(); i++)
    {
        int omp_id = kernel().neuron_manager.get_neuron_thread(neurons[i]);
        omp_neuron_vec[omp_id].push_back(i);
    }

    // if (not names.empty())
    // {
    //     for (std::string name : names)
    //     {
    //         printf("%s\n", name.c_str());
    //     }
    // }
    // else
    // {
    //     printf("names empty\n");
    // }

#pragma omp parallel
    {
        int omp_id = kernel().parallelism_manager.get_thread_local_id();
        mtPtr rng  = kernel().rng_manager.get_rng(omp_id);
        std::string gc_model;
        bool gc_model_set;
        GCPtr gc_ptr;

        for (size_t i : omp_neuron_vec[omp_id])
        {
            NeuronPtr neuron = kernel().neuron_manager.get_neuron(neurons[i]);
            statusMap status = params[i];

            size_t existing_neurites = neuron->get_num_neurites();
            bool has_axon            = neuron->has_axon();

            gc_model_set = get_param(status, "growth_cone_model", gc_model);

            if (not gc_model_set)
            {
                gc_model = neuron->get_gc_model();
            }

            gc_ptr = kernel().model_manager.get_model(gc_model);

            if (neurite_types.empty())
            {
                for (size_t j=0; j < num_neurites; j++)
                {
                    if ((names.empty() and has_axon
                        and existing_neurites + j == 0) or
                        (not names.empty() and names[j] == "axon"))
                    {
                        neuron->new_neurite("axon", "axon", gc_ptr ,rng);
                        neuron->set_neurite_status("axon", status);
                    }
                    else
                    {
                        std::string name =
                            names.empty() ?
                            "dendrite_" + std::to_string(existing_neurites + j)
                            : names[j];
                        neuron->new_neurite(name, "dendrite", gc_ptr, rng);
                        neuron->set_neurite_status(name, status);
                    }
                }
            }
            else
            {
                for (size_t j=0; j < num_neurites; j++)
                {
                    if (neurite_types[j] == "axon")
                    {
                        if (not names.empty() and names[j] != "axon")
                        {
                            exit(INVALID_AXON_NAME);
                        }

                        neuron->new_neurite("axon", "axon", gc_ptr, rng);
                        neuron->set_neurite_status("axon", status);
                    }
                    else
                    {
                        std::string name =
                            names.empty() ?
                            "dendrite_" + std::to_string(existing_neurites + j)
                            : names[j];
                        neuron->new_neurite(name, "dendrite", gc_ptr, rng);
                        neuron->set_neurite_status(name, status);
                    }
                }
            }
        }
    }
}


void delete_neurons_(const std::vector<size_t> &gids)
{
    if (gids.empty())
    {
        kernel().neuron_manager.finalize();
        kernel().record_manager.finalize();

        kernel().neuron_manager.initialize();
        kernel().record_manager.initialize();
    }
    else
    {
        // first delete the pointers to the neurons
        kernel().record_manager.neurons_deleted(gids);
        // then delete the neurons
        kernel().neuron_manager.delete_neurons(gids);
    }
}


void delete_neurites_(const std::vector<size_t> &gids,
                      const std::vector<std::string> &names)
{
    if (gids.empty())
    {
        std::vector<NeuronPtr> neurons;

        kernel().neuron_manager.get_all_neurons(neurons);

        for (NeuronPtr neuron : neurons)
        {
            neuron->delete_neurites(names);
        }
    }
    else
    {
        for (size_t neuron : gids)
        {
            NeuronPtr nptr = kernel().neuron_manager.get_neuron(neuron);

            nptr->delete_neurites(names);
        }
    }
}


/*
 * Init/Finalize and simulation functions
 */

void init_growth_(int *argc, char **argv[])
{
    KernelManager::create_kernel_manager();
    kernel().parallelism_manager.mpi_init(argc, argv);
    kernel().initialize();
}


void finalize_growth_()
{
    kernel().parallelism_manager.mpi_finalize();
    KernelManager::destroy_kernel_manager();
}


void simulate_(const Time &simtime)
{
    try
    {
        kernel().simulation_manager.simulate(simtime);
    }
    catch (...)
    {
        std::string message = "Error occurred in `simulate`.\n";
        get_backtrace_(message, 0);

        throw std::runtime_error(message);
    }
}


void reset_kernel_() { kernel().reset(); }


std::string get_simulation_id_() { return kernel().get_simulation_ID(); }


/*
 * Environment
 */

void set_environment_(
    GEOSGeometry * environment, const std::vector<GEOSGeometry *> &areas,
    std::vector<double> heights, const std::vector<std::string> &names,
    const std::vector<std::unordered_map<std::string, double>> &properties)
{
    kernel().space_manager.set_environment(environment, areas, heights, names,
                                           properties);
}


void get_environment_(
    GEOSGeometry *&environment, std::vector<GEOSGeometry *> &areas,
    std::vector<double> &heights, std::vector<std::string> &names,
    std::vector<std::unordered_map<std::string, double>> &properties)
{
    kernel().space_manager.get_environment(environment, areas, heights, names,
                                           properties);
}


/*
 * Statuses
 */


void set_kernel_status_(const statusMap &status_dict, std::string simulation_ID)
{
    kernel().set_simulation_ID(simulation_ID);
    kernel().set_status(status_dict);
}


void set_status_(size_t gid, statusMap neuron_status, statusMap axon_status,
                 statusMap dendrites_status)
{
    // @todo: do this in parallel
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
    neuron->set_neurite_status("dendrite", local_params);

    // update max_resolution for simulation
    kernel().simulation_manager.set_max_resolution();
}


void set_neurite_status_(size_t neuron, std::string neurite, statusMap status)
{
    NeuronPtr nptr = kernel().neuron_manager.get_neuron(neuron);

    nptr->set_neurite_status(neurite, status);
}


/*
 * Models
 */


void get_models_(std::unordered_map<std::string, std::string> &models,
                 bool abbrev)
{
    kernel().model_manager.get_models(models, abbrev);
}


std::vector<std::string> get_elongation_types_()
{
    return kernel().model_manager.get_extension_types();
}


std::vector<std::string> get_steering_methods_()
{
    return kernel().model_manager.get_steering_methods();
}


std::vector<std::string> get_direction_selection_methods_()
{
    return kernel().model_manager.get_direction_selection_methods();
}


/*
 * Getter functions
 */

const Time get_current_time_()
{
    return kernel().simulation_manager.get_time();
}


statusMap get_kernel_status_() { return kernel().get_status(); }


statusMap get_status_(size_t gid)
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
        throw std::runtime_error("Only neurons an recorders are supported so "
                                 "far.");
    }
}


double get_state_(size_t gid, const std::string& level,
                  const std::string& variable)
{
    const char* cvar = variable.c_str();

    if (level == "neuron")
    {
        return kernel().neuron_manager.get_neuron(gid)->get_state(cvar);
    }
    else
    {
        return kernel().neuron_manager.get_neuron(gid)->get_neurite(level)
                   .lock()->get_state(cvar);
    }
}


size_t get_num_objects_() { return kernel().get_num_objects(); }


size_t get_num_created_objects_() { return kernel().get_num_created_objects(); }


statusMap get_neurite_status_(size_t gid, const std::string &neurite,
                             const std::string &level)
{
    return kernel().neuron_manager.get_neurite_status(gid, neurite, level);
}


std::string get_default_model_()
{
    return kernel().model_manager.default_model;
}


void get_defaults_(const std::string &object_name,
                   const std::string &object_type, const std::string &gc_model,
                   bool detailed, statusMap &status)
{
    GCPtr gc = nullptr;

    if (object_type == "growth_cone")
    {
        gc = kernel().model_manager.get_model(object_name);
        gc->get_status(status);
    }
    else if (object_type == "neuron" || object_type == "neurite")
    {
        kernel().neuron_manager.get_defaults(status, object_name, detailed);
        if (detailed)
        {
            std::string model = (gc_model == "") ? "default" : gc_model;
            gc = kernel().model_manager.get_model(model);
            gc->get_status(status);
        }
    }
    else if (object_type == "recorder")
    {
        // return default NeuronContinuousRecorder
        kernel().record_manager.get_defaults(status);
    }
}


bool is_neuron_(size_t gid)
{
    return kernel().neuron_manager.is_neuron(gid);
}


bool is_neurite_(size_t gid, const std::string& neurite)
{
    return kernel().neuron_manager.get_neuron(gid)->is_neurite(neurite);
}


std::string object_type_(size_t gid)
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


// ------------------------------------------------------------------------- //
// Neuron/structure related

std::vector<size_t> get_neurons_()
{
    return kernel().neuron_manager.get_gids();
}


std::vector<std::string> get_neurites_(size_t gid)
{
    std::vector<std::string> neurite_names;

    NeuronPtr n = kernel().neuron_manager.get_neuron(gid);

    std::unordered_map<std::string, NeuritePtr>::const_iterator it =
        n->neurite_cbegin();
    std::unordered_map<std::string, NeuritePtr>::const_iterator it_last =
        n->neurite_cend();

    while (it != it_last)
    {
        neurite_names.push_back(it->first);
        it++;
    }

    return neurite_names;
}


bool neuron_has_axon_(size_t gid)
{
    return kernel().neuron_manager.get_neuron(gid)->has_axon();
}


/*
 * Morphology data
 */

void get_branches_data_(size_t neuron, const std::string &neurite_name,
                        std::vector<std::vector<std::vector<double>>> &points,
                        std::vector<double> &diameters,
                        std::vector<int> &parents, std::vector<size_t> &nodes,
                        size_t start_point)
{
    NeuronPtr n            = kernel().neuron_manager.get_neuron(neuron);
    NeuriteWeakPtr neurite = n->get_neurite(neurite_name);

    size_t idx = 0;

    auto node_it  = neurite.lock()->nodes_cbegin();
    auto node_end = neurite.lock()->nodes_cend();

    while (node_it != node_end)
    {
        std::vector<std::vector<double>> points_tmp;
        BranchPtr b = node_it->second->get_branch();


        if (b->size() != 0)
        {
            const std::vector<double> &xx = b->get_xlist();
            const std::vector<double> &yy = b->get_ylist();

            switch (start_point)
            {
            case 0:
                points_tmp.push_back(xx);
                points_tmp.push_back(yy);
                break;
            default:
                std::vector<double> row_x, row_y;
                // x
                auto it  = xx.cbegin();
                auto end = xx.cend();
                row_x.insert(row_x.end(), it + start_point, end);
                points_tmp.push_back(row_x);
                // y
                it  = yy.cbegin();
                end = yy.cend();
                row_y.insert(row_y.end(), it + start_point, end);
                points_tmp.push_back(row_y);
            }

            parents.push_back(node_it->second->get_parent().lock()->get_node_id());
            nodes.push_back(node_it->second->get_node_id());

            points.push_back(points_tmp);
            diameters.push_back(node_it->second->get_diameter());
        }

        node_it++;
    }

    for (auto gc : neurite.lock()->gc_range())
    {
        std::vector<std::vector<double>> points_tmp;
        BranchPtr b = gc.second->get_branch();

        if (b->size() != 0)
        {
            points_tmp.push_back(b->get_xlist());
            points_tmp.push_back(b->get_ylist());

            size_t parent_id = gc.second->get_parent().lock()->get_node_id();

            parents.push_back(parent_id);
            nodes.push_back(gc.second->get_node_id());

            points.push_back(points_tmp);

            diameters.push_back(gc.second->get_diameter());
        }
    }
}


void get_skeleton_(SkelNeurite &axon, SkelNeurite &dendrites,
                   SkelNeurite &nodes, SkelNeurite &growth_cones,
                   SkelSomas &somas, std::vector<size_t> gids,
                   unsigned int resolution)
{
    std::vector<NeuronPtr> neurons_vector;
    for (const auto &neuron_gid : gids)
    {
        neurons_vector.push_back(
            kernel().neuron_manager.get_neuron(neuron_gid));
    }

    for (auto const &neuron : neurons_vector)
    {
        Skeleton neuron_skel = Skeleton(neuron.get(), resolution);

        // fill the neurites
        _fill_skel(neuron_skel.axon, axon, true);
        _fill_skel(neuron_skel.dendrites, dendrites, true);
        _fill_skel(neuron_skel.branching_points, nodes, false);
        _fill_skel(neuron_skel.growth_cones, growth_cones, false);

        somas[0].push_back(neuron_skel.soma_position.x());
        somas[1].push_back(neuron_skel.soma_position.y());
        somas[2].push_back(neuron_skel.soma_radius);
    }
#ifndef NDEBUG
    printf(" %lu neurons has been imported for visualization \n"
           " the size of neurites vector is: %lu \n"
           " the size of axon vector is :    %lu \n"
           " the size of soma vector is :    %lu \n"
           " the size of growth_cones is:    %lu \n",
           neurons_vector.size(), dendrites.first.size(), axon.first.size(),
           somas[0].size(), growth_cones.first.size());
#endif
}


void get_geom_skeleton_(std::vector<size_t> gids,
                        std::vector<GEOSGeometry*>& axons,
                        std::vector<GEOSGeometry*>& dendrites,
                        std::vector<size_t>& dendrite_gids,
                        std::vector< std::vector<double> >& somas)
{
    std::vector<GEOSGeometry *> vec;
    //~ std::vector<BPolygon> vec_tmp;
    //~ std::vector<BMultiPolygon> vec_geom;
    std::vector<BPolygon> vec_tmp, vec_geom;

    GEOSContextHandle_t ch = kernel().space_manager.get_context_handler();
    GEOSWKTReader *reader  = GEOSWKTReader_create_r(ch);
    GEOSGeometry *geom_tmp, *geom_union;
    std::stringstream s;
    BMultiPolygon mp;
    std::string wkt;
    size_t num_poly;

    for (size_t gid : gids)
    {
        NeuronPtr neuron = kernel().neuron_manager.get_neuron(gid);
        auto neurite_it  = neuron->neurite_cbegin();
        auto neurite_end = neuron->neurite_cend();

        while (neurite_it != neurite_end)
        {
            auto node_it  = neurite_it->second->nodes_cbegin();
            auto node_end = neurite_it->second->nodes_cend();

            vec.clear();
            vec_geom.clear();

            while (node_it != node_end)
            {
                for (BPolygonPtr p : node_it->second->segment_range())
                {
                    s.str("");
                    s << std::setprecision(12) << bg::wkt(*(p.get()));
                    geom_tmp = GEOSWKTReader_read_r(ch, reader, s.str().c_str());
                    vec.push_back(geom_tmp);
                }

                if (vec_tmp.size())
                {
                    vec_geom.push_back(vec_tmp.back());
                }
                node_it++;
            }

            for (auto gc : neurite_it->second->gc_range())
            {
                for (BPolygonPtr p : gc.second->segment_range())
                {
                    s.str("");
                    s << std::setprecision(12) << bg::wkt(*(p.get()));
                    geom_tmp = GEOSWKTReader_read_r(ch, reader, s.str().c_str());
                    vec.push_back(geom_tmp);
                }

                // to nicely finish the neurite, we add a disk to mark the
                // growth cone position
                BPolygon disk = kernel().space_manager.make_disk(
                    gc.second->get_position(), 0.5*gc.second->get_diameter());
                
                s.str("");
                s << std::setprecision(12) << bg::wkt(disk);
                geom_tmp = GEOSWKTReader_read_r(ch, reader, s.str().c_str());
                vec.push_back(geom_tmp);
            }

            //~ // create the stupid collection to make the union
            geom_tmp = GEOSGeom_createCollection_r(ch, GEOS_MULTIPOLYGON,
                                                   vec.data(), vec.size());
            geom_union = GEOSUnaryUnion_r(ch, geom_tmp);

            GEOSGeom_destroy_r(ch, geom_tmp);

            if (neurite_it->second->get_type() == "axon")
            {
                axons.push_back(geom_union);
                //~ axons.insert(axons.end(), vec.begin(), vec.end());
            }
            else
            {
                dendrites.push_back(geom_union);
                dendrite_gids.push_back(gid);
                //~ dendrites.insert(dendrites.end(), vec.begin(), vec.end());
            }

            vec.clear();

            neurite_it++;
        }

        BPoint soma = neuron->get_position();
        somas.push_back({soma.x(), soma.y(), neuron->get_soma_radius()});
    }
}


void generate_synapses_(
  bool crossings_only, double density, bool only_new_syn, bool autapse_allowed,
  const std::set<size_t> &presyn_pop, const std::set<size_t> &postsyn_pop,
  std::vector<size_t> &presyn_neurons, std::vector<size_t> &postsyn_neurons,
  std::vector<std::string> &presyn_neurites,
  std::vector<std::string> &postsyn_neurites,
  std::vector<size_t> &presyn_nodes, std::vector<size_t> &postsyn_nodes,
  std::vector<size_t> &presyn_segments, std::vector<size_t> &postsyn_segments,
  std::vector<double> &pre_syn_x, std::vector<double> &pre_syn_y,
  std::vector<double> &post_syn_x, std::vector<double> &post_syn_y)
{
    if (crossings_only)
    {
        kernel().space_manager.generate_synapses_crossings(
            density, only_new_syn, autapse_allowed, presyn_pop, postsyn_pop,
            presyn_neurons, postsyn_neurons, presyn_neurites, postsyn_neurites,
            presyn_nodes, postsyn_nodes, presyn_segments, postsyn_segments,
            pre_syn_x, pre_syn_y);
        
        post_syn_x = pre_syn_x;
        post_syn_y = pre_syn_y;
    }
    else
    {
        kernel().space_manager.generate_synapses_all(
            density, only_new_syn, autapse_allowed, presyn_pop, postsyn_pop,
            presyn_neurons, postsyn_neurons, presyn_neurites, postsyn_neurites,
            presyn_nodes, postsyn_nodes, presyn_segments, postsyn_segments,
            pre_syn_x, pre_syn_y, post_syn_x, post_syn_y);
    }
}


void get_distances_(size_t gid, const std::string &neurite_name, size_t node,
                    size_t segment, double &dist_to_parent,
                    double &dist_to_soma)
{
    NeuronPtr neuron = kernel().neuron_manager.get_neuron(gid);

    if (neurite_name.empty())
    {
        // we got the soma
        dist_to_parent = 0.;
        dist_to_soma   = 0.;
    }
    else
    {
        NeuriteWeakPtr neurite = neuron->get_neurite(neurite_name);

        neurite.lock()->get_distances(node, segment, dist_to_parent,
                                      dist_to_soma);
    }
}


void get_swc_(std::string output_file, std::vector<size_t> gids,
              unsigned int resolution)
{
    std::sort(gids.begin(), gids.end());
    Swc swc(output_file, resolution);

    for (const auto &neuron_gid : gids)
    {
        const NeuronPtr neuron = kernel().neuron_manager.get_neuron(neuron_gid);
        // @todo: pass non default arguments
        swc.to_swc(neuron.get(), neuron_gid);
    }

    swc.close_file();
}


// ------------------------------------------------------------------------- //
// Recorder related

bool get_next_recording_(size_t gid, std::vector<Property> &ids,
                         std::vector<double> &values)
{
    return kernel().record_manager.get_next_recording(gid, ids, values);
}


bool get_next_time_(size_t gid, std::vector<Property> &ids,
                    std::vector<double> &values, const std::string &time_units)
{
    return kernel().record_manager.get_next_time(gid, ids, values, time_units);
}


void get_recorder_type_(size_t gid, std::string &level, std::string &event_type)
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


bool walk_neurite_tree_(size_t neuron, std::string neurite, NodeProp& np)
{
    NeuronPtr n                = kernel().neuron_manager.get_neuron(neuron);
    NeuriteWeakPtr neurite_ptr = n->get_neurite(neurite);

    return neurite_ptr.lock()->walk_tree(np);
}


void get_backtrace_(std::string &msg, int depth = 0)
{
    try
    {
        throw;
    }
    catch (const std::exception &e)
    {
        // handle exceptions of known type
        msg += std::to_string(depth) + ": [" + std::string(typeid(e).name()) +
               "] " + std::string(e.what()) + "\n";
        try
        {
            std::rethrow_if_nested(e);
        }
        catch (...)
        {
            get_backtrace_(msg, ++depth);
        }
    }
    catch (const std::nested_exception &ne)
    {
        // Not all nesting exceptions will be of a known type, but if they use
        // the mixin type std::nested_exception, then we can at least handle
        // them enough to get the nested exception:
        msg += std::to_string(depth) + ": Unknown nested exception\n";

        try
        {
            ne.rethrow_nested();
        }
        catch (...)
        {
            get_backtrace_(msg, ++depth);
        }
    }
    catch (...)
    {
        // Exception nesting works through inheritance, hence if the type cannot
        // be inherited std::nested_exception will not work.
        // Trying something like std::throw_with_nested( int{10} ) will hit this
        // final catch block.
        msg += std::to_string(depth) + ": Unknown nested exception\n";
    }
}


void test_random_generator_(Random_vecs &values, size_t size)
{
    kernel().simulation_manager.test_random_generator(values, size);
    printf("%lu number generated from rng\n", values[0].size());
}

} /* namespace */
