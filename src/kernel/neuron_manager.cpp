/*
 * neuron_manager.cpp
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

#include "neuron_manager.hpp"

#include <cmath>
#include <cstring>
#include <limits>
#include <mutex>

#include "config.hpp"
#include "config_impl.hpp"
#include "kernel_manager.hpp"

#include "Neuron.hpp"


namespace growth
{

NeuronManager::NeuronManager() {}


void NeuronManager::initialize()
{
    // create default neuron with two neurites (axon + dendrite)
    // set growth cone model to resource-based to get these parameters
    num_created_neurons_ = 0;
    statusMap empty_params;
    statusMap params({{names::num_neurites, Property(2, "")},
                      {names::growth_cone_model,
                       Property("resource-based_pull-only_run-and-tumble", "")},
                      {"x", Property(0., "micrometer")},
                      {"y", Property(0., "micrometer")}});

    // set their status to use all possible parameters to have them all
    // when using GetDefaults
    std::vector<std::string> options({"use_van_pelt", "use_uniform_split",
                                      "use_uniform_branching",
                                      "use_flpl_branching"});

    for (std::string opt : options)
    {
        Property p  = Property(true, "");
        params[opt] = p;
    }

    // get random engine
    mtPtr rnd_ptr = kernel().rng_manager.get_rng(0);

    // create default neuron
    model_neuron_ = std::make_shared<Neuron>(0);
    model_neuron_->init_status(params, empty_params, empty_params, rnd_ptr);

    // remove information about angles
    model_neuron_->neurite_angles_.clear();
}


void NeuronManager::finalize()
{
    neurons_.clear();
    neurons_on_thread_.clear();
    model_map_.clear();
    num_created_neurons_ = 0;
}


/**
 * @brief Create new neurons with custom parameters.
 *
 * @return Number of objects created.
 */
stype NeuronManager::create_neurons(
    const std::vector<statusMap> &neuron_params,
    const std::unordered_map<std::string, std::vector<statusMap>> &neurite_params)
{
    stype first_id             = kernel().get_num_created_objects();
    stype previous_num_neurons = neurons_.size();

    // put the neurons on the thread list they belong to
    stype num_omp = kernel().parallelism_manager.get_num_local_threads();
    int omp_id    = kernel().parallelism_manager.get_thread_local_id();
    std::vector<std::vector<stype>> thread_neurons(num_omp);

    for (stype i = 0; i < neuron_params.size(); i++)
    {
        double x, y;
        get_param(neuron_params[i], "x", x);
        get_param(neuron_params[i], "y", y);

        if (kernel().space_manager.has_environment())
        {
            // printf("check point %f, %f in the environment \n",x,y);
            if (not kernel().space_manager.env_contains(BPoint(x, y)))
            {
                throw std::runtime_error(
                    " a Neuron was positioned outside the environment\n");
            }
        }

        // @TODO change temporary round-robin for neuron assignement
        int omp_id      = (first_id + i) % num_omp;
        stype neuron_id = first_id + i;
        thread_neurons[omp_id].push_back(neuron_id);
        thread_of_neuron_[neuron_id] = omp_id;
    }

    // exception_capture_flag "guards" captured_exception. std::called_once()
    // guarantees that will only execute any of its Callable(s) ONCE for each
    // unique std::once_flag. See C++11 Standard Library documentation
    // (``<mutex>``). These tools together ensure that we can capture exceptions
    // from OpenMP parallel regions in a thread-safe way.
    std::once_flag exception_capture_flag;
    // pointer-like object that manages an exception captured with
    // std::capture_exception(). We use this to capture exceptions thrown from
    // the OpenMP parallel region.
    std::exception_ptr captured_exception;

// create the neurons on the respective threads
#pragma omp parallel
    {
        std::vector<NeuronPtr> local_neurons;
        int omp_id       = kernel().parallelism_manager.get_thread_local_id();
        mtPtr rnd_engine = kernel().rng_manager.get_rng(omp_id);
        std::vector<stype> gids(thread_neurons[omp_id]);
        std::unordered_map<std::string, statusMap> neurite_status;

        for (stype gid : gids)
        {
            stype idx        = gid - first_id;
            NeuronPtr neuron = std::make_shared<Neuron>(gid);

            for (auto entry : neurite_params)
            {
                neurite_status[entry.first] = entry.second[idx]
            }

            try
            {
                neuron->init_status(neuron_params[idx], neurite_status,
                                    rnd_engine);
            }
            catch (const std::exception &except)
            {
                std::call_once(exception_capture_flag, [&captured_exception]() {
                    captured_exception = std::current_exception();
                });
            }

            local_neurons.push_back(neuron);
        }

#pragma omp critical
        {
            for (stype i = 0; i < gids.size(); i++)
            {
                neurons_.insert({gids[i], local_neurons[i]});
                neurons_on_thread_[omp_id].push_back(local_neurons[i]);
            }
        }
    }

    // check if an exception was thrown there
    if (captured_exception != nullptr)
    {
        // rethrowing nullptr is illegal
        std::rethrow_exception(captured_exception);
    }

    // tell the kernel manager to update the number of objects
    stype num_created = neurons_.size() - previous_num_neurons;
    kernel().update_num_objects(num_created);

    return num_created;
}


void NeuronManager::delete_neurons(const std::vector<stype> &gids)
{
    for (stype neuron : gids)
    {
        auto it = neurons_.find(neuron);

        NeuronPtr n = neurons_[neuron];

        if (it != neurons_.end())
        {
            for (auto v : neurons_on_thread_)
            {
                for (stype i = 0; i < v.size(); i++)
                {
                    if (v[i] == n)
                    {
                        v.erase(v.begin() + i);
                        break;
                    }
                }
            }

            thread_of_neuron_.erase(neuron);
            neurons_.erase(it);
        }
    }
}


void NeuronManager::init_neurons_on_thread(unsigned int num_local_threads)
{
    assert(neurons_.size() == 0); // no changes once neurons exist
    neurons_on_thread_ = std::vector<std::vector<NeuronPtr>>(num_local_threads);
    max_resolutions_   = std::vector<std::unordered_map<stype, double>>(
        num_local_threads, std::unordered_map<stype, double>());
}


void NeuronManager::update_kernel_variables()
{
    // exception_capture_flag "guards" captured_exception. std::called_once()
    // guarantees that will only execute any of its Callable(s) ONCE for each
    // unique std::once_flag. See C++11 Standard Library documentation
    // (``<mutex>``). These tools together ensure that we can capture exceptions
    // from OpenMP parallel regions in a thread-safe way.
    std::once_flag exception_capture_flag;
    // pointer-like object that manages an exception captured with
    // std::capture_exception(). We use this to capture exceptions thrown from
    // the OpenMP parallel region.
    std::exception_ptr captured_exception;

#pragma omp parallel
    {
        try
        {
            int omp_id = kernel().parallelism_manager.get_thread_local_id();

            gidNeuronMap local_neurons = get_local_neurons(omp_id);

            for (auto &neuron : local_neurons)
            {
                neuron.second->update_kernel_variables();
            }
        }
        catch (const std::exception &except)
        {
            std::call_once(exception_capture_flag, [&captured_exception]() {
                captured_exception = std::current_exception();
            });
        }
    }

    // check if an exception was thrown there
    if (captured_exception != nullptr)
    {
        // rethrowing nullptr is illegal
        std::rethrow_exception(captured_exception);
    }
}


// getters

NeuronPtr NeuronManager::get_neuron(stype gid) { return neurons_[gid]; }


gidNeuronMap NeuronManager::get_local_neurons(int local_thread_id)
{
    gidNeuronMap local_neurons;

    for (auto &n : neurons_on_thread_[local_thread_id])
    {
        local_neurons[n->get_gid()] = n;
    }

    return local_neurons;
}


void NeuronManager::get_all_neurons(std::vector<NeuronPtr> &neuron_ptr_vec)
{
    for (const auto &neuron : neurons_)
    {
        neuron_ptr_vec.push_back(neuron.second);
    }
}


void NeuronManager::get_defaults(statusMap &status, const std::string &object,
                                 bool detailed) const
{
    if (object == "neuron")
    {
        if (detailed)
        {
            model_neuron_->get_neurite_status(status, "dendrite", "neurite");
            model_neuron_->get_neurite_status(status, "axon", "neurite");
        }
        model_neuron_->get_status(status);
    }
    else if (object == "axon")
    {
        model_neuron_->get_neurite_status(status, "axon", "neurite");
    }
    else if (object == "dendrite" || object == "neurite")
    {
        model_neuron_->get_neurite_status(status, "dendrite", "neurite");
    }
}


const statusMap NeuronManager::get_neuron_status(stype gid) const
{
    NeuronPtr n = neurons_.at(gid);
    statusMap neuron_status;

    n->get_status(neuron_status);

    return neuron_status;
}


std::vector<stype> NeuronManager::get_gids() const
{
    std::vector<stype> gids;

    for (auto it : neurons_)
    {
        gids.push_back(it.first);
    }

    return gids;
}


const statusMap
NeuronManager::get_neurite_status(stype gid, const std::string &neurite,
                                  const std::string &level) const
{
    statusMap status;
    neurons_.at(gid)->get_neurite_status(status, neurite, level);
    return status;
}


bool NeuronManager::is_neuron(stype gid) const
{
    auto it = neurons_.find(gid);

    if (it != neurons_.end())
    {
        return true;
    }

    return false;
}


gidNeuronMap::const_iterator NeuronManager::iter_neurons()
{
    return neurons_.begin();
}


stype NeuronManager::num_neurons() const { return neurons_.size(); }


int NeuronManager::get_neuron_thread(stype gid) const
{
#ifdef MSVC
    // insane error with MSVC requires unsigned long instead of stype
    unsigned long ul_gid = gid;
    return thread_of_neuron_.at(ul_gid);
#else
    return thread_of_neuron_.at(gid);
#endif
}


void NeuronManager::register_model(std::string model_name, GCPtr model_ptr)
{
    model_map_[model_name] = model_ptr;
}


void NeuronManager::set_max_resol(stype neuron, double max_resol)
{
    max_resolutions_[thread_of_neuron_[neuron]][neuron] = max_resol;
}


double NeuronManager::get_max_resol() const
{
    double max_resol = std::numeric_limits<double>::max();

    for (const auto &map : max_resolutions_)
    {
        for (const auto &pair : map)
        {
            max_resol = std::min(max_resol, pair.second);
        }
    }

    return max_resol;
}

} // namespace growth
