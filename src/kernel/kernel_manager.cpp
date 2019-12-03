/*
 * kernel_manager.cpp
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

#include "kernel_manager.hpp"

#include <stdio.h>

// kernel includes
#include "neuron_manager.hpp"
#include "parallelism_manager.hpp"
#include "record_manager.hpp"
#include "rng_manager.hpp"
#include "simulation_manager.hpp"
#include "space_manager.hpp"

// models include
#include "models_manager.hpp"

// spatial include
#include "Environment.hpp"


namespace growth
{

KernelManager *KernelManager::kernel_manager_instance_ = 0;


/*
 * Instantiate KernelManager
 */
void KernelManager::create_kernel_manager()
{
// only one thread at a time must enter this, which means only one
// will indeed be created (prevents race condition)
#pragma omp critical(create_kernel_manager)
    {
        if (kernel_manager_instance_ == 0)
        {
            kernel_manager_instance_ = new KernelManager();
            assert(kernel_manager_instance_);
        }
    }
}


/*
 * Remove KernelManager
 */
void KernelManager::destroy_kernel_manager()
{
    kernel_manager_instance_->finalize();
    delete kernel_manager_instance_;
}


/*
 * Constructor and destructor
 */
KernelManager::KernelManager()
    // managers
    : parallelism_manager()
    , rng_manager()
    , simulation_manager()
    , space_manager()
    , record_manager()
    , model_manager()
    , neuron_manager()
    // status
    , initialized_(false)
    // settings
    , simulation_ID_("Kernel_ID_00000000")
    , record_enabled_(false)
    , env_required_(false)
    , num_objects_(0)
    , num_created_objects_(0)
    , adaptive_timestep_(-1.)
    , version_("0.1.0")
{
}


KernelManager::~KernelManager() {}


/*
 * Init, finalize and reset
 */
void KernelManager::initialize()
{
    // parralelism comes first
    parallelism_manager.initialize();

    // then RNG
    rng_manager.initialize();

    // then the rest
    simulation_manager.initialize();
    space_manager.initialize();
    record_manager.initialize();

    // models_manager init_models() must come before neuron_manager
    model_manager.init_models();
    neuron_manager.initialize();

    num_objects_         = 0;
    num_created_objects_ = 0;
    initialized_         = true;
}


void KernelManager::finalize()
{
    initialized_ = false;
    neuron_manager.finalize();
    record_manager.finalize();
    space_manager.finalize();
    simulation_manager.finalize();
    rng_manager.finalize();
    parallelism_manager.finalize();
}


void KernelManager::reset()
{
    finalize();
    initialize();
}


// Getters

/*
 * Return full config
 */
const statusMap KernelManager::get_status() const
{
    assert(is_initialized());

    statusMap status;
    // local
    set_param(status, "version", version_, "");
    set_param(status, "environment_required", env_required_, "");
    set_param(status, "record_enabled", record_enabled_, "");
    set_param(status, "adaptive_timestep", adaptive_timestep_, "");

    // set_param(status, "simulation_ID", simulation_ID_);
    /*
     * delegate the rest; no set_status for "neuron_manager", "record_manager"
     */
    parallelism_manager.get_status(status);
    rng_manager.get_status(status);
    space_manager.get_status(status);
    simulation_manager.get_status(status);

    return status;
}


void KernelManager::set_simulation_ID(std::string simulation_ID)
{
    simulation_ID_ = simulation_ID;
}


std::string KernelManager::get_simulation_ID() const { return simulation_ID_; }


/**
 * Returns the number of objects
 */
stype KernelManager::get_num_objects() const { return num_objects_; }


stype KernelManager::get_num_created_objects() const
{
    return num_created_objects_;
}


double KernelManager::get_adaptive_timestep() const
{
    return adaptive_timestep_;
}


bool KernelManager::using_environment() const { return env_required_; }


void KernelManager::update_num_objects(stype num_new_objects)
{
    num_objects_ = neuron_manager.num_neurons();
    num_objects_ += record_manager.num_recorders();

    num_created_objects_ += num_new_objects;
}


// Setters

void KernelManager::set_status(const statusMap &status)
{
    // local
    get_param(status, "version", version_);
    get_param(status, "record_enabled", record_enabled_);
    get_param(status, "simulation_ID", simulation_ID_);

    double at_old = adaptive_timestep_;
    bool at_updated =
        get_param(status, "adaptive_timestep", adaptive_timestep_);

    bool env_required_old = env_required_;
    bool env_updated = get_param(status, "environment_required", env_required_);

    /*
     * delegate the rest; no set_status for:
     * - rng_manager
     * - neuron_manager
     * - record_manager
     */
    parallelism_manager.set_status(status);
    space_manager.set_status(status);
    double old_resol = simulation_manager.get_resolution();
    simulation_manager.set_status(status);

    // update the objects
    env_updated *= (env_required_old != env_required_);
    at_updated *= (at_old != adaptive_timestep_);

    if (env_updated or at_updated)
    {
        neuron_manager.update_kernel_variables();
    }
}
} // namespace growth
