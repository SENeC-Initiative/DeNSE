/*
 * kernel_manager.hpp
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

#ifndef KERNEL_M_H
#define KERNEL_M_H

#include <cassert>
#include <string>
#include <vector>

#include "config_impl.hpp"
#include "exceptions.hpp"

#include "models_manager.hpp"

#include "models_manager.hpp"
#include "neuron_manager.hpp"
#include "parallelism_manager.hpp"
#include "record_manager.hpp"
#include "rng_manager.hpp"
#include "simulation_manager.hpp"
#include "space_manager.hpp"


namespace growth
{

class KernelManager
{
  private:
    KernelManager();
    ~KernelManager();
    static KernelManager *kernel_manager_instance_;

    KernelManager(KernelManager const &);  // do not implement
    void operator=(KernelManager const &); // do not implement

  public:
    /*
     * Create/destroy and access the KernelManager singleton.
     */
    static void create_kernel_manager();
    static void destroy_kernel_manager();
    static KernelManager &get_kernel_manager();

    /**
     * Prepare kernel for operation.
     *
     * This method calls the initialization methods of the specific
     * managers in the proper order.
     *
     * @see finalize(), reset()
     */
    void initialize();

    /**
     * Take down kernel after operation.
     *
     * This method calls the finalization methods of the specific managers
     * in the proper order, i.e., inverse to initialize().
     *
     * @see initialize(), reset()
     */
    void finalize();

    /**
     * Reset kernel.
     *
     * Resets kernel by finalizing and initalizing.
     *
     * @see initialize(), finalize()
     */
    void reset();

    //! Get full kernel configuration.
    const statusMap get_status() const;
    void set_status(const statusMap &status);
    void set_simulation_ID(std::string simulation_ID);
    std::string get_simulation_ID() const;

    //! Get the number of objects in the simulator.
    stype get_num_created_objects() const;
    stype get_num_objects() const;
    //! Update the object count
    void update_num_objects(stype num_new_objects);

    std::string simulation_ID_;
    double get_adaptive_timestep() const;
    bool using_environment() const;
    //! Returns true if kernel is initialized
    bool is_initialized() const;
    //! Space manager instance
    ParallelismManager parallelism_manager;
    RNGManager rng_manager;
    SimulationManager simulation_manager;
    SpaceManager space_manager;
    RecordManager record_manager;
    ModelManager model_manager;
    NeuronManager neuron_manager;

  private:
    statusMap config_;  //!< configuration properties
    bool env_required_; //!< true if a spatial environment is required
    bool record_enabled_;
    //!< to create neurons
    bool initialized_;          //!< true if all sub-managers initialized
    stype num_created_objects_; //!< number of objects created
    stype num_objects_;         //!< current number of objects
    std::string version_;
    double adaptive_timestep_; //! if > 1, step divider when interacting
};

KernelManager &kernel();

} // namespace growth


// Inline implementations

inline growth::KernelManager &growth::KernelManager::get_kernel_manager()
{
    assert(kernel_manager_instance_);
    return *kernel_manager_instance_;
}


inline growth::KernelManager &growth::kernel()
{
    return KernelManager::get_kernel_manager();
}


inline bool growth::KernelManager::is_initialized() const
{
    return initialized_;
}

#endif // KERNEL_M_H
