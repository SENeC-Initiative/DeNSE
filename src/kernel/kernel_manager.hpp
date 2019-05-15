#ifndef KERNEL_M_H
#define KERNEL_M_H

#include <assert.h>
#include <string>
#include <vector>

#include "config_impl.hpp"
#include "exceptions.hpp"
#include "neuron_manager.hpp"
#include "parallelism_manager.hpp"
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
    std::string get_simulation_ID();

    //! Get the number of objects in the simulator.
    size_t get_num_objects() const;
    //! Update the object count
    void update_num_objects();

    std::string simulation_ID_;
    bool angles_in_radians() const;
    double get_current_seconds() const;
    bool using_environment() const;
    //! Returns true if kernel is initialized
    bool is_initialized() const;
    //! Space manager instance
    ParallelismManager parallelism_manager;
    RNGManager rng_manager;
    SimulationManager simulation_manager;
    SpaceManager space_manager;
    NeuronManager neuron_manager;
    std::string get_simulation_ID() const;

  private:
    bool angles_in_radians_; //!< true if angles should be passed in rad
    statusMap config_;       //!< configuration properties
    bool env_required_;      //!< true if a spatial environment is required
    bool record_enabled_;
    //!< to create neurons
    bool initialized_;   //!< true if all sub-managers initialized
    size_t num_objects_; //!< number of objects created
    std::string version_;
};

KernelManager &kernel();
}

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
