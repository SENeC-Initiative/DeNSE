#include "simulation_manager.hpp"

// C includes:
#include <cmath>
#include <sys/time.h>

// Includes from kernel
#include "GrowthCone.hpp"
#include "config_impl.hpp"
#include "exceptions.hpp"
#include "growth_time.hpp"
#include "kernel_manager.hpp"

// include from elements
#include "Neuron.hpp"


namespace growth
{

SimulationManager::SimulationManager()
    : simulating_(false)          //!< true if simulation in progress
    , step_()                     //!< Current step of the simulation
    , final_step_(0L)             //!< Last step, updated once per slice
    , initial_time_(Time())       //!< Initial time (day, hour, min, sec)
    , final_time_(Time())         //!< Final time (day, hour, min, sec)
    , maximal_time_(Time())       //!< Maximal time (day, hour, min, sec)
    , terminate_(false)           //!< Terminate on signal or error
    , resolution_scale_factor_(1) //! rescale step size respct to old resolution
{
}


void SimulationManager::initialize()
{
    // set resolution
    Time::reset_resolution();
}


/**
 * Update each elements with new step and resolution.
 * This function is called at each call to simulate.
 */
void SimulationManager::initialize_simulation_(const Time &t)
{
    assert(kernel().is_initialized());

    final_time_ = initial_time_ + t;

    // set the right number of step objects
    int num_omp = kernel().parallelism_manager.get_num_local_threads();
    step_       = std::vector<Time::timeStep>(num_omp, 0L);

#pragma omp parallel
    {
        int omp_id = kernel().parallelism_manager.get_thread_local_id();

#ifndef NDEBUG
        printf("\n\n\n ######## Starting simulation ########## \n\n");
#endif

        mtPtr rnd_engine = kernel().rng_manager.get_rng(omp_id);
        std::vector<NeuronPtr> local_neurons =
            kernel().neuron_manager.get_local_neurons(omp_id);

        // first, initialize neurons
        for (auto &neuron : local_neurons)
        {
            neuron->initialize_next_event(rnd_engine, resolution_scale_factor_,
                                          final_step_);
        }
    }
}


void SimulationManager::finalize()
{
    Time::reset_resolution();

    step_.clear();
    final_step_   = 0L;
    initial_time_ = Time();
    final_time_   = Time();
}


/**
 * Reset the SimulationManager to the state at T = 0.
 */
void SimulationManager::reset_culture() {}


void SimulationManager::test_random_generator(Random_vecs &values, size_t size)
{
    std::uniform_real_distribution<> uniform_(0, 1.);
#pragma omp parallel
    {
        int omp_id       = kernel().parallelism_manager.get_thread_local_id();
        mtPtr rnd_engine = kernel().rng_manager.get_rng(omp_id);
        std::vector<double> randoms;
        for (size_t n = 0; n < size; n++)
        {
            randoms.push_back(uniform_(*(rnd_engine).get()));
        }
#pragma omp critical
        {
            values.push_back(randoms);
        }
    }
}


/**
 * @brief Take change in number of threads into account
 */
void SimulationManager::num_threads_changed(int num_omp)
{
    Time::timeStep old_step = (step_.size() > 0) ? step_.front() : 0L;
    step_.clear();
    step_ = std::vector<Time::timeStep>(num_omp, old_step);
}


/**
 * Simulate for the given time .
 * This function performs the following steps
 * 1. set the new simulation time
 * 2. call prepare_simulation()
 * 3. call resume()
 * 4. call finalize_simulation()
 */
void SimulationManager::simulate(const Time &t)
{
    initialize_simulation_(t);

    Time::timeStep steps = Time::to_steps(t);

#pragma omp parallel
    {
        int omp_id = kernel().parallelism_manager.get_thread_local_id();

        mtPtr rnd_engine = kernel().rng_manager.get_rng(omp_id);
        std::vector<NeuronPtr> local_neurons =
            kernel().neuron_manager.get_local_neurons(omp_id);

        // then run the simulation
        while (step_[omp_id] < steps)
        {
            step_[omp_id]++;

            // update neurons
            for (auto &neuron : local_neurons)
            {
                neuron->grow(rnd_engine);
            }

            // record their state
            kernel().record_manager.record(omp_id);

#ifndef NDEBUG
            if (step_[omp_id] % 50 == 0)
            {
                printf("##simulation step is %lu \n", step_[omp_id]);
                printf("##simulation seconds is %f \n", get_current_seconds());
            }
#endif /* NDEBUG */
        }
    }

    finalize_simulation_(t);
}


void SimulationManager::finalize_simulation_(const Time &t)
{
    Time::timeStep steps = Time::to_steps(t);
    for (size_t s : step_)
    {
        assert(s == steps);
    }
    //    assert(initial_time_.get_total_seconds() ==
    //    final_time_.get_total_seconds());
    initial_time_.update(step_.front());

    // finalize recorder times
    kernel().record_manager.finalize_simulation(steps);

    final_step_ = step_.front();
}


//###################################################
//              Getter/setter functions
//###################################################

void SimulationManager::set_status(const statusMap &status)
{
    if (status.find("resolution") != status.end())
    {
        double resolution;
        get_param(status, "resolution", resolution);
        resolution_scale_factor_ = Time::RESOLUTION / resolution;
        Time::set_resolution(resolution);
    }
}


void SimulationManager::get_status(statusMap &status) const
{
    set_param(status, "resolution", Time::RESOLUTION);

    Time time = Time(initial_time_, final_step_);

    set_param(status, "second", time.get_sec());
    set_param(status, "minute", time.get_min());
    set_param(status, "hour", time.get_hour());
    set_param(status, "day", time.get_day());
}


double SimulationManager::get_current_seconds() const
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    return step_.at(omp_id) * Time::RESOLUTION +
           initial_time_.get_total_seconds();
}


Time SimulationManager::get_time() const
{
    assert(not simulating_);
    return initial_time_;
}


//! Get the timestep in second
double SimulationManager::get_resolution() const { return Time::RESOLUTION; }


size_t SimulationManager::get_current_step() const
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    return step_.at(omp_id);
}

} // namespace
