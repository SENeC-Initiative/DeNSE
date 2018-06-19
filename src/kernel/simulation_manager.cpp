#include "simulation_manager.hpp"

// C includes:
#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
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
    : simulating_(false)    //!< true if simulation in progress
    , step_()               //!< Current step of the simulation
    , substep_()            //!< Precise time inside current step
    , initial_step_(0)      //!< first step, updated once per slice
    , final_step_(0)        //!< Last step, updated once per slice
    , initial_time_(Time()) //!< Initial time (day, hour, min, sec)
    , final_time_(Time())   //!< Final time (day, hour, min, sec)
    , maximal_time_(Time()) //!< Maximal time (day, hour, min, sec)
    , terminate_(false)     //!< Terminate on signal or error
    , previous_resolution_(Time::RESOLUTION)
    , resolution_scale_factor_(1) //! rescale step size respct to old resolution
    , max_resol_(std::numeric_limits<double>::max())
{
}


void SimulationManager::initialize()
{
    // set resolution
    Time::reset_resolution();
}


void SimulationManager::finalize()
{
    Time::reset_resolution();

    step_.clear();
    substep_.clear();

    initial_step_ = 0;
    final_step_   = 0;
    initial_time_ = Time();
    final_time_   = Time();
    max_resol_    = std::numeric_limits<double>::max();
}


/**
 * Reset the SimulationManager to the state at T = 0.
 */
void SimulationManager::reset_culture() {}


void SimulationManager::test_random_generator(Random_vecs &values, size_t size)
{
    std::uniform_real_distribution<> uniform_(0, 1.);

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


/**
 * @brief Take change in number of threads into account
 */
void SimulationManager::num_threads_changed(int num_omp)
{
    Time::timeStep old_step = (step_.size() > 0) ? step_.front() : 0L;

    step_.clear();
    substep_.clear();

    step_    = std::vector<Time::timeStep>(num_omp, old_step);
    substep_ = std::vector<double>(num_omp, 0.);
}


/**
 * Update each elements with new step and resolution.
 * This function is called at each call to simulate.
 */
void SimulationManager::initialize_simulation_(const Time &t)
{
    assert(kernel().is_initialized());
    resolution_scale_factor_ = previous_resolution_ / Time::RESOLUTION;

    // @todo: remove final_step and reset step_ to zero everytime, use Time
    // objects for discrete events
    final_time_ = initial_time_ + t;
    final_step_ += Time::to_steps(t);

    // set the right number of step objects
    int num_omp = kernel().parallelism_manager.get_num_local_threads();
    step_       = std::vector<Time::timeStep>(num_omp, initial_step_);
    substep_    = std::vector<double>(num_omp, 0.);

    // reset branching events
    branching_ev_.clear();
    branching_ev_tmp_.clear();

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

#ifndef NDEBUG
            if (omp_id == 0)
                printf("\n\n\n ######## Starting simulation ########## \n\n");
#endif

            mtPtr rnd_engine = kernel().rng_manager.get_rng(omp_id);
            gidNeuronMap local_neurons =
                kernel().neuron_manager.get_local_neurons(omp_id);

            // first, initialize neurons
            for (auto &neuron : local_neurons)
            {
                neuron.second->initialize_next_event(
                    rnd_engine, resolution_scale_factor_, initial_step_);
            }

            // record the first step if time is zero
            if (initial_time_ == Time())
            {
                kernel().record_manager.record(omp_id);
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


void SimulationManager::finalize_simulation_()
{
    for (size_t s : step_)
    {
        assert(s == final_step_);
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

    // finalize neurons#pragma omp parallel
#pragma omp parallel
    {
        try
        {
            int omp_id = kernel().parallelism_manager.get_thread_local_id();

            gidNeuronMap local_neurons =
                kernel().neuron_manager.get_local_neurons(omp_id);

            for (auto &neuron : local_neurons)
            {
                neuron.second->finalize();
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

    // finalize recorder times
    kernel().record_manager.finalize_simulation(final_step_);

    //! IMPORTANT: THIS UPDATE MUST COME LAST!
    initial_time_.update(final_step_ - initial_step_);
    initial_step_ = final_step_;
}


/**
 * @brief add a branching event
 *
 * Add an event into `branching_ev_` and sort the vector.
 */
void SimulationManager::new_branching_event(const Event &ev)
{
#pragma omp critical
    branching_ev_tmp_.push_back(ev);
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
            size_t step_next_ev, current_step;
            double previous_substep = 0.;
            bool new_step           = false;
            bool branched           = false;

            mtPtr rnd_engine = kernel().rng_manager.get_rng(omp_id);
            gidNeuronMap local_neurons =
                kernel().neuron_manager.get_local_neurons(omp_id);

            // then run the simulation
            while (step_[omp_id] < final_step_)
            {
                current_step     = step_[omp_id];
                previous_substep = substep_[omp_id];
                branched         = false;

#pragma omp barrier

#pragma omp single
                {
                    if (!branching_ev_tmp_.empty())
                    {
                        branching_ev_.insert(branching_ev_.end(),
                                             branching_ev_tmp_.begin(),
                                             branching_ev_tmp_.end());

                        branching_ev_tmp_.clear();

                        std::sort(branching_ev_.begin(), branching_ev_.end(),
                                  ev_greater);
                    }
                }

#pragma omp barrier

                // check when the next event will occur
                step_next_ev = branching_ev_.empty()
                                   ? current_step + 1
                                   : std::get<0>(branching_ev_.back());

                // compute substep accordingly
                if (step_next_ev == current_step)
                {
                    substep_[omp_id] = std::get<1>(branching_ev_.back());
                    new_step         = (substep_[omp_id] == Time::RESOLUTION);
                }
                else
                {
                    substep_[omp_id] = Time::RESOLUTION;
                    new_step         = true;
                }

                if (substep_[omp_id] < 0)
                {
                    printf("current step %lu - new_step %i - substep %f\n",
                           current_step, new_step, substep_[omp_id]);
                }
                assert(substep_[omp_id] >= 0.);

                // update neurons
                for (auto &neuron : local_neurons)
                {
                    neuron.second->grow(rnd_engine, current_step,
                                        substep_[omp_id] - previous_substep);
                }

                if (step_next_ev == current_step)
                {
                    // someone has to branch
                    Event &ev            = branching_ev_.back();
                    size_t gid_branching = std::get<2>(ev);
                    auto it              = local_neurons.find(gid_branching);

                    if (it != local_neurons.end())
                    {
                        branched = local_neurons[gid_branching]->branch(
                            rnd_engine, ev);

                        // tell recorder manager
                        if (branched)
                        {
                            kernel().record_manager.new_branching_event(ev);
                        }

                        branching_ev_.pop_back();
                    }
                }

                if (new_step)
                {
                    // full step is completed, record, reset substep_, incr.
                    // step_
                    kernel().record_manager.record(omp_id);
                    substep_[omp_id] = 0.;
                    step_[omp_id]++;
                }

#ifndef NDEBUG
                if (omp_id == 0 and step_[0] % 50 == 0)
                {
                    printf("##simulation step is %lu \n", step_[omp_id]);
                    printf("##simulation seconds is %f \n",
                           get_current_seconds());
                }
#endif /* NDEBUG */
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

    finalize_simulation_();
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

        if (initial_time_ == Time(0., 0, 0, 0))
        {
            previous_resolution_ = resolution;
            // why did we decide to do that?
        }

        Time::set_resolution(resolution);
    }
}


void SimulationManager::get_status(statusMap &status) const
{
    set_param(status, "resolution", Time::RESOLUTION);
    set_param(status, "max_allowed_resolution", max_resol_);

    Time time = Time(initial_time_, final_step_);

    set_param(status, "second", time.get_sec());
    set_param(status, "minute", time.get_min());
    set_param(status, "hour", time.get_hour());
    set_param(status, "day", time.get_day());
}


double SimulationManager::get_current_seconds() const
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    return (step_.at(omp_id) - initial_step_) * Time::RESOLUTION +
           initial_time_.get_total_seconds();
}


Time SimulationManager::get_time() const
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    Time t0    = Time(initial_time_);
    t0.update(step_[omp_id] - initial_step_);
    return t0;
}


//! Get the timestep in second
double SimulationManager::get_resolution() const { return Time::RESOLUTION; }


size_t SimulationManager::get_current_step() const
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    return step_.at(omp_id);
}


double SimulationManager::get_current_substep() const
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    return substep_.at(omp_id);
}


void SimulationManager::set_max_resolution()
{
    max_resol_          = std::numeric_limits<double>::max();
    double max_modifier = std::numeric_limits<double>::max();

    // get max speed modifier from areas
    if (kernel().using_environment())
    {
        std::unordered_map<std::string, double> prop;
        for (auto name : kernel().space_manager.get_area_names())
        {
            kernel().space_manager.get_area_properties(name, prop);
            auto it = prop.find(names::speed_growth_cone);
            if (it != prop.end())
            {
                max_modifier = std::min(max_modifier, it->second);
            }
            else
            {
                max_modifier = std::min(max_modifier, 1.);
            }
        }
    }
    else
    {
        max_modifier = 1.;
    }

    // inverse to get max resolution from neuron properties
    max_modifier = 1. / max_modifier;

    max_resol_ = kernel().neuron_manager.get_max_resol() * max_modifier;
}

} // namespace growth
