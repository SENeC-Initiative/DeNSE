/*
 * simulation_manager.cpp
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

#include "simulation_manager.hpp"

// C includes:
#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>

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

/**
 * @brief compare the times of Events
 */
auto ev_greater = [](const Event &lhs, const Event &rhs) {
    return std::get<edata::TIME>(lhs) > std::get<edata::TIME>(rhs);
};


SimulationManager::SimulationManager()
    : simulating_(false) //!< true if simulation in progress
    , step_()            //!< Current step of the simulation
    , substep_()         //!< Precise time inside current step
    , final_substep_(0.) //!< Last substep, updated once per simulation call
    , final_step_(0)     //!< Last step, updated once per simulation call
    , initial_time_()    //!< Initial time (day, hour, min, sec)
    , final_time_()      //!< Final time (day, hour, min, sec)
    , maximal_time_()    //!< Maximal time (day, hour, min, sec)
    , terminate_(false)  //!< Terminate on signal or error
    , previous_resolution_(Time::RESOLUTION)
    , resolution_scale_factor_(1) //! rescale step size respct to old resolution
    , max_resol_(DEFAULT_MAX_RESOL)
{
}


void SimulationManager::initialize()
{
    // set resolution
    Time::reset_resolution();

    final_step_    = 0;
    final_substep_ = 0.;
    initial_time_  = Time();
    final_time_    = Time();
    max_resol_     = DEFAULT_MAX_RESOL;
}


void SimulationManager::finalize()
{
    step_.clear();
    substep_.clear();
}


/**
 * Reset the SimulationManager to the state at T = 0.
 */
void SimulationManager::reset_culture() {}


void SimulationManager::test_random_generator(Random_vecs &values, stype size)
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
            for (stype n = 0; n < size; n++)
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

    double old_substep = final_substep_;
    Time::to_steps(t, final_step_, final_substep_);

    // update
    final_substep_ += old_substep;
    if (final_substep_ >= 1)
    {
        final_step_++;
        final_substep_ = std::max(final_substep_ - 1., 0.);
    }

    // set the right number of step objects
    int num_omp = kernel().parallelism_manager.get_num_local_threads();
    step_       = std::vector<Time::timeStep>(num_omp, 0);
    substep_    = std::vector<double>(num_omp, 0.);

    // make sure branching events are cleared if we start from t = 0
    if (initial_time_ == Time())
    {
        branching_ev_.clear();
        branching_ev_tmp_.clear();
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
                neuron.second->initialize_next_event(rnd_engine);
            }

            // if time is zero, we need to initialize recorders
            if (initial_time_ == Time())
            {
                // record the first step
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

    simulating_ = true;
}


void SimulationManager::finalize_simulation_()
{
    for (stype s : step_)
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

    // finalize neurons
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
#ifndef NDEBUG
    initial_time_.update(final_step_, final_substep_);
    assert(std::abs(initial_time_.get_total_seconds() -
                    final_time_.get_total_seconds()) < 1e-4);
#endif

    initial_time_ = final_time_;

    // do not reset final_substep_ to zero!
    final_step_ = 0;
    simulating_ = false;
}


/**
 * @brief add a branching event
 *
 * Add an event into `branching_ev_` and sort the vector.
 */
void SimulationManager::new_branching_event(const Event &ev)
{
#pragma omp critical
    {
        branching_ev_tmp_.push_back(ev);
    }
}


/**
 * simulate for the given time .
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
        int omp_id = kernel().parallelism_manager.get_thread_local_id();
        stype current_step;
        Time time_next_ev;

        double previous_substep = 0.;
        bool new_step           = false;
        bool branching          = false;

        mtPtr rnd_engine = kernel().rng_manager.get_rng(omp_id);
        gidNeuronMap local_neurons =
            kernel().neuron_manager.get_local_neurons(omp_id);

        try
        {
            // then run the simulation
            while (step_[omp_id] < final_step_ or
                   (step_[omp_id] == final_step_ and
                    substep_[omp_id] < final_substep_))
            {
                current_step     = step_[omp_id];
                previous_substep = substep_[omp_id];

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

                // -------------- //
                // EVENT HANDLING //
                // -------------- //

                // check when the next event will occur and set step/substep
                if (branching_ev_.empty())
                {
                    new_step  = true;
                    branching = false;

                    if (current_step + 1 == final_step_)
                    {
                        substep_[omp_id] = final_substep_;
                    }
                    else
                    {
                        substep_[omp_id] = Time::RESOLUTION;
                    }
                }
                else
                {
                    Time next_time = initial_time_;
                    next_time.update(current_step + 1, 0);

                    time_next_ev = std::get<edata::TIME>(branching_ev_.back());
                    branching    = false;
                    new_step     = false;

                    if (time_next_ev < next_time)
                    {
                        substep_[omp_id] = Time::RESOLUTION -
                                           (next_time.get_total_minutes() -
                                            time_next_ev.get_total_minutes());

                        if (current_step == final_step_ and
                            substep_[omp_id] > final_substep_)
                        {
                            substep_[omp_id] = final_substep_;
                        }
                        else
                        {
                            new_step  = (substep_[omp_id] == Time::RESOLUTION);
                            branching = true;
                        }
                    }
                    else if (current_step == final_step_)
                    {
                        substep_[omp_id] = final_substep_;
                    }
                    else
                    {
                        substep_[omp_id] = Time::RESOLUTION;
                        new_step         = true;
                    }
                }

                assert(substep_[omp_id] >= 0.);

                // update neurons
                for (auto &neuron : local_neurons)
                {
                    neuron.second->grow(rnd_engine, current_step,
                                        substep_[omp_id] - previous_substep);
                }

                if (branching)
                {
                    // someone has to branch
                    Event &ev           = branching_ev_.back();
                    stype gid_branching = std::get<edata::NEURON>(ev);
                    auto it             = local_neurons.find(gid_branching);

                    if (it != local_neurons.end())
                    {
                        bool branched = local_neurons[gid_branching]->branch(
                            rnd_engine, ev);

                        // tell recorder manager
                        if (branched)
                        {
                            kernel().record_manager.new_branching_event(ev);
                        }
                    }

                    // wait for everyone to check, then remove event
#pragma omp barrier
#pragma omp single
                    {
                        branching_ev_.pop_back();
                    }
                }

#pragma omp barrier
                // update the R-tree
                kernel().space_manager.update_rtree();
                // wait for the update
#pragma omp barrier

                if (new_step)
                {
                    // full step is completed, record
                    // reset substep_, increment step_
                    kernel().record_manager.record(omp_id);
                    substep_[omp_id] = 0.;
                    step_[omp_id]++;
                }

#ifndef NDEBUG
                if (omp_id == 0 and step_[0] % 50 == 0)
                {
                    printf("##simulation step is %lu \n", step_[omp_id]);
                    printf("##simulated minutes: %f \n", get_current_minutes());
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
        get_param(status, names::max_allowed_resolution, max_resol_);
        double resolution;
        get_param(status, names::resolution, resolution);

        if (resolution > max_resol_)
        {
            throw std::invalid_argument(
                "`" + names::resolution + "` must be smaller or equal to `" +
                names::max_allowed_resolution + "`, which is currently " +
                std::to_string(max_resol_) + " minute.");
        }

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
    set_param(status, names::resolution, Time::RESOLUTION, "minute");
    set_param(status, names::max_allowed_resolution, max_resol_, "minute");

    // initial time is always up to date
    set_param(status, "second", initial_time_.get_sec(), "second");
    set_param(status, "minute", initial_time_.get_min(), "minute");
    set_param(status, "hour", initial_time_.get_hour(), "hour");
    set_param(status, "day", initial_time_.get_day(), "day");
}


/**
 * Gives current time in minutes.
 */
double SimulationManager::get_current_minutes() const
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    return (step_.at(omp_id)) * Time::RESOLUTION +
           initial_time_.get_total_minutes();
}


/**
 * Gives initial time.
 */
const Time &SimulationManager::get_initial_time() const
{
    return initial_time_;
}


Time SimulationManager::get_time() const
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    Time t0    = Time(initial_time_);
    //~ printf("initial time: %lu days %i hours %i minutes %f seconds\n",
    //~ t0.get_day(), t0.get_hour(), t0.get_min(), t0.get_sec());
    t0.update(step_[omp_id], substep_[omp_id]);
    //~ printf("final time: %lu days %i hours %i minutes %f seconds\n",
    //~ t0.get_day(), t0.get_hour(), t0.get_min(), t0.get_sec());
    return t0;
}


//! Get the timestep in minutes
double SimulationManager::get_resolution() const { return Time::RESOLUTION; }


stype SimulationManager::get_current_step() const
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


bool SimulationManager::simulating() const { return simulating_; }

} // namespace growth
