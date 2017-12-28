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
    , substep_()                  //!< Precise time inside current step
    , final_step_(0L)             //!< Last step, updated once per slice
    , initial_time_(Time())       //!< Initial time (day, hour, min, sec)
    , final_time_(Time())         //!< Final time (day, hour, min, sec)
    , maximal_time_(Time())       //!< Maximal time (day, hour, min, sec)
    , terminate_(false)           //!< Terminate on signal or error
    , previous_resolution_(Time::RESOLUTION)
    , resolution_scale_factor_(1) //! rescale step size respct to old resolution
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
    substep_.clear();
    step_ = std::vector<Time::timeStep>(num_omp, old_step);
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
    final_time_  = initial_time_ + t;
    size_t initial_step = final_step_;
    final_step_ += Time::to_steps(t);

    // set the right number of step objects
    int num_omp = kernel().parallelism_manager.get_num_local_threads();
    step_       = std::vector<Time::timeStep>(num_omp, initial_step);
    substep_    = std::vector<double>(num_omp, 0.);

    // reset branching events
    branching_ev_.clear();
    branching_ev_tmp_.clear();

#pragma omp parallel
    {
        int omp_id = kernel().parallelism_manager.get_thread_local_id();

#ifndef NDEBUG
        printf("\n\n\n ######## Starting simulation ########## \n\n");
#endif

        mtPtr rnd_engine = kernel().rng_manager.get_rng(omp_id);
        gidNeuronMap local_neurons =
            kernel().neuron_manager.get_local_neurons(omp_id);

        // first, initialize neurons
        for (auto &neuron : local_neurons)
        {
            neuron.second->initialize_next_event(
                rnd_engine, resolution_scale_factor_, initial_step);
        }
    }
}


void SimulationManager::finalize_simulation_()
{
    for (size_t s : step_)
    {
        assert(s == final_step_);
    }
    //    assert(initial_time_.get_total_seconds() ==
    //    final_time_.get_total_seconds());
    initial_time_.update(step_.front());

    // finalize neurons#pragma omp parallel
#pragma omp parallel
    {
        int omp_id = kernel().parallelism_manager.get_thread_local_id();

        gidNeuronMap local_neurons =
            kernel().neuron_manager.get_local_neurons(omp_id);

        for (auto &neuron : local_neurons)
        {
            neuron.second->finalize();
        }
    }

    // finalize recorder times
    kernel().record_manager.finalize_simulation(final_step_);
}


/**
 * @brief add a branching event
 *
 * Add an event into `branching_ev_` and sort the vector.
 */
void SimulationManager::new_branching_event(const Event& ev)
{
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

#pragma omp parallel
    {
        int omp_id = kernel().parallelism_manager.get_thread_local_id();
        size_t step_next_ev, current_step;
        double previous_substep = 0.;
        bool new_step = false;
        bool branched = false;

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
                if (! branching_ev_tmp_.empty())
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
                new_step = (substep_[omp_id] == Time::RESOLUTION);
            }
            else
            {
                substep_[omp_id] = Time::RESOLUTION;
                new_step = true;
            }

            assert(substep_[omp_id] >= 0.);

            // update neurons
            for (auto &neuron : local_neurons)
            {
                neuron.second->grow(
                    rnd_engine, current_step,
                    substep_[omp_id] - previous_substep);
            }

            if (step_next_ev == current_step)
            {
                // someone has to branch
                Event& ev = branching_ev_.back();
                size_t gid_branching = std::get<2>(ev);
                branched = local_neurons[gid_branching]->branch(rnd_engine, ev);
                // tell recorder manager
                if (branched)
                {
                    kernel().record_manager.new_branching_event(ev);
                }
                //~ std::string nname = std::get<3>(ev);
                //~ printf("Branching event ar step: %lu and substep %f\n"
                       //~ "Neurite %s of neuron %lu is branching at %lu:%f\n",
                       //~ current_step, substep_[omp_id], nname.c_str(),
                       //~ gid_branching, std::get<0>(ev), std::get<1>(ev));
                branching_ev_.pop_back();
            }

            if (new_step)
            {
                // full step is completed, record, reset substep_, incr. step_
                kernel().record_manager.record(omp_id);
                substep_[omp_id] = 0.;
                step_[omp_id]++;
            }

#ifndef NDEBUG
            if (step_[omp_id] % 50 == 0)
            {
                printf("##simulation step is %lu \n", step_[omp_id]);
                printf("##simulation seconds is %f \n", get_current_seconds());
            }
#endif /* NDEBUG */
        }
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
        }
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
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    Time t0 = Time(initial_time_);
    Time t = Time(substep_[omp_id], 0, 0, 0);
    t0.update(step_[omp_id]);
    return t0 + t;
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

} // namespace
