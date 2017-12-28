#ifndef SIMULATION_M_H
#define SIMULATION_M_H

#include <random>
#include <vector>

#include "config.hpp"
#include "elements_types.hpp"
#include "growth_time.hpp"
#include "manager_interface.hpp"

namespace growth
{
class Neuron;

class SimulationManager : public ManagerInterface
{
  public:
    SimulationManager();

    // Init/finalize functions
    virtual void initialize();
    virtual void finalize();

    // Simulation functions
    void simulate(const Time &t);
    void terminate();
    void reset_culture(); // @TODO

    // getter/setter functions
    virtual void set_status(const statusMap &);
    virtual void get_status(statusMap &) const;
    void num_threads_changed(int num_omp);
    void new_branching_event(const Event& ev);

    Time get_time() const;
    double get_resolution() const;
    double get_current_seconds() const;
    size_t get_current_step() const;
    double get_current_substep() const;
    void test_random_generator(Random_vecs &values, size_t size);

  private:
    void initialize_simulation_(const Time &t);
    void finalize_simulation_();

    bool simulating_;
    double previous_resolution_;
    std::vector<Time::timeStep> step_;
    std::vector<double> substep_;
    std::vector<Event> branching_ev_;
    std::vector<Event> branching_ev_tmp_;
    Time::timeStep final_step_;
    Time initial_time_;
    Time final_time_;
    Time maximal_time_;
    bool terminate_;
    double resolution_scale_factor_;
};


/**
 * @brief compare the times of Events
 */
auto ev_greater = [](const Event& lhs, const Event& rhs)
{
    return std::tie(std::get<0>(lhs), std::get<1>(lhs))
           > std::tie(std::get<0>(rhs), std::get<1>(rhs));
};


/**
 * Terminate the simulation after the time-slice is finished.
 */
inline void SimulationManager::terminate() { terminate_ = true; }


//~ inline Time const&
//~ SimulationManager::get_initial_time() const
//~ {
//~ return initial_time_;
//~ }

//~ inline Time const
//~ SimulationManager::get_time() const
//~ {
//~ assert( not simulating_ );
//~ return initial_time_ + Time::time( step_ );
//~ }

//~ inline size_t
//~ SimulationManager::get_slice() const
//~ {
//~ return slice_;
//~ }

//~ inline Time const&
//~ SimulationManager::get_clock() const
//~ {
//~ return clock_;
//~ }
}

#endif /* SIMULATION_M_H */
