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
    void new_branching_event(const Event &ev);
    bool simulating() const;

    Time get_time() const;
    const Time& get_initial_time() const;
    double get_resolution() const;
    double get_current_minutes() const;
    size_t get_current_step() const;
    double get_current_substep() const;

    void set_max_resolution();
    void push_max_resolution(int omp_id, double max_resol,
                             double old_max_resol);

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
    double final_substep_;
    Time::timeStep final_step_;
    Time initial_time_;
    Time final_time_;
    Time maximal_time_;
    bool terminate_;
    double resolution_scale_factor_;
    double max_resol_; // maximum allowed resolution
};


/**
 * Terminate the simulation after the time-slice is finished.
 */
inline void SimulationManager::terminate() { terminate_ = true; }

} // namespace growth

#endif /* SIMULATION_M_H */
