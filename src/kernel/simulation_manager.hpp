/*
 * simulation_manager.hpp
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
    const Time &get_initial_time() const;
    double get_resolution() const;
    double get_current_minutes() const;
    stype get_current_step() const;
    double get_current_substep() const;

    void set_max_resolution();
    void push_max_resolution(int omp_id, double max_resol,
                             double old_max_resol);

    void test_random_generator(Random_vecs &values, stype size);

  private:
    void initialize_simulation_(const Time &t);
    void finalize_simulation_();

    bool simulating_;
    bool print_time_;
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
