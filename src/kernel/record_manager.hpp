/*
 * record_manager.hpp
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

#ifndef RECORD_M_H
#define RECORD_M_H

// C++ includes
#include <fstream>
#include <iostream>
#include <vector>

// kernel includes
#include "manager_interface.hpp"
#include "recorders.hpp"

// libgrowth includes
#include "config.hpp"
#include "elements_types.hpp"


namespace growth
{

class RecordManager : public ManagerInterface
{
  public:
    RecordManager();

    virtual void initialize();
    virtual void finalize();

    stype create_recorder(const std::vector<statusMap> &obj_params);

    void record(int omp_id);
    void finalize_simulation(stype steps);

    void get_status(statusMap &status) const;
    void get_defaults(statusMap &status) const;
    bool is_recorder(stype gid) const;
    stype num_recorders() const;
    void set_status(const statusMap &status);
    void num_threads_changed(int num_omp);
    void new_branching_event(const Event &ev);

    void neurons_deleted(const std::vector<stype> &gids);
    void new_neurite(stype neuron, const std::string &neurite);
    void gc_died(stype neuron, const std::string &neurite, stype gc_id);

    statusMap get_recorder_status(stype gid) const;
    void get_recorder_type(stype gid, std::string &level,
                           std::string &event_type) const;
    bool get_next_recording(stype gid, std::vector<Property> &ids,
                            std::vector<double> &values);
    bool get_next_time(stype gid, std::vector<Property> &ids,
                       std::vector<double> &values,
                       const std::string &time_units);

  private:
    std::unordered_map<stype, std::shared_ptr<BaseRecorder>> c_recorders_;
    std::unordered_map<stype, std::shared_ptr<BaseRecorder>> d_recorders_;
    std::unordered_map<stype, std::vector<stype>> neuron_to_d_recorder_;
    std::unordered_map<stype, std::vector<stype>> neuron_to_c_recorder_;
    std::vector<std::vector<stype>> omp_id_crec_;
    std::vector<std::vector<stype>> omp_id_drec_;
    stype last_omp_id_;
    int num_threads_;
    std::vector<Event> events_;
};

} // namespace growth

#endif
