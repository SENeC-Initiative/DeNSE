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

    size_t create_recorder(const std::vector<statusMap> &obj_params);

    void record(int omp_id);
    void finalize_simulation(size_t steps);

    void get_status(statusMap &status) const;
    void get_defaults(statusMap &status) const;
    bool is_recorder(size_t gid) const;
    size_t num_recorders() const;
    void set_status(const statusMap &status);
    void num_threads_changed(int num_omp);
    void new_branching_event(const Event &ev);

    statusMap get_recorder_status(size_t gid) const;
    void get_recorder_type(size_t gid, std::string &level,
                           std::string &event_type) const;
    bool get_next_recording(size_t gid, std::vector<Property> &ids,
                            std::vector<double> &values);
    bool get_next_time(size_t gid, std::vector<Property> &ids,
                       std::vector<double> &values,
                       const std::string &time_units);

  private:
    std::unordered_map<size_t, std::shared_ptr<BaseRecorder>> c_recorders_;
    std::unordered_map<size_t, std::shared_ptr<BaseRecorder>> d_recorders_;
    std::unordered_map<size_t, std::vector<size_t>> neuron_to_d_recorder_;
    std::unordered_map<size_t, std::vector<size_t>> neuron_to_c_recorder_;
    std::vector<std::vector<size_t>> omp_id_crec_;
    std::vector<std::vector<size_t>> omp_id_drec_;
    size_t last_omp_id_;
    int num_threads_;
    std::vector<Event> events_;
};

} /* namespace */

#endif
