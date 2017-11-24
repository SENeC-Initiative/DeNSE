#ifndef RECORD_M_H
#define RECORD_M_H

// C++ includes
#include <vector>
#include <fstream>
#include <iostream>

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

    bool is_recorder(size_t gid) const;
    void get_status(statusMap &status) const;
    size_t num_recorders() const;
    void set_status(const statusMap &status);
    void num_threads_changed(int num_omp);

    statusMap get_recorder_status(size_t gid) const;
    void get_recorder_type(size_t gid, std::string& level,
                           std::string& event_type) const;
    bool get_next_recording(size_t gid, std::vector<Property>& ids,
                            std::vector<double>& values);
    bool get_next_time(size_t gid, std::vector<Property>& ids,
                       std::vector<double>& values);

  private:
    std::unordered_map< size_t, std::shared_ptr<BaseRecorder> > recorders_;
    std::vector< std::vector<size_t> > omp_id_rec_;
    size_t last_omp_id_;
    int num_threads_;
};

} /* namespace */

#endif