/*
 * parallelism_manager.hpp
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

#ifndef PARALLEL_M_H
#define PARALLEL_M_H

#include "config.hpp"
#include "exceptions.hpp"
#include "manager_interface.hpp"

#ifdef WITH_MPI
#include <mpi.h>
#endif
#ifdef WITH_OMP
#include <omp.h>
#endif

namespace growth
{

class ParallelismManager : public ManagerInterface
{
  public:
    ParallelismManager();

    void mpi_init(int *argc, char **argv[]);
    void mpi_finalize();

    virtual void initialize();
    virtual void finalize();

    void set_status(const statusMap &status);
    void get_status(statusMap &status) const;

    /*
     * Get number of MPI processes
     */
    int get_num_mpi_proc() const;

    /*
     * Get number of virtual processes
     */
    stype get_num_virtual_processes() const;

    /*
     * Get rank of MPI process
     */
    int get_mpi_rank() const;

    /*
     * Return the process id for a given virtual process. The real process'
     * id of a virtual process is defined by the relation: t = (vp mod T),
     * where T is the total number of processes.
     */
    int get_thread_local_id() const;
    int get_thread_local_id(stype vp) const;

    int get_num_local_threads() const;
    void set_num_local_threads(int n_threads);

  private:
    // MPI
    int num_mpi_;
    int rank_;
    bool using_mpi_;
    int send_buffer_size_;
    int recv_buffer_size_;
#ifdef WITH_MPI
    MPI_Comm comm;
#endif
    // OpenMP
    int num_omp_;
    bool force_singlethread_;
};

// MPI implementations

inline int ParallelismManager::get_num_mpi_proc() const { return num_mpi_; }

inline int ParallelismManager::get_mpi_rank() const { return rank_; }

// OpenMP implementations

inline int ParallelismManager::get_thread_local_id() const
{
#ifdef WITH_OMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

inline int ParallelismManager::get_num_local_threads() const
{
    return num_omp_;
}

// Mixed implementations

inline int ParallelismManager::get_thread_local_id(stype vp) const
{
    return vp % num_mpi_;
}

inline stype ParallelismManager::get_num_virtual_processes() const
{
    return num_omp_ * num_mpi_;
}
} // namespace growth

#endif // PARALLEL_M_H
