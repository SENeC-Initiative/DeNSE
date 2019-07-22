/*
 * parallelism_manager.cpp
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

// C++ includes
#include <set>

// Kernel includes
#include "config_impl.hpp"
#include "kernel_manager.hpp"

#include "parallelism_manager.hpp"

namespace growth
{

ParallelismManager::ParallelismManager()
    : num_mpi_(1)
    , rank_(0)
    , send_buffer_size_(1)
    , recv_buffer_size_(1)
#ifdef WITH_MPI
    , comm(0)
#endif
    , num_omp_(1)
#ifdef WITH_OMP
    , force_singlethread_(false)
#else
    , force_singlethread_(true)
#endif
{
}


void ParallelismManager::mpi_init(int *argc, char **argv[])
{
#ifdef WITH_MPI
    int init;
    MPI_Initialized(&init);

    if (init == 0)
    {
        int provided_thread_level;
        MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED,
                        &provided_thread_level);
        comm = MPI_COMM_WORLD;
    }

    MPI_Comm_size(comm, &num_mpi_);
    MPI_Comm_rank(comm, &rank_);

    recv_buffer_size_ = send_buffer_size_ * num_mpi_;
#endif /* #ifdef HAVE_MPI */
}


/**
 * Finish off MPI routines
 */
void ParallelismManager::mpi_finalize()
{
#ifdef WITH_MPI
    int finalized;
    MPI_Finalized(&finalized);

    int initialized;
    MPI_Initialized(&initialized);
#endif /* #ifdef HAVE_MPI */
}


void ParallelismManager::initialize()
{
#ifdef WITH_OMP
    // Keeps OpenMP from automagically changing the number of threads used for
    // parallel regions.
    omp_set_dynamic(false);
#endif
    set_num_local_threads(1);
}


void ParallelismManager::finalize() {}


void ParallelismManager::set_num_local_threads(int n_threads)
{
#ifdef WITH_OMP
    // call first space_manager, then simulation manager, then neuron_manager,
    // then record_manager
    kernel().space_manager.num_threads_changed(n_threads);
    kernel().simulation_manager.num_threads_changed(n_threads);
    kernel().neuron_manager.init_neurons_on_thread(n_threads);
    kernel().record_manager.num_threads_changed(n_threads);
    omp_set_num_threads(n_threads);
#endif
    num_omp_ = n_threads;
}


void ParallelismManager::set_status(const statusMap &status)
{

    int num_omp = num_omp_;
    get_param(status, "num_local_threads", num_omp);
    bool num_omp_changed = (num_omp != num_omp_);

    if (num_omp_changed)
    {
        if (kernel().get_num_objects())
        {
            throw InvalidParameter("Cannot change the number of threads after "
                                   "objects have been created.",
                                   __FUNCTION__, __FILE__, __LINE__);
        }
        else if (num_omp != 1 && force_singlethread_)
        {
            throw InvalidParameter("Multithreading not supported.",
                                   __FUNCTION__, __FILE__, __LINE__);
        }
        else
        {
            set_num_local_threads(num_omp);
        }
    }

    // Updates in RNGManager MUST occur after set_num_local_threads
    std::vector<long> seeds = {};
    get_param(status, "seeds", seeds);
    //~ seeds.clear();

    if (!seeds.empty())
    {
        if (seeds.size() != num_mpi_ * num_omp_)
        {
            throw InvalidParameter(
                "Number of seeds must equal the number of virtual processes: "
                "expected length " +
                    std::to_string(num_mpi_ * num_omp_) + " but received " +
                    std::to_string(seeds.size()) + ".",
                __FUNCTION__, __FILE__, __LINE__);
        }
        // check if seeds are unique
        std::set<long> seedset;
        for (stype i = 0; i < seeds.size(); i++)
        {
            long s = seeds[i];
            if (!seedset.insert(s).second)
            {
                throw InvalidParameter("Seeds are not unique across threads!",
                                       __FUNCTION__, __FILE__, __LINE__);
            }
        }
        // seed RNGs [mpi_id*num_omp_, (mpi_id+1)*num_omp_)
        stype start = get_mpi_rank() * num_omp_;
        stype stop  = (get_mpi_rank() + 1) * num_omp_;
        std::vector<long> local_seeds(stop - start);
        for (stype i = start; i < stop; i++)
        {
            local_seeds[i - start] = seeds[i];
        }
        kernel().rng_manager.seed(local_seeds);
    }
    else if (num_omp_changed)
    {
        kernel().rng_manager.create_rngs_();
    }
}


void ParallelismManager::get_status(statusMap &status) const
{
    set_param(status, "num_mpi_processes", num_mpi_, "");
    set_param(status, "num_local_threads", num_omp_, "");
    set_param(status, "num_virtual_processes", num_omp_ * num_mpi_, "");
}

} // namespace growth
