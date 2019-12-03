/*
 * rng_manager.cpp
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

// Includes from kernel
#include "config_impl.hpp"
#include "kernel_manager.hpp"

#include "rng_manager.hpp"

growth::RNGManager::RNGManager()
    : rng_()
    , rng_seeds_()
{
}

void growth::RNGManager::initialize()
{
    create_rngs_();
    //~ create_grng_();
}

void growth::RNGManager::finalize() {}

void growth::RNGManager::seed(const std::vector<long> &seeds)
{

    rng_.resize(seeds.size());
    rng_seeds_.resize(seeds.size());

    for (stype i = 0; i < seeds.size(); i++)
    {
#ifndef NDEBUG
        printf(" seeding the random generator\n");
#endif
        long seed     = seeds[i];
        rng_[i]       = std::make_shared<std::mt19937>(seed);
        rng_seeds_[i] = seed;
    }
}

void growth::RNGManager::get_status(statusMap &status) const
{
    set_param(status, "seeds", rng_seeds_, "");
}

void growth::RNGManager::create_rngs_()
{
#ifndef NDEBUG
    printf("creating the random generators\n");
#endif

    if (!rng_.empty())
    {
        rng_.clear();
    }

    // initialize
    stype mpi_rank = kernel().parallelism_manager.get_mpi_rank();
    stype num_omp  = kernel().parallelism_manager.get_num_local_threads();

    rng_seeds_.resize(num_omp);

    for (stype i = 0; i < num_omp; i++)
    {
        long s = mpi_rank + i;
        /*
         * We have to ensure that each thread is provided with a different
         * stream of random numbers.  The seeding method for Knuth's LFG
         * generator guarantees that different seeds yield non-overlapping
         * random number sequences.
         *
         * We therefore have to seed with known numbers: using random
         * seeds here would run the risk of using the same seed twice.
         * For simplicity, we use 1 .. n_vps.
         */
        rng_seeds_[i] = s;
        rng_.push_back(std::make_shared<std::mt19937>(s));
    }
}
