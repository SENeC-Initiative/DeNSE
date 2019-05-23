/*
 * rng_manager.hpp
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

#ifndef RNG_MANAGER_H
#define RNG_MANAGER_H

// C++ includes:
#include <memory>
#include <random>
#include <vector>

// Kernel includes
#include "config.hpp"
#include "manager_interface.hpp"
#include "parallelism_manager.hpp"

namespace growth
{

/**
 * Responsible for generating, seeding and passing the random number generators
 * to all the objects in the simulator.
 */
class RNGManager : public ManagerInterface
{
    friend class ParallelismManager;

  public:
    RNGManager();
    virtual ~RNGManager() {}

    virtual void initialize();
    virtual void finalize();

    virtual void get_status(statusMap &status) const;

    void seed(const std::vector<long> &seeds);

    mtPtr get_rng(int t);
    // global RNG: generate same random numbers on all threads (needed for MPI)
    //~ std::mt19937 get_grng();

  private:
    void create_rngs_();
    void create_grng_();

    /**
     * Vector of random number generators for threads.
     * There must be PRECISELY one rng per thread.
     */
    std::vector<mtPtr> rng_;

    /**
     * Global random number generator.
     * This rng must be synchronized on all threads
     */
    //~ std::mt19937 grng_;

    //! The seeds of the local RNGs. These do not necessarily describe the
    //! state of the RNGs.
    std::vector<long> rng_seeds_;

    //! The seed of the global RNG, not necessarily describing the
    //! state of the GRNG.
    //~ long grng_seed_;

}; // class RNGManager

} // namespace growth

inline std::shared_ptr<std::mt19937> growth::RNGManager::get_rng(int t)
{
    return rng_[t];
}

//~ inline std::mt19937& growth::RNGManager::get_grng()
//~ {
//~ return grng_;
//~ }

#endif /* RNG_MANAGER_H */
