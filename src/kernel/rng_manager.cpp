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
    for (size_t i = 0; i < seeds.size(); i++)
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
    if (!rng_.empty())
    {
        rng_.clear();
    }

    // initialize
    size_t mpi_rank = kernel().parallelism_manager.get_mpi_rank();
    size_t num_omp  = kernel().parallelism_manager.get_num_local_threads();
    rng_seeds_.resize(num_omp);

    for (size_t i = 0; i < num_omp; i++)
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
