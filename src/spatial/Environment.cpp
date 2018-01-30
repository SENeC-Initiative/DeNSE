#include "Environment.hpp"

#include "kernel_manager.hpp"


namespace growth
{

Environment::Environment(GEOSGeom environment,
                         GEOSContextHandle_t context_handler)
    : environment_(environment)
{
    assert(0 != environment_);
    assert(GEOSisValid_r(context_handler, environment_));

    for (int i = 0; i < kernel().parallelism_manager.get_num_local_threads();
         i++)
    {
        prepared_env_.push_back(GEOSPrepare_r(context_handler, environment_));
        const GEOSGeom border = GEOSBoundary_r(context_handler, environment_);
        prepared_border_.push_back(GEOSPrepare_r(context_handler, border));
        assert(prepared_env_[i] != 0);
        assert(prepared_border_[i] != 0);
    }
}


Environment::~Environment()
{
    for (const GEOSPreparedGeometry *env : prepared_env_)
    {
        delete env;
    }
    prepared_env_.clear();

    for (const GEOSPreparedGeometry *border : prepared_border_)
    {
        delete border;
    }
    prepared_border_.clear();
}


GEOSGeom Environment::get_environment() const { return environment_; }


const GEOSPreparedGeometry *Environment::get_prepared(int omp_id) const
{
    return prepared_env_[omp_id];
}


const GEOSPreparedGeometry *Environment::get_border(int omp_id) const
{
    return prepared_border_[omp_id];
}

} /* namespace */
