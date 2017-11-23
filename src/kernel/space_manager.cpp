#include "space_manager.hpp"

// C++ includes
#include <cmath>

// kernel include
#include "kernel_manager.hpp"

// lib include
#include "config_impl.hpp"

namespace growth
{

void notice(const char *fmt, ...)
{
    //~ va_list ap;

    //~ printf( stdout, "NOTICE: ");

    //~ va_start (ap, fmt);
    //~ vfprintf( stdout, fmt, ap);
    //~ va_end(ap);
    //~ printf( stdout, "\n" );
}


void log_and_exit(const char *fmt, ...)
{
    //~ va_list ap;

    //~ printf( stdout, "ERROR: ");

    //~ va_start (ap, fmt);
    //~ vfprintf( stdout, fmt, ap);
    //~ va_end(ap);
    //~ printf( stdout, "\n" );
    //~ exit(1);
}


SpaceManager::SpaceManager()
{
    environment_initialized_ = false;
    environment_manager_     = nullptr;
    context_handler_         = nullptr;
}


void SpaceManager::initialize()
{
    initGEOS_r(notice, log_and_exit);
    context_handler_ = GEOS_init_r();
}


void SpaceManager::finalize() { finishGEOS_r(context_handler_); }


GeomPtr SpaceManager::geosline_from_points(const Point &pointA,
                                           const Point &pointB) const
{
    GEOSCoordSeq sq = GEOSCoordSeq_create_r(context_handler_, 2, 2);
    GEOSCoordSeq_setX_r(context_handler_, sq, 0, pointA.at(0));
    GEOSCoordSeq_setY_r(context_handler_, sq, 0, pointA.at(1));
    GEOSCoordSeq_setX_r(context_handler_, sq, 1, pointB.at(0));
    GEOSCoordSeq_setY_r(context_handler_, sq, 1, pointB.at(1));
    GeomPtr line = GeomPtr(GEOSGeom_createLineString_r(context_handler_, sq));
    return line;
}


GeomPtr SpaceManager::geospoint_from_point(const Point &point) const
{
    GEOSCoordSeq sq = GEOSCoordSeq_create_r(context_handler_, 1, 2);
    GEOSCoordSeq_setX_r(context_handler_, sq, 0, point.at(0));
    GEOSCoordSeq_setY_r(context_handler_, sq, 0, point.at(1));
    /*    if (dim==3)*/
    /*{GEOSCoordSeq_setZ( sequence, 2, point.at(2));}*/
    GeomPtr geospoint = GeomPtr(GEOSGeom_createPoint_r(context_handler_, sq));
    assert(GEOSisValid_r(context_handler_, geospoint.get()));
    return geospoint;
}


bool SpaceManager::sense(std::vector<double> &directions_weigths,
                         const Filopodia &filopodia, const Point &position,
                         const Move &move, const double distance,
                         const double in_value, const double out_value)
{
    //~ assert(env_contains(Point(position.at(0), position.at(1))));
    for (int n_angle = 0; n_angle < filopodia.size; n_angle++)
    {
        if (not std::isnan(directions_weigths[n_angle]))
        {
            double angle =
                move.sigma_angle * filopodia.directions[n_angle] + move.angle;
            Point target_pos = Point(position.at(0) + cos(angle) * distance,
                                     position.at(1) + sin(angle) * distance);
            if (env_intersect(geosline_from_points(position, target_pos)))
            {
                directions_weigths[n_angle] *= out_value;
            }
            else
            {
                directions_weigths[n_angle] *= in_value;
            }
        }
    }
}


bool SpaceManager::env_intersect(GeomPtr line) const
{
    if (environment_initialized_ == false)
    {
        return 0;
    }
    int omp = kernel().parallelism_manager.get_thread_local_id();
    return GEOSPreparedIntersects_r(
        context_handler_, environment_manager_->get_border(omp), line.get());
}


bool SpaceManager::env_contains(const Point &point) const
{
    if (environment_initialized_ == false)
    {
        return 1;
    }
    int omp = kernel().parallelism_manager.get_thread_local_id();
    GeomPtr geospoint = geospoint_from_point(point);
    return GEOSPreparedContains_r(context_handler_,
                                  environment_manager_->get_prepared(omp),
                                  geospoint.get());
}


bool SpaceManager::has_environment() const { return environment_initialized_; }


void SpaceManager::set_environment(GEOSGeom environment)
{
    const GEOSGeom env = static_cast<GEOSGeom>(environment);
    assert(GEOSisValid_r(context_handler_, env));
    environment_initialized_ = true;
    //@TODO  make a unique pointer<Right>
    environment_manager_ =
        std::unique_ptr<Environment>(new Environment(env, context_handler_));
}


void SpaceManager::get_environment(GEOSGeom &environment) const
{
    if (environment_manager_ != nullptr)
    {
        environment = environment_manager_->get_environment();
    }
}


void SpaceManager::init_spatial_grid(int local_num_threads)
{
    for (int i = 0; i < local_num_threads; i++)
    {
    }
    // spatial_grid.push_back(Space::Shape());
}


int SpaceManager::get_region_thread(const Point &position)
{
    return get_region_thread(position.at(0), position.at(1));
}


int SpaceManager::get_region_thread(double x, double y)
{
    // @TODO
    return 0;
}


void SpaceManager::set_status(const statusMap &config)
{
    /*    int dim = Space::DEFAULT_DIM;*/
    // get_param(config, "dimension", dim);
    /*Space::set_dimension_(dim);*/
}


void SpaceManager::get_status(statusMap &status) const
{
    set_param(status, "environment_initialized", environment_initialized_);
    // set_param(status, "dimension", Space::get_dimension());
}


Environment::Environment(GEOSGeom environment,
                         GEOSContextHandle_t context_handler)
    : environment_(environment)
{

    assert(0 != environment_);
    for (int i=0; i<kernel().parallelism_manager.get_num_local_threads(); i++)
    {
        prepared_env_.push_back(GEOSPrepare_r(context_handler, environment_));
        const GEOSGeom border = GEOSBoundary_r(context_handler, environment_);
        prepared_border_.push_back(GEOSPrepare_r(context_handler, border));
        assert(prepared_env_[i] != 0);
        assert(prepared_border_[i] != 0);
    }
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

} // namespace growth
