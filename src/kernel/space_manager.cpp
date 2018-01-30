#include "space_manager.hpp"

// C++ includes
#include <cmath>

// GEOS include
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/MultiLineString.h>

// kernel include
#include "kernel_manager.hpp"

// lib include
#include "config_impl.hpp"
#include "growth_names.hpp"

// spatial include
#include "Area.hpp"
#include "Environment.hpp"
#include "Wall.hpp"


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
  : environment_initialized_(false)
  , environment_manager_(nullptr)
  , context_handler_(nullptr)
  , wall_area_width_(5.)
{
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
    return geospoint;
}


bool SpaceManager::sense(std::vector<double> &directions_weigths,
                         std::vector<std::string> &new_pos_area,
                         const Filopodia &filopodia, const Point &position,
                         const Move &move, double distance, double substep,
                         double in_value, double out_value,
                         const std::string &area, double proba_down_move,
                         double max_height_up_move)
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();

    assert(env_contains(Point(position.at(0), position.at(1)), omp_id));

    for (int n_angle = 0; n_angle < filopodia.size; n_angle++)
    {
        if (not std::isnan(directions_weigths[n_angle]))
        {
            double angle = move.angle + move.sigma_angle * sqrt(substep) *
                                            filopodia.directions[n_angle];

            Point target_pos = Point(position.at(0) + cos(angle) * distance,
                                     position.at(1) + sin(angle) * distance);
            GeomPtr geosline = geosline_from_points(position, target_pos);


            if (env_intersect(geosline, omp_id))
            {
                // intersection with environment border?
                directions_weigths[n_angle] *= out_value;
            }
            else if (area_intersect(area, geosline, omp_id))
            {
                // check the number of intersections
                AreaPtr new_area;
                AreaPtr old_area = areas_[area];
                std::vector<Point> middle_points;
                double old_height, new_height, new_substrate_affinity, delta_h;

                size_t n_intersect = area_num_intersections(
                    area, geosline, target_pos, middle_points, omp_id);

                // intersection with area?
                old_height = old_area->get_height();

                for (unsigned int i = 0; i < n_intersect; i++)
                {
                    std::string name_new_area =
                        get_containing_area(middle_points[i], omp_id);
                    new_pos_area[n_angle] = name_new_area;
                    new_area              = areas_[name_new_area];
                    new_height            = new_area->get_height();
                    new_substrate_affinity =
                        new_area->get_property(names::substrate_affinity);
                    delta_h = abs(new_height - old_height);

                    if (old_height == new_height)
                    {
                        // same height => proba from substrate affinity
                        if (n_intersect == 1)
                        {
                            directions_weigths[n_angle] *=
                                new_substrate_affinity;
                        }
                    }
                    else if (old_height > new_height)
                    {
                        // growth cone may go down:
                        // - if filopodia can touch the bottom, the probability
                        // of
                        //   making a down move is increased or decreased
                        //   depending on the substrate affinity
                        // - if the bottom is too far, only take the probability
                        //   of making a down move
                        double h = filopodia.finger_length - delta_h;
                        double delta_proba =
                            h > 0 ? (new_substrate_affinity - 1) * h /
                                        filopodia.finger_length
                                  : 0;

                        if ((proba_down_move + delta_proba) <= 0)
                        {
                            directions_weigths[n_angle] = nan("");
                        }
                        else
                        {
                            directions_weigths[n_angle] *=
                                substep * (proba_down_move + delta_proba);
                        }
                    }
                    else
                    {
                        // growth cone may go upwards:
                        // - move is possible if the step is smaller than
                        //   `max_height_up_move`, in [0,
                        //   new_substrate_affinity]
                        // - if step is higher, move is impossible
                        directions_weigths[n_angle] *=
                            substep * new_substrate_affinity *
                            std::exp(-delta_h / max_height_up_move);
                    }
                }
            }
            else
            {
                //~ #ifndef NDEBUG
                //~ printf("on regular substrate\n");
                //~ #endif
                directions_weigths[n_angle] *= in_value;
            }

#ifndef NDEBUG
            if (directions_weigths[n_angle] < 0)
            {
                printf("negative value: %f\n", directions_weigths[n_angle]);
            }
#endif
        }
    }
    //~ #ifndef NDEBUG
    //~ printf("out\n");
    //~ #endif
}


bool SpaceManager::sense_walls(std::vector<double> &directions_weigths,
                               const Filopodia &filopodia,
                               const Point &position, const Move &move,
                               double distance, double substep,
                               double wall_affinity, const std::string &area)
{
    int omp = kernel().parallelism_manager.get_thread_local_id();

    assert(env_contains(Point(position.at(0), position.at(1)), omp));

    for (int n_angle = 0; n_angle < filopodia.size; n_angle++)
    {
        if (not std::isnan(directions_weigths[n_angle]))
        {
            double angle = move.angle + move.sigma_angle * sqrt(substep) *
                                            filopodia.directions[n_angle];

            Point target_pos = Point(position.at(0) + cos(angle) * distance,
                                     position.at(1) + sin(angle) * distance);
            GeomPtr geosline = geosline_from_points(position, target_pos);

            if (wall_area_contains(area, target_pos, omp))
            {
                directions_weigths[n_angle] *= wall_affinity;
            }
        }
    }
}


bool SpaceManager::env_intersect(GeomPtr line, int omp_id) const
{
    if (environment_initialized_ == false)
    {
        return 0;
    }

    return GEOSPreparedIntersects_r(
        context_handler_, environment_manager_->get_border(omp_id), line.get());
}


bool SpaceManager::env_contains(const Point &point, int omp_id) const
{
    if (environment_initialized_ == false)
    {
        return 1;
    }

    GeomPtr geospoint = geospoint_from_point(point);

    return GEOSPreparedContains_r(context_handler_,
                                  environment_manager_->get_prepared(omp_id),
                                  geospoint.get());
}


bool SpaceManager::wall_area_contains(const std::string &area,
                                      const Point &point, int omp_id) const
{
    if (environment_initialized_ == false)
    {
        return 1;
    }

    GeomPtr geospoint = geospoint_from_point(point);

    return GEOSPreparedContains_r(context_handler_,
                                  walls_.at(area)->get_wall_area(omp_id),
                                  geospoint.get());
}


bool SpaceManager::area_intersect(const std::string &area, GeomPtr line,
                                  int omp_id) const
{
    if (environment_initialized_ == false)
    {
        return 0;
    }

    AreaPtr a = areas_.at(area);

    return GEOSPreparedIntersects_r(context_handler_, a->get_border(omp_id),
                                    line.get());
}


size_t SpaceManager::area_num_intersections(const std::string &area,
                                            GeomPtr line, const Point &tgt,
                                            std::vector<Point> &points,
                                            int omp_id) const
{
    if (environment_initialized_ == false)
    {
        return 0;
    }

    AreaPtr a = areas_.at(area);

    //~ GeomPtr intersect = GeomPtr(GEOSIntersection_r(
    //~ context_handler_, a->get_shape(omp).get(), line.get()));
    GEOSGeometry *intersect = GEOSIntersection_r(
        context_handler_, a->get_shape(omp_id).get(), line.get());

    size_t num_intersect = GEOSGetNumGeometries_r(context_handler_, intersect);

    if (num_intersect > 1)
    {
        double x0, x1, y0, y1;

        for (unsigned int i = 0; i < num_intersect - 1; i++)
        {
            const GEOSGeometry *g0 =
                GEOSGetGeometryN_r(context_handler_, intersect, i);
            auto coords1 = GEOSGeom_getCoordSeq_r(context_handler_, g0);
            GEOSCoordSeq_getX_r(context_handler_, coords1, 0, &x0);
            GEOSCoordSeq_getY_r(context_handler_, coords1, 0, &y0);
            const GEOSGeometry *g1 =
                GEOSGetGeometryN_r(context_handler_, intersect, i + 1);
            auto coords2 = GEOSGeom_getCoordSeq_r(context_handler_, g1);
            GEOSCoordSeq_getX_r(context_handler_, coords2, 0, &x1);
            GEOSCoordSeq_getY_r(context_handler_, coords2, 0, &y1);

            points.push_back(Point(0.5 * (x0 + x1), 0.5 * (y0 + y1)));
        }
    }
    else
    {
        points.push_back(tgt);
    }

    GEOSGeom_destroy_r(context_handler_, intersect);

    return num_intersect;
}


bool SpaceManager::has_environment() const { return environment_initialized_; }


void SpaceManager::set_environment(
    GEOSGeom environment, const std::vector<GEOSGeom> &walls,
    const std::vector<GEOSGeom> &areas, const std::vector<double> &heights,
    const std::vector<std::string> &names,
    const std::vector<std::unordered_map<std::string, double>> &properties)
{
    const GEOSGeom env = static_cast<GEOSGeom>(environment);
    assert(GEOSisValid_r(context_handler_, env));
    environment_initialized_ = true;

    environment_manager_ =
        std::unique_ptr<Environment>(new Environment(env, context_handler_));

    for (size_t i = 0; i < areas.size(); i++)
    {
        walls_[names[i]] = std::make_shared<Wall>(walls[i], context_handler_);
        areas_[names[i]] = std::make_shared<Area>(
            areas[i], context_handler_, heights[i], names[i], properties[i]);
    }
}


void SpaceManager::get_environment(
    GEOSGeom &environment, std::vector<GEOSGeom> &areas,
    std::vector<double> &heights, std::vector<std::string> &names,
    std::vector<std::unordered_map<std::string, double>> &properties) const
{
    if (environment_manager_ != nullptr)
    {
        environment = environment_manager_->get_environment();
    }

    int omp_id = kernel().parallelism_manager.get_thread_local_id();

    for (auto &a : areas_)
    {
        areas.push_back(a.second->get_shape(omp_id).get());
        heights.push_back(a.second->get_height());
        names.push_back(a.second->get_name());

        std::unordered_map<std::string, double> prop;
        a.second->get_properties(prop);
        properties.push_back(prop);
    }
}


void SpaceManager::init_spatial_grid(int local_num_threads)
{
    for (int i = 0; i < local_num_threads; i++)
    {
    }
    // spatial_grid.push_back(Space::Shape());
}


int SpaceManager::get_region_thread(const Point &position) const
{
    return get_region_thread(position.at(0), position.at(1));
}


int SpaceManager::get_region_thread(double x, double y) const
{
    // @TODO
    return 0;
}


std::string SpaceManager::get_containing_area(const Point &position,
                                              int omp_id) const
{
    if (environment_initialized_ == false)
    {
        return "";
    }

    GeomPtr geospoint = geospoint_from_point(position);
    bool contained    = false;

    for (auto &area : areas_)
    {
        contained = GEOSPreparedContains_r(
            context_handler_, area.second->get_area(omp_id), geospoint.get());

        if (contained)
        {
            return area.first;
        }
    }

    return "default_area";
}


AreaPtr SpaceManager::get_area(const std::string &name) const
{
    auto it = areas_.find(name);

    if (it == areas_.end())
    {
        return nullptr;
    }

    return areas_.at(name);
}


void SpaceManager::set_status(const statusMap &config)
{
    get_param(config, "wall_area_width", wall_area_width_);
    /*    int dim = Space::DEFAULT_DIM;*/
    // get_param(config, "dimension", dim);
    /*Space::set_dimension_(dim);*/
}


void SpaceManager::get_status(statusMap &status) const
{
    set_param(status, "environment_initialized", environment_initialized_);
    set_param(status, "wall_area_width", wall_area_width_);
    // set_param(status, "dimension", Space::get_dimension());
}

} // namespace growth
