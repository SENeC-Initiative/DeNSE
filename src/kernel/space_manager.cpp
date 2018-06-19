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


bool SpaceManager::sense(std::vector<double> &directions_weights,
                         std::vector<bool> &wall_presence,
                         const Filopodia &filopodia, const Point &position,
                         const Move &move, const std::string &area,
                         double proba_down_move, double max_height_up_move,
                         double substep, double sqrt_resol,
                         unsigned int delta_filo)
{
    // Useful values
    int omp_id        = kernel().parallelism_manager.get_thread_local_id();
    double subs_afty  = filopodia.substrate_affinity;
    double len_filo   = filopodia.finger_length;
    double wall_afty  = filopodia.wall_affinity;
    AreaPtr old_area  = areas_[area];
    double old_height = old_area->get_height();

    //~ assert(env_contains(Point(position.at(0), position.at(1)), omp_id)
    //~ && printf("omp (%i) pos (%f - %f)\n", omp_id, position.at(0),
    //position.at(1)));

#ifndef NDEBUG
    if (not env_contains(Point(position.at(0), position.at(1)), omp_id))
    {
        printf("omp (%i) pos (%f - %f)\n", omp_id, position.at(0),
               position.at(1));
    }
#endif

    // compute the number of filopodia to ignore
    //~ unsigned int ignore = 0.5*((sqrt_resol - sqrt(substep)) / (sqrt_resol -
    //1) * delta_filo); ~ unsigned int n_max  = filopodia.size - ignore; ~
    //unsigned int n_min  = ignore;

    // values used locally inside loop
    double angle, new_height, new_substrate_affinity, distance;
    bool filo_wall, interacting(false);
    Point filo_pos, lamel_pos, tmp_pos;
    std::vector<Point> middle_points;
    double delta_h, h, delta_proba;
    GeomPtr filo_line, lamel_line;
    std::string name_new_area;
    unsigned int i, n_angle;
    int last_filo_wall(-1);
    size_t n_intersect;
    AreaPtr new_area;

    // test the environment for each of the filopodia's angles
    //~ for (n_angle = n_min; n_angle < n_max; n_angle++)
    for (n_angle = 0; n_angle < filopodia.size; n_angle++)
    {
        filo_wall = false;

        angle = move.angle + filopodia.directions[n_angle];

        // ================================ //
        // test filopodia/wall interactions //
        // ================================ //

        // effect of wall on previous angles is done at the wall detection step
        // effect future angles is done when their turn comes (in else block)

        // set position/line
        filo_pos  = Point(position.at(0) + cos(angle) * len_filo,
                         position.at(1) + sin(angle) * len_filo);
        filo_line = geosline_from_points(position, filo_pos);

        // weak interaction (filopodia)
        if (env_intersect(filo_line, omp_id))
        {
            interacting = true;

            // wall affinity
            directions_weights[n_angle] *= (1. + wall_afty);
            wall_presence[n_angle] = true;

            // apply partial wall affinity to previous angle if it exists
            //~ if (last_filo_wall != n_angle-1 and n_angle-1 >= n_min)
            if (last_filo_wall != n_angle - 1 and n_angle - 1 >= 0)
            {
                directions_weights[n_angle - 1] *= (1. + 0.5 * wall_afty);
            }

            filo_wall      = true;
            last_filo_wall = n_angle;
        }
        else if (area_intersect(area, filo_line, omp_id))
        {
            interacting = true;

            // check the number of intersections
            n_intersect = area_num_intersections(area, filo_line, filo_pos,
                                                 middle_points, omp_id);

            // get all other areas
            for (i = 0; i < n_intersect; i++)
            {
                tmp_pos       = middle_points[i];
                name_new_area = get_containing_area(tmp_pos, omp_id);
                distance      = sqrt((position.at(0) - tmp_pos[0]) *
                                    (position.at(0) - tmp_pos[0]) +
                                (position.at(1) - tmp_pos[1]) *
                                    (position.at(1) - tmp_pos[1]));

                new_area   = areas_[name_new_area];
                new_height = new_area->get_height();
                new_substrate_affinity =
                    new_area->get_property(names::substrate_affinity);

                delta_h = abs(new_height - old_height);

                if (old_height == new_height)
                {
                    // same height => proba from substrate affinity
                    if (n_intersect == 1)
                    {
                        directions_weights[n_angle] *= new_substrate_affinity;
                    }
                }
                else if (old_height > new_height)
                {
                    // growth cone may go down:
                    // - if filopodia can touch the bottom, the probability of
                    //   making a down move is increased or decreased
                    //   depending on the substrate affinity
                    // - if the bottom is too far, only take the probability
                    //   of making a down move
                    h = filopodia.finger_length - distance - delta_h;

                    delta_proba = h > 0 ? new_substrate_affinity * h /
                                              filopodia.finger_length
                                        : 0.;

                    delta_proba += proba_down_move;
                    //~ delta_proba = proba_down_move/substep;
                    //~ delta_proba = proba_down_move/sqrt(substep);
                    //~ delta_proba *=
                    //proba_down_move*substep*sin(0.5*M_PI/substep); ~
                    //delta_proba *= proba_down_move*sin(0.5*M_PI/substep); ~
                    //delta_proba *=
                    //2*proba_down_move*(0.5+cos(0.5*M_PI/substep))/sqrt(substep);
                    //~ delta_proba += proba_down_move/pow(substep, 1./3.);
                    //~ delta_proba +=
                    //proba_down_move*sqrt(1-cos(0.5*M_PI/sqrt(substep))); ~
                    //delta_proba +=
                    //proba_down_move*pow(1-cos(0.5*M_PI/sqrt(substep)), 1/3.);
                    //~ delta_proba +=
                    //proba_down_move*sqrt(sin(0.5*M_PI/sqrt(substep))); ~
                    //delta_proba +=
                    //proba_down_move/(1+sin(0.5*M_PI/sqrt(substep)));
                    directions_weights[n_angle] *= delta_proba;
                    //~ directions_weights[n_angle] *= new_substrate_affinity;
                }
                else
                {
                    // growth cone may go upwards:
                    // - move is possible if the step is smaller than
                    //   `max_height_up_move`, in [0, new_substrate_affinity]
                    // - if step is higher, move is impossible
                    directions_weights[n_angle] *=
                        new_substrate_affinity *
                        std::exp(-delta_h / max_height_up_move);
                }
            }

            // set wall affinity if necessary
            std::string name_new_area =
                get_containing_area(middle_points[0], omp_id);
            new_area   = areas_[name_new_area];
            new_height = new_area->get_height();

            if (old_height < new_height)
            {
                // apply partial wall affinity to previous angle if it exists
                //~ if (last_filo_wall != n_angle-1 and n_angle-1 >= n_min)
                if (last_filo_wall != n_angle - 1 and n_angle - 1 >= 0)
                {
                    directions_weights[n_angle - 1] *= (1. + 0.5 * wall_afty);
                }

                directions_weights[n_angle] *= (1. + wall_afty);
                wall_presence[n_angle] = true;

                filo_wall      = true;
                last_filo_wall = n_angle;
            }
        }
        else
        {
            directions_weights[n_angle] *= subs_afty;
        }

        // ================================== //
        // test lamelipodia/wall interactions //
        // ================================== //

        // interaction with lamelipodia is twice stronger compared to filo/wall
        if (filo_wall)
        {
            lamel_pos  = Point(position.at(0) + 0.5 * cos(angle) * len_filo,
                              position.at(1) + 0.5 * sin(angle) * len_filo);
            lamel_line = geosline_from_points(position, lamel_pos);

            if (env_intersect(lamel_line, omp_id))
            {
                directions_weights[n_angle] *= 2;
            }
            else if (area_intersect(area, lamel_line, omp_id))
            {
                // everything was already computed for filopodia
                directions_weights[n_angle] *= 2;
            }
        }
        else
        {
            // apply partial wall affinity to next angle if relevant
            //~ if (n_angle != n_min and last_filo_wall == n_angle-1)
            if (n_angle != 0 and last_filo_wall == n_angle - 1)
            {
                directions_weights[n_angle] *= (1. + 0.5 * wall_afty);
            }
        }
    }

    return interacting;
}


void SpaceManager::move_possibility(std::vector<double> &directions_weights,
                                    std::vector<std::string> &new_pos_area,
                                    const Filopodia &filopodia,
                                    const Point &position, const Move &move,
                                    double substep, double sqrt_resol,
                                    unsigned int delta_filo)
{
    // Useful values
    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    unsigned int n_angle;
    Point target_pos, starting_pos;
    double angle;

    // compute the number of filopodia to ignore
    //~ unsigned int ignore = 0.5*((sqrt_resol - sqrt(substep)) / (sqrt_resol -
    //1) * delta_filo); ~ unsigned int n_max  = filopodia.size - ignore; ~
    //unsigned int n_min  = ignore;

    //~ // set ignored direction weights to nan
    //~ for (n_angle = 0; n_angle < n_min; n_angle++)
    //~ {
    //~ directions_weights[n_angle] = std::nan("");
    //~ }
    //~ for (n_angle = n_max; n_angle < filopodia.size; n_angle++)
    //~ {
    //~ directions_weights[n_angle] = std::nan("");
    //~ }

    starting_pos = Point(position.at(0), position.at(1));
    assert(env_contains(starting_pos, omp_id));

    // test the environment for each of the filopodia's angles
    //~ for (n_angle = n_min; n_angle < n_max; n_angle++)
    for (n_angle = 0; n_angle < filopodia.size; n_angle++)
    {
        angle = move.angle + filopodia.directions.at(n_angle);

        target_pos = Point(position.at(0) + cos(angle) * move.module,
                           position.at(1) + sin(angle) * move.module);

        GeomPtr line = geosline_from_points(starting_pos, target_pos);

        if (env_intersect(line, omp_id))
        {
            directions_weights[n_angle] = std::nan(""); // cannot escape env
        }
        else
        {
            new_pos_area[n_angle] = get_containing_area(target_pos, omp_id);
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


double SpaceManager::get_wall_distance(const Point &position, int omp_id) const
{
    GeomPtr geospoint = geospoint_from_point(position);

    GEOSGeometry *boundary;
    double distance;

    boundary = GEOSBoundary_r(context_handler_,
                              environment_manager_->get_environment(omp_id));

    GEOSDistance_r(context_handler_, boundary, geospoint.get(), &distance);

    GEOSGeom_destroy_r(context_handler_, boundary);

    return distance;
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

    // always add last point
    points.push_back(tgt);

    GEOSGeom_destroy_r(context_handler_, intersect);

    return num_intersect;
}


/**
 * @brief returns the absolute value of the angle widening necessary to unstuck
 */
double SpaceManager::unstuck_angle(const Point &position, double current_angle,
                                   double radius, const std::string &area,
                                   int omp_id)
{
    AreaPtr a = areas_.at(area);

    GEOSGeom p;

    GEOSCoordSeq sq = GEOSCoordSeq_create_r(context_handler_, 1, 2);
    GEOSCoordSeq_setX_r(context_handler_, sq, 0, position.at(0));
    GEOSCoordSeq_setY_r(context_handler_, sq, 0, position.at(1));
    p = GEOSGeom_createPoint_r(context_handler_, sq);

    GEOSGeom disk   = GEOSBuffer_r(context_handler_, p, radius, 50);
    GEOSGeom circle = GEOSBoundary_r(context_handler_, disk);
    GEOSGeom a_border =
        GEOSBoundary_r(context_handler_, a->get_shape(omp_id).get());

    GEOSGeom intersect = GEOSIntersection_r(context_handler_, a_border, circle);

    size_t num_intersect = GEOSGetNumGeometries_r(context_handler_, intersect);
    double x0, x1, y0, y1, angle_first(0.), angle_last(3.15);

    if (num_intersect > 0)
    {
        const GEOSGeometry *g0 =
            GEOSGetGeometryN_r(context_handler_, intersect, 0);
        auto coords1 = GEOSGeom_getCoordSeq_r(context_handler_, g0);
        GEOSCoordSeq_getX_r(context_handler_, coords1, 0, &x0);
        GEOSCoordSeq_getY_r(context_handler_, coords1, 0, &y0);

        angle_first = fmod(abs(atan2(y0 - position.at(1), x0 - position.at(0)) -
                               current_angle),
                           M_PI);

        if (num_intersect > 1)
        {
            const GEOSGeometry *g1 = GEOSGetGeometryN_r(
                context_handler_, intersect, num_intersect - 1);
            auto coords2 = GEOSGeom_getCoordSeq_r(context_handler_, g1);
            GEOSCoordSeq_getX_r(context_handler_, coords2, 0, &x1);
            GEOSCoordSeq_getY_r(context_handler_, coords2, 0, &y1);

            angle_last =
                fmod(abs(atan2(y1 - position.at(1), x1 - position.at(0)) -
                         current_angle),
                     M_PI);
            if (angle_last < angle_first)
            {
                angle_first = angle_last;
            }
        }
    }

    GEOSGeom_destroy_r(context_handler_, p);
    GEOSGeom_destroy_r(context_handler_, circle);
    GEOSGeom_destroy_r(context_handler_, disk);
    GEOSGeom_destroy_r(context_handler_, a_border);

    return angle_first;
}


bool SpaceManager::has_environment() const { return environment_initialized_; }


void SpaceManager::set_environment(
    GEOSGeom environment, const std::vector<GEOSGeom> &areas,
    const std::vector<double> &heights, const std::vector<std::string> &names,
    const std::vector<std::unordered_map<std::string, double>> &properties)
{
    const GEOSGeom env = static_cast<GEOSGeom>(environment);
    assert(GEOSisValid_r(context_handler_, env));
    environment_initialized_ = true;

    environment_manager_ =
        std::unique_ptr<Environment>(new Environment(env, context_handler_));

    for (size_t i = 0; i < areas.size(); i++)
    {
        areas_[names[i]] = std::make_shared<Area>(
            areas[i], context_handler_, heights[i], names[i], properties[i]);
    }
}


void SpaceManager::get_environment(
    GEOSGeom &environment, std::vector<GEOSGeom> &areas,
    std::vector<double> &heights, std::vector<std::string> &names,
    std::vector<std::unordered_map<std::string, double>> &properties) const
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();

    if (environment_manager_ != nullptr)
    {
        environment = environment_manager_->get_environment(omp_id);
    }

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


std::vector<std::string> SpaceManager::get_area_names() const
{
    std::vector<std::string> names;

    for (const auto a : areas_)
    {
        names.push_back(a.first);
    }

    return names;
}


void SpaceManager::get_area_properties(
    const std::string &area,
    std::unordered_map<std::string, double> &prop) const
{
    areas_.at(area)->get_properties(prop);
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

} // namespace growth
