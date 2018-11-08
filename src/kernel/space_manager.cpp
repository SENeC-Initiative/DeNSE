#include "space_manager.hpp"

// C++ includes
#include <cmath>
#include <deque>

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

const double epsilon = 1e-4;


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
    double lamel_factor = 2.;

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

    // values used locally inside loop
    double angle, distance, min_wall_dist;
    bool filo_wall, interacting(false), keep_point, intsct;
    Point filo_pos, tmp_pos;
    // intersections of the filopodia
    double x0(position.at(0)), y0(position.at(1)), x, y, point_angle;
    std::unordered_map<std::string, AreaPtr> tested_areas;
    std::vector<double> distances, affinities;
    std::vector<Point> intersections;
    //~ std::vector<AreaPtr> new_areas;
    GeomPtr filo_line, tmp_line;
    std::vector<size_t> indices;
    std::vector<bool> is_wall;
    std::string name_new_area;
    size_t n_intersect;
    AreaPtr new_area;
    GEOSGeom border;
    // recursive loop
    std::deque<AreaPtr> area_deque;
    std::deque<Point> area_points;
    // computation of the affinity
    double affinity, old_dist, current_dist;
    // loop indices
    unsigned int i, n_angle, s_i;

    // test the environment for each of the filopodia's angles
    //~ for (n_angle = n_min; n_angle < n_max; n_angle++)
    for (n_angle = 0; n_angle < filopodia.size; n_angle++)
    {
        filo_wall = false;

        angle = move.angle + filopodia.directions[n_angle];

        // clear all containers
        distances.clear();
        affinities.clear();
        is_wall.clear();
        tested_areas.clear();
        area_deque.clear();
        area_points.clear();
        //~ new_areas.clear();

        // set position/line
        filo_pos  = Point(position.at(0) + cos(angle) * len_filo,
                         position.at(1) + sin(angle) * len_filo);
        filo_line = geosline_from_points(position, filo_pos);

        // set first (current) affinity
        affinities.push_back(subs_afty);

        // ============================== //
        // Find lamelipodia intersections //
        // ============================== //

        // in a first time, we get all the intersections between the filopodia
        // and the environment (walls + areas)

        min_wall_dist = std::numeric_limits<double>::infinity();

        // intersections with areas (and indirectly with walls) "recursively"
        if (area_intersects(area, filo_line, omp_id))
        {
            interacting = true;

            area_points.push_back(position);
            area_deque.push_back(old_area);

            while (not area_deque.empty())
            {
                // get new area and its associated point
                new_area      = area_deque.front();
                tmp_pos       = area_points.front();
                name_new_area = new_area->get_name();
                tmp_line      = geosline_from_points(tmp_pos, filo_pos);

                auto it = tested_areas.find(name_new_area);
                intsct  = it == tested_areas.end();

                // remove from the deque
                area_deque.pop_front();
                area_points.pop_front();

                // test intersections only if not already done
                if (intsct and area_intersects(name_new_area, tmp_line, omp_id))
                {
                    // add area into tested areas
                    tested_areas[name_new_area] = new_area;
                    // check intersections
                    intersections.clear();
                    border = GEOSBoundary_r(
                        context_handler_, new_area->get_shape(omp_id).get());
                    get_intersections(
                        tmp_line, border, intersections);

                    for (auto p : intersections)
                    {
                        // compute distance
                        x = p[0];
                        y = p[1];
                        distance = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));

                        keep_point  = std::abs(x-tmp_pos[0]) > epsilon or
                                      std::abs(y-tmp_pos[1]) > epsilon;
                        keep_point *= (min_wall_dist - distance > epsilon);

                        // intersection must be before the wall, we also exclude
                        // the wall intersection itself because we did it in
                        // a previous iteration
                        if (keep_point)
                        {
                            // we have an intersection before the wall
                            // compute the point just after the intersection
                            point_angle = atan2(y-y0, x-x0);
                            tmp_pos     = Point(
                                x0 + (distance + epsilon)*cos(point_angle),
                                y0 + (distance + epsilon)*sin(point_angle));

                            // save distance
                            distances.push_back(distance);

                            // get containing area
                            name_new_area = get_containing_area(tmp_pos, omp_id);
                            it            = tested_areas.find(name_new_area);

                            if (name_new_area.empty())
                            {
                                filo_wall = true;
                                //~ new_areas.push_back(nullptr);
                                is_wall.push_back(true);
                                affinities.push_back(wall_afty);
                                // compute distance of first wall
                                min_wall_dist = std::min(distance, min_wall_dist);
                            }
                            else
                            {
                                //~ new_areas.push_back(areas_[name_new_area]);
                                is_wall.push_back(false);
                                affinities.push_back(
                                    areas_[name_new_area]->get_property(
                                        names::substrate_affinity));
                            }

                            // store the new areas with which we should check
                            // for intersections in the deque
                            if (not name_new_area.empty()
                                and it == tested_areas.end())
                            {
                                area_deque.push_back(areas_[name_new_area]);
                                area_points.push_back(tmp_pos);
                            }
                        }
                    }
                }
            }
        }

        // ===================== //
        // Sort them by distance //
        // ===================== //

        // sort distances and get indices
        indices.resize(distances.size());
        std::iota(indices.begin(), indices.end(), 0);
        auto comparator = [&distances](int a, int b){ return distances[a] < distances[b]; };
        std::sort(indices.begin(), indices.end(), comparator);

        // we are going to add a last entry, push back indices
        indices.push_back(distances.size());
        // add the last point
        distances.push_back(len_filo);
        is_wall.push_back(false);
        // note that we don't need to push-back affinities because they got
        // the additional entry from the start (the first affinity where the
        // filopodia initiates)

        // ============================ //
        // Compute the average affinity //
        // ============================ //

        old_dist = 0.;
        affinity = indices.empty() ? len_filo*subs_afty : 0.;

        for (i=0; i < indices.size(); i++)
        {
            // sorted index
            s_i          = indices[i];
            current_dist = distances[s_i];

            // are we compressing the lamellipodia (half its width)?
            if (current_dist < 0.25*len_filo)
            {
                //~ printf("%i - lamel begin %f vs %f\n", i, current_dist, old_dist);
                // if its a wall, then affinity goes to NaN and we break
                if (is_wall[s_i])
                {
                    affinity = std::nan("");
                    break;
                }
                else
                {
                    //~ printf("affinity gets %f = %f * %f * %f\n", lamel_factor*affinities[s_i]*(current_dist - old_dist), lamel_factor, affinities[s_i], (current_dist - old_dist));
                    affinity += lamel_factor*affinities[s_i]*(current_dist - old_dist);
                }
            }  // are we dealing with the lamellipodia?
            else if (current_dist < 0.5*len_filo)
            {
                //~ printf("%i - lamel %f vs %f\n", i, current_dist, old_dist);
                affinity += lamel_factor*affinities[s_i]*(current_dist - old_dist);
            }  // are we in-between
            else if (old_dist < 0.5*len_filo)
            {
                //~ printf("%i - mixed filo/lamel %f vs %f\n", i, current_dist, old_dist);
                // part for lamellipodia
                affinity += lamel_factor*affinities[s_i]*(0.5*len_filo - old_dist);
                // part for filopodia
                affinity += affinities[s_i]*(current_dist - 0.5*len_filo);
            }
            else
            {
                //~ printf("%i - filo %f vs %f\n", i, current_dist, old_dist);
                affinity += affinities[s_i]*(current_dist - old_dist);
            }

            // check where we stand: did we reach a wall or do we keep going?
            if (is_wall[s_i])
            {
                // we will break the loop so we jump directly to the last
                // point (the filopodia length)
                affinity += wall_afty*(len_filo - current_dist);
                break;
            }
            else
            {
                old_dist = current_dist;
            }
        }
        affinity = std::max(affinity / len_filo, 0.);
        directions_weights[n_angle] = affinity;

        //~ std::string s = "angle " + std::to_string(n_angle) + " got affinity " + std::to_string(affinity) + " from";

        //~ for (size_t j=0; j < distances.size(); j++)
        //~ {
            //~ s += " d " + std::to_string(distances[j]) + " a " + std::to_string(affinities[j]);
        //~ }

        //~ printf("%s\n", s.c_str());

        // question: what about the transitions for up/down?
        // answer: move to compute_accessibility! (noob)
    }

    //~ std::string weights("");
    //~ for (auto w : directions_weights)
    //~ {
        //~ weights += std::to_string(w) + " ";
    //~ }
    //~ printf("%s\n", weights.c_str());

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

        if (env_intersects(line, omp_id))
        {
            directions_weights[n_angle] = std::nan(""); // cannot escape env
        }
        else
        {
            new_pos_area[n_angle] = get_containing_area(target_pos, omp_id);
        }
    }
}


bool SpaceManager::env_intersects(GeomPtr line, int omp_id) const
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


bool SpaceManager::area_intersects(const std::string &area, GeomPtr line,
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


void SpaceManager::get_intersections(GeomPtr line, GEOSGeom geometry,
                                     std::vector<Point> &points) const
{
    if (environment_initialized_)
    {
        GEOSGeometry *intersect = GEOSIntersection_r(
            context_handler_, geometry, line.get());

        size_t num_intersect = GEOSGetNumGeometries_r(context_handler_,
                                                      intersect);

        double x, y;

        for (unsigned int i = 0; i < num_intersect; i++)
        {
            const GEOSGeometry *g0 =
                GEOSGetGeometryN_r(context_handler_, intersect, i);
            auto coords = GEOSGeom_getCoordSeq_r(context_handler_, g0);
            GEOSCoordSeq_getX_r(context_handler_, coords, 0, &x);
            GEOSCoordSeq_getY_r(context_handler_, coords, 0, &y);
            points.push_back(Point(x, y));
        }

        GEOSGeom_destroy_r(context_handler_, intersect);
    }
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

        angle_first = fmod(std::abs(atan2(y0 - position.at(1), x0 - position.at(0)) -
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
                fmod(std::abs(atan2(y1 - position.at(1), x1 - position.at(0)) -
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

    // if we get there, it means that this point is out of the environment
    return "";
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
    set_param(status, "environment_initialized", environment_initialized_, "");
    // set_param(status, "dimension", Space::get_dimension());
}

} // namespace growth
