/*
 * space_manager.hpp
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

#ifndef SPACE_M_H
#define SPACE_M_H
#ifndef GEOS_USE_ONLY_R_API
#define GEOS_USE_ONLY_R_API
#endif

#define _USE_MATH_DEFINES

// C++ include
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

// libgrowth include
#include "config.hpp"
#include "elements_types.hpp"
#include "growth_space.hpp"
#include "spatial_types.hpp"

// kernel includes
#include "manager_interface.hpp"


namespace growth
{

// forward declare spatial classes
class Environment;


class SpaceManager : public ManagerInterface
{
  public:
    SpaceManager();

    virtual void initialize();
    virtual void finalize();

    void init_spatial_grid(int local_num_threads);

    void set_status(const statusMap &config);
    void get_status(statusMap &status) const;

    void add_object(BMultiPolygonPtr geom, const ObjectInfo &info, int omp_id);
    void add_object(const BPoint &start, const BPoint &stop, double diam,
                    double length, double taper, const ObjectInfo &info,
                    BranchPtr b, int omp_id);
    void update_objects_branching(TNodePtr old_node, NodePtr new_node,
                                  stype branching_point, stype neuron,
                                  const std::string &neurite, int omp_id);
    void remove_object(const BBox &box, const ObjectInfo &info, int omp_id);
    void get_objects_in_range(const BPoint &p, double radius,
                              std::vector<ObjectInfo> &v) const;
    void get_intersected_objects(const BPoint &start, const BPoint &stop,
                                 std::vector<ObjectInfo> &v) const;
    void get_intersected_objects(const BPoint &start, const BPoint &stop,
                                 std::vector<ObjectInfo> &vi,
                                 std::vector<BPolygonPtr> &vn) const;
    void update_rtree();

    inline BLineString line_from_points(const BPoint &pointA,
                                        const BPoint &pointB) const;
    bool env_contains(const BPoint &point) const;
    bool is_inside(const BPoint &point, stype neuron,
                   const std::string &neurite, double radius,
                   BPolygon &polygon) const;

    inline bool is_close(const BPoint &p1, const BPoint &p2) const;

    double get_wall_distance(const BPoint &position, int omp_id) const;

    bool intersects(const std::string &object_name,
                    const BLineString &line) const;
    bool intersects(const BPolygon &object, const BLineString &line) const;
    void get_intersections(const BLineString &line, const BPolygon &geometry,
                           BMultiPoint &p) const;
    void get_intersections(const BLineString &line,
                           const BMultiLineString &geometry,
                           BMultiPoint &p) const;

    bool get_point_at_distance(const BLineString &line,
                               const std::string &geom_name, double radius,
                               BPoint &p, double &distance);
    bool get_point_at_distance(const BLineString &line,
                               const BPolygonPtr polygon, double radius,
                               BPoint &p, double &distance);

    double unstuck_angle(const BPoint &position, double current_angle,
                         double radius, const std::string &area, int omp_id);

    bool sense(std::vector<double> &directions_weights,
               std::vector<bool> &wall_presence, const Filopodia &filopodia,
               const BPoint &position, const Move &move,
               const std::string &area, double proba_down_move,
               double max_height_up_move, Affinities aff_values, double substep,
               double radius, GCPtr gc_ptr, space_tree_map &neighbors);

    void check_accessibility(std::vector<double> &directions_weights,
                             const Filopodia &filopodia, const BPoint &position,
                             const Move &move, const BPolygonPtr last_segment);

    void add_putative_synapse(
        int omp_id, stype other_neuron, const std::string& other_neurite,
        stype neuron_id, const std::string& neurite_name);

    void check_synaptic_site(const BPoint &position, double distance,
                             stype neuron_id, const std::string &neurite_name,
                             stype other_neuron,
                             const std::string &other_neurite,
                             BPolygonPtr poly);
    void generate_synapses_crossings(
        double synapse_density, bool only_new_syn, bool autapse_allowed,
        const std::set<stype> &presyn_pop, const std::set<stype> &postsyn_pop,
        std::vector<stype> &presyn_neurons, std::vector<stype> &postsyn_neurons,
        std::vector<std::string> &presyn_neurites,
        std::vector<std::string> &postsyn_neurites,
        std::vector<stype> &presyn_nodes, std::vector<stype> &postsyn_nodes,
        std::vector<stype> &presyn_segments,
        std::vector<stype> &postsyn_segments, std::vector<double> &pre_syn_x,
        std::vector<double> &pre_syn_y);
    void generate_synapses_all(
        double spine_density, bool only_new_syn, bool autapse_allowed,
        const std::set<stype> &presyn_pop, const std::set<stype> &postsyn_pop,
        std::vector<stype> &presyn_neurons, std::vector<stype> &postsyn_neurons,
        std::vector<std::string> &presyn_neurites,
        std::vector<std::string> &postsyn_neurites,
        std::vector<stype> &presyn_nodes, std::vector<stype> &postsyn_nodes,
        std::vector<stype> &presyn_segments,
        std::vector<stype> &postsyn_segments, std::vector<double> &pre_syn_x,
        std::vector<double> &pre_syn_y, std::vector<double> &post_syn_x,
        std::vector<double> &post_syn_y);

    void set_environment(
        GEOSGeom environment, const std::vector<GEOSGeom> &areas,
        const std::vector<double> &heights,
        const std::vector<std::string> &names,
        const std::vector<std::unordered_map<std::string, double>> &properties);
    void get_environment(
        GEOSGeom &environment, std::vector<GEOSGeom> &areas,
        std::vector<double> &heights, std::vector<std::string> &names,
        std::vector<std::unordered_map<std::string, double>> &properties) const;
    const BRing &get_env_border(int omp_id) const;
    GEOSContextHandle_t get_context_handler() const;
    void destroy_geom(GEOSGeom geom) const;
    void copy_polygon(BMultiPolygonPtr copy, const BPolygon &p);
    void copy_polygon(BMultiPolygonPtr copy, const BMultiPolygon &p);
    BPolygon make_disk(BPoint position, double radius) const;

    int get_region_thread(const BPoint &position) const;
    int get_region_thread(double x, double y) const;

    std::vector<std::string> get_area_names() const;
    std::string get_containing_area(const BPoint &position) const;
    AreaPtr get_area(const std::string &name) const;
    void
    get_area_properties(const std::string &area,
                        std::unordered_map<std::string, double> &prop) const;

    void num_threads_changed(int num_omp);

    bool has_environment() const;
    bool interactions_on() const;

  private:
    // buffer strategies
    int points_per_circle_;
    bg::strategy::buffer::join_round join_strategy_;
    boost::geometry::strategy::buffer::end_flat end_strategy_;
    bg::strategy::buffer::point_circle circle_strategy_;
    bg::strategy::buffer::side_straight side_strategy_;
    bool initialized_;
    bool environment_initialized_;
    bool interactions_; // whether neurites interact together
    GEOSContextHandle_t context_handler_;
    std::unique_ptr<Environment> environment_manager_;
    std::unordered_map<std::string, AreaPtr> areas_;
    bgi::rtree<RtreeValue, bgi::quadratic<16>> rtree_;
    space_tree_map map_geom_;
    std::vector<space_tree_map> geom_add_buffer_;
    std::vector<std::vector<box_tree_tuple>> box_buffer_;
    // potential synaptic sites
    double max_syn_distance_;
    BMultiPolygon known_synaptic_sites_;
    std::uniform_real_distribution<double> uniform_;
    std::vector<BPoint> old_potential_synapse_crossing_;
    std::vector<BPoint> new_potential_synapse_crossing_;
    std::vector<BPoint> old_potential_synapse_near_;
    std::vector<BPoint> new_potential_synapse_near_;
};


inline BLineString SpaceManager::line_from_points(const BPoint &pointA,
                                                  const BPoint &pointB) const
{
    return BLineString({pointA, pointB});
}


inline bool SpaceManager::is_close(const BPoint &p1, const BPoint &p2) const
{
    return (std::abs(p1.x() - p2.x()) < 1e-6 and
            std::abs(p1.y() - p2.y()) < 1e-6);
}

} // namespace growth

#endif // SPACE_M_H
