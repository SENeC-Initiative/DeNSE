#ifndef SPACE_M_H
#define SPACE_M_H
#ifndef GEOS_USE_ONLY_R_API
#define GEOS_USE_ONLY_R_API
#endif

// C++ include
#include <memory>
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

    GeomPtr geospoint_from_point(const Point &point) const;
    GeomPtr geosline_from_points(const Point &pointA,
                                 const Point &pointB) const;
    bool env_contains(const Point &point, int omp_id) const;

    double get_wall_distance(const Point &position, int omp_id) const;

    bool env_intersects(GeomPtr line, int omp_id) const;
    bool area_intersects(const std::string &area, GeomPtr line,
                        int omp_id) const;
    void get_intersections(GeomPtr line, GEOSGeom geometry,
                           std::vector<Point> &p) const;

    double unstuck_angle(const Point &position, double current_angle,
                         double radius, const std::string &area, int omp_id);

    bool sense(std::vector<double> &directions_weigths,
               std::vector<bool> &wall_presence, const Filopodia &filopodia,
               const Point &position, const Move &move, const std::string &area,
               double down_move_proba, double max_height_up_move,
               double substep, double radius);

    void move_possibility(std::vector<double> &directions_weights,
                          std::vector<std::string> &new_pos_area,
                          const Filopodia &filopodia, const Point &position,
                          const Move &move, double substep, double sqrt_resol,
                          unsigned int delta_filo);

    void set_environment(
        GEOSGeom environment, const std::vector<GEOSGeom> &areas,
        const std::vector<double> &heights,
        const std::vector<std::string> &names,
        const std::vector<std::unordered_map<std::string, double>> &properties);
    void get_environment(
        GEOSGeom &environment, std::vector<GEOSGeom> &areas,
        std::vector<double> &heights, std::vector<std::string> &names,
        std::vector<std::unordered_map<std::string, double>> &properties) const;
    GEOSGeom get_env_border(int omp_id) const;
    void destroy_geom(GEOSGeom geom) const;

    int get_region_thread(const Point &position) const;
    int get_region_thread(double x, double y) const;

    std::vector<std::string> get_area_names() const;
    std::string get_containing_area(const Point &position, int omp_id) const;
    AreaPtr get_area(const std::string &name) const;
    void
    get_area_properties(const std::string &area,
                        std::unordered_map<std::string, double> &prop) const;

    bool has_environment() const;

  private:
    // std::vector<Space::Shape> spatial_grid;
    bool environment_initialized_;
    GEOSContextHandle_t context_handler_;
    std::unique_ptr<Environment> environment_manager_;
    std::unordered_map<std::string, AreaPtr> areas_;
};

} // namespace growth

#endif // SPACE_M_H
