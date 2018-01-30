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
    bool wall_area_contains(const std::string &area, const Point &point,
                            int omp_id) const;
    bool env_intersect(GeomPtr line, int omp_id) const;
    bool area_intersect(const std::string &area, GeomPtr line,
                        int omp_id) const;
    size_t area_num_intersections(const std::string &area, GeomPtr line,
                                  const Point &tgt, std::vector<Point> &p,
                                  int omp_id) const;

    int get_region_thread(const Point &position) const;
    int get_region_thread(double x, double y) const;

    std::string get_containing_area(const Point &position, int omp_id) const;
    AreaPtr get_area(const std::string &name) const;

    bool sense(std::vector<double> &directions_weigths,
               std::vector<std::string> &new_pos_area,
               const Filopodia &filopodia, const Point &position,
               const Move &move, double distance, double substep,
               double in_value, double out_value, const std::string &area,
               double down_move_proba, double max_height_up_move);

    bool sense_walls(std::vector<double> &directions_weigths,
                     const Filopodia &filopodia, const Point &position,
                     const Move &move, double distance, double substep,
                     double wall_affinity, const std::string &area);

    void set_environment(
        GEOSGeom environment, const std::vector<GEOSGeom> &walls,
        const std::vector<GEOSGeom> &areas, const std::vector<double> &heights,
        const std::vector<std::string> &names,
        const std::vector<std::unordered_map<std::string, double>> &properties);
    void get_environment(
        GEOSGeom &environment, std::vector<GEOSGeom> &areas,
        std::vector<double> &heights, std::vector<std::string> &names,
        std::vector<std::unordered_map<std::string, double>> &properties) const;
    bool has_environment() const;

  private:
    // std::vector<Space::Shape> spatial_grid;
    bool environment_initialized_;
    GEOSContextHandle_t context_handler_;
    std::unique_ptr<Environment> environment_manager_;
    std::unordered_map<std::string, AreaPtr> areas_;
    std::unordered_map<std::string, WallPtr> walls_;
};

} /* namespace */

#endif // SPACE_M_H
