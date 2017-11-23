#ifndef SPACE_M_H
#define SPACE_M_H
#ifndef GEOS_USE_ONLY_R_API
#define GEOS_USE_ONLY_R_API
#endif

// C++ include
#include <memory>
#include <vector>

#include "config.hpp"
#include "elements_types.hpp"
#include "growth_space.hpp"
#include "manager_interface.hpp"
#include "spatial_types.hpp"

namespace growth
{

class Environment
{
    friend class SpaceManager;

  private:
    Environment(GEOSGeom environment, GEOSContextHandle_t handler);
    friend class SpaceManager;
    GEOSGeom environment_;
    std::vector<const GEOSPreparedGeometry*> prepared_env_;
    std::vector<const GEOSPreparedGeometry*> prepared_border_;

  public:
    GEOSGeom get_environment() const;
    const GEOSPreparedGeometry *get_prepared(int omp_id) const;
    const GEOSPreparedGeometry *get_border(int omp_id) const;
};


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
    bool env_contains(const Point &point) const;
    bool env_intersect(GeomPtr line) const;

    int get_region_thread(const Point &position);
    int get_region_thread(double x, double y);

    bool sense(std::vector<double> &directions_weigths,
               const Filopodia &filopodia, const Point &position,
               const Move &move, const double distance, const double in_value,
               const double out_value);

    void set_environment(GEOSGeom environment);
    void get_environment(GEOSGeom &environment) const;
    bool has_environment() const;

  private:
    // std::vector<Space::Shape> spatial_grid;
    bool environment_initialized_;
    GEOSContextHandle_t context_handler_;
    std::unique_ptr<Environment> environment_manager_;
};
}

#endif // SPACE_M_H
