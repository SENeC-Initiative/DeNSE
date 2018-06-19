#ifndef WALL_H
#define WALL_H
#ifndef GEOS_USE_ONLY_R_API
#define GEOS_USE_ONLY_R_API
#endif

// C++ include
#include <string>
#include <unordered_map>
#include <vector>

// libgrowth include
#include "spatial_types.hpp"

namespace growth
{

class Wall
{
  public:
    Wall(GEOSGeom area, GEOSContextHandle_t handler);
    ~Wall();

    const GEOSPreparedGeometry *get_wall_area(int omp_id) const;
    const GEOSPreparedGeometry *get_border(int omp_id) const;

  private:
    std::vector<GeomPtr> shape_;
    std::vector<const GEOSPreparedGeometry *> prepared_area_;
    std::vector<const GEOSPreparedGeometry *> prepared_border_;
};
} // namespace growth

#endif /* AREAS_H */
