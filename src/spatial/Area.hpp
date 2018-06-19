#ifndef AREA_H
#define AREA_H
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

class Area
{
  public:
    Area(GEOSGeom area, GEOSContextHandle_t handler, double height,
         const std::string &name,
         std::unordered_map<std::string, double> properties);
    ~Area();

    const GEOSPreparedGeometry *get_area(int omp_id) const;
    const GEOSPreparedGeometry *get_border(int omp_id) const;
    double get_property(const std::string &name) const;
    double get_height() const;
    std::string get_name() const;
    void
    get_properties(std::unordered_map<std::string, double> &properties) const;
    GeomPtr get_shape(int omp_id) const;

  private:
    std::vector<GeomPtr> shape_;
    std::string name_;
    std::vector<const GEOSPreparedGeometry *> prepared_area_;
    std::vector<const GEOSPreparedGeometry *> prepared_border_;
    double height_;
    std::unordered_map<std::string, double> properties_;
};
} // namespace growth

#endif /* AREAS_H */
