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
    Area(BMultiPolygonPtr area, double height, const std::string &name,
         std::unordered_map<std::string, double> properties);

    const BMultiPolygonPtr get_area() const;
    double get_property(const std::string &name) const;
    double get_height() const;
    std::string get_name() const;
    void
    get_properties(std::unordered_map<std::string, double> &properties) const;
    const BMultiLineString& get_boundary() const;

  private:
    BMultiPolygonPtr shape_;
    BMultiLineString boundary_;
    std::string name_;
    double height_;
    std::unordered_map<std::string, double> properties_;
};
} // namespace growth

#endif /* AREAS_H */
