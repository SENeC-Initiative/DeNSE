#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#ifndef GEOS_USE_ONLY_R_API
#define GEOS_USE_ONLY_R_API
#endif

// C++ include
#include <vector>

// libgrowth include
#include "config.hpp"
#include "growth_space.hpp"
#include "spatial_types.hpp"

// spatial include
#include "Area.hpp"


namespace growth
{

// forward declare friend classs
class SpaceManager;


class Environment
{
    friend class SpaceManager;

  private:
    Environment(BMultiPolygonPtr environment);
    BMultiPolygonPtr environment_;
    BMultiLineString boundary_;

  public:
    const BMultiPolygonPtr get_environment() const;
    const BMultiLineString& get_boundary() const;
};

} // namespace growth

#endif /* ENVIRONMENT_H */
