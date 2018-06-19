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
    Environment(GEOSGeom environment, GEOSContextHandle_t handler);
    std::vector<GEOSGeom> environment_;
    std::vector<const GEOSPreparedGeometry *> prepared_env_;
    std::vector<const GEOSPreparedGeometry *> prepared_border_;

  public:
    ~Environment();
    GEOSGeom get_environment(int omp_id) const;
    const GEOSPreparedGeometry *get_prepared(int omp_id) const;
    const GEOSPreparedGeometry *get_border(int omp_id) const;
};

} // namespace growth

#endif /* ENVIRONMENT_H */
