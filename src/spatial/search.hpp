#ifndef SEARCH_H
#define SEARCH_H
#ifndef GEOS_USE_ONLY_R_API
#define GEOS_USE_ONLY_R_API
#endif

#include "elements_types.hpp"
#include "spatial_types.hpp"

#include <geos_c.h>


namespace growth
{

void locate_from_distance(Point &xy, double &angle, const BranchPtr branch,
                          double distanceToNode);
void locate_from_idx(Point &xy, double &angle, double &distance,
                     const BranchPtr branch, size_t id_x);
} // namespace growth

#endif /* SEARCH_H */
