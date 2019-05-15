#ifndef SEARCH_H
#define SEARCH_H

#include "elements_types.hpp"
#include "spatial_types.hpp"

#include <geos_c.h>


namespace growth
{

void locate_from_distance(Point &xy, double &angle, const BranchPtr branch,
                          double distanceToNode);
void locate_from_idx(Point &xy, double &angle, const BranchPtr branch,
                     size_t id_x);
}

#endif /* SEARCH_H */
