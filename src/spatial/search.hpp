/*
 * search.hpp
 *
 * This file is part of DeNSE.
 *
 * Copyright (C) 2019 SeNEC Initiative
 *
 * DeNSE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * DeNSE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DeNSE. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SEARCH_H
#define SEARCH_H
#ifndef GEOS_USE_ONLY_R_API
#define GEOS_USE_ONLY_R_API
#endif

#define _USE_MATH_DEFINES

#include "elements_types.hpp"
#include "spatial_types.hpp"

#include <geos_c.h>


namespace growth
{

void locate_from_distance(BPoint &xy, double &angle, const BranchPtr branch,
                          double distanceToNode);
void locate_from_idx(BPoint &xy, double &angle, double &distance,
                     const BranchPtr branch, size_t id_x);
} // namespace growth

#endif /* SEARCH_H */
