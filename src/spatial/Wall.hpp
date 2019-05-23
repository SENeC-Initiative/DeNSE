/*
 * Wall.hpp
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
