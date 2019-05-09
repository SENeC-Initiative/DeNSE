/*
 * Area.hpp
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
