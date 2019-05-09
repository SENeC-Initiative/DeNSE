/*
 * Environment.hpp
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
