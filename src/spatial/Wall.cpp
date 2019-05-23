/*
 * Wall.cpp
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

#include "Wall.hpp"

// kernel include
#include "kernel_manager.hpp"


namespace growth
{

Wall::Wall(GEOSGeom wall, GEOSContextHandle_t handler)
{
    assert(0 != wall);
    assert(GEOSisValid_r(handler, wall));

    for (int i = 0; i < kernel().parallelism_manager.get_num_local_threads();
         i++)
    {
        shape_.push_back(GeomPtr(GEOSGeom_clone_r(handler, wall)));
        prepared_area_.push_back(GEOSPrepare_r(handler, wall));
        const GEOSGeom border = GEOSBoundary_r(handler, wall);
        prepared_border_.push_back(GEOSPrepare_r(handler, border));
        assert(prepared_area_[i] != 0);
        assert(prepared_border_[i] != 0);
    }
}


Wall::~Wall()
{
    for (const GEOSPreparedGeometry *shape : prepared_area_)
    {
        delete shape;
    }
    prepared_area_.clear();

    for (const GEOSPreparedGeometry *border : prepared_border_)
    {
        delete border;
    }
    prepared_border_.clear();
}

const GEOSPreparedGeometry *Wall::get_wall_area(int omp_id) const
{
    return prepared_area_[omp_id];
}


const GEOSPreparedGeometry *Wall::get_border(int omp_id) const
{
    return prepared_border_[omp_id];
}

} // namespace growth
