/*
 * search.cpp
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

#include "search.hpp"

#define _USE_MATH_DEFINES

#include <cassert>
#include <cmath>

#include "Branch.hpp"

namespace growth
{

void locate_from_distance(BPoint &xy, double &angle, const BranchPtr branch,
                          double distanceToNode)
{
    /*    //implemnet BINARY SEARCH TREE*/
    // size_t size     = branch->size();
    // while (distance <= distanceToTarget_)
    //{
    // id_x--;
    // distance = targetNode_->get_branch()->at(id_x)[2];
    //}
    // Point xy_0 = Point(targetNode_->get_branch()->at(id_x)[0],
    // targetNode_->get_branch()->at(id_x)[1]);
    // Point xy_1 = Point(targetNode_->get_branch()->at(id_x + 1)[0],
    // targetNode_->get_branch()->at(id_x + 1)[1]);
    // double direction;
    // if (xy_1.y() - xy_0.y())
    //{
    // direction = atan((xy_1.x() - xy_0.x()) / (xy_1.y() - xy_0.y()));
    //}
    // else
    //{
    // direction = M_PI / 2.;
    //}
    // xy = Point(targetNode_->get_branch()->at(id_x)[0],
    // targetNode_->get_branch()->at(id_x)[1]);
    /*angle = get_angle(rnd_engine, direction);*/
}

void locate_from_idx(BPoint &xy, double &angle, double &distance,
                     const BranchPtr branch, size_t id_x)
{
    xy       = branch->xy_at(id_x);
    distance = branch->at(id_x)[2];

    BPoint xy_1;      // second point to get the local direction of the branch
    double sign = 0; // correct the direction if xy_1 before xy

    // all branches must always be of size at least 2 (branching events occur
    // after the steps).
    assert(branch->size() > 1);
    double delta_x, delta_y;
    if (id_x == branch->size())
    {
        delta_y = xy.y() - branch->at(id_x - 1)[1];
        delta_x = xy.x() - branch->at(id_x - 1)[0];
    }
    else
    {
        delta_y = xy.y() - branch->at(id_x + 1)[1];
        delta_x = xy.x() - branch->at(id_x + 1)[0];
        sign    = M_PI;
    }
    // get the principal argument (between -Pi and Pi)
    angle = -sign + atan2(delta_y, delta_x);
}

} // namespace growth
