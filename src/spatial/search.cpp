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

#include <cassert>
#include <cmath>

#include "Node.hpp"
#include "Branch.hpp"

namespace growth
{

void locate_from_distance(BPoint &xy, double &angle, const BranchPtr branch,
                          double distanceToNode)
{
    /*    //implemnet BINARY SEARCH TREE*/
    // stype size     = branch->size();
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


/**
 * Get point closest to a given distance along the branch
 */
stype get_closest_point(TNodePtr branching_node, double branching_dist)
{
    BranchPtr branch = branching_node->get_branch();

    // remove distance to soma
    double initial_distance = branch->initial_distance_to_soma();

    stype left  = 0;
    stype right = branching_node->get_branch_size() - 1;
    stype mid   = 0.5*right;

    PointArray ppl, ppr, ppm;
    double distl, distr, distm;

    while (left != right - 1)
    {
        ppm = branch->at(mid);
        distm = ppm[2] - initial_distance - branching_dist;

        if (distm < 0)
        {
            left = mid;
            mid  = 0.5*(left + right);
        }
        else
        {
            right = mid;
            mid   = 0.5*(left + right);
        }
    }

    ppl = branch->at(left);
    ppr = branch->at(right);

    distl = std::abs(ppl[2] - initial_distance - branching_dist);
    distr = std::abs(ppr[2] - initial_distance - branching_dist);

    if (distl < distr)
    {
        return left;
    }

    return right;
}


void locate_from_idx(BPoint &xy, double &angle, double &distance,
                     const BranchPtr branch, stype id_x)
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
