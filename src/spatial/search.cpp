#include "search.hpp"

#define _USE_MATH_DEFINES

#include <cassert>
#include <cmath>

#include "Branch.hpp"

namespace growth
{

void locate_from_distance(Point &xy, double &angle, const BranchPtr branch,
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
    // if (xy_1.at(1) - xy_0.at(1))
    //{
    // direction = atan((xy_1.at(0) - xy_0.at(0)) / (xy_1.at(1) - xy_0.at(1)));
    //}
    // else
    //{
    // direction = M_PI / 2.;
    //}
    // xy = Point(targetNode_->get_branch()->at(id_x)[0],
    // targetNode_->get_branch()->at(id_x)[1]);
    /*angle = get_angle(rnd_engine, direction);*/
}

void locate_from_idx(Point &xy, double &angle, double &distance,
                     const BranchPtr branch, size_t id_x)
{
    xy       = branch->xy_at(id_x);
    distance = branch->at(id_x)[2];

    Point xy_1;      // second point to get the local direction of the branch
    double sign = 0; // correct the direction if xy_1 before xy

    // all branches must always be of size at least 2 (branching events occur
    // after the steps).
    assert(branch->size() > 1);
    double delta_x, delta_y;
    if (id_x == branch->size())
    {
        delta_y = xy.at(1) - branch->at(id_x - 1)[1];
        delta_x = xy.at(0) - branch->at(id_x - 1)[0];
    }
    else
    {
        delta_y = xy.at(1) - branch->at(id_x + 1)[1];
        delta_x = xy.at(0) - branch->at(id_x + 1)[0];
        sign    = M_PI;
    }
    // get the principal argument (between -Pi and Pi)
    angle = -sign + atan2(delta_y, delta_x);
    /*    printf("The arc tangent for (x=%f, y=%f) is %f degrees\n", delta_x,
     * delta_y,*/
    /*angle * 180 / M_PI);*/
}
} // namespace growth
