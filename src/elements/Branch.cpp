#include "Branch.hpp"

#include <cmath>
#define _USE_MATH_DEFINES
#include <assert.h>
#include <iostream>

namespace growth
{

Branch::Branch()
{
    std::vector<double> x, y, l;
    points = {{x, y, l}};
}

Branch::Branch(const Branch &copy)
{
    points[0].insert(points[0].end(), copy.points[0].begin(),
                     copy.points[0].end());
    points[1].insert(points[1].end(), copy.points[0].begin(),
                     copy.points[1].end());
    points[2].insert(points[2].end(), copy.points[2].begin(),
                     copy.points[2].end());
}


Branch::Branch(const Point &initial_position, double initial_length)
    : Branch()
{
    points[0].push_back(initial_position.at(0));
    points[1].push_back(initial_position.at(1));
    assert(!points.empty());
    points[2].push_back(initial_length);
}

Branch::~Branch() {}

void Branch::add_point(const Point &p, double length)
{
    points[0].push_back(p.at(0));
    points[1].push_back(p.at(1));
    points[2].push_back(points[2].back() + length);


    // printf("x, y, l : %f, %f, %f (%f) \n", points[0].back(),
    // points[1].back(),
    // points[2].back(), length);
    // printf(" branch size is %lu \n", size());
}

double Branch::module_from_points(const Point &p)
{
    return sqrt(pow(p.at(0) - points[0].back(), 2) +
                pow(p.at(1) - points[1].back(), 2));
}


void Branch::set_first_point(const Point p, double length)
{

    assert(!points.empty());
    points[0][0] = p.at(0);
    points[1][0] = p.at(1);
    points[2][0] = length;
}


void Branch::retract()
{
    points[0].pop_back();
    points[1].pop_back();
    points[2].pop_back();
}
void Branch::resize_tail(size_t new_size)
{
    points[0].resize(new_size);
    points[1].resize(new_size);
    points[2].resize(new_size);
}
/**
 * @brief Resize the head of the Branch, last segment stays, first segment is
 * cut out
 *
 * @param id_x  first element of new Branch
 *
 * @return Branch object of size Branch.size - id_x
 */
BranchPtr Branch::resize_head(size_t id_x) const
{
    assert(size() > id_x);
    BranchPtr new_branch = std::make_shared<Branch>();
    new_branch->points[0].insert(new_branch->points[0].end(),
                                 points[0].cbegin() + id_x, points[0].cend());
    new_branch->points[1].insert(new_branch->points[1].end(),
                                 points[1].cbegin() + id_x, points[1].cend());
    new_branch->points[2].insert(new_branch->points[2].end(),
                                 points[2].cbegin() + id_x, points[2].cend());
    return new_branch;
}


void Branch::append_branch(BranchPtr appended_branch)
{
    size_t total_size = 0;
    total_size += appended_branch->size() + size();

    // auto id_end = appended_branch->size();
    points[0].insert(points[0].end(), appended_branch->points[0].cbegin(),
                     appended_branch->points[0].cend());
    points[1].insert(points[1].end(), appended_branch->points[1].cbegin(),
                     appended_branch->points[1].cend());
    points[2].insert(points[2].end(), appended_branch->points[2].cbegin(),
                     appended_branch->points[2].cend());

    assert(points[2].size() == total_size);
}
PointArray Branch::get_last_point() const
{
    return {{points[0].back(), points[1].back(), points[2].back()}};
}

Point Branch::get_last_xy() const
{
    return Point(points[0].back(), points[1].back());
}

double Branch::get_distance_to_soma() const { return points[2].back(); }

PointArray Branch::at(size_t idx) const
{
    return {{points[0].at(idx), points[1].at(idx), points[2].at(idx)}};
}

size_t Branch::size() const { return points[0].size(); }
}
