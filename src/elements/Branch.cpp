#include "Branch.hpp"

#define _USE_MATH_DEFINES

#include <cmath>
#include <cassert>
#include <iostream>


namespace growth
{

Branch::Branch()
    : initial_point_(0, 0)
{
    std::vector<double> x, y, l;
    points_ = {{x, y, l}};
}


Branch::Branch(const Branch &copy)
    : initial_point_(copy.initial_point_)
{
    points_[0].insert(points_[0].end(), copy.points_[0].begin(),
                     copy.points_[0].end());
    points_[1].insert(points_[1].end(), copy.points_[0].begin(),
                     copy.points_[1].end());
    points_[2].insert(points_[2].end(), copy.points_[2].begin(),
                     copy.points_[2].end());
}


/**
 * Create a new Branch with an initial point and a given distance to soma.
 */
Branch::Branch(const Point &initial_position, double initial_distance_to_soma)
    : Branch()
{
    initial_point_ = initial_position;
    points_[0].push_back(initial_position.at(0));
    points_[1].push_back(initial_position.at(1));
    points_[2].push_back(initial_distance_to_soma);
    assert(!points_.empty());
}


Branch::~Branch() {}


void Branch::add_point(const Point &p, double length)
{
    points_[0].push_back(p.at(0));
    points_[1].push_back(p.at(1));
    points_[2].push_back(points_[2].back() + length);
}


double Branch::module_from_points(const Point &p)
{
    return sqrt(pow(p.at(0) - points_[0].back(), 2) +
                pow(p.at(1) - points_[1].back(), 2));
}


void Branch::set_first_point(const Point &p, double length)
{
    initial_point_ = p;

    if (points_[0].size() == 0)
    {
        points_[0].push_back(p.at(0));
        points_[1].push_back(p.at(1));
        points_[2].push_back(length);
    }
    else
    {
        points_[0][0] = p.at(0);
        points_[1][0] = p.at(1);
        points_[2][0] = length;
    }
}


void Branch::retract()
{
    points_[0].pop_back();
    points_[1].pop_back();
    points_[2].pop_back();
}


void Branch::resize_tail(size_t new_size)
{
    assert(new_size <= size());

    points_[0].resize(new_size);
    points_[1].resize(new_size);
    points_[2].resize(new_size);
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
    new_branch->points_[0].insert(new_branch->points_[0].end(),
                                 points_[0].cbegin() + id_x, points_[0].cend());
    new_branch->points_[1].insert(new_branch->points_[1].end(),
                                 points_[1].cbegin() + id_x, points_[1].cend());
    new_branch->points_[2].insert(new_branch->points_[2].end(),
                                 points_[2].cbegin() + id_x, points_[2].cend());

    // done in neurite branching
    //~ Point new_init_point = Point(points_[0].at(id_x), points_[1].at(id_x));
    //~ new_branch->set_first_point(new_init_point, points_[2].at(id_x));

    return new_branch;
}


/**
 * Append branch (inplace operation)
 */
void Branch::append_branch(BranchPtr appended_branch)
{
    size_t total_size = points_[0].size() + appended_branch->size();

    // auto id_end = appended_branch->size();
    points_[0].insert(points_[0].end(), appended_branch->points_[0].cbegin(),
                     appended_branch->points_[0].cend());
    points_[1].insert(points_[1].end(), appended_branch->points_[1].cbegin(),
                     appended_branch->points_[1].cend());
    points_[2].insert(points_[2].end(), appended_branch->points_[2].cbegin(),
                     appended_branch->points_[2].cend());

    assert(points_[2].size() == total_size);
}


PointArray Branch::get_last_point() const
{
    return {{points_[0].back(), points_[1].back(), points_[2].back()}};
}


Point Branch::get_last_xy() const
{
    if (points_[0].size() == 0)
    {
        return initial_point_;
    }
    else
    {
        return Point(points_[0].back(), points_[1].back());
    }
}


double Branch::initial_distance_to_soma() const { return points_[2].front(); }


double Branch::final_distance_to_soma() const { return points_[2].back(); }


double Branch::get_length() const
{
    if (not points_[2].empty())
    {
        return points_[2].back() - points_[2].front();
    }

    return 0.;
}


double Branch::get_segment_length_at(size_t idx) const
{
    if (idx == 0)
    {
        return 0;
    }

    return points_.at(2).at(idx) - points_.at(2).at(idx-1);
}


double Branch::get_last_segment_length() const
{
    if (points_[2].size() >= 2)
    {
        return points_.at(2).back() - points_.at(2).at(size()-2);
    }
    else if (points_[2].size() == 1)
    {
        double dx = points_[0].back() - initial_point_.at(0);
        double dy = points_[1].back() - initial_point_.at(1);
        return sqrt(dx*dx + dy*dy);
    }

    return 0.;
}


const std::vector<double>& Branch::get_xlist() const
{
    return points_.at(0);
}


const std::vector<double>& Branch::get_ylist() const
{
    return points_.at(1);
}


PointArray Branch::at(size_t idx) const
{
    return {{points_[0].at(idx), points_[1].at(idx), points_[2].at(idx)}};
}


Point Branch::xy_at(size_t idx) const
{
    return Point(points_[0].at(idx), points_[1].at(idx));
}


size_t Branch::size() const { return points_[0].size(); }

} // namespace growth
