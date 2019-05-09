/*
 * Branch.cpp
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
    , last_points_(copy.last_points_)
    , segments_(copy.segments_)
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
Branch::Branch(const BPoint &initial_position, double initial_distance_to_soma)
    : Branch()
{
    initial_point_ = initial_position;

    points_[0].push_back(initial_position.x());
    points_[1].push_back(initial_position.y());
    points_[2].push_back(initial_distance_to_soma);

    assert(!points_.empty());
}


Branch::~Branch() {}


/**
 * @brief Add a point to the branch
 *
 * \param BPoint& point in x,y to add to the vector
 * \param double l is required for an easy computation.
 * \param BPolygonPtr poly is the full polygon representing the segment.
 */
void Branch::add_point(const BPoint &p, double length, BPolygonPtr poly,
                       const BPoint &lp1, const BPoint &lp2)
{
    points_[0].push_back(p.x());
    points_[1].push_back(p.y());
    points_[2].push_back(points_[2].back() + length);

    last_points_.first  = lp1;
    last_points_.second = lp2;

    segments_.push_back(poly);
}


const std::pair<BPoint, BPoint>& Branch::get_last_points() const
{
    return last_points_;
}


double Branch::module_from_points(const BPoint &p)
{
    return sqrt(pow(p.x() - points_[0].back(), 2) +
                pow(p.y() - points_[1].back(), 2));
}


void Branch::set_first_point(const BPoint &p, double length)
{
    initial_point_ = p;

    if (points_[0].size() == 0)
    {
        points_[0].push_back(p.x());
        points_[1].push_back(p.y());
        points_[2].push_back(length);
    }
    else
    {
        points_[0][0] = p.x();
        points_[1][0] = p.y();
        points_[2][0] = length;
    }
}


void Branch::retract()
{
    points_[0].pop_back();
    points_[1].pop_back();
    points_[2].pop_back();

    // recover the new last points
    if (not segments_.empty())
    {
        std::vector<BPoint> new_lp;

        BPolygonPtr last_poly = segments_.back();
        BRing ring            = last_poly->outer();

        segments_.pop_back();

        if (not segments_.empty())
        {
            BPolygonPtr new_last_poly = segments_.back();

            for (size_t i=0; i < ring.size() - 1; i++)
            {
                if (bg::covered_by(ring[i], *(new_last_poly.get())))
                {
                    new_lp.push_back(ring[i]);
                }
            }

            assert(new_lp.size() == 2);

            last_points_.first  = new_lp[0];
            last_points_.second = new_lp[1];
        }
    }
}


void Branch::resize_tail(size_t new_size)
{
    assert(new_size <= size());

    points_[0].resize(new_size);
    points_[1].resize(new_size);
    points_[2].resize(new_size);

    // for segment, set new size
    size_t size_seg = new_size == 0 ? 0 : new_size - 1;

    // recover the new last points
    if (size_seg > 0 and size_seg < segments_.size())
    {
        std::vector<BPoint> new_lp;

        BPolygonPtr end_poly = segments_[size_seg];
        BPolygonPtr last_poly = segments_[size_seg - 1];

        auto ring = end_poly->outer();

        for (size_t i=0; i < ring.size() - 1; i++)
        {
            if (bg::covered_by(ring[i], *(last_poly.get())))
            {
                new_lp.push_back(ring[i]);
            }
        }

        assert(new_lp.size() == 2);

        last_points_.first  = new_lp[0];
        last_points_.second = new_lp[1];
    }

    segments_.resize(size_seg);
}


/**
 * @brief Resize the head of the Branch, last segment stay, first segments are
 * cut out.
 *
 * @param id_x  first element of new Branch
 *
 * @return Branch object of size Branch.size - id_x
 */
BranchPtr Branch::resize_head(size_t id_x) const
{
    assert(size() > id_x);
    BranchPtr new_branch = std::make_shared<Branch>();

    new_branch->points_[0].insert(
        new_branch->points_[0].end(), points_[0].cbegin() + id_x,
        points_[0].cend());

    new_branch->points_[1].insert(
        new_branch->points_[1].end(), points_[1].cbegin() + id_x,
        points_[1].cend());

    new_branch->points_[2].insert(
        new_branch->points_[2].end(), points_[2].cbegin() + id_x,
        points_[2].cend());

    // resize segments
    new_branch->segments_.insert(
        new_branch->segments_.end(), segments_.cbegin() + id_x,
        segments_.cend());

    // set points
    new_branch->initial_point_ = BPoint(new_branch->points_[0][0],
                                        new_branch->points_[0][1]);

    new_branch->last_points_ = last_points_;

    return new_branch;
}


/**
 * Append branch (inplace operation)
 */
void Branch::append_branch(BranchPtr appended_branch)
{
    size_t total_size = points_[0].size() + appended_branch->size();

    const BPoint lp = get_last_xy();
    const BPoint ip = appended_branch->xy_at(0);
    int dsize       = 0;

    if (std::abs(lp.x() - ip.x()) < 1e-6 and std::abs(lp.y() - ip.y()) < 1e-6)
    {
        dsize = 1;
    }

    total_size -= dsize;

    // auto id_end = appended_branch->size();
    points_[0].insert(points_[0].end(),
                      appended_branch->points_[0].cbegin() + dsize,
                      appended_branch->points_[0].cend());

    points_[1].insert(points_[1].end(),
                      appended_branch->points_[1].cbegin() + dsize,
                      appended_branch->points_[1].cend());

    points_[2].insert(points_[2].end(),
                      appended_branch->points_[2].cbegin() + dsize,
                      appended_branch->points_[2].cend());

    segments_.insert(segments_.end(),
                     appended_branch->segments_.cbegin(),
                     appended_branch->segments_.cend());

    assert(points_[2].size() == total_size);
}


// double Branch::get_normal_direction() const
// {
//     // normal direction to the 
//     double x1, x2, y1, y2;

//     x1 = new_branch->last_points_.first.x();
//     y1 = new_branch->last_points_.first.y();

//     x2 = new_branch->last_points_.second.x();
//     y2 = new_branch->last_points_.second.y();

//     double theta = atan2(-(x2-x1), y2-y1);
//     double current_dir = atan2();
// }


PointArray Branch::get_last_point() const
{
    return {{points_[0].back(), points_[1].back(), points_[2].back()}};
}


BPoint Branch::get_last_xy() const
{
    if (points_[0].size() == 0)
    {
        return initial_point_;
    }
    else
    {
        return BPoint(points_[0].back(), points_[1].back());
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
        double dx = points_[0].back() - initial_point_.x();
        double dy = points_[1].back() - initial_point_.y();
        return sqrt(dx*dx + dy*dy);
    }

    return 0.;
}


const BPolygonPtr Branch::get_last_segment() const
{
    if (segments_.empty())
    {
        return nullptr;
    }

    return segments_.back();
}


const BPolygonPtr Branch::get_segment_at(size_t idx) const
{
    return segments_.at(idx);
}


seg_range Branch::segment_range() const
{
    return boost::make_iterator_range(segments_);
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


BPoint Branch::xy_at(size_t idx) const
{
    return BPoint(points_[0].at(idx), points_[1].at(idx));
}


size_t Branch::size() const
{
    return points_[0].size();
}

} // namespace growth
