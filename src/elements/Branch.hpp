/*
 * Branch.hpp
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

// Copyright (c) 2017 Copyright Holder All Rights Reserved.

#include <vector>

#ifndef BRANCH_H
#define BRANCH_H

#include "elements_types.hpp"
#include "spatial_types.hpp"

namespace growth
{

/**
The Branch class is responsible for storing the spatial points representing
the neurite in space.
Each Branch instance is associated to a unique TopologicalNode and thus stores
the continuous set of points that defines its trajectory over time.
 */
class Branch
{
  private:
    BPoint initial_point_;
    PointsArray points_;
    std::vector<BPolygonPtr> segments_;
    std::pair<BPoint, BPoint> last_points_;

  public:
    Branch(const Branch &copy);
    //! Create a branch with initial position and if necesary an initial length
    Branch(const BPoint &origin);
    Branch(const BPoint &origin, double initial_length);
    //! default constructor
    Branch();
    ~Branch();

    void add_point(const BPoint &pos, double length, BPolygonPtr poly,
                   const BPoint &lp1, const BPoint &lp2);

    /**
     * @brief Get the module of length between the point and the last in the
     * array.
     *
     * \param BPoint&
     */
    double module_from_points(const BPoint &);

    /**
     * @brief set first element of Branch container
     *
     * @param xy the points where the branch starts
     * @param distanceToSoma the distance from soma at branch start.
     */
    void set_first_point(const BPoint &xy, double distanceToSoma);

    /**
     * @brief resize the head of the branch: first point invariate
     *
     * @param id_x the length of the new branch container
     */
    void resize_tail(stype new_size);

    /**
     * @brief Create a new branch from tail of this Branch: tail invariate.
     *
     *
     * @param id_x the number of elements to delete from the begin of the
     * container
     *
     * @return
     */
    BranchPtr resize_head(stype id_x) const;

    /**
     * @brief Remove one point from the container.
     */
    void retract();

    /**
     * @brief Append a Branch to existing one
     *
     * @param appended_branch Branch to add after the last element of the
     * existing
     * Branch
     *
     */
    void append_branch(BranchPtr appended_branch);

    // Getter Functions
    double final_distance_to_soma() const;
    double initial_distance_to_soma() const;
    double get_segment_length_at(stype idx) const;
    double get_last_segment_length() const;
    double get_length() const;
    const std::pair<BPoint, BPoint> &get_last_points() const;
    const std::vector<double> &get_xlist() const;
    const std::vector<double> &get_ylist() const;
    PointArray get_last_point() const;
    const BPolygonPtr get_last_segment() const;
    const BPolygonPtr get_segment_at(stype idx) const;
    seg_range segment_range() const;
    // double get_normal_direction() const;
    //! return the last xy
    BPoint get_last_xy() const;
    /**
    Return a BPoint object from the 'idx' element of the Branch
    */
    PointArray at(stype idx) const;
    BPoint xy_at(stype idx) const;
    stype size() const;
};
} // namespace growth

#endif // BRANCH_H
