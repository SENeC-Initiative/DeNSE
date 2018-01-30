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
    Point initial_point;

  public:
    PointsArray points;
    Branch(const Branch &copy);
    //! Create a branch with initial position and if necesary an initial length
    Branch(const Point &);
    Branch(const Point &, double);
    //! default constructor
    Branch();
    ~Branch();

    /**
     * @brief Add a point to the branch
     *
     * \param Point& point in x,y to add to the vector
     * \param double l is optional but required for an easy computation.
     */
    void add_point(const Point &, double);

    /**
     * @brief Get the module of length between the point and the last in the
     * array.
     *
     * \param Point&
     */
    double module_from_points(const Point &);

    /**
     * @brief set first element of Branch container
     *
     * @param xy the points where the branch starts
     * @param distanceToSoma the distance from soma at branch start.
     */
    void set_first_point(const Point &xy, double distanceToSoma);


    /**
     * @brief resize the head of the branch: first point invariate
     *
     * @param id_x the length of the new branch container
     */
    void resize_tail(size_t new_size);

    /**
     * @brief Create a new branch from tail of this Branch: tail invariate.
     *
     *
     * @param id_x the number of elements to delete from the begin of the
     * container
     *
     * @return
     */
    BranchPtr resize_head(size_t id_x) const;

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
    double get_distance_to_soma() const;
    //! return a Point object from last point of the Branch
    PointArray get_last_point() const;
    //! return the last xy
    Point get_last_xy() const;
    /**
    Return a Point object from the 'idx' element of the Branch
    */
    PointArray at(size_t idx) const;
    size_t size() const;
};
}

#endif // BRANCH_H
