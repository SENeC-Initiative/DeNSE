/*
 * Skeleton.hpp
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

#ifndef SKELETON_H
#define SKELETON_H
// libgrowth includes
#include "elements_types.hpp"
#include "spatial_types.hpp"


namespace growth
{
class Node;
/*!
    Implementation of the skeleton of a ``Neuron`` for easy graphical
    representation.
*/
class Skeleton
{
  public:
    //! Center of mass of the soma
    BPoint soma_position;
    double soma_radius;
    //! Container for the dendrites representation
    SkelNeurite dendrites;
    //! Container for the axon representation
    SkelNeurite axon;
    //! List of the points where a neurite branches
    SkelNeurite branching_points;
    //! List of active growth Cones
    SkelNeurite growth_cones;
    //! Default constructor
    Skeleton();
    //! Constructor to make the representation of a neuron
    Skeleton(const Neuron *, unsigned int resolution);
};
} // namespace growth

#endif /* SKELETON_H */
