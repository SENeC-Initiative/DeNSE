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
    Point soma_position;
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
    Skeleton(const Neuron *);
};
}

#endif /* SKELETON_H */
