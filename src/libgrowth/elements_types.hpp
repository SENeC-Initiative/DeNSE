#ifndef ELT_TYPE_H
#define ELT_TYPE_H

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "config.hpp"
#include "growth_time.hpp"
#include "spatial_types.hpp"

namespace growth
{

/*
 * The main element classes
 */

class Neuron;
class Neurite;
class BaseNode;
class TopologicalNode;
class Node;
class GrowthCone;
class Branch;
class Branching;
class ActinWave;


/*
 * Handy typedefs
 */
typedef std::vector<std::vector<double>> Random_vecs;
typedef std::unordered_map<std::string, double> Param;
typedef std::array<std::vector<double>, 3> PointsArray;
typedef std::array<double, 3> PointArray;
typedef std::tuple<size_t, double, size_t, std::string, signed char> Event;

typedef struct NodeProp
{
    NodeProp() : n_id(0), p_id(0), diameter(0), dist_to_parent(0) {};
    NodeProp(size_t n, size_t p, double diam, double dtp,
             std::vector<double> pos)
        : n_id(n), p_id(p), diameter(diam), dist_to_parent(dtp), position(pos)
    {};

    size_t n_id;
    size_t p_id;
    double diameter;
    double dist_to_parent;
    std::vector<double> position;
} NodeProp;


/*
 * Typedef for shared and weak pointers
 */
typedef std::shared_ptr<BaseNode> BaseNodePtr;
typedef std::shared_ptr<TopologicalNode> TNodePtr;
typedef std::shared_ptr<Node> NodePtr;
typedef std::shared_ptr<GrowthCone> GCPtr;

typedef std::weak_ptr<Node> NodeWeakPtr;
typedef std::weak_ptr<TopologicalNode> TNodeWeakPtr;
typedef std::weak_ptr<BaseNode> BaseWeakNodePtr;
typedef std::weak_ptr<const GrowthCone> GCWeakPtr;

typedef std::shared_ptr<Branch> BranchPtr;

typedef std::shared_ptr<Neurite> NeuritePtr;
typedef std::weak_ptr<const Neurite> NeuriteWeakPtr;

typedef std::shared_ptr<Branching> BranchingPtr;
typedef std::weak_ptr<const Branching> BranchingWeakPtr;

typedef std::shared_ptr<Neuron> NeuronPtr;
typedef std::weak_ptr<const Neuron> NeuronWeakPtr;

typedef std::shared_ptr<ActinWave> ActinPtr;
typedef std::weak_ptr<const ActinWave> ActinWeakPtr;

typedef std::shared_ptr<std::mt19937> mtPtr;


/*
 *typedef for visualization of Neuron neurites
 */
typedef std::unordered_map<std::string, NeuritePtr> NeuriteMap;
//! SkelNeurite is a vectyor of 2D points
typedef std::vector<std::vector<double>> SkelSomas;
typedef std::pair<std::vector<double>, std::vector<double>> SkelNeurite;
typedef std::vector<std::tuple<int, int, double, double, double, int>>
    SwcNeurite;
typedef std::shared_ptr<SkelNeurite> SkelPtr;
} // namespace growth

#endif /* ELT_TYPE_H */
