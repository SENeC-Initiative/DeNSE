#ifndef ELT_TYPE_H
#define ELT_TYPE_H

#include <array>
#include <memory>
#include <random>
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


// Event type, contains (event_time, neuron, neurite, gc, event_type)
typedef std::tuple<Time, size_t, std::string, int, signed char> Event;

namespace edata
{
  enum edata
  {
    TIME,
    NEURON,
    NEURITE,
    GC,
    EV_TYPE
  };
}


typedef struct NodeProp
{
    NodeProp() : n_id(0), p_id(0), diameter(0), dist_to_parent(0) {};
    NodeProp(size_t n, size_t p, double diam, double dtp,
             const std::vector<double>& pos)
        : n_id(n), p_id(p), diameter(diam), dist_to_parent(dtp), position(pos)
    {};

    size_t n_id;
    size_t p_id;
    double diameter;
    double dist_to_parent;
    std::vector<double> position;
} NodeProp;


typedef struct Affinities
{
    Affinities()
      : affinity_self(std::nan(""))
      , affinity_axon_same_neuron(std::nan(""))
      , affinity_axon_other_neuron(5.)
      , affinity_dendrite_same_neuron(std::nan(""))
      , affinity_dendrite_other_neuron(2.)
      , affinity_soma_same_neuron(std::nan(""))
      , affinity_soma_other_neuron(0.)
    {};

    Affinities(double self, double asn, double aon, double dsn, double don,
               double ssn, double son)
      : affinity_self(self)
      , affinity_axon_same_neuron(asn)
      , affinity_axon_other_neuron(aon)
      , affinity_dendrite_same_neuron(dsn)
      , affinity_dendrite_other_neuron(don)
      , affinity_soma_same_neuron(ssn)
      , affinity_soma_other_neuron(son)
    {};

    Affinities(const Affinities &copy)
      : affinity_self(copy.affinity_self)
      , affinity_axon_same_neuron(copy.affinity_axon_same_neuron)
      , affinity_axon_other_neuron(copy.affinity_axon_other_neuron)
      , affinity_dendrite_same_neuron(copy.affinity_dendrite_same_neuron)
      , affinity_dendrite_other_neuron(copy.affinity_dendrite_other_neuron)
      , affinity_soma_same_neuron(copy.affinity_soma_same_neuron)
      , affinity_soma_other_neuron(copy.affinity_soma_other_neuron)
    {};

    double affinity_self;
    double affinity_axon_same_neuron;
    double affinity_axon_other_neuron;
    double affinity_dendrite_same_neuron;
    double affinity_dendrite_other_neuron;
    double affinity_soma_same_neuron;
    double affinity_soma_other_neuron;
} Affinities;


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
