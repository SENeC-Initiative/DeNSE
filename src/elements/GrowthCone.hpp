#ifndef GROWTHCONE_H
#define GROWTHCONE_H

#include <random>

#include "Neuron.hpp"
#include "Node.hpp"

#include "config_impl.hpp"
#include "elements_types.hpp"
#include "spatial_types.hpp"


namespace growth
{

/*
 * Growth Cone is an abstract Class
 * It is a base model which s overloaded by more detailed models.
 *
 * Growthcone is a (biological) relevant class, it computes the forces acting on
 * the growth cone from environment or intrinsic phenomena; them it actualizes
 * the step, storing points into its branch.
 *
 * A GrowthCone is always at the head of a certain branch.
 * The object is never deleted: when neurites branch, a new GC is created,
 * the other GC skips the branching point and keeps on growing.
 * GrowthCone inherits from Node.
 */
class GrowthCone : public TopologicalNode,
                   public std::enable_shared_from_this<GrowthCone>
{

    friend class Neurite;
    friend class Branching;

  protected:
    // store the state of the kernel to avoid repeted calls to kernel_manager
    bool using_environment_; // whether we're embedded in space
    double timestep_;        // the time resolution of a step in seconds
    size_t gc_ID_;           // unique number for growth cones
    bool stuck_;

    // motion-related data
    double delta_angle_;
    double rw_sensing_angle_;
    double speed_growth_cone_;
    double speed_variance_;
    double average_speed_;

    Filopodia filopodia_;
    Move move_;

    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;

  public:
    GrowthCone();
    ~GrowthCone();
    //! Copy constructor for GrowthCone
    GrowthCone(const GrowthCone &copy);

    void update_topology(BaseWeakNodePtr parent, NeuritePtr ownNeurite,
                         float distanceToParent, const std::string &binaryID,
                         Point position, double angle);

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle);

    // growth
    void grow(mtPtr rnd_engine, size_t cone_n);
    void step(mtPtr rnd_engine);
    void retraction(double module);
    void prune(size_t cone_n);

    // compute direction
    void compute_pull_and_accessibility(std::vector<double> &directions_weights,
                                        mtPtr rnd_engine);
    void compute_intrinsic_direction(std::vector<double> &directions_weights);
    void choose_pull_direction(std::vector<double> &directions_weights,
                               mtPtr rnd_engine);
    virtual void compute_new_direction(mtPtr rnd_engine);

    // elongation
    void compute_module();
    virtual void compute_speed(mtPtr rnd_engine);

    double init_filopodia();
    void set_cone_ID();
    size_t get_cone_ID() const;

    void set_angle(double angle);
    virtual void prepare_for_split();
    virtual void after_split();
    virtual double compute_CR_demand(mtPtr rnd_engine);
    virtual void reset_CR_demand();

    // get functions
    double get_module() const;
    virtual double get_CR_received() const;
    virtual double get_CR_left() const;
    virtual double get_CR_used() const;
    virtual double get_growth_cone_speed() const;
    virtual double get_CR_speed_factor() const;
    virtual double get_CR_topo_coeff() const;

    // status and kernel-related functions
    virtual void set_status(const statusMap &status);
    virtual void get_status(statusMap &status) const;
    void update_kernel_variables();
};
}

#endif
