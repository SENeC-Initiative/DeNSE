#ifndef GROWTHCONE_H
#define GROWTHCONE_H

#include <random>

#include "Neuron.hpp"
#include "Node.hpp"

// lib include
#include "config_impl.hpp"
#include "cttrie.hpp"
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
    bool using_environment_;   // whether we're embedded in space
    size_t gc_ID_;             // unique number for growth cones
    std::string current_area_; // name of the area where the GC is
    bool stuck_;
    std::vector<std::string> observables_;

    // motion-related data
    double delta_angle_;
    double sensing_angle_;
    double avg_speed_;
    double speed_variance_;
    double duration_retraction_; // duration of a retraction period (seconds)
    double max_sensing_angle_;
    double proba_retraction_; // proba of retracting when stuck
    double retracting_todo_;  // duration left to retract
    double speed_ratio_retraction_;
    double proba_down_move_; // proba of going down if bottom out of reach
    double scale_up_move_;   // maximal height that GC can cross upwards

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
                         const Point &position, double angle);

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle);

    // growth
    void grow(mtPtr rnd_engine, size_t cone_n, double substep);
    void step(mtPtr rnd_engine);
    void retraction(double module);
    void prune(size_t cone_n);

    // compute direction
    void
    compute_pull_and_accessibility(std::vector<double> &directions_weights,
                                   std::vector<std::string> &change_to_new_area,
                                   mtPtr rnd_engine, double substep);
    void compute_intrinsic_direction(std::vector<double> &directions_weights);
    int choose_pull_direction(const std::vector<double> &directions_weights,
                              const std::vector<std::string> &new_pos_area,
                              mtPtr rnd_engine);
    virtual void compute_new_direction(mtPtr rnd_engine, double substep);
    void widen_sensing_angle();
    bool nonmax_sensing_angle();

    // elongation
    void compute_module(double substep);
    virtual void compute_speed(mtPtr rnd_engine, double substep);

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
    virtual double get_state(const char *observable) const;
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
    void update_growth_properties(const std::string &area_name);
};


inline bool GrowthCone::nonmax_sensing_angle()
{
    return abs(move_.sigma_angle - max_sensing_angle_) < 1e-6;
}

inline void GrowthCone::widen_sensing_angle()
{
    move_.sigma_angle = std::min(1.5 * move_.sigma_angle, max_sensing_angle_);
}
}

#endif
