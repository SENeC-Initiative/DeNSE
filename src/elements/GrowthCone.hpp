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
    double sqrt_resol_;
    size_t gc_ID_;             // unique number for growth cones
    std::string current_area_; // name of the area where the GC is
    bool stuck_;
    bool stopped_;
    bool update_filopodia_;
    std::vector<std::string> observables_;

    // motion-related data
    double delta_angle_;
    double sensing_angle_;
    double avg_speed_;
    double local_avg_speed_;
    double speed_variance_;
    double local_speed_variance_;
    double duration_retraction_; // duration of a retraction period (seconds)
    double proba_retraction_; // proba of retracting when stuck
    double retracting_todo_;  // duration left to retract
    double speed_ratio_retraction_;
    double proba_down_move_; // proba of going down if bottom out of reach
    double scale_up_move_;   // maximal height that GC can cross upwards
    double retraction_time_;

    double base_sensing_angle_;
    double max_sensing_angle_;
    size_t min_filopodia_;   // minimal number of filopodia
    size_t num_filopodia_;   // minimal number of filopodia

    Filopodia filopodia_;
    Move move_;

    double total_proba_;  // integrated probability of all possible moves
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    std::exponential_distribution<double> exponential_;

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
    void retraction();
    void prune(size_t cone_n);

    // compute direction
    bool compute_pull(std::vector<double> &directions_weights,
                        std::vector<bool> &wall_presence, double substep,
                        mtPtr rnd_engine);
    void compute_accessibility(std::vector<double> &directions_weights,
                               std::vector<std::string> &change_to_new_area);
    void compute_intrinsic_direction(std::vector<double> &directions_weights,
                                     double substep);
    void make_move(const std::vector<double> &directions_weights,
                   const std::vector<std::string> &new_pos_area,
                   double substep, mtPtr rnd_engine, int omp_id);
    virtual Point compute_new_position(
        const std::vector<double> &directions_weights, mtPtr rnd_engine,
        double substep, double frac, int n, int omp_id);
    double check_retraction(double substep, mtPtr rnd_engine);
    void change_sensing_angle(double angle);

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
    void update_filopodia(double substep);
};

}

#endif
