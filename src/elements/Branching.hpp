#ifndef BRANCHING_H
#define BRANCHING_H

// std c++
#include <random>
// elements includes
#include "Branch.hpp"
#include "Node.hpp"
#include "growth_names.hpp"


// libgrowth includes
#include "elements_types.hpp"
#include "spatial_types.hpp"


namespace growth
{

class Neurite;

class Branching
{
    friend class Neurite;

  private:
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    std::poisson_distribution<> poisson_;
    std::exponential_distribution<> exponential_;
    std::exponential_distribution<> exponential_uniform_;
    std::exponential_distribution<> exponential_flpl_;

    NeuritePtr neurite_;

    // CR_params
    bool use_critical_resource_;
    double CR_amount_;
    double initial_CR_amount_;
    double CR_tot_demand_;
    double CR_split_th_;
    double CR_split_tau_;

    // variables for van Pelt branching model
    bool use_van_pelt_;
    double van_pelt_norm_;
    Event next_vanpelt_event_;
    double B_;
    double E_;
    double S_;
    double T_;

    // variables for uniform branching
    bool use_uniform_branching_;
    Event next_uniform_event_;
    double uniform_branching_rate_;

    // variables for non-uniform lateral branching
    bool use_flpl_branching_;
    Event next_flpl_event_;
    double flpl_branching_rate_;

  public:
    Branching(NeuritePtr neurite);
    Branching();
    Branching(const Branching &cpy);
    // event handlers functions
    void compute_next_event(mtPtr rnd_engine);
    bool branching_event(mtPtr rnd_engine, const Event &ev);

    void update_growth_cones(mtPtr rnd_engine);

    // van Pelt branching functions
    bool vanpelt_new_branch(mtPtr rnd_engine);
    void compute_vanpelt_event(mtPtr rnd_engine);

    // uniform branching functions
    void compute_uniform_event(mtPtr rnd_engine);
    bool uniform_new_branch(mtPtr rnd_engine);

    // powerlaw branching functions
    void compute_flpl_event(mtPtr rnd_engine);
    bool flpl_new_branch(mtPtr rnd_engine);

    // critical_resource functions
    double get_CR_quotient() const;
    double get_CR_amount() const;
    void CR_new_branch(mtPtr rnd_engine, GCPtr splitting_cone);
    void initialize_next_event(mtPtr rnd_engine, double new_resolution,
                               size_t previous_step);

    // Get/set functions
    void set_status(const statusMap &);
    //~ void get_status(statusMap &) const;
    void get_status(statusMap &) const;
};
}
#endif /* BRANCHING_H */
