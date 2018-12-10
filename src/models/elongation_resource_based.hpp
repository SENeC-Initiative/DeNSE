#ifndef RB_EL_H
#define RB_EL_H

#include "elongation_interface.hpp"


namespace growth
{

class ResourceBasedElongationModel : public virtual ElongationModel
{
    friend class Neurite;

  protected:
    // parameters
    double leakage_;
    double use_ratio_;
    double consumption_rate_;

    // resource
    double stochastic_tmp_;
    double received_;
    double stored_;

    // noise
    double variance_;
    double correlation_;
    double sqrt_corr_;
    double noise_;

    // weights
    double weight_diameter_;
    double weight_centrifugal_;

    // speed
    double elongation_factor_;
    double elongation_th_;
    double retraction_factor_;
    double retraction_th_;
    double branching_th_;
    double branching_proba_;

    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;

  public:
    ResourceBasedElongationModel(GCPtr gc, NeuritePtr neurite);
    ResourceBasedElongationModel(const ResourceBasedElongationModel &) = delete;
    ResourceBasedElongationModel(const ResourceBasedElongationModel &copy, GCPtr gc, NeuritePtr neurite);

    void initialize_CR();
    void prepare_for_split() override;
    void after_split() override;
    void reset_CR_demand();

    double compute_speed(mtPtr rnd_engine, double substep) override;

    void compute_CR_received(double substep);
    double compute_CR(mtPtr rnd_engine, double substep, double step_length,
                      bool stuck);

    // getter functions
    void printinfo() const;
    double get_CR_demand();
    double get_CR_received() const;
    double get_CR_speed_factor() const;
    double get_speed() const;

    // status
    void set_status(const statusMap &) override;
    void get_status(statusMap &) const override;
    virtual double get_state(const char *observable) const override;
};

} // namespace growth
#endif
