#ifndef GC_TUB_H

#define GC_TUB_H

#include "GrowthCone.hpp"
#include "growth_names.hpp"

namespace growth
{

class GrowthCone_Critical : public virtual GrowthCone
{
    friend class Neurite;

  private:
    std::normal_distribution<double> normal_;

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

  public:
    GrowthCone_Critical();
    GrowthCone_Critical(const GrowthCone_Critical &);

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle) override;

    void initialize_CR();
    void prepare_for_split() override;
    void after_split() override;
    void reset_CR_demand();

    void compute_speed(mtPtr rnd_engine, double substep) override;
    double compute_cr_speed(mtPtr rnd_engine, double substep);

    void compute_CR_received(double substep);
    double compute_CR(mtPtr rnd_engine, double substep);

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
