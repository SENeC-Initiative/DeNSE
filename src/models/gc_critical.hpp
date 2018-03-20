#ifndef GC_TUB_H

#define GC_TUB_H

#include "GrowthCone.hpp"
#include "growth_names.hpp"

namespace growth
{


typedef struct CR_Demand
{
    double weight;
    double consumption_rate;
    double demand;
} CR_Demand;

typedef struct Critical_Speed
{
    double elongation_factor;
    double elongation_th;
    double retraction_factor;
    double retraction_th;
} CR_Speed;

typedef struct Critical_Resource
{
    double stochastic_tmp;
    double received;
    double sqrt_corr;
    double stored;
    double noise;

    double leakage;
    double use_ratio;
    double correlation;
    double variance;


    //Critical_Resource(double received, double demand, double initial_demand,
                      //double leakage, double retraction_th,
                      //double elongation_th, double speed_factor,
                      //double use_ratio)
        //: received(received)
        //, demand(demand)
        //, leakage(leakage)
        //, initial_demand(initial_demand)
        //, retraction_th(retraction_th)
        //, elongation_th(elongation_th)
        //, speed_factor(speed_factor)
        //, use_ratio(use_ratio)
    //{
    //}
} Critical_Resource;


class GrowthCone_Critical : public virtual GrowthCone
{

  private:
    std::normal_distribution<double> normal_;

  protected:
    Critical_Resource critical_;
    CR_Demand demand_;
    CR_Speed cr_speed_;
    const Branching *neurite_dyn;
    double timestep_;
    double euler_step_;

  public:
    GrowthCone_Critical();
    GrowthCone_Critical(const GrowthCone_Critical &);

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle) override;

    void initialize_CR();
    void prepare_for_split() override;
    void after_split() override;
    void reset_CR_demand() override;

    void compute_speed(mtPtr rnd_engine, double substep) override;

    void compute_CR_demand(mtPtr rnd_engine) override;
    void compute_CR_received();
    void compute_CR(mtPtr rnd_engine);

    // getter functions
    void printinfo() const;
    double get_CR_received() const override;
    double get_CR_speed_factor() const override;
    double get_CR_demand() const override;
    double get_speed() const;
    double CR_differential(double a, double k,  double rate, double received, double dt);

    // status
    void set_status(const statusMap &) override;
    void get_status(statusMap &) const override;
    virtual double get_state(const char *observable) const override;
};

}
#endif
