#ifndef GC_TUB_H

#define GC_TUB_H

#include "GrowthCone.hpp"
#include "growth_names.hpp"

namespace growth
{

typedef struct Demand
{
    double local_demand;
    double correlation;
    double sqrt_corr;
    double std_dev;
    double mean;
    Demand(double local_demand, double correlation, double sqrt_corr,
           double std_dev, double mean)
        : local_demand(local_demand)
        , correlation(correlation)
        , sqrt_corr(sqrt_corr)
        , std_dev(std_dev)
        , mean(mean)
    {
        printf("mean is %f, correlation is %f", mean, correlation);
    }
} Demand;


typedef struct Critical_Resource
{
    double left;
    double used;
    double received;
    double demand;
    double leakage;
    double initial_demand;
    double retraction_th;
    double elongation_th;
    double speed_factor;
    double use_ratio;

    Critical_Resource(double received, double demand, double initial_demand,
                      double leakage, double retraction_th,
                      double elongation_th, double speed_factor,
                      double use_ratio)
        : received(received)
        , demand(demand)
        , leakage(leakage)
        , initial_demand(initial_demand)
        , retraction_th(retraction_th)
        , elongation_th(elongation_th)
        , speed_factor(speed_factor)
        , use_ratio(use_ratio)
    {
    }
} Critical_Resource;


class GrowthCone_Critical : public virtual GrowthCone
{

  private:
    std::normal_distribution<double> normal_;

  protected:
    Critical_Resource critical_;
    Demand demand_;
    const Branching *neurite_dyn;

  public:
    GrowthCone_Critical();

    GrowthCone_Critical(const GrowthCone_Critical &);

    // ~GrowthCone_Critical();

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle) override;


    void initialize_CR();
    void prepare_for_split() override;
    void after_split() override;
    void reset_CR_demand() override;

    void compute_speed(mtPtr rnd_engine, double substep) override;

    double compute_CR_demand(mtPtr rnd_engine) override;
    void compute_CR_received();
    virtual void compute_CR();

    // getter functions
    void printinfo() const;
    double get_CR_received() const override;
    double get_CR_speed_factor() const override;
    double get_CR_demand() const;
    double get_CR_left() const override;
    double get_CR_used() const override;

    // status
    void set_status(const statusMap &) override;
    void get_status(statusMap &) const override;
    virtual double get_state(const char* observable) const override;
};


class GrowthCone_Critical_Langevin : public virtual GrowthCone_Critical
{
  public:
    GrowthCone_Critical_Langevin();
    GrowthCone_Critical_Langevin(const GrowthCone_Critical &);
    void compute_CR() override;
};

class GrowthCone_Critical_Gaussian : public virtual GrowthCone_Critical
{
  public:
    GrowthCone_Critical_Gaussian();
    GrowthCone_Critical_Gaussian(const GrowthCone_Critical &);
    void compute_CR() override;
};
class GrowthCone_Critical_Lurd : public virtual GrowthCone_Critical
{
  public:
    GrowthCone_Critical_Lurd();
    GrowthCone_Critical_Lurd(const GrowthCone_Critical &);
    void compute_CR() override;
};
}
#endif
