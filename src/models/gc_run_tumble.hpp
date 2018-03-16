#ifndef GC_RT_H

#define GC_RT_H
#include "GrowthCone.hpp"
#include "elements_types.hpp"
#include "kernel_manager.hpp"
namespace growth
{

class GrowthCone_RunTumble : public virtual GrowthCone
{
    // typedef void (GrowthCone_RunTumble::*model_function)(void);

  public:
    // public constructor required for models (consider moving to class?)
    GrowthCone_RunTumble();
    GrowthCone_RunTumble(const GrowthCone_RunTumble &);

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, std::string binaryID,
                        const Point &position, double angle) override;

    void prepare_for_split() override;
    void after_split() override;
    void initialize_RT();

    virtual Point compute_target_position(
        const std::vector<double> &directions_weights, mtPtr rnd_engine,
        double &substep, double &new_angle) override;

    virtual void set_status(const statusMap &status) override;
    virtual void get_status(statusMap &status) const override;
    
  private:
    double deterministic_angle_;
    double persistence_length_;
    double tau_;
    double next_tumble_;
    std::exponential_distribution<double> exponential_;
};
}
#endif
