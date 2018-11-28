#ifndef CST_EL_H
#define CST_EL_H

// elements
#include "elongation_interface.hpp"

// libgrowth
#include "config_impl.hpp"


namespace growth
{

class CstElongationModel
  : public virtual ElongationModel
  , public std::enable_shared_from_this<CstElongationModel>
{
  protected:
    double speed_growth_cone_;

  public:
    CstElongationModel(GCPtr gc, NeuritePtr neurite)
      : ElongationModel(gc, neurite), speed_growth_cone_(SPEED_GROWTH_CONE) {};

    CstElongationModel(const CstElongationModel& copy) = delete;

    CstElongationModel(const CstElongationModel& copy, GCPtr gc, NeuritePtr neurite)
      : ElongationModel(copy, gc, neurite), speed_growth_cone_(copy.speed_growth_cone_) {};

    double compute_speed(mtPtr rnd_engine, double substep) override final
    { return speed_growth_cone_; };

    virtual void set_status(const statusMap &status) override final
    { get_param(status, names::speed_growth_cone, speed_growth_cone_); };

    virtual void get_status(statusMap &status) const override final
    {
        set_param(status, names::speed_growth_cone, speed_growth_cone_,
                  "micrometer / minute");
    };
};

} // namespace growth

#endif
