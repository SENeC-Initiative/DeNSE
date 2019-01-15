#ifndef CST_EL_H
#define CST_EL_H

// elements
#include "extension_interface.hpp"

// libgrowth
#include "config_impl.hpp"


namespace growth
{

class CstExtensionModel
  : public virtual ExtensionModel
{
  protected:
    double speed_growth_cone_;

  public:
    CstExtensionModel(GCPtr gc, NeuritePtr neurite)
      : ExtensionModel(gc, neurite), speed_growth_cone_(SPEED_GROWTH_CONE) {};

    CstExtensionModel(const CstExtensionModel& copy) = delete;

    CstExtensionModel(const CstExtensionModel& copy, GCPtr gc, NeuritePtr neurite)
      : ExtensionModel(copy, gc, neurite), speed_growth_cone_(copy.speed_growth_cone_) {};

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
