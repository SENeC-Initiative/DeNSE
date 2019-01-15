#ifndef GFLUCT_EL_H
#define GFLUCT_EL_H

// elements
#include "extension_interface.hpp"


namespace growth
{

class GFluctExtensionModel : public virtual ExtensionModel
{
  protected:
    double speed_gc_avg_;
    double speed_gc_std_;

    std::normal_distribution<double> normal_;

  public:
    GFluctExtensionModel(GCPtr gc, NeuritePtr neurite);
    GFluctExtensionModel(const GFluctExtensionModel &copy) = delete;
    GFluctExtensionModel(const GFluctExtensionModel &copy, GCPtr gc, NeuritePtr neurite);

    double compute_speed(mtPtr rnd_engine, double substep) override final;

    virtual void set_status(const statusMap &status) override final;

    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
