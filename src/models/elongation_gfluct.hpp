#ifndef GFLUCT_EL_H
#define GFLUCT_EL_H

// elements
#include "elongation_interface.hpp"


namespace growth
{

class GFluctElongationModel : public virtual ElongationModel
{
  protected:
    double speed_gc_avg_;
    double speed_gc_std_;

    std::normal_distribution<double> normal_;

  public:
    GFluctElongationModel(GCPtr gc, NeuritePtr neurite);
    GFluctElongationModel(const GFluctElongationModel &copy) = delete;
    GFluctElongationModel(const GFluctElongationModel &copy, GCPtr gc, NeuritePtr neurite);

    double compute_speed(mtPtr rnd_engine, double substep) override final;

    virtual void set_status(const statusMap &status) override final;

    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
