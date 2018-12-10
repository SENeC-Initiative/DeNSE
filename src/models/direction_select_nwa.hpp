#ifndef NWA_DIRSEL_H
#define NWA_DIRSEL_H

#include "direction_select_interface.hpp"


namespace growth
{

class NWADirectionSelector: public virtual DirectionSelectModel
{
  private:
    double persistence_length_;
    double noise_amplitude_;

    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    
  public:
    NWADirectionSelector() = delete;
    NWADirectionSelector(const NWADirectionSelector& copy) = delete;

    NWADirectionSelector(GCPtr gc, NeuritePtr neurite);
    NWADirectionSelector(const NWADirectionSelector& copy, GCPtr gc, NeuritePtr neurite);

    virtual void
    select_direction(
        const std::vector<double> &directions_weights, const Filopodia &filo,
        mtPtr rnd_engine, double total_proba, bool interacting,
        double old_angle, double &substep, double &step_length,
        double &new_angle, bool &stopped,
        size_t &default_direction) override final;

    virtual void set_status(const statusMap &status) override final;
    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
