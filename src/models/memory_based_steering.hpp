#ifndef MEM_STEER_H
#define MEM_STEER_H

// elements
#include "steering_interface.hpp"


namespace growth
{

class MemBasedSteeringModel : public virtual SteeringModel
{
  private:
    double memory_angle_;     // priviledged direction
    double rigidity_factor_;  // "strength" of the memory's influence
    double decay_factor_;     // decay of a segment's influence after 1 um

  public:
    MemBasedSteeringModel(GCPtr gc, NeuritePtr neurite);
    MemBasedSteeringModel(const MemBasedSteeringModel &copy) = delete;
    MemBasedSteeringModel(const MemBasedSteeringModel &copy, GCPtr gc, NeuritePtr neurite);

    virtual void
    compute_intrinsic_direction(std::vector<double> &directions_weights,
                                const Filopodia& filo, double substep,
                                double &total_proba, bool &stuck) override final;

    virtual void set_status(const statusMap &status) override final;

    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
