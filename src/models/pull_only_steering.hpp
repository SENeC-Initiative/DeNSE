#ifndef PO_STEER_H
#define PO_STEER_H

// elements
#include "steering_interface.hpp"


namespace growth
{

class PullOnlySteeringModel : public virtual SteeringModel
{
  public:
    PullOnlySteeringModel(GCPtr gc, NeuritePtr neurite) : SteeringModel(gc, neurite) {};
    PullOnlySteeringModel(const PullOnlySteeringModel &copy) = delete;
    PullOnlySteeringModel(const PullOnlySteeringModel &copy, GCPtr gc, NeuritePtr neurite)
    : SteeringModel(copy, gc, neurite) {};

    virtual void
    compute_intrinsic_direction(std::vector<double> &directions_weights,
                                const Filopodia& filo, double substep,
                                double &total_proba, bool &stuck) override final
    {
        total_proba = 0.;
        stuck       = true;

        for (double weight : directions_weights)
        {
            if (not std::isnan(weight))
            {
                total_proba += weight;
                stuck        = false;
            }
        }
    };

    virtual void set_status(const statusMap &status) override final {};

    virtual void get_status(statusMap &status) const override final {};
};

} // namespace growth

#endif
