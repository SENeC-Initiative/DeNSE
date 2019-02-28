#ifndef SRF_STEER_H
#define SRF_STEER_H

// elements
#include "steering_interface.hpp"


namespace growth
{

enum smode { window, sine };

class SrfSteeringModel : public virtual SteeringModel
{
  private:
    double rigidity_factor_;        // strength of the "forward" force
    double somatropic_factor_;      // strength of the somatropic force
    double somatropic_scale_;       // typical scale of somatropism (microns)
    double self_avoidance_factor_;  // strength of the self avoidance force
    smode somatropic_mode_;         // how somatropism is computed
    double self_avoidance_scale_;   // typical scale of self avoidance (microns)
    BPoint soma_;

  public:
    SrfSteeringModel(GCPtr gc, NeuritePtr neurite);
    SrfSteeringModel(const SrfSteeringModel &copy) = delete;
    SrfSteeringModel(const SrfSteeringModel &copy, GCPtr gc,
                     NeuritePtr neurite);

    virtual void
    compute_direction_probabilities(
        std::vector<double> &directions_weights, const Filopodia& filo,
        double substep, double &total_proba, bool &stuck) override final;

    virtual void set_status(const statusMap &status) override final;

    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
