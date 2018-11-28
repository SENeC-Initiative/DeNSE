#ifndef RT_DIRSEL_H
#define RT_DIRSEL_H

#include "direction_select_interface.hpp"


namespace growth
{

class RTDirectionSelector: public virtual DirectionSelectModel
{
  private:
    double persistence_length_;  // typical distance on which gc goes straight
    double critical_pull_;       // pull necessary to trigger tumble with proba 1 if applied orthogonally
    double sensing_angle_;       // duplicate of GC entry to convenience
    double tau_;                 // tumbling rate
    double next_tumble_;         // distance left to next tumble
    double tumbling_;            // whether a tumble is happening
    double p_tumble_on_stop_;    // proba of tumbling if stopped
    size_t num_tumbles_;         // keep track of tumbling occurrences

    std::exponential_distribution<double> exponential_rt_;
    std::uniform_real_distribution<double> uniform_;
    
  public:
    RTDirectionSelector() = delete;
    RTDirectionSelector(GCPtr gc, NeuritePtr neurite);
    RTDirectionSelector(const RTDirectionSelector& copy) = delete;
    RTDirectionSelector(const RTDirectionSelector& copy, GCPtr gc, NeuritePtr neurite);

    virtual void
    compute_target_angle(
        const std::vector<double> &directions_weights, const Filopodia &filo,
        mtPtr rnd_engine, double total_proba, bool interacting,
        double old_angle, double &substep, double &step_length,
        double &new_angle, bool &stopped) override final;

    void initialize_rt();

    virtual double get_state(const char *observable) const override final;
    virtual void set_status(const statusMap &status) override final;
    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
