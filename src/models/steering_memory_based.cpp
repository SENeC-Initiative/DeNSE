#include "steering_memory_based.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

#include "config.hpp"

#include "GrowthCone.hpp"


namespace growth
{

MemBasedSteeringModel::MemBasedSteeringModel(GCPtr gc, NeuritePtr neurite)
  : SteeringModel(gc, neurite)
  , memory_angle_(fmod(gc->get_state("angle"), 2*M_PI))
  , rigidity_factor_(1.)
  , decay_factor_(0.9)
{
    if (memory_angle_ < -M_PI)
    {
        memory_angle_ += 2*M_PI;
    }
    else if (memory_angle_ > M_PI)
    {
        memory_angle_ -= 2*M_PI;
    }
}


MemBasedSteeringModel::MemBasedSteeringModel(const MemBasedSteeringModel &copy, GCPtr gc, NeuritePtr neurite)
  : SteeringModel(copy, gc, neurite)
  , memory_angle_(fmod(gc->get_state("angle"), 2*M_PI))
  , rigidity_factor_(copy.rigidity_factor_)
  , decay_factor_(copy.decay_factor_)
{
    if (memory_angle_ < -M_PI)
    {
        memory_angle_ += 2*M_PI;
    }
    else if (memory_angle_ > M_PI)
    {
        memory_angle_ -= 2*M_PI;
    }
}


void MemBasedSteeringModel::compute_direction_probabilities(
    std::vector<double> &directions_weights, const Filopodia& filo, double substep,
    double &total_proba, bool &stuck)
{
    double weight, angle, previous_angle, da, previous_da;

    stuck       = true;
    total_proba = 0.;

    double current_angle = fmod(gc_weakptr_.lock()->get_state("angle"), 2*M_PI);

    // update memory angle:
    // - first get the last segment's length
    double rodlen = gc_weakptr_.lock()->get_branch()->get_last_segment_length();
    // - second, get mean radius from taper rate (current radius is smaller end)
    double taper  = neurite_ptr_->get_taper_rate();
    double radius = 0.5*(gc_weakptr_.lock()->get_diameter() + 0.5*rodlen*taper);
    // - third, compute segment volume
    double volume = M_PI*radius*radius*rodlen;
    // - then update through volume-weighted algorithm
    double decay  = pow(decay_factor_, rodlen);
    memory_angle_ = (volume*current_angle + decay*memory_angle_)
                    / (volume + decay);

    // loop over the angles, get total probability, and add memory contribution
    size_t n_max = filo.directions.size() - 1;
    bool do_test = true;

    for (unsigned int n = 0; n < filo.directions.size(); n++)
    {
        angle  = filo.directions[n] + current_angle;

        // check for the direction of the memory angle
        if (do_test and n == 0 and angle > memory_angle_)
        {
            da                    = angle - memory_angle_;
            directions_weights[0] = std::max(
                directions_weights[0] + cos(da) * rigidity_factor_, 0.);

            do_test = false;
        }
        else if (do_test and angle > memory_angle_)
        {
            previous_angle = filo.directions[n - 1];
            da             = angle - memory_angle_;
            previous_da    = previous_angle - memory_angle_;

            directions_weights[n] = std::max(
                directions_weights[n] + da / std::max(-previous_da, da)
                * cos(da) * rigidity_factor_, 0.);

            directions_weights[n-1] = std::max(
                directions_weights[n-1] - previous_da
                / std::max(-previous_da, da) * cos(previous_da)
                * rigidity_factor_, 0.);

            do_test = false;
        }
        else if (n == n_max and angle < memory_angle_)
        {
            da                        = memory_angle_ - angle;
            directions_weights[n_max] = std::max(
                directions_weights[n_max] + cos(da) * rigidity_factor_, 0.);
        }

        weight = directions_weights[n];

        if (not std::isnan(weight))
        {
            total_proba += weight;
            stuck        = false;
        }
    }

    if (memory_angle_ < -2*M_PI)
    {
        memory_angle_ += 2*M_PI;
    }
    else if (memory_angle_ > 2*M_PI)
    {
        memory_angle_ -= 2*M_PI;
    }
}


void MemBasedSteeringModel::set_status(const statusMap &status)
{
    double rf, md;
    bool b;

    b = get_param(status, names::rigidity_factor, rf);
    if (b)
    {
        if (rf < 0)
        {
            throw std::invalid_argument(
                "`" + names::rigidity_factor + "` must be positive.");
        }

        rigidity_factor_ = rf;
    }

    b = get_param(status, names::decay_factor, md);
    if (b)
    {
        if (md < 0)
        {
            throw std::invalid_argument(
                "`" + names::decay_factor + "` must be positive.");
        }

        decay_factor_ = md;
    }
}


void MemBasedSteeringModel::get_status(statusMap &status) const
{
    set_param(status, names::rigidity_factor, rigidity_factor_, "");
    set_param(status, names::decay_factor, decay_factor_, "");
}

}
