#include "nm_direction_selector.hpp"

#include "kernel_manager.hpp"

#include "GrowthCone.hpp"


namespace growth
{

NMDirectionSelector::NMDirectionSelector(GCPtr gc, NeuritePtr neurite)
  : DirectionSelectModel(gc, neurite)
  , persistence_length_(PERSISTENCE_LENGTH)
{
    noise_amplitude_ = sqrt(2*SPEED_GROWTH_CONE / persistence_length_);
    normal_          = std::normal_distribution<double>(0, 1);
    uniform_         = std::uniform_real_distribution<double>(0., 1.);
}


NMDirectionSelector::NMDirectionSelector(const NMDirectionSelector& copy,
                                           GCPtr gc, NeuritePtr neurite)
  : DirectionSelectModel(copy, gc, neurite)
  , persistence_length_(copy.persistence_length_)
  , noise_amplitude_(copy.noise_amplitude_)
{
    normal_  = std::normal_distribution<double>(0, 1);
    uniform_ = std::uniform_real_distribution<double>(0., 1.);
}


void NMDirectionSelector::compute_target_angle(
  const std::vector<double> &directions_weights, const Filopodia &filo,
  mtPtr rnd_engine, double total_proba, bool interacting, double old_angle,
  double &substep, double &step_length, double &new_angle, bool &stopped)
{
    auto it = std::adjacent_find(
        directions_weights.begin(), directions_weights.end(),
        std::not_equal_to<double>());

    // get max only if some weights differ
    if (it != directions_weights.end())
    {
        auto it_max  = std::max_element(
            directions_weights.begin(), directions_weights.end());
        // get its index and the associated angle
        double n_max = std::distance(directions_weights.begin(), it_max);

        // set angle + add previous angle plus random rotation
        new_angle    =
            filo.directions[n_max] + old_angle +
            noise_amplitude_*sqrt(substep)*normal_(*(rnd_engine.get()));
    }
    else
    {
        // "keep straight"
        new_angle = old_angle +
                    noise_amplitude_*sqrt(substep)*normal_(*(rnd_engine.get()));
    }
}


void NMDirectionSelector::set_status(const statusMap &status)
{
    double pl, na;
    bool b;

    b = get_param(status, names::persistence_length, pl);
    if (b)
    {
        if (pl < 0)
        {
            throw std::invalid_argument(
                "`" + names::persistence_length + "` must be positive.");
        }

        persistence_length_ = pl;
    }

    b = get_param(status, names::noise_amplitude, na);
    if (b)
    {
        if (na < 0)
        {
            throw std::invalid_argument(
                "`" + names::noise_amplitude + "` must be positive.");
        }

        noise_amplitude_ = na;
    }
}


void NMDirectionSelector::get_status(statusMap &status) const
{
    set_param(status, names::persistence_length, persistence_length_, "micrometer");
    set_param(status, names::noise_amplitude, noise_amplitude_, "radian");
}

} // namespace growth
