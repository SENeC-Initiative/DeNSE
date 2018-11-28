#include "nwa_direction_selector.hpp"

#include "kernel_manager.hpp"

#include "GrowthCone.hpp"


namespace growth
{

NWADirectionSelector::NWADirectionSelector(GCPtr gc, NeuritePtr neurite)
  : DirectionSelectModel(gc, neurite)
  , persistence_length_(PERSISTENCE_LENGTH)
{
    noise_amplitude_ = sqrt(2*SPEED_GROWTH_CONE / persistence_length_);
    normal_          = std::normal_distribution<double>(0, 1);
    uniform_         = std::uniform_real_distribution<double>(0., 1.);
}


NWADirectionSelector::NWADirectionSelector(const NWADirectionSelector& copy,
                                           GCPtr gc, NeuritePtr neurite)
  : DirectionSelectModel(copy, gc, neurite)
  , persistence_length_(copy.persistence_length_)
  , noise_amplitude_(copy.noise_amplitude_)
{
    normal_  = std::normal_distribution<double>(0, 1);
    uniform_ = std::uniform_real_distribution<double>(0., 1.);
}


void NWADirectionSelector::compute_target_angle(
  const std::vector<double> &directions_weights, const Filopodia &filo,
  mtPtr rnd_engine, double total_proba, bool interacting, double old_angle,
  double &substep, double &step_length, double &new_angle, bool &stopped)
{
    double weight;

    new_angle = 0.;

    // get weighted values
    for (size_t n=0; n < directions_weights.size(); n++)
    {
        weight = directions_weights[n];

        if (not std::isnan(weight))
        {
            new_angle += weight*filo.directions[n];
        }
    }

    // normalize
    new_angle /= total_proba;

    // add previous angle plus random rotation
    new_angle += old_angle +
                 noise_amplitude_*sqrt(substep)*normal_(*(rnd_engine.get()));
}


void NWADirectionSelector::set_status(const statusMap &status)
{
    double pl, na, gc_speed;
    bool bl, bn;

    // get speed (all possible versions)
    get_param(status, names::speed_growth_cone, gc_speed);

    bl = get_param(status, names::persistence_length, pl);
    bn = get_param(status, names::noise_amplitude, na);

    if (bl)
    {
        if (pl < 0)
        {
            throw std::invalid_argument(
                "`" + names::persistence_length + "` must be positive.");
        }

        persistence_length_ = pl;
    }

    if (bn)
    {
        if (na < 0)
        {
            throw std::invalid_argument(
                "`" + names::noise_amplitude + "` must be positive.");
        }

        noise_amplitude_ = na;
    }
    else
    {
        noise_amplitude_ = sqrt(2*gc_speed / persistence_length_);
    }

}


void NWADirectionSelector::get_status(statusMap &status) const
{
    set_param(status, names::persistence_length, persistence_length_, "micrometer");
    set_param(status, names::noise_amplitude, noise_amplitude_, "radian");
}

} // namespace growth
