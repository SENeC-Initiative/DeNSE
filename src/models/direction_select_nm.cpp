#include "direction_select_nm.hpp"

// kernel includes
#include "kernel_manager.hpp"

// element includes
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


void NMDirectionSelector::select_direction(
  const std::vector<double> &directions_weights, const Filopodia &filo,
  mtPtr rnd_engine, double total_proba, bool interacting, double old_angle,
  double &substep, double &step_length, double &new_angle, bool &stopped,
  size_t &default_direction)
{
    auto it = std::adjacent_find(
        directions_weights.begin(), directions_weights.end(),
        std::not_equal_to<double>());

    // get max only if some weights differ
    if (it != directions_weights.end())
    {
        double w, w_max = std::numeric_limits<double>::lowest();

        for (size_t n = 0; n < directions_weights.size(); n++)
        {
            w = directions_weights[n];

            if (not std::isnan(w) and w > w_max)
            {
                w_max             = w;
                default_direction = n;
            }
        }

        // set angle + add previous angle plus random rotation
        new_angle    =
            filo.directions[default_direction] + old_angle +
            noise_amplitude_*sqrt(substep)*normal_(*(rnd_engine.get()));
    }
    else
    {
        // "keep straight"
        // default angle is closest to zero
        double dist, min_dist(std::numeric_limits<double>::max());

        for (size_t n=0; n < directions_weights.size(); n++)
        {
            if (not std::isnan(directions_weights[n]))
            {
                dist = std::abs(filo.directions[n]);
                if (dist < min_dist)
                {
                    default_direction = n;
                    min_dist          = dist;
                }
            }
        }

        new_angle = old_angle +
                    noise_amplitude_*sqrt(substep)*normal_(*(rnd_engine.get()));
    }
}


void NMDirectionSelector::set_status(const statusMap &status)
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


void NMDirectionSelector::get_status(statusMap &status) const
{
    set_param(status, names::persistence_length, persistence_length_, "micrometer");
    set_param(status, names::noise_amplitude, noise_amplitude_, "radian");
}

} // namespace growth
