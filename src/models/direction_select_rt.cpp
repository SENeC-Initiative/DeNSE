/*
 * direction_select_rt.cpp
 *
 * This file is part of DeNSE.
 *
 * Copyright (C) 2019 SeNEC Initiative
 *
 * DeNSE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * DeNSE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DeNSE. If not, see <http://www.gnu.org/licenses/>.
 */

#include "direction_select_rt.hpp"

#include "config.hpp"

#include "kernel_manager.hpp"


namespace growth
{

RTDirectionSelector::RTDirectionSelector(GCPtr gc, NeuritePtr neurite)
  : DirectionSelectModel(gc, neurite)
  , persistence_length_(PERSISTENCE_LENGTH)
  , critical_pull_(2*FILOPODIA_WALL_AFFINITY)
  , sensing_angle_(SENSING_ANGLE)
  , p_tumble_on_stop_(0.2)
  , tumbling_(false)
  , num_tumbles_(0)
{
    observables_.push_back("num_tumbles");

    uniform_ = std::uniform_real_distribution<double>(0., 1.);

    initialize_rt();
}


RTDirectionSelector::RTDirectionSelector(const RTDirectionSelector& copy,
                                         GCPtr gc, NeuritePtr neurite)
  : DirectionSelectModel(copy, gc, neurite)
  , persistence_length_(copy.persistence_length_)
  , critical_pull_(copy.critical_pull_)
  , sensing_angle_(copy.sensing_angle_)
  , p_tumble_on_stop_(copy.p_tumble_on_stop_)
  , tumbling_(false)
  , num_tumbles_(0)
{
    uniform_ = std::uniform_real_distribution<double>(0., 1.);

    initialize_rt();
}


void RTDirectionSelector::initialize_rt()
{
    // this renormalization of the "tumbling rate" is necessary to obtain
    // the correct persistence length
    tau_ = 24. / (sensing_angle_ * sensing_angle_ * persistence_length_);

    exponential_rt_ = std::exponential_distribution<double>(tau_);

    int omp_id = kernel().parallelism_manager.get_thread_local_id();
    mtPtr rng  = kernel().rng_manager.get_rng(omp_id);

    next_tumble_ = exponential_rt_(*(rng).get());
}


void RTDirectionSelector::select_direction(
  const std::vector<double> &directions_weights, const Filopodia &filo,
  mtPtr rnd_engine, double total_proba, bool interacting, double old_angle,
  double &substep, double &step_length, double &new_angle, bool &stopped,
  size_t &default_direction)
{
    new_angle = old_angle;

    // test if out of retraction, stopped, or interacting
    if (gc_weakptr_.lock()->just_retracted())
    {
        tumbling_ = true;
    }
    else if (interacting)
    {
        // increased probability of tumbling
        // get max pul direction
        auto it_max = std::max_element(
            directions_weights.begin(), directions_weights.end());
        // get its index and the associated angle
        size_t n_max     = std::distance(directions_weights.begin(), it_max);
        double max_angle = filo.directions[n_max];

        // compute influence
        double influence = directions_weights[n_max]*abs(sin(new_angle - max_angle)) / critical_pull_;

        // test if enough to trigger tumble
        double x = uniform_(*(rnd_engine.get()));

        if (x < influence*substep)
        {
            tumbling_    = true;
        }
    }

    // if not currently tumbling, check whether we're running into a wall
    if (not tumbling_)
    {
        size_t middle = 0.5*directions_weights.size();
        double weight = directions_weights[middle];

        bool wall = std::isnan(weight) or weight == 0.;

        // if even number of filopodia, check also the other one
        if (not (directions_weights.size() % 2))
        {
            weight = directions_weights[middle+1];
            wall  *= (std::isnan(weight) or weight == 0.);
        }

        if (wall)
        {
            double x = uniform_(*(rnd_engine.get()));
            if (x < p_tumble_on_stop_*substep)
            {
                tumbling_ = true;
            }
        }
    }

    // randomly pick the angle if tumbling, otherwise it is left unchanged
    if (tumbling_)
    {
        stopped = false;
        // make a tumble (leave delta_angle_ unchanged)
        num_tumbles_++;
        // compute next tumble and reset tumbling
        next_tumble_ = exponential_rt_(*(rnd_engine).get());
        tumbling_ = false;

        // weighted random choice for the new angle
        //~ new_angle += sensing_angle_*(uniform_(*(rnd_engine.get())) - 0.5);
        double weight, cumulated_weight(0);
        double x = total_proba * uniform_(*(rnd_engine.get()));

        for (size_t i = 0; i < directions_weights.size(); i++)
        {
            weight = directions_weights[i];

            if (not std::isnan(weight))
            {
                cumulated_weight += weight;

                if (x < cumulated_weight)
                {
                    new_angle += filo.directions[i];
                    break;
                }
            }
        }

        // default angle is closest to new_angle
        double dist, min_dist(std::numeric_limits<double>::max());

        for (size_t n=0; n < directions_weights.size(); n++)
        {
            if (not std::isnan(directions_weights[n]))
            {
                dist = std::abs(new_angle - filo.directions[n]);
                if (dist < min_dist)
                {
                    default_direction = n;
                    min_dist          = dist;
                }
            }
        }
    }
    else
    {
        // keep straight: default angle is closest to zero
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
    }

    // test whether a new tumble will happen before the end of the substep
    if (next_tumble_ - step_length <= 0)
    {
        // we will tumble before the full module, change distance done
        double distance_done = next_tumble_;
        // compute how much time remains until tumble and update the current
        // substep accordingly.
        substep = substep * (distance_done / step_length);
        // we are still running, so set delta_angle_ to zero

        // IMPORTANT: update step_length and set tumbling for next substep
        step_length = distance_done;
        tumbling_   = true;
    }
    else
    {
        // no spontaneous tumble during substep, keep straight
        next_tumble_ -= step_length;
    }
}


double RTDirectionSelector::get_state(const char *observable) const
{
    double value = 0.;

    TRIE(observable)
    CASE("num_tumbles")
    value = num_tumbles_;
    ENDTRIE;

    return value;
}


void RTDirectionSelector::set_status(const statusMap &status)
{
    double pl, cp;
    bool b;

    // this one was already checked in GrowthCone, no need to double-check
    get_param(status, names::sensing_angle, sensing_angle_);

    // persistence length
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

    // critical pull
    b = get_param(status, names::critical_pull, cp);
    if (b)
    {
        if (cp < 0)
        {
            throw std::invalid_argument(
                "`" + names::critical_pull + "` must be positive.");
        }

        critical_pull_ = cp;
    }

    initialize_rt();
}


void RTDirectionSelector::get_status(statusMap &status) const
{
    set_param(status, names::persistence_length, persistence_length_, "micrometer");
    set_param(status, names::critical_pull, critical_pull_, "");
}

} // namespace growth
