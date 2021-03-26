/*
 * steering_memory_based.cpp
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

#include "steering_memory_based.hpp"

#include <cmath>

#include "config.hpp"

#include "GrowthCone.hpp"


namespace growth
{

MemBasedSteeringModel::MemBasedSteeringModel(GCPtr gc, NeuritePtr neurite)
    : SteeringModel(gc, neurite)
    , memory_angle_(fmod(gc->get_state("angle"), 2 * M_PI))
    , memory_influence_(1.)
    , memory_dist_exp_(2.)
    , memory_dist_cut_(100.)
    , dist_cut_set_(false)
{}


MemBasedSteeringModel::MemBasedSteeringModel(const MemBasedSteeringModel &copy,
                                             GCPtr gc, NeuritePtr neurite)
    : SteeringModel(copy, gc, neurite)
    , memory_angle_(fmod(gc->get_state("angle"), 2 * M_PI))
    , memory_influence_(copy.memory_influence_)
    , memory_dist_exp_(copy.memory_dist_exp_)
    , memory_dist_cut_(copy.memory_dist_cut_)
    , dist_cut_set_(false)
{}


void MemBasedSteeringModel::compute_direction_probabilities(
    std::vector<double> &directions_weights, const Filopodia &filo,
    double substep, double &total_proba, bool &stuck)
{
    double current_angle =
        fmod(gc_weakptr_.lock()->get_state("angle"), 2 * M_PI);

    stuck = true;

    // get branch data
    BranchPtr branch = gc_weakptr_.lock()->get_branch();

    stype start = branch->size() == 0 ? 0 : branch->size() - 1;

    double dmax = branch->final_distance_to_soma();
    double dmin = branch->initial_distance_to_soma();

    double taper = neurite_ptr_->get_taper_rate();
    double rinit = 0.5*gc_weakptr_.lock()->get_diameter();

    double dcumul = 0.;

    BPoint pos = gc_weakptr_.lock()->get_position();
    double x(pos.x()), y(pos.y());

    // compute the memory vector
    PointArray p;
    double radius, dist, volume, inv_norm, dx, dy;
    double xdir(0.), ydir(0.), old_dist(0.);

    if (start > 1)
    {
        for (stype i=start - 1; i > 0; i--)
        {
            p = branch->at(i);

            dist = dmax - p[2];

            radius = rinit + 0.25*taper*(dist + old_dist);
            volume = M_PI*radius*radius*(dist - old_dist);

            dx = (p[0] - x);
            dy = (p[1] - y);

            inv_norm = 1. / sqrt(dx*dx + dy*dy);

            xdir += volume / pow(0.5*(dist + old_dist), memory_dist_exp_)
                    * dx * inv_norm;

            ydir += volume / pow(0.5*(dist + old_dist), memory_dist_exp_)
                    * dy * inv_norm;

            if (dist > memory_dist_cut_)
            {
                break;
            }

            x = p[0];
            y = p[1];

            old_dist = dist;
        }
    

        // compute the angle (from target to source since the vector was
        // computed from the tip to the rear)
        memory_angle_ = std::atan2(-ydir, -xdir);

        if (memory_angle_ < 0)
        {
            memory_angle_ += 2*M_PI;
        }
    }
    else
    {
        memory_angle_ = current_angle;
    }

    // get all angles that are not NaN
    std::vector<unsigned int> valid_directions;
    std::vector<double> distances;

    double angle, adist, weight;

    for (unsigned int n = 0; n < filo.directions.size(); n++)
    {
        weight = directions_weights[n];

        if (not std::isnan(weight))
        {
            total_proba += weight;
            stuck = false;

            // store index and compute angluar distance to memory angle
            valid_directions.push_back(n);

            angle = fmod(filo.directions[n] + current_angle, 2*M_PI);

            adist = angle - memory_angle_;

            if (adist < - M_PI)
            {
                adist += 2*M_PI;
            }
            else if (adist > M_PI)
            {
                adist -= 2*M_PI;
            }

            distances.push_back(std::abs(adist));
        }
    }

    // then: find the closest angle (if not stuck)
    if (not stuck)
    {
        unsigned int idx = 
            std::min_element(distances.begin(), distances.end())
            - distances.begin();

        unsigned int chosen = valid_directions[idx];

        directions_weights[chosen] += memory_influence_;
    }
}


void MemBasedSteeringModel::set_status(const statusMap &status)
{
    double minfl, dexp, dcut;
    bool b;

    b = get_param(status, names::memory_influence, minfl);
    if (b)
    {
        if (minfl <= 0)
        {
            throw std::invalid_argument("`" + names::memory_influence +
                                        "` must be strictly positive.");
        }

        memory_influence_ = minfl;
    }

    b = get_param(status, names::memory_dist_exp, dexp);
    if (b)
    {
        memory_dist_exp_ = dexp;
    }

    b = get_param(status, names::memory_dist_cut, dcut);
    if (b)
    {
        if (dcut <= 0)
        {
            throw std::invalid_argument("`" + names::memory_dist_cut +
                                        "` must be strictly positive.");
        }

        memory_dist_cut_ = dcut;

        dist_cut_set_ = true;
    }

    if (not dist_cut_set_)
    {
        // if distance cut has not been set it default to 10^4/p
        memory_dist_cut_ = pow(10, 4/memory_dist_exp_);
    }
}


void MemBasedSteeringModel::get_status(statusMap &status) const
{
    set_param(status, names::memory_influence, memory_influence_, "");
    set_param(status, names::memory_dist_exp, memory_dist_exp_, "");
    set_param(status, names::memory_dist_cut, memory_dist_cut_,
              "micrometer");
}

} // namespace growth
