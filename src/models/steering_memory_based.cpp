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
    , memory_dist_exp_(1.2)
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
    // get branch data
    BranchPtr branch = gc_weakptr_.lock()->get_branch();

    stype start = branch->size();

    double dmax = branch->final_distance_to_soma();
    double dmin = branch->initial_distance_to_soma();

    double taper = gc_weakptr_.lock()->get_taper_rate();
    double rinit = 0.5*gc_weakptr_.lock()->get_diameter();

    double dcumul = 0.;

    BPoint pos = gc_weakptr_.lock()->get_position();
    double x(pos.x()), y(pos.y());

    // compute the memory vector
    PointArray p;
    double radius, dist, xdir, ydir, volume, old_dist(0.);

    if (start > 1)
    {
        for (stype i=start - 2; i > 0; i--)
        {
            p = branch->at(i);

            dist = dmax - p[2];

            radius = rinit + taper*dist;
            volume = M_PI*radius*radius*(dist - old_dist);

            xdir += volume / power(dist, memory_dist_exp_) * p[0];
            ydir += volume / power(dist, memory_dist_exp_) * p[1];

            if (dist > memory_dist_cut_)
            {
                break;
            }
        }
    }

    // compute the angle (from target to source since the vector was computed
    // from the tip to the rear)
    memory_angle_ = std::atan2(-ydir, -x) + M_PI;
}


void MemBasedSteeringModel::set_status(const statusMap &status)
{
    double minfl, dexp, dcut;
    bool b;

    b = get_param(status, names::memory_dist_exp_, minfl);
    if (b)
    {
        if (minfl <= 0)
        {
            throw std::invalid_argument("`" + names::memory_dist_exp +
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
        // if distance cut has not been set it defaults to persistence length
        double plength;

        b = get_param(status, names::persistence_length, plength);

        if (b)
        {
            memory_dist_cut_ = 2*plength;
        }
        else
        {
            statusMap status;
            memory_dist_cut_ = gc_weakptr_.lock()->get_status(status);

            get_param(status, names::persistence_length, plength);
            memory_dist_cut_ = 2*plength;
        }
    }
}


void MemBasedSteeringModel::get_status(statusMap &status) const
{
    set_param(status, names::rigidity_factor, rigidity_factor_, "");
    set_param(status, names::memory_decay_factor, memory_decay_factor_, "");
}

} // namespace growth
