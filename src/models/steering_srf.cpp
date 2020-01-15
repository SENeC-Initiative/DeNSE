/*
 * steering_srf.cpp
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

#include "steering_srf.hpp"

// standard includes
#define _USE_MATH_DEFINES
#include <cmath>

// lib includes
#include "config.hpp"
#include "exceptions.hpp"

// kernel includes
#include "kernel_manager.hpp"

// element includes
#include "GrowthCone.hpp"


namespace growth
{

SrfSteeringModel::SrfSteeringModel(GCPtr gc, NeuritePtr neurite)
    : SteeringModel(gc, neurite)
    , rigidity_factor_(SRF_RIGIDITY_FACTOR)
    , somatropic_factor_(SRF_SOMATROPIC_FACTOR)
    , somatropic_scale_(SRF_SOMATROPIC_SCALE)
    , somatropic_mode_(sine)
    , self_avoidance_factor_(SRF_AVOIDANCE_FACTOR)
    , self_avoidance_scale_(SRF_AVOIDANCE_SCALE)
    , soma_()
{
}


SrfSteeringModel::SrfSteeringModel(const SrfSteeringModel &copy, GCPtr gc,
                                   NeuritePtr neurite)
    : SteeringModel(copy, gc, neurite)
    , rigidity_factor_(copy.rigidity_factor_)
    , somatropic_factor_(copy.somatropic_factor_)
    , somatropic_scale_(copy.somatropic_scale_)
    , somatropic_mode_(copy.somatropic_mode_)
    , self_avoidance_factor_(copy.self_avoidance_factor_)
    , self_avoidance_scale_(copy.self_avoidance_scale_)
    , soma_(neurite->get_parent_neuron().lock()->get_position())
{
}


void SrfSteeringModel::compute_direction_probabilities(
    std::vector<double> &directions_weights, const Filopodia &filo,
    double substep, double &total_proba, bool &stuck)
{
    stuck       = true;
    total_proba = 0.;

    // in this scheme, the notion of "rigidity", modeling the tendency of the
    // growth cone to continue in a straight lineis modeled by a Gaussian
    // modulation of the affinities of each of the filopodia.
    // This is because a) modeling it as a force is wrong b) this should not
    // modify the total probability, i.e. the probability of extending or
    // retracting.
    double sigma = 0.5 * (filo.directions.front() - filo.directions.back());
    double gexp_factor = rigidity_factor_ / (2 * sigma * sigma * substep);

    // the somatropism is computed as a defavorable influence which reduces the
    // probability of going towards the soma, this is given by a the power of
    // a cosine function: cos(0.5*(\theta_s - \theta))^somatropic_factor
    double current_angle =
        fmod(gc_weakptr_.lock()->get_state("angle"), 2 * M_PI);
    if (current_angle < 0)
    {
        current_angle += 2 * M_PI;
    }

    BPoint current_pos = gc_weakptr_.lock()->get_position();

    double dx(current_pos.x() - soma_.x()), dy(current_pos.y() - soma_.y());

    double somatropic_angle = atan2(dy, dx);
    if (somatropic_angle < 0)
    {
        somatropic_angle += 2 * M_PI;
    }

    double soma_distance = sqrt(dx * dx + dy * dy);
    double sexp_factor   = exp(-soma_distance / somatropic_scale_);
    double spower        = somatropic_factor_ * substep;

    // the interaction with the neighbors will be tested inside the loop
    std::vector<ObjectInfo> neighbors_info;
    std::vector<BPolygonPtr> neighbors;
    BPolygonPtr neighbor;
    ObjectInfo info;

    double max_dist     = 5 * self_avoidance_scale_;
    stype neuron_id     = neurite_ptr_->get_parent_neuron().lock()->get_gid();
    std::string neurite = neurite_ptr_->get_name();
    stype gc_id         = gc_weakptr_.lock()->get_node_id();
    double radius       = 0.5 * gc_weakptr_.lock()->get_diameter();

    // loop over the angles, get total probability, and add memory contribution
    double weight, abs_angle, angle, gfactor, sangle, sfactor, distance, tmp;
    double inv_pi = 1. / M_PI;
    std::string nneurite;
    stype nneuron, ngc;
    BPoint target_pos;

    for (unsigned int n = 0; n < filo.directions.size(); n++)
    {
        weight = directions_weights[n];

        if (not std::isnan(weight))
        {
            angle = filo.directions[n];

            // compute the rigidity effect (multiplication by Gaussian factor)
            gfactor = exp(-gexp_factor * angle * angle);
            weight *= gfactor;

            // compute the somatropic effect
            // 1) multiplication by window:
            // filopodia sense negative gradients, i.e. direction probability
            // decreases more or less sharply as we pass the 90° limit and
            // start going "towards" the soma
            // 2) smooth effect via sine function
            abs_angle = fmod(current_angle + angle, 2 * M_PI);

            if (abs_angle < 0)
            {
                abs_angle += 2 * M_PI;
            }

            sangle = std::abs(somatropic_angle - abs_angle);
            // printf("sangle %f\n", sangle*180./M_PI);
            switch (somatropic_mode_)
            {
            case window:
                sfactor = (sangle <= 0.5 * M_PI)
                              ? 1
                              : pow(1 - (2 * sangle * inv_pi - 1) * sexp_factor,
                                    spower);
                break;
            case sine:
                sfactor = pow(1 - sin(0.5 * sangle) * sexp_factor, spower);
                break;
            }

            weight *= sfactor;

            // compute the interactions with the neighbors
            target_pos = BPoint(current_pos.x() + max_dist * cos(abs_angle),
                                current_pos.y() + max_dist * sin(abs_angle));

            kernel().space_manager.get_intersected_objects(
                current_pos, target_pos, neighbors_info, neighbors);

            distance = std::nan("");

            for (stype i = 0; i < neighbors_info.size(); i++)
            {
                info     = neighbors_info[i];
                neighbor = neighbors[i];
                nneuron  = std::get<0>(info);
                nneurite = std::get<1>(info);
                ngc      = std::get<2>(info);

                // SRF avoidance is only interactions with other growth cones
                // //// of the same neuron
                if (nneuron == neuron_id and not nneurite.empty() and
                    not(nneurite == neurite and ngc == gc_id))
                {
                    // get distance, note that contact or close-to-contact
                    // interactions are ignored here, since they are handled by
                    // the affinities
                    tmp = bg::distance(current_pos, *(neighbor.get()));

                    if (std::isnan(distance) and tmp > radius)
                    {
                        distance = tmp - radius;
                    }
                    else if (tmp > radius)
                    {
                        distance = std::min(tmp - radius, distance);
                    }
                }
            }

            if (not std::isnan(distance))
            {
                weight *=
                    std::max(1 - self_avoidance_factor_ *
                                     exp(-distance / self_avoidance_scale_),
                             0.);
            }

            directions_weights[n] = weight;

            if (not std::isnan(weight))
            {
                stuck = false;
                total_proba += weight;
            }
        }
    }
}


void SrfSteeringModel::set_status(const statusMap &status)
{
    get_param(status, names::rigidity_factor, rigidity_factor_);
    get_param(status, names::somatropic_factor, somatropic_factor_);
    get_param(status, names::somatropic_scale, somatropic_scale_);

    bool has_mode;
    std::string mode;
    has_mode = get_param(status, names::somatropic_mode, mode);

    if (has_mode)
    {
        if (mode == "window")
        {
            somatropic_mode_ = window;
        }
        else if (mode == "sine")
        {
            somatropic_mode_ = sine;
        }
        else
        {
            throw InvalidArg("`somatropic_mode` must be either 'window' or "
                             "'sine' (got '" +
                                 mode + "').",
                             __FUNCTION__, __FILE__, __LINE__);
        }
    }

    get_param(status, names::self_avoidance_factor, self_avoidance_factor_);
    get_param(status, names::self_avoidance_scale, self_avoidance_scale_);
}


void SrfSteeringModel::get_status(statusMap &status) const
{
    set_param(status, names::rigidity_factor, rigidity_factor_, "");
    set_param(status, names::somatropic_factor, somatropic_factor_, "");
    set_param(status, names::somatropic_scale, somatropic_scale_, "micrometer");

    std::string mode = (somatropic_mode_ == window) ? "window" : "sine";
    set_param(status, names::somatropic_mode, mode, "");

    set_param(status, names::self_avoidance_factor, self_avoidance_factor_, "");
    set_param(status, names::self_avoidance_scale, self_avoidance_scale_,
              "micrometer");
}

} // namespace growth
