/*
 * steering_pull_only.hpp
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

#ifndef PO_STEER_H
#define PO_STEER_H

// elements
#include "steering_interface.hpp"


namespace growth
{

class PullOnlySteeringModel : public virtual SteeringModel
{
  public:
    PullOnlySteeringModel(GCPtr gc, NeuritePtr neurite) : SteeringModel(gc, neurite) {};
    PullOnlySteeringModel(const PullOnlySteeringModel &copy) = delete;
    PullOnlySteeringModel(const PullOnlySteeringModel &copy, GCPtr gc, NeuritePtr neurite)
    : SteeringModel(copy, gc, neurite) {};

    virtual void
    compute_direction_probabilities(
      std::vector<double> &directions_weights, const Filopodia& filo,
      double substep, double &total_proba, bool &stuck) override final
    {
        total_proba = 0.;
        stuck       = true;

        for (double weight : directions_weights)
        {
            if (not std::isnan(weight))
            {
                total_proba += weight;
                stuck        = false;
            }
        }
    };

    virtual void set_status(const statusMap &status) override final {};

    virtual void get_status(statusMap &status) const override final {};
};

} // namespace growth

#endif
