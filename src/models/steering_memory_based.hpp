/*
 * steering_memory_based.hpp
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

#ifndef MEM_STEER_H
#define MEM_STEER_H

// elements
#include "steering_interface.hpp"


namespace growth
{

class MemBasedSteeringModel : public virtual SteeringModel
{
  private:
    double memory_angle_;    // priviledged direction
    double rigidity_factor_; // "strength" of the memory's influence
    double decay_factor_;    // decay of a segment's influence after 1 um

  public:
    MemBasedSteeringModel(GCPtr gc, NeuritePtr neurite);
    MemBasedSteeringModel(const MemBasedSteeringModel &copy) = delete;
    MemBasedSteeringModel(const MemBasedSteeringModel &copy, GCPtr gc,
                          NeuritePtr neurite);

    virtual void compute_direction_probabilities(
        std::vector<double> &directions_weights, const Filopodia &filo,
        double substep, double &total_proba, bool &stuck) override final;

    virtual void set_status(const statusMap &status) override final;

    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
