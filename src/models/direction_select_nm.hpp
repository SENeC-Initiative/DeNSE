/*
 * direction_select_nm.hpp
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

#ifndef NM_DIRSEL_H
#define NM_DIRSEL_H

#include "direction_select_interface.hpp"


namespace growth
{

class NMDirectionSelector: public virtual DirectionSelectModel
{
  private:
    double persistence_length_;
    double noise_amplitude_;

    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    
  public:
    NMDirectionSelector() = delete;
    NMDirectionSelector(const NMDirectionSelector& copy) = delete;

    NMDirectionSelector(GCPtr gc, NeuritePtr neurite);
    NMDirectionSelector(const NMDirectionSelector& copy, GCPtr gc, NeuritePtr neurite);

    virtual void
    select_direction(
        const std::vector<double> &directions_weights, const Filopodia &filo,
        mtPtr rnd_engine, double total_proba, bool interacting,
        double old_angle, double &substep, double &step_length,
        double &new_angle, bool &stopped,
        stype &default_direction) override final;

    virtual void set_status(const statusMap &status) override final;
    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
