/*
 * direction_select_rt.hpp
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

#ifndef RT_DIRSEL_H
#define RT_DIRSEL_H

#include "direction_select_interface.hpp"


namespace growth
{

class RTDirectionSelector : public virtual DirectionSelectModel
{
  private:
    double persistence_length_; // typical distance on which gc goes straight
    double critical_pull_; // pull necessary to trigger tumble with proba 1 if
                           // applied orthogonally
    double sensing_angle_; // duplicate of GC entry to convenience
    double tau_;           // tumbling rate
    double next_tumble_;   // distance left to next tumble
    double tumbling_;      // whether a tumble is happening
    double p_tumble_on_stop_; // proba of tumbling if stopped
    stype num_tumbles_;       // keep track of tumbling occurrences

    std::exponential_distribution<double> exponential_rt_;
    std::uniform_real_distribution<double> uniform_;

  public:
    RTDirectionSelector() = delete;
    RTDirectionSelector(GCPtr gc, NeuritePtr neurite);
    RTDirectionSelector(const RTDirectionSelector &copy) = delete;
    RTDirectionSelector(const RTDirectionSelector &copy, GCPtr gc,
                        NeuritePtr neurite);

    virtual void select_direction(const std::vector<double> &directions_weights,
                                  const Filopodia &filo, mtPtr rnd_engine,
                                  double total_proba, bool interacting,
                                  double old_angle, double &substep,
                                  double &step_length, double &new_angle,
                                  bool &stopped,
                                  stype &default_direction) override final;

    void initialize_rt();

    virtual double
    get_state(const std::string &observable) const override final;
    virtual void set_status(const statusMap &status) override final;
    virtual void get_status(statusMap &status) const override final;
};

} // namespace growth

#endif
