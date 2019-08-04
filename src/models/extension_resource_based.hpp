/*
 * extension_resource_based.hpp
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

#ifndef RB_EL_H
#define RB_EL_H

#include "extension_interface.hpp"


namespace growth
{

class ResourceBasedExtensionModel : public virtual ExtensionModel
{
    friend class Neurite;

  protected:
    // parameters
    double leakage_;
    double use_ratio_;
    double consumption_rate_;

    // resource
    double stochastic_tmp_;
    double received_;
    double stored_;

    // noise
    double variance_;
    double correlation_;
    double sqrt_corr_;
    double noise_;

    // weights
    double weight_diameter_;
    double weight_centrifugal_;

    // speed
    double elongation_factor_;
    double elongation_th_;
    double retraction_factor_;
    double retraction_th_;
    double branching_th_;
    double branching_proba_;

    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;

  public:
    ResourceBasedExtensionModel(GCPtr gc, NeuritePtr neurite);
    ResourceBasedExtensionModel(const ResourceBasedExtensionModel &) = delete;
    ResourceBasedExtensionModel(const ResourceBasedExtensionModel &copy, GCPtr gc, NeuritePtr neurite);

    void initialize_CR();
    void prepare_for_split() override;
    void after_split() override;
    void reset_res_demand();

    double compute_speed(mtPtr rnd_engine, double substep) override;

    void compute_res_received(double substep);
    double compute_CR(mtPtr rnd_engine, double substep, double step_length,
                      bool stuck);

    // getter functions
    void printinfo() const;
    double get_res_demand();
    double get_res_received() const;
    double get_res_speed_factor() const;
    double get_speed() const;

    // status
    void set_status(const statusMap &) override;
    void get_status(statusMap &) const override;
    virtual double get_state(const std::string& observable) const override;
};

} // namespace growth
#endif
