/*
 * extension_resource_based.cpp
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

#include "extension_resource_based.hpp"

// c++ includes
#include <cmath>

// kernel include
#include "kernel_manager.hpp"

// lib includes
#include "config.hpp"


namespace growth
{
ResourceBasedExtensionModel::ResourceBasedExtensionModel(GCPtr gc, NeuritePtr neurite)
    : ExtensionModel(gc, neurite)
    , stochastic_tmp_(0)
    , received_(0)
    , sqrt_corr_(0)
    , noise_(0)
    , stored_(CRITICAL_ELONGATION_TH)
    , leakage_(CRITICAL_LEAKAGE)
    , use_ratio_(CRITICAL_USE_RATIO)
    , correlation_(CRITICAL_CORRELATION)
    , variance_(CRITICAL_VARIANCE)
    , elongation_factor_(CRITICAL_ELONGATION_FACTOR)
    , elongation_th_(CRITICAL_ELONGATION_TH)
    , retraction_factor_(CRITICAL_RETRACTION_FACTOR)
    , retraction_th_(CRITICAL_RETRACTION_TH)
    , branching_th_(CRITICAL_BRANCHING_TH)
    , branching_proba_(CRITICAL_BRANCHING_PROBA)
    , weight_diameter_(CRITICAL_WEIGHT_DIAMETER)
    , weight_centrifugal_(CRITICAL_WEIGHT_CENTRIFUGAL)
{
    observables_.push_back("resource");
    consumption_rate_ = use_ratio_ + 1. / leakage_;

    normal_          = std::normal_distribution<double>(0, 1);
    uniform_         = std::uniform_real_distribution<double>(0., 1.);
}


ResourceBasedExtensionModel::ResourceBasedExtensionModel(const ResourceBasedExtensionModel &copy, GCPtr gc, NeuritePtr neurite)
    : ExtensionModel(copy, gc, neurite)
    , stochastic_tmp_(0)
    , received_(copy.received_)
    , sqrt_corr_(copy.sqrt_corr_)
    , noise_(copy.noise_)
    , stored_(copy.stored_)
    , leakage_(copy.leakage_)
    , use_ratio_(copy.use_ratio_)
    , consumption_rate_(copy.consumption_rate_)
    , correlation_(copy.correlation_)
    , variance_(copy.variance_)
    , elongation_factor_(copy.elongation_factor_)
    , elongation_th_(copy.elongation_th_)
    , retraction_factor_(copy.retraction_factor_)
    , retraction_th_(copy.retraction_th_)
    , branching_th_(copy.branching_th_)
    , branching_proba_(copy.branching_proba_)
    , weight_diameter_(copy.weight_diameter_)
    , weight_centrifugal_(copy.weight_centrifugal_)
{
    normal_          = std::normal_distribution<double>(0, 1);
    uniform_         = std::uniform_real_distribution<double>(0., 1.);
}


void ResourceBasedExtensionModel::initialize_CR()
{
    sqrt_corr_ = sqrt(1 - correlation_ * correlation_);
    stored_    = (stored_ == 0) ? 1.5 * elongation_th_ : stored_;
}


/**
 * @brief General purpose function run by the
 * Neurite before the split
 *
 * This function is overwritten by each mode.
 */
void ResourceBasedExtensionModel::prepare_for_split()
{
    stored_ *= 0.5;
    received_ *= 0.5;
}


/**
 * @brief General purpose function run by the
 * Neurite after the split
 *
 * This function is overwritten by each model.
 */
void ResourceBasedExtensionModel::after_split() {}


void ResourceBasedExtensionModel::compute_res_received(double substep)
{
    received_ = neurite_ptr_->get_available_cr() *
                // local demand (weighted)
                get_res_demand() *
                // normalization factor for the neurite
                neurite_ptr_->get_quotient_cr();
}


/**
 * @brief Compute the demand of CR at the actual step
 *
 * @param rnd_engine
 *
 * @return res_demand
 */
double ResourceBasedExtensionModel::get_res_demand()
{
    double current_demand = stored_ * consumption_rate_;
    // weight by centrifugal order and diameter if required
    if (weight_centrifugal_ > 0)
    {
        current_demand *= powf(2, -weight_centrifugal_*gc_weakptr_.lock()->get_centrifugal_order());
    }
    if (weight_diameter_ > 0)
    {
        double diameter = gc_weakptr_.lock()->get_diameter();
        current_demand *= (1 + weight_diameter_*diameter*diameter);
    }

    return current_demand;
}


/**
 * @brief Compute the speed of next step
 *
 * NOTE: for the critical model, the speed is computed by the neurite as the
 * average speed over the substep, i.e. it has be set previously using the
 * results of compute_cr_speed and does not need to be recomputed here.
 */
double ResourceBasedExtensionModel::compute_speed(mtPtr rnd_engine, double substep)
{
    double speed = 0;

    // if it's over the elongation threshold the neurite will extend
    if (stored_ < retraction_th_)
    {
        speed =
            retraction_factor_ * (stored_ - retraction_th_) / retraction_th_;
    }
    else if (stored_ >= elongation_th_)
    {
        speed = elongation_factor_ * (stored_ - elongation_th_) /
                      (stored_ + elongation_th_);
    }

    return speed;
}


double ResourceBasedExtensionModel::compute_CR(
  mtPtr rnd_engine, double substep, double step_length, bool stuck)
{
    if (not stuck)
    {
        // compute received CR from the soma, with respect to other GC.
        compute_res_received(substep);

        // correlated gaussian (unit standard deviation)
        noise_ =
            noise_ * correlation_ + sqrt_corr_ * normal_(*(rnd_engine).get());

        // then use Euler for deterministic term and corrected Wiener term
        stored_ = stored_ +
                  substep * (received_ - stored_ * consumption_rate_) +
                  sqrt(substep) * variance_ * noise_;

        double branch_length = gc_weakptr_.lock()->get_branch()->get_length();

        // amount of molecule cannot be negative
        if (stored_ < 0.)
        {
            stored_ = 0.;
        }
        else if (stored_ > branching_th_ and branch_length > step_length)
        {
            double rnd_throw = uniform_(*(rnd_engine).get());
            double threshold = branching_proba_ * (stored_ - branching_th_)
                               / (stored_ + branching_th_);

            if (rnd_throw < substep*threshold)
            {
                // create an event so the growth cone will split at the next
                // step
                Time ev_time = kernel().simulation_manager.get_time();
                ev_time.update(1UL, 0.);

                auto neuron         = neurite_ptr_->get_parent_neuron().lock();
                size_t neuron_gid   = neuron->get_gid();
                std::string neurite = neurite_ptr_->get_name();
                int cone_id         = gc_weakptr_.lock()->get_node_id();

                Event ev = std::make_tuple(
                    ev_time, neuron_gid, neurite, cone_id, names::gc_splitting);

                kernel().simulation_manager.new_branching_event(ev);
            }
        }
    }

    return stored_;
}

void ResourceBasedExtensionModel::reset_res_demand()
{
    //.demand = critical_.initial_demand;
}

double ResourceBasedExtensionModel::get_speed() const { return elongation_factor_; }


double ResourceBasedExtensionModel::get_res_received() const { return received_; }


double ResourceBasedExtensionModel::get_res_speed_factor() const
{
    return elongation_factor_;
}


void ResourceBasedExtensionModel::printinfo() const
{
    printf("################ \n");
    printf("CR stored %f \n", stored_);
    printf("CR received %f \n", received_);
    printf("CR elongation_th:  %f \n", elongation_th_);
    printf("CR elongation_th:  %f \n", retraction_th_);
    printf("CR elongation_factor:  %f \n", elongation_factor_);
    printf("CR retraction_factor:  %f \n", retraction_factor_);
    printf("CR use_ratio:  %f \n", use_ratio_);
    //~ printf("CR available:  %f \n",
    //neurite_ptr_->get_available_cr(substep));

    printf("################ \n");
}


void ResourceBasedExtensionModel::set_status(const statusMap &status)
{
    // state parameters
    get_param(status, names::resource, stored_);

    // speed-related stuff
    get_param(status, names::res_elongation_factor, elongation_factor_);
    get_param(status, names::res_retraction_factor, retraction_factor_);
    get_param(status, names::res_elongation_threshold, elongation_th_);
    get_param(status, names::res_retraction_threshold, retraction_th_);

    get_param(status, names::res_branching_threshold, branching_th_);
    get_param(status, names::res_branching_proba, branching_proba_);

    // use and leakage
    get_param(status, names::res_use_ratio, use_ratio_);
    get_param(status, names::res_leakage, leakage_);
    get_param(status, names::res_correlation, correlation_);
    get_param(status, names::res_variance, variance_);
    get_param(status, names::res_weight_diameter, weight_diameter_);
    get_param(status, names::res_weight_centrifugal, weight_centrifugal_);

    consumption_rate_ = use_ratio_ + 1. / leakage_;

    initialize_CR();

#ifndef NDEBUG
    printinfo();
#endif
}


void ResourceBasedExtensionModel::get_status(statusMap &status) const
{
    // state parameters
    set_param(status, names::resource, stored_, "micromole / liter");

    // speed-related
    set_param(status, names::res_elongation_factor, elongation_factor_, "micrometer / minute");
    set_param(status, names::res_retraction_factor, retraction_factor_, "micrometer / minute");
    set_param(status, names::res_elongation_threshold, elongation_th_, "micromole / liter");
    set_param(status, names::res_retraction_threshold, retraction_th_, "micromole / liter");

    set_param(status, names::res_branching_threshold, branching_th_, "micromole / liter");
    set_param(status, names::res_branching_proba, branching_proba_, "");

    // use and leakage
    set_param(status, names::res_use_ratio, use_ratio_, "1 / minute");
    set_param(status, names::res_leakage, leakage_, "minute");
    set_param(status, names::res_correlation, correlation_, "");
    set_param(status, names::res_variance, variance_, "micromole / liter / minute**0.5");
    set_param(status, names::res_weight_diameter, weight_diameter_, "");
    set_param(status, names::res_weight_centrifugal, weight_centrifugal_, "");
}


/**
 * @brief Get the current value of one of the observables
 */
double ResourceBasedExtensionModel::get_state(const char *observable) const
{
    double value = 0.;

    TRIE(observable)
    CASE("resource")
    value = stored_;
    ENDTRIE;

    return value;
}

} // namespace growth
