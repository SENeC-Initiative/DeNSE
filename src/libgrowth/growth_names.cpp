/*
 * growth_names.cpp
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

#include "growth_names.hpp"

namespace growth
{

namespace names
{

const std::string actin_content("actin_content");
const std::string actin_content_tau("actin_content_tau");
const std::string actin_freq("actin_freq");
const std::string actin_wave_speed("actin_wave_speed");
const std::string active("active");
const std::string
    affinity_axon_axon_other_neuron("affinity_axon_axon_other_neuron");
const std::string
    affinity_axon_dendrite_other_neuron("affinity_axon_dendrite_other_neuron");
const std::string
    affinity_axon_dendrite_same_neuron("affinity_axon_dendrite_same_neuron");
const std::string affinity_axon_self("affinity_axon_self");
const std::string
    affinity_axon_soma_other_neuron("affinity_axon_soma_other_neuron");
const std::string
    affinity_axon_soma_same_neuron("affinity_axon_soma_same_neuron");
const std::string
    affinity_dendrite_axon_other_neuron("affinity_dendrite_axon_other_neuron");
const std::string
    affinity_dendrite_axon_same_neuron("affinity_dendrite_axon_same_neuron");
const std::string affinity_dendrite_dendrite_other_neuron(
    "affinity_dendrite_dendrite_other_neuron");
const std::string affinity_dendrite_dendrite_same_neuron(
    "affinity_dendrite_dendrite_same_neuron");
const std::string affinity_dendrite_self("affinity_dendrite_self");
const std::string
    affinity_dendrite_soma_other_neuron("affinity_dendrite_soma_other_neuron");
const std::string
    affinity_dendrite_soma_same_neuron("affinity_dendrite_soma_same_neuron");
const std::string axon_angle("axon_angle");
const std::string axon_polarization_weight("axon_polarization_weight");

const std::string B("B");
const std::string branching_proba_default("branching_proba_default");

const std::string critical_pull("critical_pull");

const std::string memory_decay_factor("memory_decay_factor");
const std::string dendrite_angles("dendrite_angles");
const std::string description("description");
const std::string diameter_eta_exp("diameter_eta_exp");
const std::string diameter_fraction_lb("diameter_fraction_lb");
const std::string diameter_ratio_avg("diameter_ratio_avg");
const std::string diameter_ratio_std("diameter_ratio_std");
const std::string duration_retraction("duration_of_retraction");

const std::string E("E");

const std::string filopodia_angular_resolution("filopodia_angular_resolution");
const std::string filopodia_finger_length("filopodia_finger_length");
const std::string filopodia_min_number("filopodia_min_number");
const std::string filopodia_wall_affinity("filopodia_wall_affinity");
const std::string flpl_branching_rate("flpl_branching_rate");

const std::string gamma_coeff("gamma_coeff");
const std::string gc_split_angle_mean("gc_split_angle_mean");
const std::string gc_split_angle_std("gc_split_angle_std");
const std::string growth_cone_model("growth_cone_model");

const std::string has_axon("has_axon");

const std::string initial_branch_lenght("initial_branch_lenght");
const std::string initial_diameter("initial_diameter");
const std::string interactions("interactions");

const std::string lateral_branching_angle_mean("lateral_branching_angle_mean");
const std::string lateral_branching_angle_std("lateral_branching_angle_std");

const std::string max_allowed_resolution("max_allowed_resolution");
const std::string max_arbor_length("max_arbor_length");
const std::string max_gc_number("max_gc_number");
const std::string max_sensing_angle("max_sensing_angle");
const std::string max_synaptic_distance("max_synaptic_distance");
const std::string min_branching_distance("min_branching_distance");

const std::string neurite_angles("neurite_angles");
const std::string neurite_names("neurite_names");
const std::string neurite_type("neurite_type");
const std::string noise_amplitude("noise_amplitude");
const std::string num_neurites("num_neurites");
const std::string num_growth_cones("num_growth_cones");

const std::string persistence_length("persistence_length");
const std::string points_per_circle("points_per_circle");
const std::string polarization_strength("polarization_strength");
const std::string print_time("print_time");
const std::string proba_down_move("proba_down_move");
const std::string proba_retraction("retraction_probability");

const std::string random_rotation_angles("random_rotation_angles");
const std::string res_branching_proba("res_branching_proba");
const std::string res_branching_threshold("res_branching_threshold");
const std::string res_correlation("res_correlation");
const std::string res_elongation_factor("res_elongation_factor");
const std::string res_elongation_threshold("res_elongation_threshold");
const std::string res_increase_slope("res_increase_slope");
const std::string res_leakage("res_leakage");
const std::string res_neurite_split_threshold("res_neurite_split_threshold");
const std::string res_neurite_available("res_neurite_available");
const std::string res_neurite_variance("res_neurite_variance");
const std::string res_neurite_generated("res_neurite_generated");
const std::string res_neurite_generated_tau("res_neurite_generated_tau");
const std::string res_neurite_delivery_tau("res_neurite_delivery_tau");
const std::string res_retraction_factor("res_retraction_factor");
const std::string res_retraction_threshold("res_retraction_threshold");
const std::string res_typical_gc_support("res_typical_gc_support");
const std::string res_use_ratio("res_use_ratio");
const std::string res_variance("res_variance");
const std::string res_weight_centrifugal("res_weight_centrifugal");
const std::string res_weight_diameter("res_weight_diameter");
const std::string resolution("resolution");
const std::string resource("resource");
const std::string rigidity_factor("rigidity_factor");

const std::string S("S");
const std::string scale_up_move("scale_up_move");
const std::string self_avoidance_factor("self_avoidance_factor");
const std::string self_avoidance_scale("self_avoidance_scale");
const std::string sensing_angle("sensing_angle");
const std::string soma_radius("soma_radius");
const std::string somatropic_factor("somatropic_factor");
const std::string somatropic_mode("somatropic_mode");
const std::string somatropic_scale("somatropic_scale");
const std::string speed_decay_factor("speed_decay_factor");
const std::string speed_growth_cone("speed_growth_cone");
const std::string speed_ratio_retraction("speed_ratio_retraction");
const std::string speed_variance("speed_variance");
const std::string srf_avoidance_decay("srf_avoidance_decay");
const std::string srf_avoidance_force("srf_avoidance_force");
const std::string srf_inertial_force("srf_inertial_force");
const std::string srf_inertial_decay("srf_inertial_decay");
const std::string srf_somatropic_force("srf_somatropic_force");
const std::string srf_somatropic_decay("srf_somatropic_decay");
const std::string substrate_affinity("substrate_affinity");

const std::string T("T");
const std::string taper_rate("taper_rate");

const std::string uniform_branching_rate("uniform_branching_rate");
const std::string uniform_split_rate("uniform_split_rate");
const std::string use_actin_waves("use_actin_waves");
const std::string use_critical_resource("use_critical_resource");
const std::string use_flpl_branching("use_flpl_branching");
const std::string use_uniform_branching("use_uniform_branching");
const std::string use_uniform_split("use_uniform_split");
const std::string use_van_pelt("use_van_pelt");

// recorders

const std::string event_type("event_type");
const std::string interval("interval");
const std::string level("level");
const std::string observable("observable");
const std::string observables("observables");
const std::string record_to("record_to");
const std::string restrict_to("restrict_to");
const std::string targets("targets");

const signed char lateral_branching(0);
const signed char gc_splitting(1);
const signed char gc_deletion(2);
} // namespace names

} // namespace growth
