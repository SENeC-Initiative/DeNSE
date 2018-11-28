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
const std::string axon_angle("axon_angle");
const std::string axon_diameter("axon_diameter");
const std::string axon_polarization_weight("axon_polarization_weight");

const std::string B("B");
const std::string branching_proba_default("branching_proba_default");

const std::string CR_branching_proba("CR_branching_proba");
const std::string CR_branching_th("CR_branching_th");
const std::string CR_correlation("CR_correlation");
const std::string CR_elongation_factor("CR_elongation_factor");
const std::string CR_elongation_th("CR_elongation_th");
const std::string CR_increase_slope("CR_increase_slope");
const std::string CR_leakage("CR_leakage");
const std::string CR_neurite_split_th("CR_neurite_split_th");
const std::string CR_neurite_available("CR_neurite_available");
const std::string CR_neurite_variance("CR_neurite_variance");
const std::string CR_neurite_generated("CR_neurite_generated");
const std::string CR_neurite_generated_tau("CR_neurite_generated_tau");
const std::string CR_neurite_delivery_tau("CR_neurite_delivery_tau");
const std::string CR_retraction_factor("CR_retraction_factor");
const std::string CR_retraction_th("CR_retraction_th");
const std::string CR_typical_gc_support("CR_typical_gc_support");
const std::string CR_use_ratio("CR_use_ratio");
const std::string CR_variance("CR_variance");
const std::string CR_weight_centrifugal("CR_weight_centrifugal");
const std::string CR_weight_diameter("CR_weight_diameter");
const std::string critical_pull("critical_pull");

const std::string decay_factor("decay_factor");
const std::string dendrite_angles("dendrite_angles");
const std::string dendrite_diameter("dendrite_diameter");
const std::string description("description");
const std::string diameter_eta_exp("diameter_eta_exp");
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

const std::string lateral_branching_angle_mean("lateral_branching_angle_mean");
const std::string lateral_branching_angle_std("lateral_branching_angle_std");

const std::string max_arbor_length("max_arbor_length");
const std::string max_gc_number("max_gc_number");
const std::string max_sensing_angle("max_sensing_angle");

const std::string neurite_angles("neurite_angles");
const std::string noise_amplitude("noise_amplitude");
const std::string num_neurites("num_neurites");
const std::string num_growth_cones("num_growth_cones");

const std::string persistence_length("persistence_length");
const std::string polarization_strength("polarization_strength");
const std::string proba_down_move("proba_down_move");
const std::string proba_retraction("retraction_probability");

const std::string random_walk_submodel("random_walk_submodel");
const std::string resource("resource");
const std::string rigidity_factor("rigidity_factor");
const std::string rw_memory_tau("rw_memory_tau");
const std::string rw_delta_corr("rw_delta_corr");

const std::string S("S");
const std::string scale_up_move("scale_up_move");
const std::string sensing_angle("sensing_angle");
const std::string soma_radius("soma_radius");
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

const std::string random_rotation_angles("random_rotation_angles");

const std::string uniform_branching_rate("uniform_branching_rate");
const std::string use_actin_waves("use_actin_waves");
const std::string use_critical_resource("use_critical_resource");
const std::string use_flpl_branching("use_flpl_branching");
const std::string use_uniform_branching("use_uniform_branching");
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
