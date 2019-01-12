#ifndef NAMES_MODELS_H
#define NAMES_MODELS_H

#include <string>
#include <limits>

#include "config.hpp"


namespace growth
{

namespace names
{

/*
 * Neuron, neurite and growth cone properties
 */

//! neuron description [string]
extern const std::string description;
//! growth con model for          'Rw'   [string]
extern const std::string growth_cone_model;
//! number of neurites for created neuron [int]
extern const std::string num_neurites;
//! radius of soma for neurite starts  1    [micrometers]
extern const std::string soma_radius;
extern const std::string has_axon;
extern const std::string axon_polarization_weight;
extern const std::string neurite_angles;
extern const std::string polarization_strength;
extern const std::string random_rotation_angles;

extern const std::string max_gc_number;
extern const std::string max_arbor_length;

extern const std::string active;

extern const std::string axon_diameter;
extern const std::string dendrite_diameter;
//! axon angle set to perform experiments
extern const std::string axon_angle;
//! initial branching lenght        0     [micrometers]
extern const std::string initial_branch_lenght;
//! default branching probability
extern const std::string branching_proba_default;

extern const std::string taper_rate;
extern const std::string diameter_ratio_avg;
extern const std::string diameter_ratio_std;
extern const std::string diameter_eta_exp;


#define BRANCHING_PROBA_DEFAULT 0.05
#define AXON_DIAMETER 6.
#define DENDRITE_DIAMETER 6.
#define SOMA_RADIUS 8.
#define THINNING_RATIO 0.005  // lose 1 micrometer every 200 micrometers
#define MIN_DIAMETER 0.1      // diameter when a gc stops growing [micrometers]
#define DIAMETER_RATIO_AVG 1.
#define DIAM_RATIO_STD 0.02   // deviation from identical diameters on split
#define DIAMETER_ETA_EXP 2.67
#define POLA_STRENGTH 5.
#define AXON_POLA_WEIGHT 2.
#define MAX_GC_NUM 500
#define MAX_ARBOR_LENGTH 5000.


/*
 * SPACE AND SENSING
 */

extern const std::string affinity_axon_self;
extern const std::string affinity_axon_dendrite_same_neuron;
extern const std::string affinity_axon_soma_same_neuron;
extern const std::string affinity_axon_axon_other_neuron;
extern const std::string affinity_axon_dendrite_other_neuron;
extern const std::string affinity_axon_soma_other_neuron;
extern const std::string affinity_dendrite_self;
extern const std::string affinity_dendrite_axon_same_neuron;
extern const std::string affinity_dendrite_dendrite_same_neuron;
extern const std::string affinity_dendrite_soma_same_neuron;
extern const std::string affinity_dendrite_axon_other_neuron;
extern const std::string affinity_dendrite_dendrite_other_neuron;
extern const std::string affinity_dendrite_soma_other_neuron;

extern const std::string duration_retraction;
extern const std::string filopodia_min_number;
extern const std::string filopodia_finger_length;
extern const std::string filopodia_wall_affinity;
extern const std::string max_sensing_angle;
extern const std::string proba_down_move;
extern const std::string proba_retraction;
extern const std::string scale_up_move;
extern const std::string sensing_angle;
extern const std::string speed_ratio_retraction;
extern const std::string substrate_affinity;
extern const std::string points_per_circle;

#define DURATION_RETRACTION 4.
#define FILOPODIA_MIN_NUM 24
#define FILOPODIA_FINGER_LENGTH 12.
#define FILOPODIA_SUBSTRATE_AFINITY 0.1
#define FILOPODIA_WALL_AFFINITY 2.
#define MAX_SENSING_ANGLE 1.5707963267948966 // 180 degrees max for 1 s resol
#define ONE_DEGREE 0.017453292519943295
#define PROBA_RETRACTION 0.001
#define PROBA_DOWN_MOVE 0.008
#define PERSISTENCE_LENGTH 500.
#define SCALE_UP_MOVE 20.
#define SENSING_ANGLE 1.2217 // approximately 70 deg
#define SPEED_RATIO_RETRACTION 0.2
#define SPEED_GROWTH_CONE 0.05  // um/min
#define WALL_AFNTY_DECAY_CST 19.098593171027442 // inverse of 3 deg in radians
#define DEFAULT_POINTS_PER_CIRCLE 12 // number of points used to create a circle
#define MIN_FILOPODIA_FINGER_LENGTH 5. // WARNING: THIS VALUE *MUST* BE EQUAL TO
                                       // MAX_MAX_SYN_DIST


/*
 * Steering models
 */

// memory-based steering
extern const std::string rigidity_factor;
extern const std::string decay_factor;


/*
 * Direction selection models
 */
extern const std::string noise_amplitude;
extern const std::string critical_pull;


/*
 * Steering models
 */

// memory-based steering
extern const std::string rigidity_factor;
extern const std::string decay_factor;


/*
 * Direction selection models
 */
extern const std::string noise_amplitude;
extern const std::string critical_pull;


/*
 * CRITICAL MODELS
 */

extern const std::string resource;

//! @param tub_topo_coefficient   0.1    [natural]
extern const std::string use_critical_resource;
#define USE_CRITICAL false
extern const std::string CR_use_ratio;
#define CRITICAL_USE_RATIO 1.
extern const std::string CR_leakage;
#define CRITICAL_LEAKAGE 6.
extern const std::string CR_correlation;
#define CRITICAL_CORRELATION 0.
extern const std::string CR_variance;
#define CRITICAL_VARIANCE 0.1 //
extern const std::string CR_weight_diameter;
#define CRITICAL_WEIGHT_DIAMETER 1.
extern const std::string CR_weight_centrifugal;
#define CRITICAL_WEIGHT_CENTRIFUGAL 0.

extern const std::string CR_elongation_factor;
#define CRITICAL_ELONGATION_FACTOR 0.5
extern const std::string CR_elongation_th;
#define CRITICAL_ELONGATION_TH 0.35
extern const std::string CR_retraction_factor;
#define CRITICAL_RETRACTION_FACTOR 0.1
extern const std::string CR_retraction_th;
#define CRITICAL_RETRACTION_TH 0.15
extern const std::string CR_branching_th;
#define CRITICAL_BRANCHING_TH std::numeric_limits<double>::infinity()
extern const std::string CR_branching_proba;
#define CRITICAL_BRANCHING_PROBA 0.1

extern const std::string CR_neurite_available;
#define CRITICAL_AVAILABLE .
extern const std::string CR_neurite_variance;
#define CRITICAL_GEN_VAR 5.
extern const std::string CR_neurite_generated;
#define CRITICAL_GENERATED 150.
extern const std::string CR_neurite_generated_tau;
#define CRITICAL_GEN_TAU 100.
extern const std::string CR_neurite_delivery_tau;
#define CRITICAL_DEL_TAU 50.
#define CRITICAL_GEN_CORR 0.
extern const std::string CR_increase_slope;
#define CRITICAL_SLOPE 1.
extern const std::string CR_typical_gc_support;
#define CRITICAL_GC_SUPPORT 8.


/*
 * RANDOM WALK MODEL
 */

extern const std::string random_walk_submodel;
//! @param speed_growth_cone      10     [micormeter/second]
extern const std::string speed_growth_cone;
extern const std::string speed_variance;
//! @param persistence_length      500  [micrometer]
extern const std::string persistence_length;
//@param sensing angle is choosen from experimental
// data and it's 8.2 degrees

extern const std::string rw_memory_tau;
extern const std::string rw_delta_corr;

#define RW_DELTA_CORR 100.
#define RW_MEMORY_TAU 100.


// SELF REFERENTIAL MODEL
//
#define SRF_AVOIDANCE_FORCE 1
extern const std::string srf_avoidance_force;
#define SRF_AVOIDANCE_DECAY 2
extern const std::string srf_avoidance_decay;
#define SRF_INERTIAL_FORCE 1
extern const std::string srf_inertial_force;
#define SRF_INERTIAL_DECAY 2
extern const std::string srf_inertial_decay;
#define SRF_SOMATROPIC_FORCE 1
extern const std::string srf_somatropic_force;
#define SRF_SOMATROPIC_DECAY 2
extern const std::string srf_somatropic_decay;


/*
 * GROWTH CONE SPLITTING PARAMETERS
 */

extern const std::string gc_split_angle_mean;
extern const std::string gc_split_angle_std;
//! @param van_pelt model for branching probability and direction default: True
extern const std::string use_van_pelt;
//! Van_Pelt BEST model parameters
extern const std::string B;
extern const std::string E;
extern const std::string S;
extern const std::string T;

#define GC_SPLIT_ANGLE_MEAN 98.0 / 180 * 3.14
#define GC_SPLIT_ANGLE_STD 10. / 180 * 3.14
#define USE_VAN_PELT true
#define VP_B 5.
#define VP_E 0.05
#define VP_S 1.
#define VP_T 0.01


/*
 * LATERAL BRANCHING PARAMETERS
 */

extern const std::string use_flpl_branching;
extern const std::string flpl_branching_rate;
extern const std::string use_uniform_branching;
extern const std::string uniform_branching_rate;
extern const std::string lateral_branching_angle_mean;
extern const std::string lateral_branching_angle_std;

#define ANGLE_IN_DEGREES true
#define LATERAL_BRANCHING_ANGLE_MEAN 90 * 3.14 / 180
#define LATERAL_BRANCHING_ANGLE_STD 1. / 180 * 3.14
#define UNIFORM_BRANCHING_RATE 0.001


/*
 * ACTIN WAVE MODEL
 */

//! actin wave trigger or not     False  [bool]
extern const std::string use_actin_waves;
#define USE_ACTIN_WAVES false
//! Actin Waves model parameters
extern const std::string actin_content;
#define ACTIN_CONTENT 0.
extern const std::string actin_content_tau;
#define ACTIN_CONTENT_TAU -1.
extern const std::string actin_wave_speed;
#define ACTIN_WAVE_SPEED 150.
extern const std::string actin_freq;
#define AW_GENERATION_STEP -1.


/*
 * RECORDERS
 */

extern const std::string event_type;
extern const std::string interval;
extern const std::string level;
extern const std::string observable;
extern const std::string observables;
extern const std::string record_to;
extern const std::string restrict_to;
extern const std::string targets;

extern const signed char lateral_branching;
extern const signed char gc_splitting;
extern const signed char gc_deletion;

extern const std::string num_growth_cones;


/*
 * Kernel and space
 */

extern const std::string interactions;
extern const std::string max_allowed_resolution;
extern const std::string max_synaptic_distance;
extern const std::string resolution;

#define DEFAULT_MAX_RESOL 30.
#define MAX_MAX_SYN_DIST 5.


} // namespace names

} // namespace growth

#endif /* NAMES_MODELS_H */
