#ifndef NAMES_MODELS_H
#define NAMES_MODELS_H
#include "config.hpp"
#include <string>


namespace growth
{

namespace names
{

/*
 * Neuron, neurite and growth cone properties
 */

//! growth con model for          'Rw'   [string]
extern const std::string growth_cone_model;
//! number of neurites for created neuron [int]
extern const std::string num_neurites;
//! radius of soma for neurite starts  1    [micrometers]
extern const std::string soma_radius;

extern const std::string axon_diameter;
extern const std::string dendrite_diameter;
//! axon angle set to perform experiments
extern const std::string axon_angle;
//! initial branching lenght        0     [micrometers]
extern const std::string initial_branch_lenght;
//! default branching probability
extern const std::string branching_proba_default;

#define BRANCHING_PROBA_DEFAULT 0.05
#define AXON_DIAMETER 6.
#define DENDRITE_DIAMETER 6.
#define SOMA_RADIUS 20.


/*
 * SPACE SENSING
 */

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

#define DURATION_RETRACTION 200.
#define FILOPODIA_MIN_NUM 24
#define FILOPODIA_FINGER_LENGTH 50.
#define FILOPODIA_SUBSTRATE_AFINITY 0.1
#define FILOPODIA_WALL_AFFINITY 2.
#define MAX_SENSING_ANGLE 1.5707963267948966  // 100 degrees max for 1 s resol
#define ONE_DEGREE 0.017453292519943295
#define PROBA_RETRACTION 0.001
#define PROBA_DOWN_MOVE 0.008
#define RW_DELTA_CORR 100.
#define RW_MEMORY_TAU 100.
#define RW_PERSISTENCE_LENGTH 10.
#define SCALE_UP_MOVE 20.
#define SENSING_ANGLE 0.1433
#define SPEED_RATIO_RETRACTION 0.2
#define SPEED_GROWTH_CONE 1.
#define WALL_AFNTY_DECAY_CST 19.098593171027442 // inverse of 3 deg in radians


/*
 * CRITICAL MODELS
 */

//! @param tub_topo_coefficient   0.1    [natural]
extern const std::string use_critical_resource;
#define USE_CRITICAL false
// extern const std::string CR_topo_coeff;
//#define CRITICAL_TOPO_COEFF 1.
// extern const std::string CR_geom_coeff;
//#define CRITICAL_GEOM_COEFF 1.
extern const std::string CR_retraction_th;
#define CRITICAL_RETRACTION_TH 0.02
extern const std::string CR_elongation_th;
#define CRITICAL_ELONGATION_TH 0.2
extern const std::string CR_initial_demand;
#define CRITICAL_INITIAL_DEMAND 1.
extern const std::string CR_amount;
#define CRITICAL_AMOUNT 1.
extern const std::string CR_speed_factor;
#define CRITICAL_SPEED_FACTOR 1.
extern const std::string CR_split_th;
#define CRITICAL_SPLIT_TH 8.
extern const std::string CR_use_ratio;
#define CRITICAL_USE_RATIO 0.2
extern const std::string CR_leakage;
#define CRITICAL_LEAKAGE 0.2
// extern const std::string CR_use_ratio;
//#define CRITICAL_USE_RATIO 0.2
// extern const std::string CR_use_ratio;
/*#define CRITICAL_USE_RATIO 0.2*/
extern const std::string CR_demand_correlation;
#define CR_DEMAND_CORRELATION 0.1
extern const std::string CR_demand_stddev;
#define CR_DEMAND_STDDEV 0.1
extern const std::string CR_demand_mean;
#define CR_DEMAND_MEAN 1


/*
 * RANDOM WALK MODEL
 */

extern const std::string random_walk_submodel;
//! @param speed_growth_cone      10     [micormeter/second]
extern const std::string speed_growth_cone;
extern const std::string speed_variance;
//! @param persistenc_length      2000  [micrometer]
extern const std::string rw_persistence_length;
extern const std::string rw_memory_tau;
extern const std::string rw_delta_corr;
//@param sensing angle is choosen from experimental
// data and it's 8.2 degrees


/*
 * VAN PELT BEST MODEL
 */

//! @param van_pelt model for branching probability and direction default: True
extern const std::string use_van_pelt;
#define USE_VAN_PELT true
//! Van_Pelt BEST model parameters
extern const std::string B;
#define VP_B 5.
extern const std::string E;
#define VP_E 0.05
extern const std::string S;
#define VP_S 1.
extern const std::string T;
#define VP_T 0.01


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
 * Uniform branching event parameters
 */

extern const std::string use_lateral_branching;
extern const std::string uniform_branching_rate;
#define UNIFORM_BRANCHING_RATE 0.1
extern const std::string lateral_branching_angle_mean;
#define LATERAL_BRANCHING_ANGLE_MEAN 90 * 3.14 / 180
extern const std::string lateral_branching_angle_std;
#define LATERAL_BRANCHING_ANGLE_STD 1. / 180 * 3.14
extern const std::string gc_split_angle_mean;
#define GC_SPLIT_ANGLE_MEAN 98.0 / 180 * 3.14
extern const std::string gc_split_angle_std;
#define GC_SPLIT_ANGLE_STD 10. / 180 * 3.14
extern const std::string angle_in_degrees;
#define ANGLE_IN_DEGREES true
extern const std::string diameter_variance;
#define DIAMETER_VARIANCE 0.1
extern const std::string diameter_eta_exp;
#define DIAMETER_ETA_EXP 1.5


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
}
}
#endif
