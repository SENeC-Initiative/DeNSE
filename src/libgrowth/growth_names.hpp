#ifndef NAMES_MODELS_H
#define NAMES_MODELS_H
#include "config.hpp"
#include <string>


namespace growth
{
namespace names
{


//! growth con model for          'Rw'   [string]
extern const std::string growth_cone_model;
//! number of neurites for created neuron [int]
extern const std::string num_neurites;
//! radius of soma for neurite starts  1    [micrometers]
extern const std::string soma_radius;
#define SOMA_RADIUS 20.

extern const std::string dendrite_diameter;
#define DENDRITE_DIAMETER 6.
extern const std::string axon_diameter;
#define AXON_DIAMETER 6.
//! axon angle set to perfermo experiments
extern const std::string axon_angle;
//! initial branching lenght        0     [micrometers]
extern const std::string initial_branch_lenght;
//! default branching probability
extern const std::string branching_proba_default;
#define BRANCHING_PROBA_DEFAULT 0.05

//! SPACE SENSING
extern const std::string filopodia_angular_resolution;
#define FILOPODIA_ANGULAR_RES 24
extern const std::string filopodia_finger_length;
#define FILOPODIA_FINGER_LENGTH 50.
extern const std::string filopodia_wall_affinity;
#define FILOPODIA_WALL_AFFINITY 0.2


//! CRITICAL MODELS
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

//! RANDOM WALK MODEL
extern const std::string random_walk_submodel;
//! @param speed_growth_cone      10     [micormeter/second]
extern const std::string speed_growth_cone;
extern const std::string speed_variance;
#define RW_SPEED_GROWTH_CONE 1.
//! @param persistenc_length      2000  [micrometer]
extern const std::string rw_persistence_length;
#define RW_PERSISTENCE_LENGTH 10.
extern const std::string rw_memory_tau;
#define RW_MEMORY_TAU 100.
extern const std::string rw_delta_corr;
#define RW_DELTA_CORR 100.
//@param sensing angle is choosen from experimental
// data and it's 8.2 degrees
extern const std::string rw_sensing_angle;
#define RW_SENSING_ANGLE 0.1433


//! VAN PELT BEST MODEL
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


//! ACTIN WAVE MODEL
//! actin wave trigger or not     False  [bool]
extern const std::string use_actin_waves;
#define USE_ACTIN_WAVES false
//! Actin Waves model parameters
extern const std::string actin_content;
#define ACTIN_CONTENT 0 /
extern const std::string actin_content_tau;
#define ACTIN_CONTENT_TAU -1
extern const std::string actin_wave_speed;
#define ACTIN_WAVE_SPEED 150
extern const std::string actin_freq;
#define AW_GENERATION_STEP -1

//! uniform branching event parameters
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
}
}
#endif
