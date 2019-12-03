/*
 * GrowthCone.hpp
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

#ifndef GROWTHCONE_H
#define GROWTHCONE_H

#define _USE_MATH_DEFINES

#include <random>

#include "Neuron.hpp"
#include "Node.hpp"

// lib include
#include "config_impl.hpp"
#include "elements_types.hpp"
#include "spatial_types.hpp"


namespace growth
{

/*
 * Growth Cone is an abstract class
 * It is a base model which is overloaded by more detailed models.
 *
 * Growthcone is a (biological) relevant class, it computes the forces acting on
 * the growth cone from environment or intrinsic phenomena; them it actualizes
 * the step, storing points into its branch.
 *
 * A GrowthCone is always at the head of a certain branch.
 * The object is never deleted: when neurites branch, a new GC is created,
 * the other GC skips the branching point and keeps on growing.
 * GrowthCone inherits from Node.
 */
class GrowthCone : public TopologicalNode,
                   public std::enable_shared_from_this<GrowthCone>
{

    friend class Neurite;
    friend class Branching;

  protected:
    stype neuron_id_;
    std::string neurite_name_;
    const std::string model_;
    bool using_environment_; // whether we're embedded in space
    bool sensing_required_;
    double resol_;
    double sqrt_resol_;
    double adaptive_timestep_;
    double timestep_divider_;
    std::string current_area_; // name of the area where the GC is
    unsigned char max_stop_;
    unsigned char current_stop_;
    bool stuck_;
    bool stopped_;
    bool interacting_;
    bool update_filopodia_;
    std::vector<std::string> observables_;

    // affinity
    Affinities aff_;

    // motion-related data
    bool active_;
    bool just_retracted_;
    char turning_;
    double turned_;
    double delta_angle_;
    double sensing_angle_;
    bool sensing_angle_set_;
    double avg_speed_;
    double local_avg_speed_;
    double speed_variance_;
    double local_speed_variance_;
    double duration_retraction_; // duration of a retraction period (seconds)
    double proba_retraction_;    // proba of retracting when stuck
    double retracting_todo_;     // duration left to retract
    double speed_ratio_retraction_;
    double proba_down_move_; // proba of going down if bottom out of reach
    double scale_up_move_;   // maximal height that GC can cross upwards
    double retraction_time_;

    space_tree_map current_neighbors_;

    double max_sensing_angle_;
    stype min_filopodia_; // minimal number of filopodia
    stype num_filopodia_; // minimal number of filopodia

    Filopodia filopodia_;
    Move move_;

    double total_proba_; // integrated probability of all possible moves
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    std::exponential_distribution<double> exponential_;

  public:
    // public status required for models and model manager
    GrowthCone(const std::string &model);
    ~GrowthCone();
    //! Copy constructor for GrowthCone
    GrowthCone(const GrowthCone &copy);

    virtual GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                        double distanceToParent, const BPoint &position,
                        double angle) = 0;

    void update_topology(BaseWeakNodePtr parent, NeuritePtr ownNeurite,
                         double distanc_to_parent,
                         const BPoint &position, double angle);

    // growth
    void grow(mtPtr rnd_engine, stype cone_n, double substep);
    void retraction(double distance, stype cone_n, int omp_id);

    virtual void prune(stype cone_n);

    // compute direction
    bool sense_surroundings(std::vector<double> &directions_weights,
                            std::vector<bool> &wall_presence, double substep,
                            mtPtr rnd_engine);

    virtual void
    compute_direction_probabilities(std::vector<double> &directions_weights,
                                    double substep) = 0;
    void make_move(const std::vector<double> &directions_weights,
                   const std::vector<std::string> &new_pos_area,
                   double &substep, mtPtr rnd_engine, int omp_id);
    virtual void select_direction(const std::vector<double> &directions_weights,
                                  mtPtr rnd_engine, double &substep,
                                  double &new_angle,
                                  stype &default_direction) = 0;

    double check_retraction(double substep, mtPtr rnd_engine);
    void change_sensing_angle(double angle);

    // extension
    void compute_module(double substep);
    virtual void compute_speed(mtPtr rnd_engine, double substep) = 0;

    void init_filopodia();

    void set_angle(double angle);
    virtual void prepare_for_split() = 0;
    virtual void after_split()       = 0;

    virtual void set_position(const BPoint &pos) override final;

    // get functions
    double get_module() const;
    virtual double get_state(const std::string &observable) const;
    virtual double get_growth_cone_speed() const;
    bool just_retracted() const;
    stype get_neuron_id() const;
    const std::string &get_neurite_name() const;
    const std::string &get_model_name() const;
    const BPolygonPtr get_last_segment() const;
    bool is_active() const;
    double get_self_affinity() const;

    // status and kernel-related functions
    virtual void set_status(const statusMap &status);
    virtual void get_status(statusMap &status) const;
    void update_kernel_variables();
    void update_growth_properties(const std::string &area_name);
    void update_filopodia();
};

} // namespace growth

#endif
