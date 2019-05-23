/*
 * ActinWave.hpp
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

#ifndef ACTINWAVE_H
#define ACTINWAVE_H

#define _USE_MATH_DEFINES

#include <random>

#include "elements_types.hpp"
#include "growth_names.hpp"
#include "spatial_types.hpp"


namespace growth
{

class Node;

class ActinWave : public std::enable_shared_from_this<ActinWave>
{
    friend class Neurite;

  private:
    NeuritePtr ownNeurite_;
    TNodePtr targetNode_;
    float distanceToTarget_;
    // models variable
    double actin_content_;
    double actin_content_tau_;
    double actin_wave_speed_;
    double lateral_branching_proba_;
    std::uniform_real_distribution<double> branch_distribution_;
    std::normal_distribution<double> normal_distribution;

  public:
    /**
     * @brief ActinWave is a single ActinWave element propagating in the
     * neurite.
     * ActinWave requires the 'Neuron' supports the AW models
     * Each 'ActinWave', for each step can:
     *
     * - propagate and keep walking up the neurite Branch
     * - reinforce a GrowthCone
     *
     * @param NodePtr Node is the Node leading the segment of neurite where AW
     * is.
     * @param double DistanceToTarget is the distance to such a node.
     */
    ActinWave(TNodePtr targetNode, double actin_content, NeuritePtr);
    ActinWave(const ActinWave &);
    //~ ~ActinWave();

    void step(mtPtr rnd_engine, double substep);
    void initialize_AW_distributions();
    void actin_on_node();
    void actin_on_growth_cone();
    void actin_make_branch(mtPtr);
    void get_geometry(BPoint &xy, double &angle, mtPtr rnd_engine) const;

    virtual void set_status(const statusMap &status);
    virtual void get_status(statusMap &status) const;
};

double get_angle(mtPtr, double);
} // namespace growth

#endif
