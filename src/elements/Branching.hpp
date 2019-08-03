/*
 * Branching.hpp
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

#ifndef BRANCHING_H
#define BRANCHING_H

// std c++
#include <random>
// elements includes
#include "Branch.hpp"
#include "Node.hpp"
#include "growth_names.hpp"


// libgrowth includes
#include "elements_types.hpp"
#include "spatial_types.hpp"


namespace growth
{
class Neurite;


class Branching
{
    friend class Neurite;

  private:
    std::uniform_real_distribution<double> uniform_;
    std::exponential_distribution<> exponential_;
    std::exponential_distribution<> exponential_uniform_;
    std::exponential_distribution<> exponential_flpl_;
    std::exponential_distribution<> exponential_usplit_;

    NeuritePtr neurite_;

    // variables for van Pelt branching model
    bool use_van_pelt_;
    double van_pelt_norm_;
    Event next_vanpelt_event_;
    double B_;
    double E_;
    double S_;
    double T_;

    bool use_critical_resource_;

    // variables for uniform branching
    bool use_uniform_branching_;
    Event next_uniform_event_;
    double uniform_branching_rate_;

    // variables for non-uniform lateral branching
    bool use_flpl_branching_;
    Event next_flpl_event_;
    double flpl_branching_rate_;

    // variables for non-uniform lateral branching
    bool use_uniform_split_;
    Event next_usplit_event_;
    double uniform_split_rate_;

  public:
    Branching(NeuritePtr neurite);
    Branching();
    Branching(const Branching &cpy);
    // event handlers functions
    void compute_next_event(mtPtr rnd_engine);
    void set_branching_event(Event &ev, signed char ev_type,
                             double duration);
    bool branching_event(mtPtr rnd_engine, const Event &ev);
    void update_splitting_cones(TNodePtr branching_cone,
                                GCPtr second_cone, NodePtr new_node);

    // van Pelt branching functions
    bool vanpelt_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                            stype &branching_point, mtPtr rnd_engine,
                            GCPtr &second_cone);
    void compute_vanpelt_event(mtPtr rnd_engine);

    // uniform split functions
    bool usplit_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                          size_t &branching_point, mtPtr rnd_engine,
                          GCPtr &second_cone);
    void compute_usplit_event(mtPtr rnd_engine);

    // uniform branching functions
    bool uniform_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                            stype &branching_point, mtPtr rnd_engine);
    void compute_uniform_event(mtPtr rnd_engine);

    // powerlaw branching functions
    bool flpl_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                         stype &branching_point, mtPtr rnd_engine);
    void compute_flpl_event(mtPtr rnd_engine);

    // critical_resource functions
    bool res_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                        stype &branching_point, mtPtr rnd_engine,
                        GCPtr& second_cone, const Event &ev);

    void initialize_next_event(mtPtr rnd_engine);

    // Get/set functions
    void set_status(const statusMap &);
    //~ void get_status(statusMap &) const;
    void get_status(statusMap &) const;
};
} // namespace growth
#endif /* BRANCHING_H */
