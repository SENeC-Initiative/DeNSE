#include "Branching.hpp"

// C++ includes
#include <functional>
#include <cmath>

// kernel includes
#include "kernel_manager.hpp"

// elements include
#include "GrowthCone.hpp"
#include "Neurite.hpp"
#include "Node.hpp"


namespace growth
{

const Event invalid_ev(
    std::make_tuple(0, std::nan(""), 0, std::string(""), -1));


Branching::Branching()
    : Branching(nullptr)
{}


Branching::Branching(NeuritePtr neurite)
    : neurite_(neurite)
    // parameters for van Pelt model
    , use_van_pelt_(false)
    , next_vanpelt_event_(invalid_ev)
    , B_(VP_B)
    , E_(VP_S)
    , S_(VP_E)
    , T_(VP_T)
    // parameters for critical_resource-driven growth
    , use_critical_resource_(false)
    , CR_amount_(CRITICAL_AMOUNT)
    , CR_tot_demand_(0)
    , CR_split_th_(CRITICAL_SPLIT_TH)
    // parameters for uniform branching
    , use_lateral_branching_(false)
    , uniform_branching_rate_(UNIFORM_BRANCHING_RATE)
    , next_lateral_event_(invalid_ev)
{
}


/**
 * @brief Initialize the scheduler mechanism for event generation
 * @details The event generation needs to be initialized for each branching
 * model. The initialization take care of reconnecting the simulation with the
 * previous one
 * even if the resolution of the simulation has changed
 *
 * @param rnd_engine
 *
 * @param new_resolution_ratio ratio between previous resolution and new
 * resolution
 *
 * @param last_simulation_n_steps steps of simulation in the previous session
 */
void Branching::initialize_next_event(mtPtr rnd_engine,
                                      double new_resolution_ratio,
                                      size_t last_simulation_n_steps)
{
    double resol = kernel().simulation_manager.get_resolution();
    if (use_lateral_branching_)
    {
        if (std::isnan(std::get<1>(next_lateral_event_)))
        {
            compute_lateral_event(rnd_engine);
        }
        else
        {
            size_t old_steps   = std::get<0>(next_lateral_event_);
            double old_substep = std::get<1>(next_lateral_event_);
            size_t new_steps   =
                old_steps*new_resolution_ratio + old_substep/resol;
            double new_substep = old_steps*resol
                               + old_substep - new_steps*resol;

            next_lateral_event_ = std::make_tuple(
                new_steps, new_substep,
                std::get<2>(next_lateral_event_),
                std::get<3>(next_lateral_event_),
                std::get<4>(next_lateral_event_)
            );

            // send it to the simulation and recorder managers
            kernel().simulation_manager.new_branching_event(
                next_lateral_event_);
        }
    }
    if (use_van_pelt_)
    {
        if (std::isnan(std::get<1>(next_vanpelt_event_)))
        {
            compute_vanpelt_event(rnd_engine);
        }
        else
        {
            size_t old_steps   = std::get<0>(next_vanpelt_event_);
            double old_substep = std::get<1>(next_vanpelt_event_);
            size_t new_steps   =
                old_steps*new_resolution_ratio + old_substep/resol;
            double new_substep = old_steps*resol
                               + old_substep - new_steps*resol;

            next_vanpelt_event_ = std::make_tuple(
                new_steps, new_substep,
                std::get<2>(next_vanpelt_event_),
                std::get<3>(next_vanpelt_event_),
                std::get<4>(next_vanpelt_event_)
            );

            // send it to the simulation and recorder managers
            kernel().simulation_manager.new_branching_event(
                next_vanpelt_event_);
        }
    }
}


/**
 * @brief Update each growth cone as the user defined in the branching model
 * @details Each model of neurite, like critical_resource has it's own
 * parameters to update, this function will be overriden by neurite's models
 *
 * @param rnd_engine
 */
void Branching::update_growth_cones(mtPtr rnd_engine)
{
    // if using critical_resource model it's necessary to recompute the amount
    // of critical_resource
    // required from each growth cone.
    //
    if (use_critical_resource_)
    {
        CR_tot_demand_ = 0;
        for (auto &gc : neurite_->growth_cones_)
        {
            //@TODO move this into the neurite
            CR_tot_demand_ += gc.second->compute_CR_demand(rnd_engine);
        }
        if (CR_tot_demand_ != 0)
        {
            CR_tot_demand_ = 1. / CR_tot_demand_;
        }
        CR_amount_ = CR_amount_ * (1 + log(1 + neurite_->growth_cones_.size()));
    }
}


/**
 * @brief Verify and actuate if a branching event is scheduled for the present
 * step
 *
 * @param step present step of the simulator
 * @param rnd_engine
 *
 * @return whether the branching event was successful or not.
 */
bool Branching::branching_event(mtPtr rnd_engine, const Event& ev)
{
    // uniform_branching_event
    bool uniform_occurence =
        (std::get<0>(next_lateral_event_) ==
         std::get<0>(ev)) &&
         (std::abs(std::get<1>(next_lateral_event_)
                   - std::get<1>(ev)) < 1e-8);
    bool van_pelt_occurence =
        (std::get<0>(next_vanpelt_event_) ==
         std::get<0>(ev)) &&
        (std::abs(std::get<1>(next_vanpelt_event_)
                  - std::get<1>(ev)) < 1e-8);
    
    if (use_lateral_branching_ and uniform_occurence)
    {
        return uniform_new_branch(rnd_engine);
    }
    // verify vanpelt event
    if (use_van_pelt_ and van_pelt_occurence)
    {
        return vanpelt_new_branch(rnd_engine);
    }
    // verify CR for each growth cone
    if (use_critical_resource_)
    {
        for (auto &gc : neurite_->growth_cones_)
        {
            if (gc.second->get_CR_received() > CR_split_th_)
            {
                gc.second->reset_CR_demand();
                CR_new_branch(rnd_engine, gc.second);
                return true;
            }

            // no new growth cone was created
            return false;
        }
    }
}


//#######################################################
//              Uniform Branching Model
//#######################################################

/**
 * @brief Compute the next lateral branching event
 * @details
 * The segment where next branching event will happen is computed
 * with an uniform distribution respect to the length of segments considered.
 * The time (in step unit) of next lateral branching event is computed with an
 * exponential distribution
 * whose rate is set with Branching::set_status(...)
 * Values are stored in the Branching instance itself.
 *
 * @param rnd_engine
 */
void Branching::compute_lateral_event(mtPtr rnd_engine)
{
    // here we compute next event with exponential distribution
    // we add one to save the case in which the distance between
    // two branching event is so short to appear zero.
    double duration = exponential_uniform_(*(rnd_engine).get());
    double current_substep = kernel().simulation_manager.get_current_substep();
    size_t current_step = kernel().simulation_manager.get_current_step();

    size_t ev_step;
    double ev_substep;
    double resol = kernel().simulation_manager.get_resolution();

    if (duration < resol - current_substep)
    {
        ev_step = current_step;
        ev_substep = current_substep + duration;
    }
    else
    {
        // remove what's left until next step
        duration -= current_substep;
        ev_step = duration / resol;
        ev_substep = duration - ev_step*resol;
        ev_step += current_step + 1;
    }

    auto n = neurite_->get_parent_neuron().lock();
    size_t neuron_gid = n->get_gid();
    std::string neurite_name = neurite_->get_name();
    next_lateral_event_ = std::make_tuple(
        ev_step, ev_substep, neuron_gid, neurite_name,
        names::lateral_branching);

    // send it to the simulation and recorder managers
    kernel().simulation_manager.new_branching_event(next_lateral_event_);

#ifndef NDEBUG
    printf("after lateral event, next lateral event in %lu, I m in %lu \n",
           std::get<0>(next_lateral_event_), kernel().simulation_manager.get_current_step());
#endif
}


/**
 * @brief Branch the neurite in a uniformly randoim chosen node
 *
 * This function is enable with the flag 'use_lateral_branching'
 * This function implement the branchin throught
 * @function Neurite::lateral_branching
 *
 * @param rnd_engine
 */
bool Branching::uniform_new_branch(mtPtr rnd_engine)
{
#ifndef NDEBUG
    printf("@@@@@@@ Lateral branching @@@@@@@@\n");
#endif
    TNodePtr branching_node;
    GCPtr branching_cone;
    double new_length = 0.;  // cone created exactly on the branch to make sure
                             // it is not created outside the environment.

    // select a random node of the tree, excluding the firstNode_
    // This is a reservoir sampling algorithm
    double max = 0;
    for (auto& cone : neurite_->growth_cones_)
    {
        // check if the elected cone is not dead and waiting for removal!
        if (not cone.second->is_dead() and cone.second->get_branch_size() > 3)
        {
            double key = powf(uniform_(*(rnd_engine).get()),
                              1. / cone.second->get_branch_size());
            if (key > max)
            {
                max            = key;
                branching_cone = cone.second;
            }
        }
    }

    branching_node = branching_cone;

    for (auto& node : neurite_->nodes_)
    {
        if (node.second->get_branch()->size() > 3)
        {
            double key = powf(uniform_(*(rnd_engine).get()),
                              1. / node.second->get_branch()->size());
            if (key > max)
            {
                max            = key;
                branching_node = node.second;
            }
        }
    }

    // if no node was suited for lateral branching skip the branching.
    if (max > 0)
    {
        // choose the point uniformly on the branch, except for firts 2 and last
        // 2 points.
        size_t branch_point = uniform_(*(rnd_engine).get()) *
                                  (branching_node->get_branch()->size() - 4) +
                              2;

        // actuate lateral branching on the elected node through the NEURITE.
        neurite_->lateral_branching(branching_node, branch_point, new_length,
                                    rnd_engine);
        next_lateral_event_ = invalid_ev;
    }

    compute_lateral_event(rnd_engine);

    if (max == 0)
    {
        return false;
    }
    return true;
}

//###################################################
//                 Van Pelt Model
//###################################################

/**
 * @brief Compute the next vanpelt branching event
 * @details
 * The time of branching is computed with the statistical branching algorithm
 * from Van Pelt
 * (see article for details) the relevant parameters are B, E, T
 * The growth cone who is going to branch is computed with a reservoir
 * sample algorithm; the weight are defined in the Van Pelt model and
 * the model parameter S is set in the Branching::set_status
 *
 * This distribution is verified to reproduce the Van Pelt results.
 *
 * @param rnd_engine
 */
void Branching::compute_vanpelt_event(mtPtr rnd_engine)
{
    // get the current second to compute the time-dependent
    // exponential decreasing probability of having a branch.
    double t_0 = kernel().get_current_seconds();
    if (t_0 > 1. / T_ * log(B_ * powf(neurite_->growth_cones_.size(), 1 - E_)))
    {
        // when previous condition is verified the expected time goes to
        // infinity, and the the branching is deactivated.
        // since it's the exponential leading this phenomena, it's not going to
        // branch again.
        use_van_pelt_ = false;
#ifndef NDEBUG
        printf("deactivating van pelt for divergence in probability\n");
#endif
    }
    else
    {
        // effectively compute the time of next branch event.
        double delta =
            -t_0 - log(exp(-T_ * t_0) -
                       powf(neurite_->growth_cones_.size(), E_ - 1) / B_) /
                       T_;
        exponential_ = std::exponential_distribution<double>(1. / delta);

        double duration = exponential_(*(rnd_engine).get());
        double current_substep =
            kernel().simulation_manager.get_current_substep();
        size_t current_step = kernel().simulation_manager.get_current_step();

        size_t ev_step;
        double ev_substep;
        double resol = kernel().simulation_manager.get_resolution();

        if (duration < resol - current_substep)
        {
            ev_step = current_step;
            ev_substep = current_substep + duration;
        }
        else
        {
            duration -= current_substep;  // remove what's left until next step
            ev_step = duration / resol;
            ev_substep = duration - ev_step*resol;
            ev_step += current_step + 1;
        }

        auto n = neurite_->get_parent_neuron().lock();
        size_t neuron_gid = n->get_gid();
        std::string neurite_name = neurite_->get_name();
        next_vanpelt_event_ = std::make_tuple(
            ev_step, ev_substep, neuron_gid, neurite_name, names::gc_splitting);

        // send it to the simulation and recorder managers
        kernel().simulation_manager.new_branching_event(next_vanpelt_event_);

#ifndef NDEBUG
        printf("after vanpelt event, next vanpelt event in %lu:%f \n",
               ev_step, ev_substep);
#endif
    }
}


/**
 * @brief Compute the Van Pelt normalization factor and update GrowthCones
 * This function will be run every time a branching happen if the
 * 'use_van_pelt' flag is set True.
 * This function implement the branchin throught
 * @function growth_cone_split(...)
 *
 */
bool Branching::vanpelt_new_branch(mtPtr rnd_engine)
{
    // simple implementation of a reservoir sampling algorithm for weigthed
    // choice
    GCPtr next_vanpelt_cone_;
    double max = 0;
    for (auto& cone : neurite_->growth_cones_)
    {
        // printf("this cone %s \n", cone->get_treeID().c_str());

        auto key = powf(uniform_(*(rnd_engine).get()),
                                 cone.second->get_centrifugal_order());
        // printf("this key %f\n", key);
        if (key > max)
        {
            max                = key;
            next_vanpelt_cone_ = cone.second;
        }
    }
#ifndef NDEBUG
    printf(" order%i, size %lu of selected\n",
           next_vanpelt_cone_->get_centrifugal_order(),
           next_vanpelt_cone_->get_branch()->size());
#endif
    //@TODO define new legth for van pelt
    double new_length = next_vanpelt_cone_->get_module();
    double new_angle, old_angle;
    double old_diameter = next_vanpelt_cone_->get_diameter();
    double new_diameter = old_diameter;
    neurite_->gc_split_angles_diameter(rnd_engine, new_angle, old_angle,
                                       new_diameter, old_diameter);
    bool success = neurite_->growth_cone_split(
        next_vanpelt_cone_, new_length, new_angle, old_angle, new_diameter,
        old_diameter);

    next_vanpelt_event_ = invalid_ev;
    compute_vanpelt_event(rnd_engine);
    return success;
}


//###################################################
//                  Tubuline model
//###################################################
/**
 * @brief Growth Cone split when critical_resource value exceed the treshold
 * @details
 * This function is enable with the flag 'use_critical_resource'
 * The GC is selected in Branching::branching_event with the chosen condition
 * This function implement the branchin throught
 * @function growth_cone_split(...)
 * The splitting angles are computed with a experimentally derived distribution
 * whose
 * parameters are set in Branching::set_status(...)
 *
 * @param rnd_engine
 * @param splitting_cone the cone who is going to split
 */
void Branching::CR_new_branch(mtPtr rnd_engine, GCPtr splitting_cone)
{
    double new_length = splitting_cone->get_CR_speed_factor() *
                        splitting_cone->get_CR_received();
    double new_angle, old_angle;
    double old_diameter = splitting_cone->get_diameter();
    double new_diameter = old_diameter;
    neurite_->gc_split_angles_diameter(rnd_engine, new_angle, old_angle,
                                       new_diameter, old_diameter);
    neurite_->growth_cone_split(splitting_cone, new_length, new_angle,
                                old_angle, new_diameter, old_diameter);
    printf("splitting cone angle is %f\n ",
           splitting_cone->move_.angle * 180 / 3.14);
}


double Branching::get_CR_quotient() const { return CR_tot_demand_; }


double Branching::get_CR_amount() const { return CR_amount_; }


void Branching::set_status(const statusMap &status)
{

    //                 Tubuline Params
    //###################################################
    get_param(status, names::use_critical_resource, use_critical_resource_);
    if (use_critical_resource_)
    {
        get_param(status, names::CR_amount, CR_amount_);
        get_param(status, names::CR_split_th, CR_split_th_);
#ifndef NDEBUG
        printf("\n"
               " CRITICAL RESOURCE BRANCHING \n"
               "%s : %f \n"
               "%s : %f \n"
               "%s : %f \n"
               "%s : %f \n",
               names::CR_amount.c_str(), CR_amount_, names::CR_split_th.c_str(),
               CR_split_th_, names::gc_split_angle_mean.c_str(),
               neurite_->gc_split_angle_mean_ * 180 / M_PI,
               names::gc_split_angle_std.c_str(),
               neurite_->gc_split_angle_std_ * 180 / M_PI);
#endif
    }

    //                 Van Pelt Params
    //###################################################
    get_param(status, names::use_van_pelt, use_van_pelt_);

    if (use_van_pelt_)
    {
        get_param(status, names::B, B_);
        get_param(status, names::E, E_);
        get_param(status, names::S, S_);
        get_param(status, names::T, T_);
#ifndef NDEBUG
        printf("############################# \n"
               " Van Pelt branching resume \n"
               "%s : %f \n"
               "%s : %f \n"
               "%s : %f \n"
               "%s : %f \n"
               "%s : %f \n"
               "%s : %f \n",
               names::B.c_str(), B_, names::E.c_str(), E_, names::S.c_str(), S_,
               names::T.c_str(), T_, names::gc_split_angle_mean.c_str(),
               neurite_->gc_split_angle_mean_ * 180 / M_PI,
               names::gc_split_angle_std.c_str(),
               neurite_->gc_split_angle_std_ * 180 / M_PI);
#endif
    }
    else
    {
        next_vanpelt_event_ = invalid_ev;
    }

    //                 Uniform_branching Params
    //###################################################
    get_param(status, names::use_lateral_branching, use_lateral_branching_);

    if (use_lateral_branching_)
    {
        get_param(status, names::uniform_branching_rate,
                  uniform_branching_rate_);

        exponential_uniform_ =
            std::exponential_distribution<double>(uniform_branching_rate_);
#ifndef NDEBUG
        printf("############################# \n"
               " Uniform branching parameters resume \n"
               "%s : %d \n"
               "%s : %f \n"
               "%s : %f \n"
               "%s : %f \n",

               names::use_lateral_branching.c_str(), use_lateral_branching_,
               names::uniform_branching_rate.c_str(), uniform_branching_rate_,
               names::lateral_branching_angle_mean.c_str(),
               neurite_->lateral_branching_angle_mean_ * 180 / M_PI,
               names::lateral_branching_angle_std.c_str(),
               neurite_->lateral_branching_angle_std_ * 180 / M_PI);
#endif
    }
    else
    {
        next_lateral_event_ = invalid_ev;
    }

    //                 Set Growth COnes Params
    //###################################################
}

void Branching::get_status(statusMap &status) const
{
    set_param(status, names::use_critical_resource, use_critical_resource_);
    if (use_critical_resource_)
    {
        set_param(status, names::CR_amount, CR_amount_);
        set_param(status, names::CR_split_th, CR_split_th_);
    }

    set_param(status, names::use_van_pelt, use_van_pelt_);
    if (use_van_pelt_)
    {
        set_param(status, names::B, B_);
        set_param(status, names::E, E_);
        set_param(status, names::S, S_);
        set_param(status, names::T, T_);
    }

    set_param(status, names::use_lateral_branching, use_lateral_branching_);
    if (use_lateral_branching_)
    {
        set_param(status, names::uniform_branching_rate,
                  uniform_branching_rate_);
    }
}
}
