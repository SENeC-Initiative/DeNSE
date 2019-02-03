#include "Branching.hpp"

// C++ includes
#include <cmath>
#include <functional>
#include <memory>

// kernel includes
#include "kernel_manager.hpp"

// elements include
#include "GrowthCone.hpp"
#include "Neurite.hpp"
#include "Node.hpp"

// models includes
//~ #include "gc_critical.hpp"


// Integrate the differential equation
// dA =  - (A_m-A)/tau
double f_deterministic(double A_m, double A, double tau, double dt,
                       double stochastic_tmp);


namespace growth
{

const Event invalid_ev(std::make_tuple(Time(), 0, std::string(""), -1));


Branching::Branching()
    : Branching(nullptr)
{
}

Branching::Branching(NeuritePtr neurite)
    : neurite_(neurite)
    // parameters for van Pelt model
    , use_van_pelt_(false)
    , next_vanpelt_event_(invalid_ev)
    , B_(VP_B)
    , E_(VP_S)
    , S_(VP_E)
    , T_(VP_T)
    // parameters for uniform branching
    , use_uniform_branching_(false)
    , uniform_branching_rate_(UNIFORM_BRANCHING_RATE)
    , next_uniform_event_(invalid_ev)
    // parameters for front lateral branching
    , use_flpl_branching_(false)
    , flpl_branching_rate_(UNIFORM_BRANCHING_RATE)
    , next_flpl_event_(invalid_ev)
{
    exponential_uniform_ =
        std::exponential_distribution<double>(uniform_branching_rate_);
    exponential_flpl_ =
        std::exponential_distribution<double>(flpl_branching_rate_);
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
void Branching::initialize_next_event(mtPtr rnd_engine)
{
    if (use_uniform_branching_ and next_uniform_event_ == invalid_ev)
    {
        compute_uniform_event(rnd_engine);
    }
    if (use_flpl_branching_ and next_flpl_event_ == invalid_ev)
    {
        compute_flpl_event(rnd_engine);
    }
    if (use_van_pelt_ and next_vanpelt_event_ == invalid_ev)
    {
        compute_vanpelt_event(rnd_engine);
    }
}


/**
 * Set the event content based on the duration between the current time and
 * the next event
 */
void Branching::set_branching_event(Event &ev, signed char ev_type,
                                    double duration)
{
    Time ev_time = kernel().simulation_manager.get_time();

    // separate duration into days, hours, minutes, seconds
    double total_hours = std::floor(duration / 60.);

    double seconds = (duration - std::floor(duration))*60;
    ev_time.set_sec(ev_time.get_sec() + seconds);

    unsigned char minutes = duration - 60*total_hours;
    ev_time.set_min(ev_time.get_min() + minutes);

    size_t days = total_hours / 24.;
    ev_time.set_day(ev_time.get_day() + days);

    unsigned char hours = total_hours - 24*days;
    ev_time.set_hour(ev_time.get_hour() + hours);

    // set the informations of the event
    auto n                   = neurite_->get_parent_neuron().lock();
    size_t neuron_gid        = n->get_gid();
    std::string neurite_name = neurite_->get_name();

    ev = std::make_tuple(ev_time, neuron_gid, neurite_name, ev_type);
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
bool Branching::branching_event(mtPtr rnd_engine, const Event &ev)
{
    // uniform_branching_event
    bool uniform_occurence =
        std::get<edata::TIME>(next_uniform_event_) == std::get<edata::TIME>(ev);
    bool flpl_occurence =
        std::get<edata::TIME>(next_flpl_event_) == std::get<edata::TIME>(ev);
    bool van_pelt_occurence =
        std::get<edata::TIME>(next_vanpelt_event_) == std::get<edata::TIME>(ev);

    // prepare pointers to nodes and branching_point
    TNodePtr branching_node = nullptr;
    NodePtr new_node = nullptr;
    size_t branching_point;

    bool success(false);

    // check uniform event
    if (use_uniform_branching_ and uniform_occurence)
    {
        success = uniform_new_branch(branching_node, new_node, branching_point,
                                     rnd_engine);
    }
    // check flpl event
    if (use_flpl_branching_ and flpl_occurence)
    {
        success = flpl_new_branch(branching_node, new_node, branching_point,
                                  rnd_engine);
    }
    // verify vanpelt event
    if (use_van_pelt_ and van_pelt_occurence)
    {
        success = vanpelt_new_branch(branching_node, new_node, branching_point,
                                     rnd_engine);
    }

    if (success)
    {
        // update the R-tree: all points from branch_point on need to be
        // removed and reassigned to the new point.
        int omp_id = kernel().parallelism_manager.get_thread_local_id();

        kernel().space_manager.update_objects_branching(
            branching_node, new_node, branching_point,
            neurite_->get_parent_neuron().lock()->get_gid(),
            neurite_->get_name(), omp_id);
    }

    // note that critical resource splitting is instantaneous and unplanned
    // thus it is never communicated through branching events
    // @todo change that

    return success;
}


//#######################################################
//              Uniform Branching Model
//#######################################################

/**
 * @brief Compute the next uniform branching event
 * @details
 * The segment where next branching event will happen is computed
 * with an uniform distribution with respect to the length of segments
 * considered. The time (in step unit) of next lateral branching event is
 * computed with an exponential distribution whose rate is set with
 * Branching::set_status(...) Values are stored in the Branching instance
 * itself.
 *
 * @param rnd_engine
 */
void Branching::compute_uniform_event(mtPtr rnd_engine)
{
    // here we compute next event with exponential distribution
    // we add one to save the case in which the distance between
    // two branching event is so short to appear zero.
    if (not neurite_->growth_cones_.empty())
    {
        double duration = exponential_uniform_(*(rnd_engine).get());

        set_branching_event(next_uniform_event_, names::lateral_branching, duration);

        // send it to the simulation and recorder managers
        kernel().simulation_manager.new_branching_event(next_uniform_event_);
    }
}


/**
 * @brief Compute the next non-uniform flpl branching event
 * @details
 * The segment where next branching event will happen is computed
 * with an power law distribution with respect to the distance from the branch
 * tip.
 * The time (in step unit) of next flpl branching event is computed with an
 * exponential distribution
 * whose rate is set with Branching::set_status(...)
 * Values are stored in the Branching instance itself.
 *
 * @param rnd_engine
 */
void Branching::compute_flpl_event(mtPtr rnd_engine)
{
    if (not neurite_->growth_cones_.empty())
    {
        // here we compute next event with exponential distribution
        double duration = exponential_flpl_(*(rnd_engine).get());

        set_branching_event(next_flpl_event_, names::lateral_branching, duration);

        // send it to the simulation and recorder managers
        kernel().simulation_manager.new_branching_event(next_flpl_event_);
    }
}


/**
 * @brief Branch on a random node following power law distribution from the tips
 *
 * This function is enabled through the flag 'use_flpl_branching'
 * This function implements the branchin throught
 * @function Neurite::lateral_branching
 *
 * @param rnd_engine
 */
bool Branching::flpl_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                                size_t &branching_point, mtPtr rnd_engine)
{
#ifndef NDEBUG
    printf("@@@@@@@ Lateral branching (FLPL) @@@@@@@@\n");
#endif
    branching_node = nullptr;

    if (not neurite_->growth_cones_.empty())
    {
        GCPtr branching_cone;
        bool success = false;

        // select a random node of the tree, excluding the firstNode_
        // This is a reservoir sampling algorithm
        double max = 0;
        for (auto &cone : neurite_->growth_cones_)
        {
            // check if the elected cone is not dead and waiting for removal!
            if (not cone.second->is_dead()
                and cone.second->get_branch_size() > 3)
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

        for (auto &node : neurite_->nodes_)
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
            // choose the point with a power law distribution over the branch
            // length,
            // where y is a uniform variate, n is the distribution power,
            // x0 and x1 define the range of the distribution, and x is your
            // power-law distributed variate.
            double y = uniform_(*(rnd_engine).get());
            int x_0  = 2;
            int x_1  = branching_node->get_branch()->size() - 2;
            // TODO check where this 2 come from.
            int n = 2;
            double a =
                (powf(x_1, (n + 1)) - powf(x_0, (n + 1))) * y + powf(x_0, (n + 1));
            int x = (int)powf(a, 1. / (n + 1));

            branching_point = x;

            // lateral branching on the elected node through the NEURITE.
            success = neurite_->lateral_branching(
                branching_node, branching_point, new_node, rnd_engine);

            next_flpl_event_ = invalid_ev;
        }

        compute_flpl_event(rnd_engine);

        return success;
    }

    next_flpl_event_ = invalid_ev;

    return false;
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
bool Branching::uniform_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                                   size_t &branching_point, mtPtr rnd_engine)
{
#ifndef NDEBUG
    printf("@@@@@@@ Lateral branching @@@@@@@@\n");
#endif
    branching_node = nullptr;
    new_node       = nullptr;

    if (not neurite_->growth_cones_.empty())
    {
        GCPtr branching_cone;
        bool success = false;

        // select a random node of the tree, excluding the firstNode_
        // This is a reservoir sampling algorithm
        double max = 0;
        for (auto &cone : neurite_->growth_cones_)
        {
            // check if the elected cone is not dead and waiting for removal!
            if (not cone.second->is_dead()
                and cone.second->get_branch_size() > 3)
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

        for (auto &node : neurite_->nodes_)
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
            // choose the point uniformly on the branch, except for firts 2 and
            // last 2 points.
            branching_point = uniform_(*(rnd_engine).get()) *
                                  (branching_node->get_branch()->size() - 4) +
                              2;

            // actuate lateral branching on the elected node through the
            // NEURITE.
            success = neurite_->lateral_branching(
                branching_node, branching_point, new_node, rnd_engine);
            next_uniform_event_ = invalid_ev;
        }

        compute_uniform_event(rnd_engine);

        return success;
    }

    next_uniform_event_ = invalid_ev;

    return false;
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
    if (not neurite_->growth_cones_.empty())
    {
        // get the current second to compute the time-dependent
        // exponential decreasing probability of having a branch.
        double t_0     = kernel().simulation_manager.get_current_minutes();
        size_t num_gcs = neurite_->growth_cones_.size() +
                         neurite_->growth_cones_inactive_.size();

        double delta = exp((t_0 + 1) / T_) * T_ / B_ * powf(num_gcs, E_);

        //-t_0 - log(exp(-T_ * t_0) -
        // powf(neurite_->growth_cones_.size(), E_ - 1) / B_) /T_;
        exponential_ = std::exponential_distribution<double>(1. / delta);

        double duration = exponential_(*(rnd_engine).get());

        if (duration > std::numeric_limits<size_t>::max())
        {
            next_vanpelt_event_ = invalid_ev;
        }
        else
        {
            set_branching_event(next_vanpelt_event_, names::gc_splitting, duration);

            // send it to the simulation and recorder managers
            kernel().simulation_manager.new_branching_event(next_vanpelt_event_);
        }
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
bool Branching::vanpelt_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                                   size_t &branching_point, mtPtr rnd_engine)
{
    branching_node = nullptr;
    new_node       = nullptr;

    if (not neurite_->growth_cones_.empty())
    {
        // simple implementation of a reservoir sampling algorithm for weigthed
        // choice
        GCPtr next_vanpelt_cone_;
        double weight, total_weight(0);
        std::unordered_map<size_t, double> weights;

        for (auto cone : neurite_->growth_cones_)
        {
            weight = powf(2, -cone.second->get_centrifugal_order() * S_);
            weights[cone.first] = weight;
            total_weight += weight;
        }
        double extracted = uniform_(*(rnd_engine).get()) * total_weight;

        total_weight = 0;
        for (auto cone : neurite_->growth_cones_)
        {
            total_weight      += weights[cone.first];
            next_vanpelt_cone_ = cone.second;
            if (total_weight >= extracted)
            {
                break;
            }
        }

        // set branching node and point
        branching_node  = next_vanpelt_cone_;
        branching_point = branching_node->get_branch()->size() - 1;

        //@TODO define new legth for van pelt
        double new_length = 0.5*next_vanpelt_cone_->get_diameter();
        double new_angle, old_angle;
        double old_diameter = next_vanpelt_cone_->get_diameter();
        double new_diameter = old_diameter;
        neurite_->gc_split_angles_diameter(rnd_engine, new_angle, old_angle,
                                           new_diameter, old_diameter);
        bool success =
            neurite_->growth_cone_split(next_vanpelt_cone_, new_length,
                                        new_angle, old_angle, new_diameter,
                                        old_diameter, new_node);

        next_vanpelt_event_ = invalid_ev;
        compute_vanpelt_event(rnd_engine);
#ifndef NDEBUG
        printf("VP branching on %lu, %s, %lu at %lu\n", neurite_->get_parent_neuron().lock()->get_gid(), neurite_->get_name().c_str(), branching_node->get_nodeID(), branching_point);
#endif
        return success;
    }

    next_vanpelt_event_ = invalid_ev;

    return false;
}


//###################################################
//                  Tubuline model
//###################################################
/**
 * @brief Growth Cone split when critical_resource value exceed the treshold
 * @details
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
void Branching::res_new_branch(TNodePtr &branching_node, NodePtr &new_node,
                              size_t &branching_point, mtPtr rnd_engine)
{
    //~ std::shared_ptr<GrowthCone_Critical> gcc =
        //~ std::dynamic_pointer_cast<GrowthCone_Critical>(splitting_cone);

    //~ double new_length = gcc->get_res_speed_factor() * gcc->get_res_received();
    //~ double new_angle, old_angle;
    //~ double old_diameter = splitting_cone->get_diameter();
    //~ double new_diameter = old_diameter;
    //~ neurite_->gc_split_angles_diameter(rnd_engine, new_angle, old_angle,
                                       //~ new_diameter, old_diameter);
    //~ neurite_->growth_cone_split(splitting_cone, new_length, new_angle,
                                //~ old_angle, new_diameter, old_diameter);
}


void Branching::set_status(const statusMap &status)
{

    //                 Van Pelt Params
    //###################################################
    get_param(status, names::use_van_pelt, use_van_pelt_);

    get_param(status, names::B, B_);
    get_param(status, names::E, E_);
    get_param(status, names::S, S_);
    get_param(status, names::T, T_);

    if (not use_van_pelt_)
    {
        next_vanpelt_event_ = invalid_ev;
    }

    //                 Uniform_branching Params
    //###################################################
    get_param(status, names::use_uniform_branching, use_uniform_branching_);

    get_param(status, names::uniform_branching_rate,
              uniform_branching_rate_);

    exponential_uniform_ =
        std::exponential_distribution<double>(uniform_branching_rate_);

    if (not use_uniform_branching_)
    {
        next_uniform_event_ = invalid_ev;
    }

    //                 FLPL branching Params
    //###################################################

    get_param(status, names::use_flpl_branching, use_flpl_branching_);

    get_param(status, names::flpl_branching_rate, flpl_branching_rate_);

    exponential_flpl_ =
        std::exponential_distribution<double>(flpl_branching_rate_);

    if (not use_flpl_branching_)
    {
        next_flpl_event_ = invalid_ev;
    }
}


void Branching::get_status(statusMap &status) const
{
    set_param(status, names::use_van_pelt, use_van_pelt_, "");
    if (use_van_pelt_)
    {
        set_param(status, names::B, B_, "count / minute");
        set_param(status, names::E, E_, "");
        set_param(status, names::S, S_, "");
        set_param(status, names::T, T_, "minute");
    }

    set_param(status, names::use_uniform_branching, use_uniform_branching_, "");
    if (use_uniform_branching_)
    {
        set_param(status, names::uniform_branching_rate,
                  uniform_branching_rate_, "count / minute");
    }

    set_param(status, names::use_flpl_branching, use_flpl_branching_, "");
    if (use_flpl_branching_)
    {
        set_param(status, names::flpl_branching_rate, flpl_branching_rate_,
                  "count / minute");
    }
}
} // namespace growth
