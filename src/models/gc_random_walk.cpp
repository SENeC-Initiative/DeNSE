#include "gc_random_walk.hpp"

// C++ includes
#include <assert.h>
#include <cmath>
#include <memory>

// elements includes
#include "Neurite.hpp"

// libgrowth includes
#include "config_impl.hpp"
#include "growth_names.hpp"
#include "growth_time.hpp"
#include "kernel_manager.hpp"


namespace growth
{

/**
 * This class provides a Random walk model with:
 *
 * - an history weighted direction, which is parametrized into memory_
 * - an exponential decaying correlation between gaussian, which is
 *   parametrized into corr_rw_
 * - a tunable sensing angle: the variance of the RW which determines the
 *   environment interaction.
 *
 * If any more sophisticated model for elongation isn't offered in
 * growth::set_status then a gaussian distribution for the speed is adopted.
 */
GrowthCone_RandomWalk::GrowthCone_RandomWalk()
    : GrowthCone()
    /* Model variables:: will be passed to statusMap
     defaults are initialized:
    */
    , deterministic_angle_(0)
    , persistence_length_(0)
    , corr_rw_{RW_DELTA_CORR, 0, 0, 0, 0}
    , memory_{RW_MEMORY_TAU, 1, 0, 0}
{
}


GrowthCone_RandomWalk::GrowthCone_RandomWalk(const GrowthCone_RandomWalk &copy)
    : GrowthCone(copy)
    , deterministic_angle_(0)
    , corr_rw_(copy.corr_rw_)
    , memory_(copy.memory_)
{
    memory_.effective_angle = move_.angle;
}


// Clone function

GCPtr GrowthCone_RandomWalk::clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                                   double distanceToParent,
                                   std::string binaryID, const Point &position,
                                   double angle)
{
#ifndef NDEBUG
    printf(" It's calling RandomWalk->clone! with direction %f\n", angle);
#endif
    auto newCone = std::make_shared<GrowthCone_RandomWalk>(*this);
    newCone->update_topology(parent, neurite, distanceToParent, binaryID,
                             position, angle);
    return newCone;
}


/**
 * @brief Initialize the Random Walker
 * Initialize memory and gaussian correlation algorithms.
 * The dimensional parameters set by the user are translated into coefficient
 * with respect to the average speed and the
 * time resolution of the simulator
 *
 *
 */
void GrowthCone_RandomWalk::initialize_RW()
{
    double average_step = timestep_ * average_speed_;
    // set the memory parameter.

    { // initializing memory random walk
        if (memory_.tau > 0)
        {
            memory_.alpha_coeff = exp(-average_step / memory_.tau);
        }
        else
        {
            memory_.alpha_coeff = 1;
        }
    }

    // initializing persistance random walk
    if (corr_rw_.tau > 0)
    {
        corr_rw_.f_coeff = exp(-average_step / corr_rw_.tau);
    }
    else
    {
        corr_rw_.f_coeff = 0;
    }
    corr_rw_.sqrt_f_coeff = sqrt(1 - corr_rw_.f_coeff * corr_rw_.f_coeff);

    memory_.effective_angle = move_.angle;

// setting the diffusion process for the persistence length;
#ifndef NDEBUG
    printf("diffusive_sigma: %f\n"
           "correlation f coefficient: %f\n"
           "memory alpha coefficient: %f\n",
           rw_sensing_angle_, corr_rw_.f_coeff, memory_.alpha_coeff);
#endif
}


/**
 * @brief Compute the new angle with persistence and memory length
 *
 * The persistence length is a random walk with coloured noise.
 * the correlation coefficient 'f' of the noise
 *
 * Given that 'g_i' is a normal distributed value with mean zero and unitary
 * sigma.
 * the recursive algorithm is:
 * r_n+1 = f*r_n + sqrt(1-f*f) * g_n
 * r_0   = g_0
 * In the code, the r variable is stored in ``corr_rw_.det_delta``.
 *
 * It's possible to obtain a gaussian with mean mu and variance sigma with:
 * x_n = mu + sigma* r_n
 * This process is markovian in the phasespace of x_n and r_n; it is stored in
 * ``move_.angle`` in the code.
 *
 * The memory length is a another type of persistence length and consists of
 * changing the mean of the diffusive process at each step, in respect to the
 * history
 * of the angle.
 * The weight of each angle decreases exponentially with a coefficient alpha,
 * this coefficient
 * is a simple memory kernel which the user can choose the carachteristic
 * length.
 * the recursive algorithm is:
 * y_n+1 = (A_n *alpha * y_n +x_n )/(alpha*A_n +1)
 * x_n+1 = g_n+1 + y_n+1
 * A_n+1 = alpha *A_n +1
 * the algortithm requires an initialization:
 * A_0 = 1
 * x_0 = move_angle
 * y_0 = move_angle
 *
 * whit alpha = 0 and f =0 it retrieves a diffusive process.
 *
 */
void GrowthCone_RandomWalk::compute_new_direction(mtPtr rnd_engine)
{
    // update the memory length algorithm
    memory_.effective_angle =
        (move_.angle + memory_.alpha_temp_sum * memory_.alpha_coeff *
                           memory_.effective_angle) /
        (memory_.alpha_coeff * memory_.alpha_temp_sum + 1);

    memory_.alpha_temp_sum = memory_.alpha_coeff * memory_.alpha_temp_sum + 1;

    // compute the persistence length angle deviation based on the random angle
    // delta_angle of the current step and its previous value.
    corr_rw_.det_delta = corr_rw_.sqrt_f_coeff * delta_angle_ // random part
                         + corr_rw_.f_coeff * corr_rw_.det_delta;

    // choose the new angle
    double default_angle = move_.angle + move_.sigma_angle * delta_angle_;
    double new_angle =
        memory_.effective_angle + move_.sigma_angle * corr_rw_.det_delta;

    // check that this step is allowed, otherwise move towards default_angle
    Point target_pos =
        Point(geometry_.position.at(0) + cos(new_angle) * move_.module,
              geometry_.position.at(1) + sin(new_angle) * move_.module);

    while (not kernel().space_manager.env_contains(target_pos))
    {
        new_angle = 0.5 * (new_angle + default_angle);
        target_pos =
            Point(geometry_.position.at(0) + cos(new_angle) * move_.module,
                  geometry_.position.at(1) + sin(new_angle) * move_.module);
    }
    move_.angle = new_angle;
}


void GrowthCone_RandomWalk::prepare_for_split()
{
    // memory_.alpha_temp_sum=0;
}


void GrowthCone_RandomWalk::after_split()
{
    memory_.effective_angle = move_.angle;
}


//#############################################
//          Get Set status
//#############################################
void GrowthCone_RandomWalk::set_status(const statusMap &status)
{
    GrowthCone::set_status(status);

    get_param(status, names::rw_memory_tau, memory_.tau);
    get_param(status, names::rw_delta_corr, corr_rw_.tau);
    if (get_param(status, names::rw_persistence_length, persistence_length_))
    {
        _persistence_set_ = true;
        memory_.tau       = persistence_length_;
        corr_rw_.tau      = persistence_length_;
    }
    else
    {
        persistence_length_ = -1;
    }
    assert(persistence_length_ != 0);

    initialize_RW();
}


void GrowthCone_RandomWalk::get_status(statusMap &status) const
{
    GrowthCone::get_status(status);
    set_param(status, names::rw_delta_corr, corr_rw_.tau);
    set_param(status, names::rw_memory_tau, memory_.tau);
}
}
