/*
 * Neuron.cpp
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

#include "Neuron.hpp"

// c++ includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_set>

// elements includes
#include "GrowthCone.hpp"
#include "Node.hpp"

// lib includes
#include "growth_names.hpp"
#include "spatial_types.hpp"
#include "tools.hpp"

// kernel includes
#include "Skeleton.hpp"
#include "Swc.hpp"
#include "config_impl.hpp"
#include "kernel_manager.hpp"
#include "neuron_manager.hpp"


namespace growth
{

//--------------------------------------------------------------------------
// Neuron class

// Constructors and destructor

//! Default constructor
Neuron::Neuron(stype gid)
    : gid_(gid)
    , description_("standard_neuron")
    , soma_radius_(SOMA_RADIUS)
    , observables_(
          {"length", "speed", "num_growth_cones", "retraction_time", "stopped"})
    , use_actin_waves_(false)
    , aw_generation_step_(-1)
    , actin_content_(0)
    , next_actin_event_(0)
    , axon_angle_(0)
    , axon_angle_set_(false)
    , has_axon_(true)
    , polarization_strength_(POLA_STRENGTH)
    , axon_polarization_weight_(AXON_POLA_WEIGHT)
    , random_rotation_angles_(false)
    , rnd_angle_(0.)
    , uniform_(0., 1.)
    , normal_(0., 1.)
{
    // initialize the soma
    soma_ = std::make_shared<BaseNode>();
}


Neuron::~Neuron()
{
    // we need to delete the branching pointer of the neurites to let them die
    for (auto &neurite : neurites_)
    {
        // break circular dependency of neurite/branching model
        neurite.second->branching_model_ = nullptr;
        // clear all neurite containers to free growth cones and nodes
        neurite.second->nodes_.clear();
        neurite.second->growth_cones_.clear();
        neurite.second->growth_cones_tmp_.clear();
        neurite.second->growth_cones_inactive_.clear();
        neurite.second->growth_cones_inactive_tmp_.clear();
    }

    neurites_.clear();
}


//###################################################
//                  Init functions
//###################################################

void Neuron::init_status(
  const statusMap &status,
  const std::unordered_map<std::string, statusMap> &neurite_statuses,
  mtPtr rnd_engine)
{
    // set growth cone model and actin wave
    set_status(status);

    int num_neurites = 0;

    get_param(status, names::num_neurites, num_neurites);
    get_param(status, names::growth_cone_model, growth_cone_model_);
    get_param(status, names::polarization_strength,
              polarization_strength_);
    get_param(status, names::axon_polarization_weight,
              axon_polarization_weight_);
    get_param(status, names::has_axon, has_axon_);

    // send the soma to the space_manager
    int omp_id = kernel().parallelism_manager.get_thread_local_id();

    try
    {
        kernel().space_manager.add_object(
            get_position(), get_position(), 2 * soma_radius_, 0.,
            0., std::make_tuple(gid_, std::string(""), 0UL, 0UL),
            nullptr, omp_id);
    }
    catch (...)
    {
        std::throw_with_nested(
            std::runtime_error("Passed from `Neuron::init_status`."));
    }

    // prepare the neurites
    std::unordered_set<std::string> neurite_names;

    bool names_set = get_param(status, names::neurite_names,
                               neurite_names);

    if (names_set)
    {
        if (neurite_names.size() != num_neurites)
        {
            throw InvalidParameter("`neurite_names` must contain one "
                                   "entry per neurite.",
                                   __FUNCTION__, __FILE__, __LINE__);
        }

        if (not has_axon_ and
            neurite_names.find("axon") != neurite_names.end())
        {
            throw InvalidParameter("`neurite_names` contains an 'axon' "
                                   "entry but 'has_axon' is set to "
                                   "False.", __FUNCTION__, __FILE__,
                                   __LINE__);
        }

        if (neurite_statuses.find("dendrites") == neurite_statuses.end())
        {
            for (auto p : neurite_statuses)
            {
                if (neurite_names.find(p.first) == neurite_names.end())
                {
                    throw InvalidParameter(
                        "`neurite_names` and the parameters contain "
                        "different  neurite names.", __FUNCTION__,
                        __FILE__, __LINE__);
                }
            }
        }
    }

    std::unordered_map<std::string, double> nas;
    bool angles_set = get_param(status, names::neurite_angles, nas);

    if (angles_set)
    {
        if (nas.size() != num_neurites)
        {
            throw InvalidParameter("`neurite_angles` must contain one "
                                   "entry per neurite, got " +
                                   std::to_string(nas.size()) + ".",
                                   __FUNCTION__, __FILE__, __LINE__);
        }

        if (neurite_statuses.find("dendrites") == neurite_statuses.end())
        {
            for (auto p : neurite_statuses)
            {
                if (nas.find(p.first) == nas.end())
                {
                    throw InvalidParameter(
                        "`neurite_angles` contain invalid neurite names.",
                        __FUNCTION__, __FILE__, __LINE__);
                }
            }
        }

        neurite_angles_ = nas;
    }

    GCPtr gc_model;
    // get growth cone model from the model manager;
    if (growth_cone_model_ != "")
    {
        gc_model = kernel().model_manager.get_model(growth_cone_model_);
    }
    else
    {
        gc_model = kernel().model_manager.get_default_model();
    }

    // initialize the soma giving the position of neuron
    //@TODO move the neuron in a new position with all the branches
    auto pos_x = status.find("x");
    if (pos_x != status.end())
    {
        double x, y;
        get_param(status, "x", x);
        get_param(status, "y", y);
        soma_->set_position(BPoint(x, y));
    }
    else
    {
        throw InvalidParameter("Position was not set.", __FUNCTION__,
                               __FILE__, __LINE__);
    }

    // prepare random angle and cg models if necessary
    if (random_rotation_angles_)
    {
        rnd_angle_ = 2 * M_PI * uniform_(*(rnd_engine.get()));
    }

    GCPtr axon_gc     = gc_model;
    GCPtr dendrite_gc = gc_model;

    std::string model_name;
    statusMap astatus, dstatus;

    // create axon
    if (has_axon_ and num_neurites > 0)
    {
        astatus = status;
        auto it_status = neurite_statuses.find("axon");

        if (it_status != neurite_statuses.end())
        {
            for (auto entry : it_status->second)
            {
                astatus[entry.first] = entry.second;
            }
        }

        auto it2 = astatus.find(names::growth_cone_model);
        if (it2 != astatus.end())
        {
            get_param(astatus, names::growth_cone_model, model_name);
            axon_gc = kernel().model_manager.get_model(model_name);
        }

        new_neurite("axon", "axon", axon_gc, rnd_engine);
        set_neurite_status("axon", astatus);
    }

    // create dendrites
    for (auto name : neurite_names)
    {
        dstatus = status;

        if (name != "axon")
        {
            // check for placeholder params
            auto it_dendrites = neurite_statuses.find("dendrites");

            if (it_dendrites != neurite_statuses.end())
            {
                for (auto p : it_dendrites->second)
                {
                    dstatus[p.first] = p.second;
                }
            }

            // get for specific dendrite parameters
            auto it_status = neurite_statuses.find(name);

            if (it_status != neurite_statuses.end())
            {
                for (auto entry : it_status->second)
                {
                    dstatus[entry.first] = entry.second;
                }
            }

            // get dendrite growth cone model
            auto it2 = dstatus.find(names::growth_cone_model);

            if (it2 != dstatus.end())
            {
                get_param(dstatus, names::growth_cone_model, model_name);
                dendrite_gc = kernel().model_manager.get_model(model_name);
            }

            new_neurite(name, "dendrite", dendrite_gc, rnd_engine);
            set_neurite_status(name, dstatus);
        }
    }
}


void Neuron::update_angles(const std::unordered_map<std::string, double> &angles)
{
    for (auto entry : angles)
    {
        neurite_angles_[entry.first] = rnd_angle_ + entry.second;
    }
}


void Neuron::initialize_next_event(mtPtr rnd_engine)
{
    if (use_actin_waves_ and next_actin_event_ == 0)
    {
        next_actin_event(rnd_engine);
    }
    for (auto &neurite : neurites_)
    {
        neurite.second->branching_model_->initialize_next_event(rnd_engine);
    }
}


void Neuron::finalize()
{
    for (auto &neurite : neurites_)
    {
        neurite.second->finalize();
    }
}


//###################################################
//                  Growth
//###################################################

/**
 * Function that simulates the growth of the neuron.
 */
void Neuron::grow(mtPtr rnd_engine, stype current_step, double substep)
{
    // if we use actin waves, tell each neurite to update them
    if (use_actin_waves_)
    {
        // new actin wave generation
        // @todo: update this to check substep
        if (current_step == aw_generation_step_)
        {
            next_actin_event(rnd_engine);
            // pick random neurite
            stype idx_neurite =
                neurites_.size() * uniform_(*(rnd_engine).get());
            auto rnd_neurite = std::next(std::begin(neurites_), idx_neurite);
            rnd_neurite->second->start_actin_wave(actin_content_);
        }
        // actin wave update
        for (const auto &neurite : neurites_)
        {
            neurite.second->update_actin_waves(rnd_engine, substep);
        }
    }

    // For each neurite of the neuron, apply growth
    for (auto &neurite : neurites_)
    {
        if (neurite.second->active_)
        {
            try
            {
                neurite.second->grow(rnd_engine, current_step, substep);
            }
            catch (...)
            {
                std::throw_with_nested(
                    std::runtime_error("Passed from `Neuron::grow`."));
            }
        }
    }
}


/**
 * @brief branch according to the event
 *
 * Try to create a new growth cone, either through splitting or lateral
 * branching and return whether the branching was sucessful.
 */
bool Neuron::branch(mtPtr rnd_engine, const Event &ev)
{
    std::string name_neurite = std::get<edata::NEURITE>(ev);
    NeuritePtr neurite       = neurites_[name_neurite];

    try
    {
        return neurite->branching_model_->branching_event(rnd_engine, ev);
    }
    catch (...)
    {
        std::throw_with_nested(
            std::runtime_error("Passed from `Neuron::branch`."));
    }
}


// Neurite

/*! Function creating a new ``Neurite`` object for the ``Neuron``.
 *
 *  \param name: string - name of the neurite.
 *  \param neurite_type: string, optional (default: "dendrite") -
 *      type of neurite (either "dendrite" or "axon").
 *  \return string name of the created ``Neurite`` object.
 *  the neurite constructor will initialize the ''Growth Cone'' with
 *  the proper model.
 */
std::string Neuron::new_neurite(const std::string &name,
                                const std::string &neurite_type,
                                const GCPtr gc_model, mtPtr rnd_engine)
{
    int omp_id = kernel().parallelism_manager.get_thread_local_id();

    // update the neurite in the neuron vector
    NeuronWeakPtr my_weak_ptr =
        std::dynamic_pointer_cast<Neuron>(shared_from_this());
    neurites_.insert(
        {name, std::make_shared<Neurite>(name, neurite_type, growth_cone_model_,
                                         my_weak_ptr)});

    //#####################################
    // add first growth cone to the neurite
    //#####################################
    // initialize parameters
    double angle = 0;
    BPoint cone_start_point;
    bool contained = false;

    // a minimal trophism approach: set the neurite on the other side of the
    // neuron
    // then ensure that the first cone it's not created outside the environment.
    if (name == "axon")
    {
        if (neurite_angles_.find(name) != neurite_angles_.end())
        {
            angle           = neurite_angles_[name] + rnd_angle_;
            BPoint position = get_position();
            cone_start_point =
                BPoint(position.x() + soma_radius_ * cos(angle),
                       position.y() + soma_radius_ * sin(angle));

            contained = kernel().space_manager.env_contains(cone_start_point);

            int num_test = 0;

            while (not contained and num_test < 1000)
            {
                angle += sgn(uniform_(*(rnd_engine).get()) - 0.5) * 0.1;
                axon_angle_             = angle;
                neurite_angles_["axon"] = axon_angle_;
                cone_start_point =
                    BPoint(position.x() + soma_radius_ * cos(angle),
                           position.y() + soma_radius_ * sin(angle));
                contained =
                    kernel().space_manager.env_contains(cone_start_point);

                num_test++;
            }

            if (num_test == 1000)
            {
                throw std::runtime_error(
                    "Could not initialize neurite, please "
                    "that your neuron does not lay outside "
                    "of the environment.");
            }
        }
        else
        {
            do
            {
                if (axon_angle_set_)
                {
                    angle = fmod(axon_angle_, 2 * M_PI);
                }
                else
                {
                    angle = uniform_(*(rnd_engine).get()) * 2 * M_PI;
                }
                axon_angle_             = angle;
                neurite_angles_["axon"] = axon_angle_;
                BPoint position         = get_position();
                cone_start_point =
                    BPoint(position.x() + soma_radius_ * cos(angle),
                           position.y() + soma_radius_ * sin(angle));
                contained =
                    kernel().space_manager.env_contains(cone_start_point);
                if (axon_angle_set_ and not contained)
                {
                    throw InvalidParameter(
                        "Invalid axon angle for neuron: growth "
                        "cone is outside the environment.",
                        __FUNCTION__, __FILE__, __LINE__);
                }
            } while (not contained);
        }

        neurites_[name]->init_first_node(soma_, get_position(), name,
                                         soma_radius_);
    }
    else
    {
        if (neurite_angles_.find(name) != neurite_angles_.end())
        {
            angle           = neurite_angles_[name] + rnd_angle_;
            BPoint position = get_position();
            cone_start_point =
                BPoint(position.x() + soma_radius_ * cos(angle),
                       position.y() + soma_radius_ * sin(angle));

            contained = kernel().space_manager.env_contains(cone_start_point);

            char sgn = uniform_(*(rnd_engine).get()) > 0.5 ? 1 : -1;

            while (not contained)
            {
                angle += sgn * 0.1;
                neurite_angles_[name] = angle;
                cone_start_point =
                    BPoint(position.x() + soma_radius_ * cos(angle),
                           position.y() + soma_radius_ * sin(angle));
                contained =
                    kernel().space_manager.env_contains(cone_start_point);
            }
        }
        else
        {
            double other_angle, largest_dtheta, dtheta, polarization_angle;
            std::vector<double> angles;

            for (auto na : neurite_angles_)
            {
                angles.push_back(na.second);
            }

            std::sort(angles.begin(), angles.end());

            // we iterate over the angular apertures between the existing
            // angles until we find one where we can insert the new neurite
            do
            {
                polarization_angle = 0.;
                largest_dtheta     = 0.;

                if (angles.size() == 1)
                {
                    polarization_angle = angles[0];
                    largest_dtheta     = 2 * M_PI;
                }
                else
                {
                    for (unsigned int i = 0; i < angles.size(); i++)
                    {
                        if (i != angles.size() - 1)
                        {
                            other_angle = angles[i + 1];
                        }
                        else
                        {
                            other_angle = angles[0] + 2 * M_PI;
                        }

                        dtheta = std::abs(other_angle - angles[i]);
                        // Be more explicit in this function: the axon has
                        // different direction of resting neurites, which is the
                        // weigth scale of axon polarization weigt
                        if (has_axon_ and
                            (std::abs(other_angle - neurite_angles_["axon"]) <
                                 0.0001 or
                             std::abs(angles[i] - neurite_angles_["axon"]) <
                                 0.0001))
                        {
                            if (dtheta / (1 + axon_polarization_weight_) >
                                largest_dtheta)
                            {
                                polarization_angle = angles[i];
                                largest_dtheta     = dtheta;
                            }
                        }
                        else if (dtheta > largest_dtheta)
                        {
                            polarization_angle = angles[i];
                            largest_dtheta     = dtheta;
                        }
                    }
                }

                angle =
                    fmod(polarization_angle +
                             largest_dtheta *
                                 (0.5 + (uniform_(*(rnd_engine).get()) - 0.5) /
                                            polarization_strength_),
                         2 * M_PI);

                neurite_angles_[name] = angle;

                BPoint position = get_position();
                cone_start_point =
                    BPoint(position.x() + soma_radius_ * cos(angle),
                           position.y() + soma_radius_ * sin(angle));

            } while (not kernel().space_manager.env_contains(cone_start_point));
        }

        neurites_[name]->init_first_node(soma_, get_position(), name,
                                         soma_radius_);
    }

    // initialize the neurite firstNode, a copy of the soma with child!
    // @todo: clean up this mess!!

    // eventually create the growth cone, cloning the default model.
    neurites_[name]->growth_cone_model_ = gc_model->get_model_name();
    GCPtr first_gc =
        gc_model->clone(neurites_[name]->get_first_node(), neurites_[name],
                        soma_radius_, cone_start_point, angle);

    neurites_[name]->add_cone(first_gc);

    // set initial and current diameter
    first_gc->set_diameter(neurites_[name]->get_first_node()->diameter_);

    // then reset the branch with the right initial position
    first_gc->set_first_point(cone_start_point, soma_radius_);
    // and adjust the topology of the new growth cone
    neurites_[name]->nodes_[0]->children_.push_back(first_gc);
    neurites_[name]->update_tree_structure(neurites_[name]->get_first_node());

    if (name == "axon")
    {
        has_axon_ = true;
    }

    return name;
}


void Neuron::delete_neurites(const std::vector<std::string> &names)
{
    if (names.empty())
    {
        neurites_.clear();
        has_axon_ = false;
    }
    else
    {
        for (const std::string &neurite_name : names)
        {
            auto it = neurites_.find(neurite_name);

            if (it != neurites_.end())
            {
                NeuritePtr neurite = neurites_[neurite_name];
                // break circular dependency of neurite/branching model
                neurite->branching_model_ = nullptr;
                // clear all neurite containers to free growth cones and nodes
                neurite->nodes_.clear();
                neurite->growth_cones_.clear();
                neurite->growth_cones_tmp_.clear();
                neurite->growth_cones_inactive_.clear();
                neurite->growth_cones_inactive_tmp_.clear();

                // remove from neurites
                neurites_.erase(it);

                if (neurite_name == "axon")
                {
                    has_axon_ = false;
                }
            }
            else
            {
                throw std::runtime_error("Neurite '" + neurite_name +
                                         "' does "
                                         "not exist.");
            }
        }
    }
}


void Neuron::next_actin_event(mtPtr rnd_engine)
{
    // aw_generation_step_ += (stype)1. / actin_freq_;
    // printf( "next actin event at %lu", aw_generation_step_);
}


//###################################################
//                  Getter/setter functions
//###################################################

void Neuron::set_status(const statusMap &status)
{
    double sr, ad, dd, aa;
    bool bsr, bad, bdd, baa;

    bool has_neurites = neurites_.size() > 0;
    bool sim_started  = (kernel().simulation_manager.get_time() != Time());

    bsr = get_param(status, names::soma_radius, sr);
    if (bsr and sr != soma_radius_)
    {
        if (has_neurites)
        {
            throw InvalidArg("Cannot change the soma radius for neuron " +
                                 std::to_string(gid_) +
                                 " after neurites have been created.",
                             __FUNCTION__, __FILE__, __LINE__);
        }
        else
        {
            soma_radius_ = sr;
        }
    }

    //@TODO change the growth cone_model during the growth
    std::string gc_model;
    bool gc_change = get_param(status, names::growth_cone_model, gc_model);
    if (gc_change and gc_model != growth_cone_model_)
    {
        if (has_neurites)
        {
            throw InvalidArg("Cannot change the growth cone model for neuron " +
                                 std::to_string(gid_) +
                                 " after neurites have been created.",
                             __FUNCTION__, __FILE__, __LINE__);
        }
        else
        {
            growth_cone_model_ = gc_model;
        }
    }

    get_param(status, names::random_rotation_angles, random_rotation_angles_);
    get_param(status, names::polarization_strength, polarization_strength_);
    get_param(status, names::axon_polarization_weight,
              axon_polarization_weight_);

    auto nit = neurites_.find("axon");
    baa      = get_param(status, names::axon_angle, aa);
    baa *= (aa != axon_angle_);
    if (baa)
    {
        if (nit != neurites_.end())
        {
            throw InvalidArg("Cannot change the axon angle for neuron " +
                                 std::to_string(gid_) +
                                 " after axon has been created.",
                             __FUNCTION__, __FILE__, __LINE__);
        }
        else
        {
            axon_angle_     = aa;
            axon_angle_set_ = true;
        }
    }

    auto it = status.find("x");
    if (it != status.end())
    {
        double x, y;
        get_param(status, "x", x);
        get_param(status, "y", y);

        BPoint p = BPoint(x, y);

        if (has_neurites and
            not kernel().space_manager.is_close(p, get_position()))
        {
            throw InvalidArg("Cannot change the position of neuron " +
                                 std::to_string(gid_) +
                                 " after neurites have been created.",
                             __FUNCTION__, __FILE__, __LINE__);
        }
        else
        {
            soma_->set_position(p);
        }
    }

    //@TODO implement actin waves
    // auto use_actin_waves_old = use_actin_waves_;
    // get_param(status, names::use_actin_waves, use_actin_waves_);
    // if (use_actin_waves_ != use_actin_waves_old)
    //{
    // next_actin_event_ = 0;
    //}
    // generate time for first actin wave it it wasn't set before
}


void Neuron::set_neurite_status(const std::string &neurite,
                                const statusMap &status)
{
    auto it = neurites_.find(neurite);

    if (it != neurites_.end())
    {
        // we are setting the properties of a specific neurite
        it->second->set_status(status);
    }
    else
    {
        // set the properties of all neurites for which the type is `neurite`
        for (auto &neur : neurites_)
        {
            if (neur.second->neurite_type_ == neurite)
            {
                neur.second->set_status(status);
            }
        }
    }

    // update max_resol
    double max_resol = std::numeric_limits<double>::max();

    for (auto &neurite : neurites_)
    {
        max_resol = std::min(max_resol, neurite.second->get_max_resol());
    }

    kernel().neuron_manager.set_max_resol(gid_, max_resol);
}


/**
 * @brief Get the current value of one of the observables
 */
double Neuron::get_state(const std::string &observable) const
{
    double value = 0.;

    for (const auto &neurite : neurites_)
    {
        value += neurite.second->get_state(observable);
    }

    return value;
}


/**
 * @brief Get the current value of one of the observables
 */
double Neuron::get_state(const std::string &observable,
                         std::string &unit) const
{
    double value = 0.;

    for (const auto &neurite : neurites_)
    {
        value += neurite.second->get_state(observable, unit);
    }

    return value;
}


void Neuron::get_status(statusMap &status) const
{
    set_param(status, names::soma_radius, soma_radius_, "micrometer");
    set_param(status, names::observables, observables_, "");
    set_param(status, names::axon_angle, axon_angle_, "rad");
    set_param(status, names::description, description_, "");
    set_param(status, names::has_axon, has_axon_, "");
    set_param(status, names::polarization_strength, polarization_strength_, "");
    set_param(status, names::axon_polarization_weight,
              axon_polarization_weight_, "");
    set_param(status, names::neurite_angles, neurite_angles_, "rad");
    set_param(status, names::random_rotation_angles, random_rotation_angles_,
              "");
    set_param(status, names::growth_cone_model, growth_cone_model_, "");

    // neurite names
    std::vector<std::string> nnames;

    for (auto it : neurites_)
    {
        nnames.push_back(it.first);
    }

    set_param(status, names::neurite_names, nnames, "");

    // set position
    BPoint pos = soma_->get_position();
    set_param(status, "x", pos.x(), "micrometer");
    set_param(status, "y", pos.y(), "micrometer");
}


void Neuron::get_neurite_status(statusMap &status, std::string neurite,
                                const std::string &level)
{
    auto it = neurites_.find(neurite);

    if (it == neurites_.end())
    {
        if (neurite == "dendrite" or neurite == "dendrites")
        {
            for (const auto &neurite : neurites_)
            {
                if (neurite.second->neurite_type_ == "dendrite")
                {
                    neurite.second->get_status(status, level);
                    break;
                }
            }
        }
        else
        {
            throw InvalidArg("Unknown neurite `" + neurite + "`.",
                             __FUNCTION__, __FILE__, __LINE__);
        }
    }
    else
    {
        it->second->get_status(status, level);
    }
}


void Neuron::update_kernel_variables()
{
    for (const auto &neurite : neurites_)
    {
        neurite.second->update_kernel_variables();
    }
}


double Neuron::get_soma_radius() const { return soma_radius_; }


BaseNodePtr Neuron::get_soma() const { return soma_; }


int Neuron::get_num_neurites() const { return neurites_.size(); }


std::string Neuron::get_gc_model() const { return growth_cone_model_; }


BPoint Neuron::get_position() const { return soma_->get_position(); }


stype Neuron::get_gid() const { return gid_; }


bool Neuron::has_axon() const { return has_axon_; };


bool Neuron::is_neurite(const std::string &neurite)
{
    return neurites_.find(neurite) != neurites_.end();
}


NeuriteWeakPtr Neuron::get_neurite(const std::string &name) const
{
    return NeuriteWeakPtr(neurites_.at(name));
}

} // namespace growth
