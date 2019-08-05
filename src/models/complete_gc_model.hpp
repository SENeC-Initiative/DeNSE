/*
 * complete_gc_model.hpp
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

#ifndef CMPLT_GC_MODEL_H
#define CMPLT_GC_MODEL_H

// elements
#include "GrowthCone.hpp"

// models
#include "extension_interface.hpp"
#include "steering_interface.hpp"
#include "direction_select_interface.hpp"


namespace growth
{

/**
 * @brief Full growth cone model template class
 *
 * The GrowthConeModel class inherits from :cpp:class:`GrowthCone`. It overrides
 * the pure virtual functions responsible for:
 *
 * - extension through the `ElType` class (`elongator_` member, inherits from
 *   :cpp:class:`ExtensionModel` interface),
 * - steering through the `SteerMethod` class (`steerer_` member, inherits from
 *   :cpp:class:`SteeringModel` interface),
 * - direction selection through the `DirSelMethod` class (`dir_selector_`
 *   member, inherits from :cpp:class:`DirectionSelectModel` interface).
 *
 * Most of the work is still performed in the :cpp:class:`GrowthCone` as
 * GrowthConeModel simply overrides :cpp:func:`compute_speed` and
 * :cpp:func:`compute_direction_probabilities` (both called in
 * :cpp:func:`GrowthCone::grow`), as well as :cpp:func:`select_direction`
 * (called in :cpp:func:`GrowthCone::make_move`).
 */
template <class ElType, class SteerMethod, class DirSelMethod>
class GrowthConeModel
  : public virtual GrowthCone
{
  private:
    std::shared_ptr<ElType> elongator_;
    std::shared_ptr<SteerMethod> steerer_;
    std::shared_ptr<DirSelMethod> dir_selector_;

  public:
    GrowthConeModel() = delete;
    // "default" constructor for model manager initial generation
    GrowthConeModel(const std::string &model);
    GrowthConeModel(const GrowthConeModel& copy);

    // factory function
    static std::shared_ptr<GrowthConeModel<ElType, SteerMethod, DirSelMethod>>
    create_gc_model(const std::string &model);

    GCPtr clone(BaseWeakNodePtr parent, NeuritePtr neurite,
                double distanceToParent, const BPoint &position,
                double angle) override final;

    void compute_speed(mtPtr rnd_engine, double substep) override final;

    void
    compute_direction_probabilities(std::vector<double> &directions_weights,
                                    double substep) override final;

    void
    select_direction(const std::vector<double> &directions_weights,
                     mtPtr rnd_engine, double &substep, double &new_angle,
                     stype &default_direction) override final;

    void prepare_for_split() override final;
    void after_split() override final;

    void set_status(const statusMap &) override final;
    void get_status(statusMap &) const override final;
    double get_state(const std::string& observable) const override final;
};


/*******************************************************************************
 * Implementation
 ******************************************************************************/

/**
 * @brief Initial constructor for the ghost models (in model_manager)
 *
 * This constructor must not be used directly, use :cpp:func:create_gc_model`.
 * @todo make private.
 */
template <class ElType, class SteerMethod, class DirSelMethod>
GrowthConeModel<ElType, SteerMethod, DirSelMethod>::GrowthConeModel(const std::string &model)
  : GrowthCone(model)
{}


/**
 * @brief Copy constructor to instantiate proper GrowthCone models
 *
 * This constructor must not be used directly, use :cpp:func:clone`.
 * @todo make private.
 */
template <class ElType, class SteerMethod, class DirSelMethod>
GrowthConeModel<ElType, SteerMethod, DirSelMethod>::GrowthConeModel(const GrowthConeModel<ElType, SteerMethod, DirSelMethod> &copy)
  : GrowthCone(copy)
{}


/**
 * @brief Factory function for the ghost models (in model_manager)
 *
 * Create the GrowthConeModel and initialize its members.
 */
template <class ElType, class SteerMethod, class DirSelMethod>
std::shared_ptr<GrowthConeModel<ElType, SteerMethod, DirSelMethod>>
GrowthConeModel<ElType, SteerMethod, DirSelMethod>::create_gc_model(
  const std::string &model)
{
    std::shared_ptr<GrowthConeModel<ElType, SteerMethod, DirSelMethod>> gc =
        std::make_shared<GrowthConeModel<ElType, SteerMethod, DirSelMethod>>(model);

    gc->elongator_ = std::make_shared<ElType>(gc, nullptr);

    gc->steerer_ = std::make_shared<SteerMethod>(gc, nullptr);

    gc->dir_selector_ = std::make_shared<DirSelMethod>(gc, nullptr);

    return gc;
}


/**
 * @brief Function to make a proper GrowthCone.
 *
 * Create the GrowthConeModel and initialize its members.
 */
template <class ElType, class SteerMethod, class DirSelMethod>
GCPtr GrowthConeModel<ElType, SteerMethod, DirSelMethod>::clone(
    BaseWeakNodePtr parent, NeuritePtr neurite, double distanceToParent,
    const BPoint &position, double angle)
{
    auto new_cone = std::make_shared<GrowthConeModel<ElType, SteerMethod, DirSelMethod>>(*this);
    
    // update topology
    int omp_id   = kernel().parallelism_manager.get_thread_local_id();
    new_cone->update_topology(parent, neurite, distanceToParent,
                              position, angle);

    // init the compontents
    new_cone->elongator_ = std::make_shared<ElType>(
        *(elongator_.get()), new_cone, neurite);

    new_cone->steerer_ = std::make_shared<SteerMethod>(
        *(steerer_.get()), new_cone, neurite);

    new_cone->dir_selector_ = std::make_shared<DirSelMethod>(
        *(dir_selector_.get()), new_cone, neurite);

    // update containing area
    new_cone->current_area_ =
        new_cone->using_environment_
        ? kernel().space_manager.get_containing_area(position)
        : "";

    if (new_cone->using_environment_)
    {
        new_cone->update_growth_properties(new_cone->current_area_);
    }

    return new_cone;
}


/*******************************************************************************
 * Model functions
 ******************************************************************************/

/**
 * @brief use the `elongator_` member to compute the speed.
 */
template <class ElType, class SteerMethod, class DirSelMethod>
void GrowthConeModel<ElType, SteerMethod, DirSelMethod>::compute_speed(
  mtPtr rnd_engine, double substep)
{
    move_.speed = elongator_->compute_speed(rnd_engine, substep);
}


/**
 * @brief use the `steerer_` member to evaluate the probability of each angle.
 */
template <class ElType, class SteerMethod, class DirSelMethod>
void GrowthConeModel<ElType, SteerMethod, DirSelMethod>::compute_direction_probabilities(
  std::vector<double> &directions_weights, double substep)
{
    steerer_->compute_direction_probabilities(
        directions_weights, filopodia_, substep, total_proba_, stuck_);
}


/**
 * @brief use the `dir_selector_` member to chose the next direction.
 */
template <class ElType, class SteerMethod, class DirSelMethod>
void GrowthConeModel<ElType, SteerMethod, DirSelMethod>::select_direction(
  const std::vector<double> &directions_weights, mtPtr rnd_engine,
  double &substep, double &new_angle, stype &default_direction)
{
    dir_selector_->select_direction(directions_weights, filopodia_,
        rnd_engine, total_proba_, interacting_, move_.angle, substep,
        move_.module, new_angle, stopped_, default_direction);
}


/**
 * @brief prepare the growth cone for a split.
 * @todo
 */
template <class ElType, class SteerMethod, class DirSelMethod>
void GrowthConeModel<ElType, SteerMethod, DirSelMethod>::prepare_for_split()
{
    elongator_->prepare_for_split();
    steerer_->prepare_for_split();
    dir_selector_->prepare_for_split();
}


/**
 * @brief clean up the growth cone after a split.
 * @todo
 */
template <class ElType, class SteerMethod, class DirSelMethod>
void GrowthConeModel<ElType, SteerMethod, DirSelMethod>::after_split()
{
    elongator_->after_split();
    steerer_->after_split();
    dir_selector_->after_split();
}


/**
 * @brief prepare the growth cone for a split.
 * @todo
 */
template <class ElType, class SteerMethod, class DirSelMethod>
void GrowthConeModel<ElType, SteerMethod, DirSelMethod>::set_status(
  const statusMap &status)
{
    // set the default parameters
    GrowthCone::set_status(status);

    // then set the models
    elongator_->set_status(status);
    steerer_->set_status(status);

    // direction selector may need to access all the properties
    statusMap base_status;
    GrowthCone::get_status(base_status);
    elongator_->get_status(base_status);
    steerer_->get_status(base_status);

    // update the status and set dire_selector
    base_status.insert(status.begin(), status.end());
    dir_selector_->set_status(base_status);
}


template <class ElType, class SteerMethod, class DirSelMethod>
void GrowthConeModel<ElType, SteerMethod, DirSelMethod>::get_status(
  statusMap &status) const
{
    GrowthCone::get_status(status);
    elongator_->get_status(status);
    steerer_->get_status(status);
    dir_selector_->get_status(status);

    // update observables
    std::vector<std::string> tmp;
    get_param(status, names::observables, tmp);
    elongator_->get_observables(tmp);
    steerer_->get_observables(tmp);
    dir_selector_->get_observables(tmp);
    // use set to keep only unique values
    std::set<std::string> obs(tmp.cbegin(), tmp.cend());
    obs.insert(observables_.cbegin(), observables_.cend());
    // set the full vector and update status
    tmp = std::vector<std::string>(obs.begin(), obs.end());
    set_param(status, names::observables, tmp, "");
}


template <class ElType, class SteerMethod, class DirSelMethod>
double GrowthConeModel<ElType, SteerMethod, DirSelMethod>::get_state(
  const std::string& observable) const
{
    double value = GrowthCone::get_state(observable);

    if (std::isnan(value))
    {
        value = elongator_->get_state(observable);
    }

    if (std::isnan(value))
    {
        value = steerer_->get_state(observable);
    }

    if (std::isnan(value))
    {
        value = dir_selector_->get_state(observable);
    }

    return value;
}

} // namespace growth
#endif
