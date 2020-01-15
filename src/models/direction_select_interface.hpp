/*
 * direction_select_interface.hpp
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

#ifndef DIRSELMODEL_H
#define DIRSELMODEL_H

#include <string>
#include <vector>

// lib includes
#include "elements_types.hpp"

// elements include
#include "GrowthCone.hpp"
#include "Neurite.hpp"


namespace growth
{

/**
 * @brief Direction selection interface, from which all submodels must inherit.
 *
 * This class serves as a guide for developers wanting to create a new
 * direction selection model:
 *
 * - Default (copy) constructors must not be used and are therefore deleted.
 * - Correct constructor must get a valid instance of
 *   :cpp:class:`growth::GrowthCone` and :cpp:class:`growth::Neurite`
 *   shared pointer to store in `gc_weakptr_` and `neurite_ptr_`.
 * - If necessary, the correct constructor should update the `observables_` and
 *   override the :cpp:func:`growth::DirectionSelectModel::get_state` function.
 * - There is only one method that has to be implemented for a direction
 *   model to be valid, which is
 *   :cpp:func:`growth::DirectionSelectModel::select_direction`.
 *
 * See also
 * --------
 * :cpp:class:`growth::NMDirectionSelector` or
 * :cpp:class:`growth::RTDirectionSelector` for complete model
 * declarations.
 */
class DirectionSelectModel;
typedef std::shared_ptr<DirectionSelectModel> DSPtr;


class DirectionSelectModel
{
  protected:
    // GrowthCone must be destructed upon release by Neurite (don't keep it
    // alive) and we don't need to modify it (done through template)
    std::weak_ptr<GrowthCone> gc_weakptr_;
    // We must be able to modify Neurite (branching)
    std::shared_ptr<Neurite> neurite_ptr_;

    std::vector<std::string> observables_;

  public:
    DirectionSelectModel()                                 = delete;
    DirectionSelectModel(const DirectionSelectModel &copy) = delete;

    DirectionSelectModel(GCPtr gc, NeuritePtr neurite)
        : gc_weakptr_(GCPtr(gc))
        , neurite_ptr_(neurite){};

    DirectionSelectModel(const DirectionSelectModel &copy, GCPtr gc,
                         NeuritePtr neurite)
        : gc_weakptr_(gc)
        , neurite_ptr_(neurite)
        , observables_(copy.observables_){};

    virtual void select_direction(const std::vector<double> &directions_weights,
                                  const Filopodia &filo, mtPtr rnd_engine,
                                  double total_proba, bool interacting,
                                  double old_angle, double &substep,
                                  double &step_length, double &new_angle,
                                  bool &stopped, stype &default_direction) = 0;

    void get_observables(std::vector<std::string> &obs) const
    {
        obs.insert(obs.end(), observables_.begin(), observables_.end());
    }

    virtual void prepare_for_split(){};
    virtual void after_split(){};
    virtual void kernel_updated(){};
    virtual double get_state(const std::string &observable) const
    {
        return std::nan("");
    };
    virtual void set_status(const statusMap &status) = 0;
    virtual void get_status(statusMap &status) const = 0;
};

} // namespace growth

#endif
