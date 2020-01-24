/*
 * extension_interface.hpp
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

#ifndef ELMODEL_H
#define ELMODEL_H

#include <string>
#include <vector>

// lib includes
#include "elements_types.hpp"

// elements include
#include "GrowthCone.hpp"
#include "Neurite.hpp"


namespace growth
{

class ExtensionModel;
typedef std::shared_ptr<ExtensionModel> EMPtr;


/**
 * @brief Extension interface, from which all elongation models must inherit.
 *
 * This class serves as a guide for developers wanting to create a new
 * elongation model:
 *
 * - Default (copy) constructors must not be used and are therefore deleted.
 * - Correct constructor must get a valid instance of
 *   :cpp:class:`growth::GrowthCone` and :cpp:class:`growth::Neurite`
 *   shared pointer to store in `gc_weakptr_` and `neurite_ptr_`.
 * - If necessary, the correct constructor should update the `observables_` and
 *   override the :cpp:func:`growth::ExtensionModel::get_state` function.
 * - There is only one method that has to be implemented for an elongation
 *   model to be valid, which is
 *   :cpp:func:`growth::ExtensionModel::compute_speed`.
 *
 * See also
 * --------
 * :cpp:class:`growth::CstExtensionModel` or
 * :cpp:class:`growth::ResourceBasedExtensionModel` for complete model
 * declarations.
 */
class ExtensionModel
{
  protected:
    // GrowthCone must be destructed upon release by Neurite (don't keep it
    // alive) and we don't need to modify it (done through template)
    std::weak_ptr<GrowthCone> gc_weakptr_;
    // We must be able to modify Neurite (branching)
    std::shared_ptr<Neurite> neurite_ptr_;
    std::vector<std::string> observables_;

  public:
    ExtensionModel()                           = delete;
    ExtensionModel(const ExtensionModel &copy) = delete;

    ExtensionModel(GCPtr gc, NeuritePtr neurite)
        : gc_weakptr_(GCPtr(gc))
        , neurite_ptr_(neurite){};

    ExtensionModel(const ExtensionModel &copy, GCPtr gc, NeuritePtr neurite)
        : gc_weakptr_(gc)
        , neurite_ptr_(neurite)
        , observables_(copy.observables_){};

    virtual double compute_speed(mtPtr rnd_engine, double substep) = 0;

    virtual void update_speed(double speed_factor) = 0;
    virtual void update_local_speed(double area_factor) = 0;
    virtual double get_max_speed() const = 0;

    void get_observables(std::vector<std::string> &obs) const
    {
        obs.insert(obs.end(), observables_.begin(), observables_.end());
    }

    virtual void prepare_for_split(){};
    virtual void after_split(){};

    virtual double get_state(const std::string &observable) const
    {
        return std::nan("");
    };
    virtual bool set_status(const statusMap &status) = 0;
    virtual void get_status(statusMap &status) const = 0;
};

} // namespace growth

#endif
