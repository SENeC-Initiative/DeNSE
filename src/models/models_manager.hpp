/*
 * models_manager.hpp
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

#ifndef MODELS_MANAGER_H
#define MODELS_MANAGER_H

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "GrowthCone.hpp"


namespace growth
{

class GrowthCone;
typedef std::shared_ptr<GrowthCone> GCPtr;


class ModelManager
{
  public:
    const std::vector<std::string> extension_types;
    const std::vector<std::string> steering_methods;
    const std::vector<std::string> direction_selection;

    const std::unordered_map<std::string, std::string> specials;

    const std::unordered_map<std::string, std::string> abbrev_to_full;
    const std::unordered_map<std::string, std::string> full_to_abbrev;

    const std::string default_model;

  private:
    std::unordered_map<std::string, GCPtr> models_;

  public:
    ModelManager();

    void init_models();

    std::vector<std::string> get_extension_types();
    std::vector<std::string> get_steering_methods();
    std::vector<std::string> get_direction_selection_methods();

    void get_models(std::unordered_map<std::string, std::string> &models, bool abbrev);
    GCPtr get_model(const std::string& name);
    GCPtr get_default_model();
};

}
#endif /* MODELS_MANAGER_H */
