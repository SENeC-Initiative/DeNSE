/*
 * Area.cpp
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

#include "Area.hpp"


namespace growth
{

Area::Area(BMultiPolygonPtr area, double height, const std::string &name,
           std::unordered_map<std::string, double> properties)
  : height_(height)
  , name_(name)
  , properties_(properties)
  , shape_(area)
{
    for (const auto& polygon : *(area.get()))
    {
        boundary_.push_back(BLineString());

        auto &ls = boundary_.back();
        ls.insert(ls.end(), polygon.outer().begin(), polygon.outer().end());

        for (const auto& inner : polygon.inners())
        {
            boundary_.push_back(BLineString());

            auto &l = boundary_.back();
            l.insert(l.end(), inner.begin(), inner.end());
        }
    }
}


/**
 * @brief return the value of a property modifier for the substrate of the area
 *
 * @param name : name of a property
 *
 * @returns 1. if the name is not contained in the `Area.properties_` (will not
 * change the property) or a positive real that will modify the behavior of the
 * growth on this area.
 */
double Area::get_property(const std::string &name) const
{
    auto it = properties_.find(name);

    if (it == properties_.end())
    {
        return 1.;
    }

    return properties_.at(name);
}


double Area::get_height() const { return height_; }


std::string Area::get_name() const { return name_; }


void Area::get_properties(
    std::unordered_map<std::string, double> &properties) const
{
    properties.insert(properties_.cbegin(), properties_.cend());
}


const BMultiPolygonPtr Area::get_area() const
{
    return shape_;
}


const BMultiLineString& Area::get_boundary() const
{
    return boundary_;
}

} // namespace growth
