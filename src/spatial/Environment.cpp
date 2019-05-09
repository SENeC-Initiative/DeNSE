/*
 * Environment.cpp
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

#include "Environment.hpp"


namespace growth
{

Environment::Environment(BMultiPolygonPtr environment)
  : environment_(environment)
{
    for (const auto& polygon : *(environment_.get()))
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


const BMultiPolygonPtr Environment::get_environment() const
{
    return environment_;
}


const BMultiLineString& Environment::get_boundary() const
{
    return boundary_;
}

} // namespace growth
