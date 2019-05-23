/*
 * growth_space.cpp
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

#include "growth_space.hpp"

namespace growth
{

const int Space::DEFAULT_DIM = 2;
int Space::dim_              = Space::DEFAULT_DIM;

Space::Space() {}

int Space::get_dimension() { return dim_; }

void Space::set_dimension_(int dim) { dim_ = dim; }

Space::Shape::Shape() {}
} // namespace growth
