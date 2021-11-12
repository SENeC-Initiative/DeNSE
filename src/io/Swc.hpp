/*
 * Swc.hpp
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

#ifndef SWC_H
#define SWC_H

#include <fstream>
#include <iostream>
/// libgrowth includes
#include "elements_types.hpp"
#include "spatial_types.hpp"

namespace growth
{

class Swc
{
  private:
    static const int undefinedID = 0;
    static const int somaID      = 1;
    static const int axonID      = 2;
    // all dendrites are define basal for now.
    static const int b_dendriteID = 3;
    static const int a_dendriteID = 4;
    static const int forkID       = 5;
    // growth cone
    static const int endID = 6;
    // unused for now
    static const int customID = 7;
    std::ofstream swc_file_;
    unsigned int resolution_;

  public:
    Swc();
    Swc(std::string output_file, unsigned int resolution);
    Swc(Swc&& rhs);
    Swc& operator=(Swc &&rhs);
    void to_swc(const Neuron *, stype gid);
    void close_file();
    ~Swc();
};
} // namespace growth

#endif /* SWC_H */
