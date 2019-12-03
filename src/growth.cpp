/*
 * growth.cpp
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

#include "config.hpp"
#include "module.hpp"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char **argv)
{
    int num_neurons = 6;

    growth::init_growth_(&argc, &argv);
    growth::reset_kernel_();

    growth::statusMap kernelMap;
    kernelMap.insert({"seeds", growth::Property({33.}, "")});
    kernelMap.insert({"num_local_threads", growth::Property(1, "")});

    growth::set_kernel_status_(kernelMap, "ID");


    growth::statusMap mainMap;
    mainMap.insert({"x", growth::Property(0.0, "micrometer")});
    mainMap.insert({"y", growth::Property(0.0, "micrometer")});
    mainMap.insert({"growth_cone_model", growth::Property("default", "")});
    mainMap.insert({"use_uniform_branching", growth::Property(true, "")});
    mainMap.insert({"use_van_pelt", growth::Property(true, "")});
    mainMap.insert({"B", growth::Property(4.2, "count / minute")});
    mainMap.insert({"E", growth::Property(0.05, "")});
    mainMap.insert({"S", growth::Property(2., "")});
    mainMap.insert({"T", growth::Property(0.1, "minute")});
    mainMap.insert({"num_neurites", growth::Property(2, "")});
    mainMap.insert(
        {"speed_growth_cone", growth::Property(20., "micrometer / minute")});

    vector<growth::statusMap> emptyMap = vector<growth::statusMap>(num_neurons);
    vector<growth::statusMap> mainVec =
        vector<growth::statusMap>(num_neurons, mainMap);

    growth::create_neurons_(mainVec, emptyMap, emptyMap);

    growth::simulate_(growth::Time(200, 0, 0, 0));

    growth::SkelNeurite axon, dendrites, nodes, growth_cones, somas;
    std::vector<growth::stype> gid = {0};
    // growth::get_skeleton(axon, dendrites, nodes, growth_cones, somas, gid);
    growth::get_swc_("myswc", gid, 10);

    return 0;
}
