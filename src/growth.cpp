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
    mainMap.insert({"speed_growth_cone", growth::Property(20., "micrometer / minute")});

    vector<growth::statusMap> emptyMap = vector<growth::statusMap>(num_neurons);
    vector<growth::statusMap> mainVec =
        vector<growth::statusMap>(num_neurons, mainMap);

    growth::create_neurons_(mainVec, emptyMap, emptyMap);

    growth::simulate_(growth::Time(200, 0, 0, 0));

    growth::SkelNeurite axon, dendrites, nodes, growth_cones, somas;
    std::vector<size_t> gid = {0};
    // growth::get_skeleton(axon, dendrites, nodes, growth_cones, somas, gid);
    growth::get_swc_("myswc", gid, 10);

    return 0;
}
