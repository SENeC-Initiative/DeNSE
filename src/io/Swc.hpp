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
    Swc(std::string output_file, unsigned int resolution);
    void to_swc(const Neuron *, size_t gid);
    void close_file();
    ~Swc();
};
}

#endif /* SWC_H */
