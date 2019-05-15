#include "models_manager.hpp"

// Include from kernel
#include "kernel_manager.hpp"

// Include from libgrowth
#include "elements_types.hpp"

// Include from elements
#include "GrowthCone.hpp"

// Include from models
#include "gc_dir_el.cpp"
#include "gc_random_walk.hpp"


namespace growth
{

void init_models()
{
    kernel().neuron_manager.register_model(
        "default",
        std::dynamic_pointer_cast<GrowthCone>(std::make_shared<GrowthCone>()));

    kernel().neuron_manager.register_model(
        "random_walk_Langevin",
        std::dynamic_pointer_cast<GrowthCone>(
            std::make_shared<GrowthCone_Elongation_Direction<
                GrowthCone_Critical_Langevin, GrowthCone_RandomWalk>>()));

    kernel().neuron_manager.register_model(
        "random_walk_Gaussian",
        std::dynamic_pointer_cast<GrowthCone>(
            std::make_shared<GrowthCone_Elongation_Direction<
                GrowthCone_Critical_Gaussian, GrowthCone_RandomWalk>>()));

    kernel().neuron_manager.register_model(
        "random_walk_Lurd",
        std::dynamic_pointer_cast<GrowthCone>(
            std::make_shared<GrowthCone_Elongation_Direction<
                GrowthCone_Critical_Lurd, GrowthCone_RandomWalk>>()));

    kernel().neuron_manager.register_model(
        "random_walk", std::dynamic_pointer_cast<GrowthCone>(
                           std::make_shared<GrowthCone_RandomWalk>()));

#ifndef NDEBUG
    for (auto member : kernel().neuron_manager.model_map_)
    {
        printf(" model_map key %s\n", member.first.c_str());
    }
#endif
}
}
