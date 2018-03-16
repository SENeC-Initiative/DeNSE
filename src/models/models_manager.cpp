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
#include "gc_self_referential_forces.hpp"
#include "gc_run_tumble.hpp"


namespace growth
{

void init_models()
{
    // create all model neurons
    kernel().neuron_manager.register_model(
        "default", std::dynamic_pointer_cast<GrowthCone>(
                           std::make_shared<GrowthCone_RunTumble>()));

    kernel().neuron_manager.register_model(
        "simple_random_walk",
        std::dynamic_pointer_cast<GrowthCone>(std::make_shared<GrowthCone>("simple_random_walk")));

    kernel().neuron_manager.register_model(
        "persistent_random_walk", std::dynamic_pointer_cast<GrowthCone>(
                           std::make_shared<GrowthCone_RandomWalk>()));

    kernel().neuron_manager.register_model(
        "self_referential_forces", std::dynamic_pointer_cast<GrowthCone>(
                           std::make_shared<GrowthCone_SelfReferentialForces>()));

    kernel().neuron_manager.register_model(
        "run_tumble", std::dynamic_pointer_cast<GrowthCone>(
                           std::make_shared<GrowthCone_RunTumble>()));

#ifndef NDEBUG
    for (auto member : kernel().neuron_manager.model_map_)
    {
        printf(" model_map key %s\n", member.first.c_str());
    }
#endif
}
}
