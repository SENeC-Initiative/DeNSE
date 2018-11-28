#!/usr/bin/env python
#-*- coding:utf-8 -*-

from collections import namedtuple

''' Generate the complete models_manager.cpp '''


in_file  = "models_manager.cpp.in"
out_file = "models_manager.cpp"
mc_args  = ("method", "filename", "classname") 


class ModelComponent(namedtuple("component", mc_args)):

    def __new__(cls, *args, **kwargs):
        method = ""
        if args:
            method = args[0]
        method = kwargs.get("method", method)

        assert "_" not in method, "`method` name must not contain underscores."

        obj = super(ModelComponent, cls).__new__(cls, *args, **kwargs)
        return obj


'''
List of algorithms
------------------
Each of the three lists (elongation, steering, direction selection) contains
all the possible variants for one component.
Each method variant is described by a `ModelComponent` object containing:
- the complete name of the method (`method` argument),
- the name of the header file (.hpp) where the method is declared (`filename`),
- the name of the class, as it appears inside the header file (`classname`).

Note
----
The `method` argument must not contain underscores.
'''

# list of elongation methods
elongation_methods = [
    ModelComponent(method="constant",
                   filename="cst_elongation.hpp",
                   classname="CstElongationModel"),

    ModelComponent(method="gaussian-fluctuations",
                   filename="gfluct_elongation.hpp",
                   classname="GFluctElongationModel"),

    ModelComponent(method="resource-based",
                   filename="resource_based_elongation.hpp",
                   classname="ResourceBasedElongationModel"),
]

# list of steering methods
steering_methods = [
    ModelComponent(method="pull-only",
                   filename="pull_only_steering.hpp",
                   classname="PullOnlySteeringModel"),

    ModelComponent(method="memory-based",
                   filename="memory_based_steering.hpp",
                   classname="MemBasedSteeringModel"),
]

# list of direction selection methods
direction_selection_methods = [
    ModelComponent(method="noisy-maximum",
                   filename="nm_direction_selector.hpp",
                   classname="NMDirectionSelector"),

    ModelComponent(method="noisy-weighted-average",
                   filename="nwa_direction_selector.hpp",
                   classname="NWADirectionSelector"),

    ModelComponent(method="run-and-tumble",
                   filename="rt_direction_selector.hpp",
                   classname="RTDirectionSelector"),
]

# dictionary of abbreviations (complete to abbrev)
abbrev = {
    "constant": "cst",
    "gaussian-fluctuations": "gf",
    "resource-based": "res",
    "pull-only": "po",
    "memory-based": "mem",
    "noisy-maximum": "nm",
    "noisy-weighted-average": "nwa",
    "run-and-tumble": "rt",
}

# special models with a simplified name
specials = {
    "simple-random-walk": "constant_pull-only_noisy-weighted-average",
    "run-and-tumble": "constant_pull-only_run-and-tumble",
    "netmorph-like": "constant_memory-based_noisy-maximum",
    # ~ "neuromac-like": "constant_self-referential-forces_noisy-weighted-average" @todo
}


''' Strings to fill the C++ file '''

include       = '#include "{}"\n'
translate_str = '      {{"{}", "{}"}},\n'
# ~ model_str    = '''
    # ~ kernel().neuron_manager.register_model(
        # ~ "{model_name}",
        # ~ std::dynamic_pointer_cast<GrowthCone>(
            # ~ std::make_shared<GrowthConeModel<{el}, {steer}, {dirsel}>>("{model_name}")));
# ~ '''
model_str    = '''
    models_.insert({{
        "{model_name}",
        GrowthConeModel<{el}, {steer}, {dirsel}>::create_gc_model("{model_name}")
    }});
'''


def generate_models_manager():
    mm_input = ""

    with open(in_file, "r") as f:
        mm_input = f.read()

    default_model            = '"run-and-tumble"'
    model_includes           = "// elongation methods\n"
    elongation_list          = "\n"
    steering_list            = "\n"
    direction_selection_list = "\n"
    abbrev_to_full           = "\n"
    full_to_abbrev           = "\n"
    special_map              = "\n"
    init_str                 = ""

    # fill lists and includes
    for el in elongation_methods:
        elongation_list += "      " + '"' + el.method + '"' + ",\n"
        model_includes  += include.format(el.filename)

    model_includes += "// steering methods\n"

    for steer in steering_methods:
        steering_list  += "      " + '"' + steer.method + '"' + ",\n"
        model_includes += include.format(steer.filename)

    model_includes += "// direction selection methods\n"

    for dirsel in direction_selection_methods:
        direction_selection_list += "      " + '"' + dirsel.method + '"' + ",\n"
        model_includes           += include.format(dirsel.filename)

    # fill abbreviations and specials
    for k, v in abbrev.items():
        abbrev_to_full += translate_str.format(v, k)
        full_to_abbrev += translate_str.format(k, v)

    for k, v in specials.items():
        special_map += translate_str.format(k, v)

    # fill 'init_models' function
    for el in elongation_methods:
        for steer in steering_methods:
            for dirsel in direction_selection_methods:
                # make model name
                model = el.method + "_" + steer.method + "_" + dirsel.method
                # fill string
                init_str += model_str.format(
                    model_name=model,
                    el=el.classname,
                    steer=steer.classname,
                    dirsel=dirsel.classname
                )

    # format file string
    mm_input = mm_input.format(
        model_includes=model_includes,
        elongation_list=elongation_list,
        steering_list=steering_list,
        direction_selection_list=direction_selection_list,
        abbrev_to_full=abbrev_to_full,
        full_to_abbrev=full_to_abbrev,
        specials=special_map,
        default_model=default_model,
        init_models=init_str,
    )

    # write file string
    with open(out_file, "w") as f:
        f.write(mm_input)


''' Make the file '''

generate_models_manager()
