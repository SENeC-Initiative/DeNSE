#!/usr/bin/python3

import dense as ds
import numpy as np
from arg_parser import getparser
import os
import json
# import datetime

neuron_params = {
    # "growth_cone_model": "persistent_random_walk",
    # "growth_cone_model": "self_referential_forces",
    "growth_cone_model": "run_tumble",
    "axon_angle": 90.,

    "persistence_length": 99000050.,

    "srf_inertial_force": 0.9,
    "sfr_inertial_decay": 2.,
    "sfr_somatropic_force": 0.5,
    "sfr_somatropic_decay": .5,

    "axon_diameter" : 1.,
    "persistence_length": -1.,
    "angular_resolution": 30.0,
    "speed_growth_cone": 10.,
    "sensing_angle": 0.1195,

    "filopodia_wall_affinity": 1.5,
    "filopodia_finger_length": 20.0,
    "use_lateral_branching": False,
    "uniform_branching_rate": 0.0008,


    "use_van_pelt": True,
    "B": 4.,
    "E": 1.0,
    "S": 0.001,
    "T": 50.,


    "gc_split_angle_mean": 30.3,

    "use_critical_resource": False,
    "res_speed_factor": 0.10,
    "res_amount": 1.,
    "res_leakage": 0.05,
    "res_retraction_threshold": 0.30,
    "res_elongation_threshold": 0.50,
    "res_split_th": 0.80,
    "res_demand_correlation": 0.9910,
    "res_demand_stddev": 0.2,
    "res_demand_mean": 1.,
    "res_use_ratio": 0.7
}

dendrite_params = {
    "growth_cone_model": "run_tumble",
    "persistence_length": 10.0,
    "angular_resolution": 30.0,
    "speed_growth_cone": 5.,
    "sensing_angle": 0.1195,

    "filopodia_wall_affinity": 2.0,
    "filopodia_finger_length": 50.0,
    "use_lateral_branching": False,
    "uniform_branching_rate": 0.0002,

    "use_van_pelt": True,
    "B": 10.0,
    "E": 0.905,
    "S": 1.0,
    "T": 0.001,

    "gc_split_angle_mean": 10.3,

    "res_speed_factor": 0.1,
    "res_amount": 1.,
    "res_leakage": 0.05,
    "res_retraction_threshold": 0.30,
    "res_elongation_threshold": 0.50,
    "res_split_th": 0.80,
    "res_demand_correlation": 0.9910,
    "res_demand_stddev": 0.2,
    "res_demand_mean": 1.,
    "res_use_ratio": 0.7
}


def step(n, loop_n, plot=True):
    ds.simulate(n)
    if plot:
        ds.plot_neurons(
            show_nodes=True, save_path="examples/last_run_5")


if '__main__' is __name__:

    args = getparser()

    if args.json_dict:
        with open(args.json_dict, "r") as params:
            ext_params = json.load(params)
        axon_params = ext_params["axon"]
        dendrite_params = ext_params["dendrite"]
        neuron_params = ext_params["neuron"]
    else:
        axon_params = neuron_params
        dendrite_params = dendrite_params

    np.random.seed(args.seeds[0])

    if args.critical_resource_model is None:
        neuron_params["use_critical_resource"] = False
    else:
        print(args.critical_resource_model)
        neuron_params['growth_cone_model'] = \
            neuron_params['growth_cone_model'] + args.critical_resource_model
        neuron_params["use_critical_resource"] = True
        dendrite_params["use_critical_resource"] = True

    if args.vanpelt_on:
        neuron_params["use_van_pelt"] = True
    else:
        neuron_params["use_van_pelt"] = False

    if args.lateral_on:
        neuron_params["use_lateral_branching"] = True
    else:
        neuron_params["use_lateral_branching"] = False

    kernel_dict = {"num_local_threads": args.proc, "seeds": args.seeds}
    simulation_ID = ds.generate_simulation_id(
        kernel_dict, neuron_params, dendrite_params)
    kernel_dict["record_enabled"] = args.record_enable
    if args.record_enable and not os.path.exists(simulation_ID):
        os.makedirs(simulation_ID)

    # record_file_=simulation_ID+"/record.dat"
    # swc_file_=simulation_ID+"/morphology.swc"
    # json_file=simulation_ID+"/info.json"

    gids = None
    if not args.d_params:
        dendrite_params = axon_params

    if not args.not_environment:
        ds.get_kernel_status(kernel_dict, simulation_ID)
        if args.culture_file is not None:
            culture = ds.CreateEnvironment(
                args.culture_file, min_x=0, max_x=1000)
        else:
            culture = ds.CreateEnvironment(
                os.getcwd() + '/culture_from_filled_polygons.svg', min_x=0,
                max_x=100)
        gids = ds.create_neurons(culture=culture,
                                       n=args.neurons,
                                       params=neuron_params,
                                       dendrites_params=dendrite_params,
                                       axon_params=axon_params,
                                       num_neurites=1 + args.dendrites)
    else:
        kernel_dict["environment_required"] = False
        kernel_dict["resolution"] = 1.

        ds.get_kernel_status(kernel_dict, simulation_ID)
        neuron_params["position"] = np.random.uniform(
            -500, 500, (args.neurons, 2))
        gids = ds.create_neurons(n=args.neurons,
                                       params=neuron_params,
                                       dendrites_params=dendrite_params,
                                       axon_params=axon_params,
                                       num_neurites=1 + args.dendrites)

    # ds.plot_environment(culture)

    if args.test_random:
        vals = ds.TestRandomGen(1000000)
        vals = np.array(vals)
        np.savetxt('random_test.txt', vals)

    if args.simulate:
        step(100, 1, False)
        neuron_params["use_van_pelt"] = True
        ds.set_object_properties(gids,
                            axon_params = neuron_params,
                            dendrites_params= neuron_params,
                            params = neuron_params)
        for loop_n in range(args.n_loops):
            step(args.step_size, loop_n, args.plot)

        kernel_ID = ds.get_simulation_id()

        if args.save:
            save_path = simulation_ID
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            # save to default path: "simulation_ID/..."
            swc_file = ds.SaveSwc(
                filepath=save_path, swc_resolution=args.swc_resolution)
            json_file = ds.save_json_info(filepath=save_path)
        if args.bt_visualize:
            ds.BtmorphVisualize(save_path)
        if args.dyn_analyze:
            ds.GrowthConeDynamicsAnalyzer()
