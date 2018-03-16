#!/usr/bin/python3

import NetGrowth
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

    "rt_persistence_length": 99000050.,

    "srf_inertial_force": 0.9,
    "sfr_inertial_decay": 2.,
    "sfr_somatropic_force": 0.5,
    "sfr_somatropic_decay": .5,

    "axon_diameter" : 1.,
    "rw_persistence_length": -1.,
    "rw_memory_tau": -1.,
    "rw_delta_corr": 1.8,
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
    "CR_speed_factor": 0.10,
    "CR_amount": 1.,
    "CR_leakage": 0.05,
    "CR_retraction_th": 0.30,
    "CR_elongation_th": 0.50,
    "CR_split_th": 0.80,
    "CR_demand_correlation": 0.9910,
    "CR_demand_stddev": 0.2,
    "CR_demand_mean": 1.,
    "CR_use_ratio": 0.7
}

dendrite_params = {
    "growth_cone_model": "run_tumble",
    "rt_persistence_length": 10.0,
    "rw_memory_tau": 200.,
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

    "CR_speed_factor": 0.1,
    "CR_amount": 1.,
    "CR_leakage": 0.05,
    "CR_retraction_th": 0.30,
    "CR_elongation_th": 0.50,
    "CR_split_th": 0.80,
    "CR_demand_correlation": 0.9910,
    "CR_demand_stddev": 0.2,
    "CR_demand_mean": 1.,
    "CR_use_ratio": 0.7
}


def step(n, loop_n, plot=True):
    NetGrowth.Simulate(n)
    if plot:
        NetGrowth.PlotNeuron(
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
    simulation_ID = NetGrowth.GenerateSimulationID(
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
        NetGrowth.SetKernelStatus(kernel_dict, simulation_ID)
        if args.culture_file is not None:
            culture = NetGrowth.CreateEnvironment(
                args.culture_file, min_x=0, max_x=1000)
        else:
            culture = NetGrowth.CreateEnvironment(
                os.getcwd() + '/culture_from_filled_polygons.svg', min_x=0,
                max_x=100)
        gids = NetGrowth.CreateNeurons(culture=culture,
                                       n=args.neurons,
                                       params=neuron_params,
                                       dendrites_params=dendrite_params,
                                       axon_params=axon_params,
                                       num_neurites=1 + args.dendrites)
    else:
        kernel_dict["environment_required"] = False
        kernel_dict["resolution"] = 1.

        NetGrowth.SetKernelStatus(kernel_dict, simulation_ID)
        neuron_params["position"] = np.random.uniform(
            -500, 500, (args.neurons, 2))
        gids = NetGrowth.CreateNeurons(n=args.neurons,
                                       params=neuron_params,
                                       dendrites_params=dendrite_params,
                                       axon_params=axon_params,
                                       num_neurites=1 + args.dendrites)

    # NetGrowth.PlotEnvironment(culture)

    if args.test_random:
        vals = NetGrowth.TestRandomGen(1000000)
        vals = np.array(vals)
        np.savetxt('random_test.txt', vals)

    if args.simulate:
        step(100, 1, False)
        neuron_params["use_van_pelt"] = True
        NetGrowth.SetStatus(gids,
                            axon_params = neuron_params,
                            dendrites_params= neuron_params,
                            params = neuron_params)
        for loop_n in range(args.n_loops):
            step(args.step_size, loop_n, args.plot)

        kernel_ID = NetGrowth.GetSimulationID()

        if args.save:
            save_path = simulation_ID
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            # save to default path: "simulation_ID/..."
            swc_file = NetGrowth.SaveSwc(
                filepath=save_path, swc_resolution=args.swc_resolution)
            json_file = NetGrowth.SaveJson(filepath=save_path)
        if args.bt_visualize:
            NetGrowth.BtmorphVisualize(save_path)
        if args.dyn_analyze:
            NetGrowth.GrowthConeDynamicsAnalyzer()
