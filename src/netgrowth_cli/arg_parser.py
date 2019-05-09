# -*- coding: utf-8 -*-
#
# arg_parser.py
#
# This file is part of DeNSE.
#
# Copyright (C) 2019 SeNEC Initiative
#
# DeNSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# DeNSE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DeNSE. If not, see <http://www.gnu.org/licenses/>.


import argparse
import numpy


def getparser():
        parser = argparse.ArgumentParser(description='test DeNSE simulator,set Simulation, Neurons, Record data, Process, Environment')

        ##Simulation
        parser.add_argument('-S','--simulate', action='store_true',default=False,
                            help='run a the simulation')
        parser.add_argument('-Ss','--step_size', type=int, default=10,
                            help='duration of step in seconds')
        parser.add_argument('-Sn','--n_loops', type=int, default=10,
                            help='loops of dimension step_size')
        parser.add_argument('-Sj','--json_dict',  type=str, default=None,
                                help="path to the json dict to load, json needs: neuron_params, axon_params and dendrite_params")

        ##Neurons and Neurite numbers
        parser.add_argument('-Nn','--neurons', type=int, default=1,
                            help='number of simulated neurons')
        parser.add_argument('-Nd','--dendrites', type=int, default=0,
                            help='number of dendrites for neuron')
        parser.add_argument('-Np','--d_params',  action='store_true',default=False,
                            help='use different params for dendrites')

        ##Data IO
        parser.add_argument('-Rp','--plot',  action='store_true',default=False,
                            help="perform the plot at the end of the simulation")
        parser.add_argument('-Rb','--bt_visualize',  action='store_true',default=False,
                            help="import the simulation with btmorph and visualize neurons")
        parser.add_argument('-Rd','--dyn_analyze',  action='store_true',default=False,
                            help="import the record file and visualize CR and length for each GC")
        parser.add_argument('-Rs','--save',  action='store_true', default=False,
                            help="save the swc file of the simulation")
        parser.add_argument('-Re','--record_enable', action='store_true', default=False,
                            help="record data during the simulation")
        parser.add_argument('-Rr','--swc_resolution', type=int, default=40,
                            help='swc resolution')

        ##multithreading and random
        parser.add_argument('-Pp','--proc', type = int, default =1,
                            help='set the number of used processors')
        parser.add_argument('-Pt','--test_random', action='store_true',default=False,
                            help='run a test over the rnd_gen instead of simulation')
        parser.add_argument('-Ps','--seeds', type=int, default=[numpy.random.randint(1000)], nargs='+',
                            help='seeds for the rundom generator')

        ##environment
        parser.add_argument('-En','--not_environment', action='store_true', default=False,
                            help='switch off the environment')
        parser.add_argument('-Ec','--culture_file', type=str, default=None,
                            help='filepath to the culture file')

        #neuron model
        parser.add_argument('-Mc','--critical_resource_model',type=str, default=None,
                help='set critical resource model: \n"Langevin"\n "Lurd"\n "Gaussian", else use Gaussian elongation')
        parser.add_argument('-Ml','--lateral_on',action='store_true', default=False,
                            help='switch on lateral branching')
        parser.add_argument('-Mv','--vanpelt_on',action='store_true', default=False,
                            help='switch on vanpelt branching')

        return parser.parse_args()
