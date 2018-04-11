# !/usr/bin/env python
# -*- coding:utf-8 -*-
#
# This software is part of the NetGrowth project and the SENEC initiative

import numpy as np
from scipy import optimize

from .. import _pygrowth as _pg
from .algorithms import msd_1D, cosine_correlation, msd_2D, tortuosity_local, contraction
from .containers import Population


class EnsembleRW(Population):
    """
    Store all the neurons in an unique object.
    Create the ensemble adding a population to it.
    Once a population is ready perform analysis on it.
    """

    def __init__(self, population, name=None):
        # Population.__init__(self, iption)
        # Single neurons properties
        self.population = population
        self.name = self.population.name if name is None else name

        self.theta = None
        self.r = None
        self.xy = None
        self.tortuosity_local = None
        self.msd_1D = None
        self.cosine = None
        self.msd_2D = None
        self.contraction = None
        self.effective_length = None

        self.swc_resolution = None
        self.kernel_resolution = _pg.GetKernelStatus("resolution")

    def set_swc_resolution(self, swc_resolution):
        self.swc_resolution = swc_resolution

    def set_kernel_resolution(self, kernel_resolution):
        self.kernel_resolution = kernel_resolution

    def normalize_paths(self, neurite_type):
        # measure the shortest path to make all path equal length
        shapes=[]
        for neuron in self.population.neurons:
            if neurite_type is "axon":
                for branch in neuron.axon.branches:
                    shapes.append(branch.r.shape[0])
            if neurite_type is "dendrite":
                for dendrite in neuron.dendrites:
                    for branch in dendrite.branches:
                        shapes.append(branch.r.shape[0])
        max_len = np.min(shapes)
        # for neuron in self.neurons:
            # if neuron.axon.xy.shape[1] < min_shap:
                # min_shap = neuron.axon.xy.shape[1]
        # min_shap = min_shap-2
        self.theta = []
        self.r     = []
        self.xy    = []
        for neuron in self.population.neurons:
            if neurite_type is "axon":
                for branch in neuron.axon.branches:
                    self.theta.append(branch.theta[:max_len])
                    self.r.append(branch.r[:max_len])
                    self.xy.append(branch.xy[:max_len])

            if neurite_type is "dendrite":
                for dendrite in neuron.dendrites:
                    for branch in dendrite.branches:
                        self.theta.append(branch.theta[:max_len])
                        self.r.append(branch.r[:max_len])
                        self.xy.append(branch.xy[:max_len])
        self.theta = np.array(self.theta )
        self.r     = np.array(self.r     )
        self.xy    = np.array(self.xy    )
        assert(self.theta.shape == self.r.shape)
        assert(self.xy.shape[0] == self.r.shape[0])


    def characterizeRW(self, neurite_type, first=1):
        """
        Run characterization algorithm on the random_walk ensemble:
        contraction  : measure the tortuosity length over distance ratio between curvilinear distance and euclidean
        tortuosity_local : measure the correlation between successive delta.
        msd          : measure the mean square displacement of angles
        space_msd    : measure the mean square displacement of cartesian coordinates.

        all the measure are performed with a single istance of random walk, the average is perfermed
        over block with a distance longer than the correlation length.

        Params:
        ------
        max_len: maximal relative distance to measure for the different algorithms.
        decorrelation_length: distance of correlation to resize the path in uncorrelated segments.
        """

        if first < 1:
            raise Exception(
                "first value has to be greater than zero or division by zero occurs")
        self.normalize_paths(neurite_type)
        self.ensemble_shapes = self.theta.shape
        # compute some statistical measures
        # self.tortuosity_loca = tortuosity_local(self.theta, self.r, first)
        # self.tortuosity_local= self.tortuosity_loca[:-10,:];
        self.msd_1D = msd_1D(self.theta, first)
        self.cosine = cosine_correlation(self.theta, first)
        self.msd_2D = msd_2D(self.xy, first)
        self.contraction = contraction(self.r, self.xy, first)


        # import matplotlib.pyplot as plt

        # print(self.msd_1D.shape)
        # print(self.msd_2D.shape)
        # print(self.tortuosity_local.shape)
        # print(self.cosine.shape)

        self.effective_length = np.zeros((self.ensemble_shapes[1]))
        self.effective_length[:]
        length = np.mean(self.r, axis=0)
        for shift in range(first, self.ensemble_shapes[1]):
            self.effective_length[shift-first] = np.sum(length[first:shift])

        # import pdb; pdb.set_trace()  # XXX BREAKPOINT

        self.msd_1D             = self.msd_1D[:-10,:];
        self.cosine             = self.cosine[:-10,:];
        self.msd_2D             = self.msd_2D[:-10,:];
        self.contraction        = self.contraction[:-10,:];
        self.effective_length   = self.effective_length[:-10];
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT


    def fit(self, first=0):
        """
        Fit the averaged quantities
        """
        def constant(x, k):
            return k

        def linear(x, m, b):
            return x*m+b

        def quadratic(x, a, b, c):
            return x**2*a+c

        def exponential(x, tau, a, b):
            return a*np.exp(-x/tau)+b

        def get_tau(array):
            array - 0.367879441
            for n, x in enumerate(array):
                if x < 0:
                    return n

        def get_transient(theta):
            for n, x in enumerate(theta):
                if x > 1.5:
                    break
            return n
            # popt, pcov = optimize.curve_fit(exp_low, xdata, ydata)
            # if _plt is None:
            # _plt = plt
            # _plt.plot(xdata, exp_low(xdata,*popt))
            # _plt.plot(xdata, ydata)
            # return popt[0], popt[1]

        self.fits = {}
        # print (self.msd_1D[:,0])
        theta_max = get_transient(self.msd_1D[:, 0])
        self.fits['msd_1D_ramp'] = optimize.curve_fit(linear, self.effective_length[:theta_max],
                                                    self.msd_1D[:theta_max, 0])
                                                    # bounds=([0.0001,0], [1., 1. ]))
                                                    # sigma=np.abs(self.msd_1D[:theta_max, 1]-self.msd_1D[:theta_max, 2]))



        self.fits['cosine'] = optimize.curve_fit(exponential, self.effective_length[:],
                                                self.cosine[:, 0], bounds=([2, -1, -1], [np.inf, 1, 1]))

        self.fits['msd_2D_quad'] = optimize.curve_fit(quadratic, self.effective_length,
                                                    self.msd_2D[:, 0])
        self.fits['msd_2D_lin'] = optimize.curve_fit(linear, self.effective_length,
                                                    self.msd_2D[:, 0])

        self.fits['contraction'] = optimize.curve_fit(exponential, self.effective_length[:],
                                                      self.contraction[:,0], check_finite=True)

        # not_Nan = np.where(self.tortuosity_local is not np.nan)
        # self.tortuosity_local=self.tortuosity_local[not_Nan]
        # self.fits['tortuosity_local'] = optimize.curve_fit(constant, self.effective_length,self.tortuosity_local[:,0])

        # if np.sum(np.diag(self.fits['msd_2D_quad'][1])**2) > np.sum(np.diag(self.fits['msd_2D_lin'][1])**2):
            # self.msd_2D_fit_wrong = 'msd_2D_quad'
            # self.fits.pop(self.msd_2D_fit_wrong)
            # self.fits['msd_2D'] = self.fits['msd_2D_lin']
            # self.fits.pop('msd_2D_lin')
        # else:
            # self.msd_2D_fit_wrong = 'msd_2D_lin'
            # self.fits.pop(self.msd_2D_fit_wrong)
            # self.fits['msd_2D'] = self.fits['msd_2D_quad']
            # self.fits.pop('msd_2D_quad')
        self.results = {}
        DEBUG_FIT = False
        if DEBUG_FIT:
            import matplotlib.pyplot as plt
            print(self.fits['cosine'])
            plt.plot(self.effective_length[:], self.cosine[:, 0])
            plt.plot(self.effective_length[:], exponential(
                self.effective_length[:], *self.fits['cosine'][0]))
            plt.show()
        for key in self.fits:
            self.results[key] = {"values": {"a"+str(n): x for n, x in enumerate(self.fits[key][0])},
                                 "errors": {"a"+str(n): x for n, x in enumerate(np.diag(self.fits[key][1]))}}
        return self.results
