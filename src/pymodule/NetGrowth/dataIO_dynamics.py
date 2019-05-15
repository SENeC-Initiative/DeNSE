#!/usr/bin/env cython
#-*- coding:utf-8 -*-
#cython: boundscheck=False, wraparound=False, initializedcheck=False
#cython: cdivision=True, embedsignature=True

from . import _pygrowth as _pg
from .dataIO import ImportRecordFile
import numpy as np
import matplotlib.pyplot as plt
import os


def GrowthConeDynamicsAnalyzer(record_file_="default"):
    if record_file_=="default":
        record_file_ = os.path.join(_pg.GetSimulationID(),"/record.dat")
    events,steps = ImportRecordFile(record_file_)
    gc_list = _step_data_array(steps)
    plot_dynamic_data(gc_list)

def _step_data_array(steps):
    """
    Assign each growth cone a matrix. and store them in a list.

    Parameters:
    ----------

    steps : list
    the list with 'step' information retrieved in ImportRecordFile

    Return:
    ----------
    gc_list, list of np.array with data.

    """
    steps=np.array(steps, dtype=np.float)
    gc_columns = int (np.max(steps[:,1]))
    gc_list=[]
    for gc in range(1, gc_columns+1):
        gc_list.append(steps[np.where(steps[:,1]==gc)])
    return gc_list


def plot_dynamic_data(gc_list):
    """
    Plot the data on growth cone dynamics, each GC with different color

    Parameters:
    gc_list: list of np.array
    """
    fig, ((ax1, ax2),(ax3, ax4))  = plt.subplots(2,2,sharex=True)
    ax3.set_title("Critical_resource received")
    ax2.set_title("Distance from soma over time")
    ax1.set_title("Critical resource used")
    ax4.set_title("Critical resource left")
    ax1.set_xlabel("time (step)")
    ax3.set_ylabel("CR Demand")
    ax2.set_ylabel("CR Used")
    ax1.set_ylabel("CR Left")
    ax2.set_ylabel("distance from soma ('um')")
    for gc in gc_list:
        ax1.plot(gc[:,0],gc[:,4], ls='-')
        ax2.plot(gc[:,0],gc[:,2], ls='-')
        ax3.plot(gc[:,0],gc[:,3], ls='-')
        ax4.plot(gc[:,0],gc[:,5], ls='-')

    # ax1.plot([0], ls='-', label="CR received",c='k')
    # for gc in gc_list:
        # ax1.plot(gc[:,0],gc[:,4], ls='--')
    # ax1.plot([0], ls='--', label="CR used", c='k')
    # ax1.set_prop_cycle(None)
    # for gc in gc_list:
        # ax1.plot(gc[:,0],gc[:,5], ls=':')
    # ax1.plot([0], ls=':', label="CR left", c='k')
    # ax1.set_prop_cycle(None)
    # ax1.legend(loc='upper center', shadow=True)

    plt.show()

