# -*- coding:utf-8 -*-
# This software is part of the DeNSE project and the SENEC initiative
#
# The file contains algorithms to compute the following characterization

# cosine_correlation
# msd_1D
# msd_2D
# msd_fft
# contraction Euclidean distance over curvilinear absciss
# tortuosity_local Measures the average variation of the angle

import numpy as np
from scipy.special import digamma

from .. import _pygrowth as _pg


# ============================
#  Correlation futions
# ============================

def local_tortuosity(theta, rho, first=1):
    """
    Measure the local tortuosity as defined in u
    'Computation of tortuosity vessels' 10.1109/ICAPR.2015.7050711

    Parameters
    ----------
    theta: array of size N
        Relative angle between subsequent segments
    rho:   array of size N
        Length of each segment
    first: int, optional (default: 1)
        First element to process, skip previous.

    Returns
    -------
    tortuosity : array of size N - first
    """
    theta   = np.abs(theta[first:] - theta[:-first])
    rho     = rho[:, :]
    max_len = rho.shape[1]

    tortuosity_local = np.zeros((max_len, len(percentiles)))
    length_local = np.zeros((max_len, len(percentiles)))

    # since distance needs to be greater thean 1, first elements are skipped
    for shift in range(first, max_len):
        somma = np.sum(theta[:, first:shift], axis =1)
        tortuosity_local[shift-first,1:]= np.percentile(
            somma, q=percentiles, axis=0, interpolation='midpoint')
        tortuosity_local[shift-first,0] = np.mean(somma)
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT
        somma = np.sum(rho[:, first:shift], axis=1)
        length_local[shift-first,1:] = np.percentile(somma, q=[75, 25], axis=0, interpolation='midpoint')
        length_local[shift-first,0] =np.mean(somma, axis=0)

    length_local[0]     = np.array([1, 1, 1])
    tortuosity_local    = tortuosity_local / length_local
    tortuosity_local[0] = np.array([1., 1.2, 0.8])
    return tortuosity_local


def tortuosity_local(theta, rho, percentiles=(50, 75, 25), first=1):
    """
    Measure the local tortuosity as defined in u
    'Computation of tortuosity vessels' 10.1109/ICAPR.2015.7050711

    Parameters
    ----------
    theta: np.array (N realizations, L size of realization )
        Relative angle between subsequent segments
    rho:   np.array (N realizations, L size of realization )
        Length of each segment
    first: int, optional (default: 1)
        First element to process, skip previous.

    Returns
    -------
    np.array (L size of realization, 3 percentile values)
        Percentile distribution with [50%, 75%, 25%]
    """
    theta   = np.abs(theta[:, first:] - theta[:, :-first])
    rho     = rho[:, :]
    max_len = rho.shape[1]

    tortuosity_local = np.zeros((max_len, len(percentiles)))
    length_local = np.zeros((max_len, len(percentiles)))

    # since distance needs to be greater thean 1, first elements are skipped
    for shift in range(first, max_len):
        somma = np.sum(theta[:, first:shift], axis =1)
        tortuosity_local[shift-first,1:]= np.percentile(
            somma, q=percentiles, axis=0, interpolation='midpoint')
        tortuosity_local[shift-first,0] = np.mean(somma)
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT
        somma = np.sum(rho[:, first:shift], axis=1)
        length_local[shift-first,1:] = np.percentile(somma, q=[75, 25], axis=0, interpolation='midpoint')
        length_local[shift-first,0] =np.mean(somma, axis=0)

    length_local[0]     = np.array([1, 1, 1])
    tortuosity_local    = tortuosity_local / length_local
    tortuosity_local[0] = np.array([1., 1.2, 0.8])
    return tortuosity_local


def contraction(curvilinear, xy, first=1):
    """
    Measure the  ratio between euclidean distance from the origin and
    curvilinear distance over the path. It is always smaller than 1.

    Parameters
    ----------
    xy:    np.array (N realizations, L size of realization, 2 [x,y] )
        Position of the origin of each segment
    rho:   np.array (N realizations, L size of realization )
        Length of each segment
    first: int, optional (default: 1)
        First element to process, skip previous.

    Returns
    -------
    np.array (L size of realization, 3 percentile values)
        Percentile distribution with [50%, 75%, 25%]
    """
    dx = (xy[:, :, 0].transpose() - xy[:, first, 0]).transpose()
    dy = (xy[:, :, 1].transpose() - xy[:, first, 1]).transpose()

    max_len = dx.shape[1]
    ratio = np.zeros((max_len, 3))

    for shift in range(first, max_len):
        r = np.sqrt((dx[:, shift])**2 + (dy[:, shift])**2)
        length = np.sum(np.abs(curvilinear[:, first:shift]), axis=1)
        fraction = r/length
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT

        ratio[shift-first,1:] = np.percentile(
            fraction, axis=0, q=[75, 25], interpolation='midpoint')
        ratio[shift-first,0] = np.mean(fraction, axis=0)

        ratio[0] = np.array([1., 1.01, 0.09])
    return ratio


def cosine_correlation(array, first=1):
    """
    Compute the mean value of <cos(origin-value_i)>, where origin is the
    value of `array`'s first element.

    Parameters
    ----------
    array: np.array (N realizations, L size of realization )
    first: first element to be processed, skip the previous.

    Returns
    -------
    2D array, the percentile distribution with [50%, 75%, 25%]
    """

    theta = (array[:, :].transpose() - array[:, first]).transpose()
    max_len = theta.shape[1]
    cosine = np.zeros((max_len, 3))
    for shift in range(first, max_len):
        cosine[shift-first, 0:] = np.percentile(
            np.cos(theta[:, shift]), q=[50,60,40], axis=0,
            interpolation='midpoint')
        # cosine[shift-first,0] = np.mean(np.cos(theta[:, shift]), axis = 0)
        cosine[0] = np.array([1., 1.5, 0.5])
    return cosine


def msd_1D(array, first=1):
    """
    Compute the mean square displacement.
    Delta_x^2(n) = <(x_n - x_0)^2 >

    Parameters
    ----------
    array: np.array (N realizations, L size of realization)
           the array to compute the MSD
    first: first element to be processed, skip the previous.

    Returns
    -------
    2D array, the percentile distribution with [50%, 75%, 25%]
    """
    theta = (array[:, :].transpose() - array[:, first]).transpose()
    max_len = theta.shape[1]
    msd = np.zeros((max_len, 3))
    for shift in range(first, max_len):
        theta_sum = theta[:, shift]**2
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT
        msd[shift-first, :] = np.percentile(
            theta_sum, q=[50,75,25], axis=0, interpolation='midpoint')
        # msd[shift-first,0] = np.mean(theta_sum, axis=0)
    msd[0][1] = 0.01
    msd[0][2] = 0
    return msd


def msd_2D(xy, first=1):
    """
    Compute the mean square displacement of 2 dimensional vector.\
    Delta^2(n) = <(z_n - z_0)^2 >
    where z = \sqrt(x^2+y^2)

    Parameters
    ----------
    xy: np.array (N realizations, L size of realization, 2 [x,y])
             the array to compute the MSD
    first:   first element to be processed, skip the previous.

    Returns
    -------
    numpy 2D array, the percentile distribution with [50%, 75%, 25%]
    """

    dx = (xy[:, :, 0].transpose() - xy[:, first, 0]).transpose()
    dy = (xy[:, :, 1].transpose() - xy[:, first, 1]).transpose()
    max_len = dx.shape[1]
    msd = np.zeros((max_len, 3))
    for shift in range(first, max_len):
        delta = dx[:, shift]**2+dy[:, shift]** 2
        msd[shift-first, 1:] = np.percentile(
            delta, q=[75, 25], axis=0, interpolation='midpoint')
        msd[shift-first, 0] = np.mean(delta, axis=0)
    return msd


def tree_asymmetry(neurite, neuron=None):
    '''
    Return the tree-assymetry of a neurite.

    Parameters
    ----------
    neurite : Neurite object or str
        Neurite to analyze.
    neuron : int, optional (default: None)
        GID of the neuron if the neurite was passed as a str and not as an
        object.
    '''
    import neurom as nm

    try:
        tree = neurite.get_tree().neurom_tree()
    except AttributeError:
        nrn = _pg.get_neurons(neuron)[0]
        nrt = nrn.get_neurite(neurite)
        tree = neurite.get_tree().neurom_tree()
    asyms = nm.get("partition_asymmetry", tree)
    return np.average(asyms) / _max_asym(len(neurite.get_tree().tips))


def _max_asym(n):
    return 1. - (digamma(n) - digamma(1.))/(n-1.)
