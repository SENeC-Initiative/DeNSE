# !/usr/bin/env python
# -*- coding:utf-8 -*-
# This software is part of the NetGrowth project and the SENEC initiative
#
# The file contains algorithms to compute the following characterization

# cosine_correlation
# msd_1D
# msd_2D
# msd_fft
# contraction Euclidean distance over curvilinear absciss
# tortuosity_local Measures the average variation of the angle

import numpy as np
# import uncertainties as un
# from uncertainties import unumpy


# ============================
#  Correlation futions
# ============================


def tortuosity_local(theta, rho, first=1):
    """
    Measure the local tortuosity as defined in u
    'Computation of tortuosity vessels' 10.1109/ICAPR.2015.7050711

    Params:
    -------
    theta: np.array (N realizations, L size of realization )
           it is the relative angle between 2 subsequent\
           segments
    rho:   np.array (N realizations, L size of realization )
           it is the length of the segment
    first: first element to be processed, skip the previous.

    Returns:
    --------
    np.array (L size of realization, 3 percentile values)
                the percentile distribution with [50%, 75%, 25%]
    """

    # theta =  (arrays[:,:].transpose()- arrays[:,0]).transpose()
    # differential measure reduce size -1
    # print(theta.shape)
    theta = np.abs(theta[:, first:] - theta[:, :-first])
    rho = rho[:, :]
    max_len = rho.shape[1]
    tortuosity_local = np.zeros((max_len, 3))
    length_local = np.zeros((max_len, 3))
    # since distance needs to be greater thean 1, first element is jumped!
    for shift in range(first, max_len):
        somma = np.sum(theta[:, first:shift], axis =1)
        tortuosity_local[shift-first,1:]= np.percentile(somma,q=[75, 25], axis=0, interpolation='midpoint')
        tortuosity_local[shift-first,0] = np.mean(somma)
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT
        somma = np.sum(rho[:, first:shift], axis=1)
        length_local[shift-first,1:] = np.percentile(somma, q=[75, 25], axis=0, interpolation='midpoint')
        length_local[shift-first,0] =np.mean(somma, axis=0)
    length_local[0] = np.array([1, 1, 1])
    tortuosity_local = tortuosity_local/length_local
    tortuosity_local[0] = np.array([1., 1.2, 0.8])
    return tortuosity_local


def contraction(curvilinear, xy, first=1):
    """
    Measure the  ratio between euclidean distance from the origin and \
        curvilinear distance over the path.
        It is always less than 1.

    Params:
    -------
    xy:    np.array (N realizations, L size of realization, 2 [x,y] )
           The position of the origin of each segment
    rho:   np.array (N realizations, L size of realization )
           it is the length of the segment
    first: first element to be processed, skip the previous.

    Returns:
    --------
    np.array (L size of realization, 3 percentile values)
                the percentile distribution with [50%, 75%, 25%]
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

        ratio[shift-first,1:] = np.percentile(fraction, axis = 0,
                                           q=[75, 25], interpolation='midpoint')
        ratio[shift-first,0] = np.mean(fraction,\
                                        axis =0)

        ratio[0] = np.array([1., 1.01, 0.09])
    return ratio


def cosine_correlation(array, first=1):
    """
    Compute the mean value of <cos(origin-value_i)>, where origin is the \
        value of array's first element.

    Params:
    -------
    array: np.array (N realizations, L size of realization )
    first: first element to be processed, skip the previous.

    Returns:
    2D array, the percentile distribution with [50%, 75%, 25%]
    """

    theta = (array[:, :].transpose() - array[:, first]).transpose()
    max_len = theta.shape[1]
    cosine = np.zeros((max_len, 3))
    for shift in range(first, max_len):
        cosine[shift-first,0:] = np.percentile(np.cos(theta[:, shift]), q=[50,60,40], axis=0, interpolation='midpoint')
        # cosine[shift-first,0] = np.mean(np.cos(theta[:, shift]), axis = 0)
        cosine[0] = np.array([1., 1.5, 0.5])
    return cosine


def msd_1D(array, first=1):
    """
    Compute the mean square displacement.\
    Delta_x^2(n) = <(x_n - x_0)^2 >

       Params:
    -------
    array: np.array (N realizations, L size of realization)
           the array to compute the MSD
    first: first element to be processed, skip the previous.

    Returns:
    --------
    2D array, the percentile distribution with [50%, 75%, 25%]
    """

    print (array.shape)
    theta = (array[:, :].transpose() - array[:, first]).transpose()
    max_len = theta.shape[1]
    msd = np.zeros((max_len, 3))
    for shift in range(first, max_len):
        theta_sum = theta[:, shift]**2
        # import pdb; pdb.set_trace()  # XXX BREAKPOINT

        msd[shift-first,:] = np.percentile(theta_sum, q=[50,75,25], axis=0,\
                 interpolation='midpoint')
        # msd[shift-first,0] = np.mean(theta_sum, axis=0)
    msd[0][1] = 0.01
    msd[0][2] = 0
    return msd


def msd_2D(xy, first=1):
    """
    Compute the mean square displacement of 2 dimensional vector.\
    Delta^2(n) = <(z_n - z_0)^2 >
    where z = \sqrt(x^2+y^2)

    Params:
    -------
    2Darray: np.array (N realizations, L size of realization, 2 [x,y])
             the array to compute the MSD
    first:   first element to be processed, skip the previous.

    Returns:
    --------
    numpy 2D array, the percentile distribution with [50%, 75%, 25%]
    """

    dx = (xy[:, :, 0].transpose() - xy[:, first, 0]).transpose()
    dy = (xy[:, :, 1].transpose() - xy[:, first, 1]).transpose()
    max_len = dx.shape[1]
    msd = np.zeros((max_len, 3))
    for shift in range(first, max_len):
        delta = dx[:, shift]**2+dy[:, shift]** 2
        msd[shift-first,1:] = np.percentile(delta, q=[75, 25], axis=0, interpolation='midpoint')
        msd[shift-first,0] = np.mean(delta, axis = 0)
    return msd


# def msd_fft(r):
    # """
    # Compute the mean square displacement for the sequence
    # """
    # def autocorrFFT(x):
        # N = len(x)
        # F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
        # PSD = F * F.conjugate()
        # res = np.fft.ifft(PSD)
        # res = (res[:N]).real  # now we have the autocorrelation in convention B
        # n = N*np.ones(N)-np.arange(0, N)  # divide res(m) by (N-m)
        # return res/n  # this is the autocorrelation in convention A
    # N = len(r)
    # if len(r.shape) > 1:
        # D = np.square(r).sum(axis=1)
        # D = np.append(D, 0)
        # S2 = sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    # else:
        # D = np.square(r)
        # D = np.append(D, 0)
        # S2 = autocorrFFT(r)
    # Q = 2*D.sum()
    # S1 = np.zeros(N)
    # for m in range(N):
        # Q = Q-D[m-1]-D[N-m]
        # S1[m] = Q/(N-m)
    # return S1-2*S2
