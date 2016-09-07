#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
CONSTANd
Normalizes the data matrix <data> by raking the m by n matrix such that
the row mean and column mean equals to 1/n. Missing information needs to
be presented as NaN values and not as zero values because CONSTANd
employs the Matlab functionality 'nanmean' that is able to ignore
NaN-values when calculating the mean. The variable <nbIterations> is an
integer value that denotes the number of raking cycles. The variable <h>
defines the stopping criteria based on the L1-norm as defined by
Friedrich Pukelsheim, Bruno Simeone in "On the Iterative Proportional
Fitting Procedure: Structure of Accumulation Points and L1-Error
Analysis"

© Dirk Valkenborg & Jef Hooyberghs, 2014
Ported to Python by Joris Van Houtven, 2016
"""

import numpy as np
from numpy import nan as NaN


def constand(data, accuracy, maxIterations):
    assert isinstance(data, np.matrix)
    assert data.dtype is np.dtype('float64')
    assert accuracy > 0
    assert maxIterations > 0

    # initialize variables
    Nrows, Ncols = data.shape
    convergenceTrail = [NaN]*(2*maxIterations)
    convergence = np.inf
    dataT = data.T
    normalizedData = data

    i = 0 # number of current iteration
    # without reshaping the ndarrays, they have shape (x,) (no second value) and the procedure fails.
    Rd = np.ones((Nrows,1)) # R matrix diagonal from RAS
    Sd = np.ones((Ncols,1)) # S matrix diagonal from RAS
    # main loop; iterates until convergence is reached (i.e., L1-norm below variable <h>) or the maximum number of
    # iteration cycles is surpassed.
    while convergence > accuracy and i < maxIterations:

        # fit the rows
        Rd = Rd * 1/Ncols * np.asarray(1/np.nanmean(normalizedData,1))
        normalizedData = (dataT * Rd).T *Sd
        # calculate deviation from column marginals; row deviation is zero at even indices. (index start = 0)
        convergenceTrail[2*i] = Nrows * 0.5 * np.nansum(np.abs(np.nanmean(normalizedData,0) - 1/Ncols))

        # fit the columns
        Sd = Sd * 1/Ncols * np.asarray(1/np.nanmean(normalizedData, 0))
        normalizedData = (dataT * Rd).T *Sd
        # calculate deviation from row marginals; column deviation is zero at odd indices. (index start = 0)
        convergenceTrail[2*i+1] = Ncols * 0.5 * np.nansum(np.abs(np.nanmean(normalizedData, 1) - 1/Ncols))

        convergence = convergenceTrail[2*i]
        i += 1

    return normalizedData, convergenceTrail, Rd, Sd

