#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
CONSTANd
Normalizes the data matrix <data> by raking the m by n matrix such that
the row mean and column mean equals to 1/n. Missing information needs to
be presented as nan values and not as zero values because CONSTANd
employs the Matlab/NumPy functionality 'nanmean' that is able to ignore
nan-values when calculating the mean. The variable <nbIterations> is an
integer value that denotes the number of raking cycles. The variable <h>
defines the stopping criteria based on the L1-norm as defined by
Friedrich Pukelsheim, Bruno Simeone in "On the Iterative Proportional
Fitting Procedure: Structure of Accumulation Points and L1-Error
Analysis"

© Dirk Valkenborg & Jef Hooyberghs, 2014
Ported to Python by Joris Van Houtven, 2016
"""

import numpy as np
from numpy import nan


def constand(data, accuracy=1e-5, maxIterations=50):
	"""
	Return the normalized version of the input data (matrix) as an ndarray, as well as the convergence trail (residual
	error after each iteration) and the row and column multipliers R and S.
	:param data:                np.ndArray  (N,6) absolute intensities
	:param accuracy:            float       combined allowed deviation (residual error) of col and row means from 1/6
	:param maxIterations:       int         maximum amount of iterations (1x row and 1x col per iteration)
	:return normalizedData:     np.ndArray  (N,6) normalized intensities
	:return convergenceTrail:   np.ndArray  list of the residual error after each iteration
	:return R:                  np.ndArray  (N,) row multipliers
	:return S:                  np.ndArray  (6,) column multipliers
	"""
	assert isinstance(data, np.ndarray) and data.dtype in ['float64', 'int64']
	assert accuracy > 0
	assert maxIterations > 0

	# initialize variables
	Nrows, Ncols = data.shape
	convergenceTrail = np.asarray([nan]*(2*maxIterations))
	convergence = np.inf
	normalizedData = data

	i = 0  # number of current iteration
	# without reshaping the ndarrays, they have shape (x,) (no second value) and the procedure fails.
	R = np.ones((Nrows))  # R matrix diagonal from RAS
	S = np.ones((Ncols))  # S matrix diagonal from RAS
	# main loop; iterates until convergence is reached (i.e., L1-norm below variable <h>) or the maximum number of
	# iteration cycles is surpassed.
	while convergence > accuracy and i < maxIterations:

		# fit the rows
		Ri = 1/Ncols * np.asarray(1/np.nanmean(normalizedData, 1)).reshape(Nrows,)
		normalizedData = (normalizedData.T * Ri).T
		R *= Ri
		# calculate deviation from column marginals; row deviation is zero at even indices. (index start = 0)
		convergenceTrail[2*i] = Nrows * 0.5 * np.nansum(np.abs(np.nanmean(normalizedData, 0) - 1/Ncols))

		# fit the columns
		Si = 1/Ncols * np.asarray(1/np.nanmean(normalizedData, 0)).reshape(Ncols,)
		normalizedData *= Si
		S *= Si
		# calculate deviation from row marginals; column deviation is zero at odd indices. (index start = 0)
		convergenceTrail[2*i+1] = Ncols * 0.5 * np.nansum(np.abs(np.nanmean(normalizedData, 1) - 1/Ncols))

		convergence = convergenceTrail[2*i+1]
		i += 1

	return normalizedData, convergenceTrail[~np.isnan(convergenceTrail)], R, S

