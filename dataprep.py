#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that prepare the data before it can be normalized by CONSTANd. Includes:
* adding columns for retaining valuable information after collapsing
* removing high isolation interference cases
* removing redundancy due to:
	* different peptide spectrum match (PSM) algorithms
	* different retention time (RT) values
	* different charges
* correct for isotopic impurities in the reporters
* get/set the intensity matrix of a dataFrame
"""

import numpy as np


def addColumns(df, bools):
	""" Add extra columns to the dataFrame that will contain condensed information from to-be-deleted rows. """
	# if bools is None: add ALL.
	# allColumnNames = ['foo', 'bar', 'baz']
	# df.addcolumns(allColumnNames[bools])
	return df


def removeIsolationInterference(df, threshold):
	# TODO: retain deleted info in compact way (peptide sequence, scannumber, isolation interference?)
	# remove all rows with value > threshold
	return df


def collapsePSMAlgo(df, master, exclusive):
	# TODO: retain deleted info in compact way
	# df = df.drop('column_name', 0) # 0 = row, 1 = column
	# if exclusive: do not select those detected only by the slave
	return df


def collapseRT(df, centerMeasure_channels='mean', centerMeasure_intensities='mean', maxRelativeChannelVariance=None):
	# TODO: retain deleted info in compact way
	# what if the peptides resulting from the PSM do not agree between RT's? -> within-algorithm disagreement doesn't occur.
	# TODO: second switch: what if user wants not  a peak as high as the highest peak, but as high as the mean/median?
	# todo: check that the max RELATIVE variance on the channel intensities do not exceed given value.
	""" Averages over all peaks for each channel/reporter, and then rescales the magnitude of the resulting peak to
	match the magnitude of the largest constituent peak. In this way, the absolute intensity is still that of the
	largest peak, but the within-peak relative intensities are the average of all the constituent peaks. """
	# setIntensities(df, intensities, location)
	return df


def collapseCharge(df):
	# TODO: retain deleted info in compact way
	return df


def isotopicCorrection(df, correctionsMatrix):
	# solve the linear system
	# Observed(6,1) = correctionMatrix(6,6) * Real(6,1)
	# if Det(cM) = 0 no solution can be found.
	# returns corrected intensities also INSIDE df if user is paying # TODO this should be more clear in main()
	payingUser = False
	correctedIntensities = getIntensities(df)
	if payingUser:
		df = setIntensities(df, correctedIntensities)
	return df, correctedIntensities


def getIntensities(df):
	"""
	Extracts the (absolute) intensity matrix from the dataFrame.
	:param df:              pd.dataFrame    Pandas dataFrame from which to extract the intensities
	:return intensities:    np.ndArray      matrix with the intensities
	"""
	return np.asarray(df[['126', '127', '128', '129', '130', '131']])


def setIntensities(df, intensities, location=[0,-1,0,-1]):
	""" Sets the intensities of the dataFrame at the specified location equal to the ndArray of given intensities."""
	#df[['126', '127', '128', '129', '130', '131']].iloc() = [intensities[:,1], intensities[:,2], intensities[:,3],
	#                                                  intensities[:,4], intensities[:,5], intensities[:,6]]
	return df
