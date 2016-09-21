#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that prepare the data before it can be normalized by CONSTANd. Includes:
* removing redundancy due to:
	* different peptide spectrum match (PSM) algorithms
	* different retention time (RT) values
	* different charges
* correct for isotopic impurities in the reporters
"""

import numpy as np


def addColumns(df, bools=None):
	""" Add extra columns to the dataFrame that will contain condensed information from to-be-deleted rows. """
	# if bools is None: add ALL.
	# allColumnNames = ['foo', 'bar', 'baz']
	# df.addcolumns(allColumnNames[bools])
	return df


def collapsePSMAlgo(df, master='mascot'):
	# TODO: retain deleted info in compact way
	# df = df.drop('column_name', 0) # 0 = row, 1 = column
	return df


def collapseRT(df, centerMeasure='mean'):
	# TODO: retain deleted info in compact way
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
	return df


def getIntensities(df):
	"""
	Extracts the (absolute) intensity matrix from the dataFrame.
	:param df:              pd.dataFrame    Pandas dataFrame from which to extract the intensities
	:return intensities:    np.ndArray      matrix with the intensities
	"""
	intensities = np.asarray(df[['126', '127', '128', '129', '130', '131']])


def setIntensities(df, intensities, location=[0,-1,0,-1]):
	""" Sets the intensities of the dataFrame at the specified location equal to the ndArray of given intensities."""
	#df[['126', '127', '128', '129', '130', '131']].iloc() = [intensities[:,1], intensities[:,2], intensities[:,3],
	#                                                  intensities[:,4], intensities[:,5], intensities[:,6]]
	return df
