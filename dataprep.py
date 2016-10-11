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


def removeIsolationInterference(df, threshold):
	"""
	Remove the data where there is too much isolation interference (above threshold) and return the remaining dataFrame
	along with info about the deletions.
	:param df:              pd.dataFrame    unfiltered data
	:param threshold:       float           remove all data with isolation interference above this value
	:return df:             pd.dataFrame    filtered data
	:return removedData:    pd.dataFrame    basic info about the removed values
	"""
	colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'Isolation Interference [%]']
	toDelete = df[df['Isolation Interference [%]'] > threshold].index # indices of rows to delete
	removedData = df.iloc[toDelete][colsToSave]
	df.drop(toDelete, inplace=True)
	return df, removedData


def collapsePSMAlgo(df, master, exclusive):
	"""
	Removes redundant data due to different PSM algorithms producing the same peptide match. The 'master' algorithm
	values are taken over the 'slave' algorithm values, the latter whom are removed and have their basic information
	saved in removedData. If exclusive=true, this function only keeps master data (and saves slave(s) basic info).
	:param df:              pd.dataFrame    unfiltered data
	:param master:          string          master PSM algorithm (master/slave relation)
	:param exclusive:       bool            save master data exclusively or include slave data where necessary?
	:return df:             pd.dataFrame    collapsed data
	:return removedData:    pd.dataFrame    basic info about the removed values
	"""
	if master == 'mascot':
		colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'XCorr']
		if exclusive:
			toDelete = df[df['Identifying Node'] == 'Sequest HT (A2)'].index
		else:
			toDelete = df[(df['Identifying Node'] == 'Sequest HT (A2)') * (df['Quan Info'] == 'Redundant')].index
	elif master == 'sequest':
		colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'Ions Score']
		if exclusive:
			toDelete = df[df['Identifying Node'] == 'Mascot (A6)'].index
		else:
			toDelete = df[(df['Identifying Node'] == 'Mascot (A6)') * (df['Quan Info'] == 'Redundant')].index
	else:
		raise Exception(
			"Invalid master PSM algorithm: '" + master + "'. Please pick 'mascot' or 'sequest'.")
	removedData = ('master: '+master, df.iloc[toDelete][colsToSave])
	dflen=df.shape[0] # TEST
	df.drop(toDelete, inplace=True)
	assert(dflen == df.shape[0]+removedData[1].shape[0]) # TEST

	return df, removedData


def collapseRT(df, centerMeasure_channels='mean', centerMeasure_intensities='mean', maxRelativeChannelVariance=None):
	# TODO: retain deleted info in compact way
	# what if the peptides resulting from the PSM do not agree between RT's? -> within-algorithm disagreement doesn't occur.
	# TODO: second switch: what if user wants not  a peak as high as the highest peak, but as high as the mean/median?
	# todo: check that the max RELATIVE variance on the channel intensities do not exceed given value. (better: read below)
	# todo: report when RT differences exceed a certain threshold
	""" Averages over all peaks for each channel/reporter, and then rescales the magnitude of the resulting peak to
	match the magnitude of the largest constituent peak. In this way, the absolute intensity is still that of the
	largest peak, but the within-peak relative intensities are the average of all the constituent peaks. """
	# setIntensities(df, intensities, location)
	return df


def collapseCharge(df):
	# TODO: retain deleted info in compact way
	return df


def isotopicCorrection(intensities, correctionsMatrix):
	# solve the linear system
	# Observed(6,1) = correctionMatrix(6,6) * Real(6,1)
	# if Det(cM) = 0 no solution can be found.
	correctedIntensities = intensities
	return correctedIntensities


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
