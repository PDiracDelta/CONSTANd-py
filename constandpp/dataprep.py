#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that prepare the data before it can be normalized by CONSTANd.
Includes:
* removing unnecessary variables/columns
* removing detections with missing values that are essential
* removing high isolation interference cases
* removing redundancy due to different peptide spectrum match (PSM) algorithms
* correct for isotopic impurities in the reporters
* get/set the intensity matrix of a dataFrame
Excludes (see collapse.py):
* removing redundancy due to:
	* different retention time (RT) values
	* different charges
	* different (post-translational) modifications (PTMs)
"""

import numpy as np

intensityColumns = None


def setIntensityColumns(ics):
	"""
	Sets the value of the global variable intensityColumns for use in the module functions.
	:param intensityColumns: list   names of the columns that contain the MS2 intensities
	"""
	globals()['intensityColumns'] = ics
	print(intensityColumns)


def selectEssentialColumns(df, essentialColumns):
	"""
	Returns a dataFrame with only the specified columns of the input dataFrame.
	:param df:                  pd.dataFrame    input dataFrame
	:param essentialColumns:    list            specified columns
	:return:                    pd.dataFrame    dataFrame with only the specified columns of the input dataFrame
	"""
	return df[essentialColumns]


def removeBadConfidence(df, minimum):
	# remove detections with confidence < minimum
	# map strings to integers
	# WATCH OUT FOR CAPITALIZATION
	return df # TODO


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
	removedData = df.loc[toDelete][colsToSave]
	df.drop(toDelete, inplace=True)
	return df, removedData


def undoublePSMAlgo(df, master, exclusive):
	"""
	Removes redundant data due to different PSM algorithms producing the same peptide match. The 'master' algorithm
	values are taken over the 'slave' algorithm values, the latter whom are removed and have their basic information
	saved in removedData. If exclusive=true, this function only keeps master data (and saves slave(s) basic info).
	:param df:              pd.dataFrame    unfiltered data
	:param master:          string          master PSM algorithm (master/slave relation)
	:param exclusive:       bool            save master data exclusively or include slave data where necessary?
	:return df:             pd.dataFrame    collapsed data
	:return removedData:    pd.dataFrame    basic info about the removed entries
	"""
	if master == 'mascot':
		colsToSave = ['First Scan', 'Annotated Sequence', 'Master Protein Accessions', 'XCorr']
		if exclusive:
			toDelete = df[df['Identifying Node'] == 'Sequest HT (A2)'].index
		else:
			toDelete = df[(df['Identifying Node'] == 'Sequest HT (A2)') * (df['Quan Info'] == 'Redundant')].index
	elif master == 'sequest':
		colsToSave = ['First Scan', 'Annotated Sequence', 'Master Protein Accessions', 'Ions Score']
		if exclusive:
			toDelete = df[df['Identifying Node'] == 'Mascot (A6)'].index
		else:
			toDelete = df[(df['Identifying Node'] == 'Mascot (A6)') * (df['Quan Info'] == 'Redundant')].index
	else:
		raise Exception(
			"Invalid master PSM algorithm: '" + master + "'. Please pick 'mascot' or 'sequest'.")
	removedData = ('master: '+master, df.loc[toDelete][colsToSave])
	dflen=df.shape[0] # TEST
	df.drop(toDelete, inplace=True)
	assert(dflen == df.shape[0]+removedData[1].shape[0]) # TEST

	return df, removedData


def isotopicCorrection(intensities, correctionsMatrix):
	"""
	Corrects isotopic impurities in the intensities using a given corrections matrix by solving the linear system:
	Observed(6,1) = correctionMatrix(6,6) * Real(6,1)
	:param intensities:         np.ndarray  array of arrays with the uncorrected intensities
	:param correctionsMatrix:   np.matrix   matrix with the isotopic corrections
	:return:                    np.ndarray  array of arrays with the corrected intensities
	"""
	correctedIntensities = []
	for row in intensities:
		correctedIntensities.append(np.linalg.solve(correctionsMatrix, row))
	return np.asarray(correctedIntensities)


def getIntensities(df):
	"""
	Extracts the (absolute) intensity matrix from the dataFrame.
	:param df:              pd.dataFrame    Pandas dataFrame from which to extract the intensities
	:return intensities:    np.ndArray      matrix with the intensities
	"""
	return np.asarray(df[intensityColumns])


def setIntensities(df, intensitiesDict):
	"""
	Sets the intensities of the dataFrame at the specified location equal to the ndArray of given intensities.
	:param df:              pd.dataFrame    input dataFrame
	:param intensitiesDict: dict            dict {index:[values]} with index and values of all df entries to be modified
	:return df:             pd.dataFrame    output dataFrame with updated intensities
	"""
	for index in intensitiesDict.keys():
		df.loc[index][intensityColumns] = intensitiesDict[index]
	return df
