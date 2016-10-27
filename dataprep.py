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
from warnings import warn

intensityColumns = None


def setIntensityColumns(intensityColumns):
	"""
	Sets the value of the global variable intensityColumns for use in the module functions.
	:param intensityColumns: list   names of the columns that contain the MS2 intensities
	"""
	globals()['intensityColumns'] = intensityColumns


def selectRequiredColumns(df, requiredColumns):
	"""
	Returns a dataFrame with only the specified columns of the input dataFrame.
	:param df:                  pd.dataFrame    input dataFrame
	:param requiredColumns:    list            specified columns
	:return:                    pd.dataFrame    dataFrame with only the specified columns of the input dataFrame
	"""
	return df[requiredColumns]


def removeMissing(df):
	"""
	Removes detections for which entries in essential columns is missing, or which have no quan values or labels.
	:param df:  pd.dataFrame    with missing values
	:return df: pd.dataFrame    without missing values
	"""
	toDelete = []
	for column in ['First Scan', 'Annotated Sequence', 'Identifying Node', 'Charge', 'Modifications']:
		# delete all detections that have a missing value in this column
		toDelete.extend(df[df[column].isnull()].index)
	# delete all detections that have a missing value in both columns: XCorr and Ions Score
	toDelete.extend(df[[x and y for x, y in zip(df['XCorr'].isnull(), df['Ions Score'].isnull())]].index)
	# delete all detections which have no quan values or no quan labels
	toDelete.extend(df[df['Quan Info'] == 'NoQuanValues'].index)
	toDelete.extend(df[df['Quan Info'] == 'NoQuanLabels'].index)
	toDelete = np.unique(toDelete)
	removedData = df.loc[toDelete]
	if toDelete.size > 0:
		warn("Some detections have been removed before due to missing values: see removedData['missing'].")
	df.drop(toDelete, inplace=True)
	return df, removedData


def removeBadConfidence(df, minimum):
	"""
	Removes detections from the input dataFrame if they have a confidence level worse than the given minimum. Saves some
	info about data with lower than minimum confidence levels in removedData.
	:param df:              pd.dataFrame    data with all confidence levels
	:param minimum:         str             minimum confidence level
	:return df:             pd.dataFrame    data with confidence levels > minimum
	:return removedData:    pd.dataFrame    data with confidence levels < minimum
	"""
	columnsToSave = ['First Scan', 'Annotated Sequence', 'Identifying Node', 'Master Protein Accessions', 'Confidence']
	conf2int = {'Low': 1, 'Medium': 2, 'High': 3}
	try:
		minimum = conf2int[minimum]
		badConfidences = [conf2int[x] < minimum for x in df['Confidence']]
	except KeyError:
		raise KeyError("Illegal Confidence values (allowed: Low, Medium, High). Watch out for capitalization.")
	toDelete = df[badConfidences].index  # indices of rows to delete
	removedData = df.loc[toDelete,columnsToSave]
	df.drop(toDelete, inplace=True)
	return df, removedData


def removeIsolationInterference(df, threshold):
	"""
	Remove the data where there is too much isolation interference (above threshold) and return the remaining dataFrame
	along with info about the deletions.
	:param df:              pd.dataFrame    unfiltered data
	:param threshold:       float           remove all data with isolation interference above this value
	:return df:             pd.dataFrame    filtered data
	:return removedData:    pd.dataFrame    basic info about the removed values
	"""
	columnsToSave = ['First Scan', 'Annotated Sequence', 'Identifying Node', 'Master Protein Accessions', 'Isolation Interference [%]']
	toDelete = df[df['Isolation Interference [%]'] > threshold].index # indices of rows to delete
	removedData = df.loc[toDelete,columnsToSave]
	df.drop(toDelete, inplace=True)
	return df, removedData


def undoublePSMAlgo(df, master, exclusive):
	"""
	Removes redundant data due to different PSM algorithms producing the same peptide match. The 'master' algorithm
	values are preferred over the 'slave' algorithm values, the latter whom are removed and have their basic information
	saved in removedData. If exclusive=true, this function only keeps master data (and saves slave(s) basic info).
	:param df:              pd.dataFrame    data with double First Scan numbers due to PSMAlgo redundancy
	:param master:          string          master PSM algorithm (master/slave relation)
	:param exclusive:       bool            save master data exclusively or include slave data where necessary?
	:return df:             pd.dataFrame    data without double First Scan numbers due to PSMAlgo redundancy
	:return removedData:    pd.dataFrame    basic info about the removed entries
	"""
	byFirstScanDict = df.groupby('Identifying Node').groups # {Identifying Node : [list of indices]}
	if master == 'mascot':
		columnsToSave = ['First Scan', 'Annotated Sequence', 'Master Protein Accessions', 'XCorr']
		toDelete = set(df.index.values).difference(set(byFirstScanDict['Mascot (A6)'])) # all indices not discovered by Mascot
		if not exclusive:
			toDelete = toDelete.difference(byFirstScanDict['Sequest HT (A2)']) # indices not discovered by Sequest either
	elif master == 'sequest':
		columnsToSave = ['First Scan', 'Annotated Sequence', 'Master Protein Accessions', 'Ions Score']
		toDelete = set(df.index.values).difference(set(byFirstScanDict['Sequest HT (A2)''Mascot (A6)']))  # all indices not discovered by Sequest
		if not exclusive:
			toDelete = toDelete.difference(byFirstScanDict['Mascot (A6)'])  # indices not discovered by Mascot either

	removedData = df.loc[toDelete,columnsToSave]
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
		df.loc[index,intensityColumns] = intensitiesDict[index]
	return df
