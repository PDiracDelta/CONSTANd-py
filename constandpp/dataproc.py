#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that process the data before it can be normalized by CONSTANd.
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
Removed data is always saved into a removedData dataFrame.
"""

import numpy as np
from warnings import warn
from pandas import Series

intensityColumns = None
removalColumnsToSave = None
noMissingValuesColumns = None


def setProcessingGlobals(intensityColumns, removalColumnsToSave, noMissingValuesColumns):
	"""
	Sets the value of the global variable intensityColumns for use in the module functions.
	:param intensityColumns: list   names of the columns that contain the MS2 intensities
	"""
	globals()['intensityColumns'] = intensityColumns
	globals()['removalColumnsToSave'] = removalColumnsToSave
	globals()['noMissingValuesColumns'] = noMissingValuesColumns


def removeObsoleteColumns(df, wantedColumns):
	"""
	Returns a dataFrame with only the specified columns of the input dataFrame.
	:param df:                  pd.dataFrame    input dataFrame
	:param wantedColumns:       list            specified columns
	:return:                    pd.dataFrame    dataFrame with only the specified columns of the input dataFrame
	"""
	obsolete = set(df.columns.values).difference(set(wantedColumns))
	return df.drop(list(obsolete), axis=1)


def removeMissing(df):
	"""
	Removes detections for which entries in essential columns is missing, or which have no quan values or labels.
	:param df:  pd.dataFrame    with missing values
	:return df: pd.dataFrame    without missing values
	"""
	toDelete = []
	for column in noMissingValuesColumns:
		# delete all detections that have a missing value in this column
		toDelete.extend(df.loc[df.loc[:, column].isnull(), :].index)
	# delete all detections that have a missing value in both columns: XCorr and Ions Score
	toDelete.extend(df.loc[[x and y for x, y in zip(df['XCorr'].isnull(), df['Ions Score'].isnull())]].index)
	# delete all detections which have no quan values or no quan labels
	toDelete.extend(df.loc[df['Quan Info'] == 'NoQuanValues'].index)
	toDelete.extend(df.loc[df['Quan Info'] == 'NoQuanLabels'].index)
	toDelete = np.unique(toDelete)
	removedData = df.loc[toDelete]
	if toDelete.size > 0:
		warn("Some detections have been removed from the workflow due to missing values: see removedData['missing'].")
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
	columnsToSave = ['Confidence'] + removalColumnsToSave
	conf2int = {'Low': 1, 'Medium': 2, 'High': 3}
	try:
		minimum = conf2int[minimum]
		badConfidences = [conf2int[x] < minimum for x in df.loc[:, 'Confidence']]
	except KeyError:
		raise KeyError("Illegal Confidence values (allowed: Low, Medium, High). Watch out for capitalization.")
	toDelete = df.loc[badConfidences, :].index  # indices of rows to delete
	removedData = df.loc[toDelete, columnsToSave]
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
	columnsToSave = ['Isolation Interference [%]'] + removalColumnsToSave
	toDelete = df.loc[df['Isolation Interference [%]'] > threshold].index # indices of rows to delete
	removedData = df.loc[toDelete, columnsToSave]
	df.drop(toDelete, inplace=True)
	return df, removedData


def setMasterProteinDescriptions(df):
	"""
	Takes a dataframe and removes all non-master protein accessions (entire column) and descriptions (selective).
	:param df:  pd.dataFrame     dataFrame with all descriptions and accessions
	:return df: pd.dataFrame     dataFrame with only Master Protein descriptions and accessions
	"""
	if 'Protein Accessions' in df.columns.names:
		# [[master Proteins] per peptide]
		masterProteinsLists = df.loc[:, 'Master Protein Accessions'].astype(str).apply(lambda x: x.split('; '))
		# [[all Proteins] per peptide]
		proteinsLists = df.loc[:, 'Protein Accessions'].astype(str).apply(lambda x: x.split('; '))
		# [[all descriptions] per peptide]
		descriptionsLists = df.loc[:, 'Protein Descriptions'].astype(str).apply(lambda x: x.split('; '))
		# [[indices of master descriptions with respect to list of all descriptions] per peptide]
		correctIndicesLists = [[proteins.index(masterProtein) for masterProtein in list(filter(lambda x: x.lower()!='nan', masterProteins))]
			              for masterProteins, proteins in zip(masterProteinsLists, proteinsLists)]
		# [[master descriptions] per peptide]
		df.loc[:, 'Protein Descriptions'] = ['; '.join([descriptionsList[i] for i in correctIndices]) for (descriptionsList, correctIndices) in zip(descriptionsLists, correctIndicesLists)]
		df.drop('Protein Accessions', axis=1, inplace=True)
	else:
		warn("No 'Protein Accessions' column found; leaving all 'Protein Descriptions' intact.")
	return df


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
	byIdentifyingNodeDict = df.groupby('Identifying Node').groups # {Identifying Node : [list of indices]}
	if master == 'mascot':
		masterName = 'Mascot (A6)'
		slaveScoreName = 'XCorr'
	elif master == 'sequest':
		masterName = 'Sequest HT (A2)'
		slaveScoreName = 'Ions Score'
	columnsToSave = [slaveScoreName] + removalColumnsToSave
	masterIndices = set(byIdentifyingNodeDict[masterName])
	toDelete = set(df.index.values).difference(masterIndices)  # all indices of detections not done by MASTER
	if not exclusive:  # remove unique SLAVE scans from the toDelete list
		byFirstScanDict = df.groupby('First Scan').groups
		singles = set(map(lambda e: e[0], filter(lambda e: len(e) == 1,
		                                         byFirstScanDict.values())))  # indices of detections done by only 1 PSMAlgo
		singlesNotByMasterIndices = singles.difference(masterIndices)
		toDelete = toDelete.difference(singlesNotByMasterIndices)  # keep only indices not discovered by SLAVE
	removedData = df.loc[toDelete, columnsToSave]
	dflen=df.shape[0] # TEST
	df.drop(toDelete, inplace=True)
	assert(dflen == df.shape[0]+removedData.shape[0]) # TEST

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
	noCorrectionIndices = []
	warnedYet = False
	try:
		for row in intensities:
			if not np.isnan(row).any():
				correctedIntensities.append(np.linalg.solve(correctionsMatrix, row))
			else:
				correctedIntensities.append(row)
				noCorrectionIndices.append(np.where(intensities == row)[0][0]) # np.where()[0][0] is numpy equivalent van .index()
				if not warnedYet:
					warn("Cannot correct isotope impurities for detections with NaN reporter intensities; skipping those.")
	except TypeError:
		pass
	return np.asarray(correctedIntensities), noCorrectionIndices


def getIntensities(df, indices=None):
	"""
	Extracts the (absolute) intensity matrix from the dataFrame.
	:param df:              pd.dataFrame or pd.Series   Pandas dataFrame/Series from which to extract the intensities
	:param indices:         list                        indices of the detections for which to obtain the intensities
	:return intensities:    np.ndArray                  matrix with the intensities
	"""
	if isinstance(df, Series): # this is a dataframe with only 1 entry: indexing [:, cols] doesnt work.
		return np.asarray(df.loc[intensityColumns])
	if indices is None:
		#return np.asarray(df.loc(axis=1)[intensityColumns])
		return np.asarray(df.loc[:, intensityColumns])
	else:
		try: # TEST
			return np.asarray(df.loc[indices, intensityColumns])
		except:
			pass


def setIntensities(df, intensities):
	"""
	Sets the intensities of the dataFrame at the specified location equal to the array of given intensities, at the
	specified locations if a dict is provided instead of an array.
	:param df:              pd.dataFrame    input dataFrame
	:param intensities:     np.ndarray      matrix with MS2 intensities
							dict            dict {index:[values]} with index and values of all df entries to be modified
	:return df:             pd.dataFrame    output dataFrame with updated intensities
	"""
	if isinstance(intensities, np.ndarray):
		assert df.loc[:, intensityColumns].shape == intensities.shape
		df.loc[:, intensityColumns] = intensities
	elif isinstance(intensities, dict):
		for index in intensities.keys():
			df.loc[index, intensityColumns] = intensities[index]
	return df
