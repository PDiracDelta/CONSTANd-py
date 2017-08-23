#!/usr/bin/python3
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

from constandpp.tools import getIntensities
import numpy as np
import logging
from pandas import DataFrame, Series


def getAllPresentProteins(df):
	""" Returns the set of all master proteins appearing in at least one PSM, regardless of the PSM usefulness. """
	from constandpp.tools import partition, unnest
	if 'Master Protein Accessions' in df.columns.values:
		allMPAStrings = df.loc[:, 'Master Protein Accessions'].astype(str)
		allMPAListsAndItems = [x.split('; ') for x in allMPAStrings]
		# allMPALists, allMPAElements = partition(lambda x: len(x) > 1, allMPAListsAndItems)  # sort into
		# allMPAListElements = unnest(allMPALists)
		# return set(unnest(allMPAElements).extend(allMPAListElements))
		return set(unnest(allMPAListsAndItems))
	else:
		logging.warning("No master protein accessions found! You will not get a useful differential expression result.")
		return None


def removeObsoleteColumns(df, wantedColumns):
	"""
	Returns a dataFrame with only the specified columns of the input dataFrame, IF they exist.
	:param df:                  pd.dataFrame    input dataFrame
	:param wantedColumns:       list            specified columns
	:return:                    pd.dataFrame    dataFrame with only the specified columns of the input dataFrame
	"""
	obsolete = set(df.columns.values).difference(set(wantedColumns))
	return df.drop(list(obsolete), axis=1)


def removeMissing(df, noMissingValuesColumns, quanColumns):
	"""
	Removes detections for which entries in essential columns is missing, or which have no quan values or labels.
	:param df:  					pd.dataFrame    data with missing values
	:param noMissingValuesColumns:	list			columns which may not contain missing values
	:param quanColumns:		list			columns that contain the quantification values
	:return df:						pd.dataFrame    data without missing values
	"""
	toDelete = []
	try:
		for column in noMissingValuesColumns:
			# delete all detections that have a missing value in this column
			toDelete.extend(df.loc[df.loc[:, column].isnull(), :].index)
		# delete all detections that have a missing value in both columns: XCorr and Ions Score
		toDelete.extend(df.loc[[x and y for x, y in zip(df['XCorr'].isnull(), df['Ions Score'].isnull())]].index)
	except KeyError as e:
		raise KeyError("Required column '" + str(e.args[0]) + "' was not found.")
	# get the indices of all detections which have no quan values at all (those have their nansum equal to zero)
	noIntensitiesBool = np.nansum(getIntensities(df=df, quanColumns=quanColumns), axis=1) == 0.
	toDelete.extend(df.index[noIntensitiesBool])
	
	toDelete = np.unique(toDelete)
	removedData = df.loc[toDelete]
	if toDelete.size > 0:
		logging.warning(
			"Some detections have been removed from the workflow due to missing values: see removedData['missing'].")
	df.drop(toDelete, inplace=True)
	return df, removedData


def removeBadConfidence(df, minimum, removalColumnsToSave):
	"""
	Removes detections from the input dataFrame if they have a confidence level worse than the given minimum. Saves some
	info about data with lower than minimum confidence levels in removedData.
	Information of removed entries is saved according to removalColumnsToSave.
	:param df:              		pd.dataFrame    data with all confidence levels
	:param minimum:         		str             minimum confidence level
	:param removalColumnsToSave:	list			fields to save if an entry gets removed
	:return df:             		pd.dataFrame    data with confidence levels > minimum
	:return removedData:			pd.dataFrame    data with confidence levels < minimum
	"""
	if 'Confidence' not in df.columns:
		logging.warning("No 'Confidence' column found: did NOT filter on Confidence level!")
		return df, DataFrame()
	columnsToSave = ['Confidence'] + removalColumnsToSave
	allConfidenceLevels = ('low', 'Low', 'medium', 'Medium', 'high', 'High')
	minIndex = allConfidenceLevels.index(minimum)
	allowedConfidenceLevels = allConfidenceLevels[int(minIndex-np.ceil(minIndex % 2)):]
	badConfidences = [x not in allowedConfidenceLevels for x in df.loc[:, 'Confidence']]  # True if confidence level is bad
	toDelete = df.loc[badConfidences, :].index  # indices of rows to delete
	if len(set(df.loc[toDelete, 'Confidence']).difference(set(allConfidenceLevels))) > 0:  # illegal values detected
		logging.warning("Either the Confidence column is missing or it contains illegal values (allowed: Low, Medium, High). I am removing all PSMs with illegal values.")
	try:
		removedData = df.loc[toDelete, columnsToSave]
	except KeyError as e:
		removedData = DataFrame()
		logging.warning("Could not save removedData in removeBadConfidence step because a data column is missing: " + str(e.args[0]))
	df.drop(toDelete, inplace=True)
	return df, removedData


def removeIsolationInterference(df, threshold, removalColumnsToSave):
	"""
	Remove the data where there is too much isolation interference (above threshold) and return the remaining dataFrame
	along with info about the deletions.
	Information of removed entries is saved according to removalColumnsToSave.
	:param df:              		pd.dataFrame    unfiltered data
	:param threshold:       		float           remove all data with isolation interference above this value
	:param removalColumnsToSave:	list			fields to save if an entry gets removed
	:return df:             		pd.dataFrame    filtered data
	:return removedData:    		pd.dataFrame    basic info about the removed values
	"""
	if 'Isolation Interference [%]' not in df.columns:
		logging.warning("No 'Isolation Interference' column found: did NOT filter on Isolation Interference!")
		return df, DataFrame()
	columnsToSave = ['Isolation Interference [%]'] + removalColumnsToSave
	toDelete = df.loc[df['Isolation Interference [%]'] > threshold].index  # indices of rows to delete
	try:
		removedData = df.loc[toDelete, columnsToSave]
	except KeyError as e:
		removedData = DataFrame()
		logging.warning("Could not save removedData in removeIsolationInterference step because a data column is missing: " + str(e.args[0]))
	df.drop(toDelete, inplace=True)
	return df, removedData


def setMasterProteinDescriptions(df):
	"""
	Takes a dataframe and removes all non-master protein accessions (entire column) and descriptions (selective).
	:param df:  pd.dataFrame     dataFrame with all descriptions and accessions
	:return df: pd.dataFrame     dataFrame with only Master Protein descriptions and accessions
	"""
	try:
		masterProteinsLists = df.loc[:, 'Master Protein Accessions'].astype(str).apply(lambda x: x.split('; '))
		# [[all Proteins] per peptide]
		proteinsLists = df.loc[:, 'Protein Accessions'].astype(str).apply(lambda x: x.split('; '))
		# [[all descriptions] per peptide]
		descriptionsLists = df.loc[:, 'Protein Descriptions'].astype(str).apply(lambda x: x.split('; '))
		# [[indices of master descriptions with respect to list of all descriptions] per peptide]
		correctIndicesLists = [[proteins.index(masterProtein) for masterProtein in
								list(filter(lambda x: x.lower() != 'nan', masterProteins))]
							   for masterProteins, proteins in zip(masterProteinsLists, proteinsLists)]
		# [[master descriptions] per peptide]
		df.loc[:, 'Protein Descriptions'] = ['; '.join([descriptionsList[i] for i in correctIndices]) for
											 (descriptionsList, correctIndices) in
											 zip(descriptionsLists, correctIndicesLists)]
		df.drop('Protein Accessions', axis=1, inplace=True)
	except KeyError as e:
		logging.warning("Not all necessary columns found (see below); Adding 'Protein Descriptions' column with all empty strings.\n"+str(e))
		# add a descriptions column that has all empty strings
		df.loc[:, 'Protein Descriptions'] = Series(['', ]*len(df)).astype(str)
	return df


def undoublePSMAlgo(df, identifyingNodes, exclusive, quanColumns, removalColumnsToSave):
	"""
	Removes redundant data due to different PSM algorithms producing the same peptide match. The 'master' algorithm
	values are preferred over the 'slave' algorithm values, the latter whom are removed and have their basic information
	saved in removedData. If exclusive=True, this function only keeps master data and never slave data.
	If no slave algorithm was specified, the algorithm proceeds as if exclusive==True.
	:param df:              	pd.dataFrame    data with double First Scan numbers due to PSMAlgo redundancy
	:param identifyingNodes:		dict            master-slave PSM algorithm specifier
	:param exclusive:       		bool            save master data exclusively or include slave data where necessary?
	:param quanColumns:				list			columns that contain the quantification values
	:param removalColumnsToSave:	list			fields to save if an entry gets removed
	:return df:             		pd.dataFrame    data without double First Scan numbers due to PSMAlgo redundancy
	:return removedData:    		pd.dataFrame    basic info about the removed entries
	"""
	masterName = identifyingNodes['master'][0]
	byIdentifyingNodeDict = df.groupby('Identifying Node Type').groups  # {Identifying Node Type : [list of indices]}
	masterIndices = set(byIdentifyingNodeDict[masterName])
	toDelete = set(df.index.values).difference(masterIndices)  # all indices of detections not done by MASTER
	columnsToSave = removalColumnsToSave + quanColumns
	# check whether some slave data should be kept.
	if len(identifyingNodes['slaves']) != 0 and not exclusive:  # a SLAVE was specified: save its information in
		# removedData. ALso, KEEP UNIQUE SLAVE scans by removing them from the toDelete list
		byFirstScanDict = df.groupby('First Scan').groups
		singles = set(map(lambda e: e[0], filter(lambda e: len(e) == 1,
												 byFirstScanDict.values())))  # indices of detections done by only 1 PSMAlgo
		singlesNotByMasterIndices = singles.difference(masterIndices)
		toDelete = toDelete.difference(singlesNotByMasterIndices)  # keep only indices not discovered by SLAVE
		slaveScoreName = identifyingNodes['slaves'][0][1]
		columnsToSave.append(slaveScoreName)
	removedData = df.loc[toDelete, columnsToSave]
	df.drop(toDelete, inplace=True)  # drop all or some non-master data.
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
	for row in intensities:
		if not np.isnan(row).any():
			correctedIntensities.append(np.linalg.solve(correctionsMatrix, row))
		else:
			correctedIntensities.append(row)
			noCorrectionIndices.append(
				np.where(intensities == row)[0][0])  # np.where()[0][0] is numpy equivalent van .index()
			if not warnedYet:
				logging.warning(
					"Cannot correct isotope impurities for detections with NaN reporter intensities; skipping those.")
				warnedYet = True
	return np.asarray(correctedIntensities), noCorrectionIndices
