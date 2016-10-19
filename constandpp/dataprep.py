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

intensityColumns = ['126', '127', '128', '129', '130', '131']


def getDuplicates(df, indices, checkTrueDuplicates):
	"""
	Takes a list of indices of candidate-duplicates (all df entries with identical annotated sequence) and returns
	a dict of first occurrences and their true duplicates due to charge difference, as well as the corresponding
	data extracted from the original dataFrame df. First occurrences without duplicates do not appear in the dict.
	:param df:                  pd.dataFrame    data which is to be checked for duplicates
	:param indices:             list            indices of the locations of candidate-duplicates in the original dataFrame df.
	:param checkTrueDuplicates: function        function which returns True if two entries given as arguments are true duplicates.
	:return duplicatesDict:     dict            {firstOccurrenceIndex:[duplicateIndices]}
	:return duplicatesDf:       pd.dataFrame    data of only the entries involved in true duplication due to charge
	"""
	import pandas as pd

	candidatesDf = df.loc[indices]  # create new dataframe with only the possible duplicates to reduce overhead.
	candidatesDf['index'] = pd.Series(indices)  # add the original indices as a column 'index'
	candidatesDf.set_index('index')  # set the index to the original indices
	duplicatesDict = {}  # keep a dict of which indices are duplicates of which
	stillMatchable = indices  # keep a dict of indices that can still be checked for having duplicates
	for i in indices:
		if i in stillMatchable:  # if i has already been matched as a duplicate
			stillMatchable.remove(i)
			duplicatesDict[i] = []  # (see above) keep a dict of which indices are duplicates of which
			stillMatchableTemp = stillMatchable  # cannot modify the variable while its being iterated over -> keep temp
			for j in stillMatchable:
				if checkTrueDuplicates(candidatesDf.loc[i], candidatesDf.loc[j]):
					duplicatesDict[i].append(j)  # mark index of rowj as a duplicate of index of rowi
					stillMatchableTemp.remove(j)
			stillMatchable = stillMatchableTemp  # now that iteration is done, modify.
	duplicatesDict = dict((x, y) for x, y in duplicatesDict.items() if y)  # remove empty lists of duplicates
	duplicatesDf = candidatesDf.loc[
		list(duplicatesDict.keys()) + list(duplicatesDict.values())]  # df of only the duplicates
	return duplicatesDict, duplicatesDf


def combineDetections(duplicatesDf, centerMeasure):
	if centerMeasure == 'mean':
		pass
	if centerMeasure == 'median':
		pass
	if centerMeasure == 'weighted':
		pass
	return newIntensities


def getNewIntensities(duplicatesDf, duplicatesDict, method, centerMeasure, maxRelativeReporterVariance):
	"""
	Combines the true duplicates' intensities into one new entry per first occurrence, conform the duplicatesDict structure.
	:param duplicatesDict:          dict            {firstOccurrenceIndex:[duplicateIndices]}
	:param duplicatesDf:            pd.dataFrame    data of only the first occurrences and duplicates
	:return weightedMS2Intensities: dict            {firstOccurrenceIndex:np.array(newIntensities)}
	"""
	import warnings
	weightedMS2Intensities = {}  # dict with the new MS2 intensities for each firstOccurrence
	if method == 'bestMatch':
		pass
	elif method == 'mostIntense':
		pass
	elif method == 'centerMeasure':
		newIntensities = combineDetections(duplicatesDf, centerMeasure)
	# TODO the next section is obsolete if you use combineDetections
	for firstOccurrence, duplicates in duplicatesDict:  # TODO flag PTM differences.
		totalMS1Intensity = sum(duplicatesDf.loc[[firstOccurrence] + duplicates]['Intensity'])
		allWeights = duplicatesDf.loc[[firstOccurrence] + duplicates][
			             'Intensity'] / totalMS1Intensity  # TODO this is very probably NOT correct: you are weighting absolute MS2 intensities by MS1 intensity
		allMS2Intensities = getIntensities(duplicatesDf.loc[[firstOccurrence] + duplicates])  # np.array
		weightedMS2Intensities[firstOccurrence] = np.sum((allMS2Intensities.T * allWeights).T,
		                                                 0)  # TODO check if the dimension are correct
		if np.any(np.var(allMS2Intensities,
		                 0) > maxRelativeReporterVariance):  # TODO this can only be consistent if executed on RELATIVE intensities.
			warnings.warn(
				"maxRelativeReporterVariance too high for peptide with index " + firstOccurrence + ".")  # TODO this shouldnt just warn, you should also decide what to do.
	return weightedMS2Intensities  # update the intensities


def collapse(df, checkTrueDuplicates, colsToSave, method, maxRelativeReporterVariance):
	"""
	Generic collapse function. Looks for duplicate 'Annotated Sequence' values in the dataFrame and verifies
	true duplication using checkTrueDuplicates function. Modifies df according to true duplicates and newly acquired
	intensities (via getNewIntensities function): remove all duplicates and enter one replacement detection.
	Returns removedData according to the colsToSave list.
	:param df:                          pd.dataFrame    with sequence duplicates due to difference in certain variables/columns.
	:param checkTrueDuplicates:         function        returns true if two detections are true duplicates in current context
	:param colsToSave:                  list            list of variables to be saved for detections that ought to be removed
	:param method:                      str             defines how the new detection is to be selected/constructed
	:param maxRelativeReporterVariance: float           UNUSED value that restricts reporter variance
	:return df:                         pd.dataFrame    without sequence duplicates according to to checkTrueDuplicates.
	:return removedData:                dict            {firstOccurrenceIndex : [annotated_sequence, [other, values, to, be, saved] for each duplicate]}
	"""
	allSequences = df.groupby('Annotated Sequence').groups  # dict of SEQUENCE:[INDICES]
	allDuplicatesHierarchy = {}  # {firstOccurrence:[duplicates]}
	for sequence, indices in allSequences.items():
		if len(indices) > 1:  # only treat duplicated sequences
			# dict with duplicates per first occurrence, dataFrame with df indices but with only the duplicates
			duplicatesDict, duplicatesDf = getDuplicates(df, indices, checkTrueDuplicates)
			if False:  # TODO flag isolated peaks
				pass
			# get the new intensities per first occurrence index (df index)
			intensitiesDict = getNewIntensities(duplicatesDf, duplicatesDict, method, centerMeasure, maxRelativeReporterVariance)
			allDuplicatesHierarchy.update(duplicatesDict)
	setIntensities(df, intensitiesDict)
	toDelete = list(allDuplicatesHierarchy.values())
	# save as {firstOccurrenceIndex : [annotated_sequence, [values, to, be, saved] for each duplicate]}
	removedData = dict((firstOccurrence, [df.loc[firstOccurrence][colsToSave[0]],
	                                      df.loc[allDuplicatesHierarchy[firstOccurrence]][colsToSave[1:]]])
	                   for firstOccurrence in allDuplicatesHierarchy.keys())
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
	colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'Isolation Interference [%]']
	toDelete = df[df['Isolation Interference [%]'] > threshold].index # indices of rows to delete
	removedData = df.loc[toDelete][colsToSave]
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
	:return removedData:    pd.dataFrame    basic info about the removed entries
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
	removedData = ('master: '+master, df.loc[toDelete][colsToSave])
	dflen=df.shape[0] # TEST
	df.drop(toDelete, inplace=True)
	assert(dflen == df.shape[0]+removedData[1].shape[0]) # TEST

	return df, removedData


def collapseRT(df, method, maxRelativeReporterVariance):
	# todo: check that the max RELATIVE variance on the channel intensities do not exceed given value. (better: read below)
	# todo: report when RT differences exceed a certain threshold
	"""
	Combines detections in the dataFrame that differ only in retention time but may have the same PTMs and Charge into a
	new detection. The duplicates are removed and replaced by this new detection. Essential info about removed data is
	also returned as removedData.
	:param df:                          pd.dataFrame    with sequence duplicates due to difference in certain variables/columns.
	:param method:                      str             defines how the new detection is to be selected/constructed
	:param maxRelativeReporterVariance: float           UNUSED value that restricts reporter variance
	:return df:                         pd.dataFrame    without sequence duplicates solely due to difference in RT.
	:return removedData:                dict            {firstOccurrenceIndex : [annotated_sequence, [other, values, to, be, saved] for each duplicate]}
	"""
	# setIntensities(df, intensities, location)
	import warnings

	def checkTrueDuplicates(x, y):
		"""
		Checks whether dataFrame entries x and y are truly duplicates only due to RT difference.
		:param x:   pd.Sequence candidate firstOccurrence data
		:param y:   pd.Sequence candidate duplicate data
		"""
		if x['First Scan'] == y['First Scan']: # identical peptides have identical RT (you didn't do collapsePSMAlgo?)
			return False
		if x['Charge'] != y['Charge']:  # different charge = different ("modified, non-redundant") peptide
			return False
		if x['Modifications'] != y['Modifications']: # different PTM = different sequence (before collapsePTM() anyway).
			return False
		return True

	colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'RT [min]', 'MS2Intensity', 'PSMscore']
	df, removedData = collapse(df, checkTrueDuplicates=checkTrueDuplicates, colsToSave=colsToSave, method=method, maxRelativeReporterVariance=maxRelativeReporterVariance)
	return df, removedData


def collapseCharge(df, method, maxRelativeReporterVariance): # no method parameter because always method=centerMeasure (see below)
	"""
	This function should always be preceeded by collapseRT(). Combines detections in the dataFrame that differ only in
	Charge but may have the same PTMs into a new detection. The duplicates are removed and replaced by this new detection.
	Essential info about removed data is also returned as removedData.
	:param method:                      str             defines how the new detection is to be selected/constructed
	:param maxRelativeReporterVariance: float           UNUSED value that restricts reporter variance
	:param df:                          pd.dataFrame    with sequence duplicates due to difference in Charge.
	:return df:                         pd.dataFrame    without sequence duplicates solely due to difference in Charge.
	:return removedData:                dict            {firstOccurrenceIndex : [annotated_sequence, [other, values, to, be, saved] for each duplicate]}
	"""

	def checkTrueDuplicates(x, y):
		"""
		Checks whether dataFrame entries x and y are truly duplicates only due to charge difference.
		:param x:   pd.Sequence candidate firstOccurrence data
		:param y:   pd.Sequence candidate duplicate data
		"""
		if x['First Scan'] == y['First Scan']:  # identical peptides have identical RT (you didn't do collapsePSMAlgo?)
			assert False # THIS SHOULD NOT BE REACHABLE ***UNLESS*** YOU DIDNT COLLAPSEPSMALGO() # TEST
			return False
		if x['Charge'] == y['Charge']:  # well obviously they should duplicate due to charge difference...
			return False
		if x['Modifications'] != y['Modifications']: # different PTM = different sequence (before collapsePTM() anyway).
			return False
		return True

	colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'Charge']
	df, removedData = collapse(df, checkTrueDuplicates=checkTrueDuplicates, colsToSave=colsToSave, method=method, maxRelativeReporterVariance=maxRelativeReporterVariance)
	# method=centerMeasure because after RT collapse the PSM algorithm scores may be gone
	return df, removedData


def collapsePTM(df, method, maxRelativeReporterVariance):
	"""
	Combines detections in the dataFrame that differ only in PTMs ('Modification') but may have the same RT or Charge
	into a new detection. The duplicates are removed and replaced by this new detection. Essential info about removed
	data is	also returned as removedData.
	:param maxRelativeReporterVariance:
	:param method:
	:param df:              pd.dataFrame    with sequence duplicates due to difference in PTMs.
	:return df:             pd.dataFrame    without sequence duplicates solely due to difference in PTMs.
	:return removedData:    dict            {firstOccurrenceIndex : [annotated_sequence, [other, values, to, be, saved] for each duplicate]}
	"""

	def checkTrueDuplicates(x, y):
		"""
		Checks whether dataFrame entries x and y are truly duplicates only due to PTM difference.
		:param x:   pd.Sequence candidate firstOccurrence data
		:param y:   pd.Sequence candidate duplicate data
		"""
		if x['First Scan'] == y['First Scan']:  # identical peptides have identical RT (you didn't do collapsePSMAlgo?)
			assert False # THIS SHOULD NOT BE REACHABLE ***UNLESS*** YOU DIDNT COLLAPSEPSMALGO() # TEST
			return False
		if x['Charge'] != y['Charge']:  # well obviously they should duplicate due to charge difference... # TODO should Charge variable still exist at this point?
			assert False  # THIS SHOULD NOT BE REACHABLE ***UNLESS*** YOU DIDNT COLLAPSECHARGE() # TEST
			return False
		if x['Modifications'] == y['Modifications']:
			assert False # THIS SHOULD NOT BE REACHABLE ***UNLESS*** YOU DIDNT PERFORM ALL PREVIOUS COLLAPSES # TEST
			return False
		return True

	colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'Modifications']
	df, removedData = collapse(df, checkTrueDuplicates=checkTrueDuplicates, colsToSave=colsToSave, method=method, maxRelativeReporterVariance=maxRelativeReporterVariance)
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
		df.loc[index][[intensityColumns]] = intensitiesDict[index]
	return df
