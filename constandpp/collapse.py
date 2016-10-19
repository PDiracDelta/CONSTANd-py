#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that help removing redundancy in the data due to detections which have:
* different retention time (RT) values
* different charges
* different (post-translational) modifications (PTMs)
and replaces the duplicates with one representative detection and a combination/summary/selection of their intensities.
"""

import numpy as np
from dataprep import setIntensities, getIntensities

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
	if centerMeasure == 'geometricMedian':
		pass
	if centerMeasure == 'weighted':
		pass
	return newIntensities # TODO


def getRepresentative(duplicatesDf, duplicatesDict, representativeMethod='bestmatch'):
	# get the detection with the best PSM match
	# values of BEST PSM detection in all duplicates (master/slave-wise best)
	return detection # TODO


def getNewIntensities(duplicatesDf, duplicatesDict, method, maxRelativeReporterVariance):
	"""
	Combines the true duplicates' intensities into one new entry per first occurrence, conform the duplicatesDict structure.
	:param duplicatesDict:          dict            {firstOccurrenceIndex:[duplicateIndices]}
	:param duplicatesDf:            pd.dataFrame    data of only the first occurrences and duplicates
	:return weightedMS2Intensities: dict            {firstOccurrenceIndex:np.array(newIntensities)}
	"""
	import warnings
	weightedMS2Intensities = {}  # dict with the new MS2 intensities for each firstOccurrence
	if method == 'bestMatch':
		newIntensities = None
		representative = getRepresentative(duplicatesDf, duplicatesDict)
		setIntensities(representative, newIntensities)
		pass # TODO
	elif method == 'mostIntense':
		newIntensities = None
		representative = getRepresentative(duplicatesDf, duplicatesDict)
		setIntensities(newIntensities, newIntensities)
		pass # TODO
	else:
		newIntensities = combineDetections(duplicatesDf, centerMeasure=method)
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


def collapse(toCollapse, df, colsToSave, method, maxRelativeReporterVariance):
	"""
	Generic collapse function. Looks for duplicate 'Annotated Sequence' values in the dataFrame and verifies
	true duplication using checkTrueDuplicates function. Modifies df according to true duplicates and newly acquired
	intensities (via getNewIntensities function): remove all duplicates and enter one replacement detection.
	Returns removedData according to the colsToSave list.
	:param toCollapse:                  str             variable of which true duplicates are to be collapsed.
	:param df:                          pd.dataFrame    with sequence duplicates due to difference in certain variables/columns.
	:param colsToSave:                  list            list of variables to be saved for detections that ought to be removed
	:param method:                      str             defines how the new detection is to be selected/constructed
	:param maxRelativeReporterVariance: float           UNUSED value that restricts reporter variance
	:return df:                         pd.dataFrame    without sequence duplicates according to to checkTrueDuplicates.
	:return removedData:                dict            {firstOccurrenceIndex : [annotated_sequence, [other, values, to, be, saved] for each duplicate]}
	"""

	if toCollapse == 'RT':
		def checkTrueDuplicates(x, y):
			"""
			Checks whether dataFrame entries x and y are truly duplicates only due to RT difference.
			:param x:   pd.Sequence candidate firstOccurrence data
			:param y:   pd.Sequence candidate duplicate data
			"""
			if x['First Scan'] == y[
				'First Scan']:  # identical peptides have identical RT (you didn't do collapsePSMAlgo?)
				return False
			if x['Charge'] != y['Charge']:  # different charge = different ("modified, non-redundant") peptide
				return False
			if x['Modifications'] != y[
				'Modifications']:  # different PTM = different sequence (before collapsePTM() anyway).
				return False
			return True

		colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'RT [min]', 'MS2Intensity',
		              'PSMscore']

	elif toCollapse == 'Charge':
		def checkTrueDuplicates(x, y):
			"""
			Checks whether dataFrame entries x and y are truly duplicates only due to charge difference.
			:param x:   pd.Sequence candidate firstOccurrence data
			:param y:   pd.Sequence candidate duplicate data
			"""
			if x['First Scan'] == y[
				'First Scan']:  # identical peptides have identical RT (you didn't do collapsePSMAlgo?)
				assert False  # THIS SHOULD NOT BE REACHABLE ***UNLESS*** YOU DIDNT COLLAPSEPSMALGO() # TEST
				return False
			if x['Charge'] == y['Charge']:  # well obviously they should duplicate due to charge difference...
				return False
			if x['Modifications'] != y[
				'Modifications']:  # different PTM = different sequence (before collapsePTM() anyway).
				return False
			return True
		colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'Charge']

	elif toCollapse == 'PTM':
		def checkTrueDuplicates(x, y):
			"""
			Checks whether dataFrame entries x and y are truly duplicates only due to PTM difference.
			:param x:   pd.Sequence candidate firstOccurrence data
			:param y:   pd.Sequence candidate duplicate data
			"""
			if x['First Scan'] == y[
				'First Scan']:  # identical peptides have identical RT (you didn't do collapsePSMAlgo?)
				assert False  # THIS SHOULD NOT BE REACHABLE ***UNLESS*** YOU DIDNT COLLAPSEPSMALGO() # TEST
				return False
			if x['Charge'] != y[
				'Charge']:  # well obviously they should duplicate due to charge difference... # TODO should Charge variable still exist at this point?
				assert False  # THIS SHOULD NOT BE REACHABLE ***UNLESS*** YOU DIDNT COLLAPSECHARGE() # TEST
				return False
			if x['Modifications'] == y['Modifications']:
				assert False  # THIS SHOULD NOT BE REACHABLE ***UNLESS*** YOU DIDNT PERFORM ALL PREVIOUS COLLAPSES # TEST
				return False
			return True
		colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'Modifications']

	allSequences = df.groupby('Annotated Sequence').groups  # dict of SEQUENCE:[INDICES]
	allDuplicatesHierarchy = {}  # {firstOccurrence:[duplicates]}
	for sequence, indices in allSequences.items():
		if len(indices) > 1:  # only treat duplicated sequences
			# dict with duplicates per first occurrence, dataFrame with df indices but with only the duplicates
			duplicatesDict, duplicatesDf = getDuplicates(df, indices, checkTrueDuplicates)
			if False:  # TODO flag isolated peaks
				pass
			# get the new intensities per first occurrence index (df index)
			intensitiesDict = getNewIntensities(duplicatesDf, duplicatesDict, method, maxRelativeReporterVariance)
			allDuplicatesHierarchy.update(duplicatesDict)
	setIntensities(df, intensitiesDict)
	toDelete = list(allDuplicatesHierarchy.values())
	# save as {firstOccurrenceIndex : [annotated_sequence, [values, to, be, saved] for each duplicate]}
	removedData = dict((firstOccurrence, [df.loc[firstOccurrence][colsToSave[0]],
	                                      df.loc[allDuplicatesHierarchy[firstOccurrence]][colsToSave[1:]]])
	                   for firstOccurrence in allDuplicatesHierarchy.keys())
	df.drop(toDelete, inplace=True)

	return df, removedData
