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
	#TODO: implement the use of this function so that it should not be needlessly ran again by each collapse function.
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

	candidatesDf = df.iloc[indices]  # create new dataframe with only the possible duplicates to reduce overhead.
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
				if checkTrueDuplicates(candidatesDf.iloc[i], candidatesDf.iloc[j]):
					duplicatesDict[i].append(j)  # mark index of rowj as a duplicate of index of rowi
					stillMatchableTemp.remove(j)
			stillMatchable = stillMatchableTemp  # now that iteration is done, modify.
	duplicatesDict = dict((x, y) for x, y in duplicatesDict.items() if y)  # remove empty lists of duplicates
	duplicatesDf = candidatesDf.iloc[
		list(duplicatesDict.keys()) + list(duplicatesDict.values())]  # df of only the duplicates
	return duplicatesDict, duplicatesDf


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
	removedData = ('master: '+master, df.iloc[toDelete][colsToSave])
	dflen=df.shape[0] # TEST
	df.drop(toDelete, inplace=True)
	assert(dflen == df.shape[0]+removedData[1].shape[0]) # TEST

	return df, removedData


def collapseRT(df, centerMeasure_channels='mean', centerMeasure_intensities='mean', maxRelativeReporterVariance=np.inf):
	# TODO: retain deleted info in compact way
	# what if the peptides resulting from the PSM do not agree between RT's? -> within-algorithm disagreement doesn't occur.
	# TODO: second switch: what if user wants not  a peak as high as the highest peak, but as high as the mean/median?
	# todo: check that the max RELATIVE variance on the channel intensities do not exceed given value. (better: read below)
	# todo: report when RT differences exceed a certain threshold
	""" Averages over all peaks for each channel/reporter, and then rescales the magnitude of the resulting peak to
	match the magnitude of the largest constituent peak. In this way, the absolute intensity is still that of the
	largest peak, but the within-peak relative intensities are the average of all the constituent peaks. """
	# setIntensities(df, intensities, location)
	import warnings

	def checkTrueDuplicates(x, y):  # RT close?
		"""
		Checks whether dataFrame entries x and y are truly duplicates only due to RT difference.
		:param x:   pd.Sequence candidate firstOccurrence data
		:param y:   pd.Sequence candidate duplicate data
		"""
		if x['First Scan'] == y['First Scan']: # identical peptides are always in the same MS1 sweep
			return False
		if x['RT [min]'] == y['RT [min]']:  # HOW CLOSE SHOULD THIS BE? HOW ARE THE MS2 SCANS BINNED BELONGING TO THE SAME MS1 SWEEP?
			return False  # the difference should be due to RT difference (different ms1 sweep) TODO this stuff ^^^
		return True

	def getNewIntensities(duplicatesDf, duplicatesDict): # sum intensities in a weighted way
		"""
		Combines the duplicates' intensities into one new entry per first occurrence, conform the duplicatesDict structure.
		:param duplicatesDict:          dict            {firstOccurrenceIndex:[duplicateIndices]}
		:param duplicatesDf:            pd.dataFrame    data of only the first occurrences and duplicates
		:return weightedMS2Intensities: dict            {firstOccurrenceIndex:np.array(newIntensities)}
		"""
		# TODO combine this function with the getNewIntensities from collapseCharge? Although you might not want to weigh by MS1 intensity here, you might just want the BEST MS1 peak height.
		weightedMS2Intensities = {} # dict with the new MS2 intensities for each firstOccurrence
		if centerMeasure_channels and centerMeasure_intensities and False: # TODO switch for the center measures.
			pass
		for firstOccurrence,duplicates in duplicatesDict:
			totalMS1Intensity = sum(duplicatesDf.iloc[[firstOccurrence]+duplicates]['Intensity'])
			allWeights = duplicatesDf.iloc[[firstOccurrence] + duplicates]['Intensity'] / totalMS1Intensity # TODO this is very probably NOT correct: you are weighting absolute MS2 intensities by MS1 intensity
			allMS2Intensities = getIntensities(duplicatesDf.iloc[[firstOccurrence]+duplicates]) # np.array
			weightedMS2Intensities[firstOccurrence] = np.sum((allMS2Intensities.T*allWeights).T,0) # TODO check if the dimension are correct
			if np.any(np.var(allMS2Intensities,0) > maxRelativeReporterVariance): # TODO this can only be consistent if executed on RELATIVE intensities.
				warnings.warn("maxRelativeReporterVariance too high for peptide with index "+firstOccurrence+".") #TODO this shouldnt just warn, you should also decide what to do.
		return weightedMS2Intensities # update the intensities

	colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'RT [min]', 'Intensity']
	allSequences = df.groupby('Annotated Sequence').groups  # dict of SEQUENCE:[INDICES]
	allDuplicatesHierarchy = {}  # {firstOccurrence:[duplicates]}
	for sequence,indices in allSequences.items():
		if len(indices)>1: # only treat duplicated sequences
			# dict with duplicates per first occurrence, dataFrame with df indices but with only the duplicates
			duplicatesDict, duplicatesDf = getDuplicates(df, indices, checkTrueDuplicates)
			if False:  # x['RT [min]'] is too far from y['RT [min]']
				pass  # keep track of this isolated peak # TODO implement something for isolated peaks?
			# get the new intensities per first occurrence index (df index)
			intensitiesDict = getNewIntensities(duplicatesDf, duplicatesDict)
			allDuplicatesHierarchy.update(duplicatesDict)
	setIntensities(df, intensitiesDict)
	toDelete = list(allDuplicatesHierarchy.values())
	# save as {firstOccurrenceIndex : [annotated_sequence, [values, to, be, saved] for each duplicate]}
	removedData = dict(firstOccurrence=[df.iloc[firstOccurrence][colsToSave[0]],
	                                    df.iloc[allDuplicatesHierarchy[firstOccurrence]][colsToSave[1:]]] for
	                   firstOccurrence in allDuplicatesHierarchy.keys())
	df.drop(toDelete, inplace=True)

	return df, removedData


def collapseCharge(df):
	"""
	Replaces identical sequence entries with a different charge by one containing the sum of their intensities. The
	'Charge' column is then deleted from the dataFrame. Deleted charges, scan numbers, ... are saved in removedData,
	each associated with its first occurrence duplicate index.
	:param df:              pd.dataFrame    with sequence duplicates due to difference in Charge.
	:return df:             pd.dataFrame    without sequence duplicates due to difference in Charge.
	:return removedData:    dict            {firstOccurrenceIndex : [annotated_sequence, [other, values, to, be, saved] for each duplicate]}
	"""

	def checkTrueDuplicates(x, y):  # RT close? same PSMAlgo->same first scan number->same charge?
		"""
		Checks whether dataFrame entries x and y are truly duplicates only due to charge difference.
		:param x:   pd.Sequence candidate firstOccurrence data
		:param y:   pd.Sequence candidate duplicate data
		"""
		if x['Charge'] == y['Charge']:  # well obviously they should duplicate due to charge difference...
			return False
		if x['RT [min]'] != y['RT [min]']:  # HOW CLOSE SHOULD THIS BE? HOW ARE THE MS2 SCANS BINNED BELONGING TO THE SAME MS1 SWEEP?
			return False  # the difference should not be due to RT difference (different ms1 sweep) TODO this stuff ^^^
		return True

	def getNewIntensities(duplicatesDf, duplicatesDict): # sum intensities in a weighted way
		"""
		Combines the duplicates' intensities into one new entry per first occurrence, conform the duplicatesDict structure.
		:param duplicatesDict:          dict            {firstOccurrenceIndex:[duplicateIndices]}
		:param duplicatesDf:            pd.dataFrame    data of only the first occurrences and duplicates
		:return weightedMS2Intensities: dict            {firstOccurrenceIndex:np.array(newIntensities)}
		"""
		weightedMS2Intensities = {} # dict with the new MS2 intensities for each firstOccurrence
		for firstOccurrence,duplicates in duplicatesDict:
			totalMS1Intensity = sum(duplicatesDf.iloc[[firstOccurrence]+duplicates]['Intensity'])
			allWeights = duplicatesDf.iloc[[firstOccurrence] + duplicates]['Intensity'] / totalMS1Intensity # TODO this is very probably NOT correct: you are weighting absolute MS2 intensities by MS1 intensity
			allMS2Intensities = getIntensities(duplicatesDf.iloc[[firstOccurrence]+duplicates]) # np.array
			weightedMS2Intensities[firstOccurrence] = np.sum((allMS2Intensities.T*allWeights).T,0) # TODO check if the dimension are correct
		return weightedMS2Intensities # update the intensities

	colsToSave = ['Annotated Sequence', 'Master Protein Accessions', 'First Scan', 'Charge', 'Intensity']
	allSequences = df.groupby('Annotated Sequence').groups  # dict of SEQUENCE:[INDICES]
	allDuplicatesHierarchy = {} # {firstOccurrence:[duplicates]}
	for sequence,indices in allSequences.items():
		if len(indices)>1: # only treat duplicated sequences
			# dict with duplicates per first occurrence, dataFrame with df indices but with only the duplicates
			duplicatesDict, duplicatesDf = getDuplicates(df, indices, checkTrueDuplicates)
			# get the new intensities per first occurrence index (df index)
			intensitiesDict = getNewIntensities(duplicatesDf, duplicatesDict)
			allDuplicatesHierarchy.update(duplicatesDict)
	setIntensities(df, intensitiesDict)
	toDelete = list(allDuplicatesHierarchy.values())
	# save as {firstOccurrenceIndex : [[values, to, be, saved] for each duplicate]}
	removedData = dict(firstOccurrence=[df.iloc[firstOccurrence][colsToSave[0]],
	                                    df.iloc[allDuplicatesHierarchy[firstOccurrence]][colsToSave[1:]]] for
	                   firstOccurrence in allDuplicatesHierarchy.keys())
	df.drop(toDelete, inplace=True)
	return df, removedData


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
	return np.asarray(df[intensityColumns])


def setIntensities(df, intensitiesDict):
	"""
	Sets the intensities of the dataFrame at the specified location equal to the ndArray of given intensities.
	:param df:              pd.dataFrame    input dataFrame
	:param intensitiesDict: dict            dict {index:[values]} with index and values of all df entries to be modified
	:return df:             pd.dataFrame    output dataFrame with updated intensities
	"""
	for index in intensitiesDict.keys():
		df.iloc[index][[intensityColumns]] = intensitiesDict[index]
	return df
