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
from dataprep import intensityColumns, setIntensities, getIntensities
from warnings import warn

columnsToSave = None


def setCollapseColumnsToSave(columnsToSave):
	"""
	Sets the value of the global variable columnsToSave for use in the module functions.
	:param columnsToSave: list   names of the columns that ought to be saved when removing data in a collapse.
	"""
	globals()['columnsToSave'] = columnsToSave


def collapse(toCollapse, df, method, maxRelativeReporterVariance, masterPSMAlgo, undoublePSMAlgo_bool): #
	"""
	Generic collapse function. Looks for duplicate 'Annotated Sequence' values in the dataFrame and verifies
	true duplication using checkTrueDuplicates function. Modifies df according to true duplicates and newly acquired
	intensities (via getNewIntensities function): remove all duplicates and enter one replacement detection.
	Adds a 'Degeneracy' column to the dataFrame if it didn't exist already: this contains the number of peptides that
	have been collapsed onto that (synthetic) detection.
	Returns removedData according to the columnsToSave list.
	:param toCollapse:                  str             variable of which true duplicates are to be collapsed.
	:param df:                          pd.dataFrame    with sequence duplicates due to difference in certain variables/columns.
	:param columnsToSave:                  list            list of variables to be saved for detections that ought to be removed
	:param method:                      str             defines how the new detection is to be selected/constructed
	:param maxRelativeReporterVariance: float           UNUSED value that restricts reporter variance
	:return df:                         pd.dataFrame    without sequence duplicates according to to checkTrueDuplicates.
	:return removedData:                dict            {firstOccurrenceIndex : [annotated_sequence, [other, values, to, be, saved] for each duplicate]}
	"""

	def getDuplicates():
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
		# todo outdated documentation

		duplicateLists = []  # [list of [list of duplicate indices] for each duplicate]

		def groupByIdenticalProperties(byPropDict, remainingProperties):
			# todo: if the code inside this function doesnt work, use the one outside this function instead
			if remainingProperties:
				for prop, byPropIndices in byPropDict:
					if len(byPropIndices) > 1:  # only if there are duplicates
						## SELECT IDENTICAL <NEXTPROPERTY> ##
						groupByIdenticalProperties(df[byPropIndices].groupby(remainingProperties[0]),
						                           remainingProperties[1:])
			else:
				duplicateLists.extend(byPropDict.values)
			return duplicateLists

		youreFeelingLucky = True  # todo: disable this if the code above doesnt work (TRIGGERS CODE IN FUNCTION ABOVE)
		if youreFeelingLucky:
			properties = []
			if not undoublePSMAlgo_bool:  # only if you didn't undoublePSMAlgo
				## SELECT IDENTICAL PSMALGO (i.e. different First Scan) ##
				byFirstPropDict = df.groupby('Identifying Node').groups
				properties.append('Annotated Sequence')
			else:
				## SELECT IDENTICAL SEQUENCE ##
				byFirstPropDict = df.groupby('Annotated Sequence').groups
			if toCollapse == 'RT':
				groupByIdenticalProperties(byFirstPropDict, properties + ['Charge', 'Modifications'])
				return duplicateLists
			elif toCollapse == 'Charge':
				groupByIdenticalProperties(byFirstPropDict, properties + ['Modifications'])
			elif toCollapse == 'PTM':
				groupByIdenticalProperties(byFirstPropDict, properties + ['Charge'])
			return duplicateLists

		elif not youreFeelingLucky:
			## SELECT IDENTICAL SEQUENCE ##
			bySequenceDict = df.groupby('Annotated Sequence').groups
			if toCollapse == 'RT':
				for sequence, bySequenceIndices in bySequenceDict:
					if len(bySequenceIndices) > 1:  # only if there are duplicates
						## SELECT IDENTICAL CHARGE ##
						byChargeBySequenceDict = df[bySequenceIndices].groupby('Charge').groups
						for charge, byChargeIndices in byChargeBySequenceDict:
							if len(byChargeIndices) > 1:  # only if there are duplicates
								## SELECT IDENTICAL PTM ##
								byPTMByChargeBySequenceDict = df[byChargeIndices].groupby('Modifications').groups
								if not undoublePSMAlgo_bool:  # only if you didn't undoublePSMAlgo
									for PTM, byPTMIndices in byPTMByChargeBySequenceDict:
										if len(byPTMIndices) > 1:  # only if there are duplicates
											## SELECT IDENTICAL PSMALGO (i.e. different First Scan) ##
											byPSMAlgoByPTMByChargeBySequenceDict = df[byPTMIndices].groupby(
												'Identifying Node').groups
											duplicateLists.extend(byPSMAlgoByPTMByChargeBySequenceDict.values)
								else:
									duplicateLists.extend(byPTMByChargeBySequenceDict.values)

			elif toCollapse == 'Charge':
				for sequence, bySequenceIndices in bySequenceDict:
					if len(bySequenceIndices) > 1:  # only if there are duplicates
						## SELECT IDENTICAL PTM ##
						byPTMBySequenceDict = df[bySequenceIndices].groupby('Modifications').groups
						if not undoublePSMAlgo_bool:  # only if you didn't undoublePSMAlgo
							for PTM, byPTMIndices in byPTMBySequenceDict:
								if len(byPTMIndices) > 1:  # only if there are duplicates
									## SELECT IDENTICAL PSMALGO (i.e. different First Scan) ##
									byPSMAlgoByPTMBySequenceDict = df[byPTMIndices].groupby('Identifying Node').groups
									duplicateLists.extend(byPSMAlgoByPTMBySequenceDict.values)
						else:  # you did undoublePSMAlgo? Great, you're done.
							duplicateLists.extend(byPTMBySequenceDict.values)
						## SANITY CHECK ##
						for PTM, byPTMIndices in byPTMBySequenceDict:  # TEST
							if len(byPTMIndices) > 1:  # only if there are duplicates
								## SELECT IDENTICAL CHARGE ##
								byChargeByPTMBySequenceDict = df[byPTMIndices].groupby('Charge').groups
								for charge, byChargeIndices in byChargeByPTMBySequenceDict:
									assert len(
										byChargeIndices) < 2  # if same Sequence and same PTM, Charge cannot be the same because it would have been RT-collapsed.

			elif toCollapse == 'PTM':
				for sequence, bySequenceIndices in bySequenceDict:
					if len(bySequenceIndices) > 1:  # only if there are duplicates
						## SELECT IDENTICAL CHARGE ##
						byChargeBySequenceDict = df[bySequenceIndices].groupby('Charge').groups
						if not undoublePSMAlgo_bool:  # only if you didn't undoublePSMAlgo
							for charge, byChargeIndices in byChargeBySequenceDict:
								if len(byChargeIndices) > 1:  # only if there are duplicates
									## SELECT IDENTICAL PSMALGO ##
									byPSMAlgoByChargeBySequenceDict = df[byChargeIndices].groupby(
										'Identifying Node').groups
									duplicateLists.extend(byPSMAlgoByChargeBySequenceDict.values)
						else:  # you did undoublePSMAlgo? Great, you're done.
							duplicateLists.extend(byChargeBySequenceDict.values)
						## SANITY CHECK ##
						for charge, byChargeIndices in byChargeBySequenceDict:  # TEST
							if len(byChargeIndices) > 1:  # only if there are duplicates
								## SELECT IDENTICAL PTM ##
								byPTMByChargeBySequenceDict = df[byChargeIndices].groupby('Modifications').groups
								for PTM, byPTMIndices in byPTMByChargeBySequenceDict:
									assert len(
										byPTMIndices) < 2  # if same Sequence and same Charge, PTM cannot be the same because it would have been RT-collapsed.
			return duplicateLists

	def combineDetections(duplicateLists, centerMeasure):
		for duplicatesList in duplicateLists:
			# calculate the total MS2 intensities for each duplicate
			allMS2Intensities = np.asarray(df.loc[duplicatesList][intensityColumns])

		flagMaxReporterVariance = False  # TODO flag when maxReporterVariance is exceeded
		if flagMaxReporterVariance: # TODO flag when maxReporterVariance is exceeded
			# this can only be consistent if executed on RELATIVE intensities.
			Ri = 1 / allMS2Intensities.shape[1] * np.asarray(1 / np.nanmean(allMS2Intensities, 1)).reshape(allMS2Intensities.shape[0], )
			relativeIntensities = (allMS2Intensities.T * Ri).T
			if np.any(np.var(relativeIntensities,axis=0) > maxRelativeReporterVariance):
				# TODO this shouldnt just warn, you should also decide what to do.
				warn("maxRelativeReporterVariance too high for duplicates with indices: " + str(duplicatesList) + ".")

		if centerMeasure == 'mean':
			pass
		if centerMeasure == 'geometricMedian':
			pass
		if centerMeasure == 'weighted':
			pass
		return newIntensities  # TODO

	def getNewIntensities(duplicatesDf, duplicatesDict, masterPSMAlgo):
		"""
		Combines the true duplicates' intensities into one new entry per first occurrence, conform the duplicatesDict structure.
		:param duplicatesDict:          dict            {firstOccurrenceIndex:[duplicateIndices]}
		:param duplicatesDf:            pd.dataFrame    data of only the first occurrences and duplicates
		:return newIntensitiesDict:     dict            {firstOccurrenceIndex:np.array(newIntensities)}
		"""

		return detection  # TODO

	def getBestIndices(duplicateLists):
		# get the detection with the best PSM match
		# values of BEST PSM detection in all duplicates (master/slave-wise best)
		# dont forget to increase Degeneracy # todo
		bestIndices = []
		if masterPSMAlgo == 'mascot':
			for duplicatesList in duplicateLists:
				bestIndex = df['Ions Score'].loc[duplicatesList].idxmax(axis=0,skipna=True)
				if np.isnan(bestIndex): # no Mascot scores found --> take best Sequest
					bestIndex = df['XCorr'].loc[duplicatesList].idxmax(axis=0, skipna=True)
					assert not np.isnan(bestIndex)
					bestIndices.append(bestIndex)
		elif masterPSMAlgo == 'sequest':
			for duplicatesList in duplicateLists:
				bestIndex = df['XCorr'].loc[duplicatesList].idxmax(axis=0,skipna=True)
				if np.isnan(bestIndex): # no Sequest scores found --> take best Mascot
					bestIndex = df['Ions Score'].loc[duplicatesList].idxmax(axis=0, skipna=True)
					assert not np.isnan(bestIndex)
					bestIndices.append(bestIndex)
		return bestIndices

	def getIntenseIndices(duplicateLists):
		intenseIndices = []
		for duplicatesList in duplicateLists:
			# calculate the total MS2 intensities for each duplicate
			totalIntensities = np.sum(np.asarray(df.loc[duplicatesList][intensityColumns]),axis=1)
			# get the most intense duplicate
			intenseIndex = duplicatesList[np.argmax(totalIntensities)]
			assert not np.isnan(intenseIndex)
			intenseIndices.append(intenseIndex)
		return intenseIndices

	def getRepresentativesDf(bestIndices):
		representativesDf = df.loc[bestIndices]
		if method == 'bestMatch':
			pass
		elif method == 'mostIntense':
			intenseIndices = getIntenseIndices(duplicateLists)
			# generate { bestIndex : [mostIntense intensities] }
			intensitiesDict = dict(zip(list(bestIndices), getIntensities(df.loc[intenseIndices])))
			# set the representative intensities to be the most intense intensities
			representativesDf = setIntensities(representativesDf, intensitiesDict)
		else: # method == 'centerMeasure'
			newIntensities = combineDetections(duplicateLists, centerMeasure=method)
			# generate { intenseIndex : [intensities] }
			intensitiesDict = dict(zip(list(bestIndices), newIntensities))
			# set the representative intensities to be the most intense intensities
			representativesDf = setIntensities(representativesDf, intensitiesDict)

		# reindex representativesDf so it can be concatenated properly with new indices
		representativesDf.index = list(range(df.index[-1], df.index[-1] + len(representativesDf.index)))
		return representativesDf

		i = len(df.index)
		for representative in representatives:
			df.loc[i] = representative
			i += 1

	if 'Degeneracy' not in df.columns:
		# contains the number of peptides that have been collapsed onto each (synthetic) detection.
		df['Degeneracy'] = [1,]*len(df.index)

	# get a nested list of duplicates according to toCollapse. [[duplicates1], [duplicates2], ...]
	duplicateLists = getDuplicates()
	if method == 'RT':
		pass # TODO flag isolated RT peaks
	elif method == 'PTM':
		pass # TODO flag PTM differences.
	# get the new intensities per first occurrence index (df index)
	bestIndices = getBestIndices(duplicateLists) # {bestIndex : [duplicates]}
	# add the new representative detections to the dataFrame
	representativesDf = getRepresentativesDf(bestIndices)
	df = df.append(representativesDf)
	toDelete = [item for sublist in duplicateLists for item in sublist] # unpack list of lists
	# save as {representativeIndex : df[for each duplicate,[values, to, be, saved]]}
	duplicatesDfList = [df.loc[duplicatesList,columnsToSave] for duplicatesList in duplicateLists]
	removedData = dict(zip(representativesDf.index, duplicatesDfList))
	df.drop(toDelete, inplace=True)

	return df, removedData
