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
		Takes the dataFrame df and returns a nested list of duplicates per category according to the toCollapse variable.
		Based on iterative use of the pd.DataFrame.groupby('property').groups function which returns a dict
		{ propertyValue : [duplicateIndices] }.
		:return duplicateLists:     list            [[group of duplicates] per toCollapse value in the df]
		"""
		# todo outdated documentation

		duplicateLists = []  # [list of [list of duplicate indices] for each duplicate]

		def groupByIdenticalProperties(byPropDict, remainingProperties):
			"""
			Takes a dictionary which is the result of a groupby(property) call on a dataFrame. Also takes a list of
			properties, and uses the first one to do another groupby() on each dataframe consisting of one set of
			duplicates in the input dictionary values. Then, iteratively calls itself again using that dictionary and
			the other remaining properties. The result is a nested list of duplicates which all have identical
			combination-of-properties values, but non-identical toCollapse values.
			For instance when toCollapse=='RT':
			[
				[2, 4, 8], # indices of detections with (Charge1, PTM2) and unique(RT2, RT4, RT8)==True
				[1, 5, 6], # indices of detections with (Charge1, PTM1) and unique(RT1, RT5, RT6)==True
				[3],       # indices of detections with (Charge3, PTM3) and unique(RT3)==True
				[7, 9],    # indices of detections with (Charge3, PTM2) and unique(RT7, RT9)==True
			]
			This function correctly groups by PSMAlgo when required and does not when it is prohibited.
			:param byPropDict:          dict    { propertyValue : [duplicateIndices] }
			:param remainingProperties: list    properties still to be grouped by
			:return duplicateLists:     list    [[group of duplicates] per combination-of-properties values in the dataFrame]
			"""
			# TODO: if the code inside this function doesnt work, use the one outside this function instead
			if remainingProperties:
				for propValue, byPropIndices in byPropDict.items():
					if len(byPropIndices) > 1: # only if there are duplicates
						## SELECT IDENTICAL <NEXTPROPERTY> ##
						groupByIdenticalProperties(df.loc[byPropIndices].groupby(remainingProperties[0]).groups,
						                           remainingProperties[1:]) # first pop the [0] property to both return and remove it!
			else:
				duplicateLists.extend(byPropDict.values())

		youreFeelingLucky = True  # TODO: disable this if the code above doesnt work (TRIGGERS CODE IN FUNCTION ABOVE)
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
		"""
		Takes a nested list of the indices of duplicates and combines the rows of the intensity matrix -- found in
		dataFrame df --	for each sublist into one new row of intensities for the representative that is to be their
		replacement. This function should also flag cases where the variance between the intensities (calculated per
		reporter channel) exceeds a maxRelativeReporterVariance.
		:param duplicateLists:  list        [[group of duplicates] per toCollapse value in the df]
		:param centerMeasure:   str         specifies the method of combination
		:return newIntensities: np.ndarray  new intensities of the representative detection
		"""
		flagMaxRelativeReporterVariance = False
		for duplicatesList in duplicateLists:
			# calculate the total MS2 intensities for each duplicate
			allMS2Intensities = np.asarray(df.loc[duplicatesList][intensityColumns])

			if flagMaxRelativeReporterVariance: # TODO flag when maxRelativeReporterVariance is exceeded
				# this can only be consistent if executed on RELATIVE intensities.
				Ri = 1 / allMS2Intensities.shape[1] * np.asarray(1 / np.nanmean(allMS2Intensities, 1)).reshape(allMS2Intensities.shape[0], )
				relativeIntensities = (allMS2Intensities.T * Ri).T
				if np.any(np.var(relativeIntensities,axis=0) > maxRelativeReporterVariance):
					# TODO this shouldnt just warn, you should also decide what to do.
					warn("maxRelativeReporterVariance too high for duplicates with indices: " + str(duplicatesList) + ".")

			if centerMeasure == 'mean':
				newIntensities = np.mean(allMS2Intensities, 0)
			elif centerMeasure == 'geometricMedian': # TODO
				pass
			elif centerMeasure == 'weighted': # TODO
				pass
		return newIntensities

	def getBestIndices(duplicateLists):
		"""
		For each sublist in the nested list of duplicates duplicateLists, calculates the index of the duplicate with the
		best PSM match according to dataFrame df. Does this intelligently by taking masterPSMAlgo into account.
		:param duplicateLists:  list        [[group of duplicates] per toCollapse value in the df]
		:return bestIndices:    list        [indices of detections with the best PSM score per group of duplicates]
		"""
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
		"""
		For each sublist in the nested list of duplicates duplicateLists, calculates the total MS2 intensity according
		to dataFrame df and returns the results as a list.
		:param duplicateLists:      list        [[group of duplicates] per toCollapse value in the df]
		:return intenseIndices:     list        [indices of detections with the highest total MS2 intensity per group of duplicates]
		"""
		intenseIndices = []
		for duplicatesList in duplicateLists:
			# calculate the total MS2 intensities for each duplicate
			totalIntensities = np.sum(np.asarray(df.loc[duplicatesList][intensityColumns]),axis=1)
			# get the most intense duplicate
			intenseIndex = duplicatesList[np.argmax(totalIntensities)]
			assert not np.isnan(intenseIndex)
			intenseIndices.append(intenseIndex)
		return intenseIndices

	def getRepresentativesDf(bestIndices, duplicateLists):
		"""
		Uses a list of indices of the best PSM matches bestIndices amongst each group of duplicates in the nested list
		duplicateLists, all indices with respect to dataFrame df. Based on this best PSM match, generates a
		representative detection for each group of duplicates. This is done by copying all bestMatch properties, but by
		calculating new intensities when necessary and also updating the Degeneracy parameter.
		:param bestIndices:         list            indices of the best PSM matches inside a group of duplicates (see duplicatLists)
		:param duplicateLists:      list            [[group of duplicates] per toCollapse value in the df]
		:return representativesDf:  pd.dataFrame    all representatives data that will replace the duplicate entries in the dataFrame df
		"""
		representativesDf = df.loc[bestIndices]
		# sum the degeneracies of all duplicates involved in each representative
		representativesDf['Degeneracy'] = [np.sum(np.asarray(representativesDf.loc[duplicatesList, 'Degeneracy']))
		                                   for duplicatesList in duplicateLists]
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
	representativesDf = getRepresentativesDf(bestIndices, duplicateLists)
	df = df.append(representativesDf)
	toDelete = [item for sublist in duplicateLists for item in sublist] # unpack list of lists
	# save as {representativeIndex : df[for each duplicate,[values, to, be, saved]]}
	duplicatesDfList = [df.loc[duplicatesList,columnsToSave] for duplicatesList in duplicateLists]
	removedData = dict(zip(representativesDf.index, duplicatesDfList))
	df.drop(toDelete, inplace=True)

	return df, removedData
