#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Collection of functions that help removing redundancy in the data due to detections which have:
* different retention time (RT) values
* different charges
* different (post-translational) modifications (PTMs)
and replaces the duplicates with one representative detection and a combination/summary/selection of their intensities.
"""

import numpy as np
import logging
from warnings import filterwarnings
from constandpp.tools import setIntensities, getIntensities
from scipy.spatial.distance import cdist, euclidean

columnsToSave = None


def geometricMedian(X, eps=1e-5):
	"""
	Returns the geometric median conform to the algorithm developed by Yehuda Vardi and Cun-Hui Zhang.
	(doi: 10.1073/pnas.97.4.1423) and implemented by http://stackoverflow.com/users/565635/orlp.
	:param X:   np.ndarray      multidimensional vectors
	:param eps: float64         precision
	:return:    np.ndarray      geometric median multidimensional vector
	"""
	# throw away detections with nan-values
	X = X[~np.isnan(X).any(axis=1)]

	y = np.mean(X, 0)
	while True: # as long as precision eps not reached
		D = cdist(X, [y])
		nonzeros = (D != 0)[:, 0]

		Dinv = 1 / D[nonzeros]
		Dinvs = np.sum(Dinv)
		W = Dinv / Dinvs
		T = np.sum(W * X[nonzeros], 0)

		num_zeros = len(X) - np.sum(nonzeros)
		if num_zeros == 0:
			y1 = T
		elif num_zeros == len(X):
			return y
		else:
			R = (T - y) * Dinvs
			r = np.linalg.norm(R)
			rinv = 0 if r == 0 else num_zeros/r
			y1 = max(0, 1-rinv)*T + min(1, rinv)*y

		if euclidean(y, y1) < eps:
			return y1
		y = y1


def collapse(toCollapse, df, quanColumns, method, identifyingNodes, undoublePSMAlgo_bool, columnsToSave):  #
	"""
	Generic collapse function, which removes redundancy in the data due to property toCollapse.
	Looks for duplicate 'Annotated Sequence' values in the dataFrame and further groups by other possible properties
	(PSMAlgo, Charge, Modifications) if necessary, while respecting the collapse order imposed in the main() function.
	Removes each group of duplicate from the experimental data df and replaces them by a representative entry that has
	all properties of the entry amongst them with the best PSM score, but with quantification values determined by a
	method specified through `method`.
	intensities (via getNewIntensities function): remove all duplicates and enter one replacement detection.
	Adds a 'Degeneracy' column to the dataFrame if it didn't exist already: this contains the number of peptides that
	have been collapsed onto that (synthetic) detection.
	Returns removedData according to the columnsToSave list.
	:param toCollapse:              	str             variable of which true duplicates are to be collapsed.
	:param df:                          pd.dataFrame    with sequence duplicates due to difference in certain variables/columns.
	:param quanColumns:			list			columns that contain the quantification values
	:param method:                      str             defines how the new quantification values of the representative
														are selected/constructed:
														bestMatch: those of the PSM with the best score
														mostIntsene: those of the PSM with the most intense values
														mean: the mean of all PSMs in that duplicates group
														geametricMedian: the geometric median of all PSMs in that duplicates group
	:param identifyingNodes:			dict			PSM rater algorithm names and score names. Structure:
														{"master": [NAME, SCORE_NAME], "slaves": [[NAME, SCORE_NAME], ...]}
	:param undoublePSMAlgo_bool:		bool			Has the PSM algo redundancy already been removed?
	:param columnsToSave:				list(str)		fields of to-be-removed entries that are saved into removedData
	:return df:                         pd.dataFrame    without sequence duplicates according to to checkTrueDuplicates.
	:return removedData:                pd.dataFrame    [PARENT INDEX, and, fields, to, be, saved]
	"""
	# todo docu
	def getDuplicates():
		"""
		Takes the dataFrame df from the parent scope and returns a nested list of groups of duplicates according to
		the toCollapse variable. Based on iterative use of the pd.DataFrame.groupby('property').groups function which
		returns a dict { propertyValue : [duplicateIndices] }.
		:return this_duplicateLists:     list            [[group of duplicates] per `toCollapse` value in the df]
		"""
		def groupByIdenticalProperties(byPropDict, remainingProperties):
			"""
			Takes a dictionary which is the result of a groupby(property) call on a dataFrame. Also takes a list of
			properties, and uses the first one to do another groupby() on each dataframe slice consisting of one set of
			duplicates in the input dictionary values. Then, iteratively calls itself again using that dictionary and
			the other remaining properties. The result is a nested list of duplicates which all have identical
			combination-of-properties values, but non-identical toCollapse values.
			For instance when toCollapse=='RT':
			[
				[2, 4, 8], # indices of detections with (Charge2, PTM2) and unique(RT2, RT4, RT8)==True
				[1, 5, 6], # indices of detections with (Charge1, PTM1) and unique(RT1, RT5, RT6)==True
				[3],       # indices of detections with (Charge3, PTM3) and unique(RT3)==True
				[7, 9],    # indices of detections with (Charge7, PTM7) and unique(RT7, RT9)==True
			]
			This function correctly groups by PSMAlgo when required and does not when it is prohibited.
			:param byPropDict:          	dict    { propertyValue : [duplicateIndices] }
			:param remainingProperties: 	list	properties still to be grouped by
			:return this_duplicateLists:	list    [[group of duplicates] per combination-of-properties values in the dataFrame]
			"""
			# use only indices that are not single (i.e. that have a duplicate)
			notSingleList = list(filter(lambda e: len(e) > 1, byPropDict.values()))  # [list of [lists of length > 1]]
			if remainingProperties:  # check for more identical properties before marking as duplicates
				for byPropIndices in notSingleList: # only if there actually is at least one group of dnp.asarray(df.locuplicates
					## SELECT IDENTICAL <NEXTPROPERTY> ##  # for each sublist of notSingleList
					groupByIdenticalProperties(df.loc[byPropIndices].groupby(remainingProperties[0]).groups,
											   remainingProperties[1:])  # first pop the [0] property to both return and remove it!
			else:  # no more properties to check: mark groups of indices as duplicates
				this_duplicateLists.extend(notSingleList)

		this_duplicateLists = []  # [list of [list of duplicate indices] for each duplicate]
		# list of properties
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
		elif toCollapse == 'Charge':
			groupByIdenticalProperties(byFirstPropDict, properties + ['Modifications'])
		elif toCollapse == 'PTM':
			# modifications are apparent from the sequence! Remove this dependence!
			byFirstPropDict = df.groupby(df['Annotated Sequence'].str.upper()).groups
			groupByIdenticalProperties(byFirstPropDict, properties + ['Charge'])

		return this_duplicateLists

	def combineDetections(this2_bestIndicesDict, centerMeasure):
		"""
		Takes a dict of lists of duplicate indices and combines for each list the corresponding rows of the
		quantification matrix -- found in dataFrame df (parent scope) -- into one new row of intensities for the
		representative entry that is to be their replacement in the complete dataframe. The method of combining
		quantification values is specified by centerMeasure.
		:param this2_bestIndicesDict:   dict    { bestIndex : [group of duplicate indices]}
		:param centerMeasure:           str     specifies the method of combination
		:return newIntensitiesDict:     dict    { bestIndex : new intensities of the representative detection }
		"""
		newIntensitiesDict = {}
		for bestIndex, this_duplicatesList in this2_bestIndicesDict.items():
			# calculate the total MS2 intensities for each duplicate
			allMS2Intensities = getIntensities(df, quanColumns=quanColumns, indices = this_duplicatesList)

			if allMS2Intensities is None: # TEST
				print('hoi')

			if centerMeasure == 'mean':
				filterwarnings('ignore', 'Mean of empty slice') # catch warnings about empty slices
				newIntensitiesDict[bestIndex] = np.nanmean(allMS2Intensities, 0)
			elif centerMeasure == 'geometricMedian':
				# first normalize the row sums because otherwise the Median norm isn't conserved. (set it to 1 now)
				Ri = 1 / allMS2Intensities.shape[1] * np.asarray(1 / np.nanmean(allMS2Intensities, 1)).reshape(
					allMS2Intensities.shape[0], )
				relativeIntensities = (allMS2Intensities.T * Ri).T
				newIntensitiesDict[bestIndex] = geometricMedian(relativeIntensities)
			elif centerMeasure == 'weighted':  # TODO
				pass
		return newIntensitiesDict

	def getBestIndicesDict(this_duplicateLists):
		"""
		For each sublist in the nested list of duplicates duplicateLists, calculates the index of the duplicate with the
		best PSM match according to dataFrame df. Does this intelligently by taking masterPSMAlgo into account. If no
		best index can be found (due to missing score for instance) it just takes the first in this_duplicateLists.
		:param this_duplicateLists: list    [[group of duplicates] per toCollapse value in the df]
		:return this_bestIndices:   dict    { indices of detections with the best PSM score per group of duplicates : [group of duplicates] }
		"""
		isNanWarnedYet = False
		noSlavePSMAlgoWarnedYet = False
		this_bestIndicesDict = {}
		masterScoreName = identifyingNodes['master'][1]
		if len(identifyingNodes['slaves']) > 0:  # if there is at least one PSMAlgo Slave
			slaveScoreName = identifyingNodes['slaves'][0][1]
		else:  # if there are no PSMAlgo Slaves
			slaveScoreName = None

		for this_duplicatesList in this_duplicateLists:
			bestIndex = df.loc[this_duplicatesList, masterScoreName].idxmax(axis=0, skipna=True)
			if np.isnan(bestIndex) and slaveScoreName is not None:  # no MASTER scores found --> take best SLAVE (if it exists)
				try:
					bestIndex = df.loc[this_duplicatesList, slaveScoreName].idxmax(axis=0, skipna=True)
				except KeyError:  # if no slave score column is present in the data set
					if not noSlavePSMAlgoWarnedYet:  # I haven't thrown this warning yet
						logging.warning("No slave PSMAlgo score column ('"+slaveScoreName+"') present in data set. ")
						noSlavePSMAlgoWarnedYet = True
				if np.isnan(bestIndex) and not isNanWarnedYet:
					logging.warning("No best PSM score found for some lists of duplicates; first duplicate arbitrarily chosen. "
						 "First Scan numbers of first list encountered: "+str(df.loc[this_duplicatesList, 'First Scan']))
					isNanWarnedYet = True
					bestIndex = this_duplicatesList[0]
			this_bestIndicesDict[bestIndex] = this_duplicatesList

		return this_bestIndicesDict

	def getIntenseIndicesDict(this_bestIndicesDict):
		"""
		Takes (the indices of) groups of duplicates associated with the index of their 'best' member through a dict, and
		returns an analogous dict where the groups are now replaced by the index of the 'most intense' entry, according to
		the total sum of each entry's quantification values.
		:param this_bestIndicesDict:	dict	{ indices of best PSM entry per group of duplicates : [group of duplicates] }
		:return intenseIndicesDict:		dict	{ best PSM entry index per group of duplicates : indices of most intense entry }
		"""
		intenseIndicesDict = {}
		for bestIndex, this_duplicatesList in this_bestIndicesDict.items():
			# calculate the total MS2 intensities for each duplicate
			totalIntensities = np.sum(getIntensities(df, this_duplicatesList), axis=1)
			# get the most intense duplicate
			intenseIndex = this_duplicatesList[np.argmax(totalIntensities)]
			assert not np.isnan(intenseIndex)
			intenseIndicesDict[bestIndex] = intenseIndex
		return intenseIndicesDict

	def getRepresentativesDf(this_bestIndicesDict):
		"""
		Uses a list of indices of the best PSM matches bestIndices amongst each group of duplicates. Based on this best
		PSM entry, generates a representative detection for each group of duplicates. This is done by copying all
		bestMatch properties, while calculating new quantification values when necessary and also updating the Degeneracy
		parameter. Indices are taken w.r.t. df from the parent scope.
		The behaviour should not be overwritten so as to use the mostIntense indices instead, since the quantification
		values or not of importance here, but rather the PSM properties and their reliability.
		:param this_bestIndicesDict:     dict           { index of detection with best PSM score per group of duplicates : [group of duplicates] }
		:return this_representativesDf:  pd.dataFrame   all representatives data that will replace the duplicate entries in the dataFrame df
		"""
		this_bestIndices = this_bestIndicesDict.keys()
		this_duplicateLists = this_bestIndicesDict.values()

		this_representativesDf = df.loc[this_bestIndices, :].copy(deep=True)
		# sum the degeneracies of all duplicates involved in each representative
		this_representativesDf.loc[:, 'Degeneracy'] = [np.sum(np.asarray(df.loc[(this_duplicatesList, 'Degeneracy')])) for this_duplicatesList in this_duplicateLists]

		if method == 'bestMatch':
			pass
		elif method == 'mostIntense':
			intenseIndicesDict = getIntenseIndicesDict(this_bestIndicesDict)
			# generate { bestIndex : [mostIntense intensities] }
			intensitiesDict = dict((ind_best, getIntensities(df.loc[ind_intense, :])) for (ind_best, ind_intense) in intenseIndicesDict.items())
			# set the representative intensities to be the most intense intensities
			this_representativesDf = setIntensities(this_representativesDf, quanColumns=quanColumns, intensities=intensitiesDict)
		else:  # method == 'centerMeasure'
			newIntensitiesDict = combineDetections(this_bestIndicesDict, centerMeasure=method)
			# set the representative intensities to be the most intense intensities
			this_representativesDf = setIntensities(this_representativesDf, quanColumns=quanColumns, intensities=newIntensitiesDict)

		# reindex this_representativesDf so it can be concatenated properly with new indices
		this_representativesDf.index = list(range(max(df.index)+1, max(df.index)+1 + len(this_representativesDf.index)))
		# representative indices do not correspond to these indices!!!!!
		return this_representativesDf

	if 'Degeneracy' not in df.columns:
		# contains the number of peptides that have been collapsed onto each (synthetic) detection.
		df.loc[:, 'Degeneracy'] = [1, ] * len(df.index)

	# get a nested list of duplicates according to toCollapse. [[duplicates1], [duplicates2], ...]
	duplicateLists = getDuplicates()
	# get the new intensities per first occurrence index (df index)
	bestIndicesDict = getBestIndicesDict(duplicateLists)  # {bestIndex : [duplicates]}
	# add the new representative detections to the dataFrame
	representativesDf = getRepresentativesDf(bestIndicesDict)
	df = df.append(representativesDf)
	toDelete = [item for sublist in duplicateLists for item in sublist]  # unpack list of lists
	removedData = df.loc[toDelete, columnsToSave]
	# add the representative index of each collection of collapsed duplicates
	removedData.insert(loc=0, column='Representative First Scan', value=-1)
	for duplicatesList, rfs in zip(bestIndicesDict.values(),
								   representativesDf['First Scan']):  # relies on fact that order is conserved! #todo
		removedData.loc[duplicatesList, 'Representative First Scan'] = rfs
	# actually remove the toDelete detections
	df.drop(toDelete, inplace=True)

	return df, removedData
