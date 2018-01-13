#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Collection of functions involved in analyzing the data that was processed by processing.py and constand.py.
Performs a differential expression analysis on the normalized intensities as provided by CONSTANd.
"""

import numpy as np
import pandas as pd
import logging
from qcquan.tools import unnest, getOtherConditions
from collections import defaultdict
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import ttest_ind as ttest
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage


# @profile  # kernprof / line_profiler decorator
def getRTIsolationInfo(removedData_RT):
	"""
	Returns dataFrame with the mean, standard deviation, and max-min value of the RT values for each duplicate_group
	representative that is found in the removedData_RT dataframe containing the removedData for the RT aggregate.
	:param removedData_RT:  pd.DataFrame    removedData for the RT aggregate.
	:return:                pd.DataFrame    statistics 'Degeneracy', 'mean', 'std', 'max-min' about the RT values
	"""
	if 'Representative First Scan' in removedData_RT.columns.values and 'RT [min]' in removedData_RT.columns.values:
		duplicateGroups = removedData_RT.groupby('Representative First Scan').groups
		RTIsolationInfo = []
		for rfs, duplicates in duplicateGroups.items():
			RTValues = removedData_RT.loc[duplicates, 'RT [min]']
			RTIsolationInfo.append([rfs, len(duplicates), np.nanmean(RTValues), np.std(RTValues), np.ptp(RTValues)])
		return pd.DataFrame(RTIsolationInfo, columns=['Representative First Scan', 'Degeneracy', 'mean', 'std', 'max-min'])
	else:
		logging.warning("No RT isolation info could be calculated because a necessary column was missing.")
		return pd.DataFrame()


def getNoIsotopicCorrection(df, noCorrectionIndices):
	"""
	Given a dataframe and indices of PSMs that received no corrections, returns some basic info about them.
	:param df:                  pd.DataFrame	processed and possibly (if not: exception caught) isotope-corrected data frame
	:param noCorrectionIndices: list            indices of PSMs that received no isotopic correction
	:return:                    pd.DataFrame    ['First Scan', 'Identifying Node Type', 'Sequence', 'Master Protein Accessions']
	"""
	try:
		return df.loc[noCorrectionIndices, ['First Scan', 'Identifying Node Type', 'Sequence', 'Master Protein Accessions']]
	except KeyError:
		logging.warning("No NoIsotopicCorrection info could be saved because a necessary column was missing.")
		return pd.DataFrame()


def combineProcessingMetadata(metadata, perExperimentMetadata):
	""" The metadata gathered during the Processing step is still spread over multiple dict entries.
	Combine it where necessary. """
	experimentNames = perExperimentMetadata.keys()
	# Compile a list of all master proteins found at least in 1 PSM and at least in 1 experiment:
	allObservedProteins = pd.Series(list(set().union(*[perExperimentMetadata[eName]['allMasterProteins'] for eName in experimentNames])))
	# info on amount of PSMs
	metadata['numPSMs'] = pd.DataFrame(index=experimentNames, columns=['initial', 'after cleaning'])
	metadata['pctPSMsIsolInterfTooHigh'] = pd.DataFrame(columns=experimentNames)
	metadata['MS1Intensities_PSMs'] = pd.Series(index=experimentNames, dtype=object)
	metadata['injectionTimeInfo'] = pd.DataFrame(index=experimentNames, columns=['max', 'num max', 'num below'])
	metadata['deltappmStatistics'] = pd.DataFrame(index=experimentNames, columns=['max', 'mean', 'std'])
	metadata['intensityStatisticsPerExp'] = dict()  # dict because we need a whole dataframe per experiment
	metadata['relPSMScoreVsDeltaMppmPerExp'] = dict()  # dict because we need a whole dataframe per experiment
	for eName in experimentNames:
		metadata['numPSMs'].loc[eName, :] = [perExperimentMetadata[eName]['numPSMs_initial'],
											 perExperimentMetadata[eName]['numPSMs_afterCleaning']]
		metadata['pctPSMsIsolInterfTooHigh'].loc[0, eName] = perExperimentMetadata[eName]['pctPSMsIsolInterfTooHigh']
		try:
			metadata['MS1Intensities_PSMs'][eName] = perExperimentMetadata[eName]['MS1Intensities_PSMs'].tolist()
		except KeyError as e:
			logging.warning("Entry '" + str(e.args[0]) + "' was not found for experiment "+eName+". Not gathering ANY MS1 intensity QC info.")
			if 'MS1Intensities_PSMs' in metadata.keys():
				del metadata['MS1Intensities_PSMs']
		metadata['injectionTimeInfo'].loc[eName, :] = perExperimentMetadata[eName]['injectionTimeInfo'].iloc[0, :]  # there is only 1 entry
		metadata['deltappmStatistics'].loc[eName, :] = perExperimentMetadata[eName]['deltappmStatistics'].iloc[0, :]  # there is only 1 entry
		metadata['intensityStatisticsPerExp'][eName] = perExperimentMetadata[eName]['intensityStatistics']
		metadata['relPSMScoreVsDeltaMppmPerExp'][eName] = perExperimentMetadata[eName]['relPSMScoreVsDeltaMppm']
	
	return metadata, allObservedProteins


def combineExperimentDFs(dfs):
	"""
	Merge dataframes of all experiments into one multi-indexed (eName, oldIndex) dataframe, by performing an outer join.
	The intensity columns are non-identical across different dataframes and resulting empty fields are valued NaN.
	:param dfs:     [pd.DataFrame]	dataframes from multiple experiments
	:return: 		pd.DataFrame	dataframe containing an outer join of the input list of dataframes
	"""
	# how are PSMs combined with multiple charge states for instance?
	# since we are using keys= there will be NO merging of entries, because each peptide will have its own unique index
	# index = (eName, oldIndex)
	return pd.concat(dfs.values(), keys=dfs.keys(), join='outer')


def getProteinPeptidesDicts(df, fullExpression_bool):
	"""
	Returns two dicts with the peptide indices (w.r.t. dataframe df) associated with each protein in the df as a
	dictionary. One dict (min) contains only the peptide indices uniquely associated per protein, the other contains
	all peptides indices associated per protein.
	:param df:							pd.DataFrame	peptide-level data frame
	:param fullExpression_bool:			bool			also get non-injective protein-peptide relations dict
	:return minProteinPeptidesDict:		dict			{ protein : uniquely associated peptide indices (injective) }
	:return fullProteinPeptidesDict:	dict			{ protein : all associated peptide indices (non-injective) }
	:return 							pd.Series		peptide sequences without master protein accession
	"""
	# create this column if it doesn't exist. It gets removed afterwards anyway.
	if "# Protein Groups" not in df.columns.values:
		df["# Protein Groups"] = df.apply(lambda x: len(str(x['Master Protein Accessions']).split(';')), axis=1)
	# its OK that # Protein Groups is inferred only from info within the same experiment. This is how it should be.
	numProteinGroupsDict = df.groupby("# Protein Groups").groups  # { # Protein Groups : indices }
	# DEFAULTDICT doesn't return a KeyError when key not found, but rather None. !!! so you can safely .extend()
	minProteinPeptidesDict = None  # proteins get contribution only from peptides which correspond uniquely to them
	fullProteinPeptidesDict = None  # proteins get maximal contribution from all corresponding peptides even if corresponding to multiple proteins
	noMasterProteinAccession = []
	for numGroups, peptideIndices in numProteinGroupsDict.items():
		if numGroups == 0:
			logging.warning("Peptides without Master Protein Accession detected. Omitting them in the analysis.")
			noMasterProteinAccession.extend(peptideIndices)
		elif numGroups == 1:  # these have only 1 master protein accession
			# { protein : indices }
			minProteinPeptidesDict = df.loc[peptideIndices].groupby("Master Protein Accessions").groups
			fullProteinPeptidesDict = defaultdict(list, minProteinPeptidesDict.copy())
		elif fullExpression_bool:  # multiple proteins accessions per peptide: save those to fullProteinPeptidesDict only.
			# { multiple proteins : indices }
			multipleProteinPeptidesDict = df.loc[peptideIndices].groupby("Master Protein Accessions").groups
			# cast dict values from int64index to list # todo find a non-ugly fix
			fullProteinPeptidesDict = defaultdict(list, dict((k,list(v)) for k,v in fullProteinPeptidesDict.items()))  # ugly
			for multipleProteinsString, nonUniqueIndices in multipleProteinPeptidesDict.items():
				multipleProteins = multipleProteinsString.split('; ')
				for protein in multipleProteins: # extend the possibly (probably) already existing entry in the dict.
					fullProteinPeptidesDict[protein].extend(nonUniqueIndices)
	# 3rd return argument must be a dataframe!
	if not fullExpression_bool:  # should return None
		# todo also allow only max expression
		fullProteinPeptidesDict = None
	return minProteinPeptidesDict, fullProteinPeptidesDict, df.loc[noMasterProteinAccession, ['First Scan', 'Sequence']]

# @profile
def getProteinDF(df, proteinPeptidesDict, schema, referenceCondition, otherConditions):
	"""
	Transform the data index from the peptide to the protein level by using the associations in proteinPeptidesDict, and
	select only relevant information. The resulting dataframe is indexed on the Master Protein Accessions.
	Quantification values from all associated unique modified peptides are averaged per channel/sample and kept in a
	list per condition for each protein. The total (before averaging) number of peptide observations per protein is also
	kept (this is thus NOT the number of observations used in the t-test, but contains also repeated measurements).
	:param df:					pd.DataFrame				outer join of all experiments.
	:param proteinPeptidesDict:	{protein: [peptide ids]}	proteins and their associated peptide ids.
	:param schema:				dict                		schema of the experiments' hierarchy.
	:param referenceCondition:	str							reference condition for the fold change calculation.
	:param otherConditions:		list 						all non-reference conditions in the experiment
	:return proteinDF:			pd.DataFrame				transformed and selected data on the protein level.
															Structure:
															['protein', 'peptides', 'description',
															referenceCondition, condition 1, condition 2, ...]
	"""
	def uniqueMods(dfModSeries):
		""" Selects only the unique, non-TMT modifications. """
		# make set of all occurring modifications
		allMods = set(unnest(dfModSeries.values))
		# The TMT6plex ones are now removed in processing
		# # this N-terminal one is ALWAYS present -> redundant information
		# allMods.remove('TMT6plex') if 'TMT6plex' in allMods else None
		return list(allMods)
	
	def emptyProteinEntry(this_allConditions):
		""" Create an empty Protein Entry to fill. This has to be a separate function to start out completely blank. """
		this_emptyProteinEntry = pd.Series(data=[None, ] * len(proteinDFColumns), index=proteinDFColumns)
		for i in this_emptyProteinEntry.index:  # set initial values
			if '#peptides' in i:
				this_emptyProteinEntry[i] = 0
			elif i in this_allConditions:
				this_emptyProteinEntry[i] = []
		return this_emptyProteinEntry
	
	if 'Protein Descriptions' not in df.columns.values:  # in case there was no Descriptions column in the input
		df['Protein Descriptions'] = pd.Series()
	if 'Modifications' not in df.columns.values:  # in case there was no Modifications column in the input
		df['Modifications'] = pd.Series()
	proteinDFColumns = ['protein', 'peptides', 'description', 'modifications']
	allConditions = [referenceCondition]+otherConditions
	# append reference condition first, and then the other conditions
	proteinDFColumns.extend([referenceCondition, '#peptides (' + str(referenceCondition) + ')'])
	# weave #peptide columns in between condition columns
	proteinDFColumns.extend(unnest(list(zip(otherConditions, ['#peptides (' + str(c) + ')' for c in otherConditions]))))
	proteinDF = pd.DataFrame([list(proteinPeptidesDict.keys())].extend([[None], ]*len(proteinDFColumns)),
							 columns=proteinDFColumns).set_index('protein')
	# remove protein from column list because it is now the index (use this list later on)
	proteinDFColumns.remove('protein')
	
	for protein, peptideIndices in proteinPeptidesDict.items():
		# construct the new protein entry, with empty quan lists for now, and add it to the proteinDF
		# the .loc lines consume virtually all of the computing time. But I increased performance by a factor ~20 by
		# splitting the selection and manipulation operations
		proteinData = df.loc[peptideIndices, ['Sequence', 'Protein Descriptions', 'Modifications']]
		# start out with empty protein entry
		proteinEntry = emptyProteinEntry(allConditions)
		proteinEntry[['peptides', 'description', 'modifications']] = [list(set(proteinData['Sequence'])),
																	  proteinData['Protein Descriptions'][0],
																	  uniqueMods(proteinData['Modifications']), ]
		# interpret as multi index so that you can call .levels and .get_level_values()
		peptideIndices = pd.MultiIndex.from_tuples(peptideIndices)
		# per experiment, get the list of indices and append the corresponding df values to the right quanPerCondition.
		for eName in peptideIndices.levels[0]:  # peptideIndices.levels[0] is the experimentName part of the index.
			# get the indices of the current experiment
			peptideIndicesPerExperiment = peptideIndices.values[peptideIndices.get_level_values(0) == eName]
			# retrieve all quan values for these peptides
			allPeptidesQuanValues = df.loc[peptideIndicesPerExperiment, schema[eName]['allExperimentChannelAliases']]
			# for each condition in the experiment, append all its channels to the quanPerCondition Series.
			for condition in schema[eName]['allExperimentConditions']:
				for channel in schema[eName][condition]['channelAliases']:
					# variable assignment because a pd.Series is returned after the append operation
					# NP.NANMEAN TO AVOID UNDERSTIMATION OF VARIANCE BECAUSE OF REPEATED MEASUREMENTS (stay conservative)
					proteinEntry[condition].append(np.nanmean(pd.Series(allPeptidesQuanValues[channel]).dropna()))
					# keep track of the total number of OBSERVED (!=used in DEA) peptides
					proteinEntry['#peptides (' + str(condition) + ')'] += len(pd.Series(allPeptidesQuanValues[channel]).dropna())
		# fill new dataframe
		proteinDF.loc[protein, :] = proteinEntry
	
	return proteinDF


def testDifferentialExpression(this_proteinDF, alpha, referenceCondition, otherConditions):
	"""
	Perform a t-test for independent samples for each protein on its (2) associated lists of peptide quantifications,
	do a Benjamini-Hochberg correction and store the results in new "p-values" and "adjusted p-values" columns in the
	dataframe. If a test returns NaN or masked values (e.g. due to sample size==1) the corresponding protein is removed.
	:param this_proteinDF:		pd.DataFrame	data from all experiments on the protein level
	:param alpha:				float			confidence level for the t-test
	:param referenceCondition:	str				name of the condition to be used as the reference
	:param otherConditions:		list 			all non-reference conditions in the experiment
	:return this_proteinDF:		pd.DataFrame	data on protein level, with statistical test information
	:return removedData:		pd.DataFrame	protein entries removed due to invalid t-test results
	"""
	# { protein : indices of (uniquely/all) associated peptides }
	# perform t-test on the intensities lists of both conditions of each protein, assuming data is independent.
	for condition in otherConditions:
		pValueColumn = 'p-value (' + condition + ')'
		this_proteinDF[pValueColumn] = this_proteinDF.apply(lambda x: ttest(x[referenceCondition], x[condition], nan_policy='omit')[1], axis=1)
		# remove masked values
		this_proteinDF.loc[:, pValueColumn] = this_proteinDF.loc[:, pValueColumn].apply(lambda x: np.nan if x is np.ma.masked or x == 0.0 else x)
		toDeleteProteins = this_proteinDF[np.isnan(this_proteinDF[pValueColumn])].index
		removedData = this_proteinDF.loc[toDeleteProteins, :].copy()
		this_proteinDF.drop(toDeleteProteins, inplace=True)
		# Benjamini-Hochberg correction
		# is_sorted==false &&returnsorted==false makes sure that the output is in the same order as the input.
		__, this_proteinDF['adjusted '+pValueColumn], __, __ = multipletests(pvals=np.asarray(this_proteinDF.loc[:, pValueColumn]),
																	   alpha=alpha, method='fdr_bh', is_sorted=False, returnsorted=False)
	return this_proteinDF, removedData


def applyFoldChange(proteinDF, pept2protCombinationMethod, referenceCondition, otherConditions):
	"""
	Calculate the log2 fold change of the quantification values per channel for each protein according to
	pept2protCombinationMethod and add it to the new "log2 fold change (conditionName)" columns.
	:param proteinDF:					pd.DataFrame	data on the protein level with t-test results.
	:param pept2protCombinationMethod:  str				method for reducing peptide information into one figure per protein
	:param referenceCondition:	str				name of the condition to be used as the reference
	:param otherConditions:		list 			all non-reference conditions in the experiment
	:return proteinDF:					pd.DataFrame	data on the protein level, including fold changes
	"""
	for condition in otherConditions:
		if pept2protCombinationMethod == 'mean':
			proteinDF['log2 fold change ('+condition+')'] = proteinDF.apply(lambda x: np.log2(np.nanmean(x[condition])/np.nanmean(x[referenceCondition])), axis=1)
		elif pept2protCombinationMethod == 'median':
			proteinDF['log2 fold change ('+condition+')'] = proteinDF.apply(lambda x: np.log2(np.nanmedian(x[condition]) / np.nanmedian(x[referenceCondition])), axis=1)
		else:
			raise Exception("Illegal pept2protCombinationMethod '"+str(pept2protCombinationMethod)+"'.")
	return proteinDF


def applySignificance(df, otherConditions, alpha, FCThreshold):
	"""
	Adds a column with the significance level to the dataframe of proteins; specifies whether the DEA or fold change
	results or both were significant.
	:param df:          	pd.DataFrame    proteins with their DEA and FC results.
	:param otherConditions:	[ str ]			names of all non-reference conditions
	:param alpha:       	float           significance level
	:param FCThreshold: 	float           fold change threshold
	:return:            	pd.DataFrame    protein data with significance levels 'yes', 'no', 'p' or 'fc'.
	"""
	for condition in otherConditions:
		def significant(x):
			pvalueSignificant = x['adjusted p-value ('+condition+')'] < alpha
			FCSignificant = abs(x['log2 fold change ('+condition+')']) > FCThreshold
			if pvalueSignificant & FCSignificant:
				return 'yes'
			elif pvalueSignificant:
				return 'p'
			elif FCSignificant:
				return 'fc'
			else:
				return 'no'
		
		df['significant ('+condition+')'] = df.apply(significant, axis=1)
	return df


def getCommonPeptidesQuanValuesDF(dfs, schema):
	"""
	Takes a list of dataframes and selects only the sequence, modifications and intensities, then inner joins them on
	sequence and modifications.	The result is the intensity matrix as a dataframe with header with ALL experiment
	channels per peptide, for only the COMMON modified, non-redundant peptides i.e. those peptides detected in ALL
	experiments. Also returns list/DataFrame of peptides that are not common for all experiments.
	:param dfs:     [ pd.DataFrame ]	data of the experiments
	:param schema:  dict                schema of the experiments' hierarchy
	:return:        pd.DataFrame		[e1_channel1, e1_channel2, ..., eM_channel1, ..., eM_channelN] for all COMMON peptides.
	"""
	allChannelAliases = unnest([schema[eName]['allExperimentChannelAliases'] for eName in schema['allExperiments']])
	peptidesDf = pd.DataFrame()
	# join all dataframes together on the Sequence: you get ALL channels from ALL experiments as columns per peptide.
	# [peptide, e1_channel1, e1_channel2, ..., eM_channel1, ..., eM_channelN]
	allModifiedPeptides = set()
	for eName in dfs.keys():
		eChannelAliases = schema[eName]['allExperimentChannelAliases']
		# convert modifications list to strings because they need to be hashable
		dfs[eName]['Modifications'] = dfs[eName]['Modifications'].apply(';'.join)
		if peptidesDf.empty:
			peptidesDf = dfs[eName].loc[:, ['Sequence', 'Modifications'] + eChannelAliases]
		else:
			# Merge makes sure that only peptides whose Sequence appears in EACH experiment get selected
			# If a channel in eChannelAliases does not exist in this eName, then it will return a column of NaNs (good).
			# The same applies to the Modifications column: if there are none, there is no problem.
			peptidesDf = pd.merge(peptidesDf, dfs[eName].loc[:, ['Sequence', 'Modifications'] + eChannelAliases],
							 on=['Sequence', 'Modifications'])
		peptidesAndModifications = zip(dfs[eName]['Sequence'], dfs[eName]['Modifications'])
		allModifiedPeptides.update(set(peptidesAndModifications))
	allCommonModifiedPeptides = set(zip(peptidesDf['Sequence'], peptidesDf['Modifications']))
	uncommonModifiedPeptides = pd.DataFrame(list(allModifiedPeptides.difference(allCommonModifiedPeptides)))
	if len(peptidesDf) < 2:
		raise Exception("Only "+str(len(peptidesDf))+" peptides found that were common across all experiments. Cannot perform PCA nor HC.")
	return peptidesDf.loc[:, allChannelAliases], uncommonModifiedPeptides


# def getNumUnCommonModifiedPeptidesPerCondition(dfs, uncommonModifiedPeptides, schema):
# 	allModifiedPeptidesPerCondition = dict()
# 	for k in schema['allConditions']:
# 		for eName in schema['allExperiments']:
# 			if k in schema[eName]:
# 				schema[eName][k]['channelAliases']
# 				# N.B.: common peptides are defined to appear in every experiment, but in practise they appear in every
# 				# condition as well (if not, they would already have been filtered out in processing).
# 				allModifiedPeptidesPerCondition[k] = zip(dfs[eName]['Sequence'], dfs[eName]['Modifications'])
#
# 	return numUnCommonModifiedPeptidesPerCondition


def getPCA(intensities, nComponents=2):
	"""
	Returns the nComponents Principal Component scores for the transposed intensity matrix. This means the reporter
	channels are "observations" with each protein intensity as a variable/attribute. The fast randomized method by Halko
	et al. (2009) is used for calculating the SVD. NaN values are converted to zero!!!
	:param intensities: np.ndarray  MxN ndarray with intensities
	:param nComponents: int         number of PC to keep
	:return:            np.ndarray  principal component scores of the input intensities
	"""
	pca = PCA(n_components=nComponents, svd_solver='randomized')
	# assign zero so that the PCA doesn't fail! This is OK, because NaN means that either the intensity was so low that
	# it could not be detected, or it just wasn't present at all. Both cases: close to zero.
	# ALSO this doesn't affect the row sums.
	intensities[np.isnan(intensities)] = 0
	pca.fit(intensities.T)
	return pca.transform(intensities.T)


def getHC(intensities):
	"""
	Perform hierarchical clustering on the transposed intensity matrix, with nComponents principal components.
	This means the reporter channels are "observations" with each protein intensity as a variable/attribute.
	Returns the (NxN) linkage matrix describing the distances between each observation (reporter channel) according to
	the UPGMA algorithm. NaN values are converted to zero!!!
	See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html for more information.
	:param intensities: np.ndarray  MxN ndarray with intensities
	:param nClusters:   int         number of clusters we want to find (= number of conditions in the experiment(s))
	:return:            np.ndarray  Nx4 linkage matrix
	"""
	# assign zero so that the PCA doesn't fail! This is OK, because NaN means that either the intensity was so low that
	# it could not be detected, or it just wasn't present at all. Both cases: close to zero.
	# ALSO this doesn't affect the row sums.
	intensities[np.isnan(intensities)] = 0
	condensedDistanceMatrix = pdist(intensities.T) # remove nans and transpose
	return linkage(condensedDistanceMatrix, method='average') # 'average'=UPGMA


def buildHandyColumnOrder(inColumns, referenceCondition, schema):
	"""
	Reshuffles the entries of a given list of proteinDF columns into a handy order.
	:param inColumns:			[ str ]	proteinDF columns, unordered
	:param referenceCondition:	str		reference condition
	:param schema:				dict    schema of the experiments' hierarchy
	:return outColumns:			[ str ]	proteinDF columns, ordered according to: ['protein', 'description',
										{{ 'adjusted p-value', 'log2 fold change'}}, {{'#peptides'}}, {{'significant'}},
										{{'p-value'}}, {{'condition'}} , 'peptides']
	"""
	otherConditions = getOtherConditions(schema, referenceCondition)
	# ['protein', 'description', 'adjusted p-value', 'log2 fold change', '#peptides (c1, c2)', 'significant', 'p-value',
	# 'condition 1', 'condition 2', 'peptides']
	outColumns = ['description']
	for condition in otherConditions:
		outColumns.extend(['adjusted p-value ('+condition+')', 'log2 fold change ('+condition+')'])
	outColumns.append('#peptides ('+referenceCondition+')')
	for condition in otherConditions:
		outColumns.append('#peptides ('+condition+')')
	for condition in otherConditions:
		outColumns.append('significant ('+condition+')')
	for condition in otherConditions:
		outColumns.append('p-value ('+condition+')')
	outColumns.append(referenceCondition)
	for condition in otherConditions:
		outColumns.append(condition)
	outColumns.append('peptides')
	outColumns.append('modifications')
	assert len(inColumns) == len(outColumns)
	return outColumns

