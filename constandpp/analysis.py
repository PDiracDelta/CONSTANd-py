#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions involved in analyzing the data that was processed by processing.py and constand.py.
Performs a differential expression analysis on the normalized intensities as provided by CONSTANd.
"""

import numpy as np
import pandas as pd
import logging
from constandpp.tools import unnest
# from constandpp.tools import getIntensities
from collections import defaultdict
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import ttest_ind as ttest
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage


def getRTIsolationInfo(removedData_RT):
	"""
	Returns dataFrame with the mean, standard deviation, and max-min value of the RT values for each duplicate_group
	representative that is found in the removedData_RT dataframe containing the removedData for the RT collapse.
	:param removedData_RT:  pd.DataFrame    removedData for the RT collapse.
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
	Given a dataframe and indices of detections that received no corrections, returns some basic info about them.
	:param df:                  pd.DataFrame	processed and possibly (if not: exception caught) isotope-corrected data frame
	:param noCorrectionIndices: list            indices of detections that received no isotopic correction
	:return:                    pd.DataFrame    ['First Scan', 'Identifying Node Type', 'Annotated Sequence', 'Master Protein Accessions']
	"""
	return df.loc[noCorrectionIndices, ['First Scan', 'Identifying Node Type', 'Annotated Sequence', 'Master Protein Accessions']]


def combineExperimentDFs(dfs):
	"""
	Merge dataframes of all experiments into one multi-indexed (eName, oldIndex) dataframe, by performing an outer join.
	The intensity columns are non-identical across different dataframes and resulting empty fields are valued NaN.
	:param dfs:     [pd.DataFrame]	dataframes from multiple experiments
	:return: 		pd.DataFrame	dataframe containing an outer join of the input list of dataframes
	"""
	return pd.concat(dfs.values(), keys=dfs.keys())


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
		df["# Protein Groups"] = df.apply(lambda x: len(x['Master Protein Accessions'].split(';')), axis=1)
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
	if fullExpression_bool:
		return minProteinPeptidesDict, fullProteinPeptidesDict, df.loc[noMasterProteinAccession, ['First Scan', 'Annotated Sequence']]
	else:  # todo also allow only max expression
		return minProteinPeptidesDict, None, df.loc[noMasterProteinAccession, ['First Scan', 'Annotated Sequence']]


def getProteinDF(df, proteinPeptidesDict, schema):
	"""
	Transform the data index from the peptide to the protein level by using the associations in proteinPeptidesDict, and
	select only relevant information. The resulting dataframe is indexed on the Master Protein Accessions.
	Quantification values from all associated unique modified peptides are kept in a list per condition (max: 2) for
	each protein.
	:param df:					pd.DataFrame				outer join of all experiments.
	:param proteinPeptidesDict:	{protein: [peptide ids]}	proteins and their associated peptide ids.
	:param schema:				dict                		schema of the experiments' hierarchy.
	:return proteinDF:			pd.DataFrame				transformed and selected data on the protein level.
															Structure:
															['protein', 'peptides', 'description', 'condition 1', 'condition 2']
	"""
	proteinDFColumns = ['protein', 'peptides', 'description', 'condition 1', 'condition 2']
	proteinDF = pd.DataFrame([list(proteinPeptidesDict.keys())].extend([[None], ]*len(proteinDFColumns)),
							 columns=proteinDFColumns).set_index('protein')
	# define channelAliasesPerConditionDict so that it contains the channels of ALL experiments
	channelAliasesPerConditionDict = dict((eName, experiment['channelAliasesPerCondition']) for eName, experiment in schema.items())
	for protein, peptideIndices in proteinPeptidesDict.items():
		# interpret as multi index so that you can call .levels and .get_level_values()
		peptideIndices = pd.MultiIndex.from_tuples(peptideIndices)
		# combine all channels into one channel per condition.
		condition1Intensities = pd.Series()
		condition2Intensities = pd.Series()
		# per experiment, get the a list of indices per channel for both conditions, and concatenate the corresponding df values.
		for eName in peptideIndices.levels[0]: # peptideIndices.levels[0] is the experimentName part of the index.
			# get the indices of the current experiment
			peptideIndicesPerExperiment = peptideIndices.values[peptideIndices.get_level_values(0) == eName]
			# get a list of dfs, to sort the intensities per channel.
			condition1intensitiesPerChannel = [df.loc[peptideIndicesPerExperiment, channel] for channel in
											   channelAliasesPerConditionDict[eName][0]]
			condition2intensitiesPerChannel = [df.loc[peptideIndicesPerExperiment, channel] for channel in
											   channelAliasesPerConditionDict[eName][1]]
			condition1Intensities = pd.concat([condition1Intensities] + condition1intensitiesPerChannel, axis=0, ignore_index=True)
			condition2Intensities = pd.concat([condition2Intensities] + condition2intensitiesPerChannel, axis=0, ignore_index=True)
		# fill new dataframe on protein level, per condition
		if 'Protein Descriptions' not in df.columns.values:
			df['Protein Descriptions'] = pd.Series()
		proteinDF.loc[protein, :] = [df.loc[peptideIndices, 'Annotated Sequence'].tolist(),
									 df.loc[peptideIndices, 'Protein Descriptions'][0],
									 condition1Intensities.tolist(), condition2Intensities.tolist()]
	return proteinDF


def testDifferentialExpression(this_proteinDF, alpha):
	"""
	Perform a t-test for independent samples for each protein on its (2) associated lists of peptide quantifications,
	do a Benjamini-Hochberg correction and store the results in new "p-values" and "adjusted p-values" columns in the
	dataframe. If a test returns NaN or masked values (e.g. due to sample size==1) the corresponding protein is removed.
	:param this_proteinDF:	pd.DataFrame	data from all experiments on the protein level
	:param alpha:			float			confidence level for the t-test
	:return this_proteinDF:	pd.DataFrame	data on protein level, with statistical test information
	:return removedData:	pd.DataFrame	protein entries removed due to invalid t-test results
	"""
	# { protein : indices of (uniquely/all) associated peptides }
	# perform t-test on the intensities lists of both conditions of each protein, assuming data is independent.
	this_proteinDF['p-value'] = this_proteinDF.apply(lambda x: ttest(x['condition 1'], x['condition 2'], nan_policy='omit')[1], axis=1)
	# remove masked values
	this_proteinDF.loc[:, 'p-value'] = this_proteinDF.loc[:, 'p-value'].apply(lambda x: np.nan if x is np.ma.masked or x == 0.0 else x)
	toDeleteProteins = this_proteinDF[np.isnan(this_proteinDF['p-value'])].index
	removedData = this_proteinDF.loc[toDeleteProteins, :].copy()
	this_proteinDF.drop(toDeleteProteins, inplace=True)
	# Benjamini-Hochberg correction
	# is_sorted==false &&returnsorted==false makes sure that the output is in the same order as the input.
	__, this_proteinDF['adjusted p-value'], __, __ = multipletests(pvals=np.asarray(this_proteinDF.loc[:, 'p-value']),
																   alpha=alpha, method='fdr_bh', is_sorted=False, returnsorted=False)
	return this_proteinDF, removedData


def applyFoldChange(proteinDF, pept2protCombinationMethod):
	"""
	Calculate the log2 fold change of the quantification values per channel for each protein according to
	pept2protCombinationMethod and add it to the new "log2 fold change c1/c2" column.
	:param proteinDF:					pd.DataFrame	data on the protein level with t-test results.
	:param pept2protCombinationMethod:  str				method for reducing peptide information into one figure per protein
	:return proteinDF:					pd.DataFrame	data on the protein level, including fold changes
	"""
	if pept2protCombinationMethod == 'mean':
		proteinDF['log2 fold change c1/c2'] = proteinDF.apply(lambda x: np.log2(np.nanmean(x['condition 1'])/np.nanmean(x['condition 2'])), axis=1)
	elif pept2protCombinationMethod == 'median':
		proteinDF['log2 fold change c1/c2'] = proteinDF.apply(lambda x: np.log2(np.nanmedian(x['condition 1']) / np.nanmedian(x['condition 2'])), axis=1)
	return proteinDF


def applySignificance(df, alpha, FCThreshold):
	"""
	Adds a column with the significance level to the dataframe of proteins; specifies whether the DEA or fold change
	results or both were significant.
	:param df:          pd.DataFrame    proteins with their DEA and FC results.
	:param alpha:       float           significance level
	:param FCThreshold: float           fold change threshold
	:return:            pd.DataFrame    protein data with significance levels 'yes', 'no', 'p' or 'fc'.
	"""
	def significant(x):
		pvalueSignificant = x['adjusted p-value'] < alpha
		FCSignificant = abs(x['log2 fold change c1/c2']) > FCThreshold
		if pvalueSignificant & FCSignificant:
			return 'yes'
		elif pvalueSignificant:
			return 'p'
		elif FCSignificant:
			return 'fc'
		else:
			return 'no'

	df['significant'] = df.apply(significant, axis=1)
	return df


def getAllExperimentsIntensitiesPerCommonPeptide(dfs, schema):
	"""
	Takes a list of dataframes and selects only the sequence and intensities, then inner joins them on sequence.
	The result is the intensity matrix as a dataframe with header with ALL experiment channels per peptide, for only the
	COMMON peptides i.e. those peptides detected in ALL experiments. Also returns list/DataFrame of peptides that are
	not common for all experiments.
	:param dfs:     [ pd.DataFrame ]	data of the experiments
	:param schema:  dict                schema of the experiments' hierarchy
	:return:        pd.DataFrame		[e1_channel1, e1_channel2, ..., eM_channel1, ..., eM_channelN] for all COMMON peptides.
	"""
	allChannelAliases = unnest([unnest(experiments['channelAliasesPerCondition']) for experiments in schema.values()])
	peptidesDf = pd.DataFrame()
	# join all dataframes together on the Annotated Sequence: you get ALL channels from ALL experiments as columns per peptide.
	# [peptide, e1_channel1, e1_channel2, ..., eM_channel1, ..., eM_channelN]
	allPeptides = []
	for eName in dfs.keys():
		eChannelAliases = unnest(schema[eName]['channelAliasesPerCondition'])
		if peptidesDf.empty:
			peptidesDf = dfs[eName].loc[:, ['Annotated Sequence'] + eChannelAliases]
		else:
			peptidesDf = pd.merge(peptidesDf, dfs[eName].loc[:, ['Annotated Sequence'] + eChannelAliases],
							 on='Annotated Sequence')
		allPeptides.extend(dfs[eName].loc[:, 'Annotated Sequence'])
	uncommonPeptides = pd.DataFrame(list(set(allPeptides).difference(set(peptidesDf.loc[:, 'Annotated Sequence']))))
	return peptidesDf.loc[:, allChannelAliases], uncommonPeptides


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
