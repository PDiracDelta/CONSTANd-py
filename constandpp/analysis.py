#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions involved in analyzing the data that was processed by dataproc.py and constand.py.
Performs a differential expression analysis on the normalized intensities as provided by CONSTANd.
"""

import numpy as np
import pandas as pd
from warnings import warn
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
	:param removedData_RT:  pd.dataFrame    removedData for the RT collapse.
	:return:                pd.DataFrame    statistics 'Degeneracy', 'mean', 'std', 'max-min' about the RT values
	"""
	duplicateGroups = removedData_RT.groupby('Representative First Scan').groups
	RTIsolationInfo = []
	for rfs, duplicates in duplicateGroups.items():
		RTValues = removedData_RT.loc[duplicates, 'RT [min]']
		RTIsolationInfo.append([rfs, len(duplicates), np.nanmean(RTValues), np.std(RTValues), np.ptp(RTValues)])
	return pd.DataFrame(RTIsolationInfo, columns=['Representative First Scan', 'Degeneracy', 'mean', 'std', 'max-min'])


def getNoIsotopicCorrection(df, noCorrectionIndices):
	"""
	Given a dataframe and indices of detections that received no corrections, returns some basic info about them.
	:param df:                  pd.dataFrame
	:param noCorrectionIndices: list            indices of detections that received no isotopic correction
	:return:                    pd.dataFrame    ['First Scan', 'Identifying Node', 'Annotated Sequence', 'Master Protein Accessions']
	"""
	return df.loc[noCorrectionIndices, ['First Scan', 'Identifying Node', 'Annotated Sequence', 'Master Protein Accessions']]


def combineExperimentDFs(dfs): #, schema):
	# """
	# All dataframes first have their intensityColumns replaced by their aliases so that they remain distinguishable after
	# the merge. Then the dataframes are merged and attain a multi-index (experimentName, oldIndex).
	# :param dfs:     dict of pd.DataFrames   data to be merged
	# :param schema:  dict of dicts           schema describing the hierarchy of experiments and conditions and channels
	# :return:        pd.DataFrame            data of ALL experiments, with distinguishable channels and multi-indexed on
	# 										experiment name and old index
	# """
	# for eName in schema:
	# 	channelAliases = [channel for condition in schema[eName]['channelAliasesPerCondition'] for channel in condition]
	# 	channels = [channel for condition in schema[eName]['intensityColumnsPerCondition'] for channel in condition]
	# 	for channel,alias in zip(channels,channelAliases):
	# 		# replace each old intensityColumn name by its new alias
	# 		dfs[eName].columns[dfs[eName].columns.index(channel)] = alias
	return pd.concat(dfs.values(), keys=dfs.keys())


def getProteinPeptidesDicts(df):
	"""
	Returns two dicts with the peptide indices (w.r.t. dataframe df) associated with each protein in the df as a
	dictionary. One dict (min) contains only the peptide indices uniquely associated per protein, the other contains
	all peptides indices associated per protein.
	:return minProteinPeptidesDict:	dict	{ protein : uniquely associated peptide indices }
	:return maxProteinPeptidesDict:	dict	{ protein : all associated peptide indices }
	:return 						list	peptide sequences without master protein accession
	"""
	# todo make this function independent of # protein groups
	numProteinGroupsDict = df.groupby("# Protein Groups").groups  # { # Protein Groups : indices }
	# DEFAULTDICT doesn't return a KeyError when key not found, but rather None. !!! so you can safely .extend()
	minProteinPeptidesDict = None  # proteins get contribution only from peptides which correspond uniquely to them
	maxProteinPeptidesDict = None  # proteins get maximal contribution from all corresponding peptides even if corresponding to multiple proteins
	noMasterProteinAccession = []
	for numGroups, peptideIndices in numProteinGroupsDict.items():
		if numGroups == 0:
			warn("Peptides without Master Protein Accession detected. Omitting them in the analysis.")
			noMasterProteinAccession.extend(peptideIndices)
		elif numGroups == 1:  # these have only 1 master protein accession
			# { protein : indices }
			minProteinPeptidesDict = df.loc[peptideIndices].groupby("Master Protein Accessions").groups
			maxProteinPeptidesDict = defaultdict(list, minProteinPeptidesDict.copy())
		else:  # multiple proteins accessions per peptide: save those to maxProteinPeptidesDict only.
			# { multiple proteins : indices }
			multipleProteinPeptidesDict = df.loc[peptideIndices].groupby("Master Protein Accessions").groups
			# cast dict values from int64index to list # todo find a non-ugly fix
			maxProteinPeptidesDict = defaultdict(list, dict((k,list(v)) for k,v in maxProteinPeptidesDict.items()))  # ugly
			for multipleProteinsString, nonUniqueIndices in multipleProteinPeptidesDict.items():
				multipleProteins = multipleProteinsString.split('; ')
				for protein in multipleProteins: # extend the possibly (probably) already existing entry in the dict.
					maxProteinPeptidesDict[protein].extend(nonUniqueIndices)
	# 3rd return argument must be a dataframe!
	return minProteinPeptidesDict, maxProteinPeptidesDict, df.loc[noMasterProteinAccession, ['First Scan', 'Annotated Sequence']]


def getProteinDF(df, proteinPeptidesDict, schema):
	# todo docu
	#proteinDF = pd.DataFrame([list(proteinPeptidesDict.keys())].extend([[None, ]*len(proteinPeptidesDict.keys()), ]*3),
	proteinDF = pd.DataFrame(list(proteinPeptidesDict.keys()),
	                         columns=['protein', 'peptides', 'description', 'condition 1', 'condition 2']).set_index('protein')
	# define intensityColumnsPerConditionDict so that it contains the channels of ALL experiments
	channelAliasesPerConditionDict = dict((eName, experiment['channelAliasesPerCondition']) for eName, experiment in schema.items())
	for protein, peptideIndices in proteinPeptidesDict.items():
		# combine all channels into one channel per condition. peptideIndices[0] = experimentName
		condition1Intensities = pd.concat([df.loc[peptideIndices, channel] for channel in
		                                   channelAliasesPerConditionDict[peptideIndices[0]][0]], axis=0, ignore_index=True).tolist()
		condition2Intensities = pd.concat([df.loc[peptideIndices, channel] for channel in
		                                   channelAliasesPerConditionDict[peptideIndices[0]][1]], axis=0, ignore_index=True).tolist()
		# fill new dataframe on protein level, per condition
		proteinDF.loc[protein, :] = [df.loc[peptideIndices, 'Annotated Sequence'].tolist(),
		                             df.loc[peptideIndices, 'Protein Descriptions'].tolist(),
		                             condition1Intensities, condition2Intensities]
	return proteinDF


def applyDifferentialExpression(this_proteinDF, alpha):
	# todo docu
	# { protein : indices of (uniquely/all) associated peptides }
	# perform t-test on the intensities lists of both conditions of each protein, assuming data is independent.
	this_proteinDF['p-value'] = this_proteinDF.apply(lambda x: ttest(x['condition 1'], x['condition 2'], nan_policy='omit')[1], axis=1)
	# remove masked values
	this_proteinDF.loc[:, 'p-value'] = this_proteinDF.loc[:, 'p-value'].apply(lambda x: np.nan if x is np.ma.masked else x)
	# Benjamini-Hochberg correction
	# is_sorted==false &&returnsorted==false makes sure that the output is in the same order as the input.
	__, this_proteinDF['adjusted p-value'], __, __ = multipletests(pvals=np.asarray(this_proteinDF.loc[:, 'p-value']),
																   alpha=alpha, method='fdr_bh', is_sorted=False, returnsorted=False)
	return this_proteinDF


def applyFoldChange(proteinDF, pept2protCombinationMethod):
	""" Calculate the fold change for each protein (pept2protCombinationMethod) and apply it to the given protein dataframe """
	# todo proper docu
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


def getPCA(intensities, nComponents):
	"""
	Returns the nComponents Principal scores for the transposed intensity matrix. This means the reporter channels
	are "observations" with each protein intensity as a variable/attribute. The fast randomized method by Halko et al.
	(2009) is used for calculating the SVD. Missing values are thrown away.
	:param intensities: np.ndarray  MxN ndarray with intensities
	:param nComponents: int         number of PC to keep
	:return:            np.ndarray  principal component scores of the input intensities
	"""
	pca = PCA(n_components=nComponents, svd_solver='randomized')
	transposedIntensitiesNoNaN = intensities[~np.isnan(intensities).any(axis=1)].T
	pca.fit(transposedIntensitiesNoNaN)
	return pca.transform(transposedIntensitiesNoNaN)


def getHC(intensities):
	"""
	Perform hierarchical clustering on the transposed intensity matrix, with nComponents principal components.
	This means the reporter channels are "observations" with each protein intensity as a variable/attribute.
	Returns the (NxN) linkage matrix describing the distances between each observation (reporter channel) according to
	the UPGMA algorithm. Peptides with NaN values are removed from the data!!!
	See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html for more information.
	:param intensities: np.ndarray  MxN ndarray with intensities
	:param nClusters:   int         number of clusters we want to find (= number of conditions in the experiment(s))
	:return:            np.ndarray  Nx4 linkage matrix
	"""
	condensedDistanceMatrix = pdist(intensities[~np.isnan(intensities).any(axis=1)].T) # remove nans and transpose
	return linkage(condensedDistanceMatrix, method='average') # 'average'=UPGMA
