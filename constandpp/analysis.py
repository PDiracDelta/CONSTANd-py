#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions involved in analyzing the data that was processed by dataproc.py and constand.py.
Performs a differential expression analysis on the normalized intensities as provided by CONSTANd.
Includes data visualization.
"""

import numpy as np
import pandas as pd
from collections import defaultdict
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import ttest_ind as ttest


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


def getProteinPeptidesDicts(df):
	"""
	Returns two dicts with the peptide indices (w.r.t. dataframe df) associated with each protein in the df as a
	dictionary. One dict (min) contains only the peptide indices uniquely associated per protein, the other contains
	all peptides indices associated per protein.
	:return minProteinPeptidesDict:	dict	{ protein : uniquely associated peptide indices }
	:return maxProteinPeptidesDict:	dict	{ protein : all associated peptide indices }
	""" # todo make this function independent of # protein groups
	numProteinGroupsDict = df.groupby("# Protein Groups").groups  # { # Protein Groups : indices }
	# DEFAULTDICT doesn't return a KeyError when key not found, but rather None. !!! so you can safely .extend()
	minProteinPeptidesDict = defaultdict(
		list)  # proteins get contribution only from peptides which correspond uniquely to them
	maxProteinPeptidesDict = defaultdict(
		list)  # proteins get maximal contribution from all corresponding peptides even if corresponding to multiple proteins
	for numGroups in numProteinGroupsDict.keys():
		if numGroups == 1:  # these have only 1 master protein accession
			# { protein : indices }
			minProteinPeptidesDict.extend(
				df[numProteinGroupsDict[numGroups]].groupby("Master Protein Accessions").groups)
			maxProteinPeptidesDict = minProteinPeptidesDict.copy()
		else:  # multiple proteins accessions per peptide: save those to maxProteinPeptidesDict only.
			# { multiple proteins : indices }
			multipleProteinPeptidesDict = df.loc[numProteinGroupsDict[numGroups]].groupby(
				"Master Protein Accessions").groups
			for multipleProteinsString in multipleProteinPeptidesDict.keys():
				multipleProteins = multipleProteinsString.split('; ')
				for protein in multipleProteins:
					maxProteinPeptidesDict[protein].extend(multipleProteinPeptidesDict[multipleProteinsString])
	return minProteinPeptidesDict, maxProteinPeptidesDict


def getProteinDF(df, proteinPeptidesDict, intensityColumnsPerCondition):
	proteinDF = pd.DataFrame([list(proteinPeptidesDict.keys())].extend([[None, ]*len(proteinPeptidesDict.keys()), ]*3),
	                         columns=['protein', 'peptides', 'condition 1', 'condition 2']).set_index('protein')
	for protein, peptideIndices in proteinPeptidesDict.items():
		proteinDF.loc[protein] = [df.loc[peptideIndices, 'Annotated Sequence'],
		                          df.loc[peptideIndices, intensityColumnsPerCondition[0]],
		                          df.loc[peptideIndices, intensityColumnsPerCondition[1]]]
	return proteinDF


def applyDifferentialExpression(this_proteinDF, alpha):
	# TODO: careful with peptides with more than 1 master protein
	# { protein : indices of (uniquely/all) associated peptides }
	# def DE(row): # TEST obsolete code?
	# 	"""
	# 	Test the differential expression between the intensity lists of the two conditions for a certain protein.
	# 	:param row: pd.Series   row of a certain protein in the dataframe two intensity input lists for the t-test
	# 	:return:    float64     p-value of the t-test
	# 	"""
	# 	__, p = ttest(row['condition 1'], row['condition 2'])
	# 	return p
	this_proteinDF['p-value'] = [np.nan, ] * len(this_proteinDF.index)
	# perform t-test on the intensities lists of both conditions of each protein, assuming data is independent.
	this_proteinDF.loc[:, 'p-value'] = this_proteinDF.apply(
		lambda x: ttest(x['condition 1'], x['condition 2']), axis=1)
	# Benjamini-Hochberg correction
	# is_sorted==false &&returnsorted==false makes sure that the output is in the same order as the input.
	this_proteinDF['adjusted p-value'] = multipletests(
		pvals=this_proteinDF.loc[:, 'p-value'], alpha=alpha, method='fdr_bh', is_sorted=False, returnsorted=False)
	return this_proteinDF


def applyFoldChange(proteinDF, pept2protCombinationMethod):
	""" Calculate the fold change for each protein (pept2protCombinationMethod) and apply it to the given protein dataframe """
	return None


def dataVisualization(DEresults, FCThreshold):
	# TODO (if paying customer): parameter: intensity matrix on peptide or protein level?
	# TODO: only include differentials with a fold of >threshold or <1/threshold
	return None



