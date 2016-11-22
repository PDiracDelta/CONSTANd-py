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
	"""
	numProteinGroupsDict = df.groupby("# Protein Groups").groups  # { # Protein Groups : indices }
	# DEFAULTDICT doesn't return a KeyError when key not found, but rather None. !!! so you can safely .extend()
	minProteinPeptidesDict = defaultdict(
		list)  # proteins get contribution only from peptides which correspond uniquely to them
	maxProteinPeptidesDict = defaultdict(
		list)  # proteins get maximal contribution from all corresponding peptides even if corresponding to multiple proteins
	for numGroups in numProteinGroupsDict.keys():
		if numGroups == '1':  # these have only 1 master protein accession
			# { protein : indices }
			minProteinPeptidesDict.extend(
				df[numProteinGroupsDict[numGroups]].groupby("Master Protein Accessions").groups)
			maxProteinPeptidesDict = minProteinPeptidesDict.copy()
		else:  # multiple proteins accessions per peptide: save those to maxProteinPeptidesDict only.
			# { multiple proteins : indices }
			multipleProteinPeptidesDict = df[numProteinGroupsDict[numGroups]].groupby(
				"Master Protein Accessions").groups
			for multipleProteinsString in multipleProteinPeptidesDict.keys():
				multipleProteins = multipleProteinsString.split('; ')
				for protein in multipleProteins:
					maxProteinPeptidesDict[protein].extend(multipleProteinPeptidesDict[multipleProteinsString])
	return minProteinPeptidesDict, maxProteinPeptidesDict

? minProteinPeptidesDict, maxProteinPeptidesDict = getProteinPeptidesDicts(df)


def differentialExpression(normalizedIntensities, threshold=1):
	# TODO: only include differentials with a fold of >threshold or <1/threshold
	# TODO: careful with peptides with more than 1 master protein
	# { protein : indices of (uniquely/all) associated peptides }
	return None


def dataVisualization(DEresults):
	# TODO (if paying customer): parameter: intensity matrix on peptide or protein level?
	return None



