#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions involved in analyzing the data that was processed by dataprep.py and constand.py.
Performs a differential expression analysis on the normalized intensities as provided by CONSTANd.
Includes data visualization.
"""

import numpy as np
import pandas as pd


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


def differentialExpression(normalizedIntensities, threshold=1):
	# TODO: only include differentials with a fold of >threshold or <1/threshold
	# TODO: careful with peptides with more than 1 master protein
	return None


# pept2prot = mapping of peptides to proteins (1 to many) --> extract directly from df
# prot2peptMin = mapping of proteins to peptides (1 to many) --> calculate from pept2prot
# prot2peptMax = mapping of proteins to peptides (1 to many more) --> calculate from pept2prot
def mapPeptToProt(df, peptides):
	# calculate prot2peptMin/Max and alongside it protIntensitiesMin/Max:
	# Min is the case where proteins are only associated with 1 to 1 matching peptides in pept2prot,
	# Max is the case where proteins are associated with all matching peptides in pept2prot
	#return protIntensitiesMin, protIntensitiesMax, prot2peptMin, prot2peptMax
	pass

def dataVisualization(DEresults):
	# TODO (if paying customer): parameter: intensity matrix on peptide or protein level?
	return None



