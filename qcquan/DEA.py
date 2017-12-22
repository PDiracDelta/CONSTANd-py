#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Contains all steps to be undertaken in a single DEA operation:
  - map from peptide to protein level
  - perform DE test
  - calculate fold change
  - calculate significance
"""

from qcquan.analysis import *


def DEA(this_allExperimentsDF, proteinPeptidesDict, params):
	"""
	Bring the data to the protein level in the case of [minimal]{full} expression (shared peptides are [not allowed]{allowed}).
	Execute the differential expression analysis (t-test + B-H correction, compute log2 fold change and some useful
	extra columns) and gather some metadata.
	:param this_allExperimentsDF:		pd.DataFrame	dataframe containing an outer join of all the experiment dataframes
	:param proteinPeptidesDict:			dict			{ protein : all associated peptide indices (non-injective) }
	:param params:						dict			job (global) parameters
	:return proteinDF:					pd.DataFrame	Transformed and selected data on the protein level.
														Structure:
														['peptides', 'description', {{condition}}, {{'p-value', 'adjusted p-value'}},
														{{'log2 fold change'}}, {{'significant'}}, {{'#peptides'}}]
	:return singleConditionProteins:	pd.DataFrame	protein entries removed due to invalid t-test results
	:return numProteins:				int				number of proteins taken into account in the DEA
	"""
	referenceCondition = params['referenceCondition']
	# use list() so that the new variable is not an alias
	otherConditions = getOtherConditions(params['schema'], referenceCondition)
	# execute mappings to get all peptideintensities per protein, over each whole condition. Index = 'protein'
	proteinDF = getProteinDF(this_allExperimentsDF, proteinPeptidesDict, params['schema'],
							 referenceCondition=referenceCondition, otherConditions=otherConditions)
	
	# perform differential expression analysis with Benjamini-Hochberg correction. Also remove proteins that have all
	# nan values for a certain condition and keep the removed ones in metadata
	proteinDF, singleConditionProteins = testDifferentialExpression(proteinDF, params['alpha'], referenceCondition,
																	otherConditions)
	numProteins = len(proteinDF)
	
	# calculate fold changes of the average protein expression value per CONDITION/GROUP (not per channel!)
	proteinDF = applyFoldChange(proteinDF, params['pept2protCombinationMethod'], referenceCondition,
								otherConditions)
	
	# indicate significance based on given thresholds alpha and FCThreshold
	proteinDF = applySignificance(proteinDF, otherConditions, params['alpha'], params['FCThreshold'])
	
	return proteinDF, singleConditionProteins, numProteins