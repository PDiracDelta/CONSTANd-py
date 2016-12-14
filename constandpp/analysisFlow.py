#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Workflow of the analysis part of CONSTANd++.
"""

from analysis import *
from dataIO import exportData


def analyzeProcessingResult(processingResults, params, writeToDisk):
	# todo docu
	processingResultsItems = processingResults.items()
	dfs = dict((eName, result[0]) for eName, result in processingResultsItems)
	normalizedIntensitiess = dict((eName, result[1]) for eName, result in processingResultsItems)
	removedDatas = dict((eName, result[2]) for eName, result in processingResultsItems)
	noCorrectionIndicess = {}
	for eName, result in processingResultsItems:
		if len(result) > 3: # if noCorrectionIndices exists in results
			noCorrectionIndicess[eName] = result[3]

	experimentNames = processingResults.keys()
	# contains statistics and metadata (like the parameters) about the analysis.
	metadata = {}
	# record detections without isotopic correction applied applied. Multi-indexed on experiment names and old indices!
	try:
		metadata['noIsotopicCorrection'] = pd.concat([getNoIsotopicCorrection(dfs[eName], noCorrectionIndicess[eName]) for
	                                              eName in noCorrectionIndicess.keys()], keys=experimentNames)
	except ValueError:
		pass # not a single noCorrectionIndices was found. OK.
	# record RT isolation statistics. Future: flag. Multi-indexed on experiment names and old indices!
	metadata['RTIsolationInfo'] = pd.concat([getRTIsolationInfo(removedDatas[eName]['RT']) for
	                                         eName in experimentNames], keys=experimentNames)

	# merge all experiments in multi-indexed: (eName, oldIndex) dataframe # and intensityColumns are unique and distinguishable
	allExperimentsDF = combineExperimentDFs(dfs) #, params['schema'])

	nConditions = len(list(params['schema'].values())[0]['channelAliasesPerCondition'])
	# ONLY PRODUCE VOLCANO AND DEA IF CONDITIONS == 2
	if nConditions == 2:
		# get min and max protein-peptide mappings
		minProteinPeptidesDict, maxProteinPeptidesDict, metadata['noMasterProteinAccession'] = getProteinPeptidesDicts(allExperimentsDF)

		# execute mappings to get all peptideintensities per protein, over each whole condition. Index = 'protein'
		minProteinDF = getProteinDF(allExperimentsDF, minProteinPeptidesDict, params['schema'])
		fullProteinDF = getProteinDF(allExperimentsDF, maxProteinPeptidesDict, params['schema'])

		# perform differential expression analysis with Benjamini-Hochberg correction.
		minProteinDF = applyDifferentialExpression(minProteinDF, params['alpha'])
		fullProteinDF = applyDifferentialExpression(fullProteinDF, params['alpha'])

		# calculate fold changes of the average protein expression value per CONDITION/GROUP (not per channel!)
		minProteinDF = applyFoldChange(minProteinDF, params['pept2protCombinationMethod'])
		fullProteinDF = applyFoldChange(fullProteinDF, params['pept2protCombinationMethod'])

		# indicate significance based on given thresholds alpha and FCThreshold
		minProteinDF = applySignificance(minProteinDF, params['alpha'], params['FCThreshold'])
		fullProteinDF = applySignificance(fullProteinDF, params['alpha'], params['FCThreshold'])
	else:
		minProteinDF = pd.DataFrame()
		fullProteinDF = pd.DataFrame()

	# dataframe with ALL intensities per peptide: [peptide, e1_channel1, e1_channel2, ..., eM_channel1, ..., eM_channelN]
	allExperimentsIntensitiesPerCommonPeptide, metadata['uncommonPeptides'] = getAllExperimentsIntensitiesPerCommonPeptide(dfs, params['schema'])
	# save the amount of NaN values per channel for common peptides.
	metadata['commonNanValues'] = pd.DataFrame(np.sum(np.isnan(allExperimentsIntensitiesPerCommonPeptide), axis=0))
	# perform PCA
	PCAResult = getPCA(allExperimentsIntensitiesPerCommonPeptide, params['PCA_components'])
	# perform hierarchical clustering
	HCResult = getHC(allExperimentsIntensitiesPerCommonPeptide)

	# set the protein names back as columns instead of the index, and sort the columns so the df is easier to read
	handyColumnOrder = ['protein', 'significant', 'adjusted p-value', 'log2 fold change c1/c2', 'description', 'p-value', 'peptides', 'condition 1', 'condition 2']
	minProteinDF.reset_index(level=0, inplace=True)
	fullProteinDF.reset_index(level=0, inplace=True)
	minProteinDF = minProteinDF.reindex_axis(handyColumnOrder, axis=1)
	fullProteinDF = fullProteinDF.reindex_axis(handyColumnOrder, axis=1)

	""" save results """
	if writeToDisk:
		# save the protein-level dataframes
		exportData(minProteinDF, dataType='df', path_out=params['path_out'],
		           filename=params['jobname'] + '_results_minimal', delim_out=params['delim_out'])
		exportData(fullProteinDF, dataType='df', path_out=params['path_out'],
		           filename=params['jobname'] + '_results_full', delim_out=params['delim_out'])
		# save the metadata
		exportData(metadata, dataType='df', path_out=params['path_out'],
		           filename=params['jobname'] + '_metadata',
		           delim_out=params['delim_out'], inOneFile=False)
		# generate a report PDF (without the normalized intensities: behind paywall?

	return minProteinDF, fullProteinDF, PCAResult, HCResult, metadata
