#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Workflow of the analysis part of QCQuan.
"""

from constandpp.analysis import *
from constandpp.dataIO import exportData


def analyzeProcessingResult(processingResults, params, writeToDisk):
	"""
	Calls all the necessary functions to perform the analysis on a processingResults dict of objects.
	Takes the processed dataframe of each experiment in the setup, combines it into one global dataframe for all
	experiments, maps it to the protein level and performs a differential expression analysis, indicating the
	significance of each protein. The analysis can be carried out either for only injective	protein-peptide associations
	(minProteinDF) or using also the non-injective associations (fullProteinDF) or both.
	Also, on the common-peptide-level global dataframe a Principal Component Analysis is performed of which the first 2 PCs are
	saved, as well as a Hierarchical Clustering.
	Along the way, removed data as well as some metadata are saved in the corresponding variables.
	:param processingResults:	dict			for each experiment: [normalizedDf, normalizedIntensities, removedData, noCorrectionIndices]
	:param params:				dict			job (global) parameters
	:param writeToDisk:			bool			write results to disk (if not: just pass the return statement)
	:return minProteinDF:		pd.DataFrame	combined analysis data of all experiments, on the protein level, only
												for injective protein-peptide associations
	:return fullProteinDF:		pd.DataFrame	combined analysis data of all experiments, on the protein level
	:return PCAResult:			np.ndarray		PC scores for the quantification channels in the common-peptide-level global dataframe
	:return HCResult:			np.ndarray  	Nx4 linkage matrix of the quantification channels in the common-peptide-level global dataframe
	:return allExperimentsIntensitiesPerCommonPeptide:	pd.DataFrame	for all COMMON (found in all experiments) peptides:
																		[e1_channel1, e1_channel2, ..., eM_channel1, ..., eM_channelN]
	:return metadata:			pd.DataFrame	metadata: [noIsotopicCorrection, RTIsolationInfo, noMasterProteinAccession,
												minSingleConditionProteins, fullSingleConditionProteins, uncommonPeptides, commonNanValues]
	"""
	
	""" Preparation and metadata gathering """
	processingResultsItems = processingResults.items()
	dfs = dict((eName, result[0]) for eName, result in processingResultsItems)
	constandOutputs = dict((eName, result[1]) for eName, result in processingResultsItems)
	removedDatas = dict((eName, result[2]) for eName, result in processingResultsItems)
	allMasterProteinss = dict((eName, result[3]) for eName, result in processingResultsItems)
	noCorrectionIndicess = {}
	for eName, result in processingResultsItems:
		if len(result) > 4:  # if noCorrectionIndices exists in results
			noCorrectionIndicess[eName] = result[3]

	experimentNames = processingResults.keys()
	# contains statistics and metadata (like the parameters) about the analysis.
	metadata = {}
	metadata['numeric'] = pd.DataFrame()
	# Compile a list of all master proteins found at least in 1 PSM and at least in 1 experiment:
	allObservedProteins = pd.Series(list(set(unnest(allMasterProteinss.values()))))
	metadata['numeric'].loc[0, 'numObservedProteins'] = len(allObservedProteins)
	metadata['allObservedProteins'] = pd.DataFrame({'protein': allObservedProteins})
	
	# record PSMs without isotopic correction applied. Multi-indexed on experiment names and old indices!
	# This is done here instead of the processing flow because back then there was no metadata variable yet.
	try:
		metadata['noIsotopicCorrection'] = pd.concat([getNoIsotopicCorrection(dfs[eName], noCorrectionIndicess[eName]) for
												  eName in noCorrectionIndicess.keys()], keys=experimentNames)  # todo ugly
	except ValueError:
		pass  # not a single noCorrectionIndices was found. OK.
	# record RT isolation statistics. Future: flag. Multi-indexed on experiment names and old indices!
	if params['getRTIsolationInfo_bool']:
		metadata['RTIsolationInfo'] = pd.concat([getRTIsolationInfo(removedDatas[eName]['RT']) for
												 eName in experimentNames], keys=experimentNames)
	
	def DEA(this_allExperimentsDF, proteinPeptidesDict):
		"""
		Bring the data to the protein level in the case of [minimal]{full} expression (shared peptides are [not allowed]{allowed}).
		Execute the differential expression analysis (t-test + B-H correction, compute log2 fold change and some useful
		extra columns) and gather some metadata.
		:param this_allExperimentsDF:		pd.DataFrame	dataframe containing an outer join of all the experiment dataframes
		:param proteinPeptidesDict:			dict			{ protein : all associated peptide indices (non-injective) }
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
		
		# add number of peptides that represent each protein (per condition)
		proteinDF = addNumberOfRepresentingPeptides(proteinDF, referenceCondition, otherConditions)
		return proteinDF, singleConditionProteins, numProteins
	
	""" Differential Expression Analysis """
	# merge all experiments in multi-indexed: (eName, oldIndex) dataframe as an outer join
	allExperimentsDF = combineExperimentDFs(dfs)

	# nConditions = len(params['schema']['allConditions'])
	# # ONLY PRODUCE VOLCANO AND DEA IF CONDITIONS == 2
	# if nConditions == 2:
	# get min and max protein-peptide mappings
	if params['fullExpression_bool']:
		minProteinPeptidesDict, fullProteinPeptidesDict, metadata['noMasterProteinAccession'] = getProteinPeptidesDicts(allExperimentsDF, params['fullExpression_bool'])
	else:
		minProteinPeptidesDict, __, metadata['noMasterProteinAccession'] = getProteinPeptidesDicts(allExperimentsDF, params['fullExpression_bool'])
	metadata['numeric'].loc[0, 'numNoMasterProteinAccession'] = len(metadata['noMasterProteinAccession'])

	if params['minExpression_bool']:
		# Bring the data to the protein level in the case of minimal expression (no shared peptides allowed).
		# Execute the differential expression analysis and gather some metadata
		minProteinDF, metadata['minSingleConditionProteins'], metadata['numeric'].loc[0, 'minNumProteins'] = \
			DEA(allExperimentsDF, minProteinPeptidesDict)
	else:
		minProteinDF = pd.DataFrame()

	if params['fullExpression_bool']:
		# Bring the data to the protein level in the case of full expression (shared peptides allowed).
		# Execute the differential expression analysis and gather some metadata
		fullProteinDF, metadata['fullSingleConditionProteins'], metadata['numeric'].loc[0, 'fullNumProteins'] = \
			DEA(allExperimentsDF, fullProteinPeptidesDict)
	else:
		fullProteinDF = pd.DataFrame()
	
	# set the protein names back as columns instead of the index, and sort the columns so the df is easier to read
	handyColumnOrder = buildHandyColumnOrder(minProteinDF.columns, params['referenceCondition'], params['schema'])
	minProteinDF = minProteinDF.reindex_axis(handyColumnOrder, axis=1)
	fullProteinDF = fullProteinDF.reindex_axis(handyColumnOrder, axis=1)
	
	""" Quality Control """
	# dataframe with ALL intensities per peptide: [peptide, e1_channel1, e1_channel2, ..., eM_channel1, ..., eM_channelN]
	allExperimentsIntensitiesPerCommonPeptide, metadata['uncommonPeptides'] = getAllExperimentsIntensitiesPerCommonPeptide(dfs, params['schema'])
	metadata['numeric'].loc[0, 'numUnCommonPeptides'] = len(metadata['uncommonPeptides'])
	metadata['numeric'].loc[0, 'numCommonPeptides'] = len(allExperimentsIntensitiesPerCommonPeptide)
	# save the amount of NaN values per channel for common peptides.
	metadata['commonNanValues'] = pd.DataFrame(np.sum(np.isnan(allExperimentsIntensitiesPerCommonPeptide), axis=0))
	# perform PCA
	PCAResult = getPCA(allExperimentsIntensitiesPerCommonPeptide, params['PCA_components'])
	# perform hierarchical clustering
	HCResult = getHC(allExperimentsIntensitiesPerCommonPeptide)

	""" save results """
	if writeToDisk:
		# save the protein-level dataframes. Use .reset_index(level=0, inplace=True) to get the proteins as a column
		exportData(minProteinDF.reset_index(level=0, inplace=False), dataType='df', path_out=params['path_out'],
				   filename=params['jobName'] + '_results_minimal', delim_out=params['delim_out'])
		exportData(fullProteinDF.reset_index(level=0, inplace=False), dataType='df', path_out=params['path_out'],
				   filename=params['jobName'] + '_results_full', delim_out=params['delim_out'])
		# save the intensity matrix of all COMMON peptides
		exportData(allExperimentsIntensitiesPerCommonPeptide, dataType='df', path_out=params['path_out'],
				   filename=params['jobName'] + '_CommonPeptideIntensities',
				   delim_out=params['delim_out'], inOneFile=False)
		# save the metadata
		exportData(metadata, dataType='df', path_out=params['path_out'],
				   filename=params['jobName'] + '_metadata',
				   delim_out=params['delim_out'], inOneFile=False)

	return minProteinDF, fullProteinDF, PCAResult, HCResult, allExperimentsIntensitiesPerCommonPeptide, metadata
