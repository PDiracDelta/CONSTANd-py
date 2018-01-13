#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Workflow of the analysis part of QCQuan.
"""

from qcquan.analysis import *
from qcquan.dataIO import exportData
from qcquan.DEA import DEA


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
	:return commonPeptidesQuanValuesDF:	pd.DataFrame	for all COMMON (found in all experiments) peptides:
																		[e1_channel1, e1_channel2, ..., eM_channel1, ..., eM_channelN]
	:return metadata:			pd.DataFrame	metadata: [noIsotopicCorrection, RTIsolationInfo, noMasterProteinAccession,
												minSingleConditionProteins, fullSingleConditionProteins, uncommonModifiedPeptides, commonNanValues]
	:return extraOutputFiles:	list			list of full paths of extra output files to be included in the .zip file
	"""
	
	""" Preparation and metadata gathering """
	processingResultsItems = processingResults.items()
	dfs = dict((eName, result[0]) for eName, result in processingResultsItems)
	constandOutputs = dict((eName, result[1]) for eName, result in processingResultsItems)
	removedDatas = dict((eName, result[2]) for eName, result in processingResultsItems)
	metadatas = dict((eName, result[3]) for eName, result in processingResultsItems)
	processedDfFullPaths = dict((eName, result[4]) for eName, result in processingResultsItems)
	
	# contains statistics and metadata (like the parameters) about the analysis.
	metadata = dict()
	# metadataframe for all simple numeric values
	metadata['numeric'] = pd.DataFrame()
	
	# combine all metadata from each separate MS run
	metadata, allObservedProteins = combineProcessingMetadata(metadata, metadatas)
	# get MS1 intensities on the peptide level, i.e. after aggregation and cleaning.
	metadata['MS1Intensities_peptides'] = pd.Series(index=list(dfs.keys()), dtype=object)
	for eName in dfs.keys():
		# reset index because otherwise the df will get NaN values since not all MS1 intensity indices are equal
		# across all experiments
		metadata['MS1Intensities_peptides'][eName] = dfs[eName].loc[:, 'Intensity'].tolist()
	
	metadata['numeric'].loc[0, 'numObservedProteins'] = len(allObservedProteins)
	metadata['allObservedProteins'] = pd.DataFrame({'protein': allObservedProteins})
	
	# record RT isolation statistics. Future: flag. Multi-indexed on experiment names and old indices!
	if params['getRTIsolationInfo_bool']:
		experimentNames = processingResults.keys()
		metadata['RTIsolationInfo'] = pd.concat([getRTIsolationInfo(removedDatas[eName]['RT']) for
												 eName in experimentNames], keys=experimentNames)
	
	""" Differential Expression Analysis """
	# merge all experiments in multi-indexed: (eName, oldIndex) dataframe as an outer join
	allExperimentsDF = combineExperimentDFs(dfs)

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
			DEA(allExperimentsDF, minProteinPeptidesDict, params)
	else:
		minProteinDF = pd.DataFrame()
	
	metadata['numeric'].loc[0, 'numUniqueModifiedPeptides'] = len(unnest(minProteinDF.loc[:, 'peptides']))
	
	if params['fullExpression_bool']:
		# Bring the data to the protein level in the case of full expression (shared peptides allowed).
		# Execute the differential expression analysis and gather some metadata
		fullProteinDF, metadata['fullSingleConditionProteins'], metadata['numeric'].loc[0, 'fullNumProteins'] = \
			DEA(allExperimentsDF, fullProteinPeptidesDict, params)
	else:
		fullProteinDF = pd.DataFrame()
	
	# set the protein names back as columns instead of the index, and sort the columns so the df is easier to read
	handyColumnOrder = buildHandyColumnOrder(minProteinDF.columns, params['referenceCondition'], params['schema'])
	minProteinDF = minProteinDF.reindex_axis(handyColumnOrder, axis=1)
	fullProteinDF = fullProteinDF.reindex_axis(handyColumnOrder, axis=1)
	
	""" Quality Control """
	# dataframe with ALL intensities per peptide: [peptide, e1_channel1, e1_channel2, ..., eM_channel1, ..., eM_channelN]
	commonPeptidesQuanValuesDF, metadata['uncommonModifiedPeptides'] = getCommonPeptidesQuanValuesDF(dfs, params['schema'])
	metadata['numeric'].loc[0, 'numUnCommonModifiedPeptides'] = len(metadata['uncommonModifiedPeptides'])
	metadata['numeric'].loc[0, 'numCommonModifiedPeptides'] = len(commonPeptidesQuanValuesDF)
	metadata['numeric'].loc[0, 'numModifiedPeptides'] = \
		metadata['numeric'].loc[0, 'numUnCommonModifiedPeptides'] + metadata['numeric'].loc[0, 'numCommonModifiedPeptides']
	# save the amount of NaN values per channel for common peptides.
	metadata['commonNanValues'] = pd.DataFrame(np.sum(np.isnan(commonPeptidesQuanValuesDF), axis=0))
	# metadata['numUnCommonModifiedPeptidesPerCondition'] = getNumUnCommonModifiedPeptidesPerCondition(dfs, metadata['uncommonModifiedPeptides'])
	# perform PCA
	PCAResult = getPCA(commonPeptidesQuanValuesDF, params['PCA_components'])
	# perform hierarchical clustering
	HCResult = getHC(commonPeptidesQuanValuesDF)

	""" save results """
	if writeToDisk:
		# save the protein-level dataframes. Use .reset_index(level=0, inplace=True) to get the proteins as a column
		minProteinDF_fullPath = exportData(minProteinDF.reset_index(level=0, inplace=False), dataType='df', path_out=params['path_out'],
				   filename=params['jobName'] + '_results_minimal', delim_out=params['delim_out'])
		fullProteinDF_fullPath = exportData(fullProteinDF.reset_index(level=0, inplace=False), dataType='df', path_out=params['path_out'],
				   filename=params['jobName'] + '_results_full', delim_out=params['delim_out'])
		# save the intensity matrix of all COMMON peptides
		exportData(commonPeptidesQuanValuesDF, dataType='df', path_out=params['path_out'],
				   filename=params['jobName'] + '_CommonPeptideIntensities',
				   delim_out=params['delim_out'], inOneFile=False)
		# save the metadata
		exportData(metadata, dataType='df', path_out=params['path_out'],
				   filename=params['jobName'] + '_metadata',
				   delim_out=params['delim_out'], inOneFile=False)
	else:
		minProteinDF_fullPath = None
		fullProteinDF_fullPath = None
	extraOutputFiles = list(processedDfFullPaths.values()) + [minProteinDF_fullPath]
	if params['fullExpression_bool']:
		extraOutputFiles.append(fullProteinDF_fullPath)

	return minProteinDF, fullProteinDF, PCAResult, HCResult, commonPeptidesQuanValuesDF, metadata, extraOutputFiles
