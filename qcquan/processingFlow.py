#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Workflow of the processing part of CONSTANd++.
"""

from constandpp.processing import *
from constandpp.tools import setIntensities
from constandpp.dataIO import exportData
from constandpp.aggregate import aggregate
from constandpp.constand import constand


def processDf(df, params, writeToDisk, doConstand=True):
	"""
	Calls all the necessary functions to process the dataframe of one experiment and prepare the analysis input objects.
	Cleans the input data, removes redundancy due to PSM algorithm, charge (optional) and PTMs (optional), then corrects
	isotopic impurities (optional),	normalizes using the CONSTANd method and saves the output to disk if so required.
	Along the way, removed data is kept in a corresponding dict of dataframes.
	:param df:						pd.DataFrame		Peptide Spectrum Match dataframe (see documentation).
	:param params:					dict				experiment-specific processing parameters (see getConfig.py.)
	:param writeToDisk:				bool				write results to harddisk (if not: only pass via return statement).
	:return normalizedDf:			pd.DataFrame		cleaned, normalized data on the "unique modified peptide" level
	:return normalizedIntensities:	np.ndarray			matrix of quantification values on the "unique modified peptide" level
	:return removedData:			{ pd.DataFrame }	removed data: [missing, confidence, isolationInterference,
														PSMAlgo, RT, charge, modifications]
	"""
	removedData = {}  # is to contain basic info about data that will be removed during the workflow, per removal category.
	
	# remove all useless columns from the dataFrame
	df = removeObsoleteColumns(df, wantedColumns=params['wantedColumns'])
	
	# get a set of all master proteins detected in at least one PSM.
	allMasterProteins = getAllPresentProteins(df)
	
	# remove PSMs where (essential) data is missing.
	df, removedData['missing'] = removeMissing(df, params['noMissingValuesColumns'], params['quanColumns'], params['identifyingNodes'])
	
	if params['removeBadConfidence_bool']:
		df, removedData['confidence'] = removeBadConfidence(df, params['removeBadConfidence_minimum'], params['removalColumnsToSave'])
	
	if params['removeIsolationInterference_bool']:
		# remove all data with too high isolation interference
		df, removedData['isolationInterference'] = removeIsolationInterference(df, params['removeIsolationInterference_threshold'],
																			   params['removalColumnsToSave'])
		
	# remove all non-master protein accessions (entire column) and descriptions (selective).
	df = setMasterProteinDescriptions(df)
	
	if params['undoublePSMAlgo_bool']:
		# aggregate peptide list redundancy due to overlap in MASCOT/SEQUEST peptide matches
		df, removedData['PSMAlgo'] = undoublePSMAlgo(df, identifyingNodes=params['identifyingNodes'],
													 exclusive=params['undoublePSMAlgo_exclusive_bool'],
													 quanColumns=params['quanColumns'],
													 removalColumnsToSave=params['removalColumnsToSave'])
		# SANITY CHECK: no PSMs with the same scan number may exist after undoublePSMAlgo()
		assert np.prod((len(i) < 2 for (s, i) in df.groupby('First Scan').groups))
	
	if params['isotopicCorrection_bool']:
		# perform isotopic corrections and then apply them to the dataframe. No i-TRAQ copyright issues as of 2017-05-04
		# (see logboek.txt or git history)
		correctedIntensities, noCorrectionIndices = isotopicCorrection(getIntensities(df, quanColumns=params['quanColumns']),
															  correctionsMatrix=params['isotopicCorrection_matrix'])
		df = setIntensities(df, intensities=correctedIntensities, quanColumns=params['quanColumns'])
		# TEST remove the negative quan value rows
		# from scripts.tools import removeRowsWithNeg
		# df = removeRowsWithNeg(df, params['quanColumns'])
	
	# aggregate peptide list redundancy due to multiple PSMs at different RT
	df, removedData['RT'] = aggregate('RT', df, quanColumns=params['quanColumns'], method=params['aggregate_method'],
									 identifyingNodes=params['identifyingNodes'],
									 undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'], columnsToSave=params['aggregateColumnsToSave'])
	
	if params['aggregateCharge_bool']:
		# aggregate peptide list redundancy due to different charges (optional)
		df, removedData['charge'] = aggregate('Charge', df, quanColumns=params['quanColumns'], method=params['aggregate_method'],
											 identifyingNodes=params['identifyingNodes'],
											 undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'], columnsToSave=params['aggregateColumnsToSave'])
	
	if params['aggregatePTM_bool']:
		# aggregate peptide list redundancy due to different charges (optional)
		df, removedData['modifications'] = aggregate('PTM', df, quanColumns=params['quanColumns'], method=params['aggregate_method'],
													identifyingNodes=params['identifyingNodes'],
													undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'], columnsToSave=params['aggregateColumnsToSave'])

	# SANITY CHECK: there should be no more duplicates if all aggregates have been applied.
	if params['undoublePSMAlgo_bool'] and params['aggregateCharge_bool']:  # TEST
		assert np.prod((len(i) < 2 for (s, i) in df.groupby(
			'Annotated Sequence').groups))  # only 1 index vector in dict of SEQUENCE:[INDICES] for all sequences

	if doConstand:
		# perform the CONSTANd algorithm;
		intensities = getIntensities(df, quanColumns=params['quanColumns'])
		constandOutput = constand(intensities, params['precision'], params['maxIterations'])
		normalizedIntensities = constandOutput[0]
		normalizedDf = setIntensities(df, intensities=normalizedIntensities, quanColumns=params['quanColumns'])
	else:
		# TEST do NOT perform CONSTANd
		logging.warning("+++++++++++++++++++++++++++++++CONSTAND NOT PERFORMED+++++++++++++++++++++++++++++++++")
		normalizedIntensities = [0]
		normalizedDf = df

	""" save results """
	if writeToDisk:
		import os
		# save the removed data information
		# create folder first
		removedDataDir = os.path.join(params['path_out'], 'removedData')
		os.makedirs(removedDataDir)
		exportData(removedData, dataType='df', path_out=removedDataDir, filename=params['filename_out'] + '_removedData',
				   delim_out=params['delim_out'], inOneFile=params['removedDataInOneFile_bool'])
		# save the final form of the dataFrame WITHOUT normalized intensities.
		exportData(df, dataType='df', path_out=params['path_out'], filename=params['filename_out'] + '_dataFrame', delim_out=params['delim_out'])
		# save the normalized intensities obtained through CONSTANd
		exportData(normalizedIntensities, dataType='txt', path_out=params['path_out'],
				   filename=params['filename_out'] + '_normalizedIntensities', delim_out=params['delim_out'])
		# save the DE analysis results

	if params['isotopicCorrection_bool']:
		return normalizedDf, constandOutput, removedData, allMasterProteins, noCorrectionIndices  # todo find better solution than 2 returns
	else:
		return normalizedDf, constandOutput, removedData, allMasterProteins  #todo add noCorrectionIndices variable that is empty