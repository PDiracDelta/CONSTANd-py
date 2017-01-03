#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Workflow of the processing part of CONSTANd++.
"""

from processing import *
from dataIO import exportData
from collapse import collapse
from constand import constand


def processDf(df, params, writeToDisk):
	# todo docu
	removedData = {}  # is to contain basic info about data that will be removed during the workflow, per removal category.
	# remove detections where (essential) data is missing.
	df, removedData['missing'] = removeMissing(df, params['noMissingValuesColumns'], params['intensityColumns'])
	if params['removeBadConfidence_bool']:
		df, removedData['confidence'] = removeBadConfidence(df, params['removeBadConfidence_minimum'], params['removalColumnsToSave'])
	# remove all useless columns from the dataFrame
	df = removeObsoleteColumns(df, wantedColumns=params['wantedColumns'])
	if params['removeIsolationInterference_bool']:
		# remove all data with too high isolation interference
		df, removedData['isolationInterference'] = removeIsolationInterference(df, params[
			'removeIsolationInterference_threshold'], params['removalColumnsToSave'])
	# remove all non-master protein accessions (entire column) and descriptions (selective).
	df = setMasterProteinDescriptions(df)
	if params['undoublePSMAlgo_bool']:
		# collapse peptide list redundancy due to overlap in MASCOT/SEQUEST peptide matches
		df, removedData['PSMAlgo'] = undoublePSMAlgo(df, identifyingNodes=params['identifyingNodes'],
		                                             exclusive=params['undoublePSMAlgo_exclusive_bool'],
		                                             intensityColumns=params['intensityColumns'],
		                                             removalColumnsToSave=params['removalColumnsToSave'])
		# SANITY CHECK: no detections with the same scan number may exist after undoublePSMAlgo()
		assert np.prod((len(i) < 2 for (s, i) in df.groupby('First Scan').groups))
	# collapse peptide list redundancy due to multiple detections at different RT
	# TEST here the intensity columns are alraedy lost
	df, removedData['RT'] = collapse('RT', df, intensityColumns=params['intensityColumns'], method=params['collapse_method'],
	                                 identifyingNodes=params['identifyingNodes'],
	                                 undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'], columnsToSave=params['collapseColumnsToSave'])
	if params['collapseCharge_bool']:
		# collapse peptide list redundancy due to different charges (optional)
		df, removedData['charge'] = collapse('Charge', df, intensityColumns=params['intensityColumns'], method=params['collapse_method'],
		                                     identifyingNodes=params['identifyingNodes'],
		                                     undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'], columnsToSave=params['collapseColumnsToSave'])
	if params['collapsePTM_bool']:
		# collapse peptide list redundancy due to different charges (optional)
		df, removedData['modifications'] = collapse('PTM', df, intensityColumns=params['intensityColumns'], method=params['collapse_method'],
		                                            identifyingNodes=params['identifyingNodes'],
		                                            undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'], columnsToSave=params['collapseColumnsToSave'])

	# SANITY CHECK: there should be no more duplicates if all collapses have been applied.
	if params['undoublePSMAlgo_bool'] and params['collapseCharge_bool']:  # TEST
		assert np.prod((len(i) < 2 for (s, i) in df.groupby(
			'Annotated Sequence').groups))  # only 1 index vector in dict of SEQUENCE:[INDICES] for all sequences

	if params['isotopicCorrection_bool']:
		# perform isotopic corrections but do NOT apply them to df because this information is sensitive (copyright i-TRAQ)
		intensities, noCorrectionIndices = isotopicCorrection(getIntensities(df, intensityColumns=params['intensityColumns']),
		                                                      correctionsMatrix=params['isotopicCorrection_matrix'])
	else:
		intensities = getIntensities(df, intensityColumns=params['intensityColumns'])

	doConstand = True # todo # TEST
	if doConstand:
		# perform the CONSTANd algorithm;
		normalizedIntensities, convergenceTrail, R, S = constand(intensities, params['accuracy'], params['maxIterations'])
		normalizedDf = setIntensities(df, intensities=normalizedIntensities, intensityColumns=params['intensityColumns'])
	else:
		# TEST do NOT perform CONSTANd
		logging.warning("+++++++++++++++++++++++++++++++CONSTAND NOT PERFORMED+++++++++++++++++++++++++++++++++")
		normalizedIntensities=[0]
		normalizedDf = df

	""" save results """
	if writeToDisk:
		# save the removed data information
		exportData(removedData, dataType='df', path_out=params['path_out'], filename=params['filename_out'] + '_removedData',
		           delim_out=params['delim_out'], inOneFile=params['removedDataInOneFile_bool'])
		# save the final form of the dataFrame WITHOUT normalized intensities.
		exportData(df, dataType='df', path_out=params['path_out'], filename=params['filename_out'] + '_dataFrame', delim_out=params['delim_out'])
		# save the normalized intensities obtained through CONSTANd
		exportData(normalizedIntensities, dataType='txt', path_out=params['path_out'],
		           filename=params['filename_out'] + '_normalizedIntensities', delim_out=params['delim_out'])
		# save the DE analysis results

	if params['isotopicCorrection_bool']:
		return normalizedDf, normalizedIntensities, removedData, noCorrectionIndices # todo find better solution than 2 returns
	else:
		return normalizedDf, normalizedIntensities, removedData
