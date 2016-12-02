#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Python implementation of mass spectrometer protein data analysis using the CONSTANd_RAS algorithm.
"""

__author__ = "Joris Van Houtven"
__copyright__ = "Copyright ?YEAR, VITO"
__credits__ = ["Joris Van Houtven", "Dirk Valkenborg"]
# __license__ = "GPL"
# __version__ = "0.0.1"
__maintainer__ = "Joris Van Houtven"
__email__ = "vanhoutvenjoris@gmail.com"
__status__ = "Development"

import sys
# import matplotlib as mpl
# import matplotlib.pyplot as plt
from constand import constand
from os import path
from time import time
from dataIO import *
from dataproc import *
from collapse import collapse, setCollapseColumnsToSave
from analysis import *
from report import *


def performanceTest():  # remove for production # TEST
	""" Use this development method to test the performance of the CONSTANd algorithm. """
	t = []
	for i in range(100):
		params = getInput()
		df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(2*10**3, 6)), columns=list('ABCDEF'))
		start = time()
		constand(np.asarray(df), 1e-2, 50)
		stop = time()
		t.append((stop - start))
	print("average runtime: " + str(np.mean(t)))


def isotopicImpuritiesTest(): # TEST
	## test if isotopic correction is necessary:
	params = getInput()
	# get the dataframe
	df = importDataFrame(params['files_in'], delim=params['delim_in'], header=params['header_in'])
	correctedIntensities = getIntensities(df)
	normalizedIntensities, convergenceTrail, R, S = constand(correctedIntensities, params['accuracy'],
	                                                         params['maxIterations'])
	# exportData(normalizedIntensities, 'txt', path_out=params['path_out'],
	#            filename=params['filename_out'] + '_normalizedIntensities', delim_out=params['delim_out'])
	# test "impure data"
	correctedIntensities_impure = correctedIntensities
	spillover = correctedIntensities_impure[0, :] * 0.1
	correctedIntensities_impure[0, :] -= spillover
	correctedIntensities_impure[1, :] += spillover
	normalizedIntensities_impure, convergenceTrail_i, R_i, S_i = constand(correctedIntensities_impure,
	                                                                      params['accuracy'],
	                                                                      params['maxIterations'])
	diff = abs(normalizedIntensities - normalizedIntensities_impure)
	print(np.allclose(normalizedIntensities, normalizedIntensities_impure, atol=1e-3, equal_nan=True))
	print(np.nanmean(np.nanmean(diff[:, 0:1], 1)))
	print(max(np.amax(diff, 1)))
# False tot en met 1e-3 --> fouten van > 0.1%


def isotopicCorrectionsTest(params): # TEST
	if params['isotopicCorrection_bool']:
		int_in = np.array([range(6), range(6)]) + np.array([np.zeros(6), 5*np.ones(6)])
		# perform isotopic corrections but do NOT apply them to df because this information is sensitive (copyright i-TRAQ)
		icm = params['isotopicCorrection_matrix']
		icm[0,0] = 0.9
		icm[0,1] = 0.1
		int_out = isotopicCorrection(int_in, correctionsMatrix=icm)
		print(int_out)
		# M=np.eye(6); M[0,0]=0.9; M[0,1]=0.1; b=np.asarray(range(6)); c=np.asarray(range(6))+5
		# print(int_out) above should be equal to:
		# [np.linalg.solve(M, b) ; np.linalg.solve(M, c)]


def MS2IntensityDoesntMatter(df):
	I = getIntensities(df)
	r1 = constand(I, 1e-5, 50)
	I[0] *= 1e9 # this is BIG. MS2 intensity doesnt reach beyond 1e9 so if one value has magnitudeOrder 1 it's still OK.
	r2 = constand(I, 1e-5, 50)
	print(np.allclose(r1[0], r2[0], equal_nan=True))
	diff = r1[0] - r2[0]
	maxdiff = max(np.amax(diff, 1))
	print(maxdiff)


def testDataComplementarity(df):
	scannrs_init = set(df.groupby('First Scan').groups.keys())
	main(testing=False, writeToDisk=True)
	# SANITY CHECK if df + removedData scan numbers = total scan numbers.
	scannrs_final = set(df.groupby('First Scan').groups.keys())
	##### THIS IS OUTDATED SINCE COMMIT b98041f
	removedDataLoaded = pickle.load(open('../data/MB_result_removedData', 'rb'))
	for value in removedDataLoaded.values():
		scannrs_final = scannrs_final.union(set(value['First Scan']))
	print(scannrs_final == scannrs_init)


def compareIntensitySN():
	filepath1 = '../data/COON data/PSMs/BR1_a.txt'
	filepath2 = '../data/COON data/PSMs/BR1_b_SN.txt'
	df1 = importDataFrame(filepath1, delim='\t', header=0)
	df2 = importDataFrame(filepath2, delim='\t', header=0)
	intensityColumns = ["126", "127N", "127C", "128C", "129N", "129C", "130C", "131"]
	relIntensities = np.empty((len(df1.index),8), dtype='float')
	relSNs = np.empty((len(df2.index),8), dtype='float')
	for col in range(len(intensityColumns)):
		relIntensities[:, col] = 1/8*df1.loc[:, intensityColumns[col]]/np.nanmean(df1.loc[:, intensityColumns],1)
		relSNs[:, col] = 1/8*df2.loc[:, intensityColumns[col]]/np.nanmean(df2.loc[:, intensityColumns],1)
	diff = abs(relIntensities - relSNs)
	print(np.allclose(relIntensities, relSNs, atol=1e-3, equal_nan=True))
	print("mean over all values")
	print(np.nanmean(np.nanmean(diff[:, 0:6], 1)))
	print("max difference")
	print(np.nanmax(np.nanmax(diff, 1)))


def devStuff(df, params): # TEST
	# performanceTest()
	# isotopicCorrectionsTest(params)
	# MS2IntensityDoesntMatter(df)
	# testDataComplementarity(df)
	compareIntensitySN()
	pass


def processDf(df, params, writeToDisk):
	removedData = {}  # is to contain basic info about data that will be removed during the workflow, per removal category.
	# remove detections where (essential) data is missing.
	df, removedData['missing'] = removeMissing(df)
	if params['removeBadConfidence_bool']:
		df, removedData['confidence'] = removeBadConfidence(df, params['removeBadConfidence_minimum'])
	# remove all useless columns from the dataFrame
	df = removeObsoleteColumns(df, wantedColumns=params['wantedColumns'])
	if params['removeIsolationInterference_bool']:
		# remove all data with too high isolation interference
		df, removedData['isolationInterference'] = removeIsolationInterference(df, params[
			'removeIsolationInterference_threshold'])
	# remove all non-master protein accessions (entire column) and descriptions (selective).
	df = setMasterProteinDescriptions(df)
	if params['undoublePSMAlgo_bool']:
		# collapse peptide list redundancy due to overlap in MASCOT/SEQUEST peptide matches
		df, removedData['PSMAlgo'] = undoublePSMAlgo(df, master=params['masterPSMAlgo'],
		                                             exclusive=params['undoublePSMAlgo_exclusive_bool'])
		# SANITY CHECK: no detections with the same scan number may exist after undoublePSMAlgo()
		assert np.prod((len(i) < 2 for (s, i) in df.groupby('First Scan').groups))

	# collapse peptide list redundancy due to multiple detections at different RT
	# TEST here the intensity columns are alraedy lost
	df, removedData['RT'] = collapse('RT', df, method=params['collapse_method'],
	                                 maxRelativeReporterVariance=params['collapse_maxRelativeReporterVariance'],
	                                 masterPSMAlgo=params['masterPSMAlgo'],
	                                 undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'])
	if params['collapseCharge_bool']:
		# collapse peptide list redundancy due to different charges (optional)
		df, removedData['charge'] = collapse('Charge', df, method=params['collapse_method'],
		                                     maxRelativeReporterVariance=params['collapse_maxRelativeReporterVariance'],
		                                     masterPSMAlgo=params['masterPSMAlgo'],
		                                     undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'])
	if params['collapsePTM_bool']:
		# collapse peptide list redundancy due to different charges (optional)
		df, removedData['modifications'] = collapse('PTM', df, method=params['collapse_method'],
		                                            maxRelativeReporterVariance=params[
			                                            'collapse_maxRelativeReporterVariance'],
		                                            masterPSMAlgo=params['masterPSMAlgo'],
		                                            undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'])

	# SANITY CHECK: there should be no more duplicates if all collapses have been applied.
	if params['undoublePSMAlgo_bool'] and params['collapseCharge_bool']:  # TEST
		assert np.prod((len(i) < 2 for (s, i) in df.groupby(
			'Annotated Sequence').groups))  # only 1 index vector in dict of SEQUENCE:[INDICES] for all sequences

	if params['isotopicCorrection_bool']:
		# perform isotopic corrections but do NOT apply them to df because this information is sensitive (copyright i-TRAQ)
		intensities, noCorrectionIndices = isotopicCorrection(getIntensities(df),
		                                                      correctionsMatrix=params['isotopicCorrection_matrix'])
	else:
		intensities = getIntensities(df)
	# perform the CONSTANd algorithm; also do NOT include normalized intensities in df --> only for paying users.
	normalizedIntensities, convergenceTrail, R, S = constand(intensities, params['accuracy'], params['maxIterations'])
	normalizedDf = setIntensities(df, normalizedIntensities)

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

	return normalizedDf, normalizedIntensities, removedData, noCorrectionIndices


def analyzeProcessingResult(processingResults, params, writeToDisk):
	dfs = [result[0] for result in processingResults]
	normalizedIntensitiess = [result[1] for result in processingResults]
	removedDatas = [result[2] for result in processingResults]
	noCorrectionIndicess = [result[3] for result in processingResults]

	# TODO effectively implement multiple experiment analysis beyond this point
	normalizedDf = dfs[0]
	normalizedIntensities = normalizedIntensitiess[0]
	removedData = removedDatas[0]
	noCorrectionIndices = noCorrectionIndicess[0]

	# contains statistics and metadata (like the parameters) about the analysis.
	metadata = {}
	# record detections without isotopic correction applied applied
	metadata['noIsotopicCorrection'] = getNoIsotopicCorrection(normalizedDf, noCorrectionIndices)
	# record RT isolation statistics. Future: flag
	metadata['RTIsolationInfo'] = getRTIsolationInfo(removedData['RT'])

	# get min and max protein-peptide mappings
	minProteinPeptidesDict, maxProteinPeptidesDict, metadata['noMasterProteinAccession'] = getProteinPeptidesDicts(normalizedDf)
	# execute mappings to get all peptideintensities per protein, over each whole condition
	minProteinDF = getProteinDF(normalizedDf, minProteinPeptidesDict, params['intensityColumnsPerCondition'])
	fullProteinDF = getProteinDF(normalizedDf, maxProteinPeptidesDict, params['intensityColumnsPerCondition'])

	# perform differential expression analysis with Benjamini-Hochberg correction.
	minProteinDF = applyDifferentialExpression(minProteinDF, params['alpha'])
	fullProteinDF = applyDifferentialExpression(fullProteinDF, params['alpha'])

	# calculate fold changes of the average protein expression value per CONDITION/GROUP (not per channel!)
	minProteinDF = applyFoldChange(minProteinDF, params['pept2protCombinationMethod'])
	fullProteinDF = applyFoldChange(fullProteinDF, params['pept2protCombinationMethod'])

	# perform PCA
	PCAResult = getPCA(getIntensities(normalizedDf), params['PCA_components'])

	# perform hierarchical clustering
	HCResult = getHC(getIntensities(normalizedDf))

	# set the protein names back as columns instead of the index, and sort the columns so the df is easier to read
	handyColumnOrder = ['protein', 'adjusted p-value', 'fold change c1/c2', 'p-value', 'peptides', 'condition 1', 'condition 2']
	minProteinDF.reset_index(level=0, inplace=True)
	fullProteinDF.reset_index(level=0, inplace=True)
	minProteinDF = minProteinDF.reindex_axis(handyColumnOrder, axis=1)
	fullProteinDF = fullProteinDF.reindex_axis(handyColumnOrder, axis=1)

	""" save results """
	if writeToDisk:
		# save the protein-level dataframes
		exportData(minProteinDF, dataType='df', path_out=params['path_out'],
		           filename=params['filename_out'] + '_results_minimal', delim_out=params['delim_out'])
		exportData(fullProteinDF, dataType='df', path_out=params['path_out'],
		           filename=params['filename_out'] + '_results_full', delim_out=params['delim_out'])
		# save the metadata
		exportData(metadata, dataType='df', path_out=params['path_out'],
		           filename=params['filename_out'] + '_metadata',
		           delim_out=params['delim_out'], inOneFile=False)
		# generate a report PDF (without the normalized intensities: behind paywall?

	return minProteinDF, fullProteinDF, PCAResult, HCResult, metadata


def generateReport(analysisResults, params, writeToDisk):
	# todo docu
	# todo multi-experiment support
	minProteinDF = analysisResults[0]
	fullProteinDF = analysisResults[1]
	PCAResult = analysisResults[2]
	HCResult = analysisResults[3]
	metadata = analysisResults[4]

	# data visualization
	visualizationsDict = dataVisualization(minProteinDF, fullProteinDF, params['alpha'], params['FCThreshold'],
	                                       PCAResult, HCResult, params['intensityColumnsPerCondition'])
	if writeToDisk:
		# save the visualizations
		exportData(visualizationsDict, dataType='visualizationsDict', path_out=params['path_out'],
		           filename=params['filename_out'] + '_dataViz')  # TODO


def main(doProcessing, doAnalysis, doReport, writeToDisk, testing):
	start = time()
	# testing=True # TEST
	# writeToDisk=False # TEST
	"""
	For now this is just stuff for debugging and testing. Later:
	Contains and explicits the workflow of the program. Using the booleans doProcessing, doAnalysis and writeToDisk one
	can control	which parts of the workflow to perform.
	"""
	""" get all input parameters
	params:
	files_in, delim_in, header_in, intensityColumns, wantedColumns, collapseColumnsToSave, undoublePSMAlgo_bool,
	removeIsolationInterference_bool, collapse_method, collapse_maxRelativeReporterVariance,
	removeIsolationInterference_master, masterPSMAlgo, undoublePSMAlgo_exclusive_bool, collapseRT_bool,
	collapseCharge_bool, collapsePTM_bool, isotopicCorrection_bool, isotopicCorrection_matrix, accuracy, maxIterations,
	DEFoldThreshold, path_out, filename_out, delim_out
	"""
	params = getInput()
	# get the dataframes
	dfs = []
	for filepath in params['files_in']:
			dfs.append(importDataFrame(filepath, delim=params['delim_in'], header=params['header_in']))

	# define global parameters
	setProcessingGlobals(intensityColumns=params['intensityColumns'],
	                     removalColumnsToSave=params['removalColumnsToSave'],
	                     noMissingValuesColumns=params['noMissingValuesColumns'])
	setCollapseColumnsToSave(params['collapseColumnsToSave'])  # define the intensityColumns for use in dataproc.py

	if not testing:
		""" Data processing """
		processingResultsDumpFilename = path.relpath(path.join(filepath, path.pardir))+'/processingResultsDump'
		if doProcessing:
			# process every input dataframe
			processingResults = [processDf(df, params, writeToDisk) for df in dfs]
			pickle.dump(processingResults, open(processingResultsDumpFilename, 'wb')) # TEST
		elif doAnalysis:
			try:
				processingResults = pickle.load(open(processingResultsDumpFilename, 'rb'))
			except FileNotFoundError:
				raise FileNotFoundError("There is no previously processed data in this path: "+processingResultsDumpFilename)
		else:
			warn("No processing step performed nor processing file loaded!")

		""" Data analysis and visualization """
		analysisResultsDumpFilename = path.relpath(path.join(filepath, path.pardir)) + '/analysisResultsDump'
		if doAnalysis:
			analysisResults = analyzeProcessingResult(processingResults, params, writeToDisk)
			pickle.dump(analysisResults, open(analysisResultsDumpFilename, 'wb'))  # TEST
		elif doReport:
			try:
				analysisResults = pickle.load(open(analysisResultsDumpFilename, 'rb'))
			except FileNotFoundError:
				raise FileNotFoundError("There is no previously analyzed data in this path: "+analysisResultsDumpFilename)
		else:
			warn("No analysis step performed nor analysis file loaded!")

		""" generate report """
		if doReport:
			generateReport(analysisResults, params, writeToDisk)  # TODO
		else:
			warn("No report generated!")

	elif testing:
		devStuff(dfs[0], params)
	stop = time()
	print(stop - start)


if __name__ == '__main__':
	sys.exit(main(doProcessing=True, doAnalysis=True, doReport=True, testing=False, writeToDisk=True))
