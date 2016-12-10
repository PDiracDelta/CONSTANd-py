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

import sys, os, logging
from webFlow import webFlow
from getInput import getInput, getMasterInput
from constand import constand
from time import time
from dataIO import *
from dataproc import *
from collapse import collapse
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


def MAPlot(x,y):
	from matplotlib import pyplot as plt
	logx = np.log2(x)
	logy = np.log2(y)
	plt.scatter((logx+logy)*0.5, logx-logy)
	plt.title('PD2.1 Intensities versus S/N values (scaled relatively within each row/peptide)')
	plt.xlabel('A')
	plt.ylabel('M')
	plt.show()


def compareIntensitySN():
	filepath1 = '../data/COON data/PSMs/fixed_BR1_a.txt'
	filepath2 = '../data/COON data/PSMs/fixed_BR1_b_SN.txt'
	constandnorm=True
	if constandnorm:
		if path.exists('../data/compareIntensitySNProcessingResults'):
			processingResults = pickle.load(open('../data/compareIntensitySNProcessingResults', 'rb'))
		else:
			params=getInput()
			# setProcessingGlobals(intensityColumns=params['intensityColumns'],
			# 					 removalColumnsToSave=params['removalColumnsToSave'],
			# 					 noMissingValuesColumns=params['noMissingValuesColumns'])
			# setCollapseColumnsToSave(params['collapseColumnsToSave'])  # define the intensityColumns for use in dataproc.py
			dfs = []
			for filepath in [filepath1, filepath2]:
				dfs.append(importDataFrame(filepath, delim=params['delim_in'], header=params['header_in']))
			processingResults = [processDf(df, params, writeToDisk=False) for df in dfs]
			pickle.dump(processingResults, open('../data/compareIntensitySNProcessingResults', 'wb'))
		relIntensities = getIntensities(processingResults[0][0])
		relSNs = getIntensities(processingResults[1][0])
	else:
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

	# TODO nakijken of hier geen fout gebeurt
	MAPlot(relIntensities.reshape(relIntensities.size, 1), relSNs.reshape(relSNs.size, 1))


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
	                                 maxRelativeReporterVariance=params['collapse_maxRelativeReporterVariance'],
	                                 identifyingNodes=params['identifyingNodes'],
	                                 undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'], columnsToSave=params['collapseColumnsToSave'])
	if params['collapseCharge_bool']:
		# collapse peptide list redundancy due to different charges (optional)
		df, removedData['charge'] = collapse('Charge', df, intensityColumns=params['intensityColumns'], method=params['collapse_method'],
		                                     maxRelativeReporterVariance=params['collapse_maxRelativeReporterVariance'],
		                                     identifyingNodes=params['identifyingNodes'],
		                                     undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'], columnsToSave=params['collapseColumnsToSave'])
	if params['collapsePTM_bool']:
		# collapse peptide list redundancy due to different charges (optional)
		df, removedData['modifications'] = collapse('PTM', df, intensityColumns=params['intensityColumns'], method=params['collapse_method'],
		                                            maxRelativeReporterVariance=params[
			                                            'collapse_maxRelativeReporterVariance'],
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
		intensities = getIntensities(df)
	# perform the CONSTANd algorithm; also do NOT include normalized intensities in df --> only for paying users.
	normalizedIntensities, convergenceTrail, R, S = constand(intensities, params['accuracy'], params['maxIterations'])
	normalizedDf = setIntensities(df, intensities=normalizedIntensities, intensityColumns=params['intensityColumns'])

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


def analyzeProcessingResult(processingResults, params, writeToDisk):
	processingResultsItems = processingResults.items()
	dfs = dict((eName, result[0]) for eName, result in processingResultsItems)
	normalizedIntensitiess = dict((eName, result[1]) for eName, result in processingResultsItems)
	removedDatas = dict((eName, result[2]) for eName, result in processingResultsItems)
	noCorrectionIndicess = dict((eName, result[3]) for eName, result in processingResultsItems)

	# normalizedDf = dfs[0]
	# normalizedIntensities = normalizedIntensitiess[0]
	# removedData = removedDatas[0]
	# noCorrectionIndices = noCorrectionIndicess[0]

	experimentNames = processingResults.keys()
	# contains statistics and metadata (like the parameters) about the analysis.
	metadata = {}
	# record detections without isotopic correction applied applied. Multi-indexed on experiment names and old indices!
	metadata['noIsotopicCorrection'] = pd.concat([getNoIsotopicCorrection(dfs[eName], noCorrectionIndicess[eName]) for
	                                              eName in experimentNames], keys=experimentNames)
	# record RT isolation statistics. Future: flag. Multi-indexed on experiment names and old indices!
	metadata['RTIsolationInfo'] = pd.concat([getRTIsolationInfo(removedDatas[eName]['RT']) for
	                                         eName in experimentNames], keys=experimentNames)

	# merge all experiments in multi-indexed: (eName, oldIndex) dataframe # and intensityColumns are unique and distinguishable
	allExperimentsDF = combineExperimentDFs(dfs) #, params['schema'])

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
	minProteinDF = analysisResults[0]
	fullProteinDF = analysisResults[1]
	PCAResult = analysisResults[2]
	HCResult = analysisResults[3]
	metadata = analysisResults[4]

	# generate sorted (on FC) list of differentials
	minSortedDifferentialProteinsDF = getSortedDifferentialProteinsDF(minProteinDF)
	fullSortedDifferentialProteinsDF = getSortedDifferentialProteinsDF(fullProteinDF)
	minSet = set(minSortedDifferentialProteinsDF['protein'])
	fullSet = set(fullSortedDifferentialProteinsDF['protein'])
	# list( [in min but not in full], [in full but not in min] )
	metadata['diffMinFullProteins'] = [list(minSet.difference(fullSet)), list(fullSet.difference(minSet))]
	# todo combine into one

	# data visualization
	visualizationsDict = {}
	visualizationsDict['minVolcano'] = getVolcanoPlot(minProteinDF, params['alpha'], params['FCThreshold'],
	                                               params['labelVolcanoPlotAreas'])
	visualizationsDict['fullVolcano'] = getVolcanoPlot(fullProteinDF, params['alpha'], params['FCThreshold'],
	                                                  params['labelVolcanoPlotAreas'])
	visualizationsDict['pca'] = getPCAPlot(PCAResult, params['schema']) # todo multiple experimenst
	visualizationsDict['hcd'] = getHCDendrogram(HCResult, params['schema']) # todo multiple experimenst

	# generate HTML and PDF reports # todo
	htmlReport = makeHTML(minSortedDifferentialProteinsDF, fullSortedDifferentialProteinsDF, visualizationsDict, metadata)
	pdfReport = HTMLtoPDF(htmlReport)

	writeToDisk = False # TEST
	if writeToDisk:
		# save the visualizations
		exportData(visualizationsDict, dataType='viz', path_out=params['path_results'],
		           filename=params['jobname'] + '_dataViz') # TODO


def main(masterConfigFilePath, doProcessing, doAnalysis, doReport, writeToDisk, testing):
	"""
	For now this is just stuff for debugging and testing. Later:
	Contains and explicits the workflow of the program. Using the booleans doProcessing, doAnalysis and writeToDisk one
	can control	which parts of the workflow to perform.
	"""
	logFilePath = os.path.relpath(os.path.join(masterConfigFilePath, os.path.join(os.pardir, 'log.txt')))
	logging.basicConfig(filename=logFilePath, level=logging.INFO)
	start = time()
	masterParams = getMasterInput(masterConfigFilePath) # config filenames + params for the combination of experiments
	specificParams = {} # specific params for each experiment
	dfs = {}
	processingResults = {}
	experimentNames = list(masterParams['schema'].keys())
	for eName in experimentNames:
		# get all input parameters
		specificParams[eName] = getInput(masterParams['schema'][eName]['config'])
		# get the dataframes
		dfs[eName] = getData(specificParams[eName]['data'], delim=specificParams[eName]['delim_in'], header=specificParams[eName]['header_in'], wrapper=specificParams[eName]['wrapper'])
	# todo find mean of empty slices warning flood source (ABOVE this line)
	if not testing:
		for eName in experimentNames:
			# define global parameters
			# setProcessingGlobals(intensityColumns=specificParams[eName]['intensityColumns'],
			#                      removalColumnsToSave=specificParams[eName]['removalColumnsToSave'],
			#                      noMissingValuesColumns=specificParams[eName]['noMissingValuesColumns'])
			# setCollapseColumnsToSave(
			# 	specificParams[eName]['collapseColumnsToSave'])  # define the intensityColumns for use in dataproc.py
			""" Data processing """
			processing_path_out = specificParams[eName]['path_out']
			processingResultsDumpFilename = path.relpath(path.join(processing_path_out, path.pardir))+'/processingResultsDump_'+str(eName)
			if doProcessing:
				# prepare the output directory
				assert not os.path.exists(processing_path_out)  # write to empty dir
				assert os.path.exists(
					path.relpath(path.join(processing_path_out, path.pardir)))  # parent dir must exist
				os.makedirs(processing_path_out)

				# process every input dataframe
				processingResults[eName] = processDf(dfs[eName], specificParams[eName], writeToDisk)
				pickle.dump(processingResults[eName], open(processingResultsDumpFilename, 'wb')) # TEST
			elif doAnalysis:
				try:
					processingResults[eName] = pickle.load(open(processingResultsDumpFilename, 'rb'))
				except FileNotFoundError:
					raise FileNotFoundError("There is no previously processed data in this path: "+processingResultsDumpFilename)
			else:
				logging.warning("No processing step performed nor processing file loaded for experiment "+str(eName)+"!")

		""" Data analysis """
		analysis_path_out = masterParams[eName]['path_out']
		analysisResultsDumpFilename = analysis_path_out + '/analysisResultsDump'
		if doAnalysis:
			# prepare the output directory
			assert not os.path.exists(analysis_path_out)  # write to empty dir
			assert os.path.exists(path.relpath(path.join(analysis_path_out, path.pardir)))  # parent dir must exist
			os.makedirs(analysis_path_out)

			# perform analysis
			analysisResults = analyzeProcessingResult(processingResults, masterParams, writeToDisk)
			pickle.dump(analysisResults, open(analysisResultsDumpFilename, 'wb'))  # TEST
		elif doReport:
			try:
				analysisResults = pickle.load(open(analysisResultsDumpFilename, 'rb'))
			except FileNotFoundError:
				raise FileNotFoundError("There is no previously analyzed data in this path: "+analysisResultsDumpFilename)
		else:
			logging.warning("No analysis step performed nor analysis file loaded!")

		""" Visualize and generate report """
		if doReport:
			# prepare the output directory
			results_path_out = masterParams[eName]['path_results']
			assert not os.path.exists(results_path_out)  # write to empty dir
			assert os.path.exists(path.relpath(path.join(results_path_out, path.pardir)))  # parent dir must exist
			os.makedirs(results_path_out)

			# visualize and make a report
			generateReport(analysisResults, masterParams, logFilePath, writeToDisk)
		else:
			logging.warning("No report generated!")

	elif testing:
		devStuff(dfs[0], specificParams[0])
	stop = time()
	print(stop - start)


if __name__ == '__main__':
	masterConfigFilePath = 'masterConfig.ini' # TEST
	masterConfigFilePath = webFlow()
	sys.exit(main(masterConfigFilePath=masterConfigFilePath, doProcessing=True, doAnalysis=True, doReport=True, testing=False, writeToDisk=False))
