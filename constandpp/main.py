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

import sys, os, logging, datetime
from webFlow import webFlow
from getInput import getProcessingInput, getJobInput
from constand import constand
from time import time
from dataIO import *
from processing import *
from collapse import collapse
from analysis import *
from report import *


def performanceTest():  # remove for production # TEST
	""" Use this development method to test the performance of the CONSTANd algorithm. """
	t = []
	for i in range(100):
		params = getProcessingInput()
		df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(2*10**3, 6)), columns=list('ABCDEF'))
		start = time()
		constand(np.asarray(df), 1e-2, 50)
		stop = time()
		t.append((stop - start))
	print("average runtime: " + str(np.mean(t)))


def isotopicImpuritiesTest(): # TEST
	## test if isotopic correction is necessary:
	params = getProcessingInput()
	# get the dataframe
	df = importDataFrame(params['files_in'], delim=params['delim_in'], header=params['header_in'])
	correctedIntensities = getIntensities(df)
	normalizedIntensities, convergenceTrail, R, S = constand(correctedIntensities, params['accuracy'],
	                                                         params['maxIterations'])
	# exportData(normalizedIntensities, 'txt', path_out=params['path_out'],
	#            filename=params['jobname'] + '_normalizedIntensities', delim_out=params['delim_out'])
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
		if os.path.exists('../data/compareIntensitySNProcessingResults'):
			processingResults = pickle.load(open('../data/compareIntensitySNProcessingResults', 'rb'))
		else:
			params=getProcessingInput()
			# setProcessingGlobals(intensityColumns=params['intensityColumns'],
			# 					 removalColumnsToSave=params['removalColumnsToSave'],
			# 					 noMissingValuesColumns=params['noMissingValuesColumns'])
			# setCollapseColumnsToSave(params['collapseColumnsToSave'])  # define the intensityColumns for use in processing.py
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
	# todo find mean of empty slices warning flood source (ABOVE this line)
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

	# get min and max protein-peptide mappings
	minProteinPeptidesDict, maxProteinPeptidesDict, metadata['noMasterProteinAccession'] = getProteinPeptidesDicts(allExperimentsDF)

	# execute mappings to get all peptideintensities per protein, over each whole condition. Index = 'protein'
	# !!! only for 2 conditions up to now!
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
		           filename=params['jobname'] + '_results_minimal', delim_out=params['delim_out'])
		exportData(fullProteinDF, dataType='df', path_out=params['path_out'],
		           filename=params['jobname'] + '_results_full', delim_out=params['delim_out'])
		# save the metadata
		exportData(metadata, dataType='df', path_out=params['path_out'],
		           filename=params['jobname'] + '_metadata',
		           delim_out=params['delim_out'], inOneFile=False)
		# generate a report PDF (without the normalized intensities: behind paywall?

	return minProteinDF, fullProteinDF, PCAResult, HCResult, metadata


def generateReport(analysisResults, params, logFilePath, writeToDisk):
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
	minVolcanoPlot = getVolcanoPlot(minProteinDF, params['alpha'], params['FCThreshold'],
	                                               params['labelVolcanoPlotAreas'])
	if writeToDisk:
		exportData(minVolcanoPlot, dataType='fig', path_out=params['path_results'],
		           filename=params['jobname'] + '_minVolcanoPlot')
	fullVolcanoPlot = getVolcanoPlot(fullProteinDF, params['alpha'], params['FCThreshold'],
	                                                  params['labelVolcanoPlotAreas'])
	if writeToDisk:
		exportData(fullVolcanoPlot, dataType='fig', path_out=params['path_results'],
		           filename=params['jobname'] + '_fullVolcanoPlot')
	PCAPlot = getPCAPlot(PCAResult, params['schema'])
	if writeToDisk:
		exportData(PCAPlot, dataType='fig', path_out=params['path_results'],
		           filename=params['jobname'] + '_PCAPlot')
	HCDendrogram = getHCDendrogram(HCResult, params['schema'])
	if writeToDisk:
		exportData(HCDendrogram, dataType='fig', path_out=params['path_results'],
		           filename=params['jobname'] + '_HCDendrogram')

	# generate HTML and PDF reports # todo
	htmlReport = makeHTML(minSortedDifferentialProteinsDF, fullSortedDifferentialProteinsDF, minVolcanoPlot,
	                      fullVolcanoPlot, PCAPlot, HCDendrogram, metadata, logFilePath)
	pdfReport = HTMLtoPDF(htmlReport)


def main(masterConfigFilePath, doProcessing, doAnalysis, doReport, writeToDisk, testing):
	"""
	For now this is just stuff for debugging and testing. Later:
	Contains and explicits the workflow of the program. Using the booleans doProcessing, doAnalysis and writeToDisk one
	can control	which parts of the workflow to perform.
	"""
	logFilePath = os.path.abspath(os.path.join(masterConfigFilePath, os.path.join(os.pardir, 'log.txt')))
	logging.basicConfig(filename=logFilePath, level=logging.INFO)
	start = time()
	jobParams = getJobInput(masterConfigFilePath) # config filenames + params for the combination of experiments
	processingParams = {} # specific params for each experiment
	dfs = {}
	processingResults = {}
	experimentNames = list(jobParams['schema'].keys())

	if not testing:
		for eName in experimentNames:
			""" Data processing """
			# get all input parameters
			processingParams[eName] = getProcessingInput(jobParams['schema'][eName]['config'])
			# get the dataframes
			dfs[eName] = getData(processingParams[eName]['data'], delim=processingParams[eName]['delim_in'],
			                     header=processingParams[eName]['header_in'],
			                     wrapper=processingParams[eName]['wrapper'])
			processing_path_out = processingParams[eName]['path_out']
			processingResultsDumpFilename = os.path.join(processing_path_out, 'processingResultsDump_'+str(eName))
			if doProcessing:
				print('Starting processing of ' + eName + '...')
				# prepare the output directories
				if not os.path.exists(processing_path_out):  # do not overwrite dir
					assert os.path.exists(
						os.path.abspath(os.path.join(processing_path_out, os.path.pardir)))  # parent dir must exist
					os.makedirs(processing_path_out)
				else:
					raise Exception("Output path "+processing_path_out+" already exists! Aborting.")

				# process every input dataframe
				logging.info("Starting processing of experiment '" + eName + "' of job '" + jobParams['jobname'] + "' at " +
			             str(datetime.datetime.now()).split('.')[0])
				processingResults[eName] = processDf(dfs[eName], processingParams[eName], writeToDisk)
				logging.info("Finished processing of experiment '" + eName + "' of job '" + jobParams['jobname'] + "' at " +
			             str(datetime.datetime.now()).split('.')[0])
				pickle.dump(processingResults[eName], open(processingResultsDumpFilename, 'wb')) # TEST
			elif doAnalysis:
				try:
					processingResults[eName] = pickle.load(open(processingResultsDumpFilename, 'rb'))
				except FileNotFoundError:
					raise FileNotFoundError("There is no previously processed data in this path: "+processingResultsDumpFilename)
			else:
				logging.warning("No processing step performed nor processing file loaded for experiment "+str(eName)+"!")

		""" Data analysis """
		analysis_path_out = jobParams['path_out']
		analysisResultsDumpFilename = os.path.join(analysis_path_out, 'analysisResultsDump')
		if doAnalysis:
			print('Starting analysis...')
			# prepare the output directories
			if not os.path.exists(analysis_path_out):  # do not overwrite dir
				assert os.path.exists(os.path.abspath(os.path.join(analysis_path_out, os.path.pardir)))  # parent dir must exist
				os.makedirs(analysis_path_out)

			# perform analysis
			logging.info("Starting analysis of job: " + jobParams['jobname'] + "at " +
			             str(datetime.datetime.now()).split('.')[0])
			analysisResults = analyzeProcessingResult(processingResults, jobParams, writeToDisk)
			logging.info("Finished analysis of job: " + jobParams['jobname'] + "at " +
			             str(datetime.datetime.now()).split('.')[0])
			pickle.dump(analysisResults, open(analysisResultsDumpFilename, 'wb'))  # TEST
		elif doReport:
			try:
				analysisResults = pickle.load(open(analysisResultsDumpFilename, 'rb'))
			except FileNotFoundError:
				raise FileNotFoundError("There is no previously analyzed data in this path: "+analysisResultsDumpFilename)
		else:
			logging.warning("No analysis step performed nor analysis file loaded!")

		""" Visualize and generate report """
		results_path_out = jobParams['path_results']

		if doReport:
			print('Starting visualization and report generation...')
			# prepare the output directories
			if not os.path.exists(results_path_out):  # do not overwrite dir
				assert os.path.exists(os.path.abspath(os.path.join(results_path_out, os.path.pardir)))  # parent dir must exist
				os.makedirs(results_path_out)

			# visualize and make a report
			logging.info("Starting visualization end report generation of job: " + jobParams['jobname'] + "at " +
			             str(datetime.datetime.now()).split('.')[0])
			generateReport(analysisResults, jobParams, logFilePath, writeToDisk)
			logging.info("Finished visualization end report generation of job: " + jobParams['jobname'] + "at " +
			             str(datetime.datetime.now()).split('.')[0])
		else:
			logging.warning("No report generated!")

	elif testing:
		devStuff(dfs[0], processingParams[0])
	stop = time()
	print(stop - start)


if __name__ == '__main__':
	#masterConfigFilePath = 'jobConfig.ini' # TEST
	#masterConfigFilePath = webFlow(exptype='COON', previousjobdirName='2016-12-11 22:17:20.393578_COON')
	masterConfigFilePath = webFlow(exptype='COON')
	#masterConfigFilePath = webFlow(exptype='COON_SN', previousjobdirName='2016-12-11 23:56:33.683926_COON_SN')
	#masterConfigFilePath = webFlow(exptype='COON_SN')
	#masterConfigFilePath = webFlow(exptype='COON_norm') # todo constand uitzetten
	#masterConfigFilePath = webFlow(exptype='COON_norm', previousjobdirName='2016-12-12 10:10:37.693588_COON_norm')  # todo constand uitzetten
	#masterConfigFilePath = webFlow(exptype='COON_SN_norm')  # todo constand uitzetten
	#masterConfigFilePath = webFlow(exptype='COON_SN_norm', previousjobdirName='2016-12-12 10:22:41.783491_COON_SN_norm')  # todo constand uitzetten

	sys.exit(main(masterConfigFilePath=masterConfigFilePath, doProcessing=True, doAnalysis=True, doReport=True, testing=False, writeToDisk=True))
