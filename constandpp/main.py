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
from time import time
from dataIO import *
from dataprep import *
from collapse import collapse, setCollapseColumnsToSave
from analysis import *


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
	df = importDataFrame(params['file_in'], delim=params['delim_in'], header=params['header_in'])
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
		icm = params['isotopicCorrectionsMatrix']
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


def generateReport(DEresults, viz):
	pass


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


def devStuff(df, params): # TEST
	# performanceTest()
	# isotopicCorrectionsTest(params)
	# MS2IntensityDoesntMatter(df)
	# testDataComplementarity(df)
	pass


def main(testing, writeToDisk):
	start = time()
	# testing=True # TEST
	# writeToDisk=False # TEST
	"""
	For now this is just stuff for debugging and testing. Later:
	Contains and explicits the workflow of the program.
	"""
	""" get all input parameters
	params:
	file_in, delim_in, header_in, intensityColumns, requiredColumns, collapseColumnsToSave, undoublePSMAlgo_bool,
	removeIsolationInterference_bool, collapse_method, collapse_maxRelativeReporterVariance,
	removeIsolationInterference_master, masterPSMAlgo, undoublePSMAlgo_exclusive_bool, collapseRT_bool,
	collapseCharge_bool, collapsePTM_bool, isotopicCorrection_bool, isotopicCorrectionsMatrix, accuracy, maxIterations,
	DEFoldThreshold, path_out, filename_out, delim_out
	"""
	params = getInput()
	# get the dataframe
	df = importDataFrame(params['file_in'], delim=params['delim_in'], header=params['header_in'])
	if not testing:
		""" Data preparation """
		removedData={} # is to contain basic info about data that will be removed during the workflow, per removal category.
		# define global parameters
		setGlobals(intensityColumns = params['intensityColumns'], remove_ExtraColumnsToSave=params['remove_ExtraColumnsToSave'],
		           noMissingValuesColumns=params['noMissingValuesColumns'])
		setCollapseColumnsToSave(params['collapseColumnsToSave'])  # define the intensityColumns for use in dataprep.py
		# remove detections where (essential) data is missing.
		df, removedData['missing'] = removeMissing(df)
		if params['removeBadConfidence_bool']:
			df, removedData['confidence'] = removeBadConfidence(df, params['removeBadConfidence_minimum'])
		# remove all useless columns from the dataFrame
		df = selectRequiredColumns(df, requiredColumns=params['requiredColumns'])
		if params['removeIsolationInterference_bool']:
			# remove all data with too high isolation interference
			df, removedData['isolationInterference'] = removeIsolationInterference(df, params['removeIsolationInterference_threshold'])
		if params['undoublePSMAlgo_bool']:
			# collapse peptide list redundancy due to overlap in MASCOT/SEQUEST peptide matches
			df, removedData['PSMAlgo'] = undoublePSMAlgo(df, master=params['masterPSMAlgo'],
			                                             exclusive=params['undoublePSMAlgo_exclusive_bool'])
			# SANITY CHECK: no detections with the same scan number may exist after undoublePSMAlgo()
			assert np.prod((len(i) < 2 for (s, i) in df.groupby('First Scan').groups))

		# collapse peptide list redundancy due to multiple detections at different RT
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
			                                   maxRelativeReporterVariance=params['collapse_maxRelativeReporterVariance'],
			                                            masterPSMAlgo=params['masterPSMAlgo'],
			                                            undoublePSMAlgo_bool=params['undoublePSMAlgo_bool'])

		# SANITY CHECK: there should be no more duplicates if all collapses have been applied.
		if params['undoublePSMAlgo_bool'] and params['collapseRT_bool'] and params['collapseCharge_bool']: # TEST
			assert np.prod((len(i) < 2 for (s, i) in df.groupby('Annotated Sequence').groups)) # only 1 index vector in dict of SEQUENCE:[INDICES] for all sequences

		if params['isotopicCorrection_bool']:
			# perform isotopic corrections but do NOT apply them to df because this information is sensitive (copyright i-TRAQ)
			intensities = isotopicCorrection(getIntensities(df), correctionsMatrix=params['isotopicCorrectionsMatrix'])
		else:
			intensities = getIntensities(df)
		# perform the CONSTANd algorithm; also do NOT include normalized intensities in df --> only for paying users.
		normalizedIntensities, convergenceTrail, R, S = constand(intensities, params['accuracy'], params['maxIterations'])

		""" Data analysis and visualization """
		# perform differential expression analysis

		# the ID of each result is determined by the collection of propertes that were NOT collapsed (so if all are collapsed its only the sequence)
		DEresults = differentialExpression(normalizedIntensities, params['DEFoldThreshold']) # TODO
		# data visualization
		viz = dataVisualization(DEresults) # TODO

		""" Save data to disk and generate report """
		if writeToDisk:
			# save the removed data information
			exportData(removedData, dataType='df', path_out=params['path_out'], filename=params['filename_out'] + '_removedData',
			           delim_out=params['delim_out'], inOneFile=params['unusedDetectionsInOneFile_bool'])
			# save the final form of the dataFrame WITHOUT normalized intensities.
			exportData(df, dataType='df', path_out=params['path_out'], filename=params['filename_out'] + '_dataFrame', delim_out=params['delim_out'])
			# save the normalized intensities obtained through CONSTANd
			exportData(normalizedIntensities, dataType='txt', path_out=params['path_out'],
			           filename=params['filename_out'] + '_normalizedIntensities', delim_out=params['delim_out'])
			# save the DE analysis results
			exportData(DEresults, dataType='obj', path_out=params['path_out'], filename=params['filename_out'] + '_DEresults')  # TODO
			# save the visualizations
			exportData(viz, dataType='viz', path_out=params['path_out'], filename=params['filename_out']+'_dataViz') # TODO
			# generate a report PDF (without the normalized intensities: behind paywall?
			generateReport(DEresults, viz) # TODO

	elif testing:
		devStuff(df, params)
	stop = time()
	print(stop - start)


if __name__ == '__main__':
	sys.exit(main(testing=False, writeToDisk=True))
