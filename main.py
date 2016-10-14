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
from analysis import *


def performanceTest():  # remove for production
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


def isotopicImpuritiesTest():
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

def generateReport(DEresults, viz):
	pass


def devStuff():
	dummy = getInput()
	#pass
	performanceTest()
	# isotopicImpuritiesTest()


def main():
	testing=False # TEST
	writeToDisk=False # TEST
	"""
	For now this is just stuff for debugging and testing. Later:
	Contains and explicits the workflow of the program.
	"""
	""" get all input parameters
	params:
	file_in, delim_in, header_in, collapsePSMAlgo_bool, removeIsolationInterference_bool,
	removeIsolationInterference_master, collapsePSMAlgo_master, collapsePSMAlgo_exclusive_bool,	collapseRT_bool,
	collapseRT_centerMeasure_channels, collapseRT_centerMeasure_intensities, collapseRT_maxRelativeReporterVariance,
	collapseCharge_bool, isotopicCorrection_bool, isotopicCorrectionsMatrix, accuracy, maxIterations, DEFoldThreshold,
	path_out, filename_out, delim_out
	"""
	params = getInput()
	# get the dataframe
	df = importDataFrame(params['file_in'], delim=params['delim_in'], header=params['header_in'])
	if not testing:
		""" Data preparation """
		removedData={} # is to contain basic info about data that will be removed during the workflow, per removal category.
		if params['removeIsolationInterference_bool']:
			df, removedData['isolationInterference'] = removeIsolationInterference(df, params['removeIsolationInterference_threshold'])
		if params['collapsePSMAlgo_bool']:
			# collapse peptide list redundancy due to overlap in MASCOT/SEQUEST peptide matches
			df, removedData['PSMAlgo'] = collapsePSMAlgo(df, master=params['collapsePSMAlgo_master'],
			                                             exclusive=params['collapsePSMAlgo_exclusive_bool'])
		if params['collapseRT_bool']:
			# collapse peptide list redundancy due to multiple detections at different RT
			df = collapseRT(df, centerMeasure_channels=params['collapseRT_centerMeasure_channels'],
			                centerMeasure_intensities=params['collapseRT_centerMeasure_intensities'],
			                maxRelativeReporterVariance=params['collapseRT_maxRelativeReporterVariance']) # TODO
		if params['collapseCharge_bool']:
			# collapse peptide list redundancy due to different charges (optional)
			df = collapseCharge(df) # TODO
		if params['isotopicCorrection_bool']:
			# perform isotopic corrections but do NOT apply them to df because this information is sensitive (copyright i-TRAQ)
			intensities = isotopicCorrection(getIntensities(df), correctionsMatrix=params['isotopicCorrectionsMatrix'])
		else:
			intensities = getIntensities(df)
		# perform the CONSTANd algorithm; also do NOT include normalized intensities in df --> only for paying users.
		normalizedIntensities, convergenceTrail, R, S = constand(intensities, params['accuracy'], params['maxIterations'])

		""" Data analysis and visualization """
		# perform differential expression analysis
		DEresults = differentialExpression(normalizedIntensities, params['DEFoldThreshold']) # TODO
		# data visualization
		viz = dataVisualization(DEresults) # TODO

		""" Save data to disk and generate report """
		if writeToDisk:
			# save the removed data information
			exportData(removedData, path_out=params['path_out'],
			           filename=params['filename_out'] + '_removedData', delim_out=params['delim_out'])
			# save the normalized intensities obtained through CONSTANd
			exportData(normalizedIntensities, dataType='txt', path_out=params['path_out'],
			           filename=params['filename_out'] + '_normalizedIntensities', delim_out=params['delim_out'])
			# save the DE analysis results
			exportData(DEresults, path_out=params['path_out'], filename=params['filename_out'] + '_DEresults')  # TODO
			# save the visualizations
			exportData(viz, path_out=params['path_out'], filename=params['filename_out']+'_dataViz') # TODO
			# generate a report PDF (without the normalized intensities: behind paywall?
			generateReport(DEresults, viz) # TODO

	if testing:
		devStuff()

if __name__ == '__main__':
	sys.exit(main())
