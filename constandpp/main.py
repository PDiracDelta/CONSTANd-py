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
	for i in range(1000):
		params = getInput()
		df = importDataFrame(params['path'], delim=params['delim_in'])
		start = time()
		constand(getIntensities(df), 1e-2, 50)
		stop = time()
		t.append((stop - start))
	print("average runtime: " + str(np.mean(t)))


def generateReport(DEresults, viz):
	pass


def main():
	"""
	For now this is just stuff for debugging and testing. Later:
	Contains and explicits the workflow of the program.
	"""
	""" get all input parameters
	params:
	file_in, delim_in, header_in, collapsePSMAlgo_bool, collapsePSMAlgo_master, collapsePSMAlgo_bool_exclusive,
	collapseRT_bool, collapseRT_centerMeasure_channels, collapseRT_centerMeasure_intensities,
	collapseRT_maxRelativeChannelVariance, collapseCharge_bool, isotopicCorrectionsMatrix, accuracy, maxIterations,
	DEFoldThreshold, path_out, filename_out, delim_out
	"""
	params = getInput()
	# get the dataframe
	df = importDataFrame(params['file_in'], delim=params['delim_in'], header=params['header_in'])

	""" Data preparation """
	removedData={} # is to contain basic info about data that will be removed during the workflow, per removal category.
	if params['removeIsolationInterference_bool']:
		df, removedData['isolationInterference'] = removeIsolationInterference(df, params['removeIsolationInterference_threshold'])
		print(str(df.shape)+', '+str(removedData['isolationInterference'].shape)) # TEST
	if params['collapsePSMAlgo_bool']:
		# collapse peptide list redundancy due to overlap in MASCOT/SEQUEST peptide matches
		df = collapsePSMAlgo(df, master=params['collapsePSMAlgo_master'],
		                     exclusive=params['collapsePSMAlgo_bool_exclusive']) # TODO
	if params['collapseRT_bool']:
		# collapse peptide list redundancy due to multiple detections at different RT
		df = collapseRT(df, centerMeasure_channels=params['collapseRT_centerMeasure_channels'],
		                centerMeasure_intensities=params['collapseRT_centerMeasure_intensities'],
		                maxRelativeChannelVariance=params['collapseRT_maxRelativeChannelVariance']) # TODO
	if params['collapseCharge_bool']:
		# collapse peptide list redundancy due to different charges (optional)
		df = collapseCharge(df) # TODO
	# perform isotopic corrections but do NOT apply them to df because this information is sensitive (copyright i-TRAQ)
	correctedIntensities = isotopicCorrection(getIntensities(df), correctionsMatrix=params['isotopicCorrectionsMatrix']) # TODO
	# perform the CONSTANd algorithm; also do NOT include normalized intensities in df --> only for paying users.
	normalizedIntensities, convergenceTrail, R, S = constand(correctedIntensities, params['accuracy'], params['maxIterations'])

	""" Data analysis and visualization """
	# perform differential expression analysis
	DEresults = differentialExpression(normalizedIntensities, params['DEFoldThreshold']) # TODO
	# data visualization
	viz = dataVisualization(DEresults) # TODO

	""" Save data to disk and generate report """
	# save the removed data information
	exportData(removedData, path_out=params['path_out'],
	           filename=params['filename_out'] + '_removedData', delim_out=params['delim_out'])
	# save the normalized intensities obtained through CONSTANd
	exportData(normalizedIntensities, path_out=params['path_out'],
	           filename=params['filename_out'] + '_normalizedIntensities', delim_out=params['delim_out'])
	# save the DE analysis results
	exportData(DEresults, path_out=params['path_out'], filename=params['filename_out'] + '_DEresults')  # TODO
	# save the visualizations
	exportData(viz, path_out=params['path_out'], filename=params['filename_out']+'_dataViz') # TODO
	# generate a report PDF (without the normalized intensities: behind paywall?
	generateReport(DEresults, viz) # TODO

	# performanceTest()

if __name__ == '__main__':
	sys.exit(main())
