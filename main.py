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
	# get all input parameters and option switches
	""" params:
	file_in, delim_in, header_in, collapsePSMAlgo_bool, collapsePSMAlgo_master, collapsePSMAlgo_bool_exclusive, collapseRT_bool,
	collapseRT_centerMeasure, collapseCharge_bool, isotopicCorrectionsMatrix, accuracy, maxIterations, path_out, filename_out, delim_out
	 """
	params = getInput()
	# get the dataframe
	df = importDataFrame(params('file_in'), delim=params('delim_in'), header=params('header_in'))
	# add extra columns to the dataFrame for retaining condensed data after each collapse, according to bools (or not).
	addColumns(df, bools=None)
	if params['removeIsolationInterference_bool']:
		df = removeIsolationInterference(df, params['removeIsolationInterference_threshold'])
	if params['collapsePSMAlgo_bool']:
		# collapse peptide list redundancy due to overlap in MASCOT/SEQUEST peptide matches
		df = collapsePSMAlgo(df, master=params['collapsePSMAlgo_master'], exclusive=params['collapsePSMAlgo_bool_exclusive']) # TODO
	if params['collapseRT_bool']:
		# collapse peptide list redundancy due to multiple detections at different RT
		df = collapseRT(df, centerMeasure=params['collapseRT_centerMeasure']) # TODO
	if params['collapseCharge_bool']:
		# collapse peptide list redundancy due to different charges (optional)
		df = collapseCharge(df) # TODO
	# perform isotopic corrections
	intensities = isotopicCorrection(getIntensities(df), correctionsMatrix=params['isotopicCorrectionsMatrix']) # TODO
	# perform the CONSTANd algorithm
	normalizedIntensities, convergenceTrail, R, S = constand(intensities, params('accuracy'), params('maxIterations'))
	# save the normalized intensities obtained through CONSTANd
	exportData(normalizedIntensities, path=params('path_out'), filename=params['filename_out'], delim=params('delim_out'))
	# perform differential expression analysis
	DEresults = differentialExpression(normalizedIntensities) # TODO
	# save the DE analysis results
	exportData(DEresults) # TODO
	# data visualization
	viz = dataVisualization(DEresults) # TODO
	# save the visualizations
	exportData(viz) # TODO
	# generate a report PDF (without the normalized intensities: behind paywall?
	generateReport(DEresults, viz) # TODO

	# performanceTest()

if __name__ == '__main__':
	sys.exit(main())
