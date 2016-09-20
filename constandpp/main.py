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
		path, delim, accuracy, maxIterations = getInput()
		intensities, df = importData(path, delim)
		start = time()
		constand(intensities, 1e-2, 50)
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
	# get all input parameters
	path_in, delim_in, header_in, accuracy, maxIterations, path_out, delim_out = getInput()
	# get the data
	intensities, df = importData(path_in, delim_in, header_in)
	# collapse peptide list redundancy due to overlap in MASCOT/SEQUEST peptide matches
	df = collapsePSMAlgo(df) # TODO
	# collapse peptide list redundancy due to multiple detections at different RT
	df = collapseRT(df) # TODO
	# collapse peptide list redundancy due to different charges (optional)
	df = collapseCharge(df) # TODO
	# perform isotopic corrections
	intensities = isotopicCorrection(intensities) # TODO
	# perform the CONSTANd algorithm
	normalizedIntensities, convergenceTrail, R, S = constand(intensities, accuracy, maxIterations)
	# save the normalized intensities obtained through CONSTANd
	exportData(normalizedIntensities, path_out, delim_out)
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
