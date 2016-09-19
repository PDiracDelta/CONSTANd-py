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


def main():
	""" For now this is just stuff for debugging and testing. """
	path_in, delim_in, accuracy, maxIterations, path_out, delim_out = getInput()
	intensities, df = importData(path_in, delim_in)
	assert isinstance(intensities, np.ndarray)
	normalizedIntensities, convergenceTrail, R, S = constand(intensities, accuracy, maxIterations)
	# print(normalizedIntensities)
	exportData(normalizedIntensities, path_out, delim_out)
	# performanceTest()


if __name__ == '__main__':
	sys.exit(main())
