#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Python implementation of mass spectrometer protein data analysis using the CONSTANd_RAS algorithm
"""

__author__ = "Joris Van Houtven"
__copyright__ = "Copyright ?, VITO"
__credits__ = ["Joris Van Houtven", "Dirk Valkenborg"]
# __license__ = "GPL"
# __version__ = "0.0.1"
__maintainer__ = "Joris Van Houtven"
__email__ = "vanhoutvenjoris@gmail.com"
__status__ = "Development"

import sys
import pandas as pd
import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
from constand import constand
from time import time


def getInput():
	""" Get mass spec data and CONSTANd parameters """
	path=None
	delim='\t'
	accuracy=1e-2
	maxIterations=50
	return path,delim,accuracy,maxIterations


def importData(path, delim):
	# df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10, 6)), columns=list('ABCDEF'))  # TEST
	# df = pd.read_csv('../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.txt', delim='\t') # TEST
	# df = pd.DataFrame(np.arange(10*6).reshape(10,6),columns=list('ABCDEF')) # TEST
	# df['B'][0]=np.nan # TEST
	df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10**3, 6)), columns=list('ABCDEF'))  # TEST
	intensities = selectIntensities(df)
	assert isinstance(intensities, np.ndarray) # ndarray instead of matrix because this is more convenient in the calculations
	return intensities, df


def importDataFrame(path=None, filetype=None, delim=None):
	""" Get the data from disk as a Pandas DataFrame """
	assert path is not None
	if filetype is None: # set filetype equal to the file extension
		filetype = path.split('.')(-1)

	if filetype == 'xlsx':
		df = pd.read_excel(path)
	elif filetype == 'csv':
		df = pd.read_csv(path)
	elif filetype == 'tsv':
		df = pd.read_table(path)
	else:
		df = pd.read_table(path, delim=delim)

	data = np.asarray(df)
	return data


def selectIntensities(df):
	intensities=None
	return intensities


def performanceTest(): # remove for production
	""" Use this development method to test the performance of the CONSTANd algorithm """
	t = []
	for i in range(1000):
		path, delim, accuracy, maxIterations = getInput()
		data = importData(path, delim)
		start = time()
		constand(data, 1e-2, 50)
		stop = time()
		t.append((stop - start))
	print("average runtime: " + str(np.mean(t)))


def main():
	""" For now this is just stuff for debugging and testing """
	path,delim,accuracy,maxIterations = getInput()
	intensities, df = importData(path,delim)
	assert isinstance(intensities, np.ndarray)
	normalizedIntensities,convergenceTrail,R,S = constand(intensities,accuracy,maxIterations)
	#print(normalizedIntensities)
	#performanceTest()


if __name__ == '__main__':
	sys.exit(main())
