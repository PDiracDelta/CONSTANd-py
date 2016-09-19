#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Handle all I/O of data files to and from the workflow.
"""

import warnings
import pandas as pd
import numpy as np


def getInput():
	""" Get mass spec data and CONSTANd parameters from the user or from the web interface. """
	# path='../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.tsv' # TEST
	path_in = '../data/MB_noapostrophes.tsv'  # TEST
	delim_in = '\t'
	accuracy = 1e-2
	maxIterations = 50
	path_out = '../data/MB_result.tsv'  # TEST
	delim_out = '\t'
	return path_in, delim_in, accuracy, maxIterations, path_out, delim_out


def importData(path=None, delim=None):
	""" Return the intensity matrix and the dataFrame of the data specified. """
	# df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10, 6)), columns=list('ABCDEF'))  # TEST
	# df = pd.read_csv('../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.txt', delim='\t') # TEST
	# df = pd.DataFrame(np.arange(10*6).reshape(10,6),columns=list('ABCDEF')) # TEST
	# df['B'][0]=np.nan # TEST
	# df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10**3, 6)), columns=list('ABCDEF'))  # TEST
	df = importDataFrame(path, delim=delim)
	intensities = selectIntensities(df)
	assert isinstance(intensities,
	                  np.ndarray)  # ndarray instead of matrix because this is more convenient in the calculations
	return intensities, df


def exportData(data=None, path=None, delim=','):
	""" Save the results (normalized intensities) to disk. """
	assert data is not None
	assert path is not None
	assert delim is not None

	np.savetxt(path, data, delimiter=delim)


def importDataFrame(path=None, filetype=None, delim=None):
	""" Get the data from disk as a Pandas DataFrame. """
	assert path is not None
	if filetype is None:  # set filetype equal to the file extension
		filetype = path.split('.')[-1]

	if filetype == 'xlsx':
		df = pd.read_excel(path)
	elif filetype == 'csv':
		if delim is None:
			df = pd.read_csv(path, delimiter=',')
		else:
			df = pd.read_csv(path, delimiter=delim)
	elif filetype == 'tsv':
		if delim is None:
			df = pd.read_csv(path, delimiter='\t')
		else:
			df = pd.read_csv(path, delimiter=delim)
	else:
		if delim is None:
			raise Exception(
				"I don't know how to handle this data: the filetype was not recognized and no delimiter was specified.")
		warnings.warn("Did not recognize filetype: treating as delimited textfile with the delimiter you specified.")
		df = pd.read_csv(path, delimiter=delim)

	return df


def selectIntensities(df):
	""" Extracts the intensity matrix from the dataFrame. """
	intensities = np.asarray(df[['126', '127', '128', '129', '130', '131']])
	return intensities