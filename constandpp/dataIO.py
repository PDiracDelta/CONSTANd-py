#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Handle all I/O of data files to and from the workflow.
"""

import warnings
import pandas as pd
import numpy as np
from os import path


def getInput():
	"""
	Get mass spec data and CONSTANd parameters from the user or from the web interface as a dict.
	:return params:         dict    dictionary containing all paraeters mentioned below:
		:return path_in:        string  path to the input file
		:return delim_in:       char    delimiter of the data in the input file
		:return accuracy:       float   CONSTANd param: combined allowed deviation of col and row means from 1/6
		:return header_in:      integer row number containing the dataFrame header (can be None if no header)
		:return maxIterations:  int     CONSTANd param: maximum amount of iterations (1x row and 1x col per iteration)
		:return path_out:       string  path to the output file
		:return delim_out:      char    delimiter of the data in the output file
	"""
	# path_in='../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.tsv' # TEST
	path_in = '../data/MB_noapostrophes.tsv'  # TEST
	delim_in = '\t'
	header_in = 0
	accuracy = 1e-2
	maxIterations = 50
	path_out = '../data/MB_result.tsv'  # TEST
	delim_out = '\t'

	if not path.exists(path_in):
		raise FileNotFoundError("File "+path_in+" not found.")
	if path.exists(path_out):
		# raise Exception("The file "+path_out+" already exists.")
		warnings.warn("Overwriting file "+path.basename(path.normpath(path_out)))
	if not (len(delim_in) == 1 and isinstance(delim_in, str)):
		raise Exception("Delimiter of input file must be a character (string of length one).")
	if not(len(delim_out) == 1 and isinstance(delim_out, str)):
		raise Exception("Delimiter of output file must be a character (string of length one).")
	if not (accuracy > 0):
		raise Exception("Accuracy must be strictly greater than zero.")
	if not (maxIterations > 0 and isinstance(maxIterations,int)):
		raise Exception("Maximum number of iterations must be an integer strictly greater than zero.")

	params = {'path_in': path_in,
	          'delim_in': delim_in,
	          'header_in': header_in,
	          'accuracy': accuracy,
	          'maxIterations': maxIterations,
	          'path_out': path_out,
	          'delim_out': delim_out}
	return params


def importData(path_in=None, delim=None, header_in=0):
	"""
	Return the intensity matrix and the dataFrame of the data specified.
	:param path_in:         string          path to the data file
	:param delim:           char            delimiter of the data
	:return intensities:    np.ndArray      (N,6) ndarray with the absolute intensities
	:return df:             pd.dataFrame    Pandas dataframe with the contents of the data file, including the intensities
	"""
	# df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10, 6)), columns=list('ABCDEF'))  # TEST
	# df = pd.read_csv('../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.txt', delim='\t') # TEST
	# df = pd.DataFrame(np.arange(10*6).reshape(10,6),columns=list('ABCDEF')) # TEST
	# df['B'][0]=np.nan # TEST
	# df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10**3, 6)), columns=list('ABCDEF'))  # TEST
	df = importDataFrame(path_in, delim=delim, header=header_in)
	intensities = selectIntensities(df)
	assert isinstance(intensities, np.ndarray) and intensities.shape[1] == 6 and intensities.dtype == 'float64'

	return intensities, df


def exportData(data=None, path_in=None, delim=','):
	"""
	Save the results (normalized intensities) to disk.
	:param data:        obj     data object to be exported to disk
	:param path_in:     string  path+filename where data should be exported to
	:param delim:       char    delimiter of the data
	"""
	assert data is not None
	assert path.exists(path_in)

	np.savetxt(path_in, data, delimiter=delim)


def importDataFrame(path_in=None, filetype=None, delim=None, header=0):
	"""
	Get the data from disk as a Pandas DataFrame.
	:param path_in:     string          existing path to input file
	:param filetype:    string          specifier for the type of the file (file extension)
	:param delim:       char            delimiter of the data
	:return df:         pd.dataFrame    Pandas dataFrame of the file contents
	"""
	assert path.exists(path_in)

	if filetype is None:  # set filetype equal to the file extension
		filetype = path_in.split('.')[-1]

	if filetype == 'xlsx':
		df = pd.read_excel(path_in)
	elif filetype == 'csv':
		if delim is None:
			df = pd.read_csv(path_in, delimiter=',', header=header)
		else:
			df = pd.read_csv(path_in, delimiter=delim, header=header)
	elif filetype == 'tsv':
		if delim is None:
			df = pd.read_csv(path_in, delimiter='\t', header=header)
		else:
			df = pd.read_csv(path_in, delimiter=delim, header=header)
	else:
		if delim is None:
			raise Exception(
				"I don't know how to handle this data: the filetype was not recognized and no delimiter was specified.")
		warnings.warn("Did not recognize filetype: treating as delimited textfile with the delimiter you specified.")
		df = pd.read_csv(path_in, delimiter=delim, header=header)

	return df


def selectIntensities(df):
	"""
	Extracts the (absolute) intensity matrix from the dataFrame.
	:param df:              pd.dataFrame    Pandas dataFrame from which to extract the intensities
	:return intensities:    np.ndArray      matrix with the intensities
	"""
	intensities = np.asarray(df[['126', '127', '128', '129', '130', '131']])

	return intensities
