#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Handle all I/O of data files and parameters to and from both the workflow and the main dataFrame.
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
	# TODO add .lower() to all string input
	# path_in='../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.tsv' # TEST
	path_in = '../data/MB_noapostrophes.tsv'  # TEST
	delim_in = '\t'
	header_in = 0
	collapsePSMAlgo_bool = True
	collapsePSMAlgo_master = 'mascot'
	collapseRT_bool = True
	collapseRT_centermeasure = 'mean'
	collapseCharge_bool = True
	isotopicCorrectionsMatrix = np.asmatrix(np.diag(np.ones([6,])))
	accuracy = 1e-2
	maxIterations = 50
	path_out = '../data/MB_result.tsv'  # TEST
	delim_out = '\t'

	if not path.exists(path_in):
		raise FileNotFoundError("File "+path_in+" not found.")
	if not (len(delim_in) == 1 and isinstance(delim_in, str)):
		raise Exception("Delimiter of input file must be a character (string of length one).")
	if not ((isinstance(header_in, int) and header_in >= 0) or header_in is None):
		raise Exception("Header parameter of the input file must be a non-negative integer or of type None.")
	if collapsePSMAlgo_bool is None:
		raise Exception("Please indicate whether you would like to remove redundancy due to multuiple PSM Algorithms.")
	if collapsePSMAlgo_master not in ('mascot', 'sequest'):
		raise Exception("Invalid master PSM algorithm: '"+collapsePSMAlgo_master+"'. Please pick 'mascot' or 'sequest'.")
	if collapseRT_bool is None:
		raise Exception("Please indicate whether you would like to remove redundancy due to multiple retention times.")
	if collapseRT_centermeasure not in ('mean', 'median'):
		raise Exception("Invalid center measure: '"+collapseRT_centermeasure+"'. Please pick 'mean' or 'median'.")
	if collapseCharge_bool is None:
		raise Exception("Please indicate whether you would like to remove redundancy due to multiple charge states.")
	if not (isotopicCorrectionsMatrix.shape == [6,6]):
		raise Exception("Isotopic corrections matrix must have shape (6,6).")
	if not (accuracy > 0):
		raise Exception("Accuracy must be strictly greater than zero.")
	if not (maxIterations > 0 and isinstance(maxIterations,int)):
		raise Exception("Maximum number of iterations must be an integer strictly greater than zero.")
	if path.exists(path_out):
		# raise Exception("The file "+path_out+" already exists.")
		warnings.warn("Will overwrite file "+path.basename(path.normpath(path_out)))
	if not (len(delim_out) == 1 and isinstance(delim_out, str)):
		raise Exception("Delimiter of output file must be a character (string of length one).")

	params = {
		'path_in': path_in,
		'delim_in': delim_in,
		'header_in': header_in,
		'collapsePSMAlgo_bool': collapsePSMAlgo_bool,
		'collapsePSMAlgo_master': collapsePSMAlgo_master,
		'collapseRT_bool': collapseRT_bool,
		'collapseRT_centermeasure': collapseRT_centermeasure,
		'collapseCharge_bool': collapseCharge_bool,
		'isotopicCorrectionsMatrix': isotopicCorrectionsMatrix,
		'accuracy': accuracy,
		'maxIterations': maxIterations,
		'path_out': path_out,
		'delim_out': delim_out
	}
	return params


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

	# df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10, 6)), columns=list('ABCDEF'))  # TEST
	# df = pd.read_csv('../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.txt', delim='\t') # TEST
	# df = pd.DataFrame(np.arange(10*6).reshape(10,6),columns=list('ABCDEF')) # TEST
	# df['B'][0]=np.nan # TEST
	# df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10**3, 6)), columns=list('ABCDEF'))  # TEST

	return df


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
