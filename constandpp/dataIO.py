#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Handle all I/O of data files and parameters to and from both the workflow and the main dataFrame.
"""

import warnings
import pandas as pd
import numpy as np
import configparser
from os import path
from codecs import getdecoder as gd


def getInput():
	"""
	Get mass spec data and CONSTANd parameters from the user or from the web interface as a dict.
	:return params:         dict    dictionary containing all paraeters mentioned below:
		:return file_in:        string  path to the input file
		:return delim_in:       char    delimiter of the data in the input file
		:return accuracy:       float   CONSTANd param: combined allowed deviation of col and row means from 1/6
		:return header_in:      integer row number containing the dataFrame header (can be None if no header)
		:return maxIterations:  int     CONSTANd param: maximum amount of iterations (1x row and 1x col per iteration)
		:return path_out:       string  path to the output file
		:return delim_out:      char    delimiter of the data in the output file
	"""
	# TODO add all parameters in docstring
	# TODO add .lower() to all string input
	# TODO attach real input source
	# file_in='../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.tsv' # TEST

	# read the config file to obtain the defaults
	config = configparser.ConfigParser(allow_no_value=True, inline_comment_prefixes='#') # TODO split up CONFIG and PARAMS (user input files vs workflow params)
	config.optionxform = str # so that strings dont automatically get .lower()-ed
	config.read('config.ini', encoding='utf-8')
	dontForgetAnyParameters = len(config._defaults)

	# get variables from config in correct typography
	file_in = config.get('DEFAULT','file_in')
	delim_in = gd("unicode_escape")(config.get('DEFAULT','delim_in'))[0] # treat delimiters correctly: ignore first escape
	header_in = config.getint('DEFAULT','header_in')
	removeIsolationInterference_bool = config.getboolean('DEFAULT','removeIsolationInterference_bool')
	removeIsolationInterference_threshold = config.getfloat('DEFAULT','removeIsolationInterference_threshold')
	collapsePSMAlgo_bool = config.getboolean('DEFAULT','collapsePSMAlgo_bool')
	collapsePSMAlgo_master = config.get('DEFAULT','collapsePSMAlgo_master')
	collapsePSMAlgo_exclusive_bool = config.getboolean('DEFAULT','collapsePSMAlgo_exclusive_bool')
	collapseRT_bool = config.getboolean('DEFAULT','collapseRT_bool')
	collapseRT_centerMeasure_channels = config.get('DEFAULT','collapseRT_centerMeasure_channels')
	collapseRT_centerMeasure_intensities = config.get('DEFAULT','collapseRT_centerMeasure_intensities')
	collapseRT_maxRelativeReporterVariance = config.getfloat('DEFAULT','collapseRT_maxRelativeReporterVariance')
	collapseCharge_bool = config.getboolean('DEFAULT','collapseCharge_bool')
	isotopicCorrectionsMatrix = eval(config.get('DEFAULT','isotopicCorrectionsMatrix')) # TODO: separate file?
	accuracy = config.getfloat('DEFAULT','accuracy')
	maxIterations = config.getint('DEFAULT','maxIterations')
	DEFoldThreshold = config.getfloat('DEFAULT','DEFoldThreshold')
	path_out = config.get('DEFAULT','path_out')
	filename_out = config.get('DEFAULT','filename_out')
	delim_out = gd("unicode_escape")(config.get('DEFAULT','delim_in'))[0] # treat delimiters correctly: ignore first escape

	# perform checks on the validity of the parameters and raise exceptions if necessary
	# DO NOT change the value of variables here!
	# TODO the 'is None' checks are obsolete. remove them (keep the error messages for later, now).
	if not path.exists(file_in):
		raise FileNotFoundError("File "+file_in+" not found.")
	if not (len(delim_in) == 1 and isinstance(delim_in, str)):
		raise Exception("Delimiter of input file must be a character (string of length one).")
	if not ((isinstance(header_in, int) and header_in >= 0) or header_in is None):
		raise Exception("Header parameter of the input file must be a non-negative integer or of type None.")
	if removeIsolationInterference_bool is None:
		raise Exception("Please indicate whether you would like to remove high Isolation Interference cases.")
	if not (0 < removeIsolationInterference_threshold < 100 or removeIsolationInterference_bool is None):
		raise Exception("Isolation Interference Threshold should be either 'None' or between 0 and 100 (percentage).")
	if collapsePSMAlgo_bool is None:
		raise Exception("Please indicate whether you would like to remove redundancy due to multiple PSM Algorithms.")
	if collapsePSMAlgo_master not in ('mascot', 'sequest'):
		raise Exception("Invalid master PSM algorithm: '"+collapsePSMAlgo_master+"'. Please pick 'mascot' or 'sequest'.")
	if collapsePSMAlgo_exclusive_bool is None:
		raise Exception("Please indicate whether PSM Algorithm redundancy removal should be exclusive or not.")
	if collapseRT_bool is None:
		raise Exception("Please indicate whether you would like to remove redundancy due to multiple retention times.")
	if collapseRT_centerMeasure_channels not in ('mean', 'median'):
		raise Exception("Invalid center measure: '"+collapseRT_centerMeasure_channels+"'. Please pick 'mean' or 'median'.")
	if collapseRT_centerMeasure_intensities not in ('max', 'mean', 'median'):
		raise Exception("Invalid center measure: '"+collapseRT_centerMeasure_channels+"'. "
		                                                                              "Please pick 'max', 'mean' or 'median'.")
	if collapseRT_maxRelativeReporterVariance is not None:
		if not collapseRT_maxRelativeReporterVariance > 0:
			raise Exception("maxRelativeChannelVariance should be either 'None' or greater than zero.")
	if collapseCharge_bool is None:
		raise Exception("Please indicate whether you would like to remove redundancy due to multiple charge states.")
	if not (isotopicCorrectionsMatrix.shape == (6,6)):
		raise Exception("Isotopic corrections matrix must have shape (6,6).")
	if not (accuracy > 0):
		raise Exception("Accuracy must be strictly greater than zero.")
	if not (maxIterations > 0 and isinstance(maxIterations,int)):
		raise Exception("Maximum number of iterations must be an integer strictly greater than zero.")
	if not path.exists(path_out):
		raise FileNotFoundError("Path " + path_out + " not found.")
	if path.exists(path_out+'/'+filename_out):
		warnings.warn("Will overwrite file "+path.basename(path.normpath(path_out)))
	if not (len(delim_out) == 1 and isinstance(delim_out, str)):
		raise Exception("Delimiter of output file must be a character (string of length one).")

	# assign the TYPOGRAPHICALLY CORRECT values to the params dict.
	params = {
		'file_in': file_in,
		'delim_in': delim_in,
		'header_in': header_in,
		'removeIsolationInterference_bool': removeIsolationInterference_bool,
		'removeIsolationInterference_threshold': removeIsolationInterference_threshold,
		'collapsePSMAlgo_bool': collapsePSMAlgo_bool,
		'collapsePSMAlgo_master': collapsePSMAlgo_master,
		'collapsePSMAlgo_exclusive_bool': collapsePSMAlgo_exclusive_bool,
		'collapseRT_bool': collapseRT_bool,
		'collapseRT_centerMeasure_channels': collapseRT_centerMeasure_channels,
		'collapseRT_centerMeasure_intensities': collapseRT_centerMeasure_intensities,
		'collapseRT_maxRelativeReporterVariance': collapseRT_maxRelativeReporterVariance,
		'collapseCharge_bool': collapseCharge_bool,
		'isotopicCorrectionsMatrix': isotopicCorrectionsMatrix,
		'accuracy': accuracy,
		'maxIterations': maxIterations,
		'DEFoldThreshold': DEFoldThreshold,
		'path_out': path_out,
		'filename_out': filename_out,
		'delim_out': delim_out
	}

	# check if you forgot to hardcode new parameters
	if not dontForgetAnyParameters == len(params):
		raise Exception("Number of parameters in config.ini not equal to number of parameters returned by getInput().")

	return params


def importDataFrame(path_in=None, delim=None, header=0):
	"""
	Get the data from disk as a Pandas DataFrame.
	:param path_in:     string          existing path to input file
	:param filetype:    string          specifier for the type of the file (file extension)
	:param delim:       char            delimiter of the data
	:return df:         pd.dataFrame    Pandas dataFrame of the file contents
	"""
	assert path.exists(path_in)

	if delim is None:
		extension = path_in.split('.')[-1] # set extension equal to the file extension (can return None)
		delim = ext2delim(extension)
		if delim is None:
			raise Exception(
				"Cannot handle data: filetype not recognized and no delimiter specified.")

	if delim=='xlsx':
		df = pd.read_excel(path_in)
	else:
		df = pd.read_csv(path_in, delimiter=delim, header=header)

	if delim is None:
		raise Exception(
			"I don't know how to handle this data: the filetype was not recognized and no delimiter was specified.")

	# df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10, 6)), columns=list('ABCDEF'))  # TEST
	# df = pd.read_csv('../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.txt', delim='\t') # TEST
	# df = pd.DataFrame(np.arange(10*6).reshape(10,6),columns=list('ABCDEF')) # TEST
	# df['B'][0]=np.nan # TEST
	# df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(2*10**3, 6)), columns=list('ABCDEF'))  # TEST

	return df


def exportData(data, dataType, path_out, filename, delim_out=','):
	# TODO is path_in the complete path including the filename? If so, only one file can be exported (or you can choose
	# to automatically put the other files alongside it in the same root).
	"""
	Save the results (normalized intensities) to disk.
	:param data:        obj     data object to be exported to disk
	:param path_out:    string  path where data should be exported to
	:param filename:    string  filename for the data
	:param delim_out:       char    delimiter of the data
	"""
	# assert data is not None # TODO
	assert path.exists(path_out)

	extension = delim2ext(delim_out)
	fullPath = path_out + '/' + filename + extension

	if dataType == 'txt':
		np.savetxt(fullPath, data, delimiter=delim_out) # TODO
	elif dataType == 'df':
		data.to_csv(fullPath, sep=delim_out, index=False)


def delim2ext(delim):
	"""
	Returns the filename extension belonging to the specified delimiter, or returns empty string.
	"""
	if delim==',':
		return '.csv'
	elif delim=='\t':
		return '.tsv'
	else:
		return ''


def ext2delim(ext):
	"""
	Returns the delimiter belonging to the specified filename extension, or returns empty string.
	"""
	ext.lstrip('.')
	if ext=='csv':
		return ','
	elif ext=='tsv':
		return '\t'
	else:
		return None
