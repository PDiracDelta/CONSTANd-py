#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Handle all I/O of data files and parameters to and from both the workflow and the main dataFrame.
"""

import pandas as pd
import numpy as np
import configparser
import pickle
from os import path
from codecs import getdecoder as gd
from warnings import warn


def parseList(listInString):
	"""
	Takes a comma-separated list hidden in a string and returns it as a list.
	:param listInString:    str     comma-separated list hidden in a string
	:return:                list    the list that was hidden in the string
	"""
	return [x.strip(' ') for x in listInString.split(',')]


def getInput():
	"""
	Get mass spec data and CONSTANd parameters from the user or from the web interface as a dict.
	"""
	# TODO add all parameters in docstring
	# TODO add .lower() to all string input except requiredColumns and intensityColumns

	# read the config file to obtain the defaults
	config = configparser.ConfigParser(allow_no_value=True, comment_prefixes=';',
	                                   inline_comment_prefixes='@')  # TODO split up CONFIG and PARAMS (user input files vs workflow params)
	config.optionxform = str # so that strings dont automatically get .lower()-ed
	config.read('config.ini', encoding='utf-8')

	# get variables from config in correct typography
	file_in = config.get('DEFAULT','file_in')
	delim_in = gd("unicode_escape")(config.get('DEFAULT','delim_in'))[0] # treat delimiters correctly: ignore first escape
	header_in = config.getint('DEFAULT','header_in')
	removedDataInOneFile_bool = config.getboolean('DEFAULT','removedDataInOneFile_bool')
	intensityColumns = parseList(config.get('DEFAULT', 'intensityColumns'))
	requiredColumns = parseList(config.get('DEFAULT', 'requiredColumns'))
	noMissingValuesColumns = parseList(config.get('DEFAULT', 'noMissingValuesColumns'))
	remove_ExtraColumnsToSave = parseList(config.get('DEFAULT', 'remove_ExtraColumnsToSave'))
	collapseColumnsToSave = parseList(config.get('DEFAULT', 'collapseColumnsToSave'))
	removeBadConfidence_bool = config.getboolean('DEFAULT','removeBadConfidence_bool')
	removeBadConfidence_minimum = config.get('DEFAULT','removeBadConfidence_minimum')
	removeIsolationInterference_bool = config.getboolean('DEFAULT','removeIsolationInterference_bool')
	removeIsolationInterference_threshold = config.getfloat('DEFAULT','removeIsolationInterference_threshold')
	collapse_method = config.get('DEFAULT', 'collapse_method')
	collapse_maxRelativeReporterVariance = config.getfloat('DEFAULT', 'collapse_maxRelativeReporterVariance')
	undoublePSMAlgo_bool = config.getboolean('DEFAULT','undoublePSMAlgo_bool')
	masterPSMAlgo = config.get('DEFAULT','masterPSMAlgo')
	undoublePSMAlgo_exclusive_bool = config.getboolean('DEFAULT','undoublePSMAlgo_exclusive_bool')
	collapseCharge_bool = config.getboolean('DEFAULT','collapseCharge_bool')
	collapsePTM_bool = config.getboolean('DEFAULT','collapsePTM_bool')
	isotopicCorrection_bool = config.getboolean('DEFAULT','isotopicCorrection_bool')
	isotopicCorrectionsMatrix = getIsotopicCorrectionsMatrix(config.get('DEFAULT','isotopicCorrectionsMatrix'))
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
	if intensityColumns is None:
		raise Exception("Please indicate which columns contain the MS2 reporter intensities.")
	if requiredColumns is None:
		raise Exception("Please indicate which columns (in addition to the intensities) you would like to have output for.")
	if collapseColumnsToSave is None:
		raise Exception("Please indicate which columns (in addition to the intensities) to save for removed data.")
	if removeBadConfidence_bool is None:
		raise Exception("Please indicate whether you would like to remove detections with confidence lower than certain "
		                "threshold.")
	if removeBadConfidence_bool and removeBadConfidence_minimum not in ['High', 'Medium']:
		raise Exception("Invalid minimum confidence level: "+removeBadConfidence_minimum+". Must select 'Medium' or 'High'.")
	if removeIsolationInterference_bool is None:
		raise Exception("Please indicate whether you would like to remove high Isolation Interference detections.")
	if not (0 < removeIsolationInterference_threshold < 100 or removeIsolationInterference_bool is None):
		raise Exception("Isolation Interference Threshold should be either 'None' or between 0 and 100 (percentage).")
	if collapse_method not in ('bestMatch', 'mostIntense', 'mean', 'geometricMedian', 'weighted'):
		raise Exception("Invalid collapse method: '"+collapse_method+"'. Please pick 'max', 'mean' or 'median'.")
	if collapse_maxRelativeReporterVariance is not None:
		if not collapse_maxRelativeReporterVariance > 0:
			raise Exception("maxRelativeChannelVariance should be either 'None' or greater than zero.")
	if undoublePSMAlgo_bool is None:
		raise Exception("Please indicate whether you would like to remove redundancy due to multiple PSM Algorithms.")
	if masterPSMAlgo is None:
		raise Exception("You have to enter a masterPSMAlgo -- even if you don't want to undouble the PSMs -- in order to "
		                "choose the best matching representative when collapsing")
	if masterPSMAlgo not in ('mascot', 'sequest'):
		raise Exception("Invalid master PSM algorithm: '"+masterPSMAlgo+"'. Please pick 'mascot' or 'sequest'.")
	if undoublePSMAlgo_exclusive_bool is None:
		raise Exception("Please indicate whether PSM Algorithm redundancy removal should be exclusive or not.")
	if collapseCharge_bool is None:
		raise Exception("Please indicate whether you would like to remove redundancy due to multiple charge states.")
	if collapsePTM_bool is None:
		raise Exception("Please indicate whether you would like to remove redundancy due to multiple PTMs.")
	if isotopicCorrection_bool is None:
		raise Exception("Please indicate whether you would like to correct for isotopic impurities.")
	if not (isotopicCorrectionsMatrix.shape == (6,6)):
		raise Exception("Isotopic corrections matrix must have shape (6,6).")
	if not (np.allclose(np.sum(isotopicCorrectionsMatrix,0),np.ones(6),atol=1e-9)): # absolute tolerance: intensities known up to ~1e-10
		raise Exception("Isotopic corrections matrix row values do not add up to 1.")
	if np.linalg.det(isotopicCorrectionsMatrix) == 0: # if Det(cM) = 0 no solution can be found.
		raise Exception("Determinant of isotopic corrections matrix is zero; cannot solve the linear system.")
	if not (accuracy > 0):
		raise Exception("Accuracy must be strictly greater than zero.")
	if not (maxIterations > 0 and isinstance(maxIterations,int)):
		raise Exception("Maximum number of iterations must be an integer strictly greater than zero.")
	if not path.exists(path_out):
		raise FileNotFoundError("Path " + path_out + " not found.")
	if path.exists(path_out+'/'+filename_out):
		warn("Will overwrite file "+path.basename(path.normpath(path_out)))
	if not (len(delim_out) == 1 and isinstance(delim_out, str)):
		raise Exception("Delimiter of output file must be a character (string of length one).")

	# assign the TYPOGRAPHICALLY CORRECT values to the params dict and modify them if necessary.
	params = {
		'file_in': file_in,
		'delim_in': delim_in,
		'header_in': header_in,
		'removedDataInOneFile_bool': removedDataInOneFile_bool,
		'intensityColumns': intensityColumns,
		'requiredColumns': requiredColumns+intensityColumns, # needs to include intensitycolumns
		'noMissingValuesColumns': noMissingValuesColumns,
		'remove_ExtraColumnsToSave': remove_ExtraColumnsToSave+intensityColumns, # needs to include intensitycolumns
		'collapseColumnsToSave': collapseColumnsToSave+intensityColumns, # needs to include intensitycolumns
		'removeBadConfidence_bool': removeBadConfidence_bool,
		'removeBadConfidence_minimum': removeBadConfidence_minimum,
		'removeIsolationInterference_bool': removeIsolationInterference_bool,
		'removeIsolationInterference_threshold': removeIsolationInterference_threshold,
		'masterPSMAlgo': masterPSMAlgo,
		'undoublePSMAlgo_bool': undoublePSMAlgo_bool,
		'undoublePSMAlgo_exclusive_bool': undoublePSMAlgo_exclusive_bool,
		'collapse_method': collapse_method,
		'collapse_maxRelativeReporterVariance': collapse_maxRelativeReporterVariance,
		'collapseCharge_bool': collapseCharge_bool,
		'collapsePTM_bool': collapsePTM_bool,
		'isotopicCorrection_bool': isotopicCorrection_bool,
		'isotopicCorrectionsMatrix': isotopicCorrectionsMatrix,
		'accuracy': accuracy,
		'maxIterations': maxIterations,
		'DEFoldThreshold': DEFoldThreshold,
		'path_out': path_out,
		'filename_out': filename_out,
		'delim_out': delim_out
	}

	# check if you forgot to hardcode new parameters
	for param in config._defaults.keys():
		if param not in params.keys():
			raise Exception("You forgot to include "+param+" in the params dictionary.")

	return params


def getIsotopicCorrectionsMatrix(path_in='ICM_default.tsv'):
	"""
	Reads the isotopic corrections matrix from a file on disk through importDataFrame, and returns it as a matrix.
	:param path_in: str         path of the isotopic corrections matrix file
	:return icm:    pd.ndarray  isotopic corrections matrix
	"""
	return np.asmatrix(importDataFrame(path_in,delim='\t', header=None)).astype('float64') # make sure its float64


def importDataFrame(path_in, delim=None, header=0):
	"""
	Get the data from disk as a Pandas DataFrame.
	:param path_in:     string          existing path to input file
	:param delim:       char            delimiter of the data
	:param header:      int             row that contains the header of the data (None if no header)
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


def exportData(data, dataType, path_out, filename, delim_out=None, inOneFile=False):
	"""
	Save the results (normalized intensities) to disk.
	:param data:        obj     data object to be exported to disk
	:param dataType:    str     type of data; influences the way the data is to be stored
	:param path_out:    string  path where data should be exported to
	:param filename:    string  filename for the data
	:param delim_out:       char    delimiter of the data
	"""
	# assert data is not None # TODO
	assert path.exists(path_out)

	extension = delim2ext(delim_out)
	fullPath = path_out + '/' + filename + extension

	if dataType == 'txt':
		np.savetxt(fullPath, data, delimiter=delim_out)
	elif dataType == 'obj':
		pickle.dump(data, open(fullPath, 'wb'))
	elif dataType == 'df':
		if isinstance(data, dict): # there are actually multiple dataFrames
			if inOneFile: # save all removedData in one file.
				removedData=pd.DataFrame().append(list(data.values()))
				# data['missing'].columns contains all possible columns
				removedData.to_csv(path_out + '/' + filename + extension, sep=delim_out, index=False, columns=data['missing'].columns)
			else: # save all removedData in separate files per category.
				for frameName, frame in data.items():
					assert isinstance(frame, pd.DataFrame)
					frame.to_csv(path_out + '/' + filename + '_' + frameName + extension, sep=delim_out,index=False)
		else:
			assert isinstance(data, pd.DataFrame)
			data.to_csv(fullPath, sep=delim_out, index=False)
	elif dataType == 'viz': # TODO
		pass # plt.savefig('foo.png', bbox_inches='tight')


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
