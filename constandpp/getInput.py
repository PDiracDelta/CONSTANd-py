#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Get the input parameters from the config file.
"""

import configparser
import numpy as np
from os import path
from json import loads as getList
from codecs import getdecoder as gd
from dataIO import getIsotopicCorrectionsMatrix
from warnings import warn


def getInput(configFilePath):
	"""
	Get mass spec data and CONSTANd parameters from the user or from the web interface as a dict.
	"""
	# TODO add all parameters in docstring
	# TODO add .lower() to all string input except requiredColumns and intensityColumns

	# read the config file to obtain the defaults
	config = configparser.ConfigParser(allow_no_value=True, comment_prefixes=';',
	                                   inline_comment_prefixes='@')  # TODO split up CONFIG and PARAMS (user input files vs workflow params)
	config.optionxform = str # so that strings dont automatically get .lower()-ed
	config.read(configFilePath, encoding='utf-8')

	# get variables from config in correct typography
	date = config.get('DEFAULT','date')
	files_in = getList(config.get('DEFAULT','files_in'))
	delim_in = gd("unicode_escape")(config.get('DEFAULT','delim_in'))[0] # treat delimiters correctly: ignore first escape
	header_in = config.getint('DEFAULT','header_in')
	removedDataInOneFile_bool = config.getboolean('DEFAULT','removedDataInOneFile_bool')
	intensityColumnsPerCondition = getList(config.get('DEFAULT', 'intensityColumnsPerCondition'))
	wantedColumns = getList(config.get('DEFAULT', 'wantedColumns'))
	noMissingValuesColumns = getList(config.get('DEFAULT', 'noMissingValuesColumns'))
	removalColumnsToSave = getList(config.get('DEFAULT', 'removalColumnsToSave'))
	collapseColumnsToSave = getList(config.get('DEFAULT', 'collapseColumnsToSave'))
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
	isotopicCorrection_matrix = getIsotopicCorrectionsMatrix(config.get('DEFAULT','isotopicCorrection_matrix'))
	accuracy = config.getfloat('DEFAULT','accuracy')
	maxIterations = config.getint('DEFAULT','maxIterations')
	pept2protCombinationMethod = config.get('DEFAULT','pept2protCombinationMethod')
	alpha = config.getfloat('DEFAULT','alpha')
	FCThreshold = config.getfloat('DEFAULT','FCThreshold')
	labelVolcanoPlotAreas = getList(config.get('DEFAULT','labelVolcanoPlotAreas'))
	PCA_components = config.getint('DEFAULT','PCA_components')
	path_out = config.get('DEFAULT','path_out')
	filename_out = config.get('DEFAULT','filename_out')
	delim_out = gd("unicode_escape")(config.get('DEFAULT','delim_in'))[0] # treat delimiters correctly: ignore first escape

	# perform checks on the validity of the parameters and raise exceptions if necessary
	# DO NOT change the value of variables here!
	# TODO the 'is None' checks are obsolete. remove them (keep the error messages for later, now).
	if not path.exists(files_in[0]): # TODO for all files
		raise FileNotFoundError("File "+files_in+" not found.")
	if not (len(delim_in) == 1 and isinstance(delim_in, str)):
		raise Exception("Delimiter of input file must be a character (string of length one).")
	if not ((isinstance(header_in, int) and header_in >= 0) or header_in is None):
		raise Exception("Header parameter of the input file must be a non-negative integer or of type None.")
	if intensityColumnsPerCondition is None:
		raise Exception("Please indicate which columns contain the MS2 reporter intensities.")
	if wantedColumns is None:
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
	if not (isotopicCorrection_matrix.shape[0] == isotopicCorrection_matrix.shape[1]):
		raise Exception("Isotopic corrections matrix must have square shape. AND EQUAL TO NUMBER OF CHANNELS")
	if not (np.allclose(np.sum(isotopicCorrection_matrix,0),np.ones(isotopicCorrection_matrix.shape[0]),atol=1e-9)): # absolute tolerance: intensities known up to ~1e-10
		raise Exception("Isotopic corrections matrix row values do not add up to 1.")
	if np.linalg.det(isotopicCorrection_matrix) == 0: # if Det(cM) = 0 no solution can be found.
		raise Exception("Determinant of isotopic corrections matrix is zero; cannot solve the linear system.")
	if not (accuracy > 0):
		raise Exception("Accuracy must be strictly greater than zero.")
	if not (maxIterations > 0 and isinstance(maxIterations,int)):
		raise Exception("Maximum number of iterations must be an integer strictly greater than zero.")
	if PCA_components < 2:
		raise Exception("Minimum number of principal coponents is 2.")
	if not path.exists(path_out):
		raise FileNotFoundError("Path " + path_out + " not found.")
	if path.exists(path_out+'/'+filename_out):
		warn("Will overwrite file "+path.basename(path.normpath(path_out)))
	if not (len(delim_out) == 1 and isinstance(delim_out, str)):
		raise Exception("Delimiter of output file must be a character (string of length one).")

	# assign the TYPOGRAPHICALLY CORRECT values to the params dict and modify them if necessary.
	# modify
	intensityColumns = [item for sublist in intensityColumnsPerCondition for item in sublist]
	# assign
	params = {
		'date': date,
		'files_in': files_in,
		'delim_in': delim_in,
		'header_in': header_in,
		'removedDataInOneFile_bool': removedDataInOneFile_bool,
		'intensityColumnsPerCondition': intensityColumnsPerCondition,
		'intensityColumns': intensityColumns,
		'wantedColumns': wantedColumns+intensityColumns, # needs to include intensitycolumns
		'noMissingValuesColumns': noMissingValuesColumns,
		'removalColumnsToSave': removalColumnsToSave+intensityColumns, # needs to include intensitycolumns
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
		'isotopicCorrection_matrix': isotopicCorrection_matrix,
		'accuracy': accuracy,
		'maxIterations': maxIterations,
		'pept2protCombinationMethod': pept2protCombinationMethod,
		'alpha': alpha,
		'FCThreshold': FCThreshold,
		'labelVolcanoPlotAreas': labelVolcanoPlotAreas,
		'PCA_components': PCA_components,
		'path_out': path_out,
		'filename_out': filename_out,
		'delim_out': delim_out
	}

	# check if you forgot to hardcode new parameters
	for param in config._defaults.keys():
		if param not in params.keys():
			raise Exception("You forgot to include "+param+" in the params dictionary.")

	return params