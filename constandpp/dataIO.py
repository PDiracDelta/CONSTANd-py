#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Handle all I/O of data files and parameters to and from both the workflow and the main dataFrame.
"""

import pandas as pd
import numpy as np
import pickle
from os import path


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

	if delim == 'xlsx':
		df = pd.read_excel(path_in)
	else:
		df = pd.read_csv(path_in, delimiter=delim, header=header)

	if delim is None:
		raise Exception(
			"I don't know how to handle this data: the filetype was not recognized and no delimiter was specified.")
	return df


def getIsotopicCorrectionsMatrix(path_in='ICM_default.tsv'):
	"""
	Reads the isotopic corrections matrix from a file on disk through importDataFrame, and returns it as a matrix.
	:param path_in: str         path of the isotopic corrections matrix file
	:return icm:    pd.ndarray  isotopic corrections matrix
	"""
	return np.asmatrix(importDataFrame(path_in,delim='\t', header=None)).astype('float64') # make sure its float64

#
# def parseSchema(schemaPath):
# 	"""
# 	Parses the .tsv schema into a hierarchical overview with intensity columns groups per condition and experiment
# 	:param schemaPath:
# 	:return:
# 	"""
# 	schemaDF = importDataFrame(schemaPath, delim='\t', header=None)
# 	schemaDict = None
# 	return schemaDict


def fixFixableFormatMistakes(df):
	"""
	Takes a dataframe and fixes format mistakes frequently present in PD2.1 output.
	:param df:  pd.DataFrame    possibly containing a multitude of mistakes.
	:return df: pd.DataFrame    data without recognized format mistakes.
	"""
	# you've enabled "show flanking amino acids": DIRk --> [L].DIRk.[m]
	if df.sample(n=1)['Annotated Sequence'].item().count('.') == 2: # sequence contains 2 dots
		df['Annotated Sequence'] = df['Annotated Sequence'].apply(lambda x: x.split('.')[1]) # select part between dots

	return df


def applyWrapper(df, wrapper):
	# todo
	return df


def getDataFrame(path_in, delim=None, header=0, wrapper=None):
	"""
	Gets the dataframe specified by path_in from disk, applies a wrapper (optional) and fixes common format mistakes.
	:param path_in:     string          existing path to input file
	:param delim:       char            delimiter of the data
	:param header:      int             row that contains the header of the data (None if no header)
	:return:            pd.DataFrame    ready-to-use data
	"""
	df = importDataFrame(path_in, delim=delim, header=header)
	df = applyWrapper(df, wrapper) #todo
	df = fixFixableFormatMistakes(df)
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
