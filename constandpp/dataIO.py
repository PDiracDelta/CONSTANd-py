#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Handle all I/O of data files and parameters to and from both the workflow and the main dataFrame.
"""

import pandas as pd
import numpy as np
import pickle
import configparser
from json import dumps
from os import path
import re


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


def getWrapper(path_in='wrapper.tsv'):
	"""
	Reads the wrapper from a file on disk through importDataFrame, and returns it as a list of tuples.
	:param path_in: str            path of the wrapper file
	:return :       nested list    wrapper specifying column name transformations
	"""
	return list(importDataFrame(path_in, header=None).values)


def parseSchemaFile(schemaPath): #todo move to web
	"""
	Parses the .tsv schema into a hierarchical overview with intensity columns groups per condition and experiment. The
	wrapper and config entries are set to None for now.
	!!! the schema is NEVER to be changed after it has been first used !!!
		" Keys and values are iterated over in an arbitrary order which is non-random, varies across Python implementations,
		and depends on the dictionary’s history of insertions and deletions. If keys, values and items views are iterated
		over with no intervening modifications to the dictionary, the order of items will directly correspond."
	:param schemaPath:              str     path to the schema file that the user uploaded
	:return incompleteSchemaDict:   dict    schema in dict format, without config and wrapper information, in the format
											{ experiment: { channels: [[channels] per condition], aliases: [[channels] per condition] }
	"""
	# import schema file as dataframe and replace nan values by empty strings
	schemaDF = importDataFrame(schemaPath, delim='\t', header=None).replace(np.nan, '', regex=True)
	incompleteSchemaDict = {}
	assert np.mod(len(schemaDF), 2) == 0 # schema must have even number of lines
	for i in range(int(len(schemaDF)/2)):
		thisRow = schemaDF.loc[2*i, :]
		experimentName = str(thisRow[0])
		numConditions = len(thisRow)-1
		intensityColumnsPerCondition = [str(thisRow[c]).split(',') for c in range(1, numConditions+1)]
		nextRow = schemaDF.loc[2*i+1, :]
		# if the next row contains just as many values (that are not empty); aliases are provided
		if len(nextRow) == len(thisRow) and np.sum(nextRow == '') == 0:
			channelAliasesPerCondition = [str(nextRow[c]).split(',') for c in range(1, numConditions+1)]
		else: # aliases not (properly) provided
			shortExperimentName = str(nextRow[0])
			channelAliasesPerCondition = [[shortExperimentName + '_' + channelName for channelName in condition]
			                              for condition in intensityColumnsPerCondition]

		incompleteSchemaDict[experimentName] = {'intensityColumnsPerCondition': intensityColumnsPerCondition,
		                              'channelAliasesPerCondition': channelAliasesPerCondition,
		                              'config': None, 'wrapper': None}
	return incompleteSchemaDict


def constructMasterConfigContents(schemaDict, otherMasterParams): # todo move to web
	# todo docu
	def isNumeric(s):
		""" Returns whether or not the argument is numeric """
		try:
			float(s)
			return True
		except ValueError:
			return False
	contents = {}
	contents['schema'] = dumps(schemaDict)
	for k, v in otherMasterParams.items():
		if isinstance(v, str) or isNumeric(v):
			if k.split('_')[0] == 'delim': # delimiters should be saved in a visible format
				from codecs import getencoder as ge
				contents[k] = ge("unicode_escape")("\\"+v) # todo doesnt work
			else: # not a delimiter
				contents[k] = v
		else:
			contents[k] = dumps(v)

	return contents


def writeConfig(filePath, contents):
	file = open(filePath, 'w')
	config = configparser.ConfigParser()
	section = 'DEFAULT'
	# config.add_section(section)
	for k, v in contents.items():
		config.set(section, k, v)
	config.write(file)
	file.close()


def fixFixableFormatMistakes(df):
	"""
	Takes a dataframe and fixes format mistakes frequently present in PD2.1 output.
	:param df:  pd.DataFrame    possibly containing a multitude of mistakes.
	:return df: pd.DataFrame    data without recognized format mistakes.
	"""
	# you've enabled "show flanking amino acids": DIRk --> [L].DIRk.[m]
	if df.sample(n=1)['Annotated Sequence'].item().count('.') == 2: # sequence contains 2 dots
		df['Annotated Sequence'] = df['Annotated Sequence'].apply(lambda x: x.split('.')[1]) # select part between dots

	# you're using "Identifying Node" instead of "Identifying Node Type"
	if 'Identifying Node' in df.columns.values and 'Identifying Node Type' not in df.columns.values:
		# strip the everything after the last space (=node number between parentheses) + the last space.
		df['Identifying Node Type'] = df['Identifying Node'].apply(lambda s: s.rsplit(' ', 1)[0])

	return df


def applyWrapper(columns, wrapper):
	"""
	Takes a Nx2 tuple wrapper and transforms the columns in the dataframe df specified by the first column entries
	in wrapper into the new column name specified by the second colmumn entries in wrapper.
	:param columns:     pd.Index        old column names
	:param wrapper:     list(tuples)    [(oldName, newName) for some columns]
	:return newColumns: pd.Index        transformed column names
	"""
	newColumns = list(columns.values)
	for (old, new) in wrapper:
		newColumns[newColumns.index(old)] = new
	return newColumns


def getDataFrame(path_in, delim=None, header=0, wrapper=None):
	"""
	Gets the dataframe specified by path_in from disk, applies a wrapper (optional) and fixes common format mistakes.
	:param path_in:     string          existing path to input file
	:param delim:       char            delimiter of the data
	:param header:      int             row that contains the header of the data (None if no header)
	:return:            pd.DataFrame    ready-to-use data
	"""
	df = importDataFrame(path_in, delim=delim, header=header)
	df.columns = applyWrapper(df.columns, wrapper) #todo
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


def unnest(x):
	""" returns un-nested version of level 1 nested list x."""
	return [e for sublist in x for e in sublist]