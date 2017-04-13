#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Handle all I/O of data files and parameters to and from both the workflow and the main dataFrame.
"""

import os
import pandas as pd
import numpy as np
import pickle
import configparser
from json import dumps
from warnings import warn


def importDataFrame(path_in, delim=None, header=0, dtype=None):
	"""
	Get the data from disk as a Pandas DataFrame.
	:param path_in:     string          existing path to input file
	:param delim:       char            delimiter of the data fields
	:param header:      int             row that contains the header of the data (None if no header)
	:param dtype:       dict/object		specifies as which type each column should be interpreted
	:return df:         pd.dataFrame    Pandas dataFrame of the file contents
	"""
	assert os.path.exists(path_in)

	if delim is None:
		if '.' in path_in: # file has an extension
			extension = path_in.split('.')[-1] # set extension equal to the file extension (can return None)
			delim = ext2delim(extension)
		else:
			warn("No file extension nor delimiter specified; Pandas will try to automatically detect the delimiter.")

	if delim == 'xlsx':
		df = pd.read_excel(path_in)
	else: # delim is something else OR None.
		try:
			df = pd.read_csv(path_in, delimiter=delim, header=header, dtype=dtype)
		except:
			if delim is None:
				raise Exception("Data cannot be read: no delimiter specified and Pandas failed automatic recognition.")
			else:
				raise Exception("Data cannot be read: the delimiter "+str(delim)+" is not right for this file.")

	return df.dropna(how="all") # drop empty lines


def importIsotopicCorrectionsMatrix(path_in):
	"""
	Reads the isotopic corrections matrix from a file on disk through importDataFrame, and returns it as a matrix.
	:param path_in: str         path of the isotopic corrections matrix file
	:return :    	np.ndarray  isotopic corrections matrix
	"""
	return np.asmatrix(importDataFrame(path_in, delim='\t', header=None)).astype('float64') # make sure its float64


def getTMTIsotopicDistributions(path_in):
	"""
	Read the TMT isotope distributions table (IDT) in tsv format from disk as a dataframe with the correct column headers and
	labels as specified in the file. The column with the channel names and the "IDT" header is removed.
	:param path_in:	str				path to the .tsv file with the TMT IDT
	:return tmtid:	pd.DataFrame	TMT isotope distributions table. Format:
											IDT     -2  -1  +1  +2
											126     0   0   1.2 0
											127N    1.2 3.3 2.5 0.3
											127C    ...
											...
	"""
	tmtid = importDataFrame(path_in)
	if 'IDT' in tmtid.columns.values:
		tmtid.set_index('IDT', drop=True, inplace=True)
	else:
		raise Exception("Column header of the channels in the TMT isotopic distributions .tsv file should be 'IDT'.")
	if set(tmtid.columns.values) != set(['-2','-1','+1','+2']):
		raise Exception("TMT isotopic distributions .tsv file should contain columns '[-2, -1, +1, +2]'.")
	tmtid.index = tmtid.index.values.astype(str)
	return tmtid


def importWrapper(path_in='wrapper.tsv'):
	"""
	Reads the column header wrapper from a file on disk through importDataFrame, and returns it as a nested list.
	:param path_in: str		path of the wrapper file
	:return :       list	wrapper (nested list of pairs) specifying column name transformations
	"""
	return list(importDataFrame(path_in, header=None, dtype=str).values)


def parseSchemaFile(schemaPath):
	"""
	Parses the .tsv schema into a hierarchical overview with intensity columns groups per condition and experiment. The
	wrapper and config entries are set to None for now.
	!!! the schema is NEVER to be changed after it has been first used !!!
		"Keys and values are iterated over in an arbitrary order which is non-random, varies across Python implementations,
		and depends on the dictionary’s history of insertions and deletions. If keys, values and items views are iterated
		over with no intervening modifications to the dictionary, the order of items will directly correspond."
	:param schemaPath:              str     path to the schema file that the user uploaded
	:return incompleteSchemaDict:   dict    schema in dict format, without config and wrapper information, in the format
											{ experiment: { channels: [[channels] per condition], aliases: [[channels] per condition] }
	"""
	from collections import OrderedDict
	# import schema file as dataframe and replace nan values by empty strings
	schemaDF = importDataFrame(schemaPath, delim='\t', header=None, dtype=str).replace(np.nan, '', regex=True)
	incompleteSchemaDict = OrderedDict() # use ordered dict so the items() order is always the same
	assert np.mod(len(schemaDF), 2) == 0 # schema must have even number of lines
	assert len(schemaDF.columns) > 2 # at least 3 columns
	for i in range(int(len(schemaDF)/2)):
		thisRow = schemaDF.loc[2*i, :]
		experimentName = str(thisRow[0])
		numConditions = len(thisRow)-1
		channelNamesPerCondition = [str(thisRow[c]).split(',') for c in range(1, numConditions+1)]
		nextRow = schemaDF.loc[2*i+1, :]
		# if the next row contains just as many values (that are not empty); aliases are provided
		if len(nextRow) == len(thisRow) and np.sum(nextRow == '') == 0:
			channelAliasesPerCondition = [str(nextRow[c]).split(',') for c in range(1, numConditions+1)]
		else: # aliases not (properly) provided
			shortExperimentName = str(nextRow[0])
			channelAliasesPerCondition = [[shortExperimentName + '_' + channelName for channelName in condition]
										  for condition in channelNamesPerCondition]

		incompleteSchemaDict[experimentName] = {'channelNamesPerCondition': channelNamesPerCondition,
									  'channelAliasesPerCondition': channelAliasesPerCondition}
	return incompleteSchemaDict


# def writeConfig(filePath, contents): # todo remove this
# 	"""
#
# 	:param filePath:
# 	:param contents:
# 	:return:
# 	"""
# 	file = open(filePath, 'w')
# 	config = configparser.ConfigParser()
# 	section = 'DEFAULT'
# 	# config.add_section(section)
# 	for k, v in contents.items():
# 		config.set(section, k, v)
# 	config.write(file)
# 	file.close()


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
	Takes a Nx2 tuple/list wrapper and transforms the columns in the dataframe df specified by the first column entries
	in wrapper into the new column name specified by the second colmumn entries in wrapper.
	:param columns:     pd.Index        old column names
	:param wrapper:     list(tuples)    [(oldName, newName) for some columns]
	:return newColumns: pd.Index        transformed column names
	"""
	newColumns = list(columns.values)
	for (old, new) in wrapper:
		newColumns[newColumns.index(old)] = new
	return newColumns


def importExperimentData(path_in, delim=None, header=0, wrapper=None):
	"""
	Gets the experimental data specified by path_in from disk, applies a wrapper (optional) and fixes common
	format mistakes.
	:param path_in:     string          existing path to input file
	:param delim:       char            delimiter of the data
	:param header:      int             row that contains the header of the data (None if no header)
	:return df:			pd.DataFrame    ready-to-use experimental data
	"""
	df = importDataFrame(path_in, delim=delim, header=header)
	df.columns = applyWrapper(df.columns, wrapper)
	df = fixFixableFormatMistakes(df)
	return df


def exportData(data, dataType, path_out, filename, delim_out=None, inOneFile=False):
	"""
	Save any type of output to disk.
	:param data:		obj     data object to be exported to disk
	:param dataType:	str     type of data; influences the way the data is to be stored
	:param path_out:	string  path where data should be exported to
	:param filename:	string  filename for the data
	:param delim_out:	char    delimiter of the data
	"""
	# todo this function should call sub-functions per data type and optional arguments
	assert os.path.exists(path_out)

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
	elif dataType == 'fig':
		outFileFullPathNoExt = os.path.join(path_out, filename)
		with open(outFileFullPathNoExt+'.pkl', "wb") as fout:
			# save as pickle. Load again later as: ax = pickle.load(file('myplot.pickle')); then plt.show()
			pickle.dump(data, fout, protocol=4)
		from matplotlib import pyplot as plt
		#ax = data
		plt.savefig(outFileFullPathNoExt+'.png', format='png', bbox_inches='tight')
		return outFileFullPathNoExt+'.png'
	elif dataType == 'html':
		fullPath += '.html'
		with open(fullPath, 'w') as htmlFile:
			htmlFile.writelines(data)
		return fullPath


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
