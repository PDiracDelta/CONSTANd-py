#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Handle all I/O of data files and parameters to and from both the workflow and the main dataFrame.
"""

import os
import pandas as pd
import numpy as np
import pickle
import zipfile
import configparser
from json import dumps
from warnings import warn
import logging


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
		if '.' in path_in:  # file has an extension
			extension = path_in.split('.')[-1]  # set extension equal to the file extension (can return None)
			delim = ext2delim(extension)
		else:
			warn("No file extension nor delimiter specified; Pandas will try to automatically detect the delimiter.")
	
	if delim == 'xlsx':
		df = pd.read_excel(path_in)
	else:  # delim is something else OR None.
		try:
			df = pd.read_csv(path_in, delimiter=delim, header=header, dtype=dtype)
		except:
			if delim is None:
				raise Exception("Data cannot be read: no delimiter specified and Pandas failed automatic recognition.")
			else:
				raise Exception("Data cannot be read: the delimiter " + str(delim) + " is not right for this file.")
	
	return df.dropna(how="all")  # drop empty lines


def importIsotopicCorrectionsMatrix(path_in):
	"""
	Reads the isotopic corrections matrix from a file on disk through importDataFrame, and returns it as a matrix.
	:param path_in: str         path of the isotopic corrections matrix file
	:return :    	np.ndarray  isotopic corrections matrix
	"""
	# return np.asmatrix(importDataFrame(path_in, delim='\t', header=None)).astype('float64')  # make sure its float64
	return np.genfromtxt(path_in, delimiter='\t').astype('float64')  # make sure its float64


def importTMTIsotopicDistributions(path_in):
	"""
	Read the TMT isotope distributions table (IDT) in tsv format from disk as a dataframe with the correct column headers
	and labels as specified in the file. The column with the channel names and the "IDT" header is removed.
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
	if set(tmtid.columns.values) != set(['-2', '-1', '+1', '+2']):
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


def parseSchemaFile(schemaPath):  # todo either move this to web.py or redistribute file manipulation files(functions?) in web.py
	"""
	Parses the .tsv schema into a hierarchical overview with intensity columns groups per condition and experiment. The
	wrapper and config entries are set to None for now. The structure of the schema should be as follows:
		experiment_name	cond1:df_col1,df_col2:alias1,alias2		cond2:df_col3,df_col4:alias3,alias4
	IF aliases are provided for a certain row, they must all be provided. Else they are generated automatically as
	EXPERIMENTNAME_CONDITION_COLNAME.
	:param schemaPath:              str     path to the schema file that the user uploaded
	:return incompleteSchemaDict:   dict    schema in dict format, without config and wrapper information, in the format
											{ allExperiments: [experiments] ,
											  allConditions: [conditions] ,
											  experiment: {
											    allExperimentConditions: [conditions] ,
												allExperimentChannelNames: [channelNames] ,
												allExperimentChannelAliases: [channelAliases] ,
												{ condition: { channelNames: [names] , channelAliases: [aliases] } }
											  }
											}
	"""
	def extractAliases(rowElement):
		"""
		Tries to extract the comma-separated channelAliases from a schemaDF row element and returns '' if it fails.
		:param rowElement:			list	element in a schemaDF row. Structure: [ condition : channelNames : channelAliases ]
		:return channelAliasesList:	list	comma-separated list of [ channelAliases ]
		"""
		try:  # allow that no aliases are provided
			this_channelAliases = rowElement.split(':')[2]
		except IndexError:
			this_channelAliases = ''  # return '' because that's what is also returned if the substring is empty.
		return this_channelAliases
	
	""" import schema as dataframe """
	from collections import OrderedDict
	# import schema file as dataframe and replace nan values by empty strings
	schemaDF = importDataFrame(schemaPath, delim='\t', header=None, dtype=str).replace(np.nan, '', regex=True)
	incompleteSchemaDict = OrderedDict()  # use ordered dict so the items() order is always the same
	numCols = len(schemaDF.columns)
	numRows = len(schemaDF)
	if not (numCols > 2):  # at least 3 columns
		raise Exception("Schema must have at least 3 columns EXPERIMENT\\tCONDITION 1\\tCONDITION 2 (separated by tabs)")
	incompleteSchemaDict['allExperiments'] = list(schemaDF.loc[:, 0])
	# check if experiment names are unique
	if not len(set(incompleteSchemaDict['allExperiments'])) == numRows:
		raise Exception("Experiment names contain duplicates. Please provide unique names.")
	forbiddenExperimentNames = {'allExperiments', 'allConditions'}
	# check if no experiment names are forbidden
	if set(incompleteSchemaDict['allExperiments']) & forbiddenExperimentNames:
		raise Exception("Please do not use the following as experiment names: "+str(list(forbiddenExperimentNames)))
	
	""" construct schema dict """
	allConditions = set()
	allChannelNames = set()
	allChannelAliases = set()
	numAllChannelNames = 0
	numAllChannelAliases = 0
	for __, row in schemaDF.iterrows():
		experimentChannelNames = []
		experimentChannelAliases = []
		experimentName = str(row[0])
		incompleteSchemaDict[experimentName] = OrderedDict()
		conditionsList = [str(element).split(':')[0] for element in row[1:]]
		allConditions.update(conditionsList)
		try:
			channelNamesList = [str(element).split(':')[1] for element in row[1:]]
		except IndexError:
			raise Exception("Couldn't find any experiment information. Do you have a line in your schema that's (nearly) empty?")
		# check if condition names unique
		if not len(set(conditionsList)) == len(conditionsList):
			raise Exception("Same condition name used multiple times for same experiment. Please define each condition only once per experiment.")
		channelAliasesList = [extractAliases(element) for element in row[1:]]
		
		# store each condition and its channelNames in the dict
		for condition, channelNamesString, channelAliasesString in zip(conditionsList, channelNamesList, channelAliasesList):
			incompleteSchemaDict[experimentName][condition] = OrderedDict()
			channelNames = channelNamesString.split(',')
			if '' in channelNames:
				raise Exception("Channel names cannot be empty strings. Maybe you placed an extra comma somewhere?")
			incompleteSchemaDict[experimentName][condition]['channelNames'] = channelNames
			if channelAliasesString == '':  # no channelAliases provided: construct yourself
				channelAliases = [experimentName + '_' + condition + '_' + channelName for channelName in channelNames]
			else:  # channelAliases are provided.
				channelAliases = channelAliasesString.split(',')
			if '' in channelAliases:
				raise Exception("Channel aliases cannot be empty strings (if you define one alias, you should define all aliases in that condition). Maybe you placed an extra comma somewhere?")
			incompleteSchemaDict[experimentName][condition]['channelAliases'] = channelAliases
			numNames = len(channelNames)
			numAliases = len(channelAliases)
			numAllChannelNames += numNames
			numAllChannelAliases += numAliases
			experimentChannelNames += channelNames  # channel names must be tested for uniqueness per experiment
			experimentChannelAliases += channelAliases
			allChannelNames.update(channelNames)
			allChannelAliases.update(channelAliases)
			# check if there are as many names as aliases
			if not numNames == numAliases:
				raise Exception("Amount of channel names and channel aliases must either be equal, or no aliases should be provided.")
		
		incompleteSchemaDict[experimentName]['allExperimentConditions'] = conditionsList
		incompleteSchemaDict[experimentName]['allExperimentChannelNames'] = experimentChannelNames
		incompleteSchemaDict[experimentName]['allExperimentChannelAliases'] = experimentChannelAliases
		# check if channel names unique
		if not len(set(experimentChannelNames)) == len(experimentChannelNames):
			raise Exception("Same channel name used multiple times for same experiment. Please define each condition only once per experiment.")
	
	# check if channel aliases unique
	if not len(allChannelAliases) == numAllChannelAliases:
		raise Exception("Same alias name used multiple times for same experiment. Please define each condition only once per experiment.")
	
	# channelNames already checked per experiment, and are allowed to be non-unique across experiments since they get replaced anyway.
	incompleteSchemaDict['allConditions'] = list(allConditions)
	
	# check if some of the conditions provided do not actually have forbidden names (union not empty)
	forbiddenConditionNames = {'data', 'config', 'isotopicCorrection_matrix', 'wrapper', 'allExperimentConditions',
							   'protein', 'peptides', 'description', 'fold change log2(c1/c2)', 'p-value'}
	if allConditions & {'data', 'config', 'isotopicCorrection_matrix', 'wrapper'}:
		raise Exception("Please do not use as condition names any of the following: " + str(list(forbiddenConditionNames)))
	
	return incompleteSchemaDict


def fixFixableFormatMistakes(df):
	"""
	Takes a dataframe and fixes format mistakes frequently present in PD2.1 output.
	:param df:  pd.DataFrame    possibly containing a multitude of mistakes.
	:return df: pd.DataFrame    data without recognized format mistakes.
	"""
	columns = list(df.columns.values)
	# you've enabled "show flanking amino acids": DIRk --> [L].DIRk.[m]
	if df.sample(n=1)['Annotated Sequence'].item().count('.') == 2:  # sequence contains 2 dots
		df['Annotated Sequence'] = df['Annotated Sequence'].apply(lambda x: x.split('.')[1])  # select part between dots
	
	columns = list(df.columns.values)
	# you're using "Identifying Node" instead of "Identifying Node Type"
	if 'Identifying Node' in columns and 'Identifying Node Type' not in columns:
		# strip the everything after the last space (=node number between parentheses) + the last space.
		df['Identifying Node Type'] = df['Identifying Node'].apply(lambda s: s.rsplit(' ', 1)[0])
	
	columns = list(df.columns.values)
	# you're using "Isolation Interference in Percent" instead of 'Isolation Interference [%]'
	if 'Isolation Interference in Percent' in columns:
		# replace it
		newColumns = applyWrapper(columns, [('Isolation Interference in Percent', 'Isolation Interference [%]')])
		df.columns = newColumns
	
	return df


def applyWrapper(columns, wrapper):
	"""
	Takes a Nx2 tuple/list wrapper and transforms the columns in the dataframe df specified by the first column entries
	in wrapper into the new column name specified by the second colmumn entries in wrapper.
	:param columns:     pd.Index or list	old column names
	:param wrapper:     list(tuples)    	[(oldName, newName) for some columns]
	:return newColumns: list	        	transformed column names
	"""
	if isinstance(columns, pd.Index):
		newColumns = list(columns.values)
	elif isinstance(columns, list):
		newColumns = columns
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
	:return fullPath:	str		full path of where the data was saved
	:return fullPaths:	dict	full paths of where the data were saved, per data object
	"""
	# todo this function should call sub-functions per data type and optional arguments
	assert os.path.exists(path_out)
	
	extension = delim2ext(delim_out)
	fullPath = path_out + '/' + filename + extension
	
	if dataType == 'txt':
		np.savetxt(fullPath, data, delimiter=delim_out)
		return fullPath
	elif dataType == 'obj':
		pickle.dump(data, open(fullPath, 'wb'))
	elif dataType == 'df':
		if isinstance(data, dict):  # there are actually multiple dataFrames
			if inOneFile:  # save all removedData in one file.
				removedData = pd.DataFrame().append(list(data.values()))
				# data['missing'].columns contains all possible columns
				removedData.to_csv(fullPath, sep=delim_out, index=False,
								   columns=data['missing'].columns)
				return fullPath
			else:  # save all removedData in separate files per category.
				fullPaths = dict()
				for frameName, frame in data.items():
					assert isinstance(frame, pd.DataFrame)
					fullPaths[frameName] = path_out + '/' + filename + '_' + frameName + extension
					frame.to_csv(fullPaths[frameName], sep=delim_out, index=False)
				return fullPaths
		else:
			assert isinstance(data, pd.DataFrame)
			data.to_csv(fullPath, sep=delim_out, index=False)
			return fullPath
	elif dataType == 'fig':
		outFileFullPathNoExt = os.path.join(path_out, filename)
		with open(outFileFullPathNoExt + '.pkl', "wb") as fout:
			# save as pickle. Load again later as: ax = pickle.load(file('myplot.pickle')); then plt.show()
			pickle.dump(data, fout, protocol=4)
		# from matplotlib import pyplot as plt
		# ax = data
		data.savefig(outFileFullPathNoExt + '.png', format='png', bbox_inches='tight')
		return outFileFullPathNoExt + '.png'
	elif dataType == 'html':
		fullPath += '.html'
		with open(fullPath, 'w') as htmlFile:
			htmlFile.writelines(data)
		return fullPath


def delim2ext(delim):
	"""
	Returns the filename extension belonging to the specified delimiter, or returns empty string.
	"""
	if delim == ',':
		return '.csv'
	elif delim == '\t':
		return '.tsv'
	else:
		return ''


def ext2delim(ext):
	"""
	Returns the delimiter belonging to the specified filename extension, or returns empty string.
	"""
	ext.lstrip('.')
	if ext == 'csv':
		return ','
	elif ext == 'tsv':
		return '\t'
	else:
		return None


def genZip(outFilePath, inFilePaths):
	"""
	Creates a zipfile at specified ouFilePath location from a list of inFilePaths
	:param outFilePath:	str		output zip file path
	:param inFilePaths:	[ str ]	input file paths
	"""
	with zipfile.ZipFile(outFilePath, 'w') as zf:
		for f in inFilePaths:
			zf.write(f, arcname=os.path.basename(f))