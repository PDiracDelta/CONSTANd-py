#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that simulates the web interface.
"""

import os, datetime
from dataIO import parseSchemaFile, unnest
from json import dumps


def webFlow():
	# HARDCODED FILE LOCATIONS
	HC_SCHEMA = '../jobs/schema6.tsv'
	HC_DATA1 = '../jobs/MB_noapostrophes.tsv'
	HC_DATA2 = '../jobs/MB_noapostrophes_bis.tsv'
	HC_WRAPPER1 = None
	HC_WRAPPER2 = '../jobs/wrapper6_bis.tsv'
	HC_CONFIG1 = '../jobs/config.ini'
	HC_CONFIG2 = '../jobs/config_bis.ini'
	HC_ICM1 = '../jobs/ICM6_default.tsv'
	HC_ICM2 = '../jobs/ICM6_default.tsv'

	from shutil import copyfile

	def newJobDir():
		jobPath = os.path.join('../jobs', str(datetime.datetime.now()))
		os.makedirs(jobPath)
		return os.path.abspath(jobPath)

	def uploadSchema(this_job_path):
		this_schemaPath = HC_SCHEMA
		destination = os.path.join(this_job_path, os.path.basename(this_schemaPath))
		copyfile(this_schemaPath, destination)
		return os.path.abspath(destination)

	def uploadFile(this_job_path, sourceDataPath, prefix):
		if sourceDataPath is None:
			return None
		else:
			destinationData = os.path.join(this_job_path, prefix+os.path.basename(sourceDataPath))
			copyfile(sourceDataPath, destinationData)
			return os.path.basename(destinationData)

	def updateSchema(this_job_path, this_incompleteSchema):
		for eName in this_incompleteSchema:
			if eName == 'human':
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA1,
				                                                  prefix=eName+'_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER1,
				                                                     prefix=eName+'_')
				this_incompleteSchema[eName]['config'] = uploadFile(this_job_path, sourceDataPath=HC_CONFIG1,
				                                                    prefix=eName+'_')
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path, sourceDataPath=HC_ICM1,
				                                                 prefix=eName+'_')
			elif eName == 'mouse':
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA2,
				                                                  prefix=eName+'_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER2,
				                                                     prefix=eName+'_')
				this_incompleteSchema[eName]['config'] = uploadFile(this_job_path, sourceDataPath=HC_CONFIG2,
				                                                    prefix=eName+'_')
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path, sourceDataPath=HC_ICM2,
				                                                 prefix=eName+'_')
			# in case no wrapper was uploaded
			if this_incompleteSchema[eName]['wrapper'] is None:
				wrapperFileName = os.path.join(this_job_path, eName+'_wrapper.tsv')
				open(wrapperFileName, 'w').close()
				this_incompleteSchema[eName]['wrapper'] = wrapperFileName
		return this_incompleteSchema

	def getBaseConfigFile():
		return 'baseConfig.ini'

	def updateConfigs(this_job_path, this_schema):
		import fileinput
		for eName in this_schema:
			experiment = this_schema[eName]
			configFile = experiment['config']
			# open user config parameters
			with open(configFile, 'a') as fout, fileinput.input(getBaseConfigFile()) as fin:
				fout.write('\n') # so you dont accidentally append to the last line
				# write baseConfig parameters
				for line in fin:
					if '[DEFAULT]' not in line:
						fout.write(line)
				# write schema parameters
				fout.write('\n')  # so you dont accidentally append to the last line
				fout.write('data = ' + experiment['data'] + '\n')
				fout.write('wrapper = ' + experiment['wrapper'] + '\n')
					# caution! use the ALIASES, and NOT the original names (they are rewritten by the wrapper)
				fout.write('intensityColumnsPerCondition = ' + dumps(experiment['channelAliasesPerCondition']) + '\n')
				fout.write('isotopicCorrection_matrix = ' + experiment['isotopicCorrection_matrix'] + '\n')
				# write output parameters
				fout.write('path_out = '+str(os.path.join(this_job_path, eName+'_output_processing/'))+'\n')
				fout.write('filename_out = '+str(eName)+'\n')

	def updateWrappers(this_schema):
		# write the channel aliases to the wrapper
		for eName in this_schema:
			experiment = this_schema[eName]
			channelNames = unnest(experiment['intensityColumnsPerCondition'])
			channelAliases = unnest(experiment['channelAliasesPerCondition'])
			wrapperFile = experiment['wrapper']
			# open user config parameters
			with open(wrapperFile, 'a') as fout:
				fout.write('\n') # so you dont accidentally append to the last line
				for n, a in zip(channelNames, channelAliases):
					fout.write(n+'\t'+a+'\n')

	def getMasterConfig(this_job_path):
		this_masterConfigFile = uploadFile(this_job_path, sourceDataPath='../jobs/masterConfig.ini',
		                                             prefix='')
		return this_masterConfigFile

	def updateMasterConfig(this_job_path, this_masterConfigFile, this_schema):
		with open(this_masterConfigFile, 'a') as fout:
			fout.write('\n')  # so you dont accidentally append to the last line
			fout.write('path_out = ' + os.path.join(this_job_path, 'output_analysis') + '\n')
			fout.write('path_results = ' + os.path.join(this_job_path, 'results') + '\n')
			fout.write('schema = '+dumps(this_schema)+'\n')
			fout.write('date = ' + str(os.path.basename(this_job_path).split('.')[0]) + '\n')


	### STEP 1: get schema and create new job
	job_path = newJobDir()
	schemaPath = uploadSchema(job_path)
	try:
		incompleteSchema = parseSchemaFile(schemaPath)
	except Exception as e: # remove file and job dir if something went wrong
		os.remove(schemaPath)
		os.removedirs(job_path)
		raise e

	### STEP 2: upload data files, wrapper files, ICM files and config (files), while updating schema with their locations.
	schema = updateSchema(job_path, incompleteSchema)

	### STEP 3: update config files and wrapper files
	updateConfigs(job_path, schema)
	updateWrappers(schema)

	### STEP 4: get masterConfig from web and update it (add schema, date, path_out, path_results)
	masterConfigFile = getMasterConfig(job_path)
	updateMasterConfig(job_path, masterConfigFile, schema)

	return masterConfigFile