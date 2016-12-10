#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that simulates the web interface.
"""

import os, datetime
from dataIO import parseSchemaFile
from json import dumps


def webFlow():
	from shutil import copyfile
	def newJobDir():
		jobPath = os.path.join('../jobs', str(datetime.datetime.now()))
		os.makedirs(jobPath)
		return os.path.abspath(jobPath)
	def uploadSchema(this_job_path):
		this_schemaPath = '../jobs/schema6.tsv'
		destination = os.path.join(this_job_path, os.path.basename(this_schemaPath))
		copyfile(this_schemaPath, destination)
		return os.path.abspath(destination)
	def uploadFile(this_job_path, sourceDataPath, prefix):
		if sourceDataPath is None:
			return None
		else:
			destinationData = os.path.join(this_job_path, prefix+os.path.basename(sourceDataPath))
			copyfile(sourceDataPath, destinationData)
			return os.path.abspath(destinationData)
	def updateSchema(this_job_path, incompleteSchema):
		for eName in incompleteSchema:
			if eName == 'human':
				incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath='../jobs/MB_noapostrophes.tsv',
				                                             prefix=eName+'_')
				incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=None,
				                                                prefix=eName+'_')
				incompleteSchema[eName]['config'] = uploadFile(this_job_path, sourceDataPath='../jobs/config.ini',
				                                               prefix=eName+'_')
				incompleteSchema[eName]['icm'] = uploadFile(this_job_path, sourceDataPath='../jobs/ICM6_default.tsv',
				                                               prefix=eName+'_')
			elif eName == 'mouse':
				incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath='../jobs/MB_noapostrophes_bis.tsv',
				                                             prefix=eName+'_')
				incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath='../jobs/wrapper6_bis.tsv',
				                                                prefix=eName+'_')
				incompleteSchema[eName]['config'] = uploadFile(this_job_path, sourceDataPath='../jobs/config_bis.ini',
				                                               prefix=eName+'_')
				incompleteSchema[eName]['icm'] = uploadFile(this_job_path, sourceDataPath='../jobs/ICM6_default.tsv',
				                                            prefix=eName+'_')
		return incompleteSchema
	def getBaseConfigFile():
		return 'baseConfig.ini'
	def updateConfigs(this_job_path, schema):
		import fileinput
		for eName in schema:
			experiment = schema[eName]
			configFile = experiment['config']
			# open user config parameters
			with open(configFile, 'a') as fout, fileinput.input(getBaseConfigFile()) as fin:
				fout.write('\n') # so you dont accidentally append to the last line
				# write baseConfig parameters
				for line in fin:
					if '[DEFAULT]' not in line:
						fout.write(line)
				# write schema parameters
				for parameter, value in experiment.items():
					if parameter != 'config':
						fout.write(dumps(parameter, value))
				# write output parameters
				fout.write('path_out = '+dumps(os.path.join(this_job_path, 'output/'))+'\n')
				fout.write('filename_out = '+dumps(eName)+'\n')
	def getMasterConfig(this_job_path):
		this_masterConfigFile = uploadFile(this_job_path, sourceDataPath='../jobs/masterConfig.ini',
		                                             prefix='')
		return this_masterConfigFile
	def updateMasterConfig(this_job_path, this_masterConfigFile, schema):
		with open(this_masterConfigFile, 'a') as fout:
			fout.write('\n')  # so you dont accidentally append to the last line
			fout.write('schema = '+dumps(schema)+'\n')
			fout.write('date = ' + dumps(os.path.basename(this_job_path).split('.')[0]) + '\n')


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

	### STEP 3: update config files
	updateConfigs(job_path, schema)

	### STEP 4: get masterConfig from web and update it (add schema, date)
	masterConfigFile = getMasterConfig(job_path)
	updateMasterConfig(job_path, masterConfigFile, schema)

	return os.path.join(job_path, masterConfigFile)