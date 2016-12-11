#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that simulates the web interface.
"""

import os, datetime
from dataIO import parseSchemaFile, unnest
from json import dumps


def webFlow(exptype='dummy', previousjobdirName=None):
	jobsPath = '../jobs'
	jobConfigNameSuffix = '_jobConfig.ini'
	if previousjobdirName is not None: # work on previously processed/analyzed data
		previousjobdir = os.path.join(jobsPath, previousjobdirName)
		if os.path.exists(previousjobdir):
			jobName = '_'.join(os.path.basename(previousjobdir).split('_')[1:]) # haha this looks funky
			assumedjobconfigfilepath = os.path.join(previousjobdir, jobName+jobConfigNameSuffix)
			if os.path.exists(assumedjobconfigfilepath):
				return assumedjobconfigfilepath
			else:
				raise Exception("No previous job config file found. Aborting.")

	# HARDCODED FILE LOCATIONS
	if exptype == 'dummy':
		HC_JOBNAME = 'dummy'
		HC_SCHEMA = '../jobs/schema6.tsv'
		HC_ENAME1 = 'human'
		HC_ENAME2 = 'mouse'
		HC_DATA1 = '../jobs/MB_noapostrophes.tsv'
		HC_DATA2 = '../jobs/MB_noapostrophes_bis.tsv'
		HC_WRAPPER1 = None
		HC_WRAPPER2 = '../jobs/wrapper6_bis.tsv'
		HC_CONFIG1 = '../jobs/processingConfig.ini'
		HC_CONFIG2 = '../jobs/processingConfig_bis.ini'
		HC_ICM1 = '../jobs/ICM6_default.tsv'
		HC_ICM2 = '../jobs/ICM6_default.tsv'
		HC_MASTERCONFIG = '../jobs/jobConfig.ini'
	elif exptype == 'COON':
		coondatapath = '../data/COON data/PSMs/'
		datatype = '_e_ISO'
		coonconfig = '../jobs/coonProcessingConfig.ini'
		coonwrapper = None #'../jobs/coonWrapper.tsv'
		coonICM = None
		HC_JOBNAME = 'COON'
		HC_SCHEMA = '../jobs/coonSchema.tsv'
		HC_ENAME1 = 'BR1'
		HC_ENAME2 = 'BR2'
		HC_ENAME3 = 'BR3'
		HC_ENAME4 = 'BR4'
		HC_DATA1 = coondatapath+HC_ENAME1+datatype+'.txt'
		HC_DATA2 = coondatapath+HC_ENAME2+datatype+'.txt'
		HC_DATA3 = coondatapath+HC_ENAME3+datatype+'.txt'
		HC_DATA4 = coondatapath+HC_ENAME4+datatype+'.txt'
		HC_WRAPPER1 = coonwrapper
		HC_WRAPPER2 = coonwrapper
		HC_WRAPPER3 = coonwrapper
		HC_WRAPPER4 = coonwrapper
		HC_CONFIG1 = coonconfig
		HC_CONFIG2 = coonconfig
		HC_CONFIG3 = coonconfig
		HC_CONFIG4 = coonconfig
		HC_ICM1 = coonICM
		HC_ICM2 = coonICM
		HC_ICM3 = coonICM
		HC_ICM4 = coonICM
		HC_MASTERCONFIG = '../jobs/coonJobConfig.ini'
	elif exptype == 'COON_SN':
		coondatapath = '../data/COON data/PSMs/'
		datatype = '_f_ISO_SN'
		coonconfig = '../jobs/coonProcessingConfig.ini'
		coonwrapper = None  # '../jobs/coonWrapper.tsv'
		coonICM = None
		HC_JOBNAME = 'COON_SN'
		HC_SCHEMA = '../jobs/coonSchema.tsv'
		HC_ENAME1 = 'BR1'
		HC_ENAME2 = 'BR2'
		HC_ENAME3 = 'BR3'
		HC_ENAME4 = 'BR4'
		HC_DATA1 = coondatapath + HC_ENAME1 + datatype + '.txt'
		HC_DATA2 = coondatapath + HC_ENAME2 + datatype + '.txt'
		HC_DATA3 = coondatapath + HC_ENAME3 + datatype + '.txt'
		HC_DATA4 = coondatapath + HC_ENAME4 + datatype + '.txt'
		HC_WRAPPER1 = coonwrapper
		HC_WRAPPER2 = coonwrapper
		HC_WRAPPER3 = coonwrapper
		HC_WRAPPER4 = coonwrapper
		HC_CONFIG1 = coonconfig
		HC_CONFIG2 = coonconfig
		HC_CONFIG3 = coonconfig
		HC_CONFIG4 = coonconfig
		HC_ICM1 = coonICM
		HC_ICM2 = coonICM
		HC_ICM3 = coonICM
		HC_ICM4 = coonICM
		HC_MASTERCONFIG = '../jobs/coonJobConfig.ini'

	from shutil import copyfile

	def getJobName():
		return HC_JOBNAME

	def newJobDir(this_job_name):
		jobPath = os.path.join('../jobs', str(datetime.datetime.now())+'_'+this_job_name)
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
			if eName == HC_ENAME1:
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA1,
				                                                  prefix=eName+'_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER1,
				                                                     prefix=eName+'_')
				this_incompleteSchema[eName]['config'] = os.path.join(this_job_path, uploadFile(this_job_path, sourceDataPath=HC_CONFIG1,
				                                                    prefix=eName+'_')) # config needs FULL PATH!
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path, sourceDataPath=HC_ICM1,
				                                                 prefix=eName+'_')
			elif eName == HC_ENAME2:
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA2,
				                                                  prefix=eName+'_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER2,
				                                                     prefix=eName+'_')
				this_incompleteSchema[eName]['config'] = os.path.join(this_job_path, uploadFile(this_job_path, sourceDataPath=HC_CONFIG2,
				                                                    prefix=eName+'_')) # config needs FULL PATH!
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path, sourceDataPath=HC_ICM2,
				                                                 prefix=eName+'_')
			elif eName == HC_ENAME3:
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA3,
				                                                  prefix=eName+'_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER3,
				                                                     prefix=eName+'_')
				this_incompleteSchema[eName]['config'] = os.path.join(this_job_path, uploadFile(this_job_path, sourceDataPath=HC_CONFIG3,
				                                                    prefix=eName+'_')) # config needs FULL PATH!
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path, sourceDataPath=HC_ICM3,
				                                                 prefix=eName+'_')
			elif eName == HC_ENAME4:
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA4,
				                                                  prefix=eName+'_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER4,
				                                                     prefix=eName+'_')
				this_incompleteSchema[eName]['config'] = os.path.join(this_job_path, uploadFile(this_job_path, sourceDataPath=HC_CONFIG4,
				                                                    prefix=eName+'_')) # config needs FULL PATH!
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path, sourceDataPath=HC_ICM4,
				                                                 prefix=eName+'_')
			# in case no wrapper was uploaded
			if this_incompleteSchema[eName]['wrapper'] is None:
				wrapperFileName = eName+'_wrapper.tsv'
				open(os.path.join(this_job_path, wrapperFileName), 'w').close()
				this_incompleteSchema[eName]['wrapper'] = wrapperFileName
		return this_incompleteSchema

	def getBaseConfigFile():
		return 'baseProcessingConfig.ini'

	def updateConfigs(this_job_path, this_schema):
		import fileinput
		for eName in this_schema:
			experiment = this_schema[eName]
			configFile = os.path.join(this_job_path, experiment['config'])
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
				# fout.write('channelNamesPerCondition = ' + dumps(experiment['channelAliasesPerCondition']) + '\n')
				fout.write('intensityColumns = ' + dumps(unnest(experiment['channelAliasesPerCondition'])) + '\n')
				if experiment['isotopicCorrection_matrix'] is not None:
					fout.write('isotopicCorrection_matrix = ' + experiment['isotopicCorrection_matrix'] + '\n')
				else:
					fout.write('isotopicCorrection_matrix\n')
				# write output parameters
				fout.write('path_out = '+eName+'_output_processing/\n')
				fout.write('filename_out = '+str(eName)+'\n')

	def updateWrappers(this_job_path, this_schema):
		# write the channel aliases to the wrapper
		for eName in this_schema:
			experiment = this_schema[eName]
			channelNames = unnest(experiment['channelNamesPerCondition'])
			channelAliases = unnest(experiment['channelAliasesPerCondition'])
			wrapperFile = os.path.join(this_job_path, experiment['wrapper'])
			# open user config parameters
			with open(wrapperFile, 'a') as fout:
				fout.write('\n') # so you dont accidentally append to the last line
				for n, a in zip(channelNames, channelAliases):
					fout.write(n+'\t'+a+'\n')

	def getMasterConfig(this_job_path, this_job_name):
		this_masterConfigFile = uploadFile(this_job_path, sourceDataPath=HC_MASTERCONFIG,
		                                             prefix='')
		new_masterConfigFile = os.path.join(this_job_path, this_job_name+jobConfigNameSuffix)
		os.rename(os.path.join(this_job_path, this_masterConfigFile), new_masterConfigFile)
		return new_masterConfigFile

	def updateMasterConfig(this_job_path, this_masterConfigFileAbsPath, this_schema, this_jobname):
		with open(this_masterConfigFileAbsPath, 'a') as fout:
			fout.write('\n')  # so you dont accidentally append to the last line
			fout.write('jobname = ' + this_jobname+ '\n')
			fout.write('path_out = output_analysis\n')
			fout.write('path_results = results\n')
			fout.write('schema = '+dumps(this_schema)+'\n')
			fout.write('date = ' + str(os.path.basename(this_job_path).split('.')[0]) + '\n')


	### STEP 1: get schema and create new job
	job_name = getJobName()
	job_path = newJobDir(job_name)
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
	updateWrappers(job_path, schema)

	### STEP 4: get masterConfig from web and update it (add schema, date, path_out, path_results)
	masterConfigFileAbsPath = getMasterConfig(job_path, job_name)
	updateMasterConfig(job_path, masterConfigFileAbsPath, schema, job_name)

	return masterConfigFileAbsPath