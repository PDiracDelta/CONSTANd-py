#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that simulates the web interface.
"""

import os, datetime
from dataIO import parseSchemaFile, unnest
from web.web import TMT2ICM, newJobDir
from json import dumps
from shutil import copyfile
from web.web import updateConfigs


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
	if True: # collapse hard coded
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
		elif exptype == 'COON_norm':
			coondatapath = '../data/COON data/PSMs/'
			datatype = '_g_ISO_norm'
			coonconfig = '../jobs/coonProcessingConfig.ini'
			coonwrapper = None  # '../jobs/coonWrapper.tsv'
			coonICM = None
			HC_JOBNAME = 'COON_norm'
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
		elif exptype == 'COON_SN_norm':
			coondatapath = '../data/COON data/PSMs/'
			datatype = '_h_ISO_SN_norm'
			coonconfig = '../jobs/coonProcessingConfig.ini'
			coonwrapper = None  # '../jobs/coonWrapper.tsv'
			coonICM = None
			HC_JOBNAME = 'COON_SN_norm'
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
		elif exptype == 'COON_noISO':
			coondatapath = '../data/COON data/PSMs/'
			datatype = '_a'
			coonconfig = '../jobs/coonProcessingConfig_iso.ini'
			coonwrapper = None  # '../jobs/coonWrapper.tsv'
			coonICM = '../jobs/coonICM.tsv' # TEST
			HC_JOBNAME = 'COON_noISO'
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
		elif exptype == 'COON_nonormnoconstand':
			coondatapath = '../data/COON data/PSMs/'
			datatype = '_e_ISO'
			coonconfig = '../jobs/coonProcessingConfig.ini'
			coonwrapper = None #'../jobs/coonWrapper.tsv'
			coonICM = None
			HC_JOBNAME = 'COON_nonormnoconstand'
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
		elif exptype == 'COON_SN_nonormnoconstand':
			coondatapath = '../data/COON data/PSMs/'
			datatype = '_f_ISO_SN'
			coonconfig = '../jobs/coonProcessingConfig.ini'
			coonwrapper = None #'../jobs/coonWrapper.tsv'
			coonICM = None
			HC_JOBNAME = 'COON_SN_nonormnoconstand'
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

	def getJobName():
		return HC_JOBNAME

	def uploadSchema(this_job_path):
		this_schemaPath = HC_SCHEMA
		destination = os.path.join(this_job_path, os.path.basename(this_schemaPath))
		copyfile(this_schemaPath, destination)
		return os.path.abspath(destination)

	def uploadFile(this_job_path, sourceDataPath, prefix):
		if sourceDataPath is None:
			return None
		else:
			destinationData = os.path.join(this_job_path, prefix + os.path.basename(sourceDataPath))
			copyfile(sourceDataPath, destinationData)
			return os.path.basename(destinationData)

	def transformICM(filePath, this_isTMTICM, this_channelAliasesPerCondition):
		from dataIO import exportData
		if this_isTMTICM:
			from dataIO import getTMTIsotopicDistributions
			icm = TMT2ICM(getTMTIsotopicDistributions(filePath), this_channelAliasesPerCondition)
		else:
			from dataIO import getIsotopicCorrectionsMatrix
			icm = getIsotopicCorrectionsMatrix(filePath)
			raise Exception(
				"not implemented (ICM column and row order transformation from non-TMT formatted input).")  # todo
		this_path = os.path.abspath(os.path.join(filePath, os.pardir))
		this_filename = os.path.basename(filePath)
		exportData(icm, dataType='txt', path_out=this_path, filename=this_filename[0:-4],
		           delim_out='\t')  # no extention

	def updateSchema(this_job_path, this_incompleteSchema):
		for eName in this_incompleteSchema:
			if eName == HC_ENAME1:
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA1,
				                                                  prefix=eName + '_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER1,
				                                                     prefix=eName + '_')
				this_incompleteSchema[eName]['config'] = os.path.join(this_job_path, uploadFile(this_job_path,
				                                                                                sourceDataPath=HC_CONFIG1,
				                                                                                prefix=eName + '_'))  # config needs FULL PATH!
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path,
				                                                                       sourceDataPath=HC_ICM1,
				                                                                       prefix=eName + '_')
			elif eName == HC_ENAME2:
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA2,
				                                                  prefix=eName + '_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER2,
				                                                     prefix=eName + '_')
				this_incompleteSchema[eName]['config'] = os.path.join(this_job_path, uploadFile(this_job_path,
				                                                                                sourceDataPath=HC_CONFIG2,
				                                                                                prefix=eName + '_'))  # config needs FULL PATH!
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path,
				                                                                       sourceDataPath=HC_ICM2,
				                                                                       prefix=eName + '_')
			elif eName == HC_ENAME3:
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA3,
				                                                  prefix=eName + '_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER3,
				                                                     prefix=eName + '_')
				this_incompleteSchema[eName]['config'] = os.path.join(this_job_path, uploadFile(this_job_path,
				                                                                                sourceDataPath=HC_CONFIG3,
				                                                                                prefix=eName + '_'))  # config needs FULL PATH!
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path,
				                                                                       sourceDataPath=HC_ICM3,
				                                                                       prefix=eName + '_')
			elif eName == HC_ENAME4:
				this_incompleteSchema[eName]['data'] = uploadFile(this_job_path, sourceDataPath=HC_DATA4,
				                                                  prefix=eName + '_')
				this_incompleteSchema[eName]['wrapper'] = uploadFile(this_job_path, sourceDataPath=HC_WRAPPER4,
				                                                     prefix=eName + '_')
				this_incompleteSchema[eName]['config'] = os.path.join(this_job_path, uploadFile(this_job_path,
				                                                                                sourceDataPath=HC_CONFIG4,
				                                                                                prefix=eName + '_'))  # config needs FULL PATH!
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = uploadFile(this_job_path,
				                                                                       sourceDataPath=HC_ICM4,
				                                                                       prefix=eName + '_')
			# in case no wrapper was uploaded
			if this_incompleteSchema[eName]['wrapper'] is None:
				wrapperFileName = eName + '_wrapper.tsv'
				open(os.path.join(this_job_path, wrapperFileName), 'w').close()
				this_incompleteSchema[eName]['wrapper'] = wrapperFileName

		isTMTICM = True
		for eName in incompleteSchema:
			ICMFile = this_incompleteSchema[eName]['isotopicCorrection_matrix']
			if ICMFile is not None:
				transformICM(os.path.join(this_job_path, ICMFile), isTMTICM,
				             this_incompleteSchema[eName]['channelNamesPerCondition'])

		return this_incompleteSchema

	def getBaseConfigFile():
		return 'baseProcessingConfig.ini'

	def updateWrappers(this_job_path, this_schema):
		# write the channel aliases to the wrapper
		for eName in this_schema:
			experiment = this_schema[eName]
			channelNames = unnest(experiment['channelNamesPerCondition'])
			channelAliases = unnest(experiment['channelAliasesPerCondition'])
			wrapperFile = os.path.join(this_job_path, experiment['wrapper'])
			# open user config parameters
			with open(wrapperFile, 'a') as fout:
				fout.write('\n')  # so you dont accidentally append to the last line
				for n, a in zip(channelNames, channelAliases):
					fout.write(n + '\t' + a + '\n')

	def getMasterConfig(this_job_path, this_job_name):
		this_masterConfigFile = uploadFile(this_job_path, sourceDataPath=HC_MASTERCONFIG,
		                                   prefix='')
		new_masterConfigFile = os.path.join(this_job_path, this_job_name + jobConfigNameSuffix)
		os.rename(os.path.join(this_job_path, this_masterConfigFile), new_masterConfigFile)
		return new_masterConfigFile

	def updateMasterConfig(this_job_path, this_masterConfigFileAbsPath, this_schema, this_jobname):
		# allChannelAliases = unnest([unnest(experiments['channelAliasesPerCondition']) for experiments in this_schema.values()])
		with open(this_masterConfigFileAbsPath, 'a') as fout:
			fout.write('\n')  # so you dont accidentally append to the last line
			fout.write('jobname = ' + this_jobname + '\n')
			fout.write('path_out = output_analysis\n')
			fout.write('path_results = results\n')
			fout.write('PCA_components = 2\n')
			fout.write('schema = ' + dumps(this_schema) + '\n')
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
	# also make sure the ICM columns are in the right order and the ICM is transformed from TMT table if necessary.
	schema = updateSchema(job_path, incompleteSchema)

	### STEP 3: update config files and wrapper files
	updateConfigs(job_path, schema)
	updateWrappers(job_path, schema)

	### STEP 4: get masterConfig from web and update it (add schema, date, path_out, path_results)
	masterConfigFileAbsPath = getMasterConfig(job_path, job_name)
	updateMasterConfig(job_path, masterConfigFileAbsPath, schema, job_name)

	return masterConfigFileAbsPath
