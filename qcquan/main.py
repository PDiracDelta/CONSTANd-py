#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Python implementation of mass spectrometer protein data analysis using the CONSTANd_RAS algorithm.
"""

import sys
import datetime
import traceback
from constandpp.dataIO import *


def main(jobConfigFilePath, doProcessing, doAnalysis, doReport, writeToDisk):
	"""
	Contains and explicits the workflow of the program. Using the booleans doProcessing, doAnalysis, deReport and
	writeToDisk one can control	which parts of the workflow to perform.
	"""
	# todo proper docu
	from constandpp.getConfig import getProcessingConfig, getJobConfig
	from constandpp.processingFlow import processDf
	from constandpp.analysisFlow import analyzeProcessingResult
	from constandpp.reportFlow import generateReport
	from time import time
	from constandpp_web.web import DB_setJobReportRelPaths
	logFilePath = os.path.abspath(os.path.join(jobConfigFilePath, os.path.join(os.pardir, 'log.txt')))
	logging.basicConfig(filename=logFilePath, level=logging.INFO)
	start = time()
	jobParams = getJobConfig(jobConfigFilePath)  # config filenames + params for the combination of experiments
	processingParams = {}  # specific params for each experiment
	dfs = {}
	processingResults = {}
	experimentNames = jobParams['schema']['allExperiments']
	
	for eName in experimentNames:
		""" Data processing """
		# get all input parameters
		processingParams[eName] = getProcessingConfig(jobParams['schema'][eName]['config'])
		# get the dataframes
		# todo move this step to processingFlow --> NO because everything inside the Flow.py files should reside in memory, not on disk.
		dfs[eName] = importExperimentData(processingParams[eName]['data'],
										  delim=processingParams[eName]['delim_in'],
										  header=processingParams[eName]['header_in'],
										  wrapper=processingParams[eName]['wrapper'])
		processing_path_out = processingParams[eName]['path_out']
		processingResultsDumpFilename = os.path.join(processing_path_out, 'processingResultsDump_' + str(eName))
		if doProcessing:
			print('Starting processing of ' + eName + '...')
			if writeToDisk:  # prepare the output directories
				if not os.path.exists(os.path.abspath(processing_path_out)):  # do not overwrite dir
					assert os.path.exists(
						os.path.abspath(os.path.join(processing_path_out, os.path.pardir)))  # parent dir must exist
					os.makedirs(processing_path_out)
				else:
					raise Exception("Output path " + os.path.abspath(processing_path_out) + " already exists! Aborting.")
			
			# process every input dataframe
			logging.info(
				"Starting processing of experiment '" + eName + "' of job '" + jobParams['jobName'] + "' at " +
				str(datetime.datetime.utcnow()).split('.')[0])
			processingResults[eName] = processDf(dfs[eName], processingParams[eName], writeToDisk)
			logging.info(
				"Finished processing of experiment '" + eName + "' of job '" + jobParams['jobName'] + "' at " +
				str(datetime.datetime.utcnow()).split('.')[0])
			pickle.dump(processingResults[eName], open(processingResultsDumpFilename, 'wb'))  # TEST
		elif doAnalysis:
			try:
				processingResults[eName] = pickle.load(open(processingResultsDumpFilename, 'rb'))
			except FileNotFoundError:
				raise FileNotFoundError(
					"There is no previously processed data in this path: " + processingResultsDumpFilename)
		else:
			logging.warning(
				"No processing step performed nor processing file loaded for experiment " + str(eName) + "!")
	
	""" Data analysis """
	analysis_path_out = jobParams['path_out']
	analysisResultsDumpFilename = os.path.join(analysis_path_out, 'analysisResultsDump')
	if doAnalysis:
		print('Starting analysis...')
		if writeToDisk:  # prepare the output directories
			if not os.path.exists(analysis_path_out):  # do not overwrite dir
				assert os.path.exists(
					os.path.abspath(os.path.join(analysis_path_out, os.path.pardir)))  # parent dir must exist
				os.makedirs(analysis_path_out)
		
		# perform analysis
		logging.info("Starting analysis of job: " + jobParams['jobName'] + " at " +
					 str(datetime.datetime.utcnow()).split('.')[0])
		analysisResults = analyzeProcessingResult(processingResults, jobParams, writeToDisk)
		logging.info("Finished analysis of job: " + jobParams['jobName'] + " at " +
					 str(datetime.datetime.utcnow()).split('.')[0])
		pickle.dump(analysisResults, open(analysisResultsDumpFilename, 'wb'))  # TEST
	elif doReport:
		try:
			analysisResults = pickle.load(open(analysisResultsDumpFilename, 'rb'))
		except FileNotFoundError:
			raise FileNotFoundError(
				"There is no previously analyzed data in this path: " + analysisResultsDumpFilename)
	else:
		logging.warning("No analysis step performed nor analysis file loaded!")
	
	""" Visualize and generate report """
	results_path_out = jobParams['path_results']
	if doReport:
		print('Starting visualization and report generation...')
		if writeToDisk:  # prepare the output directories
			if not os.path.exists(results_path_out):  # do not overwrite dir
				assert os.path.exists(
					os.path.abspath(os.path.join(results_path_out, os.path.pardir)))  # parent dir must exist
				os.makedirs(results_path_out)
		
		# visualize and make a report
		logging.info("Starting visualization end report generation of job: " + jobParams['jobName'] + "at " +
					 str(datetime.datetime.utcnow()).split('.')[0])
		generateReport(analysisResults, jobParams, logFilePath, writeToDisk, processingParams, start)
		DB_setJobReportRelPaths(jobID=jobDirName, resultpath=jobParams['path_results'],
								jobName=jobParams['jobName'])
		logging.info("Finished visualization end report generation of job: " + jobParams['jobName'] + "at " +
					 str(datetime.datetime.utcnow()).split('.')[0])
	else:
		logging.warning("No report generated!")
	stop = time()
	print(stop - start)


if __name__ == '__main__':  # this should not execute if main.py is not the main module called by the python interpreter,
	from constandpp_web.web import DB_setJobCompleted, DB_setJobFailed
	from constandpp_web import app
	from traceback import print_exc
	
	args = sys.argv
	print(str(args))  # TEST
	if len(args) != 1:
		assert len(args) == 6  # todo use argparser
		jobConfigFilePath = args[1]
		doProcessing = (args[2] == 'True')
		doAnalysis = (args[3] == 'True')
		doReport = (args[4] == 'True')
		writeToDisk = (args[5] == 'True')
	# so if you start main.py from within web.py or something, this won't be executed
	else:  # you didn't call main.py from the command line but from pycharm
		doProcessing = False
		doAnalysis = False
		doReport = True
		writeToDisk = True
		
		from constandpp_web.config import ALLJOBSDIR
		#jobConfigFilePath = '2017-08-23 10:20:21.345488_test_updatedRemoveBadConfidence/jobConfig_test_updatedRemoveBadConfidence.ini'
		#jobConfigFilePath = '2017-08-23 17:01:13.474192_df/jobConfig_df.ini'
		# jobConfigFilePath = '2017-08-30 16:56:24.378757_ESmit_digests_v2/jobConfig_ESmit_digests_v2.ini'
		# jobConfigFilePath = '2017-09-13 17:16:18.431969_Gatto/jobConfig_Gatto.ini'
		# jobConfigFilePath = '2017-09-14 12:02:02.772295_Schmidt/jobConfig_Schmidt.ini'
		jobConfigFilePath = '2017-09-21 09:53:22.528726_test_zipfileinjobifopage/jobConfig_test_zipfileinjobifopage.ini'
		
		jobConfigFilePath = os.path.join(ALLJOBSDIR, jobConfigFilePath)
	
	with app.app_context():
		jobDirName = os.path.basename(os.path.abspath(os.path.join(jobConfigFilePath, os.pardir)))
		try:
			main(jobConfigFilePath=jobConfigFilePath, doProcessing=doProcessing, doAnalysis=doAnalysis,
				 doReport=doReport, writeToDisk=writeToDisk)
			DB_setJobCompleted(jobDirName)
		except:
			DB_setJobFailed(jobDirName)
			logging.error("\n=====\n===== A FATAL ERROR OCCURRED; PLEASE FIND THE STACK TRACE BELOW =====\n=====\n")
			logFilePath = os.path.abspath(os.path.join(jobConfigFilePath, os.path.join(os.pardir, 'log.txt')))
			traceback.print_exc(file=open(logFilePath, "a"))
			print_exc()