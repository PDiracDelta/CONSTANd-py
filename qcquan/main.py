# -*- coding: utf-8 -*-

"""
Python implementation of mass spectrometer protein data analysis using the CONSTANd_RAS algorithm.
"""

import sys
import datetime
import traceback
from qcquan.dataIO import *


def main(jobConfigFilePath, doProcessing, doAnalysis, doReport, writeToDisk):
	"""
	Contains and explicits the workflow of the program. Using the booleans doProcessing, doAnalysis, deReport and
	writeToDisk one can control	which parts of the workflow to perform.
	"""
	# todo proper docu
	from qcquan.getConfig import getProcessingConfig, getJobConfig
	from qcquan.processingFlow import processDf
	from qcquan.analysisFlow import analyzeProcessingResult
	from qcquan.reportFlow import generateReport
	from time import time
	from qcquan_web.web import DB_setJobReportRelPaths
	logFilePath = os.path.abspath(os.path.join(jobConfigFilePath, os.path.join(os.pardir, 'log.txt')))
	if os.path.exists(logFilePath):
		# clean up old log file so as not to clog it while debugging
		os.remove(logFilePath)
	logging.basicConfig(filename=logFilePath, level=logging.INFO)
	startTime = time()
	jobParams = getJobConfig(jobConfigFilePath)  # config filenames + params for the combination of MSRuns
	allProcessingParams = {}  # specific params for each MSRun
	dfs = {}
	processingResults = {}
	MSRunNames = jobParams['schema']['allMSRuns']
	
	for eName in MSRunNames:
		""" Data processing """
		# get all input parameters
		allProcessingParams[eName] = getProcessingConfig(jobParams['schema'][eName]['config'])
		# get the dataframes
		# todo move this step to processingFlow --> NO because everything inside the Flow.py files should reside in memory, not on disk.
		dfs[eName] = importMSRunData(allProcessingParams[eName]['data'],
										  delim=allProcessingParams[eName]['delim_in'],
										  header=allProcessingParams[eName]['header_in'],
										  wrapper=allProcessingParams[eName]['wrapper'])
		processing_path_out = allProcessingParams[eName]['path_out']
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
				"Starting processing of MSRun '" + eName + "' of job '" + jobParams['jobName'] + "' at " +
				str(datetime.datetime.utcnow()).split('.')[0])
			processingResults[eName] = processDf(dfs[eName], allProcessingParams[eName], writeToDisk)
			logging.info(
				"Finished processing of MSRun '" + eName + "' of job '" + jobParams['jobName'] + "' at " +
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
				"No processing step performed nor processing file loaded for MSRun " + str(eName) + "!")
	
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
		logging.info("Starting visualization and report generation of job: " + jobParams['jobName'] + " at " +
					 str(datetime.datetime.utcnow()).split('.')[0])
		generateReport(analysisResults, jobParams, logFilePath, writeToDisk, allProcessingParams, startTime)
		DB_setJobReportRelPaths(jobID=jobDirName, resultpath=jobParams['path_results'],
								jobName=jobParams['jobName'])
		logging.info("Finished visualization and report generation of job: " + jobParams['jobName'] + " at " +
					 str(datetime.datetime.utcnow()).split('.')[0])
	else:
		logging.warning("No report generated!")
	stop = time()
	print(stop - startTime)


if __name__ == '__main__':  # this should not execute if main.py is not the main module called by the python interpreter,
	# so if you start main.py from within web.py or something, this won't be executed
	from qcquan_web.web import DB_setJobCompleted, DB_setJobFailed
	from qcquan_web import app
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
	else:  # you didn't call main.py from the command line but from pycharm
		doProcessing = True
		doAnalysis = True
		doReport = True
		writeToDisk = True
		
		from qcquan_web.config import ALLJOBSDIR

		# jobConfigFilePath = '/home/pdiracdelta/Documents/UHasselt/QCQuan/jobs/2018-02-05 14:22:56.936690_Gatto_newPSMs_selfprocessed_ref129/jobConfig_Gatto_newPSMs_selfprocessed_ref129.ini'
		jobConfigFilePath = '/home/pdiracdelta/Documents/UHasselt/QCQuan/jobs/2018-01-25 15:12:07.373181_ESmit_conditions/jobConfig_ESmit_conditions.ini'
		
		# jobConfigFilePath = os.path.join(ALLJOBSDIR, jobConfigFilePath)
	
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
