#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Python implementation of mass spectrometer protein data analysis using the CONSTANd_RAS algorithm.
"""

import sys, logging, datetime
from dataIO import *


def intraInterMAPlots():
	# calculate all intra- and inter-MA plots for the Max data set (Intensities only)
	# on the PEPTIDE LEVEL
	constandFileName = '../jobs/2016-12-22 13:07:48.663095_MAX_SN_noTAM/output_analysis/MAX_SN_noTAM_CommonPeptideIntensities.tsv'
	PDFileName = '../jobs/2016-12-24 12:09:28.341520_MAX_SN_abundances_noTAM/output_analysis/MAX_SN_abundances_noTAM_CommonPeptideIntensities.tsv'
	rawFileName = '../jobs/2016-12-22 12:56:54.697920_MAX_SN_nonormnoconstand_noTAM/output_analysis/MAX_SN_nonormnoconstand_noTAM_CommonPeptideIntensities.tsv'
	cdf = importDataFrame(constandFileName, delim='\t')
	pddf = importDataFrame(PDFileName, delim='\t')
	rdf = importDataFrame(rawFileName, delim='\t')
	
	compareDIFFERENTconditions = True
	compareIDENTICALconditions = True
	diffboxplots = False  # THIS IS USELESS THE DIMENSIONS DONT MATCH
	if compareDIFFERENTconditions:
		D_intra_cmeans, D_intra_pdmeans, D_intra_rmeans = [], [], []
		D_intra_cvars, D_intra_pdvars, D_intra_rvars = [], [], []
		experiments = [[3, 4], [5, 6], [1, 2]]
		# INTRA
		for e in experiments:  # for each experiment
			for i in e:  # for each sample of BM in this experiment
				for j in e:  # compare with the samples of PM in the same experiment
					cm, cv = MA(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					D_intra_cmeans.append(cm)
					D_intra_cvars.append(cv)
					pdm, pdv = MA(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'PM' + str(j)])[
							   2:4]  # [2:4] only mean and var
					D_intra_pdmeans.append(pdm)
					D_intra_pdvars.append(pdv)
					rm, rv = MA(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					D_intra_rmeans.append(rm)
					D_intra_rvars.append(rv)
		print("+++INTRA: \n"
			  "CONSTANd means: MAD: " + str(np.mean(np.abs(D_intra_cmeans))) + "; values: " + str(D_intra_cmeans) + "\n"
																													"PD2.1 means: MAD: " + str(
			np.mean(np.abs(D_intra_pdmeans))) + "; values: " + str(D_intra_pdmeans) + "\n"
																					  "raw means: MAD: " + str(
			np.mean(np.abs(D_intra_rmeans))) + "; values: " + str(D_intra_rmeans) + "\n"
																					"CONSTANd vars: average: " + str(
			np.mean(D_intra_cvars)) + "; values: " + str(D_intra_cvars) + "\n"
																		  "PD2.1 vars: average: " + str(
			np.mean(D_intra_pdvars)) + "; values: " + str(D_intra_pdvars) + "\n"
																			"raw vars: average: " + str(
			np.mean(D_intra_rvars)) + "; values: " + str(D_intra_rvars) + "\n")
		# MAPlot(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'PM' + str(j)])
		# MAPlot(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'PM' + str(j)])
		# MAPlot(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'PM' + str(j)])
		boxPlot([D_intra_rvars, D_intra_cvars, D_intra_pdvars], labels=['raw', 'CONSTANd', 'PD2.1'], ylab='variance')
		
		# INTER
		# experiments = [[3,4],[5,6],[1,2]]
		D_inter_cmeans, D_inter_pdmeans, D_inter_rmeans = [], [], []
		D_inter_cvars, D_inter_pdvars, D_inter_rvars = [], [], []
		for l in range(len(experiments)):
			notE = unnest(experiments[:l] + experiments[l + 1:])
			e = experiments[l]
			for i in e:
				for j in notE:
					cm, cv = MA(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					D_inter_cmeans.append(cm)
					D_inter_cvars.append(cv)
					pdm, pdv = MA(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'PM' + str(j)])[
							   2:4]  # [2:4] only mean and var
					D_inter_pdmeans.append(pdm)
					D_inter_pdvars.append(pdv)
					rm, rv = MA(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					D_inter_rmeans.append(rm)
					D_inter_rvars.append(rv)
		print("+++INTER: \n"
			  "CONSTANd means: MAD: " + str(np.mean(np.abs(D_inter_cmeans))) + "; values: " + str(D_inter_cmeans) + "\n"
																													"PD2.1 means: MAD: " + str(
			np.mean(np.abs(D_inter_pdmeans))) + "; values: " + str(D_inter_pdmeans) + "\n"
																					  "raw means: MAD: " + str(
			np.mean(np.abs(D_inter_rmeans))) + "; values: " + str(D_inter_rmeans) + "\n"
																					"CONSTANd vars: average: " + str(
			np.mean(D_inter_cvars)) + "; values: " + str(D_inter_cvars) + "\n"
																		  "PD2.1 vars: average: " + str(
			np.mean(D_inter_pdvars)) + "; values: " + str(D_inter_pdvars) + "\n"
																			"raw vars: average: " + str(
			np.mean(D_inter_rvars)) + "; values: " + str(D_inter_rvars) + "\n")
		# plot one graph for each (the last one):
		# MAPlot(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'PM' + str(j)])
		# MAPlot(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'PM' + str(j)])
		# MAPlot(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'PM' + str(j)])
		boxPlot([D_inter_rvars, D_inter_cvars, D_inter_pdvars], labels=['raw', 'CONSTANd', 'PD2.1'], ylab='variance')
	
	if compareIDENTICALconditions:
		I_intra_cmeans, I_intra_pdmeans, I_intra_rmeans = [], [], []
		I_intra_cvars, I_intra_pdvars, I_intra_rvars = [], [], []
		experiments = [[3, 4], [5, 6], [1, 2]]
		# INTRA
		for e in experiments:  # for each experiment
			for l in range(len(e)):  # for each sample in (and each condition) this experiment
				notI = e[l + 1:]
				i = e[l]
				if notI:  # not empty
					for j in notI:  # compare with the other sample(s) (that are of the same condition) in this experiment
						cm, cv = MA(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'BM' + str(j)])[
								 2:4]  # [2:4] only mean and var
						I_intra_cmeans.append(cm)
						I_intra_cvars.append(cv)
						cm, cv = MA(cdf.loc[:, 'PM' + str(i)], cdf.loc[:, 'PM' + str(j)])[
								 2:4]  # [2:4] only mean and var
						I_intra_cmeans.append(cm)
						I_intra_cvars.append(cv)
						
						pdm, pdv = MA(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'BM' + str(j)])[
								   2:4]  # [2:4] only mean and var
						I_intra_pdmeans.append(pdm)
						I_intra_pdvars.append(pdv)
						pdm, pdv = MA(pddf.loc[:, 'PM' + str(i)], pddf.loc[:, 'PM' + str(j)])[
								   2:4]  # [2:4] only mean and var
						I_intra_pdmeans.append(pdm)
						I_intra_pdvars.append(pdv)
						
						rm, rv = MA(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'BM' + str(j)])[
								 2:4]  # [2:4] only mean and var
						I_intra_rmeans.append(rm)
						I_intra_rvars.append(rv)
						rm, rv = MA(rdf.loc[:, 'PM' + str(i)], rdf.loc[:, 'PM' + str(j)])[
								 2:4]  # [2:4] only mean and var
						I_intra_rmeans.append(rm)
						I_intra_rvars.append(rv)
		print("+++INTRA: \n"
			  "CONSTANd means: MAD: " + str(np.mean(np.abs(I_intra_cmeans))) + "; values: " + str(I_intra_cmeans) + "\n"
																													"PD2.1 means: MAD: " + str(
			np.mean(np.abs(I_intra_pdmeans))) + "; values: " + str(I_intra_pdmeans) + "\n"
																					  "raw means: MAD: " + str(
			np.mean(np.abs(I_intra_rmeans))) + "; values: " + str(I_intra_rmeans) + "\n"
																					"CONSTANd vars: average: " + str(
			np.mean(I_intra_cvars)) + "; values: " + str(I_intra_cvars) + "\n"
																		  "PD2.1 vars: average: " + str(
			np.mean(I_intra_pdvars)) + "; values: " + str(I_intra_pdvars) + "\n"
																			"raw vars: average: " + str(
			np.mean(I_intra_rvars)) + "; values: " + str(I_intra_rvars) + "\n")
		# MAPlot(cdf.loc[:, 'BM' + str(1)], cdf.loc[:, 'BM' + str(2)])
		# MAPlot(pddf.loc[:, 'BM' + str(1)], pddf.loc[:, 'BM' + str(2)])
		# MAPlot(rdf.loc[:, 'BM' + str(1)], rdf.loc[:, 'BM' + str(2)])
		boxPlot([I_intra_rvars, I_intra_cvars, I_intra_pdvars], labels=['raw', 'CONSTANd', 'PD2.1'], ylab='variance')
		
		# INTER
		# experiments = [[3,4],[5,6],[1,2]]
		I_inter_cmeans, I_inter_pdmeans, I_inter_rmeans = [], [], []
		I_inter_cvars, I_inter_pdvars, I_inter_rvars = [], [], []
		for l in range(len(experiments)):
			notE = unnest(experiments[l + 1:])
			e = experiments[l]
			for i in e:
				for j in notE:
					cm, cv = MA(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'BM' + str(j)])[2:4]  # [2:4] only mean and var
					I_inter_cmeans.append(cm)
					I_inter_cvars.append(cv)
					cm, cv = MA(cdf.loc[:, 'PM' + str(i)], cdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					I_inter_cmeans.append(cm)
					I_inter_cvars.append(cv)
					
					pdm, pdv = MA(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'BM' + str(j)])[
							   2:4]  # [2:4] only mean and var
					I_inter_pdmeans.append(pdm)
					I_inter_pdvars.append(pdv)
					pdm, pdv = MA(pddf.loc[:, 'PM' + str(i)], pddf.loc[:, 'PM' + str(j)])[
							   2:4]  # [2:4] only mean and var
					I_inter_pdmeans.append(pdm)
					I_inter_pdvars.append(pdv)
					
					rm, rv = MA(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'BM' + str(j)])[2:4]  # [2:4] only mean and var
					I_inter_rmeans.append(rm)
					I_inter_rvars.append(rv)
					rm, rv = MA(rdf.loc[:, 'PM' + str(i)], rdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					I_inter_rmeans.append(rm)
					I_inter_rvars.append(rv)
		print("+++INTER: \n"
			  "CONSTANd means: MAD: " + str(np.mean(np.abs(I_inter_cmeans))) + "; values: " + str(I_inter_cmeans) + "\n"
																													"PD2.1 means: MAD: " + str(
			np.mean(np.abs(I_inter_pdmeans))) + "; values: " + str(I_inter_pdmeans) + "\n"
																					  "raw means: MAD: " + str(
			np.mean(np.abs(I_inter_rmeans))) + "; values: " + str(I_inter_rmeans) + "\n"
																					"CONSTANd vars: MAD: " + str(
			np.mean(I_inter_cvars)) + "; values: " + str(I_inter_cvars) + "\n"
																		  "PD2.1 vars: MAD: " + str(
			np.mean(I_inter_pdvars)) + "; values: " + str(I_inter_pdvars) + "\n"
																			"raw vars: MAD: " + str(
			np.mean(I_inter_rvars)) + "; values: " + str(I_inter_rvars) + "\n")
		# plot one graph for each (the last one):
		# MAPlot(cdf.loc[:, 'BM' + str(3)], cdf.loc[:, 'BM' + str(6)])
		# MAPlot(pddf.loc[:, 'BM' + str(3)], pddf.loc[:, 'BM' + str(6)])
		# MAPlot(rdf.loc[:, 'BM' + str(3)], rdf.loc[:, 'BM' + str(6)])
		boxPlot([I_inter_rvars, I_inter_cvars, I_inter_pdvars], labels=['raw', 'CONSTANd', 'PD2.1'], ylab='variance')
	
	if diffboxplots:
		# THIS IS USELESS THE DIMENSIONS DONT MATCH
		# transform all to arrays
		allLists = [I_inter_rvars, I_intra_rvars, I_inter_cvars, I_intra_cvars, I_inter_pdvars, I_intra_pdvars,
					D_inter_rvars, D_intra_rvars, D_inter_cvars, D_intra_cvars, D_inter_pdvars, D_intra_pdvars]
		I_inter_rvars, I_intra_rvars, I_inter_cvars, I_intra_cvars, I_inter_pdvars, I_intra_pdvars, D_inter_rvars, D_intra_rvars, D_inter_cvars, D_intra_cvars, D_inter_pdvars, D_intra_pdvars = (
		np.array(x) for x in allLists)
		# amke boxplots of differences
		boxPlot([I_inter_rvars - I_intra_rvars, I_inter_cvars - I_intra_cvars, I_inter_pdvars - I_intra_pdvars],
				labels=['raw', 'CONSTANd', 'PD2.1'], ylab=r'\Delta V')
		boxPlot([D_inter_rvars - D_intra_rvars, D_inter_cvars - D_intra_cvars, D_inter_pdvars - D_intra_pdvars],
				labels=['raw', 'CONSTANd', 'PD2.1'], ylab=r'\Delta V')


def compareDEAresults():
	df1file = '../jobs/2016-12-26 15:44:27.523732_MAX_noTAM/output_analysis/MAX_noTAM_results_minimal.tsv'
	df2file = '../jobs/2016-12-22 13:07:48.663095_MAX_SN_noTAM/output_analysis/MAX_SN_noTAM_results_minimal.tsv'
	# df1file = '../jobs/2016-12-30 14:46:02.689189_MAX_abundances_noTAM/output_analysis/MAX_abundances_noTAM_results_minimal.tsv'
	# df2file = '../jobs/2016-12-24 12:09:28.341520_MAX_SN_abundances_noTAM/output_analysis/MAX_SN_abundances_noTAM_results_minimal.tsv'
	df1sorted = importDataFrame(df1file, delim='\t').sort('protein')
	df2sorted = importDataFrame(df2file, delim='\t').sort('protein')
	assert set(df1sorted['protein']) == set(
		df2sorted['protein'])  # else you cant compare without first matching each protein
	# fcM, fcA, fcm, fcv =
	fmean, fmax = RDHPlot(np.power(2, df1sorted.loc[:, 'log2 fold change c1/c2']),
						  np.power(2, df2sorted.loc[:, 'log2 fold change c1/c2']), quantity='log2 fold change')
	print("fc mean: " + str(fmean) + "\nmax: " + str(fmax))
	# pM, pA, pm, pv =
	pmean, pmax = RDHPlot(df1sorted.loc[:, 'adjusted p-value'], df2sorted.loc[:, 'adjusted p-value'],
						  quantity='adjusted p-value')
	print("p mean: " + str(pmean) + "\nmax: " + str(pmax))


def dataSuitabilityMA():
	# rawfile = '../jobs/2016-12-20 14:39:09.476567_COON_SN_nonormnoconstand/output_analysis/COON_SN_nonormnoconstand_CommonPeptideIntensities.tsv'
	rawfile = '../jobs/2016-12-21 16:07:19.300450_MAX_SN_nonormnoconstand/output_analysis/MAX_SN_nonormnoconstand_CommonPeptideIntensities.tsv'
	rawSNvalues = importDataFrame(rawfile, delim='\t', header=None)
	if 'MAX' in rawfile:
		plexity = 6
	elif 'COON' in rawfile:
		plexity = 8
	# i = 3; j = 1
	# MAPlot(rawSNvalues.loc[:, i], rawSNvalues.loc[:, j], title='')
	# i=3; j=1
	# MAPlot(rawSNvalues.loc[:, i], rawSNvalues.loc[:, j], title='')
	# MAX
	# i=4; j=0
	i = 2;
	j = 6
	MAPlot(rawSNvalues.loc[:, i], rawSNvalues.loc[:, j], title='')
	for i in range(plexity):
		for j in range(i):
			MAPlot(rawSNvalues.loc[:, i], rawSNvalues.loc[:, j], title=str(i + 1) + " versus " + str(j + 1))


def devStuff(df, params):  # TEST
	# performanceTest()
	# isotopicCorrectionsTest()
	# MS2IntensityDoesntMatter(df)
	# testDataComplementarity(df)
	# compareIntensitySN(None, None, title='')
	# abundancesPCAHCD()
	# compareICmethods()
	compareAbundancesIntSN()
	# intraInterMAPlots()
	# compareDEAresults()
	# dataSuitabilityMA()
	pass


def main(jobConfigFilePath, doProcessing, doAnalysis, doReport, writeToDisk, testing):
	"""
	For now this is just stuff for debugging and testing. Later:
	Contains and explicits the workflow of the program. Using the booleans doProcessing, doAnalysis and writeToDisk one
	can control	which parts of the workflow to perform.
	"""
	# todo proper docu
	from getInput import getProcessingInput, getJobInput
	from processingFlow import processDf
	from analysisFlow import analyzeProcessingResult
	from reportFlow import generateReport
	from time import time
	from web.web import DB_setJobReportRelPaths
	
	logFilePath = os.path.abspath(os.path.join(jobConfigFilePath, os.path.join(os.pardir, 'log.txt')))
	logging.basicConfig(filename=logFilePath, level=logging.INFO)
	start = time()
	jobParams = getJobInput(jobConfigFilePath)  # config filenames + params for the combination of experiments
	processingParams = {}  # specific params for each experiment
	dfs = {}
	processingResults = {}
	experimentNames = list(jobParams['schema'].keys())
	
	if not testing:
		for eName in experimentNames:
			""" Data processing """
			# get all input parameters
			processingParams[eName] = getProcessingInput(jobParams['schema'][eName]['config'])
			# get the dataframes
			dfs[eName] = importExperimentData(processingParams[eName]['data'], delim=processingParams[eName]['delim_in'],
											  header=processingParams[eName]['header_in'],
											  wrapper=processingParams[eName]['wrapper'])
			processing_path_out = processingParams[eName]['path_out']
			processingResultsDumpFilename = os.path.join(processing_path_out, 'processingResultsDump_' + str(eName))
			if doProcessing:
				print('Starting processing of ' + eName + '...')
				# prepare the output directories
				if not os.path.exists(os.path.abspath(processing_path_out)):  # do not overwrite dir
					assert os.path.exists(
						os.path.abspath(os.path.join(processing_path_out, os.path.pardir)))  # parent dir must exist
					os.makedirs(processing_path_out)
				else:
					if jobConfigFilePath != 'job/jobConfig.ini':  # TEST
						raise Exception(
							"Output path " + os.path.abspath(processing_path_out) + " already exists! Aborting.")
				
				# process every input dataframe
				logging.info(
					"Starting processing of experiment '" + eName + "' of job '" + jobParams['jobname'] + "' at " +
					str(datetime.datetime.now()).split('.')[0])
				processingResults[eName] = processDf(dfs[eName], processingParams[eName], writeToDisk)
				logging.info(
					"Finished processing of experiment '" + eName + "' of job '" + jobParams['jobname'] + "' at " +
					str(datetime.datetime.now()).split('.')[0])
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
			# prepare the output directories
			if not os.path.exists(analysis_path_out):  # do not overwrite dir
				assert os.path.exists(
					os.path.abspath(os.path.join(analysis_path_out, os.path.pardir)))  # parent dir must exist
				os.makedirs(analysis_path_out)
			
			# perform analysis
			logging.info("Starting analysis of job: " + jobParams['jobname'] + " at " +
						 str(datetime.datetime.now()).split('.')[0])
			analysisResults = analyzeProcessingResult(processingResults, jobParams, writeToDisk)
			logging.info("Finished analysis of job: " + jobParams['jobname'] + " at " +
						 str(datetime.datetime.now()).split('.')[0])
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
			# prepare the output directories
			if not os.path.exists(results_path_out):  # do not overwrite dir
				assert os.path.exists(
					os.path.abspath(os.path.join(results_path_out, os.path.pardir)))  # parent dir must exist
				os.makedirs(results_path_out)
			
			# visualize and make a report
			logging.info("Starting visualization end report generation of job: " + jobParams['jobname'] + "at " +
						 str(datetime.datetime.now()).split('.')[0])
			generateReport(analysisResults, jobParams, logFilePath, writeToDisk, processingParams, start)
			DB_setJobReportRelPaths(jobID=jobDirName, resultpath=jobParams['path_results'],
									jobName=jobParams['jobname'])
			logging.info("Finished visualization end report generation of job: " + jobParams['jobname'] + "at " +
						 str(datetime.datetime.now()).split('.')[0])
		else:
			logging.warning("No report generated!")
	
	elif testing:
		devStuff(dfs, processingParams)
	stop = time()
	print(stop - start)


if __name__ == '__main__':  # this should not execute if main.py is not the main module called by the python interpreter,
	from web.web import DB_setJobCompleted, DB_setJobFailed
	from web import app
	from traceback import print_exc
	
	args = sys.argv
	print(str(args))  # TEST
	if len(args) != 1:
		assert len(args) == 7  # todo use argparser
		jobConfigFilePath = args[1]
		doProcessing = (args[2] == 'True')
		doAnalysis = (args[3] == 'True')
		doReport = (args[4] == 'True')
		writeToDisk = (args[5] == 'True')
		testing = (args[6] == 'True')
	# so if you start main.py from within web.py or something, this won't be executed
	else:
		doProcessing = False
		doAnalysis = False
		doReport = True
		writeToDisk = True
		testing = False
		
		# jobConfigFilePath = 'job/jobConfig.ini' # TEST
		jobConfigFilePath = '../jobs/2017-01-23 02:40:26.553814_two/jobConfig_two.ini'
	# jobConfigFilePath = webFlow(exptype='COON')
	# jobConfigFilePath = webFlow(exptype='COON', previousjobdirName='2016-12-26 10:56:10.646919_COON')
	# jobConfigFilePath = webFlow(exptype='COON_SN')
	# jobConfigFilePath = webFlow(exptype='COON_SN', previousjobdirName='2016-12-20 14:21:47.786288_COON_SN')
	# jobConfigFilePath = webFlow(exptype='COON_norm') # todo constand uitzetten
	# jobConfigFilePath = webFlow(exptype='COON_norm', previousjobdirName='2016-12-12 22:43:38.030716_COON_norm')
	# jobConfigFilePath = webFlow(exptype='COON_SN_norm')  # todo constand uitzetten
	# jobConfigFilePath = webFlow(exptype='COON_SN_norm', previousjobdirName='2016-12-12 22:48:30.701250_COON_SN_norm')
	# jobConfigFilePath = webFlow(exptype='COON_nonormnoconstand')  # todo constand uitzetten
	# jobConfigFilePath = webFlow(exptype='COON_nonormnoconstand', previousjobdirName='2016-12-20 14:31:47.045927_COON_nonormnoconstand')
	# jobConfigFilePath = webFlow(exptype='COON_noISO')
	# jobConfigFilePath = webFlow(exptype='COON_noISO', previousjobdirName='2016-12-16 16:38:30.536344_COON_noISO')
	# jobConfigFilePath = webFlow(exptype='COON_SN_nonormnoconstand')  # todo constand uitzetten
	# jobConfigFilePath = webFlow(exptype='COON_SN_nonormnoconstand', previousjobdirName='2016-12-20 14:39:09.476567_COON_SN_nonormnoconstand')
		jobConfigFilePath = '/home/pdiracdelta/Documents/UHasselt/CONSTANd++/jobs/2017-04-14 17:46:37.527494_alleswerkt?/jobConfig_alleswerkt?.ini'
	
	with app.app_context():
		jobDirName = os.path.basename(os.path.abspath(os.path.join(jobConfigFilePath, os.pardir)))
		try:
			main(jobConfigFilePath=jobConfigFilePath, doProcessing=doProcessing, doAnalysis=doAnalysis,
				 doReport=doReport,
				 testing=testing, writeToDisk=writeToDisk)
			DB_setJobCompleted(jobDirName)
		except:
			DB_setJobFailed(jobDirName)
			print_exc()
