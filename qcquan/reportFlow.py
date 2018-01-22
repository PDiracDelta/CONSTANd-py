#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Workflow of the processing part of QCQuan.
"""


from qcquan.report import *
from qcquan.dataIO import exportData, genZip


def generateReport(analysisResults, params, logFilePath, writeToDisk, processingParams, startTime):
	"""
	Calls all the necessary functions to make visualizations and generate an HTML and PDF report.
	Takes an analysisResults object and uses the PCA and HC to generate a 2D PC plot as well as a HC dendrogram. If
	nConditions==2 then volcano plots for minExpression (injective) and fullExpression (non-injective, if available) are
	also produced and for each expression type a list of top differentials (according to adjusted p-value) is constructed.
	They are poured into an HTML report together with the parameters and metadata, and then converted into a PDF report.
	Because of technical reasons to do with image representation, the PDF is generated from a slightly different HTML
	file than the "public" HTML file.
	Graphs and report files are written to disk if so specified by writeToDisk.
	:param analysisResults:	list	[minProteinDF, fullProteinDF, PCAResult, HCResult,
									allExperimentsIntensitiesPerCommonPeptide, metadata, extraOutputFilesFullPaths]
	:param params:			dict	job (global) parameters
	:param logFilePath:		str		path to the log file with information about each processingFlow and analysisFlow call
	:param writeToDisk:		bool	write visualizations and reports to disk (if not: just pass the return statement)
	:param processingParams:dict	experiment-specific processing parameters (see getConfig.py.)
	:param startTime:		float	UNIX epoch timestamp at which the reportFlow was started
	"""
	import logging  # for some reason it doesnt work if you put it at the top of the file
	minProteinDF = analysisResults[0]
	fullProteinDF = analysisResults[1]
	PCAResult = analysisResults[2]
	HCResult = analysisResults[3]
	metadata = analysisResults[5]
	extraOutputFilesFullPaths = analysisResults[6]  # processedDfFullPaths, minProteinDF_fullPath, fullProteinDF_fullPath
	
	otherConditions = getOtherConditions(params['schema'], params['referenceCondition'])

	def getExpressionResults(this_proteinDF, this_schema):
		"""
		Sorts the protein dataframe, calculates the set of proteins in the results, generates a volcano plot and selects the
		top differentials by calling functions from report.py
		:param this_proteinDF:						pd.DataFrame    unsorted DE analysis results on the protein level
		:return this_sortedProteinExpressionsDF:	dict			DEA output table sorted according to adjusted p-value,
																	including proteins without DE results, per condition
		:return this_topDifferentialsDFs:			dict			top X sorted (on adjusted p-value) proteins, per condition
		:return this_volcanoPlot:					dict			plt.figure volcano plot, per condition
		:return this_set:							set				all proteins represented in the results
		"""
		# { condition: sortedProteinExpressionsDF }
		this_sortedProteinExpressionsDFs = getSortedProteinExpressionsDFs(this_proteinDF, this_schema, params['referenceCondition'])
		this_set = set()
		this_topDifferentialsDFs = dict()  # { condition: topDifferentialsDF }
		this_volcanoPlots = dict()  # { condition: volcanoPlot }
		# get the Expression results for each condition separately
		for otherCondition in otherConditions:
			this_set.update(set(this_sortedProteinExpressionsDFs[otherCondition]['protein']))
			# get top X differentials
			this_topDifferentialsDFs[otherCondition] = getTopDifferentials(this_sortedProteinExpressionsDFs[otherCondition], params['numDifferentials'])
			# data visualization
			this_volcanoPlots[otherCondition] = getVolcanoPlot(this_proteinDF, otherCondition, params['alpha'], params['FCThreshold'],
															 params['labelVolcanoPlotAreas'],
															 topIndices=list(this_topDifferentialsDFs[otherCondition]['protein']))
			# add protein IDs that were observed at least once but got removed, for completeness in the output csv.
			this_sortedProteinExpressionsDFs[otherCondition] = addMissingObservedProteins(this_sortedProteinExpressionsDFs[otherCondition],
																						  set(metadata['allObservedProteins'].loc[:, 'protein']))
		return this_sortedProteinExpressionsDFs, this_topDifferentialsDFs, this_volcanoPlots, this_set
	
	allDEResultsFullPaths = []  # paths to later pass on for mail attachments
	# include non-redundant peptide information ( = processed PSM file)
	allDEResultsFullPaths.extend(extraOutputFilesFullPaths)
	# do MINIMAL expression
	if params['minExpression_bool']:
		minSortedProteinExpressionsDFs, minTopDifferentialsDFs, minVolcanoPlots, minSet = getExpressionResults(minProteinDF, params['schema'])
		# save results
		minDEResultsFullPaths = exportData(minSortedProteinExpressionsDFs, dataType='df', path_out=params['path_results'],
										  filename=params['jobName'] + '_minSortedDifferentials',
										  delim_out=params['delim_out'])
		minVolcanoFullPaths = dict((otherCondition, exportData(minVolcanoPlots[otherCondition], dataType='fig',
														   path_out=params['path_results'],
														   filename=params['jobName'] + '_minVolcanoPlot_' + otherCondition))
							   for otherCondition in otherConditions)
		allDEResultsFullPaths.extend(list(minDEResultsFullPaths.values()))  # no need to know which path is which
		
	else:
		minTopDifferentialsDFs = None
		minVolcanoFullPaths = None
	
	# do FULL expression
	if params['fullExpression_bool']:
		fullSortedProteinExpressionsDFs, fullTopDifferentialsDFs, fullVolcanoPlots, fullSet = getExpressionResults(fullProteinDF, params['schema'])
		# save results
		fullDEResultsFullPaths = exportData(fullSortedProteinExpressionsDFs, dataType='df',
										   path_out=params['path_results'],
										   filename=params['jobName'] + '_fullSortedDifferentials',
										   delim_out=params['delim_out'])
		fullVolcanoFullPaths = dict((otherCondition, exportData(fullVolcanoPlots[otherCondition], dataType='fig',
															path_out=params['path_results'],
															filename=params['jobName'] + '_fullVolcanoPlot_' + otherCondition))
								for otherCondition in otherConditions)
		allDEResultsFullPaths.extend(list(fullDEResultsFullPaths.values()))  # no need to know which path is which
	else:
		fullTopDifferentialsDFs = None
		fullVolcanoFullPaths = None

	# metadata
	if params['minExpression_bool'] and params['fullExpression_bool']:
		# list( [in min but not in full], [in full but not in min] )
		metadata['diffMinFullProteins'] = [list(minSet.difference(fullSet)), list(fullSet.difference(minSet))]
		# todo combine into one

	# else:
	# 	minTopDifferentialsDF = pd.DataFrame()
	# 	fullTopDifferentialsDF = pd.DataFrame()
	# 	minVolcanoFullPath = None
	# 	fullVolcanoFullPath = None

	PCAPlot = getPCAPlot(PCAResult, params['schema'])
	if writeToDisk:
		PCAPlotFullPath = exportData(PCAPlot, dataType='fig', path_out=params['path_results'],
				   filename=params['jobName'] + '_PCAPlot')
	HCDendrogram = getHCDendrogram(HCResult, params['schema'])
	if writeToDisk:
		HCDendrogramFullPath = exportData(HCDendrogram, dataType='fig', path_out=params['path_results'],
				   filename=params['jobName'] + '_HCDendrogram')
	
	try:
		ScoreVsDeltaMppmScatter = getScoreVsDeltaMppmScatter(metadata['relPSMScoreVsDeltaMppmPerExp'])
		if writeToDisk:
			ScoreVsDeltaMppmScatterFullPath = exportData(ScoreVsDeltaMppmScatter, dataType='fig',
														 path_out=params['path_results'],
														 filename=params['jobName'] + '_ScoreVsDeltaMppmScatter')
	except KeyError:
		logging.warning("No relPSMScoreVsDeltaMppmPerExp QC info available. Not making MS1 calibration QC plot.")
		ScoreVsDeltaMppmScatter = None
	
	try:
		MS1IntensityHist = getMS1IntensityHist(metadata['MS1Intensities_PSMs'], metadata['MS1Intensities_peptides'])
		if writeToDisk:
			MS1IntensityHistFullPath = exportData(MS1IntensityHist, dataType='fig', path_out=params['path_results'],
					   filename=params['jobName'] + '_MS1IntensityHist')
	except KeyError as e:
		logging.warning("QC entry '" + str(e.args[0]) + "' was not found. Not producing MS1IntensityHist.")
		MS1IntensityHistFullPath = None


	if writeToDisk:
		htmlReport, pdfhtmlreport = makeHTML(jobParams=params, allProcessingParams=processingParams,
											 otherConditions=otherConditions,
											 minTopDifferentialsDFs=minTopDifferentialsDFs,
											 fullTopDifferentialsDFs=fullTopDifferentialsDFs,
											 minVolcanoFullPaths=minVolcanoFullPaths,
											 fullVolcanoFullPaths=fullVolcanoFullPaths,
											 PCAPlotFullPath=PCAPlotFullPath, HCDendrogramFullPath=HCDendrogramFullPath,
											 ScoreVsDeltaMppmScatterFullPath=ScoreVsDeltaMppmScatterFullPath,
											 MS1IntensityHistFullPath=MS1IntensityHistFullPath,
											 metadata=metadata, logFilePath=logFilePath, startTime=startTime)
		htmlFullPath = exportData(htmlReport, dataType='html', path_out=params['path_results'],
				   filename=params['jobName'] + '_report')
		pdfhtmlFullPath = exportData(pdfhtmlreport, dataType='html', path_out=params['path_results'],
				   filename=params['jobName'] + '_Report')

		pdfFullPath = HTMLtoPDF(pdfhtmlFullPath)
		# todo possibly remove need for special pdfhtml if weasyprint fetches the HTML from the web server via URL instead
		
		# zip the result files together (except the report file)
		from os.path import join
		resultsZipFullPath = join(params['path_results'], 'results.zip')
		genZip(resultsZipFullPath, allDEResultsFullPaths)
		
		from qcquan_web.web import send_mail
		### SEND JOB COMPLETED MAIL ###
		mailSuccess = send_mail(recipient=params['mailRecipient'], mailBodyFile='reportMail',
				  jobName=params['jobName'], jobID=params['jobID'], attachments=[pdfFullPath])
		if mailSuccess is not None:  # something went wrong
			import logging
			logging.error(mailSuccess)
			print(mailSuccess)
