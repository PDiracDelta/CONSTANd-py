#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Workflow of the processing part of CONSTANd++.
"""

from constandpp.report import *
from constandpp.dataIO import exportData


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
	:param analysisResults:	list	[minProteinDF, fullProteinDF, PCAResult, HCResult, allExperimentsIntensitiesPerCommonPeptide]
	:param params:			dict	job (global) parameters
	:param logFilePath:		str		path to the log file with information about each processingFlow and analysisFlow call
	:param writeToDisk:		bool	write visualizations and reports to disk (if not: just pass the return statement)
	:param processingParams:dict	experiment-specific processing parameters (see getInput.py.)
	:param startTime:		float	UNIX epoch timestamp at which the reportFlow was started
	"""
	minProteinDF = analysisResults[0]
	fullProteinDF = analysisResults[1]
	PCAResult = analysisResults[2]
	HCResult = analysisResults[3]
	allExperimentsIntensitiesPerCommonPeptide = np.asarray(analysisResults[4])
	metadata = analysisResults[5]

	### TEST
	# make MA plots for comparing conditions
	testInterexperimentalOnPeptideLevel = False
	proteinLevelMAPlots = False
	if testInterexperimentalOnPeptideLevel:
		from constandpp.tools import MAPlot
		org3data = allExperimentsIntensitiesPerCommonPeptide[:, 17]
		#MAPlot(org3data, allExperimentsIntensitiesPerCommonPeptide[:, 11],
		#       'org2 exp [1] vs org2 exp [1]')
		MAPlot(org3data, allExperimentsIntensitiesPerCommonPeptide[:, 18],
			   'org2 exp [2] vs org3 exp [2]')
		MAPlot(org3data, allExperimentsIntensitiesPerCommonPeptide[:, 26],
			   'org2 exp [2] vs org3 exp [3]')
	if proteinLevelMAPlots:
		from constandpp.tools import MAPlot
		MAPlot(minProteinDF.loc[:, '3_muscle'], minProteinDF.loc[:, '4_muscle'],
			   'org2 exp [2] vs org3 exp [2]')
		MAPlot(minProteinDF.loc[:, '3_muscle'], minProteinDF.loc[:, '4_cerebrum'],
			   'org2 exp [2] vs org3 exp [3]')

	nConditions = len(params['schema']['allConditions'])
	allDEResultsFullPaths = []  # paths to later pass on for mail attachments
	# ONLY PRODUCE VOLCANO AND DEA IF CONDITIONS == 2
	if nConditions == 2:
		# todo split off these steps in a separate function which you call for min and full analysis
		# generate sorted (on p-value) list of differentials
		if params['minExpression_bool']:
			minSortedProteinExpressionsDF = getSortedProteinExpressionsDF(minProteinDF)
			minSet = set(minSortedProteinExpressionsDF['protein'])
			# get top X differentials
			minTopDifferentialsDF = getTopDifferentials(minSortedProteinExpressionsDF, params['numDifferentials'])
			# data visualization
			minVolcanoPlot = getVolcanoPlot(minProteinDF, params['alpha'], params['FCThreshold'],
											params['labelVolcanoPlotAreas'], topIndices=minTopDifferentialsDF.index)
			# add protein IDs that were observed at least once but got removed, for completeness in the output csv.
			minSortedProteinExpressionsDF = addMissingObservedProteins(minSortedProteinExpressionsDF, metadata['allObservedProteins'].loc[:, 'protein'][0])
		else:  # todo in this case (and also for fullExpression_bool) just let the jinja template handle the None variable.
			# but don't make a fake on here and then pass it onto makeHTML() like is done now.
			minSortedProteinExpressionsDF = pd.DataFrame(columns=['protein', 'significant', 'description', 'fold change log2(c1/c2)', 'adjusted p-value'])
			minTopDifferentialsDF = pd.DataFrame(columns=minSortedProteinExpressionsDF.columns)
		
		if params['fullExpression_bool']:
			fullSortedProteinExpressionsDF = getSortedProteinExpressionsDF(fullProteinDF)
			fullSet = set(fullSortedProteinExpressionsDF['protein'])
			# get top X differentials
			fullTopDifferentialsDF = getTopDifferentials(fullSortedProteinExpressionsDF, params['numDifferentials'])
			# data visualization
			fullVolcanoPlot = getVolcanoPlot(fullProteinDF, params['alpha'], params['FCThreshold'],
											 params['labelVolcanoPlotAreas'], topIndices=fullTopDifferentialsDF.index)
			# add protein IDs that were observed at least once but got removed, for completeness in the output csv.
			fullSortedProteinExpressionsDF = addMissingObservedProteins(fullSortedProteinExpressionsDF, metadata['allObservedProteins'].loc[:, 'protein'][0])
		else:
			fullSortedProteinExpressionsDF = pd.DataFrame(columns=['protein', 'significant', 'description', 'fold change log2(c1/c2)', 'adjusted p-value'])
			fullTopDifferentialsDF = pd.DataFrame(columns=fullSortedProteinExpressionsDF.columns)
		
		if params['minExpression_bool'] and params['fullExpression_bool']:
			# list( [in min but not in full], [in full but not in min] )
			metadata['diffMinFullProteins'] = [list(minSet.difference(fullSet)), list(fullSet.difference(minSet))]
			# todo combine into one

		if writeToDisk:
			if params['minExpression_bool']:
				minDEResultsFullPath = exportData(minSortedProteinExpressionsDF, dataType='df', path_out=params['path_results'],
						   filename=params['jobName'] + '_minSortedDifferentials', delim_out=params['delim_out'])
				minVolcanoFullPath = exportData(minVolcanoPlot, dataType='fig', path_out=params['path_results'],
						   filename=params['jobName'] + '_minVolcanoPlot')
				allDEResultsFullPaths.append(minDEResultsFullPath)
			else:
				minVolcanoFullPath = None
			if params['fullExpression_bool']:
				fullDEResultsFullPath = exportData(fullSortedProteinExpressionsDF, dataType='df', path_out=params['path_results'],
						   filename=params['jobName'] + '_fullSortedDifferentials', delim_out=params['delim_out'])
				fullVolcanoFullPath = exportData(fullVolcanoPlot, dataType='fig', path_out=params['path_results'],
						   filename=params['jobName'] + '_fullVolcanoPlot')
				allDEResultsFullPaths.append(fullDEResultsFullPath)
			else:
				fullVolcanoFullPath = None
	else:
		minTopDifferentialsDF = pd.DataFrame()
		fullTopDifferentialsDF = pd.DataFrame()
		minVolcanoFullPath = None
		fullVolcanoFullPath = None

	PCAPlot = getPCAPlot(PCAResult, params['schema'])
	if writeToDisk:
		PCAPlotFullPath = exportData(PCAPlot, dataType='fig', path_out=params['path_results'],
				   filename=params['jobName'] + '_PCAPlot')
	HCDendrogram = getHCDendrogram(HCResult, params['schema'])
	if writeToDisk:
		HCDendrogramFullPath = exportData(HCDendrogram, dataType='fig', path_out=params['path_results'],
				   filename=params['jobName'] + '_HCDendrogram')

	if writeToDisk:
		htmlReport, pdfhtmlreport = makeHTML(jobParams=params, allProcessingParams=processingParams,
											 minTopDifferentialsDF=minTopDifferentialsDF,
											 fullTopDifferentialsDF=fullTopDifferentialsDF,
											 minVolcanoFullPath=minVolcanoFullPath,
											 fullVolcanoFullPath=fullVolcanoFullPath,
											 PCAPlotFullPath=PCAPlotFullPath, HCDendrogramFullPath=HCDendrogramFullPath,
											 metadata=metadata, logFilePath=logFilePath, startTime=startTime)
		htmlFullPath = exportData(htmlReport, dataType='html', path_out=params['path_results'],
				   filename=params['jobName'] + '_report')
		pdfhtmlFullPath = exportData(pdfhtmlreport, dataType='html', path_out=params['path_results'],
				   filename=params['jobName'] + '_Report')

		pdfFullPath = HTMLtoPDF(pdfhtmlFullPath)
		# todo possibly remove need for special pdfhtml if weasyprint fetches the HTML from the web server via URL instead

		from constandpp_web.web import send_mail
		### SEND JOB COMPLETED MAIL ###
		mailSuccess = send_mail(recipient=params['mailRecipient'], mailBodyFile='reportMail',
				  jobName=params['jobName'], jobID=params['jobID'], attachments=[pdfFullPath]+allDEResultsFullPaths)
		if mailSuccess is not None:  # something went wrong
			import logging
			logging.error(mailSuccess)
			print(mailSuccess)
