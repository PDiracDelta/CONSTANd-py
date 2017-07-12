#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Functions involved in generating the report that includes:
* volcano plot
* PCA plot
* HC tree
* list of differentially expressed proteins, sorted on adjusted p-value
* metadata info (warnings etc.)
* overview of parameters used in the workflow
"""
import pandas as pd
import numpy as np
from constandpp.tools import unnest, getOtherConditions
from warnings import warn
from scipy.cluster.hierarchy import dendrogram
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import markers
from matplotlib.colors import to_hex
from constandpp import fontweight, fontsize, figwidth, figheight

# from adjustText import adjust_text

# adjust font size globally
matplotlib.rcParams.update({'font.size': fontsize, 'font.weight': fontweight})
matplotlib.use('GTK3Agg')


# save matplotlib images without whitespace: savefig('foo.png', bbox_inches='tight')


def distinguishableColours(n, type='jet'):
	"""
	Generates distinguishable colours in the format of a n array.
	:param n: 			int         number of colours
	:param type:        str         specifier of which colourmap to use
	:return:            np.ndarray  nx4 array
	"""
	cmap = plt.get_cmap(type)  # brg, jet
	return cmap(np.linspace(0, 1.0, n))


def getColours(schema, hex=False):
	"""
	Returns list of colours for all the channels in all experiments (based on schema) so that the channels of the same
	condition have the same colour.
	:param schema:  			dict    schema of the experiments' hierarchy
	:return channelColoursDict:	dict    colour for each channel; a different one for each condition
	"""
	numConditions = len(schema['allConditions'])
	distColours = distinguishableColours(numConditions)
	if hex:
		distColours = [to_hex(x) for x in distColours]
	channelColoursDict = dict()
	c = 0  # colour counter
	for cond in schema['allConditions']:
		for eName in schema['allExperiments']:
			if cond in schema[eName]:
				numChannels = len(schema[eName][cond]['channelAliases'])
				channelColoursDict.update(dict(zip(schema[eName][cond]['channelAliases'], np.tile(distColours[c], (numChannels, 1)))))
		c += 1
	assert c == numConditions
	return channelColoursDict


def distinguishableMarkers(n):
	"""
	Returns a list of n (easily) distinguishable markers, or just visible markers if n is too large. If n is extremely
	large, starts repeating markers.
	:param n:   int     number of distinguishable markers
	:return:    list    (easily) (distinguishable) markers
	"""
	easilyDistinguishable = ['o', 's', 'x', '*', 'v', 'd', '+', '^']
	allMarkers = markers.MarkerStyle.markers
	visibleMarkers = allMarkers.copy()
	for k, v in list(allMarkers.items()):
		if v == 'nothing':
			del visibleMarkers[k]
	
	if n > len(easilyDistinguishable):  # not enough distinguishable markers
		if n < len(visibleMarkers):  # enough visible markers
			warn("More experiments than easily distinguishable markers; using all (visible) markers.")
			return [list(allMarkers.keys())[i] for i in range(n)]
		else:  # not enough markers at all
			warn("More experiments than markers. Using all (visible) markers with possible repetitions!")
			# number of times to re-use ALL visible markers
			nRepetitions = np.mod(len(visibleMarkers), n)
			nResidual = len(visibleMarkers) - n
			return list(visibleMarkers.keys()) * nRepetitions + list(visibleMarkers.keys())[0:nResidual]
	else:  # have enough distinguishable markers for the n experiments
		return easilyDistinguishable[0:n]


def getMarkers(schema):
	"""
	Returns list of markers for the channels in all experiments (based on schema) so that the channels of the same
	experiment have the same marker.
	:param schema:  			dict    schema of the experiments' hierarchy
	:return channelMarkersDict:	dict	marker for each channel; a different one for each experiment
	"""
	distMarkers = distinguishableMarkers(len(schema))
	channelMarkersDict = {}
	i = 0
	for eName in schema['allExperiments']:
		for alias in schema[eName]['allExperimentChannelAliases']:
			channelMarkersDict[alias] = distMarkers[i]
		i += 1
	return channelMarkersDict


def getSortedProteinExpressionsDFs(proteinDF, schema, referenceCondition):
	"""
	Returns a list of sortedProteinDFs (see getSortedProteinExpressionsDF); one for each non-reference condition.
	:param proteinDF: 					pd.DataFrame    	unsorted DE analysis results on the protein level
	:param schema: 						dict				schema of the experiments' hierarchy.
	:return sortedProteinExpressionsDFs:dict				{ condition : sortedProteinDF }
	"""
	otherConditions = getOtherConditions(schema, referenceCondition)
	sortedProteinExpressionsDFs = dict()#zip(otherConditions, [None, ]*len(otherConditions)))
	for condition in otherConditions:
		sortedProteinExpressionsDFs[condition] = getSortedProteinExpressionsDF(proteinDF, referenceCondition, condition)
	return sortedProteinExpressionsDFs


def getSortedProteinExpressionsDF(proteinDF, referenceCondition, condition):
	"""
	Sorts the differential protein data according to adjusted p-value (high p-value should also mean high DE) and resets
	the index. Returns only the columns	specified.
	Later in the workflow, the head() function will called to get the top X differentials from this df.
	:param proteinDF:			pd.DataFrame    unsorted DE analysis results on the protein level
	:param referenceCondition:	str				name of the condition to be used as the reference
	:param condition:			str				name of the non-reference condition for which to get the sorted DE
	:return:    				pd.DataFrame    sorted according to adjusted p-value and only specified columns
	"""
	reportColumns = ['protein', 'description']
	
	adjustedPValueColumn = 'adjusted p-value (' + condition + ')'
	FCColumn = 'log2 fold change ('+condition+')'
	numPeptColumn = '#peptides ('+condition+')'
	refNumPeptColumn = '#peptides ('+referenceCondition+')'
	reportColumns.extend([FCColumn, adjustedPValueColumn, numPeptColumn, refNumPeptColumn])

	return proteinDF.loc[:, reportColumns].sort_values(by=adjustedPValueColumn, ascending=True)


def addMissingObservedProteins(sortedProteinExpressionsDF, allProteinsSet):
	"""
	Add all proteins in allProteinsSet that are not present -- due to missing values or whatever -- in the proteins
	column of the DEA output, for completeness. The other columns for these entries are NaN.
	:param sortedProteinExpressionsDF:	pd.DataFrame	DEA output table with only useful entries
	:param allProteinsSet:				Set				all proteins observed in at least 1 PSM of at least 1 experiment
	:return sortedProteinExpressionsDF:	pd.DataFrame	DEA output table including proteins without DE results
	"""
	presentProteinsSet = set(sortedProteinExpressionsDF.loc[:, 'protein'])
	missingProteinsList = list(allProteinsSet.difference(presentProteinsSet))
	sortedProteinExpressionsDF.append(pd.Series({'protein': missingProteinsList}), ignore_index=True)
	return sortedProteinExpressionsDF


def getTopDifferentials(sortedDifferentialsDF, numDifferentials):
	"""
	Takes a sorted protein differentials dataframe and returns the top `numDifferentials` entries according to the order.
	:param sortedDifferentialsDF:	pd.DataFrame	all sorted proteins (normally according to adjusted p-value)
	:param numDifferentials:		int				number of top DE proteins to select
	:return:						pd.DataFrame	top X sorted proteins (normally according to adjusted p-value)
	"""
	if numDifferentials < len(sortedDifferentialsDF):
		topDifferentials = sortedDifferentialsDF.head(
			numDifferentials)  # minSortedDifferentialProteinsDF.loc[range(numDifferentials), :]
	else:
		topDifferentials = sortedDifferentialsDF  # todo log this
	return topDifferentials


def getVolcanoPlot(df, alpha, FCThreshold, labelPlot=[False, ] * 4, topIndices=None):
	"""
	Generates a volcano plot using the log2 fold changes and adjusted p-values of the differential protein data. It is
	divided in several regions according to the indicated significance level (as can be determined by alpha and FCThreshold)
	and colours each region differently. Optionally, it only shows the labels of the proteins specified by topIndices.
	The IDs of the protein data points each region may be annotated if so specified by labelPlot.
	:param df:				pd.DataFrame	DE analysis results data on the protein level.
	:param alpha:			float			significance level used in the t-test of the DE analysis
	:param FCThreshold:		float			log2 fold change threshold used in the DE analysis
	:param labelPlot:		[ bool ]		label all proteins in the different significance regions: [ yes, p, fc, no ]
	:param topIndices:		list			indices of proteins for which to show the label exclusively
	:return volcanoPlot:	plt.figure		volcano plot as a matplotlib figure object
	"""
	# todo add protein ID labels according to sorted list entry ID
	volcanoPlot = plt.figure(
		figsize=(figwidth, figheight))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	# maximize figure
	# mng = plt.get_current_fig_manager()
	# mng.full_screen_toggle()
	plt.title(r'Volcano Plot ($FC>$' + str(FCThreshold) + r'; $\alpha=$' + str(alpha) + ')', figure=volcanoPlot)
	plt.xlabel(r'log$_2$(fold change)', figure=volcanoPlot)
	plt.ylabel(r'-log$_{10}$(p-value) ', figure=volcanoPlot)
	# get indices of different levels of significance
	significantIndices_yes = df[df['significant'] == 'yes'].index
	significantIndices_p = df[df['significant'] == 'p'].index
	significantIndices_fc = df[df['significant'] == 'fc'].index
	significantIndices_no = df[df['significant'] == 'no'].index
	
	# produce scatterPlot for each category of significance
	# YES
	xdataYES = df.loc[significantIndices_yes, 'fold change log2(c1/c2)']
	ydataYES = -np.log10(df.loc[significantIndices_yes, 'adjusted p-value'])
	labelsYES = df.loc[significantIndices_yes, 'protein']
	plt.scatter(xdataYES, ydataYES, color='r', figure=volcanoPlot)
	# P
	xdataP = df.loc[significantIndices_p, 'fold change log2(c1/c2)']
	ydataP = -np.log10(df.loc[significantIndices_p, 'adjusted p-value'])
	labelsP = df.loc[significantIndices_p, 'protein']
	plt.scatter(xdataP, ydataP, color='b', figure=volcanoPlot)
	# FC
	xdataFC = df.loc[significantIndices_fc, 'fold change log2(c1/c2)']
	ydataFC = -np.log10(df.loc[significantIndices_fc, 'adjusted p-value'])
	labelsFC = df.loc[significantIndices_fc, 'protein']
	plt.scatter(xdataFC, ydataFC, color='g', figure=volcanoPlot)
	# NO
	xdataNO = df.loc[significantIndices_no, 'fold change log2(c1/c2)']
	ydataNO = -np.log10(df.loc[significantIndices_no, 'adjusted p-value'])
	labelsNO = df.loc[significantIndices_no, 'protein']
	plt.scatter(xdataNO, ydataNO, color='k', figure=volcanoPlot)
	
	if topIndices is not None:  # remove labels of non-top X differentially expressed proteins
		labelsYES.loc[~labelsYES.index.isin(topIndices)] = ''
		labelsP.loc[~labelsP.index.isin(topIndices)] = ''
		labelsFC.loc[~labelsFC.index.isin(topIndices)] = ''
		labelsNO.loc[~labelsNO.index.isin(topIndices)] = ''
	
	# annotate where requested
	for labelPlotBool, xdata, ydata, labels in zip(labelPlot, [xdataYES, xdataP, xdataFC, xdataNO],
												   [ydataYES, ydataP, ydataFC, ydataNO],
												   [labelsYES, labelsP, labelsFC, labelsNO]):
		if labelPlotBool:
			for x, y, label in zip(xdata, ydata, labels):
				plt.annotate(label, xy=(x, y), xytext=(-1, 1), textcoords='offset points', ha='right', va='bottom',
							 fontsize=20)
				
	# adjust limits
	allXdata = pd.concat([xdataFC, xdataYES, xdataNO, xdataP])
	allYdata = pd.concat([ydataFC, ydataYES, ydataNO, ydataP])
	try:
		Xabsmax = max(abs(np.floor(min(allXdata))), abs(np.ceil(max(allXdata))))
		plt.xlim([-Xabsmax, Xabsmax])
	except ValueError:
		plt.xlim([-5, 5])
	try:
		plt.ylim([0, np.ceil(max(allYdata) / 5) * 5])  # int(base * round(float(x)/base))
	except ValueError:
		plt.ylim([0, 100])  # int(base * round(float(x)/base))
	
	# plt.show() # TEST
	return volcanoPlot


def getPCAPlot(PCAResult, schema, title=None):
	"""
	Generates a 2D plot of each quantification channel's first 2 PC scores. Identical colour means identical condition,
	and identical marker means identical experiment.
	:param PCAResult:	np.ndarray		PC scores of the channels for each protein (see getPCA() in analysis.py)
	:param schema:		dict			schema of the experiments' hierarchy
	:param title:		str				title for the plot
	:return PCAPlot:	plt.figure		PCA plot as a matplotlib figure object
	"""
	# todo also show percentage variance explained per PC
	PCAPlot = plt.figure(
		figsize=(figwidth, figheight))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	# maximize figure
	if title is None:
		title = 'Principal Component scores'
	plt.title(title, figure=PCAPlot)
	plt.xlabel('First PC', figure=PCAPlot)
	plt.ylabel('Second PC', figure=PCAPlot)
	
	# labels for annotation
	allChannelAliases = unnest([schema[eName]['allExperimentChannelAliases'] for eName in schema['allExperiments']])
	# generate colors/markers so that the channels of the same condition/experiment have the same colour/markers
	channelColorsDict = getColours(schema)
	channelMarkersDict = getMarkers(schema)
	
	for (x, y, label) in zip(PCAResult[:, 0], PCAResult[:, 1], allChannelAliases):
		# produce scatterPlot of two first principal components and annotate
		plt.scatter(x, y, color=channelColorsDict[label], marker=channelMarkersDict[label], figure=PCAPlot, s=80)
		plt.annotate(label, xy=(x, y), xytext=(-1, 1),
					 textcoords='offset points', ha='right', va='bottom', fontsize=20)
	legendHandles = []
	legendStrings = []
	# look up the corresponding experiment name for each marker to construct the legend
	markersToCheck = set(channelMarkersDict.values())
	for channel, marker in channelMarkersDict.items():
		if marker in markersToCheck:
			for eName in schema['allExperiments']:
				if channel in schema[eName]['allExperimentChannelAliases']:
					handle = plt.scatter([], [], color='k', marker=marker, s=160)
					legendHandles.append(handle)
					legendStrings.append(eName)
					markersToCheck.remove(marker)
					break
	plt.legend(legendHandles, legendStrings, scatterpoints=1)  # , loc=2)
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
	# plt.show() # TEST
	return PCAPlot


def getHCDendrogram(HCResult, schema, title=None):
	"""
	Generates a hierarchical clustering dendrogram (horizontal) using the NxN linkage matrix HCResult. Each leaf
	corresponds to a quantification channel. Identical colour means identical condition, and experiment labels are shown.
	:param HCResult:		np.ndarray		NxN linkage matrix (N=#channels)
	:param schema:			dict			schema of the experiments' hierarchy
	:return HCDendrogram:	plt.figure		HC dendrogram as a matplotlib figure object
	"""
	# hierarchical clustering dendrogram
	allChannelAliases = unnest([schema[eName]['allExperimentChannelAliases'] for eName in schema['allExperiments']])
	HCDendrogram = plt.figure(
		figsize=(figwidth, figheight))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	# maximize figure
	# mng = plt.get_current_fig_manager()
	# mng.full_screen_toggle()
	if title is None:
		title = 'Hierarchical Clustering Dendrogram'
	plt.title(title, figure=HCDendrogram)
	plt.xlabel('distance', figure=HCDendrogram)
	plt.ylabel('reporter channel', figure=HCDendrogram)
	# generate colors/markers so that the channels of the same condition/experiment have the same colour/markers
	# first transform to hex code because the dendrogram() function only takes strings.
	channelColorsDict = getColours(schema, hex=False)
	dendrogram(HCResult, orientation='right', leaf_rotation=0., leaf_font_size=24, labels=allChannelAliases,
			   link_color_func=lambda x: channelColorsDict[allChannelAliases[x]] if x < len(allChannelAliases) else 'k',
			   above_threshold_color='k')
	ax = plt.gca()
	ylbls = ax.get_ymajorticklabels()
	for i in range(len(ylbls)):
		ylbls[i].set_color(channelColorsDict[ylbls[i].get_text()])
	# ylbls[i].set_color(colorsPerCondition[i])
	# plt.show()  # TEST
	return HCDendrogram


def makeHTML(jobParams, allProcessingParams, otherConditions, minTopDifferentialsDFs, fullTopDifferentialsDFs, minVolcanoFullPaths,
			 fullVolcanoFullPaths, PCAPlotFullPath, HCDendrogramFullPath, metadata, logFilePath, startTime):
	"""
	Pour all report visualizations, the list(s) of differentials, metadata and job parameters into an HTML file.
	A second HTML file used for conversion to PDF is generated slightly different from the one used	for actual HTML
	representation, for technical reasons to do with image representation.
	:param jobParams:				dict			job (global) parameters
	:param allProcessingParams:		dict			per experiment, all parameters for the processing step
	:param minTopDifferentialsDFs:	{pd.DataFrame}	top X differential protein data (minimal expression, injective)
													sorted on adjusted p-value and only specified columns, per condition
	:param fullTopDifferentialsDFs:	{pd.DataFrame}	top X differential protein data (full expression, non-injective)
													sorted on adjusted p-value and only specified columns, per condition
	:param minVolcanoFullPaths:		{str}			path to the volcano plot image (minimal expression), per condition
	:param fullVolcanoFullPaths:	{str}			path to the volcano plot image (full expression), per condition
	:param PCAPlotFullPath: 		str				path to the PCA plot image
	:param HCDendrogramFullPath: 	str				path to the HC dendrogram image
	:param metadata:				dict			[noIsotopicCorrection, RTIsolationInfo, noMasterProteinAccession,
													minSingleConditionProteins, fullSingleConditionProteins,
													uncommonPeptides, commonNanValues]
	:param logFilePath:				str				path to the log file with information about each
													processingFlow and analysisFlow call
	:param startTime:				float			UNIX epoch timestamp at which the reportFlow was started
	:return htmlReport:				str				report HTML (for actual HTML representation)
	:return pdfhtmlreport:			str				report HTML (for conversion to PDF)
	"""
	from flask import render_template
	from time import time
	from constandpp_web import app
	from os import path, pardir
	from pandas import set_option
	import logging
	
	allJobsParDir = path.abspath(path.join(app.config.get('ALLJOBSDIR'), pardir))
	set_option('display.max_colwidth', -1)  # otherwise the Description column text is truncated.
	
	def injectColumnWidthHTML(DETableHTML):
		"""
		Takes the HTML string that generates the differential proteins table and inserts a <colgroup> element before the
		<thead> element by splitting, concatenating and rejoining. The colgroup element defines the relative column widths.
		:param DETableHTML:	str		HTML table with the differential proteins
		:return:			str		HTML table with the differential proteins and an extra colgroup element.
		"""
		#columnWidthHTML = '<colgroup><col width="8%" /><col width="52%" /><col width="16%" /><col width="12%" /><col width="12%" /></colgroup>'
		#columnWidthHTML = '<colgroup><col style="min-width:10em;" /><col class="block" /><col style="min-width:10em;" /><col style="min-width:12em;" /><col style="min-width:8em;" /></colgroup>'
		columnWidthHTML = '<colgroup><col class="block" /><col class="block" /><col class="block" /><col class="block" /><col class="block" /></colgroup>'
		splitter = '<thead>'
		chunks = DETableHTML.split(splitter)
		assert len(chunks) == 2
		return str.join(splitter, [chunks[0] + columnWidthHTML, chunks[1]])
	
	def hackImagePathToSymlinkInStaticDir(old_path):
		"""
		Removes the string allJobsParDir (the full path of the jobs dir's parent dir) plus the remaining leading '/'
		from the old_path. This is done because the remaining path will be coupled to a job dir symlink in the static dir.
		:param old_path:	str		path to some child of the jobs dir
		:return:			str		tail of the input path, starting from the jobs dir.
		"""
		if old_path is not None:
			return old_path.split(allJobsParDir)[1].lstrip('/')
		else:
			return None
	
	# remove 'significant' columns
	for condition in jobParams['schema']['allConditions']:
		if 'significant' in minTopDifferentialsDFs[condition].columns:
			minTopDifferentialsDFs[condition] = minTopDifferentialsDFs[condition].drop('significant', axis=1, inplace=False)
		else:
			logging.warning("I cannot remove the 'significant' column from the minDE dataframe (it doesn't exist).")
		if 'significant' in minTopDifferentialsDFs[condition].columns:
			minTopDifferentialsDFs[condition] = minTopDifferentialsDFs[condition].drop('significant', axis=1, inplace=False)
		else:
			logging.warning("I cannot remove the 'significant' column from the fullDE dataframe (it doesn't exist).")
	
	# per condition: generate list of differentials HTML code separately because Jinja cant do this
	minTopDifferentialsHTMLDict = {(otherCondition, injectColumnWidthHTML(minTopDifferentialsDFs['otherCondition'].to_html(index=False, justify='left')))
								   for otherCondition in otherConditions}
	fullTopDifferentialsHTMLDict = {(otherCondition, injectColumnWidthHTML(fullTopDifferentialsDFs['otherCondition'].to_html(index=False, justify='left')))
								   for otherCondition in otherConditions}
	
	with open(logFilePath, 'r') as logFile:
		logContents = logFile.readlines()
	
	approxDuration = time() - startTime

	pdfhtmlreport = render_template('report.html', jobName=jobParams['jobName'], otherConditions=otherConditions,
									minVolcanoFullPathDict=minVolcanoFullPaths,
									fullVolcanoFullPathDict=fullVolcanoFullPaths,
									minExpression_bool=jobParams['minExpression_bool'],
									fullExpression_bool=jobParams['fullExpression_bool'],
									mindifferentialsdict=minTopDifferentialsHTMLDict,
									fulldifferentialsdict=fullTopDifferentialsHTMLDict, PCAFileName=PCAPlotFullPath,
									HCDFileName=HCDendrogramFullPath, metadata=metadata, date=jobParams['date'],
									duration=approxDuration, log=logContents, jobParams=jobParams,
									allProcessingParams=allProcessingParams, pdfsrc='True')#, experiments=experiments)
	# get the tails of the input paths, starting from the jobs dir, so the Jinja report template can couple it to the
	# jobs symlink in the static dir.
	for condition in otherConditions:
		minVolcanoFullPaths[condition] = hackImagePathToSymlinkInStaticDir(minVolcanoFullPaths[condition])
		fullVolcanoFullPaths[condition] = hackImagePathToSymlinkInStaticDir(fullVolcanoFullPaths[condition])
	HCDendrogramFullPath = hackImagePathToSymlinkInStaticDir(HCDendrogramFullPath)
	PCAPlotFullPath = hackImagePathToSymlinkInStaticDir(PCAPlotFullPath)
	htmlReport = render_template('report.html', jobName=jobParams['jobName'], otherConditions=otherConditions,
								 minVolcanoFullPathDict=minVolcanoFullPaths,
								 fullVolcanoFullPathDict=fullVolcanoFullPaths,
								 minExpression_bool=jobParams['minExpression_bool'],
								 fullExpression_bool=jobParams['fullExpression_bool'],
								 mindifferentialsdict=minTopDifferentialsHTMLDict,
								 fulldifferentialsdict=fullTopDifferentialsHTMLDict, PCAFileName=PCAPlotFullPath,
								 HCDFileName=HCDendrogramFullPath, metadata=metadata, date=jobParams['date'],
								 duration=approxDuration, log=logContents, jobParams=jobParams,
								 allProcessingParams=allProcessingParams)#, experiments=experiments)
	return htmlReport, pdfhtmlreport


def HTMLtoPDF(htmlReportFullPath):
	"""
	Generate a PDF file by converting the HTML report using weasyprint.
	:param htmlReportFullPath:	str		path to the HTML report file
	:return pdfReportFullPath:	str		path to the PDF report file
	"""
	# """
	# Generate a PDF file by converting the HTML report using the linux wkhtml2pdf package in a subprocess.
	# :param htmlReportFullPath:	str		path to the HTML report file
	# :return pdfReportFullPath:	str		path to the PDF report file
	# """
	from subprocess import run
	from weasyprint import HTML, CSS
	# from os import path
	# from constandpp_web import config
	
	pdfReportFullPath = htmlReportFullPath[0:-4] + 'pdf'
	# command = 'wkhtmltopdf -L 1cm -R 1cm -T 1cm -B 1cm "'+htmlReportFullPath+'" "'+pdfReportFullPath+'"'
	# command = 'weasyprint "'+htmlReportFullPath+'" "'+pdfReportFullPath+'"'# -s "'+path.abspath(config.__file__+'/../static/css/style.css')+'"'
	# run(command, shell=True
	HTML(htmlReportFullPath).write_pdf(pdfReportFullPath, stylesheets=[CSS(string='@page { size: A4; margin: 1cm; }')])
	rmcmd = 'rm -f "' + htmlReportFullPath + '"'
	run(rmcmd, shell=True)
	return pdfReportFullPath
