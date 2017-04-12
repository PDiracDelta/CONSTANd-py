#!/usr/bin/python
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
from dataIO import unnest
from warnings import warn
from scipy.cluster.hierarchy import dendrogram
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import markers
#from adjustText import adjust_text

# adjust font size globally
fontsize = 30
fontweight = 'normal'
matplotlib.rcParams.update({'font.size': fontsize, 'font.weight': fontweight})

# save matplotlib images without whitespace: savefig('foo.png', bbox_inches='tight')


def distinguishableColours(n, type='jet'):
	"""
	Generates distinguishable colours in the format of a n array.
	:param nConditions: int         number of colours
	:param type:        str         specifier of which colourmap to use
	:return:            np.ndarray  nx4 array
	"""
	cmap = plt.get_cmap(type)  # brg, jet
	return cmap(np.linspace(0, 1.0, n))


def getColours(schema, allChannelAliases):
	"""
	Returns list of colours for all the channels in all experiments (based on schema) so that the channels of the same
	condition have the same colour.
	:param schema:  			dict    schema of the experiments' hierarchy' hierarchy
	:return channelColoursDict:	dict    colour for each channel; a different one for each condition
	"""
	numConditions = len(list(schema.values())[0]['channelAliasesPerCondition'])
	distColours = distinguishableColours(numConditions)
	colours = []
	for experiment in schema.values():
		# repeat each channel's distColour as many times as there are channels in the current condition, and repeat for each experiment
		conditions = experiment['channelAliasesPerCondition']
		colours.append([np.tile(distColours[c], (len(conditions[c]),1)).tolist() for c in range(numConditions)])
	channelColoursDict = dict(zip(allChannelAliases, unnest(unnest(colours))))
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
	for k,v in list(allMarkers.items()):
		if v == 'nothing':
			del visibleMarkers[k]

	if n > len(easilyDistinguishable): # not enough distinguishable markers
		if n < len(visibleMarkers): # enough visible markers
			warn("More experiments than easily distinguishable markers; using all (visible) markers.")
			return [list(allMarkers.keys())[i] for i in range(n)]
		else: # not enough markers at all
			warn("More experiments than markers. Using all (visible) markers with possible repetitions!")
			#number of times to re-use ALL visible markers
			nRepetitions = np.mod(len(visibleMarkers), n)
			nResidual = len(visibleMarkers) - n
			return list(visibleMarkers.keys())*nRepetitions + list(visibleMarkers.keys())[0:nResidual]
	else: # have enough distinguishable markers for the n experiments
		return easilyDistinguishable[0:n]


def getMarkers(schema, allChannelAliases):
	"""
	Returns list of markers for the channels in all experiments (based on schema) so that the channels of the same
	experiment have the same marker.
	:param schema:  dict    schema of the experiments' hierarchy
	:return:        dict
	"""
	# todo docu
	distMarkers = distinguishableMarkers(len(schema))
	channelMarkersDict = {}
	i = 0
	for experiment in schema.values():
		for alias in unnest(experiment['channelAliasesPerCondition']):
			channelMarkersDict[alias] = distMarkers[i]
		i += 1
	return channelMarkersDict


def getSortedDifferentialProteinsDF(df):
	"""
	Sorts the differential protein data according to adjusted p-value and resets the index. Returns only the columns
	specified.
	:param df:  pd.DataFrame    unsorted
	:return:    pd.DataFrame    sorted according to adjusted p-value and only specified columns
	"""
	reportColumns = ['protein', 'significant', 'description', 'log2 fold change c1/c2', 'adjusted p-value']
	significantIndices = list(df[df['significant'] == 'yes'].index) + list(df[df['significant'] == 'p'].index)
	significantDf = df.loc[significantIndices, :]
	return significantDf.reindex(significantDf['adjusted p-value'].sort_values(ascending=True).index).loc[:, reportColumns]


def getVolcanoPlot(df, alpha, FCThreshold, labelPlot=[False, ] * 4):
	# todo docu
	# todo add protein ID labels according to sorted list entry ID
	volcanoPlot = plt.figure(figsize=(16, 12))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	# maximize figure
	#mng = plt.get_current_fig_manager()
	#mng.full_screen_toggle()
	plt.title(r'Volcano Plot ($FC>$' + str(FCThreshold) + r'; $\alpha=$' + str(alpha) + ')', figure=volcanoPlot)
	plt.xlabel(r'log$_2$(fold change)', figure=volcanoPlot)
	plt.ylabel(r'-log$_{10}$(p-value) ', figure=volcanoPlot)
	# get indices of different levels of significance
	significantIndices_yes = df[df['significant'] == 'yes'].index
	significantIndices_p = df[df['significant'] == 'p'].index
	significantIndices_fc = df[df['significant'] == 'fc'].index
	significantIndices_no = df[df['significant'] == 'no'].index

	# produce scatterplot for each category of significance
	# YES
	xdataYES = df.loc[significantIndices_yes, 'log2 fold change c1/c2']
	ydataYES = -np.log10(df.loc[significantIndices_yes, 'adjusted p-value'])

	labelsYES = df.loc[significantIndices_yes, 'protein']
	plt.scatter(xdataYES, ydataYES, color='r', figure=volcanoPlot)
	# P
	xdataP = df.loc[significantIndices_p, 'log2 fold change c1/c2']
	ydataP = -np.log10(df.loc[significantIndices_p, 'adjusted p-value'])
	labelsP = df.loc[significantIndices_p, 'protein']
	plt.scatter(xdataP, ydataP, color='b', figure=volcanoPlot)
	# FC
	xdataFC = df.loc[significantIndices_fc, 'log2 fold change c1/c2']
	ydataFC = -np.log10(df.loc[significantIndices_fc, 'adjusted p-value'])
	labelsFC = df.loc[significantIndices_fc, 'protein']
	plt.scatter(xdataFC, ydataFC, color='g', figure=volcanoPlot)
	# NO
	xdataNO = df.loc[significantIndices_no, 'log2 fold change c1/c2']
	ydataNO = -np.log10(df.loc[significantIndices_no, 'adjusted p-value'])
	labelsNO = df.loc[significantIndices_no, 'protein']
	plt.scatter(xdataNO, ydataNO, color='k', figure=volcanoPlot)
	
	# adjust limits
	try:
		plt.xlim([int(min(min(xdataFC), min(xdataYES))), np.ceil(max(max(xdataFC), max(xdataYES)))])
	except ValueError:
		plt.xlim([-5,5])
	try:
		plt.ylim([0, np.ceil(max(max(ydataP), max(ydataYES))/5)*5]) # int(base * round(float(x)/base))
	except ValueError:
		plt.ylim([0, 100])  # int(base * round(float(x)/base))
	# annotate where requested
	for labelPlotBool,xdata,ydata,labels in zip(labelPlot,[xdataYES, xdataP, xdataFC, xdataNO],
	                              [ydataYES, ydataP, ydataFC, ydataNO],
	                              [labelsYES, labelsP, labelsFC, labelsNO]):
		if labelPlotBool:
			for x, y, label in zip(xdata, ydata, labels):
				plt.annotate(label, xy=(x, y), xytext=(-1, 1), textcoords='offset points', ha='right', va='bottom', fontsize=20)

	#plt.show() # TEST
	return volcanoPlot


def getPCAPlot(PCAResult, schema):
	# todo docu
	PCAPlot = plt.figure(figsize=(16, 12))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	# maximize figure
	#mng = plt.get_current_fig_manager()
	#mng.full_screen_toggle()
	plt.title('Principal Component scores', figure=PCAPlot)
	plt.xlabel('First PC', figure=PCAPlot)
	plt.ylabel('Second PC', figure=PCAPlot)

	# labels for annotation
	allChannelAliases = unnest([unnest(experiments['channelAliasesPerCondition']) for experiments in schema.values()])
	#print(allChannelAliases) # TEST
	# generate colors/markers so that the channels of the same condition/experiment have the same colour/markers
	channelColorsDict = getColours(schema, allChannelAliases)
	channelMarkersDict = getMarkers(schema, allChannelAliases)

	#xmin, xmax, ymin, ymax = min(PCAResult[:, 0]), max(PCAResult[:, 0]), min(PCAResult[:, 1]), max(PCAResult[:, 1])
	for (x, y, label) in zip(PCAResult[:, 0], PCAResult[:, 1], allChannelAliases):
		# produce scatterplot of two first principal components and annotate
		plt.scatter(x, y, color=channelColorsDict[label], marker=channelMarkersDict[label], figure=PCAPlot, s=80)
		plt.annotate(label, xy=(x, y), xytext=(-1, 1),
			textcoords='offset points', ha='right', va='bottom', fontsize=20)
	legendHandles = []
	legendStrings = []
	# look for corresponding experiment name
	markersToCheck = set(channelMarkersDict.values())
	for channel, marker in channelMarkersDict.items():
		if marker in markersToCheck:
			for eName, experiment in schema.items():
				if channel in unnest(experiment['channelAliasesPerCondition']):
					handle = plt.scatter([], [], color='k', marker=marker, s=160)
					legendHandles.append(handle)
					legendStrings.append(eName)
					markersToCheck.remove(marker)
					break
	plt.legend(legendHandles, legendStrings, scatterpoints=1)#, loc=2)
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
	#plt.axhspan(xmin=xmin-0.1(xmax-xmin), ymax=ymax+0.05*(ymax-ymin))
	#plt.show() # TEST
	return PCAPlot


def getHCDendrogram(HCResult, schema):
	# todo docu
	# hierarchical clustering dendrogram
	allChannelAliases = unnest([unnest(experiments['channelAliasesPerCondition']) for experiments in schema.values()])
	HCDendrogram = plt.figure(figsize=(16, 12)) # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	# maximize figure
	#mng = plt.get_current_fig_manager()
	#mng.full_screen_toggle()
	plt.title('Hierarchical Clustering Dendrogram', figure=HCDendrogram)
	plt.xlabel('distance', figure=HCDendrogram)
	plt.ylabel('reporter channel', figure=HCDendrogram)
	# generate colors/markers so that the channels of the same condition/experiment have the same colour/markers
	channelColorsDict = getColours(schema, allChannelAliases)
	dendrogram(HCResult, orientation='right', leaf_rotation=0., leaf_font_size=24, labels=allChannelAliases,
	           link_color_func=lambda x: channelColorsDict[allChannelAliases[x]] if x < len(allChannelAliases) else 'k',
	           above_threshold_color='k')
	ax = plt.gca()
	ylbls = ax.get_ymajorticklabels()
	for i in range(len(ylbls)):
		ylbls[i].set_color(channelColorsDict[ylbls[i].get_text()])
		#ylbls[i].set_color(colorsPerCondition[i])
	#plt.show()  # TEST
	return HCDendrogram


def makeHTML(jobParams, processingParams, minSortedDifferentialProteinsDF, fullSortedDifferentialProteinsDF, minVolcanoFullPath,
             fullVolcanoFullPath, PCAPlotFullPath, HCDendrogramFullPath, metadata, logFilePath, startTime):
	"""
	Pour all report ingredients into an HTML file.
	:param minSortedDifferentialProteinsDF:     pd.DataFrame
	:param fullSortedDifferentialProteinsDF:    pd.DataFrame
	:param diffMinFullProteins:                 nested list
	:param visualizationsDict:                  dict
	:param metadata:                            dict
	:return:
	"""
	# todo docu
	from flask import render_template
	from time import time
	from web import app
	from os import path, pardir
	allJobdsParDir = path.abspath(path.join(app.config.get('ALLJOBSDIR'), pardir))

	def hackImagePathToSymlinkInStaticDir(old_path):
		if old_path is not None:
			return old_path.split(allJobdsParDir)[1].lstrip('/')
		else:
			return None

	numDifferentials = jobParams['numDifferentials']
	if jobParams['minExpression_bool']:
		if numDifferentials < len(minSortedDifferentialProteinsDF):
			minTopDifferentials = minSortedDifferentialProteinsDF.loc[range(numDifferentials), :]
		else:
			minTopDifferentials = minSortedDifferentialProteinsDF
	else:
		minTopDifferentials = pd.DataFrame(columns=minSortedDifferentialProteinsDF.columns)
	if jobParams['fullExpression_bool']:
		if numDifferentials < len(fullSortedDifferentialProteinsDF):
			fullTopDifferentials = fullSortedDifferentialProteinsDF.loc[range(numDifferentials), :]
		else:
			fullTopDifferentials = fullSortedDifferentialProteinsDF
	else:
		fullTopDifferentials = pd.DataFrame(columns=fullSortedDifferentialProteinsDF.columns)
	with open(logFilePath, 'r') as logFile:
		logContents = logFile.readlines()

	approxDuration = time() - startTime
	experiments = jobParams['schema']
	for e in experiments:
		experiments[e]['cond1Aliases'] = experiments[e]['channelAliasesPerCondition'][0]
	pdfhtmlreport = render_template('report.html', jobName=jobParams['jobname'], minVolcanoFullPath=minVolcanoFullPath,
	                             fullVolcanoFullPath=fullVolcanoFullPath, minExpression_bool=jobParams['minExpression_bool'],
	                             fullExpression_bool=jobParams['fullExpression_bool'], mindifferentials=minTopDifferentials,
	                             fulldifferentials=fullTopDifferentials, PCAFileName=PCAPlotFullPath,
	                             HCDFileName=HCDendrogramFullPath, metadata=metadata, date=jobParams['date'],
	                             duration=approxDuration, log=logContents, jobParams=jobParams,
	                             processingParams=processingParams, experiments=experiments, pdfsrc='True')
	minVolcanoFullPath = hackImagePathToSymlinkInStaticDir(minVolcanoFullPath)
	fullVolcanoFullPath = hackImagePathToSymlinkInStaticDir(fullVolcanoFullPath)
	HCDendrogramFullPath = hackImagePathToSymlinkInStaticDir(HCDendrogramFullPath)
	PCAPlotFullPath = hackImagePathToSymlinkInStaticDir(PCAPlotFullPath)
	htmlReport = render_template('report.html', jobName=jobParams['jobname'], minVolcanoFullPath=minVolcanoFullPath,
	                             fullVolcanoFullPath=fullVolcanoFullPath,
	                             minExpression_bool=jobParams['minExpression_bool'],
	                             fullExpression_bool=jobParams['fullExpression_bool'],
	                             mindifferentials=minTopDifferentials,
	                             fulldifferentials=fullTopDifferentials, PCAFileName=PCAPlotFullPath,
	                             HCDFileName=HCDendrogramFullPath, metadata=metadata, date=jobParams['date'],
	                             duration=approxDuration, log=logContents, jobParams=jobParams,
	                             processingParams=processingParams, experiments=experiments)
	return htmlReport, pdfhtmlreport


def HTMLtoPDF(htmlReportFullPath):
	# todo docu
	from subprocess import run

	pdfReportFullPath = htmlReportFullPath[0:-4]+'pdf'
	command = 'wkhtmltopdf -L 1cm -R 1cm -T 1cm -B 1cm "'+htmlReportFullPath+'" "'+pdfReportFullPath+'"'
	run(command, shell=True)
	rmcmd = 'rm -f "'+htmlReportFullPath+'"'
	run(rmcmd, shell=True)
	return pdfReportFullPath