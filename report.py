#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Functions involved in generating the report that includes:
* volcano plot
* PCA plot
* HC tree
* list of differentially expressed proteins, sorted on fold change
* metadata info (warnings etc.)
* overview of parameters used in the workflow
"""
import pandas as pd
import numpy as np
from main import unnest
from warnings import warn
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt
from matplotlib import markers
#from adjustText import adjust_text

# save matplotlib images without whitespace: savefig('foo.png', bbox_inches='tight')


def distinguishableColours(n, type='prism'):
	"""
	Generates distinguishable colours in the format of a nConditionsx4 array.
	:param nConditions: int         number of colours
	:param type:        str         specifier of which colourmap to use
	:return:            np.ndarray  nConditionsx4 array
	"""
	cmap = plt.get_cmap(type)  # brg, jet
	return cmap(np.linspace(0, 1.0, n))


def getColours(schema):
	channelsPerConditionForAllExperiments = [len(condition) for experiment in schema.values() for condition in experiment['intensityColumnsPerCondition']]
	distColours = distinguishableColours(len(channelsPerConditionForAllExperiments))
	return unnest([[c] * n for c, n in zip(distColours, channelsPerConditionForAllExperiments)])


def distinguishableMarkers(n):
	"""
	Returns a list of n (easily) distinguishable markers, or just visible markers if n is too large. If n is extremely
	large, starts repeating markers.
	:param n:   int     number of distinguishable markers
	:return:    list    (easily) (distinguishable) markers
	"""
	easilyDistinguishable = ['o','s','x','*','v','d','+','^']
	allMarkers = markers.MarkerStyle.markers
	visibleMarkers = allMarkers
	for k,v in allMarkers.items():
		if v == 'nothing':
			del visibleMarkers[k]

	if n > len(easilyDistinguishable):
		if n < len(visibleMarkers):
			warn("More experiments than easily distinguishable markers; using all (visible) markers.")
			return [visibleMarkers[list(allMarkers.keys())[i]] for i in range(n)]
		else:
			warn("More experiments than markers. Using all (visible) markers with possible repetitions!")
			#number of times to re-use ALL visible markers
			nRepetitions = np.mod(len(visibleMarkers), n)
			nResidual = len(visibleMarkers) - n
			return list(visibleMarkers.keys())*nRepetitions + list(visibleMarkers.keys())[0:nResidual]
	else:
		return [visibleMarkers[i] for i in easilyDistinguishable[0:n]]


def getMarkers(schema):
	"""
	Returns list of markers for the channels in all experiments (based on schema) so that the channels of the same
	experiment have the same marker.
	:param schema:  dict    schema of the experiments
	:return:        list    markers per channel (markers differ only across experiments)
	"""
	channelsPerExperiment = [len(unnest(experiment['intensityColumnsPerCondition'])) for experiment in schema.values()]
	distMarkers = distinguishableMarkers(len(channelsPerExperiment))
	return unnest([[m]*n for m, n in zip(distMarkers, channelsPerExperiment)])


def getSortedDifferentialProteinsDF(df):
	"""
	Sorts the differential protein data according to absolute fold change and resets the index. Returns only the columns
	specified.
	:param df:  pd.DataFrame    unsorted
	:return:    pd.DataFrame    sorted according to fold change and only specified columns
	"""
	reportColumns = ['protein', 'significant', 'description', 'log2 fold change c1/c2', 'adjusted p-value']
	significantIndices = list(df[df['significant'] == 'yes'].index) + list(df[df['significant'] == 'p'].index)
	significantDf = df.loc[significantIndices, :]
	return significantDf.reindex(significantDf['log2 fold change c1/c2'].abs().sort_values(ascending=False).index).loc[:, reportColumns]


def getVolcanoPlot(df, alpha, FCThreshold, labelPlot=[False, ] * 4):
	# todo docu
	# todo add protein ID labels according to sorted list entry ID
	volcanoPlot = plt.figure(figsize=(6, 5))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
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

	# annotate where requested
	for labelPlotBool,xdata,ydata,labels in zip(labelPlot,[xdataYES, xdataP, xdataFC, xdataNO],
	                              [ydataYES, ydataP, ydataFC, ydataNO],
	                              [labelsYES, labelsP, labelsFC, labelsNO]):
		if labelPlotBool:
			for x, y, label in zip(xdata, ydata, labels):
				plt.annotate(label, xy=(x, y), xytext=(-1, 1), textcoords='offset points', ha='right', va='bottom')

	plt.show() # TEST
	return volcanoPlot


def getPCAPlot(PCAResult, schema):
	# todo docu
	# todo add experiment labels
	allChannelAliases = unnest([unnest(experiments['channelAliasesPerCondition']) for experiments in schema.values()])
	PCAPlot = plt.figure(figsize=(6, 5))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title('Principal Component scores', figure=PCAPlot)
	plt.xlabel('First PC', figure=PCAPlot)
	plt.ylabel('Second PC', figure=PCAPlot)

	# generate colors/markers so that the channels of the same condition/experiment have the same colour/markers


	# labels for annotation
	intensityColumns = [item for sublist in schema for item in sublist]
	# produce scatterplot and annotate
	for (x, y, color, marker, label) in zip(PCAResult[:, 0], PCAResult[:, 1], colors[eName], markers[eName], intensityColumns):
		plt.scatter(x, y, color=color, marker=marker, figure=PCAPlot)  # plot first two principal components
		plt.annotate(label, xy=(x, y), xytext=(-1, 1),
			textcoords='offset points', ha='right', va='bottom')
	plt.show() # TEST
	return PCAPlot


def getHCDendrogram(HCResult, intensityColumnsPerCondition):
	# todo docu
	# hierarchical clustering dendrogram
	intensityColumns = [item for sublist in intensityColumnsPerCondition for item in sublist]
	HCDendrogram = plt.figure(figsize=(6, 5)) # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title('Hierarchical Clustering Dendrogram', figure=HCDendrogram)
	plt.xlabel('reporter channel', figure=HCDendrogram)
	plt.ylabel('distance', figure=HCDendrogram)
	dendrogram(HCResult, leaf_rotation=0., leaf_font_size=12., labels=intensityColumns)
	return HCDendrogram
	plt.show() # TEST


def makeHTML(minSortedDifferentialProteinsDF, fullSortedDifferentialProteinsDF, diffMinFullProteins,
	                      visualizationsDict, metadata):
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
	htmlReport = None
	return htmlReport


def HTMLtoPDF(htmlReport):
	# todo docu
	pdfReport = htmlReport
	return pdfReport