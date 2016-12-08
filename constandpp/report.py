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
from warnings import warn
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt
#from adjustText import adjust_text

# save matplotlib images without whitespace: savefig('foo.png', bbox_inches='tight')


def distColours(n, type='prism'):
	"""
	Generates distinguishable colours in the format of a nConditionsx4 array.
	:param nConditions: int         number of colours
	:param type:        str         specifier of which colourmap to use
	:return:            np.ndarray  nConditionsx4 array
	"""
	cmap = plt.get_cmap(type)  # brg, jet
	return cmap(np.linspace(0, 1.0, n))


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
	PCAPlot = plt.figure(figsize=(6, 5))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title('Principal Component scores', figure=PCAPlot)
	plt.xlabel('First PC', figure=PCAPlot)
	plt.ylabel('Second PC', figure=PCAPlot)
	experiments = list(schema.keys()) # names of experiments
	nConditions = len(schema[experiments[0]]['intensityColumnsPerCondition'])  # number of conditions
	# get distinguishable colours
	distinguishableColors = distColours(nConditions)
	# generate colors/shapes so that the channels of the same condition/experiment have the same colour/shape
	colors = {}
	shapes = {}
	for experiment in range(len(experiments)):
		shapes[experiment] = []
		for condition in range(nConditions):  # for each condition a different color
			for i in range(len(schema[experiment][condition])):  # add the color for each channel per condition
				colors.append(distinguishableColors[condition])

	# labels for annotation
	intensityColumns = [item for sublist in schema for item in sublist]
	# produce scatterplot and annotate
	for (x, y, color, label) in zip(PCAResult[:, 0], PCAResult[:, 1], colors, intensityColumns):
		plt.scatter(x, y, color=color, figure=PCAPlot)  # plot first two principal components
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