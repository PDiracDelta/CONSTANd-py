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
from adjustText import adjust_text

# save matplotlib images without whitespace: savefig('foo.png', bbox_inches='tight')


def getSortedDifferentials(df):
	"""
	Sorts the differential protein data according to fold change (and p-value as secondary).
	:param df:  pd.DataFrame    unsorted
	:return:    pd.DataFrame    sorted according to fold change (and p-value as secondary)
	"""
	significantIndices = list(df[df['significant'] == 'yes'].index) + list(df[df['significant'] == 'p'].index)
	return df.loc[significantIndices, :].sort_values(by=['log2 fold change c1/c2', 'adjusted p-value'],
	                                                 ascending=[False, True], axis=0)


def getVolcanoPlot(minProteinDF, fullProteinDF, alpha, FCThreshold, labelPlot=[False,]*4):
	# todo docu
	# todo add minProteinDF and combine into one plot
	# todo add protein ID labels according to sorted list entry ID
	volcanoPlot = plt.figure(figsize=(6, 5))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title(r'Volcano Plot ($FC>$' + str(FCThreshold) + r'; $\alpha=$' + str(alpha) + ')', figure=volcanoPlot)
	plt.xlabel(r'log$_2$(fold change)', figure=volcanoPlot)
	plt.ylabel(r'-log$_{10}$(p-value) ', figure=volcanoPlot)
	# get indices of different levels of significance
	significantIndices_yes = fullProteinDF[fullProteinDF['significant'] == 'yes'].index
	significantIndices_p = fullProteinDF[fullProteinDF['significant'] == 'p'].index
	significantIndices_fc = fullProteinDF[fullProteinDF['significant'] == 'fc'].index
	significantIndices_no = fullProteinDF[fullProteinDF['significant'] == 'no'].index

	# produce scatterplot for each category of significance
	# YES (also annotate with sorted list ID)
	xdataYES = fullProteinDF.loc[significantIndices_yes, 'log2 fold change c1/c2']
	ydataYES = -np.log10(fullProteinDF.loc[significantIndices_yes, 'adjusted p-value'])
	labelsYES = fullProteinDF.loc[significantIndices_yes, 'protein']
	plt.scatter(xdataYES, ydataYES, color='r', figure=volcanoPlot)
	if labelPlot[0]:
		#textsYES = []
		for x,y,label in zip(xdataYES, ydataYES, labelsYES):
			plt.annotate(label, xy=(x, y), xytext=(-1, 1), textcoords='offset points', ha='right', va='bottom')
			#textsYES.append(plt.text(x, y, label))
			#adjust_text(xdataYES, ydataYES, textsYES, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
	# P (also annotate with sorted list ID)
	xdataP = fullProteinDF.loc[significantIndices_p, 'log2 fold change c1/c2']
	ydataP = -np.log10(fullProteinDF.loc[significantIndices_p, 'adjusted p-value'])
	labelsP = fullProteinDF.loc[significantIndices_p, 'protein']
	plt.scatter(fullProteinDF.loc[significantIndices_p, 'log2 fold change c1/c2'],
	            -np.log10(fullProteinDF.loc[significantIndices_p, 'adjusted p-value']),
	            color='b', figure=volcanoPlot)
	if labelPlot[1]:
		for x,y,label in zip(xdataP, ydataP, labelsP):
			plt.annotate(label, xy=(x, y), xytext=(-1, 1), textcoords='offset points', ha='right', va='bottom')
	# FC
	xdataFC = fullProteinDF.loc[significantIndices_fc, 'log2 fold change c1/c2']
	ydataFC = -np.log10(fullProteinDF.loc[significantIndices_fc, 'adjusted p-value'])
	labelsFC = fullProteinDF.loc[significantIndices_fc, 'protein']
	plt.scatter(fullProteinDF.loc[significantIndices_fc, 'log2 fold change c1/c2'],
	            -np.log10(fullProteinDF.loc[significantIndices_fc, 'adjusted p-value']),
	            color='g', figure=volcanoPlot)
	if labelPlot[2]:
		for x,y,label in zip(xdataFC, ydataFC, labelsFC):
			plt.annotate(label, xy=(x, y), xytext=(-1, 1), textcoords='offset points', ha='right', va='bottom')
	# NO
	xdataNO = fullProteinDF.loc[significantIndices_no, 'log2 fold change c1/c2']
	ydataNO = -np.log10(fullProteinDF.loc[significantIndices_no, 'adjusted p-value'])
	labelsNO = fullProteinDF.loc[significantIndices_no, 'protein']
	plt.scatter(fullProteinDF.loc[significantIndices_no, 'log2 fold change c1/c2'],
	            -np.log10(fullProteinDF.loc[significantIndices_no, 'adjusted p-value']),
	            color='k', figure=volcanoPlot)
	if labelPlot[3]:
		for x,y,label in zip(xdataNO, ydataNO, labelsNO):
			plt.annotate(label, xy=(x, y), xytext=(-1, 1), textcoords='offset points', ha='right', va='bottom')
	plt.show()
	return volcanoPlot


def getPCAPlot(PCAResult, intensityColumnsPerCondition):
	# todo docu
	# todo add experiment labels
	PCAPlot = plt.figure(figsize=(6, 5))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title('Principal Component scores', figure=PCAPlot)
	plt.xlabel('First PC', figure=PCAPlot)
	plt.ylabel('Second PC', figure=PCAPlot)
	# get distinguishable colours
	cmap = plt.get_cmap('prism')  # brg, jet
	nConditions = len(intensityColumnsPerCondition)  # number of TMT channels
	distinguishableColors = cmap(np.linspace(0, 1.0, nConditions))
	# generate colors vector so that the channels of the same condition have the same colour
	colors = []
	for condition in range(nConditions):  # for each condition a different color
		for i in range(len(intensityColumnsPerCondition[condition])):  # add the color for each channel per condition
			colors.append(distinguishableColors[condition])
	# labels for annotation
	intensityColumns = [item for sublist in intensityColumnsPerCondition for item in sublist]
	# produce scatterplot and annotate
	for (x, y, color, label) in zip(PCAResult[:, 0], PCAResult[:, 1], colors, intensityColumns):
		plt.scatter(x, y, color=color, figure=PCAPlot)  # plot first two principal components
		plt.annotate(label, xy=(x, y), xytext=(-1, 1),
			textcoords='offset points', ha='right', va='bottom')
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
	plt.show()  # TEST
