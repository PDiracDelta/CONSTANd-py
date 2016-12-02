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

# save matplotlib images without whitespace: savefig('foo.png', bbox_inches='tight')


def dataVisualization(minProteinDF, fullProteinDF, alpha, FCThreshold, PCAResult, HCResult, intensityColumnsPerCondition):
	# TODO (if paying customer): parameter: intensity matrix on peptide or protein level?
	# TODO: only include differentials with a fold of >threshold or <1/threshold
	visualizationsDict = {}

	# volcano plot
	# todo add protein ID labels according to sorted list entry ID
	volcanoPlot = plt.figure(figsize=(6, 5))  # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title(r'Volcano Plot ($FC>$'+str(FCThreshold)+r'; $\alpha=$'+str(alpha)+')', figure=volcanoPlot)
	plt.xlabel(r'log$_2$(fold change)', figure=volcanoPlot)
	plt.ylabel(r'-log$_{10}$(p-value) ', figure=volcanoPlot)
	# get indices of different levels of significance
	significantIndices_yes = fullProteinDF[fullProteinDF['significant'] == 'yes'].index
	significantIndices_p = fullProteinDF[fullProteinDF['significant'] == 'p'].index
	significantIndices_fc = fullProteinDF[fullProteinDF['significant'] == 'fc'].index
	significantIndices_no = fullProteinDF[fullProteinDF['significant'] == 'no'].index
	# produce scatterplot for each category of significance
	plt.scatter(fullProteinDF.loc[significantIndices_yes, 'log2 fold change c1/c2'],
	            -np.log10(fullProteinDF.loc[significantIndices_yes, 'adjusted p-value']),
	            color='r', figure=volcanoPlot)
	plt.scatter(fullProteinDF.loc[significantIndices_p, 'log2 fold change c1/c2'],
	            -np.log10(fullProteinDF.loc[significantIndices_p, 'adjusted p-value']),
	            color='b', figure=volcanoPlot)
	plt.scatter(fullProteinDF.loc[significantIndices_fc, 'log2 fold change c1/c2'],
	            -np.log10(fullProteinDF.loc[significantIndices_fc, 'adjusted p-value']),
	            color='g', figure=volcanoPlot)
	plt.scatter(fullProteinDF.loc[significantIndices_no, 'log2 fold change c1/c2'],
	            -np.log10(fullProteinDF.loc[significantIndices_no, 'adjusted p-value']),
	            color='k', figure=volcanoPlot)
	visualizationsDict['pca'] = volcanoPlot

	# PCA plot
	# todo add experiment labels
	PCAPlot = plt.figure(figsize=(6, 5)) # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title('Principal Component scores', figure=PCAPlot)
	plt.xlabel('First PC', figure=PCAPlot)
	plt.ylabel('Second PC', figure=PCAPlot)
	# get distinguishable colours
	cmap = plt.get_cmap('prism') # brg, jet
	nConditions = len(intensityColumnsPerCondition) # number of TMT channels
	distinguishableColors = cmap(np.linspace(0, 1.0, nConditions))
	# generate colors vector so that the channels of the same condition have the same colour
	colors = []
	for condition in range(nConditions): # for each condition a different color
		for i in range(len(intensityColumnsPerCondition[condition])): # add the color for each channel per condition
			colors.append(distinguishableColors[condition])
	# produce scatterplot
	for (x,y,color) in zip(PCAResult[:, 0], PCAResult[:, 1], colors):
		plt.scatter(x, y, color=color, figure=PCAPlot) # plot first two principal components
	visualizationsDict['pca'] = PCAPlot

	# hierarchical clustering dendrogram
	HCDendrogram = plt.figure(figsize=(6, 5)) # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title('Hierarchical Clustering Dendrogram', figure=HCDendrogram)
	plt.xlabel('reporter channel', figure=HCDendrogram)
	plt.ylabel('distance', figure=HCDendrogram)
	dendrogram(HCResult, leaf_rotation=0., leaf_font_size=12.) # TEST
	visualizationsDict['hcd'] = HCDendrogram
	plt.show()  # TEST

	return visualizationsDict
