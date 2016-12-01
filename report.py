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
	volcano = plt.figure(figsize=(6, 5)) # size(inches wide, height); a4paper: width = 8.267in; height 11.692in

	# PCA plot
	PCAPlot = plt.figure(figsize=(6, 5)) # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title('Principal Component scores', figure=PCAPlot)
	plt.xlabel('First PC', figure=PCAPlot)
	plt.ylabel('Second PC', figure=PCAPlot)
	# get distinguishable colours
	hsv = plt.get_cmap('hsv')
	nConditions = len(intensityColumnsPerCondition) # number of TMT channels
	distinguishableColors = hsv(np.linspace(0, 1.0, nConditions))
	# generate colors vector so that the channels of the same condition have the same colour
	colors = []
	for condition in range(nConditions):
		colors.append([distinguishableColors[condition],]*len(intensityColumnsPerCondition[condition]))
	# produce scatterplot
	plt.scatter(PCAResult[:, 0], PCAResult[:, 1], c=colors, figure=PCAPlot) # plot first two principal components
	visualizationsDict['pca'] = PCAPlot
	plt.show() # TEST

	# hierarchical clustering dendrogram
	HCDendrogram = plt.figure(figsize=(6, 5)) # size(inches wide, height); a4paper: width = 8.267in; height 11.692in
	plt.title('Hierarchical Clustering Dendrogram', figure=HCDendrogram)
	plt.xlabel('reporter channel', figure=HCDendrogram)
	plt.ylabel('distance', figure=HCDendrogram)
	#dendrogram(HCResult, leaf_rotation=0., leaf_font_size=12., figure=HCDendrogram) # TEST
	visualizationsDict['hcd'] = HCDendrogram
	return visualizationsDict