#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Workflow of the processing part of CONSTANd++.
"""

from report import *
from dataIO import exportData


def generateReport(analysisResults, params, logFilePath, writeToDisk):
	# todo docu
	minProteinDF = analysisResults[0]
	fullProteinDF = analysisResults[1]
	PCAResult = analysisResults[2]
	HCResult = analysisResults[3]
	metadata = analysisResults[4]

	nConditions = len(list(params['schema'].values())[0]['channelAliasesPerCondition'])
	# ONLY PRODUCE VOLCANO AND DEA IF CONDITIONS == 2
	if nConditions == 2:
		# generate sorted (on FC) list of differentials
		minSortedDifferentialProteinsDF = getSortedDifferentialProteinsDF(minProteinDF)
		fullSortedDifferentialProteinsDF = getSortedDifferentialProteinsDF(fullProteinDF)
		minSet = set(minSortedDifferentialProteinsDF['protein'])
		fullSet = set(fullSortedDifferentialProteinsDF['protein'])
		# list( [in min but not in full], [in full but not in min] )
		metadata['diffMinFullProteins'] = [list(minSet.difference(fullSet)), list(fullSet.difference(minSet))]
		# todo combine into one

		# data visualization
		minVolcanoPlot = getVolcanoPlot(minProteinDF, params['alpha'], params['FCThreshold'],
		                                               params['labelVolcanoPlotAreas'])
		if writeToDisk:
			exportData(minVolcanoPlot, dataType='fig', path_out=params['path_results'],
			           filename=params['jobname'] + '_minVolcanoPlot')
		fullVolcanoPlot = getVolcanoPlot(fullProteinDF, params['alpha'], params['FCThreshold'],
		                                                  params['labelVolcanoPlotAreas'])
	else:
		minSortedDifferentialProteinsDF = pd.DataFrame()
		fullSortedDifferentialProteinsDF = pd.DataFrame()
		minVolcanoPlot = None
		fullVolcanoPlot = None

	if writeToDisk:
		exportData(fullVolcanoPlot, dataType='fig', path_out=params['path_results'],
		           filename=params['jobname'] + '_fullVolcanoPlot')
	PCAPlot = getPCAPlot(PCAResult, params['schema'])
	if writeToDisk:
		exportData(PCAPlot, dataType='fig', path_out=params['path_results'],
		           filename=params['jobname'] + '_PCAPlot')
	HCDendrogram = getHCDendrogram(HCResult, params['schema'])
	if writeToDisk:
		exportData(HCDendrogram, dataType='fig', path_out=params['path_results'],
		           filename=params['jobname'] + '_HCDendrogram')

	# generate HTML and PDF reports # todo
	htmlReport = makeHTML(minSortedDifferentialProteinsDF, fullSortedDifferentialProteinsDF, minVolcanoPlot,
	                      fullVolcanoPlot, PCAPlot, HCDendrogram, metadata, logFilePath)
	pdfReport = HTMLtoPDF(htmlReport)

