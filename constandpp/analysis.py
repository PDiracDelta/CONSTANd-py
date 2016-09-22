#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Performs a differential expression analysis on the normalized intensities as provided by CONSTANd.
Includes data visualization.
"""


def differentialExpression(normalizedIntensities, threshold=1):
	# TODO: only include differentials with a fold of >threshold or <1/threshold
	# TODO: careful with peptides with more than 1 master protein
	return None


# pept2prot = mapping of peptides to proteins (1 to many) --> extract directly from df
# prot2peptMin = mapping of proteins to peptides (1 to many) --> calculate from pept2prot
# prot2peptMax = mapping of proteins to peptides (1 to many more) --> calculate from pept2prot
def mapPeptToProt(df, peptides):
	# calculate prot2peptMin/Max and alongside it protIntensitiesMin/Max:
	# Min is the case where proteins are only associated with 1 to 1 matching peptides in pept2prot,
	# Max is the case where proteins are associated with all matching peptides in pept2prot
	return protIntensitiesMin, protIntensitiesMax, prot2peptMin, prot2peptMax

def dataVisualization(DEresults):
	# TODO (if paying customer): parameter: intensity matrix on peptide or protein level?
	return None



