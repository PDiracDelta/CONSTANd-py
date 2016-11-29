#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Functions involved in generating the report that includes:
* volcano plot
* PCA plot
* HC tree
* list of differentially expressed proteins, sorted on fold change or probability (p-value)
* metadata info (warnings etc.)
* overview of parameters used in the workflow
"""
import pandas as pd
from warnings import warn

def dataVisualization(minProteinDF, fullProteinDF, FCThreshold, alpha):
	# TODO (if paying customer): parameter: intensity matrix on peptide or protein level?
	# TODO: only include differentials with a fold of >threshold or <1/threshold
	return None