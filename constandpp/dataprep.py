#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Collection of functions that prepare the data before it can be normalized by CONSTANd. Includes:
* removing redundancy due to:
	* different peptide spectrum match (PSM) algorithms
	* different retention time (RT) values
	* different charges
* correct for isotopic impurities in the reporters
"""


def collapsePSMAlgo(df):
	return None


def collapseRT(df):
	return None


def collapseCharge(df):
	return None


def isotopicCorrection(intensities):
	return None
