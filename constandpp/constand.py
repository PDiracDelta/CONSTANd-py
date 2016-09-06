#!/usr/bin/python
# -*- coding: utf-8 -*-

# function [normalizedData,f,R,S] = CONSTANd_RAS(data,h,maxIterations)
"""
CONSTANd
Normalizes the data matrix <data> by raking the m by n matrix such that
the row mean and column mean equals to 1/n. Missing information needs to
be presented as NaN values and not as zero values because CONSTANd
employs the Matlab functionality 'nanmean' that is able to ignore
NaN-values when calculating the mean. The variable <nbIterations> is an
integer value that denotes the number of raking cycles. The variable <h>
defines the stopping criteria based on the L1-norm as defined by
Friedrich Pukelsheim, Bruno Simeone in "On the Iterative Proportional
Fitting Procedure: Structure of Accumulation Points and L1-Error
Analysis"

© Dirk Valkenborg & Jef Hooyberghs, 2014
Ported to Python by Joris Van Houtven, 2016
"""

