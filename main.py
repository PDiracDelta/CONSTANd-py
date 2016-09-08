#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Python implementation of mass spectrometer protein data analysis using the CONSTANd_RAS algorithm
"""

__author__ = "Joris Van Houtven"
__copyright__ = "Copyright ?, VITO"
__credits__ = ["Joris Van Houtven", "Dirk Valkenborg"]
# __license__ = "GPL"
# __version__ = "0.0.1"
__maintainer__ = "Joris Van Houtven"
__email__ = "vanhoutvenjoris@gmail.com"
__status__ = "Development"

import sys
import pandas as pd
import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
import constand as cd


def getInput():
    """Get mass spec data and CONSTANd parameters"""
    path=None
    sep='\t'
    accuracy=1e-2
    maxIterations=50
    return path,sep,accuracy,maxIterations


def importDataFrame(path=None, sep=','):
    """Get the data from disk as a Pandas DataFrame"""
    df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(10, 6)), columns=list('ABCDEF'))  # TEST
    # df = pd.read_csv('../data/MB_Bon_tmt_TargetPeptideSpectrumMatch.txt', sep='\t') # TEST
    # df = pd.read_csv(path, sep=sep)
    df = pd.DataFrame(np.arange(10*6).reshape(10,6),columns=list('ABCDEF')) # TEST
    df['B'][0]=np.nan
    print(df)
    return df


def main():
    path,sep,accuracy,maxIterations = getInput()
    df = importDataFrame(path,sep)
    assert isinstance(df, pd.DataFrame)
    data = np.asarray(df) # ndarray instead of matrix because this is more convenient in the calculations
    normalizedData,convergenceTrail,R,S = cd.constand(data,accuracy,maxIterations)
    #print(data)
    print(normalizedData)
    return normalizedData,convergenceTrail

if __name__ == '__main__':
    sys.exit(main())
