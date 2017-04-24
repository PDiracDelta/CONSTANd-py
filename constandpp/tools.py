""" Collection of handy tools used by various modules of the CONSTANd++ workflow. """


def unnest(x):
	""" returns un-nested version of level 1 nested list x."""
	return [e for sublist in x for e in sublist]


def getIntensities(df, intensityColumns, indices=None):
	"""
	Extracts the (absolute) matrix of quantification values from the dataFrame as an ndarray.
	:param df:              pd.dataFrame or pd.Series   Pandas dataFrame/Series from which to extract the intensities
	:param intensityColumns:list						columns that contain the quantification values
	:param indices:         list                        indices of the entries for which to obtain the intensities
	:return intensities:    np.ndarray                  matrix of quantification values
	"""
	import numpy as np
	import pandas.Series as Series
	if isinstance(df, Series):  # this is a dataframe with only 1 entry: indexing [:, cols] doesnt work.
		return np.asarray(df.loc[intensityColumns])
	elif indices is None:
		return np.asarray(df.loc[:, intensityColumns])
	else:
		return np.asarray(df.loc[indices, intensityColumns])


def setIntensities(df, intensities, intensityColumns):
	"""
	Sets the quantification values of the dataFrame at the specified location equal to the array of given intensities,
	at the specified locations if a dict is provided instead of an array.
	:param df:              pd.dataFrame    input dataframe
	:param intensities:     np.ndarray      matrix with MS2 intensities
							dict            dict {index:[values]} with index and values of all df entries to be modified
	:param intensityColumns:list			columns that contain the quantification values
	:return df:             pd.dataFrame    output dataframe with updated intensities
	"""
	import numpy.ndarray as ndarray
	if isinstance(intensities, ndarray):
		assert df.loc[:, intensityColumns].shape == intensities.shape
		df.loc[:, intensityColumns] = intensities
	elif isinstance(intensities, dict):
		for index in intensities.keys():
			df.loc[index, intensityColumns] = intensities[index]
	return df

