""" Collection of handy tools used by various modules of the CONSTANd++ workflow. """

import numpy as np
from constandpp import fontsize, fontweight


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
	from pandas import Series
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
	if isinstance(intensities, np.ndarray):
		assert df.loc[:, intensityColumns].shape == intensities.shape
		df.loc[:, intensityColumns] = intensities
	elif isinstance(intensities, dict):
		for index in intensities.keys():
			df.loc[index, intensityColumns] = intensities[index]
	return df


def MA(x, y):
	logx = np.log2(x)
	logy = np.log2(y)
	A = (logx + logy) * 0.5
	M = logx - logy
	m = np.mean(M[np.isfinite(M)])
	v = np.var(M[np.isfinite(M)])
	return M, A, m, v


def scatterplot(x, y, title=None, xlab=None, ylab=None):
	import matplotlib.rcParams.update
	from matplotlib import pyplot as plt
	matplotlib.rcParams.update({'font.size': fontsize, 'font.weight': fontweight})
	f = plt.figure(figsize=(12, 9))
	plt.scatter(x, y)
	if title:
		plt.title(title)
	if xlab:
		plt.xlabel(xlab)
	if ylab:
		plt.ylabel(ylab)
	return f


def MAPlot(x, y, title=None):
	M, A, m, v = MA(x, y)
	if title is None:
		title = title('PD2.1 Intensities versus S/N values (scaled relatively within each row/peptide)')
	elif title == '':
		title = title('mean(M): ' + str(m) + '; var(M):' + str(v))
	else:
		title = title(title + '; mean(M): ' + str(m) + '; var(M):' + str(v))
	fig = scatterplot(M, A, title)
	# fig.xlabel('A')
	# fig.ylabel('M')
	# plt.xlim((-10.1,0))
	# plt.ylim((-10.1, 10.1))
	fig.show()
	return M, A, m, v
