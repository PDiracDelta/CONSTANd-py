""" Collection of handy tools used by various modules of the CONSTANd++ workflow. """

import numpy as np
from constandpp import fontsize, fontweight, figheight, figwidth
import matplotlib as mpl
from matplotlib import pyplot as plt


def unnest(x):
	""" returns un-nested version of level 1 nested list x."""
	return [e for sublist in x for e in sublist]


def partition(pred, iterable):
	""" Partition an iterable into a list of trues and list of falses according to some conditional function pred. """
	trues = []
	falses = []
	for item in iterable:
		if pred(item):
			trues.append(item)
		else:
			falses.append(item)
	return trues, falses


def getIntensities(df, quanColumns, indices=None):
	"""
	Extracts the (absolute) matrix of quantification values from the dataFrame as an ndarray.
	:param df:              pd.dataFrame or pd.Series   Pandas dataFrame/Series from which to extract the intensities
	:param quanColumns:list						columns that contain the quantification values
	:param indices:         list                        indices of the entries for which to obtain the intensities
	:return intensities:    np.ndarray                  matrix of quantification values
	"""
	import numpy as np
	from pandas import Series
	if isinstance(df, Series):  # this is a dataframe with only 1 entry: indexing [:, cols] doesnt work.
		return np.asarray(df.loc[quanColumns])
	elif indices is None:
		return np.asarray(df.loc[:, quanColumns])
	else:
		return np.asarray(df.loc[indices, quanColumns])


def setIntensities(df, intensities, quanColumns):
	"""
	Sets the quantification values of the dataFrame at the specified location equal to the array of given intensities,
	at the specified locations if a dict is provided instead of an array.
	:param df:              pd.dataFrame    input dataframe
	:param intensities:     np.ndarray      matrix with MS2 intensities
							dict            dict {index:[values]} with index and values of all df entries to be modified
	:param quanColumns:list			columns that contain the quantification values
	:return df:             pd.dataFrame    output dataframe with updated intensities
	"""
	if isinstance(intensities, np.ndarray):
		assert df.loc[:, quanColumns].shape == intensities.shape
		df.loc[:, quanColumns] = intensities
	elif isinstance(intensities, dict):
		for index in intensities.keys():
			df.loc[index, quanColumns] = intensities[index]
	return df


def MA(x, y):
	"""
	Returns for (x,y) the corresponding M and A values as defined for the MA (minus-additive) plot, as well as the
	average and variance of the M values.
	:param x:	list	input values x
	:param y:	list	input values y
	:return M:	list	logx - logy
	:return A:	list	(logx + logy) * 0.5
	:return m:	float64	mean(M)
	:return v:	float64	var(M)
	"""
	logx = np.log2(x)
	logy = np.log2(y)
	A = (logx + logy) * 0.5
	M = logx - logy
	m = np.mean(M[np.isfinite(M)])
	v = np.var(M[np.isfinite(M)])
	return M, A, m, v


def scatterplot(x, y, title=None, xlab=None, ylab=None):
	"""
	Return a scatterplot for the given data, with the constandpp-specific plot settings (fonts, figsizes, ...).
	:param x:		list		x data
	:param y:		list		y data
	:param title:	str			scatterplot title
	:param xlab:	str			label for the x-axis
	:param ylab:	str			label for the y-axis
	:return f:		plt.figure	scatterplot as a matplotlib figure object
	"""
	mpl.use('GTK3Agg')
	mpl.rcParams.update({'font.size': fontsize, 'font.weight': fontweight})
	
	f = plt.figure(figsize=(figwidth, figheight))
	plt.scatter(x, y)
	if title:
		plt.title(title)
	if xlab:
		plt.xlabel(xlab)
	if ylab:
		plt.ylabel(ylab)
	return f


def MAPlot(x, y, title=None):
	"""
	Return an MA plot for the given data, with the constandpp-specific plot settings (fonts, figsizes, ...) by calling
	MA() and then scatterplot().
	:param x:		list	x data
	:param y:		list	y data
	:param title:	str		MA plot title
	:return M:		list	logx - logy
	:return A:		list	(logx + logy) * 0.5
	:return m:		float64	mean(M)
	:return v:		float64	var(M)
	"""
	M, A, m, v = MA(x, y)
	if title is None:
		title = 'PD2.1 Intensities versus S/N values (scaled relatively within each row/peptide)'
	elif title == '':
		title = 'mean(M): ' + str(m) + '; var(M):' + str(v)
	else:
		title = title + '\nmean(M): ' + str(m) + '; var(M):' + str(v)
	fig = scatterplot(A, M, title=title, xlab='A', ylab='M')
	# plt.xlim((-10.1,0))
	# plt.ylim((-10.1, 10.1))
	fig.show()
	return M, A, m, v
