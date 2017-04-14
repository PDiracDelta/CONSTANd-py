#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Python implementation of mass spectrometer protein data analysis using the CONSTANd_RAS algorithm.
"""

__author__ = "Joris Van Houtven"
__copyright__ = "Copyright ?YEAR, VITO"
__credits__ = ["Joris Van Houtven", "Dirk Valkenborg"]
# __license__ = "GPL"
# __version__ = "0.1"
__maintainer__ = "Joris Van Houtven"
__email__ = "vanhoutvenjoris@gmail.com"
__status__ = "Development"

import sys, logging, datetime
from dataIO import *

fontsize = 30
fontweight = 'normal'


def performanceTest():  # remove for production # TEST
	""" Use this development method to test the performance of the CONSTANd algorithm. """
	from constand import constand
	from getInput import getProcessingInput
	from time import time
	
	t = []
	for i in range(100):
		params = getProcessingInput()
		df = pd.DataFrame(np.random.uniform(low=10 ** 3, high=10 ** 5, size=(2 * 10 ** 3, 6)), columns=list('ABCDEF'))
		start = time()
		constand(np.asarray(df), 1e-2, 50)
		stop = time()
		t.append((stop - start))
	print("average runtime: " + str(np.mean(t)))


def isotopicImpuritiesTest():  # TEST
	from getInput import getProcessingInput
	from constand import constand
	from processing import getIntensities
	## test if isotopic correction is necessary:
	params = getProcessingInput()
	# get the dataframe
	df = importDataFrame(params['files_in'], delim=params['delim_in'], header=params['header_in'])
	correctedIntensities = getIntensities(df)
	normalizedIntensities, convergenceTrail, R, S = constand(correctedIntensities, params['accuracy'],
															 params['maxIterations'])
	# exportData(normalizedIntensities, 'txt', path_out=params['path_out'],
	#            filename=params['jobname'] + '_normalizedIntensities', delim_out=params['delim_out'])
	# test "impure data"
	correctedIntensities_impure = correctedIntensities
	spillover = correctedIntensities_impure[0, :] * 0.1
	correctedIntensities_impure[0, :] -= spillover
	correctedIntensities_impure[1, :] += spillover
	normalizedIntensities_impure, convergenceTrail_i, R_i, S_i = constand(correctedIntensities_impure,
																		  params['accuracy'],
																		  params['maxIterations'])
	diff = abs(normalizedIntensities - normalizedIntensities_impure)
	print(np.allclose(normalizedIntensities, normalizedIntensities_impure, atol=1e-3, equal_nan=True))
	print(np.nanmean(np.nanmean(diff[:, 0:1], 1)))
	print(max(np.amax(diff, 1)))


# False tot en met 1e-3 --> fouten van > 0.1%


def isotopicCorrectionsTest():  # TEST
	from processing import getIntensities
	from constand import constand
	filepath1 = '../data/COON data/PSMs/BR1_a.txt'
	filepathiso = '../data/COON data/PSMs/BR1_e_ISO.txt'
	ics = unnest([["126", "127N", "127C", "128C"], ["129N", "129C", "130C", "131"]])
	quan = np.asmatrix(
		constand(getIntensities(importDataFrame(filepath1, delim='\t'), intensityColumns=ics), accuracy=1e-5,
				 maxIterations=50)[0])
	# without constand
	# quan = np.asmatrix(getIntensities(importDataFrame(filepath1, delim='\t'), intensityColumns=ics))
	quanIso = np.asmatrix(
		constand(getIntensities(importDataFrame(filepathiso, delim='\t'), intensityColumns=ics), accuracy=1e-5,
				 maxIterations=50)[0])
	# quandf = pd.DataFrame(quan)
	# quanIsodf = pd.DataFrame(quanIso)
	# for col in range(len(ics)):
	#	quanIsodf.loc[:, col] = 1 / 8 * np.asarray(quanIsodf.loc[:, col]) / np.nanmean(np.asarray(quanIsodf), 1)
	# relSNs[:, col] = 1 / 8 * df2.loc[:, ics[col]] / np.nanmean(df2.loc[:, ics], 1)
	# quanIso = np.asmatrix(quanIsodf)
	diff = abs(quan - quanIso)
	MAPlot(quan.flatten(), quanIso.flatten(), title='')
	# RDHPlot(quan.flatten(), quanIso.flatten())
	# print(np.allclose(quan, quanIso, atol=1e-3, equal_nan=True))
	print("mean over all values")
	print(np.nanmean(np.nanmean(diff[:, 0:7], 1)))
	print("max difference")
	print(np.nanmax(np.nanmax(diff, 1)))
	print("median over all values")


# print(np.nanmean(np.nanmedian(diff[:, 0:7], 1)))
# M=np.eye(6); M[0,0]=0.9; M[0,1]=0.1; b=np.asarray(range(6)); c=np.asarray(range(6))+5
# print(int_out) above should be equal to:
# [np.linalg.solve(M, b) ; np.linalg.solve(M, c)]


def MS2IntensityDoesntMatter(df):
	from processing import getIntensities
	from constand import constand
	I = getIntensities(df)
	r1 = constand(I, 1e-5, 50)
	I[
		0] *= 1e9  # this is BIG. MS2 intensity doesnt reach beyond 1e9 so if one value has magnitudeOrder 1 it's still OK.
	r2 = constand(I, 1e-5, 50)
	print(np.allclose(r1[0], r2[0], equal_nan=True))
	diff = r1[0] - r2[0]
	maxdiff = max(np.amax(diff, 1))
	print(maxdiff)


def testDataComplementarity(df):
	scannrs_init = set(df.groupby('First Scan').groups.keys())
	main(testing=False, writeToDisk=True)
	# SANITY CHECK if df + removedData scan numbers = total scan numbers.
	scannrs_final = set(df.groupby('First Scan').groups.keys())
	##### THIS IS OUTDATED SINCE COMMIT b98041f
	removedDataLoaded = pickle.load(open('../data/MB_result_removedData', 'rb'))
	for value in removedDataLoaded.values():
		scannrs_final = scannrs_final.union(set(value['First Scan']))
	print(scannrs_final == scannrs_init)


def MA(x, y):
	logx = np.log2(x)
	logy = np.log2(y)
	A = (logx + logy) * 0.5
	M = logx - logy
	m = np.mean(M[np.isfinite(M)])
	v = np.var(M[np.isfinite(M)])
	return M, A, m, v


def MAPlot(x, y, title=None):
	import matplotlib
	from matplotlib import pyplot as plt
	matplotlib.rcParams.update({'font.size': fontsize, 'font.weight': fontweight})
	plt.figure(figsize=(16, 12))
	M, A, m, v = MA(x, y)
	plt.scatter(A, M)
	if title is None:
		plt.title('PD2.1 Intensities versus S/N values (scaled relatively within each row/peptide)')
	elif title == '':
		plt.title('mean(M): ' + str(m) + '; var(M):' + str(v))
	else:
		plt.title(title + '; mean(M): ' + str(m) + '; var(M):' + str(v))
	plt.xlabel('A')
	plt.ylabel('M')
	# plt.xlim((-10.1,0))
	# plt.ylim((-10.1, 10.1))
	plt.show()
	return M, A, m, v


def compareIntensitySN(df1, df2, title=None):
	from getInput import getProcessingInput
	from processingFlow import processDf
	from processing import getIntensities
	
	filepath1 = '../data/COON data/PSMs/BR1_e_ISO.txt'
	filepath2 = '../data/COON data/PSMs/BR1_f_ISO_SN.txt'
	intensityColumns = ["126", "127N", "127C", "128C", "129N", "129C", "130C", "131"]
	pickleFileName = 'job/compareIntensitySNProcessingResults'
	constandnorm = True
	alsoprocess = False
	if constandnorm:
		if alsoprocess and os.path.exists(pickleFileName):
			processingResults = pickle.load(open(pickleFileName, 'rb'))
		else:
			params = getProcessingInput('job/processingConfig.ini')
			dfs = []
			if df1 is None and df2 is None:
				for filepath in [filepath1, filepath2]:
					dfs.append(importDataFrame(filepath, delim=params['delim_in'], header=params['header_in']))
			else:
				dfs = [df1, df2]
			if alsoprocess:
				processingResults = [processDf(df, params, writeToDisk=False) for df in dfs]
				pickle.dump(processingResults, open(pickleFileName, 'wb'))
		if alsoprocess:
			relIntensities = getIntensities(processingResults[0][0], intensityColumns=intensityColumns)
			relSNs = getIntensities(processingResults[1][0], intensityColumns=intensityColumns)
		else:
			from constand import constand
			relIntensities, __, __, __ = constand(getIntensities(dfs[0], intensityColumns=intensityColumns), 1e-5, 50)
			relSNs, __, __, __ = constand(getIntensities(dfs[1], intensityColumns=intensityColumns), 1e-5, 50)
	else:
		if df1 is None and df2 is None:
			df1 = importDataFrame(filepath1, delim='\t', header=0)
			df2 = importDataFrame(filepath2, delim='\t', header=0)
		intensityColumns = ["126", "127N", "127C", "128C", "129N", "129C", "130C", "131"]
		relIntensities = np.empty((len(df1.index), 8), dtype='float')
		relSNs = np.empty((len(df2.index), 8), dtype='float')
		for col in range(len(intensityColumns)):
			relIntensities[:, col] = 1 / 8 * df1.loc[:, intensityColumns[col]] / np.nanmean(
				df1.loc[:, intensityColumns], 1)
			relSNs[:, col] = 1 / 8 * df2.loc[:, intensityColumns[col]] / np.nanmean(df2.loc[:, intensityColumns], 1)
	diff = abs(relIntensities - relSNs)
	print(np.allclose(relIntensities, relSNs, atol=1e-3, equal_nan=True))
	print("mean over all values")
	print(np.nanmean(np.nanmean(diff[:, 0:7], 1)))
	print("median over all values")
	print(np.nanmean(np.nanmedian(diff[:, 0:7], 1)))
	print("max difference")
	print(np.nanmax(np.nanmax(diff, 1)))
	
	return MAPlot(relIntensities.reshape(relIntensities.size, 1), relSNs.reshape(relSNs.size, 1), title)


def compareICmethods():
	PDdfFile = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/2016-12-26 10:56:10.646919_COON/BR1_output_processing/BR1_normalizedIntensities.tsv'
	CdfFile = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/2016-12-26 10:50:42.462546_COON_noISO/BR1_output_processing/BR1_normalizedIntensities.tsv'
	from dataIO import importDataFrame
	from main import compareIntensitySN
	from matplotlib import pyplot as plt
	Cdf = importDataFrame(CdfFile)
	PDdf = importDataFrame(PDdfFile)
	PDdf.columns = ["126", "127N", "127C", "128C", "129N", "129C", "130C", "131"]
	Cdf.columns = ["126", "127N", "127C", "128C", "129N", "129C", "130C", "131"]
	M, A, m, v = compareIntensitySN(PDdf, Cdf, title='')
	print("MAD: " + str(np.nanmean(np.abs(M))))
	Mfinite = M[np.isfinite(M)]
	# np.digitize(M, np.linspace(min(M),max(M),20))
	nbins = 21
	hist, bins = np.histogram(Mfinite, bins=nbins, range=(-5, 5))
	binwidth = bins[1] - bins[0]
	plt.title('')
	plt.bar(bins[0:-1], hist, width=binwidth)
	plt.xticks(np.arange(11) - 5)
	plt.xlim(-5, 5)
	plt.xlabel('M')
	plt.ylabel('amount of peptides')
	plt.show()


def abundancesPCAHCD():
	from analysisFlow import getPCA, getHC
	from analysis import getAllExperimentsIntensitiesPerCommonPeptide
	from report import getPCAPlot, getHCDendrogram
	import json
	datapath = '../data/COON data/peptidegroups/'
	wrapperpath = '../jobs/2016-12-12 22:37:48.458146_COON/'
	schema = json.loads(
		'{"BR4": {"channelNamesPerCondition": [["127C"], ["128C"], ["126"], ["131"], ["129C"], ["130C"], ["129N"], ["127N"]], "wrapper": "BR4_wrapper.tsv", "isotopicCorrection_matrix": null, "channelAliasesPerCondition": [["4_kidney"], ["4_cerebellum"], ["4_muscle"], ["4_cerebrum"], ["4_lung"], ["4_liver"], ["4_heart"], ["4_spleen"]], "data": "BR4_BR4_e_ISO.txt", "config": "/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/2016-12-12 22:37:48.458146_COON/BR4_coonProcessingConfig.ini"}, "BR3": {"channelNamesPerCondition": [["131"], ["130C"], ["128C"], ["129C"], ["126"], ["127C"], ["129N"], ["127N"]], "wrapper": "BR3_wrapper.tsv", "isotopicCorrection_matrix": null, "channelAliasesPerCondition": [["3_kidney"], ["3_cerebellum"], ["3_muscle"], ["3_cerebrum"], ["3_lung"], ["3_liver"], ["3_heart"], ["3_spleen"]], "data": "BR3_BR3_e_ISO.txt", "config": "/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/2016-12-12 22:37:48.458146_COON/BR3_coonProcessingConfig.ini"}, "BR1": {"channelNamesPerCondition": [["126"], ["127C"], ["127N"], ["128C"], ["129C"], ["129N"], ["130C"], ["131"]], "wrapper": "BR1_wrapper.tsv", "isotopicCorrection_matrix": null, "channelAliasesPerCondition": [["1_kidney"], ["1_cerebellum"], ["1_muscle"], ["1_cerebrum"], ["1_lung"], ["1_liver"], ["1_heart"], ["1_spleen"]], "data": "BR1_BR1_e_ISO.txt", "config": "/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/2016-12-12 22:37:48.458146_COON/BR1_coonProcessingConfig.ini"}, "BR2": {"channelNamesPerCondition": [["127N"], ["126"], ["128C"], ["130C"], ["129N"], ["127C"], ["129C"], ["131"]], "wrapper": "BR2_wrapper.tsv", "isotopicCorrection_matrix": null, "channelAliasesPerCondition": [["2_kidney"], ["2_cerebellum"], ["2_muscle"], ["2_cerebrum"], ["2_lung"], ["2_liver"], ["2_heart"], ["2_spleen"]], "data": "BR2_BR2_e_ISO.txt", "config": "/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/2016-12-12 22:37:48.458146_COON/BR2_coonProcessingConfig.ini"}}')
	eNames = list(schema.keys())
	datafilenames = {}
	wrapperfilenames = {}
	for e in eNames:
		# datafilenames[e] = 'P_'+e+'_h_ISO_SN_norm.txt'
		datafilenames[e] = 'P_' + e + '_g_ISO_norm.txt'
		wrapperfilenames[e] = e + '_wrapper.tsv'
	dfs = {}
	wrappers = {}
	for eName in eNames:
		df = importDataFrame(datapath + datafilenames[eName], delim='\t')
		intensityCols = unnest([["127C"], ["128C"], ["126"], ["131"], ["129C"], ["130C"], ["129N"], ["127N"]])
		oldIntensityCols = list(df.columns.values)
		newIntensityCols = oldIntensityCols.copy()
		for col in intensityCols:
			newIntensityCols[
				oldIntensityCols.index('Abundances (Normalized): F' + eName[-1] + ': ' + col + ', Sample')] = col
		df.columns = newIntensityCols
		wrappers[eName] = importWrapper(wrapperpath + wrapperfilenames[eName])
		df.columns = applyWrapper(df.columns, wrapper=wrappers[eName])
		fixFixableFormatMistakes(df)
		dfs[eName] = df
	
	AEIPCP = getAllExperimentsIntensitiesPerCommonPeptide(dfs, schema)[0].astype(np.float)
	doConstand = True
	if doConstand:
		from constand import constand
		AEIPCPcorrected = constand(AEIPCP, 1e-5, 50)[0]
	else:
		colTotals = np.nansum(AEIPCP, axis=0)
		multipliers = np.max(colTotals) / colTotals
		AEIPCPcorrected = AEIPCP * multipliers
	pca = getPCA(AEIPCPcorrected, 2)
	hcd = getHC(AEIPCPcorrected)
	PCAPlot = getPCAPlot(pca, schema)
	HCDendrogram = getHCDendrogram(hcd, schema)


def compareAbundancesIntSN():
	abundancesIntFile = '../jobs/2016-12-20 19:03:17.450992_COON_abundances/P_BR1_output_processing/P_BR1_dataFrame.tsv'
	abundancesSNFile = '../jobs/2016-12-20 19:02:48.246369_COON_SN_abundances/P_BR1_output_processing/P_BR1_dataFrame.tsv'
	abundancesInt = importDataFrame(abundancesIntFile, delim='\t')
	abundancesSN = importDataFrame(abundancesSNFile, delim='\t')
	intensityColumns = ["126", "127N", "127C", "128C", "129N", "129C", "130C", "131"]
	abundancesInt.columns = intensityColumns
	abundancesSN.columns = intensityColumns
	compareIntensitySN(abundancesInt, abundancesSN,
					   title='')  # title="Intensities versus S/N values (normalized abundances)")


def boxPlot(x, labels=None, ylab=None):
	import matplotlib
	from matplotlib import pyplot as plt
	matplotlib.rcParams.update({'font.size': fontsize, 'font.weight': fontweight})
	plt.boxplot(x, whis=[5, 95], showmeans=True, labels=labels)
	plt.ylabel(ylab)
	plt.show()


def intraInterMAPlots():
	# calculate all intra- and inter-MA plots for the Max data set (Intensities only)
	# on the PEPTIDE LEVEL
	constandFileName = '../jobs/2016-12-22 13:07:48.663095_MAX_SN_noTAM/output_analysis/MAX_SN_noTAM_CommonPeptideIntensities.tsv'
	PDFileName = '../jobs/2016-12-24 12:09:28.341520_MAX_SN_abundances_noTAM/output_analysis/MAX_SN_abundances_noTAM_CommonPeptideIntensities.tsv'
	rawFileName = '../jobs/2016-12-22 12:56:54.697920_MAX_SN_nonormnoconstand_noTAM/output_analysis/MAX_SN_nonormnoconstand_noTAM_CommonPeptideIntensities.tsv'
	cdf = importDataFrame(constandFileName, delim='\t')
	pddf = importDataFrame(PDFileName, delim='\t')
	rdf = importDataFrame(rawFileName, delim='\t')
	
	compareDIFFERENTconditions = True
	compareIDENTICALconditions = True
	diffboxplots = False  # THIS IS USELESS THE DIMENSIONS DONT MATCH
	if compareDIFFERENTconditions:
		D_intra_cmeans, D_intra_pdmeans, D_intra_rmeans = [], [], []
		D_intra_cvars, D_intra_pdvars, D_intra_rvars = [], [], []
		experiments = [[3, 4], [5, 6], [1, 2]]
		# INTRA
		for e in experiments:  # for each experiment
			for i in e:  # for each sample of BM in this experiment
				for j in e:  # compare with the samples of PM in the same experiment
					cm, cv = MA(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					D_intra_cmeans.append(cm)
					D_intra_cvars.append(cv)
					pdm, pdv = MA(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'PM' + str(j)])[
							   2:4]  # [2:4] only mean and var
					D_intra_pdmeans.append(pdm)
					D_intra_pdvars.append(pdv)
					rm, rv = MA(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					D_intra_rmeans.append(rm)
					D_intra_rvars.append(rv)
		print("+++INTRA: \n"
			  "CONSTANd means: MAD: " + str(np.mean(np.abs(D_intra_cmeans))) + "; values: " + str(D_intra_cmeans) + "\n"
																													"PD2.1 means: MAD: " + str(
			np.mean(np.abs(D_intra_pdmeans))) + "; values: " + str(D_intra_pdmeans) + "\n"
																					  "raw means: MAD: " + str(
			np.mean(np.abs(D_intra_rmeans))) + "; values: " + str(D_intra_rmeans) + "\n"
																					"CONSTANd vars: average: " + str(
			np.mean(D_intra_cvars)) + "; values: " + str(D_intra_cvars) + "\n"
																		  "PD2.1 vars: average: " + str(
			np.mean(D_intra_pdvars)) + "; values: " + str(D_intra_pdvars) + "\n"
																			"raw vars: average: " + str(
			np.mean(D_intra_rvars)) + "; values: " + str(D_intra_rvars) + "\n")
		# MAPlot(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'PM' + str(j)])
		# MAPlot(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'PM' + str(j)])
		# MAPlot(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'PM' + str(j)])
		boxPlot([D_intra_rvars, D_intra_cvars, D_intra_pdvars], labels=['raw', 'CONSTANd', 'PD2.1'], ylab='variance')
		
		# INTER
		# experiments = [[3,4],[5,6],[1,2]]
		D_inter_cmeans, D_inter_pdmeans, D_inter_rmeans = [], [], []
		D_inter_cvars, D_inter_pdvars, D_inter_rvars = [], [], []
		for l in range(len(experiments)):
			notE = unnest(experiments[:l] + experiments[l + 1:])
			e = experiments[l]
			for i in e:
				for j in notE:
					cm, cv = MA(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					D_inter_cmeans.append(cm)
					D_inter_cvars.append(cv)
					pdm, pdv = MA(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'PM' + str(j)])[
							   2:4]  # [2:4] only mean and var
					D_inter_pdmeans.append(pdm)
					D_inter_pdvars.append(pdv)
					rm, rv = MA(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					D_inter_rmeans.append(rm)
					D_inter_rvars.append(rv)
		print("+++INTER: \n"
			  "CONSTANd means: MAD: " + str(np.mean(np.abs(D_inter_cmeans))) + "; values: " + str(D_inter_cmeans) + "\n"
																													"PD2.1 means: MAD: " + str(
			np.mean(np.abs(D_inter_pdmeans))) + "; values: " + str(D_inter_pdmeans) + "\n"
																					  "raw means: MAD: " + str(
			np.mean(np.abs(D_inter_rmeans))) + "; values: " + str(D_inter_rmeans) + "\n"
																					"CONSTANd vars: average: " + str(
			np.mean(D_inter_cvars)) + "; values: " + str(D_inter_cvars) + "\n"
																		  "PD2.1 vars: average: " + str(
			np.mean(D_inter_pdvars)) + "; values: " + str(D_inter_pdvars) + "\n"
																			"raw vars: average: " + str(
			np.mean(D_inter_rvars)) + "; values: " + str(D_inter_rvars) + "\n")
		# plot one graph for each (the last one):
		# MAPlot(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'PM' + str(j)])
		# MAPlot(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'PM' + str(j)])
		# MAPlot(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'PM' + str(j)])
		boxPlot([D_inter_rvars, D_inter_cvars, D_inter_pdvars], labels=['raw', 'CONSTANd', 'PD2.1'], ylab='variance')
	
	if compareIDENTICALconditions:
		I_intra_cmeans, I_intra_pdmeans, I_intra_rmeans = [], [], []
		I_intra_cvars, I_intra_pdvars, I_intra_rvars = [], [], []
		experiments = [[3, 4], [5, 6], [1, 2]]
		# INTRA
		for e in experiments:  # for each experiment
			for l in range(len(e)):  # for each sample in (and each condition) this experiment
				notI = e[l + 1:]
				i = e[l]
				if notI:  # not empty
					for j in notI:  # compare with the other sample(s) (that are of the same condition) in this experiment
						cm, cv = MA(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'BM' + str(j)])[
								 2:4]  # [2:4] only mean and var
						I_intra_cmeans.append(cm)
						I_intra_cvars.append(cv)
						cm, cv = MA(cdf.loc[:, 'PM' + str(i)], cdf.loc[:, 'PM' + str(j)])[
								 2:4]  # [2:4] only mean and var
						I_intra_cmeans.append(cm)
						I_intra_cvars.append(cv)
						
						pdm, pdv = MA(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'BM' + str(j)])[
								   2:4]  # [2:4] only mean and var
						I_intra_pdmeans.append(pdm)
						I_intra_pdvars.append(pdv)
						pdm, pdv = MA(pddf.loc[:, 'PM' + str(i)], pddf.loc[:, 'PM' + str(j)])[
								   2:4]  # [2:4] only mean and var
						I_intra_pdmeans.append(pdm)
						I_intra_pdvars.append(pdv)
						
						rm, rv = MA(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'BM' + str(j)])[
								 2:4]  # [2:4] only mean and var
						I_intra_rmeans.append(rm)
						I_intra_rvars.append(rv)
						rm, rv = MA(rdf.loc[:, 'PM' + str(i)], rdf.loc[:, 'PM' + str(j)])[
								 2:4]  # [2:4] only mean and var
						I_intra_rmeans.append(rm)
						I_intra_rvars.append(rv)
		print("+++INTRA: \n"
			  "CONSTANd means: MAD: " + str(np.mean(np.abs(I_intra_cmeans))) + "; values: " + str(I_intra_cmeans) + "\n"
																													"PD2.1 means: MAD: " + str(
			np.mean(np.abs(I_intra_pdmeans))) + "; values: " + str(I_intra_pdmeans) + "\n"
																					  "raw means: MAD: " + str(
			np.mean(np.abs(I_intra_rmeans))) + "; values: " + str(I_intra_rmeans) + "\n"
																					"CONSTANd vars: average: " + str(
			np.mean(I_intra_cvars)) + "; values: " + str(I_intra_cvars) + "\n"
																		  "PD2.1 vars: average: " + str(
			np.mean(I_intra_pdvars)) + "; values: " + str(I_intra_pdvars) + "\n"
																			"raw vars: average: " + str(
			np.mean(I_intra_rvars)) + "; values: " + str(I_intra_rvars) + "\n")
		# MAPlot(cdf.loc[:, 'BM' + str(1)], cdf.loc[:, 'BM' + str(2)])
		# MAPlot(pddf.loc[:, 'BM' + str(1)], pddf.loc[:, 'BM' + str(2)])
		# MAPlot(rdf.loc[:, 'BM' + str(1)], rdf.loc[:, 'BM' + str(2)])
		boxPlot([I_intra_rvars, I_intra_cvars, I_intra_pdvars], labels=['raw', 'CONSTANd', 'PD2.1'], ylab='variance')
		
		# INTER
		# experiments = [[3,4],[5,6],[1,2]]
		I_inter_cmeans, I_inter_pdmeans, I_inter_rmeans = [], [], []
		I_inter_cvars, I_inter_pdvars, I_inter_rvars = [], [], []
		for l in range(len(experiments)):
			notE = unnest(experiments[l + 1:])
			e = experiments[l]
			for i in e:
				for j in notE:
					cm, cv = MA(cdf.loc[:, 'BM' + str(i)], cdf.loc[:, 'BM' + str(j)])[2:4]  # [2:4] only mean and var
					I_inter_cmeans.append(cm)
					I_inter_cvars.append(cv)
					cm, cv = MA(cdf.loc[:, 'PM' + str(i)], cdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					I_inter_cmeans.append(cm)
					I_inter_cvars.append(cv)
					
					pdm, pdv = MA(pddf.loc[:, 'BM' + str(i)], pddf.loc[:, 'BM' + str(j)])[
							   2:4]  # [2:4] only mean and var
					I_inter_pdmeans.append(pdm)
					I_inter_pdvars.append(pdv)
					pdm, pdv = MA(pddf.loc[:, 'PM' + str(i)], pddf.loc[:, 'PM' + str(j)])[
							   2:4]  # [2:4] only mean and var
					I_inter_pdmeans.append(pdm)
					I_inter_pdvars.append(pdv)
					
					rm, rv = MA(rdf.loc[:, 'BM' + str(i)], rdf.loc[:, 'BM' + str(j)])[2:4]  # [2:4] only mean and var
					I_inter_rmeans.append(rm)
					I_inter_rvars.append(rv)
					rm, rv = MA(rdf.loc[:, 'PM' + str(i)], rdf.loc[:, 'PM' + str(j)])[2:4]  # [2:4] only mean and var
					I_inter_rmeans.append(rm)
					I_inter_rvars.append(rv)
		print("+++INTER: \n"
			  "CONSTANd means: MAD: " + str(np.mean(np.abs(I_inter_cmeans))) + "; values: " + str(I_inter_cmeans) + "\n"
																													"PD2.1 means: MAD: " + str(
			np.mean(np.abs(I_inter_pdmeans))) + "; values: " + str(I_inter_pdmeans) + "\n"
																					  "raw means: MAD: " + str(
			np.mean(np.abs(I_inter_rmeans))) + "; values: " + str(I_inter_rmeans) + "\n"
																					"CONSTANd vars: MAD: " + str(
			np.mean(I_inter_cvars)) + "; values: " + str(I_inter_cvars) + "\n"
																		  "PD2.1 vars: MAD: " + str(
			np.mean(I_inter_pdvars)) + "; values: " + str(I_inter_pdvars) + "\n"
																			"raw vars: MAD: " + str(
			np.mean(I_inter_rvars)) + "; values: " + str(I_inter_rvars) + "\n")
		# plot one graph for each (the last one):
		# MAPlot(cdf.loc[:, 'BM' + str(3)], cdf.loc[:, 'BM' + str(6)])
		# MAPlot(pddf.loc[:, 'BM' + str(3)], pddf.loc[:, 'BM' + str(6)])
		# MAPlot(rdf.loc[:, 'BM' + str(3)], rdf.loc[:, 'BM' + str(6)])
		boxPlot([I_inter_rvars, I_inter_cvars, I_inter_pdvars], labels=['raw', 'CONSTANd', 'PD2.1'], ylab='variance')
	
	if diffboxplots:
		# THIS IS USELESS THE DIMENSIONS DONT MATCH
		# transform all to arrays
		allLists = [I_inter_rvars, I_intra_rvars, I_inter_cvars, I_intra_cvars, I_inter_pdvars, I_intra_pdvars,
					D_inter_rvars, D_intra_rvars, D_inter_cvars, D_intra_cvars, D_inter_pdvars, D_intra_pdvars]
		I_inter_rvars, I_intra_rvars, I_inter_cvars, I_intra_cvars, I_inter_pdvars, I_intra_pdvars, D_inter_rvars, D_intra_rvars, D_inter_cvars, D_intra_cvars, D_inter_pdvars, D_intra_pdvars = (
		np.array(x) for x in allLists)
		# amke boxplots of differences
		boxPlot([I_inter_rvars - I_intra_rvars, I_inter_cvars - I_intra_cvars, I_inter_pdvars - I_intra_pdvars],
				labels=['raw', 'CONSTANd', 'PD2.1'], ylab=r'\Delta V')
		boxPlot([D_inter_rvars - D_intra_rvars, D_inter_cvars - D_intra_cvars, D_inter_pdvars - D_intra_pdvars],
				labels=['raw', 'CONSTANd', 'PD2.1'], ylab=r'\Delta V')


def RDHPlot(x, y, quantity=None):
	import matplotlib
	import matplotlib.pyplot as plt
	matplotlib.rcParams.update({'font.size': fontsize, 'font.weight': fontweight})
	plt.figure(figsize=(12, 9))
	x = np.array(x)
	y = np.array(y)
	# relative difference histogram
	diff = np.abs(x - y)
	reldiff = diff / np.maximum(x, y)
	# finite = diff[np.isfinite(M)]
	# np.digitize(M, np.linspace(min(M),max(M),20))
	hist, bins = np.histogram(reldiff, bins=max(10, np.ceil(max(reldiff) - min(reldiff))),
							  range=(0, np.ceil(max(reldiff) * 10) / 10))
	plt.ylabel('number of proteins')
	if quantity is None:
		plt.xlabel('relative difference')
	else:
		plt.xlabel('relative difference in ' + quantity)
	binsize = (max(bins) - min(bins)) / (len(bins) - 1)
	plt.bar(bins[0:-1], hist, width=binsize)
	plt.xticks(bins)
	plt.xlim(0, min(max(bins) + binsize, 1))
	plt.show()
	meandiff = np.nanmean(reldiff)
	maxdiff = np.max(reldiff)
	return meandiff, maxdiff


def compareDEAresults():
	df1file = '../jobs/2016-12-26 15:44:27.523732_MAX_noTAM/output_analysis/MAX_noTAM_results_minimal.tsv'
	df2file = '../jobs/2016-12-22 13:07:48.663095_MAX_SN_noTAM/output_analysis/MAX_SN_noTAM_results_minimal.tsv'
	# df1file = '../jobs/2016-12-30 14:46:02.689189_MAX_abundances_noTAM/output_analysis/MAX_abundances_noTAM_results_minimal.tsv'
	# df2file = '../jobs/2016-12-24 12:09:28.341520_MAX_SN_abundances_noTAM/output_analysis/MAX_SN_abundances_noTAM_results_minimal.tsv'
	df1sorted = importDataFrame(df1file, delim='\t').sort('protein')
	df2sorted = importDataFrame(df2file, delim='\t').sort('protein')
	assert set(df1sorted['protein']) == set(
		df2sorted['protein'])  # else you cant compare without first matching each protein
	# fcM, fcA, fcm, fcv =
	fmean, fmax = RDHPlot(np.power(2, df1sorted.loc[:, 'log2 fold change c1/c2']),
						  np.power(2, df2sorted.loc[:, 'log2 fold change c1/c2']), quantity='log2 fold change')
	print("fc mean: " + str(fmean) + "\nmax: " + str(fmax))
	# pM, pA, pm, pv =
	pmean, pmax = RDHPlot(df1sorted.loc[:, 'adjusted p-value'], df2sorted.loc[:, 'adjusted p-value'],
						  quantity='adjusted p-value')
	print("p mean: " + str(pmean) + "\nmax: " + str(pmax))


def dataSuitabilityMA():
	# rawfile = '../jobs/2016-12-20 14:39:09.476567_COON_SN_nonormnoconstand/output_analysis/COON_SN_nonormnoconstand_CommonPeptideIntensities.tsv'
	rawfile = '../jobs/2016-12-21 16:07:19.300450_MAX_SN_nonormnoconstand/output_analysis/MAX_SN_nonormnoconstand_CommonPeptideIntensities.tsv'
	rawSNvalues = importDataFrame(rawfile, delim='\t', header=None)
	if 'MAX' in rawfile:
		plexity = 6
	elif 'COON' in rawfile:
		plexity = 8
	# i = 3; j = 1
	# MAPlot(rawSNvalues.loc[:, i], rawSNvalues.loc[:, j], title='')
	# i=3; j=1
	# MAPlot(rawSNvalues.loc[:, i], rawSNvalues.loc[:, j], title='')
	# MAX
	# i=4; j=0
	i = 2;
	j = 6
	MAPlot(rawSNvalues.loc[:, i], rawSNvalues.loc[:, j], title='')
	for i in range(plexity):
		for j in range(i):
			MAPlot(rawSNvalues.loc[:, i], rawSNvalues.loc[:, j], title=str(i + 1) + " versus " + str(j + 1))


def devStuff(df, params):  # TEST
	# performanceTest()
	# isotopicCorrectionsTest()
	# MS2IntensityDoesntMatter(df)
	# testDataComplementarity(df)
	# compareIntensitySN(None, None, title='')
	# abundancesPCAHCD()
	# compareICmethods()
	compareAbundancesIntSN()
	# intraInterMAPlots()
	# compareDEAresults()
	# dataSuitabilityMA()
	pass


def main(jobConfigFilePath, doProcessing, doAnalysis, doReport, writeToDisk, testing):
	"""
	For now this is just stuff for debugging and testing. Later:
	Contains and explicits the workflow of the program. Using the booleans doProcessing, doAnalysis and writeToDisk one
	can control	which parts of the workflow to perform.
	"""
	# todo proper docu
	from getInput import getProcessingInput, getJobInput
	from processingFlow import processDf
	from analysisFlow import analyzeProcessingResult
	from reportFlow import generateReport
	from time import time
	from web.web import DB_setJobReportRelPaths
	
	logFilePath = os.path.abspath(os.path.join(jobConfigFilePath, os.path.join(os.pardir, 'log.txt')))
	logging.basicConfig(filename=logFilePath, level=logging.INFO)
	start = time()
	jobParams = getJobInput(jobConfigFilePath)  # config filenames + params for the combination of experiments
	processingParams = {}  # specific params for each experiment
	dfs = {}
	processingResults = {}
	experimentNames = list(jobParams['schema'].keys())
	
	if not testing:
		for eName in experimentNames:
			""" Data processing """
			# get all input parameters
			processingParams[eName] = getProcessingInput(jobParams['schema'][eName]['config'])
			# get the dataframes
			dfs[eName] = importExperimentData(processingParams[eName]['data'], delim=processingParams[eName]['delim_in'],
											  header=processingParams[eName]['header_in'],
											  wrapper=processingParams[eName]['wrapper'])
			processing_path_out = processingParams[eName]['path_out']
			processingResultsDumpFilename = os.path.join(processing_path_out, 'processingResultsDump_' + str(eName))
			if doProcessing:
				print('Starting processing of ' + eName + '...')
				# prepare the output directories
				if not os.path.exists(os.path.abspath(processing_path_out)):  # do not overwrite dir
					assert os.path.exists(
						os.path.abspath(os.path.join(processing_path_out, os.path.pardir)))  # parent dir must exist
					os.makedirs(processing_path_out)
				else:
					if jobConfigFilePath != 'job/jobConfig.ini':  # TEST
						raise Exception(
							"Output path " + os.path.abspath(processing_path_out) + " already exists! Aborting.")
				
				# process every input dataframe
				logging.info(
					"Starting processing of experiment '" + eName + "' of job '" + jobParams['jobname'] + "' at " +
					str(datetime.datetime.now()).split('.')[0])
				processingResults[eName] = processDf(dfs[eName], processingParams[eName], writeToDisk)
				logging.info(
					"Finished processing of experiment '" + eName + "' of job '" + jobParams['jobname'] + "' at " +
					str(datetime.datetime.now()).split('.')[0])
				pickle.dump(processingResults[eName], open(processingResultsDumpFilename, 'wb'))  # TEST
			elif doAnalysis:
				try:
					processingResults[eName] = pickle.load(open(processingResultsDumpFilename, 'rb'))
				except FileNotFoundError:
					raise FileNotFoundError(
						"There is no previously processed data in this path: " + processingResultsDumpFilename)
			else:
				logging.warning(
					"No processing step performed nor processing file loaded for experiment " + str(eName) + "!")
		
		""" Data analysis """
		analysis_path_out = jobParams['path_out']
		analysisResultsDumpFilename = os.path.join(analysis_path_out, 'analysisResultsDump')
		if doAnalysis:
			print('Starting analysis...')
			# prepare the output directories
			if not os.path.exists(analysis_path_out):  # do not overwrite dir
				assert os.path.exists(
					os.path.abspath(os.path.join(analysis_path_out, os.path.pardir)))  # parent dir must exist
				os.makedirs(analysis_path_out)
			
			# perform analysis
			logging.info("Starting analysis of job: " + jobParams['jobname'] + " at " +
						 str(datetime.datetime.now()).split('.')[0])
			analysisResults = analyzeProcessingResult(processingResults, jobParams, writeToDisk)
			logging.info("Finished analysis of job: " + jobParams['jobname'] + " at " +
						 str(datetime.datetime.now()).split('.')[0])
			pickle.dump(analysisResults, open(analysisResultsDumpFilename, 'wb'))  # TEST
		elif doReport:
			try:
				analysisResults = pickle.load(open(analysisResultsDumpFilename, 'rb'))
			except FileNotFoundError:
				raise FileNotFoundError(
					"There is no previously analyzed data in this path: " + analysisResultsDumpFilename)
		else:
			logging.warning("No analysis step performed nor analysis file loaded!")
		
		""" Visualize and generate report """
		results_path_out = jobParams['path_results']
		
		if doReport:
			print('Starting visualization and report generation...')
			# prepare the output directories
			if not os.path.exists(results_path_out):  # do not overwrite dir
				assert os.path.exists(
					os.path.abspath(os.path.join(results_path_out, os.path.pardir)))  # parent dir must exist
				os.makedirs(results_path_out)
			
			# visualize and make a report
			logging.info("Starting visualization end report generation of job: " + jobParams['jobname'] + "at " +
						 str(datetime.datetime.now()).split('.')[0])
			generateReport(analysisResults, jobParams, logFilePath, writeToDisk, processingParams, start)
			DB_setJobReportRelPaths(jobDirName=jobDirName, resultpath=jobParams['path_results'],
									jobName=jobParams['jobname'])
			logging.info("Finished visualization end report generation of job: " + jobParams['jobname'] + "at " +
						 str(datetime.datetime.now()).split('.')[0])
		else:
			logging.warning("No report generated!")
	
	elif testing:
		devStuff(dfs, processingParams)
	stop = time()
	print(stop - start)


if __name__ == '__main__':  # this should not execute if main.py is not the main module called by the python interpreter,
	from web.web import DB_setJobCompleted, DB_setJobFailed
	from web import app
	from traceback import print_exc
	
	args = sys.argv
	print(str(args))  # TEST
	if len(args) != 1:
		assert len(args) == 7  # todo use argparser
		jobConfigFilePath = args[1]
		doProcessing = (args[2] == 'True')
		doAnalysis = (args[3] == 'True')
		doReport = (args[4] == 'True')
		writeToDisk = (args[5] == 'True')
		testing = (args[6] == 'True')
	# so if you start main.py from within web.py or something, this won't be executed
	else:
		doProcessing = False
		doAnalysis = False
		doReport = True
		writeToDisk = True
		testing = False
		
		# jobConfigFilePath = 'job/jobConfig.ini' # TEST
		jobConfigFilePath = '../jobs/2017-01-23 02:40:26.553814_two/jobConfig_two.ini'
	# jobConfigFilePath = webFlow(exptype='COON')
	# jobConfigFilePath = webFlow(exptype='COON', previousjobdirName='2016-12-26 10:56:10.646919_COON')
	# jobConfigFilePath = webFlow(exptype='COON_SN')
	# jobConfigFilePath = webFlow(exptype='COON_SN', previousjobdirName='2016-12-20 14:21:47.786288_COON_SN')
	# jobConfigFilePath = webFlow(exptype='COON_norm') # todo constand uitzetten
	# jobConfigFilePath = webFlow(exptype='COON_norm', previousjobdirName='2016-12-12 22:43:38.030716_COON_norm')
	# jobConfigFilePath = webFlow(exptype='COON_SN_norm')  # todo constand uitzetten
	# jobConfigFilePath = webFlow(exptype='COON_SN_norm', previousjobdirName='2016-12-12 22:48:30.701250_COON_SN_norm')
	# jobConfigFilePath = webFlow(exptype='COON_nonormnoconstand')  # todo constand uitzetten
	# jobConfigFilePath = webFlow(exptype='COON_nonormnoconstand', previousjobdirName='2016-12-20 14:31:47.045927_COON_nonormnoconstand')
	# jobConfigFilePath = webFlow(exptype='COON_noISO')
	# jobConfigFilePath = webFlow(exptype='COON_noISO', previousjobdirName='2016-12-16 16:38:30.536344_COON_noISO')
	# jobConfigFilePath = webFlow(exptype='COON_SN_nonormnoconstand')  # todo constand uitzetten
	# jobConfigFilePath = webFlow(exptype='COON_SN_nonormnoconstand', previousjobdirName='2016-12-20 14:39:09.476567_COON_SN_nonormnoconstand')
		jobConfigFilePath = '/home/pdiracdelta/Documents/UHasselt/CONSTANd++/jobs/2017-04-14 10:14:53.002433_coon2test/jobConfig_coon2test.ini'
	
	with app.app_context():
		jobDirName = os.path.basename(os.path.abspath(os.path.join(jobConfigFilePath, os.pardir)))
		try:
			main(jobConfigFilePath=jobConfigFilePath, doProcessing=doProcessing, doAnalysis=doAnalysis,
				 doReport=doReport,
				 testing=testing, writeToDisk=writeToDisk)
			DB_setJobCompleted(jobDirName)
		except:
			DB_setJobFailed(jobDirName)
			print_exc()
