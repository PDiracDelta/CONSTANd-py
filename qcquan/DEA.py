# -*- coding: utf-8 -*-

"""
Contains all steps to be undertaken in a single DEA operation:
  - map from peptide to protein level
  - perform DE test
  - calculate fold change
  - calculate significance
"""

from qcquan.analysis import *
from tools import getOtherConditions, unnest


def DEA(this_allMSRunsDF, proteinPeptidesDict, params):
	"""
	Bring the data to the protein level in the case of [minimal]{full} expression (shared peptides are [not allowed]{allowed}).
	Execute the differential expression analysis (t-test + B-H correction, compute log2 fold change and some useful
	extra columns) and gather some metadata.
	:param this_allMSRunsDF:		pd.DataFrame	dataframe containing an outer join of all the MSRun dataframes
	:param proteinPeptidesDict:			dict			{ protein : all associated peptide indices (non-injective) }
	:param params:						dict			job (global) parameters
	:return proteinDF:					pd.DataFrame	Transformed and selected data on the protein level.
														Structure:
														['peptides', 'description', {{condition}}, {{'p-value', 'adjusted p-value'}},
														{{'log2 fold change'}}, {{'significant'}}, {{'#peptides'}}]
	:return numProteins:				int				number of proteins taken into account in the DEA
	"""
	referenceCondition = params['referenceCondition']
	# use list() so that the new variable is not an alias
	otherConditions = getOtherConditions(params['schema'], referenceCondition)
	# execute mappings to get all peptideintensities per protein, over each whole condition. Index = 'protein'
	proteinDF = getProteinDF(this_allMSRunsDF, proteinPeptidesDict, params['schema'],
							 referenceCondition=referenceCondition, otherConditions=otherConditions)
	
	# perform differential expression analysis with Kai Kammer's moderated t-test
	# and do a Benjamini-Hochberg correction. Also remove proteins that have all
	# nan values for a certain condition and keep the removed ones in metadata
	proteinDF = testDifferentialExpression(proteinDF, params['alpha'], referenceCondition, otherConditions, params['schema'])
	numProteins = len(proteinDF)
	
	# calculate fold changes of the average protein expression value per CONDITION/GROUP (not per channel!)
	proteinDF = applyFoldChange(proteinDF, params['pept2protCombinationMethod'], referenceCondition,
								otherConditions)
	
	# indicate significance based on given thresholds alpha and FCThreshold
	proteinDF = applySignificance(proteinDF, otherConditions, params['alpha'], params['FCThreshold'])
	
	return proteinDF, numProteins


def get_arbitrary_column_names_by_condition(schema, referenceCondition):
	""" For each condition create as many arbitrary column names as necessary, in a dict, reference first. """
	from collections import OrderedDict
	otherConditions = getOtherConditions(schema, referenceCondition)
	allConditions = [referenceCondition, ] + otherConditions
	D_N_channels = OrderedDict.fromkeys(allConditions)
	for k in D_N_channels.keys():
		D_N_channels[k] = 0
	for c in allConditions:
		for r in schema['allMSRuns']:
			if c in schema[r]:
				# we want the NUMBER of channels, NOT actual names (we don't know which sample each value comes from)
				D_N_channels[c] += len(schema[r][c]['channelAliases'])
	D_arbitrary_column_names = {c: [c+'#'+str(i) for i in range(1, 1+n)] for c, n in D_N_channels.items()}
	return D_arbitrary_column_names


def get_design_matrix(D_arbitrary_column_names_by_condition):
	""" ANOVA-like design matrix for use in moderated_ttest, indicating group (condition) membership of each entry in
	all_channels.
	:param	D_arbitrary_column_names_by_condition	OrderedDict		All conditions as keys, referenceCondition first(!)
	"""
	otherConditions = list(D_arbitrary_column_names_by_condition.keys())[1:]
	all_channels = unnest(list(D_arbitrary_column_names_by_condition.values()))
	N_channels = len(all_channels)
	N_conditions = 1+len(otherConditions)
	design = np.zeros((N_channels, N_conditions), dtype=int)
	design[:, 0] = 1  # reference gets 1 everywhere
	for con in otherConditions:  # for each channel in each condition, put a "1" in the design matrix.
		for chan in D_arbitrary_column_names_by_condition[con]:
			design[all_channels.index(chan), 1+otherConditions.index(con)] = 1
	return design


def get_protein_intensities_as_long_format(proteinDF, D_arbitrary_column_names_by_condition):
	""" Returns one numerical value in per column (ARBITRARY sample) for each protein, in contrast with
	getProteinDF which returns a list per condition for each protein. Spans all MS runs. """
	conditions = D_arbitrary_column_names_by_condition.keys()
	all_arbitrary_column_names = unnest(D_arbitrary_column_names_by_condition.values())
	protein_intensities = pd.DataFrame(index=proteinDF.index, columns=all_arbitrary_column_names, dtype=float)
	for row in proteinDF.iterrows():
		entry = pd.Series(index=protein_intensities.columns)
		for c in conditions:
			entry[D_arbitrary_column_names_by_condition[c][:len(row[1][c])]] = row[1][c]
		protein_intensities.loc[row[0], :] = entry
	return protein_intensities


def testDifferentialExpression(this_proteinDF, alpha, referenceCondition, otherConditions, schema):
	"""
	Perform a moderated t-test (https://doi.org/10.2202/1544-6115.1055, Smyth2004.pdf) for independent samples for each
	protein on its (2) associated lists of peptide quantifications,	do a Benjamini-Hochberg correction and store the
	results in new "p-values" and "adjusted p-values" columns in the dataframe.
	This makes use of the moderated_ttest.R (which also includes the BH correction) script via rpy2.
	:param this_proteinDF:		pd.DataFrame	data from all MSRuns on the protein level
	:param alpha:				float			confidence level for the t-test
	:param referenceCondition:	str				name of the condition to be used as the reference
	:param otherConditions:		list 			all non-reference conditions in the MSRun
	:return this_proteinDF:		pd.DataFrame	data on protein level, with statistical test information
	"""
	import rpy2.robjects as ro
	from rpy2.robjects.packages import importr
	from rpy2.robjects import r, pandas2ri, numpy2ri
	from rpy2.robjects.conversion import localconverter

	FP_moderated_ttest = '/var/www/QCQuan/qcquan/qcquan/moderated_ttest.R'
	rbase = importr('base')
	importr('limma')
	rbase.source(FP_moderated_ttest)
	numpy2ri.activate()
	
	D_arbitrary_column_names_by_condition = get_arbitrary_column_names_by_condition(schema, referenceCondition)
	df_intensities = get_protein_intensities_as_long_format(this_proteinDF, D_arbitrary_column_names_by_condition)
	design = get_design_matrix(D_arbitrary_column_names_by_condition)
	# put df into R, call moderated ttest on it, and extract the result from R.
	with localconverter(ro.default_converter + pandas2ri.converter):
		r_df_intensities = ro.conversion.py2rpy(df_intensities)
	r_moderated_ttest = r['moderated_ttest']
	# assign design matrix colnames using black magic https://stackoverflow.com/a/38808519
	r_colnames = r["colnames<-"]
	r_design = r_colnames(design, list(D_arbitrary_column_names_by_condition.keys()))
	r_result_mtt = r_moderated_ttest(r_df_intensities, r_design)
	numpy2ri.deactivate()
	with localconverter(ro.default_converter + pandas2ri.converter):
		# in rpy2v3.2.0 rpy2py apparently doesn't work properly, so we need a pd.DataFrame call
		result_mtt = pd.DataFrame(ro.conversion.rpy2py(r_result_mtt))
	for condition in otherConditions:
		pValueColumn = 'p-value (' + condition + ')'
		# use .values otherwise it sets values to NaN because the indices don't match
		this_proteinDF.loc[:, pValueColumn] = result_mtt.loc[:, 'p.mod_' + condition].values
		this_proteinDF['adjusted ' + pValueColumn] = result_mtt['q.mod_' + condition].values
		# set entries where at a condition has <=1 observations to np.nan
		indices_too_few_observations = ~np.logical_and(this_proteinDF['#peptides (' + condition + ')'] > 1,
													   this_proteinDF['#peptides (' + referenceCondition + ')'] > 1)
		this_proteinDF.loc[indices_too_few_observations, [pValueColumn, 'adjusted ' + pValueColumn]] = np.nan
		# If only 1 sample per condition, then each protein has only 1 value for each condition in the t-test,
		# because the multiple values _within_ each sample have already been averaged (they are correlated).
		if all(indices_too_few_observations):
			logging.warning("No (adjusted) p-values could be calculated. You probably have only 1 sample per condition.")
	return this_proteinDF


def applyFoldChange(proteinDF, pept2protCombinationMethod, referenceCondition, otherConditions):
	"""
	Calculate the log2 fold change of the quantification values per channel for each protein according to
	pept2protCombinationMethod and add it to the new "log2 fold change (conditionName)" columns.
	:param proteinDF:					pd.DataFrame	data on the protein level with t-test results.
	:param pept2protCombinationMethod:  str				method for reducing peptide information into one figure per protein
	:param referenceCondition:	str				name of the condition to be used as the reference
	:param otherConditions:		list 			all non-reference conditions in the MSRun
	:return proteinDF:					pd.DataFrame	data on the protein level, including fold changes
	"""
	for condition in otherConditions:
		if pept2protCombinationMethod == 'mean':
			proteinDF['log2 fold change ('+condition+')'] = proteinDF.apply(lambda x: np.log2(np.nanmean(x[condition])/np.nanmean(x[referenceCondition])), axis=1)
		elif pept2protCombinationMethod == 'median':
			proteinDF['log2 fold change ('+condition+')'] = proteinDF.apply(lambda x: np.log2(np.nanmedian(x[condition]) / np.nanmedian(x[referenceCondition])), axis=1)
		else:
			raise Exception("Illegal pept2protCombinationMethod '"+str(pept2protCombinationMethod)+"'.")
	return proteinDF


def applySignificance(df, otherConditions, alpha, FCThreshold):
	"""
	Adds a column with the significance level to the dataframe of proteins; specifies whether the DEA or fold change
	results or both were significant.
	:param df:          	pd.DataFrame    proteins with their DEA and FC results.
	:param otherConditions:	[ str ]			names of all non-reference conditions
	:param alpha:       	float           significance level
	:param FCThreshold: 	float           fold change threshold
	:return:            	pd.DataFrame    protein data with significance levels 'yes', 'no', 'p' or 'fc'.
	"""
	for condition in otherConditions:
		def significant(x):
			pvalueSignificant = x['adjusted p-value ('+condition+')'] < alpha
			FCSignificant = abs(x['log2 fold change ('+condition+')']) > FCThreshold
			try:
				if pvalueSignificant & FCSignificant:
					return 'yes'
				elif pvalueSignificant:
					return 'p'
				elif FCSignificant:
					return 'fc'
				else:
					return 'no'
			except TypeError:
				print("egfs")  # TEST
		
		df['significant ('+condition+')'] = df.apply(significant, axis=1)
	return df