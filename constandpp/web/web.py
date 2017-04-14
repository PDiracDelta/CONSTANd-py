"""
Functions employed by the web interface (real one or webFlow.py) to handle:
* job preparation including input file modifications
* running the job
* database queries
* e-mailing
* misc. functions like TMT to ICM conversion (not available through web interface yet).
"""

import pandas as pd
import numpy as np
import os, datetime
from dataIO import unnest
from json import dumps
from web import app, mailer, get_db, close_connection
from subprocess import Popen, PIPE


def newJobDir(this_job_name):
	"""
	Create a new job dir (as a subdir of the jobs dir) with as name the current datetime and this_job_name.
	:param this_job_name:	str		name of the new job
	:return:				str		full path to the dir of the new job
	"""
	jobBaseDir = app.config.get('ALLJOBSDIR')
	jobPath = os.path.join(jobBaseDir, str(datetime.datetime.now()) + '_' + this_job_name)
	os.makedirs(jobPath)
	return os.path.abspath(jobPath)


def hackExperimentNamesIntoForm(form, eNames):
	"""
	Replaces the generic "experiments-x" labels in a WTForms 'experiments' field entry by the names of the experiments as
	specified in eNames, in unmodified order.
	:param form:	jobSettingsForm		Form 2 instance (see views.py) with generic 'experiments' entry labels
	:param eNames:	list(str)			names of the experiments, as extracted from job schema
	:return:		jobSettingsForm		Form 2 instance (see views.py) 'experiments' entry labels according to eNames
	"""
	import re
	
	def replace(oldstring, eName):
		newstring = re.sub(r'(e|E)xperiments-[0-9]*', eName, oldstring)
		return newstring
	#form._prefix = re.sub(r'(e|E)xperiments-[0-9]*', eName, form._prefix)
	for i in range(len(eNames)):
		form.experiments.entries[i].label.text = replace(form.experiments.entries[i].label.text, eNames[i])
		form.experiments.entries[i].name = replace(form.experiments.entries[i].name, eNames[i])
	return form


def saveFileStorage(this_path, fs, prefix=None):
	if fs.filename == '':
		return None
	else:
		fileName = fs.filename.lstrip('/')
		if prefix is not None:
			filePath = os.path.join(this_path, prefix + '_' + fileName)
		else:
			filePath = os.path.join(this_path, fileName)
		fs.save(filePath)
		return filePath


def updateSchema(this_job_path, this_incompleteSchema, form):
	#from web.forms import FileField, FormField
	for eName in this_incompleteSchema.keys():
		for experimentForm in form.experiments.entries:
			if experimentForm.name == eName:
				this_incompleteSchema[eName]['data'] = os.path.basename(saveFileStorage(this_job_path, experimentForm.data['dataFile'], prefix=eName))
				this_incompleteSchema[eName]['config'] = saveFileStorage(this_job_path, experimentForm.data['processingConfig'], prefix=eName)
				if experimentForm.data['wrapper'].filename != '':
					this_incompleteSchema[eName]['wrapper'] = os.path.basename(saveFileStorage(this_job_path, experimentForm.data['wrapper'], prefix=eName))
				else:
					wrapperFileName = eName + '_wrapper.tsv'
					open(os.path.join(this_job_path, wrapperFileName), 'w').close()
					this_incompleteSchema[eName]['wrapper'] = wrapperFileName
				this_incompleteSchema[eName]['isotopicCorrection_matrix'] = None

	return this_incompleteSchema


def DB_checkJobExist(ID):
	return get_db().execute('SELECT EXISTS(SELECT 1 FROM jobs WHERE id="' + ID + '" LIMIT 1);')


def DB_insertJob(jobDirName, jobName):
	cur = get_db().execute('INSERT INTO jobs VALUES ("' + jobDirName + '","' + jobName + '","","", 0, 0);')
	get_db().commit()
	return cur


def DB_getJobVar(ID, varName):
	return get_db().execute('SELECT '+varName+' FROM jobs WHERE id="' + ID + '" LIMIT 1;')


def updateConfigs(this_job_path, this_schema):
	import fileinput
	for eName in this_schema:
		experiment = this_schema[eName]
		configFile = os.path.join(this_job_path, experiment['config'])
		# open user config parameters
		with open(configFile, 'a') as fout, fileinput.input(app.config.get('BASEPROCESSINGCONFIG')) as fin:
			fout.write('\n')  # so you dont accidentally append to the last line
			# write baseConfig parameters
			for line in fin:
				if '[DEFAULT]' not in line:
					fout.write(line)
			# write schema parameters
			fout.write('\n')  # so you dont accidentally append to the last line
			fout.write('data = ' + experiment['data'] + '\n')
			fout.write('wrapper = ' + experiment['wrapper'] + '\n')
			# caution! use the ALIASES, and NOT the original names (they are rewritten by the wrapper)
			# fout.write('channelNamesPerCondition = ' + dumps(experiment['channelAliasesPerCondition']) + '\n')
			fout.write('intensityColumns = ' + dumps(unnest(experiment['channelAliasesPerCondition'])) + '\n')
			if experiment['isotopicCorrection_matrix'] is not None:
				fout.write('isotopicCorrection_matrix = ' + experiment['isotopicCorrection_matrix'] + '\n')
			else:
				fout.write('isotopicCorrection_matrix\n')
			# write output parameters
			fout.write('path_out = ' + eName + '_output_processing/\n')
			fout.write('filename_out = ' + str(eName) + '\n')


def updateWrappers(this_job_path, this_schema):
	# write the channel aliases to the wrapper
	for eName in this_schema:
		experiment = this_schema[eName]
		channelNames = unnest(experiment['channelNamesPerCondition'])
		channelAliases = unnest(experiment['channelAliasesPerCondition'])
		wrapperFile = os.path.join(this_job_path, experiment['wrapper'])
		# open user config parameters
		with open(wrapperFile, 'a') as fout:
			fout.write('\n')  # so you dont accidentally append to the last line
			for n, a in zip(channelNames, channelAliases):
				fout.write(n + '\t' + a + '\n')


def makeJobConfigFile(this_job_path, this_jobname, jobID, this_schema, form):
	jobConfigFullPath = os.path.join(this_job_path, 'jobConfig_'+this_jobname+'.ini')
	with open(jobConfigFullPath, 'w+') as fout:
		fout.write('[DEFAULT]\n') # todo enable
		for field in form:
			if field.name != 'csrf_token' and field.name != 'experiments':
				fout.write(field.name + ' = '+ str(field.data) + '\n')
		# HARD CODED VARIABLES
		#fout.write('\n')  # so you dont accidentally append to the last line
		fout.write('jobname = ' + this_jobname + '\n')
		fout.write('jobID = ' + jobID + '\n')
		fout.write('path_out = output_analysis\n')
		fout.write('path_results = results\n')
		fout.write('PCA_components = 2\n')
		fout.write('schema = ' + dumps(this_schema) + '\n')
		fout.write('date = ' + str(os.path.basename(this_job_path).split('.')[0]) + '\n')
	return jobConfigFullPath


def DB_setJobReportRelPaths(jobDirName, resultpath, jobName):
	get_db().execute('UPDATE jobs SET htmlreport = "' + os.path.join(resultpath, jobName+'_report.html')
					 + '", pdfreport = "' + os.path.join(resultpath, jobName+'_Report.pdf') + '" WHERE id = "' + jobDirName + '";')
	get_db().commit()


def DB_setJobCompleted(jobDirName):
	get_db().execute('UPDATE jobs SET done = 1, success = 1 WHERE id = "' + jobDirName + '";')
	get_db().commit()


def DB_setJobFailed(jobDirName):
	get_db().execute('UPDATE jobs SET done = 1, success = 0 WHERE id = "'+jobDirName+'";')
	get_db().commit()


def startJob(jobConfigFullPath):
	return Popen(['python3', app.config.get('MAIN'), jobConfigFullPath,
		 'True', # doProcessing
		 'True', # doAnalysis
		 'True', # doReport
		 'True', # writeToDisk
		 'False', # testing
		 ], shell=False, stdout=PIPE, stderr=PIPE, stdin=PIPE)


def TMT2ICM(TMTImpuritiesDF, order=None):
	"""
	Converts a dataframe of TMT-like isotopic impurities (indexed on TMT label name) into the correct isotopic
	correction matrix. Column order from the dataframe is changed according to 'order'. Rows are normalized so that
	their sums are equal to one.
	:param TMTImpuritiesDF: pd.DataFrame    TMT-like isotopic impurities
													-2  -1  +1  +2
											126     0   0   1.2 0
											127N    1.2 3.3 2.5 0.3
											127C    ...
											...
	:return ICM:            np.ndarray      isotopic corrections matrix
	"""
	# cols6plex = ['126', '127', '128', '129', '130', '131']
	# cols8plex = ['126', '127N', '127C', '128C', '129N', '129C', '130C', '131']
	# cols10plex = ['126', '127N', '127C', '128N', '128C', '129N', '129C', '130N', '130C', '131']
	Nplex = len(TMTImpuritiesDF)
	channelNames = list(TMTImpuritiesDF.index.values.astype(str))
	labelNames = ['O_'+n for n in channelNames] # O_ for Observed_
	# create empty ICM-dataframe with 100 on the diagonals and zeroes elsewhere
	icmdf = pd.DataFrame(np.eye(Nplex)*100, index=labelNames, columns=channelNames).fillna(0)

	# build the dictionary with correspondents
	if Nplex == 6: #sixplex
		correspondents = {}
		for k in range(126,132):
			correspondents[str(k)] = {'-2':str(k-2), '-1':str(k-1), '+1':str(k+1), '+2':str(k+2)}
		correspondents['126']['-2'] = 'nobody'
		correspondents['126']['-1'] = 'nobody'
		correspondents['127']['-2'] = 'nobody'
		correspondents['130']['+2'] = 'nobody'
		correspondents['131']['+1'] = 'nobody'
		correspondents['131']['+2'] = 'nobody'
	elif Nplex in [8, 10]: # 8- and 10-plex
		correspondents = {'126' : {'-2':'nobody', '-1':'nobody', '+1':'127C', '+2':'128N'},
				  '127N': {'-2': 'nobody', '-1': 'nobody', '+1': '128N', '+2': '128C'},
				  '127C': {'-2': 'nobody', '-1': '126', '+1': '128C', '+2': '129N'},
				  '128N': {'-2': 'nobody', '-1': '127N', '+1': '129N', '+2': '129C'},
				  '128C': {'-2': '126', '-1': '127C', '+1': '129C', '+2': '130N'},
				  '129N': {'-2': '127N', '-1': '128N', '+1': '130N', '+2': '130C'},
				  '129C': {'-2': '127C', '-1': '128C', '+1': '130C', '+2': '131'},
				  '130N': {'-2': '128N', '-1': '129N', '+1': '131', '+2': 'nobody'},
				  '130C': {'-2': '128C', '-1': '129C', '+1': 'nobody', '+2': 'nobody'},
				  '131': {'-2': '129N', '-1': '130N', '+1': 'nobody', '+2': 'nobody'}}
	else:
		raise Exception("Illegal plexity of your TMT labels. Only 6plex, 8plex, 10plex are supported.")

	# execute mappings
	for trueChannel in channelNames: #for each row in TMTImpurities
		for TMTisotope, observedChannel in correspondents[trueChannel].items(): # look up each isotope correspondent...
			# ... and put the TMT input value of that isotope of the true channel into the icmdf location where the
			# transfer from the observed channel to the true channel is stored. Obs = ICM * True
			if observedChannel in channelNames: # (TMT 8plex is a subset of 10plex: some observed channels don't exist.
				if observedChannel != 'nobody': # (only if the observed channel exists of course)
					icmdf.loc['O_'+observedChannel, trueChannel] = TMTImpuritiesDF.loc[trueChannel, TMTisotope]

	if order is not None: # reorder according to 'order' of the channelAliasesPerCondition
		order = unnest(order)
		assert len(order) == Nplex
		icmdf = icmdf.reindex_axis(order, axis=1).reindex_axis(['O_'+i for i in order], axis=0)
	icm = np.asmatrix(icmdf)
	rowSums = np.asarray(np.sum(icm, axis=1))
	return icm / rowSums # normalize each row so that the sum is one


def constructMasterConfigContents(schemaDict, otherMasterParams): #todo obsolete
	""" OBSOLETE """
	# todo docu
	def isNumeric(s):
		""" Returns whether or not the argument is numeric """
		try:
			float(s)
			return True
		except ValueError:
			return False
	contents = {}
	contents['schema'] = dumps(schemaDict)
	for k, v in otherMasterParams.items():
		if isinstance(v, str) or isNumeric(v):
			if k.split('_')[0] == 'delim': # delimiters should be saved in a visible format
				from codecs import getencoder as ge, getdecoder as gd
				byteDelim = ge("unicode-escape")(v)[0]
				contents[k] = gd("utf-8")(byteDelim)[0]
			else:  # not a delimiter
				contents[k] = v
		else:
			contents[k] = dumps(v)

	return contents


def send_mail(recipient, mailBodyFile, jobname, jobID, attachment):  # TODO VITO credentials
	from flask_mail import Message
	from socket import gaierror
	subject = "Your CONSTANd++ job %s" % (jobname)
	#body = "=== English version below ===\n\n"
	#body += "Beste {0} {1},\n\n"
	with open(os.path.join(app.config.root_path, 'static/'+mailBodyFile), 'r') as f:
		body = f.read()
	body = body.format(recipient, jobname, jobID)

	msg = Message(subject, recipients=[recipient], body=body, sender=("CONSTANd++", "constand@vito.be"))

	if attachment is not None:
		assert os.path.exists(attachment)
		aName = os.path.basename(attachment)
		with app.open_resource(attachment) as a:
			msg.attach(filename=aName, content_type='application/pdf', data=a.read())
	try:
		mailer.send(msg)
		return None
	except gaierror as e:  # return an error to inform the user no e-mail was sent.
		return "You will probably not receive an e-mail, because of the following internal error: "+str(e)
