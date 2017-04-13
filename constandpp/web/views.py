import os
from web import app
from flask import render_template, send_from_directory, request, redirect, url_for, session, flash
from .forms import newJobForm, experimentForm, jobSettingsForm
from werkzeug.utils import secure_filename
from web.web import updateSchema, DB_checkJobExist, DB_insertJob, DB_getJobVar, updateConfigs, updateWrappers, \
	makeJobConfigFile, startJob


#############################
# Client side
#############################

@app.route('/')
def home():
	"""
	Show the homepage with a new job form embedded into it.
	"""
	return render_template('index.html', reportFile='reportexample.html', form=newJobForm())


@app.route('/report/html')
def htmlreport():
	"""
	Download the html report.
	"""
	jobID = request.args.get('id', '')
	cur = DB_getJobVar(jobID, 'htmlreport')
	htmlreportName = os.path.basename(cur.fetchall()[0][0])
	resultsFullPath = os.path.join(app.config.get('ALLJOBSDIR'), jobID + '/results/')
	# htmlFileName = request.args.get('htmlreport', '')
	return send_from_directory(resultsFullPath, htmlreportName)


@app.route('/report/pdf')
def pdfreport():
	"""
	Download the pdf report.
	"""
	jobID = request.args.get('id', '')
	cur = DB_getJobVar(jobID, 'pdfreport')
	pdfreportName = os.path.basename(cur.fetchall()[0][0])
	resultsFullPath = os.path.join(app.config.get('ALLJOBSDIR'), jobID + '/results/')
	# pdfFileName = request.args.get('pdfreport', '')
	return send_from_directory(resultsFullPath, pdfreportName, as_attachment=True)


@app.route('/file', methods=['GET', 'POST'])
def getFile():
	"""
	Get any file from the server.
	"""
	fileFullPath = request.args.get('fileFullPath', '')
	dirFullPath = os.path.dirname(fileFullPath)
	fileName = os.path.basename(fileFullPath)
	return send_from_directory(dirFullPath, fileName, as_attachment=True)


@app.route('/docu')
def documentation():
	"""
	Show the documentation.
	"""
	return render_template('documentation.html', title="Documentation")


@app.route('/newjob', methods=['GET', 'POST'])
def newJob():
	"""
	Show the new job page containing Form 1, using as an input:
	- No data:		show fresh new job form newJobForm and render the page.
	- Form 1 data:  create a job dir on disk and save+verify the uploaded schema.
					If invalid, redirect to newJob (no data).
					If valid, redirect to jobSettings with Form 1 data.
	"""
	### STEP 1: get schema and create new job
	form = newJobForm()
	if form.validate_on_submit():  # new job had already been (partially) created
		jobName = form.jobName.data
		session['jobName'] = jobName
		from web.web import newJobDir
		jobDirPath = newJobDir(jobName)
		session['jobDirName'] = os.path.basename(jobDirPath)
		# schemaFileName = 'schema_'+secure_filename(request.files[form.schema.name])#form.schema.data.filename)
		schemaFileName = 'schema_' + secure_filename(form.schema.data.filename)
		schemaFilePath = os.path.join(jobDirPath, schemaFileName)
		form.schema.data.save(schemaFilePath)
		from dataIO import parseSchemaFile
		try:
			session['incompleteSchema'] = parseSchemaFile(schemaFilePath)
		except Exception:  # remove file and job dir if something went wrong
			os.remove(schemaFilePath)
			os.removedirs(jobDirPath)
			session.clear()
			flash("Invalid schema file format. Please refer to the documentation.")
			return redirect(url_for('newJob'))
		return redirect(url_for('jobSettings'))
	return render_template('newjob.html',
						   form=form)  # you're beginning a new request on the /newjob page instead of home.


@app.route('/jobsettings', methods=['GET', 'POST'])
def jobSettings():
	"""
	Create and show the job settings configuration page containing Form 2, using as input:
	- Form 1 data:	first populate the Form2.1 fields (upload buttons), then render and show the page.
	- Form 2 data:	update the schema, config and wrapper files using the data, make the jobConfig file, create an entry
					in the job database, RUN the job, SEND a job start email and redirect to jobInfo passing the job ID.
	"""
	incompleteSchema = session.get('incompleteSchema')
	eNames = list(incompleteSchema.keys())
	# eforms = {(eName, experimentForm()) for eName in incompleteSchema}
	# CREATE FORM
	form = jobSettingsForm()
	from web.web import hackExperimentNamesIntoForm, send_mail
	if form.validate_on_submit():  # form has been filled already
		jobID = session.get('jobDirName')
		form = hackExperimentNamesIntoForm(form, eNames)
		jobDir = os.path.join(app.config.get('ALLJOBSDIR'), jobID)
		### STEP 2: upload data files, wrapper files, ICM files and config (files), while updating schema with their locations.
		schema = updateSchema(jobDir, incompleteSchema, form)
		### STEP 3: update config files and wrapper files
		updateConfigs(jobDir, schema)
		updateWrappers(jobDir, schema)
		### STEP 4: get masterConfig from web and update it (add schema, date, path_out, path_results)
		jobConfigFullPath = makeJobConfigFile(this_job_path=jobDir, this_jobname=session.get('jobName'), jobID=jobID,
											  this_schema=schema, form=form)
		cur = DB_checkJobExist(jobID)
		if not cur.fetchall()[0][0]:  # job does not already exist --> create it in the DB and start it.
			cur = DB_insertJob(jobID, session.get('jobName'))
			cur.close()
			### RUN CONSTANd++ in independent subprocess ###
			jobProcess = startJob(jobConfigFullPath)
		### SEND JOB START MAIL ###
		send_mail(recipient=form.mailRecipient.data, mailBodyFile='jobstartedMail', jobname=session.get('jobName'),
				  jobID=jobID, attachment=None)
		return redirect(url_for('jobInfo', id=session.get('jobDirName')))
	
	elif len(form.experiments.entries) == 0:  # haven't completed the form yet: populate Form 2.1 (upload buttons)
		# form.experiments.label.text = 'experiment'
		for i in range(len(eNames)):
			form.experiments.append_entry(experimentForm(prefix=eNames[i]))
		form = hackExperimentNamesIntoForm(form, eNames)
	return render_template('jobsettings.html', jobName=session.get('jobName'), form=form)


@app.route('/jobinfo', methods=['GET', 'POST'])
def jobInfo():
	"""
	Shows the status of the job specified through the URL parameter'id'.
	"""
	jobID = request.args.get('id', '')
	if jobID:  # id has been set
		cur = DB_getJobVar(jobID, 'done')
		isDone = cur.fetchall()[0][0]
		if isDone is not None:
			if isDone:
				return render_template('jobinfo.html', done=True, jobID=jobID, jobName=session.get('jobName'))
			else:
				return render_template('jobinfo.html', done=False, jobID=jobID, jobName=session.get('jobName'),
									   autorefresh=5)
		else:
			return "Couldn't find that job (or something else went wrong)."
	else:
		return render_template('jobinfo.html')

#
# @app.route('/htmlreport/<path:jobID>', methods=['GET', 'POST'])
# def getHtmlReport(jobID):
# 	htmlFileName = request.args.get('htmlFileName', '')
# 	return send_from_directory(app.config.get('ALLJOBSDIR')+jobID, htmlFileName, as_attachment=True)
#
#
# @app.route('/pdfreport/<path:jobID>', methods=['GET', 'POST'])
# def getPdfReport(jobID):
# 	pdfFileName = request.args.get('pdfFileName', '')
# 	return send_from_directory(app.config.get('ALLJOBSDIR')+jobID, pdfFileName, as_attachment=True)
