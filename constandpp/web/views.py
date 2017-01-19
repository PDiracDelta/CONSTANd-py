import os
from web import app, mailer
from flask import render_template, send_from_directory, request, redirect, url_for, session, flash
from flask_mail import Message
from .forms import newJobForm, experimentForm, jobSettingsForm
from werkzeug.utils import secure_filename
from web import get_db
from subprocess import run

#############################
#Client side
#############################

@app.route('/')
def home():
	#from main import main
	#from web.webFlow import webFlow
	#masterConfigFilePath = webFlow(exptype='COON')
	#main(jobConfigFilePath=masterConfigFilePath, doProcessing=True, doAnalysis=True, doReport=True, testing=False, writeToDisk=True)
	return render_template('index.html', reportFile='reportexample.html', form=newJobForm())#csrf_enabled=False))


@app.route('/report/<file>')
def report(file=None):
	return render_template(file)


@app.route('/jobs/<path:filename>')
def getFile(filename):
	return send_from_directory('../../doc/figures/reportexample/', filename, as_attachment=True)


@app.route('/docu')
def documentation():
	return render_template('documentation.html', title="Documentation")


@app.route('/newjob', methods=['GET', 'POST'])
def newJob():
	#form = newJobForm(request.form)#, csrf_enabled=False)
	form = newJobForm()
	#if request.method == 'POST' and form.validate():
	if form.validate_on_submit():
		jobName = form.jobName.data
		session['jobName'] = jobName
		from web.web import newJobDir
		jobDirPath = newJobDir(jobName)
		session['jobDirName'] = os.path.basename(jobDirPath)
		#schemaFileName = 'schema_'+secure_filename(request.files[form.schema.name])#form.schema.data.filename)
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
			return redirect(url_for('newjob'))
		return redirect(url_for('jobSettings'))
	return render_template('newjob.html', form=form)


@app.route('/jobsettings', methods=['GET', 'POST'])
def jobSettings():
	incompleteSchema = session.get('incompleteSchema')
	# eforms = {(eName, experimentForm()) for eName in incompleteSchema}
	form = jobSettingsForm()
	if form.validate_on_submit():
		from web import updateSchema
		updateSchema(os.path.join(app.config.get('allJobsDir'), session['jobDirName']), incompleteSchema, form)
		cur = get_db().execute('SELECT EXISTS(SELECT 1 FROM jobs WHERE id="'+session['jobDirName']+'" LIMIT 1);')
		if cur.fetchall(): # already exists
			redirect(url_for('jobinfo'))
		else: # does not exist yet
			run('python3 '+'"/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/scripts/main.py" '
			    +' True', #doProcessing
			    +' True', #doAnalysis
			    +' True', #doReport
			    +' False', #testing
			    +' True', #writeToDisk
			    +' &',
			    shell=True) # RUN CONSTANd++ IN INDEPENDENT
			cur = get_db().execute('INSERT INTO jobs VALUES ("' + session['jobDirName'] + '","' + session['jobName'] + '","","", 0, 0);')
	else:
		numExperiments = len(incompleteSchema)
		#form.experiments.min_entries = 3#numExperiments
		for eName in incompleteSchema:
			form.experiments.append_entry(experimentForm(prefix=eName)) #{'title': session["experiments"][pif][0]}
		return render_template('jobsettings.html', jobName=session.get('jobName'), form=form)


@app.route('/jobinfo', methods=['GET', 'POST'])
def jobInfo():
	jobID = request.args.get('id', '')
	cur = get_db().execute('SELECT done FROM jobs WHERE id="' + jobID + '" LIMIT 1;')
	isDone = cur.fetchall()
	if isDone is not None:
		if isDone:
			cur = get_db().execute('SELECT htmlreport FROM jobs WHERE id="' + jobID + '" LIMIT 1;')
			htmlreportPath = cur.fetchall()
			pdfreportPath = htmlreportPath[0:-4]+'pdf'
			return render_template('jobinfo.html', done=True, html=htmlreportPath, pdf=pdfreportPath)
		else:
			return render_template('jobinfo.html', done=False)
	else:
		return "Couldn't find that job (or something else went wrong)."


@app.route('/htmlreport/<path:jobID>', methods=['GET', 'POST'])
def getHtmlReport(jobID):
	htmlFileName = request.args.get('htmlFileName', '')
	return send_from_directory(app.config.get('allJobsDir')+jobID, htmlFileName, as_attachment=True)


@app.route('/pdfreport/<path:jobID>', methods=['GET', 'POST'])
def getPdfReport(jobID):
	pdfFileName = request.args.get('pdfFileName', '')
	return send_from_directory(app.config.get('allJobsDir')+jobID, pdfFileName, as_attachment=True)


#############################
#Admin functions
#############################

def send_mail(recipient, mailBodyFile, jobname, jobID, attachment): # TODO VITO credentials
	subject = "Your CONSTANd++ job %s" % (jobname)
	#body = "=== English version below ===\n\n"
	#body += "Beste {0} {1},\n\n"
	with open(os.path.join('static', mailBodyFile), 'r') as f:
		body = f.read()
	body = body.format(recipient, jobname, jobID)

	msg = Message(subject, recipients=[recipient], body=body, sender=("ULYSSIS", "ulyssis@ulyssis.org"))
	assert os.path.exists(attachment)
	msg.attach(attachment)
	mailer.send(msg)