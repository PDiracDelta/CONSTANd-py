import os
from web import app, mailer
from flask import render_template, send_from_directory, request, redirect, url_for
from flask_mail import Message
from .forms import newJobForm
from werkzeug.utils import secure_filename


#############################
#Client side
#############################

@app.route('/')
def home():
	#from main import main
	#from web.webFlow import webFlow
	#masterConfigFilePath = webFlow(exptype='COON')
	#main(jobConfigFilePath=masterConfigFilePath, doProcessing=True, doAnalysis=True, doReport=True, testing=False, writeToDisk=True)
	return render_template('index.html', reportFile='reportexample.html', form=newJobForm(csrf_enabled=False))


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
def newjob():
	form = newJobForm(request.form, csrf_enabled=False)
	#if request.method == 'POST' and form.validate():
	if form.validate_on_submit():
		jobName = form.jobName.data
		from web.web import newJobDir
		jobDir = newJobDir(jobName)
		#schemaFileName = 'schema_'+secure_filename(request.files[form.schema.name])#form.schema.data.filename)
		schemaFileName = 'schema_' + secure_filename(request.files[request.files['schema']])
		schemaFilePath = os.path.join(jobDir, schemaFileName)
		form.schema.file.save(schemaFilePath)
		from dataIO import parseSchemaFile
		try:
			incompleteSchema = parseSchemaFile(schemaFilePath)
		except Exception as e:  # remove file and job dir if something went wrong
			os.remove(schemaFilePath)
			os.removedirs(jobDir)
			return redirect(url_for('newjob')) # todo give error message
		return redirect(url_for('jobSettings', incompleteSchema))
	return render_template('newjob.html', form=form)


@app.route('/jobSettings', methods=['GET', 'POST'])
def jobSettings(incompleteSchema):
	return str(incompleteSchema)


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