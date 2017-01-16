import os
from web import app, mailer
from flask import render_template, send_from_directory
from flask_mail import Message


#############################
#Client side
#############################

@app.route('/')
def home():
	#from main import main
	#from web.webFlow import webFlow
	#masterConfigFilePath = webFlow(exptype='COON')
	#main(jobConfigFilePath=masterConfigFilePath, doProcessing=True, doAnalysis=True, doReport=True, testing=False, writeToDisk=True)
	return render_template('index.html', reportFile='reportexample.html')


@app.route('/report/<file>')
def report(file=None):
	return render_template(file)


@app.route('/jobs/<path:filename>')
def getFile(filename):
	return send_from_directory('../../doc/figures/reportexample/', filename, as_attachment=True)


@app.route('/docu')
def documentation():
	return render_template('documentation.html', title="Documentation")


@app.route('/newjob')
def newjob():
	return render_template('newjob.html')


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