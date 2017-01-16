from runweb import app
from flask import render_template


@app.route('/')
def hello_world():
	#from main import main
	#from web.webFlow import webFlow
	#masterConfigFilePath = webFlow(exptype='COON')
	#main(jobConfigFilePath=masterConfigFilePath, doProcessing=True, doAnalysis=True, doReport=True, testing=False, writeToDisk=True)
	return render_template('index.html', reportFile='reportexample.html')


@app.route('/report/<file>')
def report(file=None):
	return render_template(file)
