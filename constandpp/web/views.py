from web import app
from flask import render_template, send_from_directory


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

@app.route('/jobs/<path:filename>')
def getImage(filename):
	return send_from_directory('../../doc/figures/reportexample/', filename, as_attachment=True)