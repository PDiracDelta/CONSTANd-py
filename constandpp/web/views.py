from runweb import app

@app.route('/')
def hello_world():
	from main import main
	from web.webFlow import webFlow
	masterConfigFilePath = webFlow(exptype='COON')
	main(jobConfigFilePath=masterConfigFilePath, doProcessing=True, doAnalysis=True, doReport=True, testing=False, writeToDisk=True)
	return 'Hello, Crazy World!'
