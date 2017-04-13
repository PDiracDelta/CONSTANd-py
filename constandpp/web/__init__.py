"""
Constructs the web app object and its config, and defines the database connection.
"""

from flask import Flask, g  # flask.g is a global object you can use to store data on. It persists between sessions and across contexts
from flask_mail import Mail
import sqlite3


app = Flask(__name__)
app.config.from_object('web.config')
app.config['allJobsDir'] = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/'
app.config['jobDB'] = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/jobs.db'

mailer = Mail(app)


def get_db():
	"""
	Initialize a new database connection or pick up an existing one.
	:return db:	sqlite3 database connection
	"""
	db = getattr(g, '_database', None)
	if db is None:
		db = g._database = sqlite3.connect(app.config.get('jobDB'))
	return db


# gets executed when the app is torn down.
@app.teardown_appcontext
def close_connection(exception):
	"""
	Close the database connection, possibly while catching an exception.
	:param exception:	Exception	???
	"""
	db = getattr(g, '_database', None)
	if db is not None:
		db.close()

# this should really be at the end of the file, otherwise imports are messed up. This is just how Flask works...
import web.views