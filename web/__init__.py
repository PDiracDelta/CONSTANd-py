from flask import Flask, g
from flask_mail import Mail
import sqlite3
#from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)
app.config.from_object('web.config')
app.config['allJobsDir'] = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/'
app.config['DB'] = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/jobs.db'
#
# jobs(id text primary key not null, jobname text, htmlreport text, pdfreport text, done integer default 0, success integer default 0);
#
#app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:////home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/jobs.db'
mailer = Mail(app)


def get_db():
	db = getattr(g, '_database', None)
	if db is None:
		db = g._database = sqlite3.connect(app.config.get('DB'))
	return db


@app.teardown_appcontext
def close_connection(exception):
	db = getattr(g, '_database', None)
	if db is not None:
		db.close()

from web import views
