"""
Config file with all parameters used in the flask web app.
"""
import os

#todo check CSRF because actually i think it is enabled anyway...
WTF_CSRF_ENABLED = False  # Enable CSRF protection
SECRET_KEY = 'azezfhbvdujhjnzeosdkohefsduiqs'

basedir = os.path.abspath(os.path.dirname(__file__))

MAIL_DEFAULT_SENDER = ('CONSTANd++', 'constand@vito.be')
MAIL_SERVER = 'smtp.vito.local'

#SERVER_NAME = 'http://localhost:5000' # THIS BREAKS THE BIND SPECIFIED IN app.run(host=...)
UPLOAD_FOLDER = '/home/pdiracdelta'  # filesystem path where uploaded posters will be stored
ALLJOBSDIR = '/home/pdiracdelta/Documents/UHasselt/CONSTANd++/jobs/'
JOBDB = os.path.join(ALLJOBSDIR, 'jobs.db')
ALLOWED_EXTENSIONS = {'tsv', 'txt', 'xlsx', 'csv', 'xls', ''}  # Allowed extensions for uploads
BASEPROCESSINGCONFIG = os.path.abspath(os.path.join(basedir, '../baseProcessingConfig.ini'))
MAIN = '/var/www/CONSTANd++/main.py'