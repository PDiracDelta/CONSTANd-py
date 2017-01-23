import os

WTF_CSRF_ENABLED = False # Enable CSRF protection
SECRET_KEY = 'azezfhbvdujhjnzeosdkohefsduiqs'

basedir = os.path.abspath(os.path.dirname(__file__))

MAIL_DEFAULT_SENDER = ('CONSTANd++', 'constand@vito.be')
MAIL_SERVER = 'smtp.vito.local'

#SERVER_NAME = 'http://localhost:5000' # THIS BREAKS THE BIND SPECIFIED IN app.run(host=...)
UPLOAD_FOLDER = '/home/pdiracdelta' # filesystem path where uploaded posters will be stored
ALLJOBSDIR = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/jobs/'
ALLOWED_EXTENSIONS = {'tsv', 'txt', 'xlsx', 'csv', 'xls', ''}  # Allowed extensions for uploads
BASEPROCESSINGCONFIG = os.path.abspath(os.path.join(basedir,'../baseProcessingConfig.ini'))
MAIN = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/scripts/main.py'