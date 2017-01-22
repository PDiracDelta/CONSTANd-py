import os

WTF_CSRF_ENABLED = False # Enable CSRF protection
SECRET_KEY = 'azezfhbvdujhjnzeosdkohefsduiqs'

basedir = os.path.abspath(os.path.dirname(__file__))

UPLOAD_FOLDER = '/home/pdiracdelta' # filesystem path where uploaded posters will be stored
ALLOWED_EXTENSIONS = {'tsv', 'txt', 'xlsx', 'csv', 'xls', ''}  # Allowed extensions for uploads
BASEPROCESSINGCONFIG = os.path.abspath(os.path.join(basedir,'../baseProcessingConfig.ini'))
MAIN = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/scripts/main.py'