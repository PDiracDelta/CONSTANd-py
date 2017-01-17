import os

WTF_CSRF_ENABLED = True # Enable CSRF protection
SECRET_KEY = 'azezfhbvdujhjnzeosdkohefsduiqs'

basedir = os.path.abspath(os.path.dirname(__file__))

UPLOAD_FOLDER = '/path/to/poster_folder' # filesystem path where uploaded posters will be stored
ALLOWED_EXTENSIONS = {'tsv', 'txt', 'xlsx', 'csv', 'xls', ''}  # Allowed extensions for uploads