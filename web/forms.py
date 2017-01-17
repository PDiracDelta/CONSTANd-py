from flask_wtf import FlaskForm
#from flask_wtf.file import FileField
from wtforms import StringField, RadioField, BooleanField, IntegerField, SelectField, FileField
from wtforms.validators import DataRequired
#from werkzeug.utils import secure_filename


class newJobForm(FlaskForm):
	jobName = StringField('job name', validators=[DataRequired()])
	schema = FileField('schema', validators=[DataRequired()])