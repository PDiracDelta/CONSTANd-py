from flask_wtf import Form
from wtforms import StringField, TextAreaField, RadioField, BooleanField, IntegerField, DateTimeField, SelectField, FileField
from wtforms.validators import DataRequired
#from werkzeug.utils import secure_filename


class newJobForm(Form):
	jobName = StringField('job name', validators=[DataRequired()])
	schema = FileField('schema')