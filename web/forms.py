from flask_wtf import FlaskForm
#from flask_wtf.file import FileField
from wtforms import StringField, RadioField, BooleanField, IntegerField, SelectField, FileField, FieldList, FormField
from wtforms.validators import DataRequired
#from werkzeug.utils import secure_filename


class newJobForm(FlaskForm):
	jobName = StringField('job name', validators=[DataRequired()])
	schema = FileField('schema', validators=[DataRequired()])


class experimentForm(FlaskForm):
	data = FileField('data', validators=[DataRequired()])
	processingConfig = FileField('processingConfig')
	wrapper = FileField('wrapper')


class jobSettingsForm(FlaskForm):
	experiments = FieldList(FormField(experimentForm))
