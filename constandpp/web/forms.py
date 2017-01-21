from flask_wtf import FlaskForm
#from flask_wtf.file import FileField
from wtforms import StringField, RadioField, BooleanField, IntegerField, SelectField, FileField, FieldList, FormField, FloatField
from wtforms.validators import DataRequired
#from werkzeug.utils import secure_filename


class newJobForm(FlaskForm):
	jobName = StringField('job name', validators=[DataRequired()])
	schema = FileField('schema', validators=[DataRequired()])


class experimentForm(FlaskForm):
	dataFile = FileField('dataFile', validators=[DataRequired()])
	processingConfig = FileField('processingConfig')
	wrapper = FileField('wrapper')


class jobSettingsForm(FlaskForm):
	experiments = FieldList(FormField(experimentForm))
	pept2protCombinationMethod = SelectField('pept2protCombinationMethod', choices=[('mean', 'mean'), ('median', 'median')], default=('mean', 'mean'))
	alpha = FloatField(default=0.05)
	FCThreshold = FloatField(default=1)
	labelVolcanoPlotAreas = StringField(default='[true, false, false, false]')
	minExpression_bool = BooleanField(default=True)
	fullExpression_bool = BooleanField(default=False)
	numDifferentials = IntegerField(default=10)