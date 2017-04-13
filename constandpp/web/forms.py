"""
Definition of all the wtforms used in the web app.
"""

from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, IntegerField, SelectField, FileField, FieldList, FormField, FloatField
from wtforms.validators import DataRequired


class newJobForm(FlaskForm):
	"""
	1: create a new job and upload a schema.
	"""
	jobName = StringField('job name', validators=[DataRequired()])
	schema = FileField('schema', validators=[DataRequired()])


class experimentForm(FlaskForm):
	"""
	2.1: upload the experiment files, config files and wrappers.
	"""
	dataFile = FileField('dataFile', validators=[DataRequired()])
	processingConfig = FileField('processingConfig')
	wrapper = FileField('wrapper')


class jobSettingsForm(FlaskForm):
	"""
	2: upload input files (see 2.1) and set the parameters for the job.
	"""
	experiments = FieldList(FormField(experimentForm))
	pept2protCombinationMethod = SelectField('pept2protCombinationMethod', choices=[('mean', 'mean'), ('median', 'median')], default=('mean', 'mean'))
	alpha = FloatField(default=0.05)
	FCThreshold = FloatField(default=1)
	labelVolcanoPlotAreas = StringField(default='[true, false, false, false]')
	minExpression_bool = BooleanField(default=True)
	fullExpression_bool = BooleanField(default=False)
	numDifferentials = IntegerField(default=10)
	delim_out = SelectField('delim_out', choices=[('\\t', 'tab'), (',', 'comma')], default=('\\t', 'tab'))
	mailRecipient = StringField(default="foo@bar.baz")