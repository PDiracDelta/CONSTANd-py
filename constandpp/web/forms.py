from flask_wtf import FlaskForm
from wtforms import StringField, TextAreaField, RadioField, BooleanField, IntegerField, DateTimeField, SelectField, FileField
from wtforms.validators import DataRequired
#from werkzeug.utils import secure_filename


class schemaForm(FlaskForm):
    schema = FileField('schema', validators=[DataRequired()])