from flask_wtf import FlaskForm
from wtforms import StringField, TextAreaField, RadioField, BooleanField, IntegerField, DateTimeField, SelectField, FileField


class SchemaForm(FlaskForm):
    name = StringField('name', validators=[DataRequired()])