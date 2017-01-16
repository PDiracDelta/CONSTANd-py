#!/usr/bin/python
# -*- coding: utf-8 -*-
__author__ = 'Joris Van Houtven'
__maintainer__ = "Joris Van Houtven"
__email__ = "vanhoutvenjoris@gmail.com"

from flask import Flask
from flask_mail import Mail

app = Flask(__name__)
#app.config.from_object('config')
#db = SQLAlchemy(app)
mail = Mail(app)

from web import views

app.run(debug=False, host='0.0.0.0')