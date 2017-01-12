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

@app.route('/')
def hello_world():
    return 'Hello, World!'

app.run(debug=True)