#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This file merely creates and runs the Flask web app object.
"""

__author__ = 'Joris Van Houtven'
__maintainer__ = "Joris Van Houtven"
__email__ = "vanhoutvenjoris@gmail.com"


from constandpp.web import app

# host=0.0.0.0 makes the server publicly available.
app.run(debug=True, host='0.0.0.0')