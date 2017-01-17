from flask import Flask
from flask_mail import Mail

app = Flask(__name__)
app.config.from_object('web.config')
#db = SQLAlchemy(app)
mailer = Mail(app)

from web import views