from runweb import app

@app.route('/')
def hello_world():
    return 'Hello, Crazy World!'
