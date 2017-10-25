indCAPS (pronounced indy-caps)
A package for finding dCAPS primers suitable for sequences with insertions or deletions rather than a single nucleotide polymorphism.


## Structure ##
wsgi.py contains the entry point for the flask framework and the function calls which analyze user input.
indCAPS.py contains all the heavy lifting functions for analyzing user input.
helperFuncs.py contains several general functions for quality-checking user input before passing to main indCAPS functions.


## Installation ##
This package was written to be implemented on an OpenShift platform using the Gunicorn WSGI server and the flask microframework. If implemented on a standalone server, the Gunicorn recommendations are to run a reverse proxy server such as nginx.

Required external packages are:
flask
bleach
gunicorn

Each package can be installed using the pip package manager for Python.


## Usage ##

If the user wishes to run indCAPS as a local web-server, the application can be started by setting the FLASK_APP environment variable to "wsgi.py", and then issuing the command "flask run" at the command line while in the indCAPS directory. After flask has started, a web server will be running on the user's local machine and accessible at http://127.0.0.1:5000/.




