#!/usr/bin/env python2
import warnings
from math import log
from flask import Flask, render_template, request
from indCAPS import hammingBool, hamming, revComp, baseNumbers, deltaG, deltaS
from indCAPS import saltAdjusted, basicTemp, nearestNeighbor, estimateTM
from indCAPS import lastSharedBase, scanUnshared, scanSequence, evaluateSites
from indCAPS import generatePrimer, crisprEdit
from enzymeList import enzymes
from helperFuncs import checkBases, removeWhitespace, evaluateInput

# Set up Flask stuff
application = Flask(__name__)
#app.config.from_object(os.environ['APP_SETTINGS'])
application.jinja_env.trim_blocks = True
		

@application.route('/', methods=['GET','POST'])
def index():
	return(render_template('index.html'))
		
# I need to change this so it recognizes WHICH button was pushed rather than just getting a button push
@application.route('/results', methods=['GET','POST'])
def results():
	seq1 = False
	seq2 = False
	hamDist = False
	errors = []
	allResults = []
	notes = None
	if request.method == "POST":
		# Get stuff the user entered
		try:
			seq1 = request.form['seq1']
			seq2 = request.form['seq2']
			hamDist = int(request.form['ham'])
		except:
			errors.append("No sequences provided.")
		if seq1 and seq2 and hamDist:
			# Evaluate the input
			inputEvaluation = evaluateInput(seq1, seq2)
			seq1 = inputEvaluation[0]
			seq2 = inputEvaluation[1]
			notes = inputEvaluation[2]
			
			# Call function to evaluate enzymes
			for eachEnzyme in enzymes:
				enzymeName = eachEnzyme
				enzymeValue = enzymes[eachEnzyme]
				enzymeResults = evaluateSites(seq1,seq2,enzymeValue,hamDist, enzymeName)
				if enzymeResults is not None:
					allResults.append(enzymeResults)
			# Call function to generate primers
			# Right now, that information is returned by evaluateSites()
			# But, to be more portable, I need to edit it so that primer returns and site returns are handled by separate functions
			
			# Assemble output using sites and primers
			
			# Display output
		return(render_template('results.html',allResults=allResults,notes=notes))
	return(render_template('index.html')) # if you tried to go to the results page on your own rather than being sent by the index, redirect to the index page
	
@application.route('/screening', methods=['GET','POST'])
def screening():
	allResults = None
	return(render_template('results.html',allResults=allResults))

if __name__ == '__main__':
	application.run()