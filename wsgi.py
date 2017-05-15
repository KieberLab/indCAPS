import warnings
from math import log
from flask import Flask, render_template, request
from indCAPS import hammingBool, hamming, revComp, baseNumbers, deltaG, deltaS
from indCAPS import saltAdjusted, basicTemp, nearestNeighbor, estimateTM
from indCAPS import lastSharedBase, scanUnshared, scanSequence, evaluateSites
from indCAPS import generatePrimer, crisprEdit
from enzymeList import enzymes
from helperFuncs import checkBases, removeWhitespace, evaluateInput, nonBasePresent
import bleach
import os

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
	notes = []
	if request.method == "POST":
		# Get stuff the user entered
		try:
			seq1 = bleach.clean(request.form['seq1'])
			seq2 = bleach.clean(request.form['seq2'])
			hamDist = int(bleach.clean(request.form['ham']))
			TM = int(bleach.clean(request.form['tm']))
		except:
			errors.append("No sequences provided or Mismatch Match missing.")
			return(render_template('results.html',allResults=[],notes=errors))
		if seq1 and seq2 and hamDist:
			if nonBasePresent(seq1) or nonBasePresent(seq2):
				return(render_template('results.html',allResults=[],notes=["Non-bases included in input."]))
			# Evaluate the input
			inputEvaluation = evaluateInput(seq1, seq2)
			seq1 = inputEvaluation[0]
			seq2 = inputEvaluation[1]
			notes = inputEvaluation[2]
			
			# Call function to evaluate enzymes
			for eachEnzyme in enzymes:
				enzymeName = eachEnzyme
				enzymeValue = enzymes[eachEnzyme]
				enzymeResults = evaluateSites(seq1,seq2,enzymeValue,hamDist,enzymeName,TM)
				if enzymeResults is not None:
					allResults.append(enzymeResults)
			# Call function to generate primers
			# Right now, that information is returned by evaluateSites()
			# But, to be more portable, I need to edit it so that primer returns and site returns are handled by separate functions
			
			# Assemble output using sites and primers
			
			# Display output
			print(allResults)
			print(notes)
			return(render_template('results.html',allResults=allResults,notes=notes))
	return(render_template('index.html')) # if you tried to go to the results page on your own rather than being sent by the index, redirect to the index page
	
@application.route('/screening', methods=['GET','POST'])
def screening():
	errors = []
	seq = ""
	hamDist = 0
	TM = 58
	cutoffPercent = 0
	organism = "Athaliana"
	if request.method == "POST":
		# Get stuff the user entered
		try:
			seq = bleach.clean(request.form['seq'])
			hamDist = int(bleach.clean(request.form['ham']))
			cutoffPercent = int(bleach.clean(request.form['cutoffPercent']))
			organism = int(bleach.clean(request.form['org']))
			TM = int(bleach.clean(request.form['tm']))
		except:
			errors.append("No sequences provided or Mismatch Match missing.")
			return(render_template('results.html',allResults=[],notes=errors))
		if seq and hamDist and TM:
			mutationResults = None#evaluateMutations(seq,altSeqs,seqThreshold,enzymeDict,hammingThreshold,enzymeName,TM)
			return(render_template('results.html',allResults=mutationResults))
	return(render_template('index.html',allResults=[]))

if __name__ == "__main__":
	application.run()