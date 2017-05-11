#!/usr/bin/env python2
import warnings
from math import log
from flask import Flask, render_template, request
import indCAPS
import helperFuncs
from enzymeList import enzymes
import bleach

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
			allowMisMatch = bleach.clean(request.form['allowMM'])
		except:
			errors.append("No sequences provided or Mismatch Match missing.")
			return(render_template('results.html',allResults=[],notes=errors))
		if seq1 and seq2 and hamDist:
			if helperFuncs.nonBasePresent(seq1) or helperFuncs.nonBasePresent(seq2):
				return(render_template('results.html',allResults=[],notes=["Non-bases included in input."]))
			# Evaluate the input
			inputEvaluation = helperFuncs.evaluateInput(seq1, seq2)
			seq1 = inputEvaluation[0]
			seq2 = inputEvaluation[1]
			notes = inputEvaluation[2]
			
			# Define the settings object
			Settings = indCAPS.SettingsObject(TM=60,ampliconLength=100,primerType='tm',primerLength=30,allowMismatch=False,hammingThreshold=3,organism=None,sodiumConc=0.05,primerConc=50*10**(-9))
			
			indCAPS.Settings = Settings
			
			# Call function to evaluate enzymes
			for eachEnzyme in enzymes:
				enzymeName = eachEnzyme
				enzymeValue = enzymes[eachEnzyme]
				enzymeResults = indCAPS.evaluateSites(seq1,seq2,enzymeValue,enzymeName)
				if enzymeResults is not None:
					allResults.append(enzymeResults)
			# Call function to generate primers
			# Right now, that information is returned by evaluateSites()
			# But, to be more portable, I need to edit it so that primer returns and site returns are handled by separate functions
			
			# Assemble output using sites and primers
			
			
			# Display output
			print(allResults)
			print(notes)
			if allResults == [] or allResults == None:
				notes.append('No primer candidates. Please consider increasing the mismatch tolerance or altering your desired amplicon length.')
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
			mutationResults = None#indCAPS.evaluateMutations(seq,altSeqs,seqThreshold,enzymeDict,hammingThreshold,enzymeName,TM)
			return(render_template('results.html',allResults=mutationResults))
	return(render_template('index.html',allResults=[]))

if __name__ == '__main__':
	application.run()