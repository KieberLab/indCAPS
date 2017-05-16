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
			allowMisMatch = bleach.clean(request.form['allowMM']) # This is now a checkbox, it will exist if it was checked and not exist if not
			sodiumConc = float(bleach.clean(request.form['sodiumConc']))
			primerConc = float(bleach.clean(request.form['allowMM']))
			ampliconLength = int(bleach.clean(request.form['ampliconLength']))
			primerType = bleach.clean(request.form['primerType'])
			primerLength = int(bleach.clean(request.form['primerLength']))
		except:
			errors.append("No sequences provided or Mismatch Match missing.")
			return(render_template('results.html',allResults=[],notes=errors))
		if seq1 and seq2:
			if helperFuncs.nonBasePresent(seq1) or helperFuncs.nonBasePresent(seq2):
				return(render_template('results.html',allResults=[],notes=["Non-bases included in input."]))
			
			# Define the settings object
			Settings = indCAPS.SettingsObject(TM=TM,ampliconLength=ampliconLength,primerType=primerType,primerLength=primerLength,allowMismatch=allowMisMatch,hammingThreshold=hamDist,organism=None,sodiumConc=sodiumConc,primerConc=primerConc*10**(-9))
			
			indCAPS.Settings = Settings
			helperFuncs.Settings = Settings
			
			# Evaluate the input
			inputEvaluation = helperFuncs.evaluateInput(seq1, seq2)
			seq1 = inputEvaluation[0]
			seq2 = inputEvaluation[1]
			notes = inputEvaluation[2]

			
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
			if allResults == [] or allResults == None:
				notes.append('No primer candidates. Please consider increasing the mismatch tolerance or altering your desired amplicon length.')
			return(render_template('results.html',allResults=allResults,notes=notes))
	return(render_template('index.html')) # if you tried to go to the results page on your own rather than being sent by the index, redirect to the index page
	
@application.route('/screening', methods=['GET','POST'])
def screening():
	# TODO: populate Settings object
	# TODO: send user to index if they try to go directly to the results page
	errors = []
	seq = ""
	hamDist = 0
	TM = 58
	cutoffPercent = 0
	organism = "Athaliana"
	if request.method == "POST":
		# Get stuff the user entered
		try:
			seq = bleach.clean(request.form['seq1'])
			targetSeq = bleach.clean(request.form['targetSeq'])
			hamDist = int(bleach.clean(request.form['ham']))
			TM = int(bleach.clean(request.form['tm']))
			allowMisMatch = bleach.clean(request.form['allowMM']) # This is now a checkbox, it will exist if it was checked and not exist if not
			sodiumConc = float(bleach.clean(request.form['sodiumConc']))
			primerConc = float(bleach.clean(request.form['allowMM']))
			ampliconLength = int(bleach.clean(request.form['ampliconLength']))
			primerType = bleach.clean(request.form['primerType'])
			primerLength = int(bleach.clean(request.form['primerLength']))
			organism = bleach.clean(request.form['organism'])
		except:
			errors.append("No sequences provided or Mismatch Match missing.")
			return(render_template('results.html',allResults=[],notes=errors))
		if seq and hamDist and TM:
			# Make settings object
			Settings = SettingsObject(TM=TM,ampliconLength=ampliconLength,primerType=primerType,primerLength=primerLength,allowMismatch=allowMisMatch,hammingThreshold=hamDist,organism=organism,sodiumConc=sodiumConc,primerConc=primerConc*10**(-9))
			
			# Populate modules with settings object
			indcaps.Settings = Settings
			helperFuncs.Settings = Settings
		
			# Call function to evaluate enzymes
			for eachEnzyme in enzymes:
				enzymeName = eachEnzyme
				enzymeValue = enzymes[eachEnzyme]
				mutationResults = indCAPS.evaluateMutations(seq,targetSeq,enzymeDict,enzymeName)
				if mutationResults is not None:
					allResults.append(enzymeResults)
			
			if allResults == [] or allResults == None:
				notes.append('No primer candidates. Please consider increasing the mismatch tolerance or altering your desired amplicon length.')
			return(render_template('screening.html',allResults=allResults,notes=notes))
	return(render_template('index.html')) # if you tried to go to the results page on your own rather than being sent by the index, redirect to the index page

if __name__ == '__main__':
	application.run()