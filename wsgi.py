import warnings
from math import log
from flask import Flask, render_template, request # 3 clause BSD
import indCAPS
import helperFuncs
from enzymeList import enzymes
import bleach # apache v2

# Set up Flask stuff
application = Flask(__name__)
application.jinja_env.trim_blocks = True

@application.route('/', methods=['GET','POST'])
def index():
	return(render_template('index.html'))
		
@application.route('/results', methods=['GET','POST'])
def results():
	seq1 = False
	seq2 = False
	hamDist = False
	errors = []
	allResults = []
	notes = []
	if request.method == "POST":
		# Try to get infomration the user entered
		try:
			seq1 = bleach.clean(request.form['seq1'])
			seq2 = bleach.clean(request.form['seq2'])
			hamDist = int(bleach.clean(request.form['ham']))
			TM = int(bleach.clean(request.form['tm']))
			allowMisMatch = bleach.clean(request.form['allowMM']) # This is now a checkbox, it will exist if it was checked and not exist if not
			sodiumConc = float(bleach.clean(request.form['sodiumConc']))
			primerConc = float(bleach.clean(request.form['primerConc']))
			ampliconLength = int(bleach.clean(request.form['ampliconLength']))
			primerType = bleach.clean(request.form['primerType'])
			primerLength = int(bleach.clean(request.form['primerLength']))
		# Throw an error if it didn't work
		except Exception as e:
			errors.append("No sequences provided or Mismatch Match missing.")
			errors.append(e)
			return(render_template('results.html',allResults=[],notes=errors))
		
		# If the user gave two sequences, check them for quality and run main routine
		if seq1 and seq2:
			if helperFuncs.nonBasePresent(seq1) or helperFuncs.nonBasePresent(seq2):
				return(render_template('results.html',allResults=[],notes=["Non-bases included in input."]))
			
			# Define the settings object
			Settings = indCAPS.SettingsObject(TM=TM,ampliconLength=ampliconLength,primerType=primerType,primerLength=primerLength,allowMismatch=allowMisMatch,hammingThreshold=hamDist,organism=None,sodiumConc=sodiumConc,primerConc=primerConc*10**(-9))
			
			# Pass the new Settings object down to the namespaces for the custom functions
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
			# To be more portable, future versions might make it so that site returns are separate from primer returns
			
			# Assemble output using sites and primers, then display
			if allResults == [] or allResults == None:
				notes.append('No primer candidates. Please consider increasing the mismatch tolerance or altering your desired amplicon length.')
			return(render_template('results.html',allResults=allResults,notes=notes))
	return(render_template('index.html')) # if you tried to go to the results page on your own rather than being sent by the index, redirect to the index page
	
@application.route('/screening', methods=['GET','POST'])
def screening():
	errors = []
	allResults = []
	notes = []
	if request.method == "POST":
		# Get stuff the user entered
		try:
			seq = bleach.clean(request.form['seq'])
			targetSeq = bleach.clean(request.form['targetSeq'])
			hamDist = int(bleach.clean(request.form['ham']))
			TM = int(bleach.clean(request.form['tm']))
			allowMisMatch = bleach.clean(request.form['allowMM']) # This is now a checkbox, it will exist if it was checked and not exist if not
			sodiumConc = float(bleach.clean(request.form['sodiumConc']))
			primerConc = float(bleach.clean(request.form['primerConc']))
			ampliconLength = int(bleach.clean(request.form['ampliconLength']))
			primerType = bleach.clean(request.form['primerType'])
			primerLength = int(bleach.clean(request.form['primerLength']))
			organism = None #bleach.clean(request.form['organism'])
			seqThreshold = float(bleach.clean(request.form['cutoffPercent']))
		except Exception as e:
			errors.append("Error in input form.")
			errors.append(e)
			return(render_template('results.html',allResults=[],notes=errors))
		if seq and hamDist and TM:
			# Make settings object
			Settings = indCAPS.SettingsObject(TM=TM,ampliconLength=ampliconLength,primerType=primerType,primerLength=primerLength,allowMismatch=allowMisMatch,hammingThreshold=hamDist,organism=organism,sodiumConc=sodiumConc,primerConc=primerConc*10**(-9),seqThreshold=seqThreshold)
			
			# Populate modules with settings object
			indCAPS.Settings = Settings
			helperFuncs.Settings = Settings
			
			# Evaluate the input
			inputEvaluation = helperFuncs.evaluateInput(seq)
			seq1 = inputEvaluation[0]
			seq2 = inputEvaluation[1]
			notes = inputEvaluation[2]
			siteMatches = helperFuncs.checkSingleSite(seq,targetSeq)
			if siteMatches == 0:
				# No matches present
				notes.append('Target does not match any region in the sequence. Please provide a new target.')
				return(render_template('screening.html',allResults=None,notes=notes))
			elif siteMatches > 1:
				# Too many matches present
				notes.append('Target matches multiple locations in provided sequence. Primers cannot be designed. Please choose a different target sequence.')
				return(render_template('screening.html',allResults=None,notes=notes))
				
			# Call function to evaluate enzymes
			for eachEnzyme in enzymes:
				enzymeName = eachEnzyme
				enzymeValue = enzymes[eachEnzyme]
				mutationResults = indCAPS.evaluateMutations(seq,targetSeq,enzymeValue,enzymeName)
				if mutationResults is not None:
					allResults.append(mutationResults)
			
			# Tell the user if the program failed
			if allResults == [] or allResults == None:
				notes.append('No primer candidates. Please consider increasing the mismatch tolerance or altering your desired amplicon length.')
				
			# Render the results, send the user to the results page
			return(render_template('screening.html',allResults=allResults,notes=notes))
	# If you tried to go to the results page on your own rather than being sent by the index, redirect the user to the index page
	return(render_template('index.html')) 

# Isogenic screening
@application.route('/isogenic', methods=['GET','POST'])
def isogenic():
	wtSeq = False
	mutSeq = False
	targetSeq = False
	errors = []
	allResults = []
	notes = []
	if request.method == "POST":
		# Get stuff the user entered
		try:
			wtSeq = bleach.clean(request.form['wtSeq'])
			mutSeq = bleach.clean(request.form['mutSeq'])
			targetSeq = bleach.clean(request.form['targetSeq'])
			hamDist = int(bleach.clean(request.form['ham']))
			TM = int(bleach.clean(request.form['tm']))
			allowMisMatch = bleach.clean(request.form['allowMM']) # This is now a checkbox, it will exist if it was checked and not exist if not
			sodiumConc = float(bleach.clean(request.form['sodiumConc']))
			primerConc = float(bleach.clean(request.form['primerConc']))
			ampliconLength = int(bleach.clean(request.form['ampliconLength']))
			primerType = bleach.clean(request.form['primerType'])
			primerLength = int(bleach.clean(request.form['primerLength']))
			organism = None #bleach.clean(request.form['organism'])
			seqThreshold = float(bleach.clean(request.form['cutoffPercent']))
		except Exception as e:
			errors.append("Error in input form: No sequences provided.")
			errors.append(e)
			return(render_template('isogenic.html',allResults=[],notes=errors))
		if wtSeq and mutSeq and targetSeq:
			if helperFuncs.nonBasePresent(wtSeq) or helperFuncs.nonBasePresent(mutSeq) or helperFuncs.nonBasePresent(targetSeq):
				return(render_template('results.html',allResults=[],notes=["Non-bases included in input."]))
			# Make settings object
			Settings = indCAPS.SettingsObject(TM=TM,ampliconLength=ampliconLength,primerType=primerType,primerLength=primerLength,allowMismatch=allowMisMatch,hammingThreshold=hamDist,organism=organism,sodiumConc=sodiumConc,primerConc=primerConc*10**(-9),seqThreshold=seqThreshold)
			
			# Populate modules with settings object
			indCAPS.Settings = Settings
			helperFuncs.Settings = Settings
			
			# Evaluate the input
			inputEvaluation = helperFuncs.evaluateInput(wtSeq)
			wtSeq = inputEvaluation[0]
			notes = inputEvaluation[2]
			siteMatches = helperFuncs.checkSingleSite(wtSeq,targetSeq)
			
			if siteMatches == 0:
				# No matches present
				notes.append("Target sequence does not match any region in the wild-type sequence. Please provide a new target.")
				return(render_template('isogenic.html',allResults=None,notes=notes))
			elif siteMatches > 1:
				# Too many matches present
				notes.append("Target sequence matches multiple locations in provided wild-type sequence. Primers cannot be designed. Please choose a different target sequence.")
				return(render_template('isogenic.html',allResults=None,notes=notes))
			
			# Call function to evaluate enzymes
			for eachEnzyme in enzymes:
				enzymeName = eachEnzyme
				enzymeValue = enzymes[eachEnzyme]
				isogenicResults = indCAPS.evaluateIsogenic(wtSeq,mutSeq,targetSeq,enzymeValue,enzymeName)
				if isogenicResults is not None:
					allResults.append(isogenicResults)
			
			# Tell the user if the program failed
			if allResults == [] or allResults == None:
				notes.append('No primer candidates. Please consider increasing the mismatch tolerance or altering your desired amplicon length.')

			# Render the results, send the user to the results page
			return(render_template('isogenic.html',allResults=allResults,notes=notes))
	# If you tried to go to the results page on your own rather than being sent by the index, redirect the user to the index page
	return(render_template('index.html'))

# Gotta show up in Google
@application.route('/google61adb2a906ee8c51.html')
def site_verification():
	return(render_template("google61adb2a906ee8c51.html"))

# Web server magic
if __name__ == '__main__':
	application.run()