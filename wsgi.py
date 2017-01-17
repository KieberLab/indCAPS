import os
import sys
import warnings
#import requests
from math import log
from flask import Flask, render_template, request
from indCAPS import hammingBool, hamming, revComp, baseNumbers, deltaG, deltaS
from indCAPS import saltAdjusted, basicTemp, nearestNeighbor, estimateTM
from indCAPS import lastSharedBase, scanUnshared, scanSequence, evaluateSites
from indCAPS import generatePrimer, compareEnzymes, makePrimer, enzymes
from enzymes import enzymes

# Set up Flask stuff
app = Flask(__name__)
#app.config.from_object(os.environ['APP_SETTINGS'])
app.jinja_env.trim_blocks = True


@app.route('/', methods=['GET','POST'])
def index():
	return(render_template('index.html'))
		
@app.route('/results', methods=['GET','POST'])
def results():
	seq1 = False
	seq2 = False
	hamDist = False
	errors = []
	allResults = []
	if request.method == "POST":
		# Get stuff the user entered
		try:
			seq1 = request.form['seq1']
			seq2 = request.form['seq2']
			hamDist = int(request.form['ham'])
		except:
			errors.append("No sequences provided.")
		if seq1 and seq2 and hamDist:
			# Call function to evaluate enzymes
			for eachEnzyme in enzymes:
				enzymeName = eachEnzyme
				enzymeValue = enzymes[eachEnzyme]
				enzymeResults = evaluateSites(seq1,seq2,enzymeValue,hamDist, enzymeName)
				if enzymeResults is not None:
					allResults.append(enzymeResults)
			# Call function to generate primers
			
			# Call function to generate output
			
			# Display output
	return(render_template('results.html',allResults=allResults))
	
if __name__ == '__main__':
	app.run()
