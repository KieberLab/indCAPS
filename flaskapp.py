import os
import sys
import warnings
#import requests
from math import log
import numpy
from flask import Flask, render_template, request
from dCAPS-new-v7 import hammingBool, hamming, revComp, baseNumbers, deltaG, deltaS
from dCAPS-new-v7 import saltAdjusted, basicTemp, nearestNeighbor, estimateTM
from dCAPS-new-v7 import lastSharedBase, scanUnshared, scanSequence, evaluateSites
from dCAPS-new-v7 import generatePrimer, compareEnzymes, makePrimer, enzymes

# Set up Flask stuff
app = Flask(__name__)
app.config.from_object(os.environ['APP_SETTINGS'])

# Import enzymes
enzymes = {}
with open("motifs.txt") as enzymeFile:
	for line in enzymeFile:
		line = line.rstrip()
		line = line.lower()
		if line[0] is "#":
			continue
		(enzyme,motif,forward,reverse) = line.split("\t")
		enzymes[enzyme] = [motif,forward,reverse]

@aap.route('/', methods=['GET','POST'])
def index():
	errors = []
	if request.method == "POST":
		# Get stuff the user entered
		try:
			seq1 = request.form['seq1','seq2','ham']
			seq2 = request.form['seq2']
			hamDist = request.form['ham']
		except:
			errors.append("No sequences provided.")
		allResults = []
		if seq1 and seq2 and hamDist:
			# Call function to evaluate enzymes
			for eachEnzyme in enzymes:
				enzymeName = eachEnzyme
				enzymeValue = enzymes[eachEnzyme]
				enzymeResults = evaluateSites(seq1,seq2,enzymeValue,hamDist, enzymeName)
				if enzymeResults is not None:
					allResults.append(enzymeReuslts)
			# Call function to generate primers
			
			# Call function to generate output
			
			# Display output
	return(render_template('index.html'))
	
if __name__ == '__main__':
	app.run()