import os
import sys
import warnings
import requests
from math import log
from flask import Flask, render_template, request
from dCAPS-new-v7 import hammingBool, hamming, revComp, baseNumbers, deltaG, deltaS
from dCAPS-new-v7 import saltAdjusted, basicTemp, nearestNeighbor, estimateTM
from dCAPS-new-v7 import lastSharedBase, scanUnshared, scanSequence, evaluateSites
from dCAPS-new-v7 import generatePrimer, compareEnzymes, makePrimer

# Set up Flask stuff
app = Flask(__name__)
app.config.from_object(os.environ['APP_SETTINGS'])

@aap.route('/', methods=['GET','POST'])
def index():
	errors = []
	if request.method == "POST":
		# Get stuff the user entered
		try:
			seq1 = request.form['seq1']
			seq2 = request.form['seq2']
			receivedRequest = requests.get() # I think i need to look up the docs on this
		except:
			errors.append("No sequences provided.")
		if receivedRequest:	
			x=1
			# Call function to evaluate enzymes
			
			# Call function to generate primers
			
			# Call function to generate output
			
			# Display output
	return(render_template('index.html'))
	
if __name__ == '__main__':
	app.run()