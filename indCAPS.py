#!/usr/bin/env python2
"""@package indCAPS
Examine multiple strings representing DNA bases for cases where 
diagnostic PCR primers can be constructed to distinguish the strings.

Currently tailored for use as a web application, but can be adapted
for standalone use.
"""
# -*- coding: utf-8 -*-
from math import log
import itertools
import warnings

## Enzymes
# Will not include nicking enzymes or intron-encoding enzymes or homing enzymes
# Also does not include two-site cutters: BaeI, BsaXI, CspCI, BcgI
from enzymeList import enzymes

## Custom Classes
class SettingsObject:
	"""
	This object should be in the global scope and hold all the 
	invariant information on things like parameters of primer design.
	When calling functions from the indCAPS package, make sure the first
	thing you do is set indCAPS.Settings = X and helperFuncs.Settings = X, 
	where X is the name of your SettingsObject() object.
	
	Nearly all functions in this package depend on this object existing.
	
	TM = float, desired melting temperature of primer in degrees C
	ampliconLength = integer, length of dCAPS amplicon
	primerType = string, choice of 'tm' or 'length'
	primerLength = integer, length of primer to design
	allowMismatch = boolean, whether to allow 3' mismatch in primrs
	hammingThreshold = integer, mismatch threshold
	organism = string, species name
	sodiumConc = float, units milliMolar
	primerConc = float, units nanoMolar
	seqThreshold = float, range from 0-100, cutoff percent for CRISPR screening
	"""
	def __init__(self,TM=55,ampliconLength=70,primerType='tm',primerLength=25,allowMismatch=False,hammingThreshold=1,organism=None,sodiumConc=0.05,primerConc=50*10**(-9),seqThreshold=0):
		self.TM = TM
		self.ampliconLength = ampliconLength
		self.primerType = primerType
		self.primerLength = primerLength
		if allowMismatch == 'yes':
			self.allowMismatch = True
		elif allowMismatch == 'no':
			self.allowMismatch = False
		else:
			self.allowMismatch = False
		self.hammingThreshold = hammingThreshold
		self.organism = organism
		self.sodiumConc = sodiumConc
		self.primerConc = primerConc
		self.seqThreshold = seqThreshold
		
		commonEdits =  {"Athaliana":[],
					"Osativa":[],
					"Hsapiens":[],
					"Mmusculus":[]}
		if self.organism in commonEdits:
			self.organismEditProfile = commonEdits[organism]
		else:
			self.organismEditProfile = [1,-1]

## Custom Functions
def proportionalDistance(mutMotif,motif):
	"""
	Computes a hamming distance between two sequences. Returns a list of 
	lists. Returned values are:
	[number of comparisons, proportion of matching comparisons, lowHamDist, highHamDist]
	
	This is specifically intended to match a hypothetical, defined DNA sequence against
	a restriction site motif. The only nonspecific base in the mutant motif should be
	the X characters introduced to represent insertions.
	
	Example: The sequence "GCAA" compared to "GCXA", representing a single
	base pair insertion, has the possible values "GCTA", "GCAA", "GCGA", "GCTA".
	Of the four possible comparisons, one is a match, and this would be returned as [4,0.25]
	"""
	# Mut will have X characters.
	# Returns assessment of most optimistic match and also proportion of all matches that it will match
	possibleMatches = {'a':['a','r','w','m','d','h','v','n','x'],
					   'g':['g','r','s','k','b','d','v','n','x'],
					   'c':['c','y','s','m','b','h','v','n','x'],
					   't':['t','y','w','k','b','d','h','n','x'],
					   'r':['a','g','n','x'],
					   'y':['c','t','n','x'],
					   's':['g','c','n','x'],
					   'w':['a','t','n','x'],
					   'k':['g','t','n','x'],
					   'm':['a','c','n','x'],
					   'b':['c','g','t','n','x'],
					   'd':['a','g','t','n','x'],
					   'h':['a','c','t','n','x'],
					   'v':['a','c','g','n','x'],
					   'n':['a','g','c','t','r','y','s','w','k','m','b','d','h','v','n','x'],
					   'x':['a','g','c','t','r','y','s','w','k','m','b','d','h','v','n','x']}
	
	# Honestly I don't know if this is necessary but it doesn't hurt.
	mutMotif = mutMotif.lower()
	motif = motif.lower()
	
	# Arrange possible comparisons of sequence and motif
	sets = [ [mutMotif,motif] , [revComp(mutMotif).lower(),motif] , [mutMotif, revComp(motif).lower()] ]
	
	setOutput = []
	
	for eachNumber in range(0,len(sets)):
		# Set up initial data
		eachSet = sets[eachNumber]
		currentMut = eachSet[0]
		currentMotif = eachSet[1]
		
		# Set up hamming distance counter
		currentHam = 0
		
		# Proportion of possible comparisons which are matches
		currentProp = 1
		
		# Number of possible comparisons (one if no X in sequence, multiple of 4 if X in sequence)
		comparisons = 0
		
		# Iterate over bases
		for eachPos in range(0,len(currentMut)):
			mutBase = currentMut[eachPos]
			motifBase = currentMotif[eachPos]
			
			# Check to see if there's a match at all
			if mutBase not in possibleMatches[motifBase]:
				currentHam += 1
			
			# If there's an X, increase comparison counter and reduce proportion variable
			if mutBase == "x":
				comparisons += 4
				if motifBase in ['g','c','a','t']:
					currentProp = currentProp * 1/4
				
				elif motifBase in ['r','y','s','w','k','m']:
					currentProp = currentProp * 1/2
				
				elif motifBase in ['b','d','h','v']:
					currentProp = currentProp * 3/4
		
		# If no Xs in sequence, then only one comparison can be made
		if comparisons == 0:
			comparisons = 1
		
		# Find worst case hamming distance and construct output
		hamDistHigh = currentHam + len([base for base in currentMut if base == "x"])
		setOutput.append([comparisons,currentProp,currentHam,hamDistHigh])
	return(setOutput)
	
def getDegenerateMatch(base):
	"""
	Given a degenerate base, returns all possible matches.
	Possible return values - a, g, c, t, n, x
	base = string, lowercase, degenerate base symbol
	"""
	base = base.lower()
	possibleMatches = {
				   'r':['a','g','n','x'],
				   'y':['c','t','n','x'],
				   's':['g','c','n','x'],
				   'w':['a','t','n','x'],
				   'k':['g','t','n','x'],
				   'm':['a','c','n','x'],
				   'b':['c','g','t','n','x'],
				   'd':['a','g','t','n','x'],
				   'h':['a','c','t','n','x'],
				   'v':['a','c','g','n','x']}
	
	# Return just the string, not the string in a list
	replacementBase = possibleMatches[base][0]
	return(replacementBase)

def hamming(seq1,seq2,allResults=False,allComparisons=True):
	"""
	Calculates the hamming distance between two equal-length sequences.
	Expects a character string, will complain if it doesn't get it.
	Returns a list of integers.
	
	Warning: Be careful how you set up sequences with N's and call this.
	Ex: hamming('gggnnn','grcnnn',allResults=True,allComparisons=True) 
		returns [4,0,0].
	
	seq1 = string
	seq2 = string
	allResults = True or False
	allComparisons = True or False
	
	If allComparisons is False, only the first output result (sequences
	compared as provided) is returned.
	If allComparisons is True, all of the output results are returned.
	If allResults is False, only the minimum hamming distance is returned.
	If allResults is True, all three hamming distances are returned.
	
	Ex: [4,0,4]
	allComparisons = True, allResults = True returns [4,0,4]
	allComparisons = False, allResults = True returns [4]
	allComparisons = True, allResults = False returns [0]
	allComparisons = False, allResults = False returns [4]
	"""
	seq1 = seq1.lower()
	seq2 = seq2.lower()
	
	# FIXME: it sends a program-killing error if I try to send a warning right now.
	if len(seq1) != len(seq2):
		warnings.warn("Comparing sequences of unequal length.")
		
	# Originally I had some fancy-pants list comprehension approach put in here but I abandoned it to have more control over the process. Keeping it around because I kinda like the elegance and might want to use it in a different project.
	# distance1 = len(seq1) - sum([x == y for (x,y) in zip(seq1,seq2)])
	# These next two are for cases of non-palindromic motifs
	# distance2 = len(seq1) - sum([x == y for (x,y) in zip(seq1,revComp(seq2))])
	# distance3 = len(seq1) - sum([x == y for (x,y) in zip(revComp(seq1),seq2)])
	# return(min(distance1,distance2,distance3))
	
	# Define permissible matches.
	# Follows IUPAC base conventions with the exception of X, which implies an inserted base following a CRISPR/Cas9 editing event. X is treated here as an N.
	possibleMatches = {'a':['a','r','w','m','d','h','v','n','x'],
				  'g':['g','r','s','k','b','d','v','n','x'],
				  'c':['c','y','s','m','b','h','v','n','x'],
				  't':['t','y','w','k','b','d','h','n','x'],
				  'r':['a','g','n','x'],
				  'y':['c','t','n','x'],
				  's':['g','c','n','x'],
				  'w':['a','t','n','x'],
				  'k':['g','t','n','x'],
				  'm':['a','c','n','x'],
				  'b':['c','g','t','n','x'],
				  'd':['a','g','t','n','x'],
				  'h':['a','c','t','n','x'],
				  'v':['a','c','g','n','x'],
				  'n':['a','g','c','t','r','y','s','w','k','m','b','d','h','v','n','x'],
				  'x':['a','g','c','t','r','y','s','w','k','m','b','d','h','v','n','x']}
	# NOTE: degenerate bases should only occur in enzyme recognition motifs
	# This should allow motifs with degenerate bases to match given bases in a genomic sequence
	# However, if this were coded for general use, it would produce incorrect results
	# eg, 'r' is returned here as a valid match for 'n', when that is not strictly true
	
	# Construct sets of sequences for comparison purposes
	sets = [ [seq1,seq2] , [revComp(seq1).lower(),seq2] , [seq1, revComp(seq2).lower()] ]
	
	setOutput = [0] * len(sets)
	
	# Iterate over each set checking its matches
	for eachSet in range(0,len(sets)):
		currentSet = sets[eachSet]
		currentSeq1 = currentSet[0]
		currentSeq2 = currentSet[1]
		
		for eachPosition in range(0,len(currentSeq1)):
			seq1Base = currentSeq1[eachPosition]
			seq2Base = currentSeq2[eachPosition]
			lookup = possibleMatches[seq1Base]
			if seq2Base not in lookup:
				setOutput[eachSet] += 1
	if not allComparisons:
		setOutput = [setOutput[0]]
		
	if allResults:
		return(setOutput)
	else:
		return(min(setOutput))


def revComp(seq1):
	"""
	Takes a string with valid DNA bases (including degenerate bases) 
	and returns the reverse complement.
	
	Input
	seq1 = string of valid DNA bases
	
	Output
	seqString = string of valid DNA bases
	"""
	seq1 = seq1.lower()
	seqNew = []
	baseDict = {'g':'c',
			 'c':'g',
			 'a':'t',
			 't':'a',
			 'n':'n',
			 'y':'r',
			 'r':'y',
			 'w':'w',
			 's':'s',
			 'k':'m',
			 'm':'k',
			 'd':'h',
			 'v':'b',
			 'h':'d',
			 'b':'v',
			 'x':'x'}
	for each in seq1:
		seqNew.append(baseDict[each])
	seqNew.reverse()
	seqString = "".join(seqNew)
	seqString = seqString.upper()
	return(seqString)

## Melting temperature calculations
def baseNumbers(seq):
	"""
	Counts the number of each base in the string. Currently only accepts g,c,a,t bases.
	
	Input
	seq = string of valid DNA bases
	
	Output
	storage = dict with lowercase base letters as keys and count of each in seq as values.
	"""
	# Expects a string
	# Returns a dict with the counts for each base
	seq = seq.lower()
	storage = {'g':0,'c':0,'a':0,'t':0}
	for eachBase in seq:
		storage[eachBase] += 1
	return(storage)

#oligo		dH	dS	dG
#dAA/dTT	-8.0	-21.9	-1.2
#dAT/dAT	-5.6	-15.2	-0.9
#dCG/dCG	-11.8	-29.0	-2.8
#dCT/dAG	-6.6	-16.4	-1.5
#dGA/dTC	-8.8	-23.5	-1.5
#dGC/dGC	-10.5	-26.4	-2.3
#dGG/dCC	-10.9	-28.4	-2.1
#dGT/dAC	-9.4	-25.5	-1.5
#dTA/dTA	-6.6	-18.4	-0.9
#dTG/dCA	-8.2	-21.0	-1.7
#initiation	0.6	-9.0	3.4
#self-comp	0.0	-1.4	0.4
#non-comp	0.0	0.0	0.0
# from 1. Sugimoto, N., Nakano, S. I., Yoneyama, M. and Honda, K. I. Improved thermodynamic parameters and helix initiation factor to predict stability of DNA duplexes. Nucleic Acids Res. 24, 4501–4505 (1996).

def deltaG(seq):
	"""
	Calcuating free energy of a sequence for use in primer melting temperature calculation.
	"""
	seq = seq.lower()
	dG = {'aa':-1.2,
		'tt':-1.2,
		'at':-0.9,
		'cg':-2.8,
		'ct':-1.5,
		'ag':-1.5,
		'ga':-1.5,
		'tc':-1.5,
		'gc':-2.3,
		'gg':-2.1,
		'cc':-2.1,
		'gt':-1.5,
		'ac':-1.5,
		'ta':-0.9,
		'tg':-1.7,
		'ca':-1.7}
	dGout = 3.4
	for eachBase in range(0,len(seq)-1):
		currentDimer = seq[eachBase:eachBase+2]
		dGout += dG[currentDimer]
	if seq == revComp(seq):
		dGout += 0.4
	return(dGout)


def deltaH(seq):
	"""
	Calculates enthalpy of a sequence for use in primer melting temperature calculation.
	Units kcal/mole
	"""
	seq = seq.lower()
	dH = {'aa':-8.0,
		'tt':-8.0,
		'at':-5.6,
		'cg':-11.8,
		'ct':-6.6,
		'ag':-6.6,
		'ga':-8.8,
		'tc':-8.8,
		'gc':-10.5,
		'gg':-10.9,
		'cc':-10.9,
		'gt':-9.4,
		'ac':-9.4,
		'ta':-6.6,
		'tg':-8.2,
		'ca':-8.2}
	dHout = 0.6
	for eachBase in range(0,len(seq)-1):
		currentDimer = seq[eachBase:eachBase+2]
		dHout += dH[currentDimer]
	if seq == revComp(seq):
		dHout += 0
	return(dHout)
		
def deltaS(seq):
	"""
	Calculates entropy of a sequence for use in primer melting temperature calculation.
	Units cal/(mole*K)
	"""
	# Units cal/(mole*K)
	seq = seq.lower()
	dS = {'aa':-21.9,
		'tt':-21.9,
		'at':-15.2,
		'cg':-29.0,
		'ct':-16.4,
		'ag':-16.4,
		'ga':-23.5,
		'tc':-23.5,
		'gc':-26.4,
		'gg':-28.4,
		'cc':-28.4,
		'gt':-25.5,
		'ac':-25.5,
		'ta':-28.4,
		'tg':-21.0,
		'ca':-21.0}
	dSout = -9.0
	for eachBase in range(0,len(seq)-1):
		currentDimer = seq[eachBase:eachBase+2]
		dSout += dS[currentDimer]
	if seq == revComp(seq):
		dSout += -1.4
	return(dSout)

def saltAdjusted(seq):
	"""
	Calculates melting temperature for an oligonucleotide using a
	salt-adjusted calculation reportedly accurate for long sequences
	"""
	bases = baseNumbers(seq)
	temp = 81.5 + (41*((bases['g']+bases['c'])/(bases['a']+bases['c']+bases['g']+bases['t']))-(500/(bases['a']+bases['c']+bases['g']+bases['t']))+16.6*log(Settings.sodiumConc,10))
	return(temp)

def basicTemp(seq):
	"""
	Calculates melting temperature for an oligonucleotide using a
	basic calculation appropriate for very short sequence lengths.
	"""
	bases = baseNumbers(seq) # dict of 'g','c','a','t'
	temp = (bases['a']+bases['t'])*2+(bases['g']+bases['t'])*4
	return(temp)

def nearestNeighbor(seq):
	"""
	Calculates DNA sequence melting temperature using the nearest neighbor method.
	"""
	# typical salt conc is [50 mM], keep it between [0.01 M] and [1 M] 
	# typical primer conc is 50 nM
	# keep it above 8 bases
	# From http://biotools.nubic.northwestern.edu/OligoCalc.html and cited sources
	R = 1.9872036 # kcal / (mole*kelvin)
	T = (1000*(-deltaH(seq))) / ((-deltaS(seq))+R*log(1/Settings.primerConc))-272.9+16.6*log(Settings.sodiumConc,10) # FIXME: i think this is messed up somewhere? dont i need a -3.4kcal/kmole? it's not drastically wrong i dont think, just slightly messed up. oh wait i think i'm fine, the 3.4 is included in the deltaH() function
	return(T)

def estimateTM(seq,func='nearestNeighbor'):
	"""
	Estimates primer melting temperature using one of several methods.
	Available methods:
	Nearest Neighbor
	Salt Adjusted
	Basic
	
	For now, only nearest neighbor is supported.
	"""
	# Default function is 'nearestNeighbor'
	# see http://biotools.nubic.northwestern.edu/OligoCalc.html
	# Find number of each base
	
	if func is 'nearestNeighbor':
		return(nearestNeighbor(seq))
	# If length between 8 and 40 bp:
	# use nearest neighbor
	
	# If length less than 8:
	
	
	# If length greater than 40:
	# use salt adjusted
	return(None)

def lastSharedBase(seq1,seq2,direction='left'):
	"""
	Examines two sequences, determines the position of the last 
	shared base between the two in the indicated direction.
	Intended usage is that the provided sequences are identical 
	on their left and right flanks with an unshared region in the middle.
	
	seq1 = string, DNA bases
	seq2 = string, DNA bases
	direction = string, "left" or "right"
	"""
	# Compares two strings, finds longest stretch of identical sequence from indicated direction
	# Sequences should be identical length
	# seq1 and seq2 are strings, direction is string of choices "left" and "right". "left" is assumed to be the default
	# returns position of last shared base from chosen direction
	# The returned position is always relative to the original provided orientation
	# Check that position string is correct
	if direction.lower() not in ["right","left"]:
		warnings.warn("Direction not specified correctly. Please use 'right' or 'left'. Using default of left.")
		direction = 'left'
	
	# Reverse strings if necessary
	if direction.lower() == "right":
		seq1 = revComp(seq1)
		seq2 = revComp(seq2)
	
	# Set variables for the loop
	exitVar = 0
	position = 0
	
	# Start the loop to compare the strings
	while exitVar is 0:
		# Compare the sequences at that position
		if seq1[position] == seq2[position]:
			position = position + 1
		else:
			exitVar = 1

		# Exit if we've exhausted the whole string
		if position > min(len(seq1),len(seq2)):
			exitVar = 1
	
	# Return the position in the original context - we revcomp'd, now we need to subtract from end
	if direction is "right":
		position = 0 - (position )

	return(position)
	
def scanUnshared(seq,currentMotif,lastShared,lastSharedReverse):
	"""
	Examines the unshared region of the given sequence 
	for matches below the given hamming threshold.
	
	seq = string, valid DNA bases
	currentMotif = string, enzyme recognition motif
	lastShared = integer, index of last shared based from left
	lastSharedReverse = integer, index of last shared base from right
	hammingThreshold = integer
	"""
	# if lastShared = 15, motifLen = 4,
	# range is [15 - 4 + 1=12 : 15 + 4 - 1 = 18]
	# border is [motifLen-1]=3, or:
	# [0, 1, 2, 3, 4, 5, 6, 7...]
	# [12, 13, 14, 15, 16...]
	# border = 3, or position 3, is base 15, the lastShared base
	
	motifLen = len(currentMotif)
	border = motifLen - 1
	seqUnshared = seq[lastShared - motifLen + 1 : lastSharedReverse + motifLen - 1]
	seqUSpos = []
	seqUSexact = []

	# Scan across seq1 unshared
	for eachPosition in range(0,len(seqUnshared)-motifLen+1):
		
		if eachPosition < border:
			partSeqLeft = seqUnshared[eachPosition:border]
			partSeqRight = seqUnshared[border:eachPosition+motifLen]
			
			partMotifLeft = currentMotif[:len(partSeqLeft)]
			partMotifRight = currentMotif[len(partSeqLeft):]	
			
			currentSeqLeft = partSeqLeft + "n"*(eachPosition+motifLen-border)
			currentSeqRight = "n"*(border-eachPosition) + partSeqRight
			currentMotifLeft = partMotifLeft + "n"*(len(currentMotif) - len(partSeqLeft))
			currentMotifRight = "n"*len(partSeqLeft) + partMotifRight
			
			leftHam = hamming(currentMotifLeft,currentSeqLeft,allComparisons=False)
			rightHam = hamming(currentMotifRight,currentSeqRight,allComparisons=False)
			
			if leftHam == 0 and rightHam == 0:
				seqUSpos.append(eachPosition+lastShared - motifLen + 1)
				seqUSexact.append(eachPosition+lastShared - motifLen + 1)
			elif leftHam <= Settings.hammingThreshold and rightHam == 0:
				seqUSpos.append(eachPosition+lastShared - motifLen + 1)

		else:
			currentSeq1 = seqUnshared[eachPosition:eachPosition+motifLen]
			currentHam = hamming(currentSeq1,currentMotif)
			
			if currentHam == 0:
				seqUSpos.append(eachPosition+lastShared - motifLen + 1)
				seqUSexact.append(eachPosition+lastShared - motifLen + 1)
	return([seqUSpos,seqUSexact])

def scanSequence(seq1,seq2,currentMotif,direction):
	"""
	This is going to scan the sequence from left to right, 
	with the assumption that left is the 5' side of the primer 
	and anything right is the downstream direction. If 'right' 
	is specified as the direction, the sequence is reversed 
	before processing.
	
	seq1 = string, valid DNA bases
	seq2 = string, valid DNA bases
	currentMotif = string, restriction enzyme motif
	direction = string, choose from "left" or "right"
	
	Output
	
	"""
	motifLen = len(currentMotif)
	lastShared = lastSharedBase(seq1,seq2,direction)
	
	possibleDirections = ["left","right"]
	directionComparison = [x in direction for x in possibleDirections]
	otherDirection = [d for d,s in zip(possibleDirections,directionComparison) if not s][0]
	lastSharedReverse = lastSharedBase(seq1,seq2,otherDirection)
	
	if direction is 'right':
		seq1 = revComp(seq1)
		seq2 = revComp(seq2)
		lastShared = 0 - lastShared
		lastSharedReverse = 0 - lastSharedReverse # FIXME: i think the indices are getting fucked on this reversal step
	
	untenablePositions = []
	suitablePositions = []
	# ===============================
	# Ranges
	# (0) to (lastShared - motifLen)
	# (lastShared - motifLen + 1) to (lastShared + motifLen)
	# (lastShared + motifLen + 1) to (end)
	
	# ===============================
	# Check the left shared region for the presence of the motif
	
	# Get a prospective primer length
	tempPrimer = putativePrimer(seq1,lastShared)
	primerCutoff = lastShared - len(tempPrimer[0])
	
	# Worried about edge cases here
	if primerCutoff < motifLen:
		primerCutoff = 0
	
	for eachPosition in range(primerCutoff,lastShared - motifLen):
		currentRegion1 = seq1[eachPosition : motifLen + eachPosition]
		currentRegion2 = seq2[eachPosition : motifLen + eachPosition]
		
		region1Ham = hamming(currentRegion1,currentMotif) # FIXME: Check this!!! Make sure i'm actually being returned an int and not a list
		region2Ham = hamming(currentRegion2,currentMotif)
		
		if region1Ham == 0 and region2Ham == 0:
			untenablePositions.append(eachPosition)
		
		if (not region1Ham or not region2Ham) and (region1Ham != region2Ham):
			suitablePositions.append(eachPosition)
		
		
	# ===============================
	# Check to see if the current motif is present on the right side of the sequence (defined to be downstream here)
	rightSeq1 = revComp(seq1[lastSharedReverse:])
	rightSeq2 = revComp(seq1[lastSharedReverse:])
	rightLastShared = abs(lastSharedReverse)
	rightSharedMotif = False
	# FIXME: if it's the right shared sequence then I don't need to subset twice at all!
	
	# Figure out limit of right side
	# TODO: check to make sure this isn't fucked
	ampliconEnd = lastShared + Settings.ampliconLength + motifLen
	if ampliconEnd > min(len(seq1),len(seq2)):
		rightSideCutoff = 0
	else:
		rightSideCutoff = min(len(seq1),len(seq2)) - ampliconEnd # TODO: make sure I don't have a fencepost error here
	
	for eachPosition in range(rightSideCutoff,rightLastShared-motifLen):
		currentRightRegion1 = rightSeq1[eachPosition : motifLen + eachPosition]
		currentRightRegion2 = rightSeq2[eachPosition : motifLen + eachPosition]

		
		rightRegion1Ham = hamming(currentRightRegion1,currentMotif)
		rightRegion2Ham = hamming(currentRightRegion2,currentMotif)
		
		if 0 in [rightRegion1Ham,rightRegion2Ham]:
			rightSharedMotif = True
	
	if rightSharedMotif == True:
		# This motif is unsuable from this direction
		return([[],[]])
	
	# ===============================
	# Check to see if the motif is found in the unshared region
	# Seq1 unshared
	seq1Results = scanUnshared(seq1,currentMotif,lastShared,lastSharedReverse)
	seq2Results = scanUnshared(seq2,currentMotif,lastShared,lastSharedReverse)
	
	seq1USexact = seq1Results[1]
	seq1USpos = seq1Results[0]
	seq2USexact = seq2Results[1]
	seq2USpos = seq2Results[0]
	
	# Make sure there's not exact matches in both
	if not (len(seq1USexact) > 0 and len(seq2USexact) > 0):
		if seq1USpos is not []:
			for eachPosition in seq1USpos:
				suitablePositions.append(eachPosition)
		if seq2USpos is not []:
			for eachPosition in seq2USpos:
				suitablePositions.append(eachPosition)
	if direction is 'right':
		newUntenable = [0 - x for x in untenablePositions]
		newSuitable = [0 - x for x in suitablePositions]
		return([newUntenable,newSuitable])
	else:
		return([untenablePositions,suitablePositions])

def scanSingle(seq,cutSite,currentMotif,direction,lastShared,lastSharedReverse,exactMatch=False):
	"""
	Modified version of scanSequence() that examines only one sequence.
	
	seq = string, valid DNA bases
	cutSite = integer, index of cut site
	currentMotif = string, restriction enzyme motif
	direction = string, choice of "left" or "right"
	hammingThreshold = integer, number of tolerable mismatches
	ampliconLength = integer, desired length of dCAPS amplicon
	TM = desired primer TM
	primerType = string, choice of "tm" or "length"
	minLength = integer, length of primer
	"""
	# consider the sequence "abcdefghijklmnop" where a cut will occur between "f" and "g"
	# the last shared base will then be "f", at index 5 (0-based)
	# the reverse last shared base is "g", at index 9 (for reversed string) and index -10 (  0-(seqLen-lastShared)+1, for forward string)
	
	motifLen = len(currentMotif)
	
	possibleDirections = ["left","right"]
	directionComparison = [x in direction for x in possibleDirections]
	otherDirection = [d for d,s in zip(possibleDirections,directionComparison) if not s][0]
		
	
	if direction is 'right':
		seq = revComp(seq)
		lastShared = 0 - lastShared
		lastSharedReverse = 0 - lastSharedReverse
	
	untenablePositions = []
	suitablePositions = []
	# ===============================
	# Ranges
	# (0) to (lastShared - motifLen)
	# (lastShared - motifLen + 1) to (lastShared + motifLen)
	# (lastShared + motifLen + 1) to (end)
	
	# ===============================
	# Check the left shared region for the presence of the motif

	tempPrimer = putativePrimer(seq,lastShared)
	primerCutoff = lastShared - len(tempPrimer[0])

	# Taking care of an edge case
	if primerCutoff < motifLen:
		primerCutoff = 0
	# TODO: Check all this, I've only done a first pass on this to remove seq1/seq2 refs
	for eachPosition in range(primerCutoff,lastShared - motifLen):
		currentRegion = seq[eachPosition : motifLen + eachPosition]
		
		regionHam = hamming(currentRegion,currentMotif)
		
		if regionHam == 0:
			untenablePositions.append(eachPosition)
		
	# ===============================
	# Check to see if the current motif is present on the right side of the sequence (defined to be downstream here)
	rightSeq = revComp(seq[lastSharedReverse:]) 
	rightLastShared = 0 - lastSharedReverse
	rightSharedMotif = False
	# Figure out limit of right side
	ampliconEnd = lastShared + Settings.ampliconLength + motifLen
	
	if ampliconEnd > len(seq):
		rightSideCutoff = 0
	else: 
		rightSideCutoff = len(seq) - ampliconEnd # TODO: make sure i dont have a fencepost error here
	
	for eachPosition in range(rightSideCutoff,rightLastShared-motifLen):
		currentRightRegion = rightSeq[eachPosition : motifLen + eachPosition]
		rightRegionHam = hamming(currentRightRegion,currentMotif)
		
		if 0 in [rightRegionHam]:
			rightSharedMotif = True
	
	if rightSharedMotif == True:
		# This motif is unsuable from this direction
		return([[],[]])

	# ===============================
	# Check to see if the motif is found in the unshared region
	# Seq1 unshared
	seqResults = scanUnshared(seq,currentMotif,lastShared,lastSharedReverse)
	
	seqUSexact = seqResults[1]
	seqUSpos = seqResults[0]
	
	# Make sure there's not exact matches in both
	if not (len(seqUSexact) > 0):
		if seqUSpos is not []:
			for eachPosition in seqUSpos:
				suitablePositions.append(eachPosition)
	if direction is 'right':
		newUntenable = [0 - x for x in untenablePositions]
		newSuitable = [0 - x for x in suitablePositions]
		return([newUntenable,newSuitable])
	else:
		return([untenablePositions,suitablePositions])
		
def evaluateMutations(seq,targetSeq,enzymeInfo,enzymeName):
	"""
	Top-level function that searches a sequence to be 
	edited with CRISPR/Cas9 methods for sites likely to be 
	usable as screening sites.
	
	Returns list of strings of output to return to html renderer.
	"""
	# ALGO: check enzyme for case where it cuts in WT. if it doesn't, skip it. if it does, check the cut position. if the cut position is flanked by Ns, skip checking. if the cut position isn't, make probable edits.

	# This calls scanSequence(), which returns a two-element list of form:
	# [untenable positions, suitable positions]
	
	# Scan targetSeq across seq to find matching position
	# I assume there is only one cut site because a check should have performed for multiple match sites before it ever got this far. Also assuming there is exactly one because of same pre-checks.
	# Initializing variables
	currentMotif = enzymeInfo[0]
	canCutLeft = False
	canCutRight = False
	directions = ["left","right"]
	currentOutput = []
	output = []
	motifLen = len(currentMotif)
	
	# Check to see where the cut site is in the given sequence
	for eachPosition in range(0,len(seq)-(len(targetSeq)-1)):
		currentSubset = seq[eachPosition:(eachPosition+len(targetSeq))]
		comparisons = hamming(currentSubset,targetSeq,True,True)
		if 0 in comparisons:
			targetStart = eachPosition
			if comparisons[1] == 0 or comparisons[2] == 0:
				targetSeq = revComp(targetSeq)
			break
	
	# Identify cut site
	# Hard assumption that the cut site is at the -3 position from the 3' end of the provided 5'->3' target sequence
	if comparisons[0] == 0:
		cutPosition = eachPosition + len(targetSeq) - 3
	elif comparisons[1] == 0 or comparisons[2] == 0:
		cutPosition = eachPosition + 3
	
	# Find last shared base on each side of putative editing
	lastSharedLeft = len(seq)
	lastSharedRight = lastSharedLeft
	tempEditedSeqs = crisprEdit(seq,cutPosition)
	if len(tempEditedSeqs) > 1:
		for eachPair in itertools.permutations(tempEditedSeqs,2):
			tempLastSharedLeft = lastSharedBase(eachPair[0],eachPair[1],'left')
			tempLastSharedRight = lastSharedBase(eachPair[0],eachPair[1],'right')
			if tempLastSharedLeft < lastSharedLeft:
				lastSharedLeft = tempLastSharedLeft
			if tempLastSharedRight > lastSharedRight:
				lastSharedRight = tempLastSharedRight
	else:
		lastSharedLeft = lastSharedBase(tempEditedSeqs[0],seq,'left') # should be positive
		lastSharedRight = lastSharedBase(tempEditedSeqs[0],seq,'right') # should be negative
	
	
	# Check if this enzyme cuts in the WT
	sitesLeft = scanSingle(seq,cutPosition,currentMotif,'left',lastSharedLeft,lastSharedRight,exactMatch=False) # should be positive
	sitesRight = scanSingle(seq,cutPosition,currentMotif,'right',lastSharedRight,lastSharedLeft,exactMatch=False) # should be negative
	
	# Count how many possible sites have rejected primers
	rejectedPrimers = 0
	
	for eachNum in [0,1]:
		# Set up some initial vairables
		eachSet = [sitesLeft,sitesRight][eachNum]
		currentDirection = ['left','right'][eachNum]
		currentSet = [[],[]]
		lastShared = [lastSharedLeft,lastSharedRight][eachNum]
		lastSharedReverse = [lastSharedRight,lastSharedLeft][eachNum]

		# If there aren't any good sites to begin with, skip this loop
		if eachSet == [[],[]]:
			continue
		elif eachSet[1] != []: # if it doesn't cut in the wild-type, it will be []. if it cuts, it won't be [].
			# Start typesetting some output
			currentOut = []
			currentOut.append("===============================")
			
			# If any values are negative, reverse the sequences and indices
			if any([value < 0 for value in eachSet[1]]): # negative values exist
				# "1234"[-1] gives "4", while "4321"[1] gives "3"
				currentSet[1] = [0-x for x in eachSet[1]]
				currentSet[0] = [0-x for x in eachSet[0]]
				currentSeq = revComp(seq)
				currentTarget = revComp(targetSeq)
				currentMotif = currentMotif
				currentOut.append("Sequences are reversed.")
				lastSharedLeft = 0 - lastShared # FIXME: Check for off-by-one errors
				lastSharedRight = 0 - lastSharedReverse # FIXME: Check for off-by-one errors
			else:
				currentSet[0] = eachSet[0]
				currentSet[1] = eachSet[1]
				currentSeq = seq
				currentTarget = targetSeq
			# Provide information on cut sites
			currentOut.append("Enzyme: " + enzymeName)
			outputstring1 = "Possible cut site found at position " + str(currentSet[1]) + "."
			outputstring2 = "Problem cut site found at position " + str(currentSet[0]) + "."
			currentOut.append(outputstring1)
			currentOut.append(outputstring2)
			
			# Provide the WT sequence and the target site
			# TODO: indicate the cut site in the output somehow
			currentOut.append("WT Sequence: " + currentSeq)
			if currentDirection == "left":
				currentOut.append("Target Sequence:" + " "*(targetStart-3) + currentTarget) 
			elif currentDirection == "right":
				currentOut.append("Target Sequence:" + " "*(len(currentSeq)-3-targetStart-len(currentTarget)) + currentTarget) # TODO: check for off-by-one errors
			
			# Get proper enzyme direction and typeset the enzyme
			# This accounts for cases of non-palindromic enzymes
			for eachIndex in currentSet[1]:
				motifDiff = hamming(currentMotif,currentSeq[eachIndex:(eachIndex+motifLen)],allResults=True)
				if motifDiff[0] > min(motifDiff[1:2]):
					currentMotif = revComp(currentMotif)
				currentOut.append(" "*(13+eachIndex) + currentMotif)
			
			# Indicate any exact cut sites in the shared regions
			for eachIndex in currentSet[0]:
				currentOut.append(" "*(13+eachIndex) + "."*len(currentMotif))
			
			# Get the sites we want to make primers for
			desiredSuitable = currentSet[1]
			untenablePositions = currentSet[0]
			# Multiple putative sites may exist, so check at each one
			for eachSite in desiredSuitable:
				# set up a counter for how many editing events cut
				usableSite = 0
				
				# Edit sequence so that the primer works
				currentSite = currentSeq[eachSite:eachSite+motifLen]
				
				if revComp(currentMotif).lower() == currentMotif.lower():
					orientedMotif = currentMotif
				else:
					# TODO: check right orientation
					orientedMotif = currentMotif
				
				orientedMotif = list(orientedMotif)
				currentSite = list(currentSite)
				
				# Change bases and reconstruct altered sequence
				for each in range(0,lastSharedLeft-eachSite-1): 
					if orientedMotif[each].lower() in ['g','c','t','a'] and orientedMotif[each].lower() != currentSite[each].lower():
						currentSite[each] = orientedMotif[each] 
					elif orientedMotif[each].lower() in ['y','r','w','s','k','m','d','v','h','b'] and hamming(orientedMotif[each].lower(),currentSite[each].lower()) != [0]:
						currentSite[each] = getDegenerateMatch(orientedMotif[each]) # TODO/FIXME: right now it includes the degenerate base in the primer, but should I have it randomly pick a compatible base?
				orientedMotif = ''.join(orientedMotif)
				currentSite = ''.join(currentSite)
				seqLeft = currentSeq[:eachSite]
				seqRight = currentSeq[(eachSite+motifLen):]
				alteredSeq = seqLeft+currentSite+seqRight
								
				# Simulate CRISPR edits
				altSeqs = crisprEdit(alteredSeq,cutPosition)
				seqNum = len(altSeqs)
				
				# Set up variables to count cuts and tests
				comparisonCount = 0
				cutCount = 0
				
				# See if the enzyme cuts the edited sites
				for eachSequence in altSeqs:
					# Find the proportional distance between the edited and motif sequences
					proportionTest = proportionalDistance(eachSequence[desiredSuitable[0]:desiredSuitable[0]+motifLen],currentMotif)
					
					# Get the set with the best match
					# setToUse = [comparisons,currentProp,currentHam,hamDistHigh]
					# There are three lists in this list, representing the three possible comparison orientations
					setToUse = [comp for comp in proportionTest if comp[0] == min([comp2[0] for comp2 in proportionTest])][0]
					#print(setToUse)
					
					# Add the number of comparisons made
					comparisonCount += setToUse[0]
					
					# setToUse[1] is the proportion of comparisons that are exact matches
					# Add the number of possible cuts if any exist
					if setToUse[2] == 0:
						cutCount += setToUse[0] * setToUse[1]
					
				# setToUse[1] is currentProp, the proportion of all possible comparisons which are exact matches
				# currentProp * comparisons gives you the number of cuts you can expect from that comparison
				
				# Examine whether the proportion of viable cuts is above the threshold for this cut site

				# YOU WANT IT NOT TO CUT, YOU WANT cutCount TO BE AS LOW AS POSSIBLE
				# THE WHOLE POINT IS IT DOESN'T CUT THE MUTANT
				# I'M YELLING BECUASE THIS IS LIKE THE THIRD TIME I'VE LOOKED AT THIS AND THOUGHT 'hey why isn't it >=' AND SPENT LIKE TWO HOURS AUDITING MY CODE AND IM SICK OF IT HAPPENING
				
				if (100*cutCount/comparisonCount) <= Settings.seqThreshold:
					# Attempt to generate a primer
					newPrimer = generatePrimer(currentSeq,untenablePositions,eachSite,lastSharedLeft,currentMotif)
					
					# Typeset the primer if it worked
					if newPrimer is not None:
						currentOut.append("Possible screening primer found.")
						currentOut.append(" "*(13+lastSharedLeft-len(newPrimer[0]))+newPrimer[0])
						output.append(currentOut)
					else:
						rejectedPrimers += 1
	if rejectedPrimers == (len(sitesLeft[1])+len(sitesRight[1])):
		return(None)
	else:
		return(output)

def evaluateSites(seq1,seq2,enzymeInfo,enzymeName):
	"""
	Top-level function accepting sequences and enzyme information. 
	Calls other functions and typesets output for rendering by the server.

	Input
	seq1 = string, DNA bases
	seq2 = string, DNA bases
	enzymeInfo = list, [motif (string), cut position upper 
		(integer), cut position lower (integer)]
	hammingThreshold = integer
	enzymeName = string
	TM = number
	
	Output
	output = list of strings to be typeset by the web page, 
		listing results of screening for diagnostic primers
	"""	
	currentMotif = enzymeInfo[0]
	# Search from Left
	sitesLeft = scanSequence(seq1,seq2,currentMotif,'left')
	# Is [ [untenable positions] , [suitable positions] ]
	# Search from Right
	sitesRight = scanSequence(seq1,seq2,currentMotif,'right')
	# Is [ [untenable positions] , [suitable positions] ]
	
	output = []
	motifLen = len(currentMotif)
	rejectedPrimers = 0
	
	for eachSet in [sitesLeft,sitesRight]:
		currentRejects = 0
		currentSet=[[],[]]
		if eachSet == [[],[]]:
			return(None)
		elif eachSet[1] != []:
			currentOut = []
			currentOut.append("===============================")
			# REMEMBER: I can't BOTH flip the sequences and indices, that just won't work.
			if any([value < 0 for value in eachSet[1]]):
				# "1234"[-1] gives "4", while "4321"[1] gives "3"
				currentSet[1] = [abs(x) for x in eachSet[1]]
				currentSet[0] = [abs(x) for x in eachSet[0]]
				currentSeq1 = revComp(seq1)
				currentSeq2 = revComp(seq2)
				currentMotif = currentMotif
				currentOut.append("Sequences are reversed.")
			else:
				currentSet[0] = eachSet[0]
				currentSet[1] = eachSet[1]
				currentSeq1 = seq1
				currentSeq2 = seq2
				currentMotif = currentMotif
			currentOut.append("Enzyme: " + enzymeName)
			outputstring1 = "Possible cut site found at position " + str(currentSet[1]) + "."
			outputstring2 = "Problem cut site found at position " + str(currentSet[0]) + "."
			currentOut.append(outputstring1)
			currentOut.append(outputstring2)
			currentOut.append("Sequence 1: " + currentSeq1)
			currentOut.append("Sequence 2: " + currentSeq2)
			
			# FIXME: I need to merge this iteration over currentSet[1] with the one in the next block because right now it might print a recognizable site that I later ignore because it's not a good primer
			for eachIndex in currentSet[1]:
				# Get direction to typeset the enzyme
				motifDiff1 = hamming(currentMotif,currentSeq1[eachIndex:(eachIndex+motifLen)], allResults=True)
				motifDiff2 = hamming(currentMotif,currentSeq2[eachIndex:(eachIndex+motifLen)], allResults=True)
				motifComp1 = motifDiff1[0] > min(motifDiff1[1:2])
				motifComp2 = motifDiff2[0] > min(motifDiff2[1:2])
				if motifComp1 or motifComp2:
					currentOut.append("reversing motif")
					currentMotif = revComp(currentMotif)
				currentOut.append(" "*(12+eachIndex) + currentMotif)
			for eachIndex in currentSet[0]:
				currentOut.append(" "*(12+eachIndex) + "."*len(currentMotif))
			
			# FIXME: Does this need to be in the scope where I define the motifDiff objects?
			if 0 in motifDiff1 or 0 in motifDiff2:
				currentOut.append("CAPS primer possible")
			
			currentLastShared = lastSharedBase(currentSeq1,currentSeq2,'left')
			for eachIndex in currentSet[1]:
				newPrimer = generatePrimer(currentSeq1,currentSet[0],eachIndex,currentLastShared,currentMotif)

				if newPrimer is not None:
					lastPrimerBase = newPrimer[0][-1] # FIXME: I don'tt think this and hte next line are ever used???
					lastSequenceBase = currentSeq1[currentLastShared-1]
					currentOut.append(" "*(12+currentLastShared-len(newPrimer[0]))+newPrimer[0])
					
					# Find out which of seq1, seq2 will be cut
					# Get amplified sequence
					amplifiedSeq1 = newPrimer[0]+currentSeq1[(currentLastShared):(currentLastShared+len(currentMotif))]
					amplifiedSeq2 = newPrimer[0]+currentSeq2[(currentLastShared):(currentLastShared+len(currentMotif))]
					amp1Cut = False
					amp2Cut = False
					for each in range(0,len(amplifiedSeq1)-motifLen):
						ampSub = amplifiedSeq1[each:(each+motifLen)]
						ampHam = hamming(currentMotif,ampSub,allResults=True)
						if any([x == 0 for x in ampHam]):
							amp1Cut = True
					if amp1Cut == False:
						for each in range(0,len(amplifiedSeq2)-motifLen):
							ampSub = amplifiedSeq2[each:(each+motifLen)]
							ampHam = hamming(currentMotif,ampSub,allResults=True)
							if any([x == 0 for x in ampHam]):
								amp2Cut = True
					cutStrand = "neither sequence"
					if amp1Cut == True:
						cutStrand = "Sequence 1"
					elif amp2Cut == True:
						cutStrand = "Sequence 2"
					
					currentOut.append("Assay for above primer cuts "+cutStrand+". The primer has Tm "+str(estimateTM(newPrimer[0]))+" degrees C.")
				else:
					currentRejects += 1
					continue
			
			if currentRejects == len(currentSet[1]):
				continue
			output.append(currentOut)
			rejectedPrimers += currentRejects
	if rejectedPrimers == len(sitesLeft[1]) + len(sitesRight[1]):
		return(None)
	return(output) 

def putativePrimer(seq,lastShared):
	"""
	Generate a mock primer based on desired TM or length 
	and end position. This is used to estimate whether an 
	exact match restriction site found in the shared region 
	of two sequences is likely to be captured by a primer 
	(rendering it necessary to modify the site or throw out 
	the enzyme) or if it can be safely ignored.
	
	Input
	seq = string of valid DNA bases
	lastShared = integer indicating last base of primer
	
	Output
	bestPrimer = string of valid DNA bases
	"""
	# type can be 'tm' or 'length'
	seq = seq.lower()
	
	# Generate primer sub-sequences
	currentStart = 0
	primerList = []
	
	while currentStart >= 0:
		currentPrimer = seq[currentStart:lastShared]
		primerList.append(currentPrimer)
		currentStart -= 1
	
	if Settings.primerType == 'tm':
		output = []
		for eachPrimer in primerList:
			output.append(estimateTM(eachPrimer))
		
		# Filter for negative values
		filterList = [x for x in range(0,len(output)) if output[x] <= 0]
		primerList = [primerList[x] for x in range(0,len(output)) if x not in filterList]
		output = [output[x] for x in range(0,len(output)) if x not in filterList]
		# Find minimum diff
		outputDiff = [abs(x - Settings.TM) for x in output]
		# Choose the best primer sub-sequence based on difference from optimum Tm
		bestPrimer = [primerList[x] for x in range(0,len(outputDiff)) if outputDiff[x] == min(outputDiff)]
		return(bestPrimer)
	elif Settings.primerType == 'length':
		# compare length of primers in list to optimum length
		positionList = list(range(0,len(primerList)))
		filterList = [abs(len(x) - Settings.minLength) for x in primerList]
		bestPrimer = [primerList[x] for x in positionList if filterList[x] == min(filterList)]
		return(bestPrimer)

def evaluateIsogenic(wtSeq,mutSeq,targetSeq,enzymeInfo,enzymeName):
	"""
	Top-level function generating primer for identifying an
	isogenic mutation in a CRISPR mutagenesis experiment.
	
	This tool will find a primer that cuts only the desired
	mutant sequence. The wild-type and undesired mutant
	sequences will not be cut by the restriction digest.
	
	Input
	wtSeq = 
	mutSeq = 
	targetSeq = 
	enzymeInfo = 
	enzymeName = 
	
	Output
	"""
	
	# Set up some variables
	currentMotif = enzymeInfo[0]
	directions = ["left","right"]
	currentOutput = []
	output = []
	motifLen = len(currentMotif)
	targetStart = None
	
	# Check wtSeq to see where the cut site is in the given sequence
	for eachPosition in range(0,len(wtSeq)-(len(targetSeq)-1)):
		currentSubset = wtSeq[eachPosition:(eachPosition+len(targetSeq))]
		comparisons = hamming(currentSubset,targetSeq,True,True)
		if 0 in comparisons:
			targetStart = eachPosition
			if comparisons[1] == 0 or comparisons[2] == 0:
				targetSeq = revComp(targetSeq)
			break
	
	# If there wasn't any cut site, just kill the function
	if targetStart == None:
		return(None)
	
	# Identify cut site
	# Hard assumption that the cut site is at the -3 position from the 3' end of the provided 5'->3' target sequence
	if comparisons[0] == 0:
		cutPosition = eachPosition + len(targetSeq) - 3
	elif comparisons[1] == 0 or comparisons[2] == 0:
		cutPosition = eachPosition + 3

	# Find last shared base on each side of putative editing
	lastSharedLeft = len(wtSeq)
	lastSharedRight = lastSharedLeft
	tempEditedSeqs = crisprEdit(wtSeq,cutPosition)
	
	# Make sure mutantSeq is in the list of edited sequences
	if mutantSeq not in tempEditedSeqs:
		tempEditedSeqs += mutantSeq
	
	# Find best shared base for set of edited sequences
	if len(tempEditedSeqs) > 1:
		for eachPair in itertools.permutations(tempEditedSeqs,2):
			tempLastSharedLeft = lastSharedBase(eachPair[0],eachPair[1],'left')
			tempLastSharedRight = lastSharedBase(eachPair[0],eachPair[1],'right')
			if tempLastSharedLeft < lastSharedLeft:
				lastSharedLeft = tempLastSharedLeft
			if tempLastSharedRight > lastSharedRight:
				lastSharedRight = tempLastSharedRight
	else:
		lastSharedLeft = lastSharedBase(tempEditedSeqs[0],wtSeq,'left') # should be positive
		lastSharedRight = lastSharedBase(tempEditedSeqs[0],wtSeq,'right') # should be negative
	
	# Check sites left and right to see if this enzyme cuts in the mutant
	sitesLeft = scanSingle() # This needs to be defined for the mutant! because all the indices are based on checking for cutting in the mutant
	sitesRight = scanSingle() # TODO: what does exactMatch need to be here...
	
	# Count rejected primers
	rejectedPrimers = 0
	
	# Iterate from left and right
	for eachNum in [0,1]:
		# Set up some initial variables
		eachSet = [sitesLeft,sitesRight][eachNum]
		currentDirection = ['left','right'][eachNum]
		currentSet = [[],[]]
		lastShared = [lastSharedLeft,lastSharedRight][eachNum]
		lastSharedReverse = [lastSharedRight,lastSharedLeft][eachNum]
		
		# If there aren't good sites, skip this loop
		if eachSet == [[],[]]:
			continue
		elif eachSet[1] != []:
			# Start typesetting output
			currentOut = []
			currentOut.append("===============================")
		
			# If any values are negative, reverse sequences and indices
			if any([value <0 for value in eachSet[1]]): # negative values exist
				currentSet[1] = [0-x for x in eachSet[1]]
				currentSet[0] = [0-x for x in eachSet[0]]
				currentWTSeq = revComp(wtSeq)
				currentMutSeq = revComp(mutSeq)
				currentTarget = revComp(targetSeq)
				currentMotif = currentMotif # Not reversing this because later on I just reverse it anyway for checks when necessary
				currentOut.append("Sequences are reversed.")
				lastSharedLeft = 0 - lastShared # FIXME: Check for off-by-one errors
				lastSharedRight = 0 - lastSharedReverse # FIXME: Check for off-by-one errors
			else:
				currentSet[0] = eachSet[0]
				currentSet[1] = eachSet[1]
				currentWTSeq = wtSeq
				currentMutSeq = mutSeq
				currentTarget = targetSeq
				currentMotif = currentMotif
			
			# Typeset information on cut sites
			# TODO: I don't need this information for this tool
			currentOut.append("Enzyme: " + enzymeName)
			outputstring1 = "Possible cut site found at position " + str(currentSet[1]) + "."
			outputstring2 = "Problem cut site found at position " + str(currentSet[0]) + "."
			currentOut.append(outputstring1)
			currentOut.append(outputstring2)
			
			# Typeset mutant sequence and the target site
			currentOut.append("WT Sequence: " + currentWTSeq)
			seqlens = [len(x) for x in [currentWTSeq,currentMutSeq]]
			
			# Generate offsets for wt/mut typesetting
			if len(currentWTSeq) > len(currentMutSeq):
				wtOffset = 0
				mutOffset = max(seqlens) - min(seqlens) # TODO: check for off-by-one errors
			else:
				wtOffset = max(seqlens) - min(seqlens) # TODO: check for off-by-one errors
				mutOffset = 0
			
			if currentDirection == "left":
				currentOut.append("Target Sequence:" + " "*(targetStart-3) + currentTarget)
			elif currentDirection == "right":
				currentOut.append("Target Sequence:" + " "*(len(currentWTSeq-3)-targetStart-len(currentTarget) + wtOffset)) # TODO: check for off-by-one errors
			
			# Get proper enzyme direction and typeset the enzyme
			# This accounts for cases of non-palindromic enzymes
			for eachIndex in currentSet[1]:
				motifDiff = hamming(currentMotif,currentMutSeq[eachIndex:(eachIndex+motifLen)],allResults=True)
				if motifDiff[0] > min(motifDiff[1:2]):
					currentMotif = revComp(currentMotif)
				currentOut.append(" "*(13+eachIndex+mutOffset)+currentMotif)
			
			# Indicate any exact cut sites in the shared region
			for eachIndex in currentSet[0]:
				currentOut.append(" ")*(13+eachIndex+mutOffset) + "."*len(currentMotif)
				
			desiredSuitable = currentSet[1]
			untenablePositions = currentSet[0]
			# Multiple putative sites may exists, so check at each one
			for eachSite in desiredSuitable:
				# Set up a counter for how many editing events cut
				usableSite = 0
				
				# Edit sequence so that the primer works
				currentSite = currentMutSeq[eachSite:eachSite+motifLen]
				
				if revComp(currentMotif).lower() == currentMotif.lower():
					orientedMotif = currentMotif
				else: # TODO: check right orientation. wait do I even need to do this? I thought I checked and covered this earlier.
					orientedMotif = currentMotif
				
				orientedMotif = list(orientedMotif)
				currentSite = list(currentSite)
				
				# Change bases and reconstruct altered sequence
				for each in range(0,lastSharedLeft - eachSite - 1):
					if orientedMotif[each].lower() in ['g','c','a','t'] and orientedMotif[each].lower() != currentSite[each].lower():
						currentSite[each] = orientedMotif[each]
					elif orientedMotif[each].lower() in ['y','r','w','s','k','m','d','v','h','b'] and hamming(orientedMotif[each].lower(),currentSite[each].lower()) != [0]:
						currentSite[each] = getDegenerateMatch(orientedMotif[each]) # TODO/FIXME: right now it includes the degenerate base in the primer, but should I have it randomly pick a compatible base?
				orientedMotif = ''.join(orientedMotif)
				currentSite = ''.join(currentSite)
				seqLeft = currentMutSeq[:eachSite]
				seqRight = currentMutSeq[(eachSite+motifLen):]
				alteredSeq = seqLeft+currentSite+seqRight
				
				# Simulate CRISPR edits
				# TODO: If mutSeq is in this set, remove it
				altSeqs = crisprEdit(alteredSeq,cutPosition)
				if currentMutSeq in altSeqs:
					altSeqs = [x for x in altSeqs if x != currentMutSeq]
				seqNum = len(altSeqs)
				
				# TODO: I need to think carefully about how indels might overlap with the edited sequences.
				
				# Set up variables to count cuts and tests
				comparisonCount = 0
				cutCount = 0
				
				# See if the enzyme cuts the edited sites
				for eachSequence in altSeqs:
					# Find proportional distance between edited and motif sequences
					proportionTest = proportionalDistnace(eachSequence[desiredSuitable[0]:desiredSuitable[0]+motifLen],currentMotif)
					
					# Get the set with the best match
					setToUse = [comp for comp in proportionTest if comp[0] == min([comp2[0] for comp2 in proportionTest])][0]
					
					# Add the number of comparisons made
					comparisonCount += setToUse[0]
					
					# Add the number of possible cuts if any exist
					if setToUse[2] == 0:
						cutCount += setToUse[0] * setToUse[1]
				
				# Examine whether the proportion of cuts is above the threshold for the site
				# TODO/FIXME: make sure the logic of this is correct! I want it to cut only mutSeq and NOT wtSeq or any other edited sequence!
				if (100*cutCount/comparisonCount) <= Settings.seqThreshold:
					# Attempt to generate a primer
					newPrimer = generatePrimer(currentSeq,untenablePositions,eachSite,lastSharedLeft,currentMotif)
					
					# Typeset the primer if it worked
					if newPrimer is not None:
						currentOut.append("Possible screening primer found.")
						currentOut.append(" "*(13+lastSharedLeft-len(newPrimer[0]))+newPrimer[0]+mutOffset)
						output.append(currentOut)
					else:
						rejectedPrimers += 1
	
	if rejectedPrimers == (len(sitesLeft[1])+len(sitesRight[1])):
		return(None)
	else:
		return(output)
		
def generatePrimer(seq,untenablePositions,desiredSuitable,lastShared,currentMotif):
	"""
	Generates primer for a given sequence and checks 
	it for errors.
	
	Input
	seq = string of valid DNA bases
	untenablePositions = list of integers, indicating 
	exact matches of the restriction enzyme motif which 
	must be killed with mismatches
	desiredSuitable = integer, starting position of the 
	restriction enzyme recognition site
	lastShared = integer, position of the last shared base 
	between two sequences. Primer will start here and move 
	backwards.
	currentMotif = string of the restriction enzyme 
		recognition motif to search for
	
	Output
	bestPrimer = None or string of valid DNA bases
	"""
	# Start the primer at the last shared base
	# If the last shared base is a mismatch to the motif, check the rules for that
	# Start extending a primer backwards
	# for each position, check melting temp
	# continue if it's in good melting temp range
	# if you have a complete upstream motif, then you have to change a base
	# Take untenablePositions so we know all the positions we need to edit
	# But this function will only return a primer for a single suitable position, 
	# So only supply one position
	
	# FIXME: Right now, it changes bases to modify the sequence before selecting a primer based on Tm
	# This can run into cases where it modifies a base that it doesn't need to.
	# Example: in one case, it modified a GGGCCC site to be GAGCCC, but the primer the function reported
	# Didn't even have a valid GGGCCC site, because the first base in the primer was partway through that sequence
	
	# TODO: Check global settings object for allowGTMM
	originalSeq = seq.lower()
	currentMotif = currentMotif.lower()
	motifLen = len(currentMotif)
	currentStart = lastShared
	primerList = []
	output = []

	# Figure out which base of the motif will have to be changed
	motifBase = 0
	while currentMotif[motifBase].lower() == "n":
		motifBase += 1
	
	# Make any necessary edits to the sequence
	# In general, edit the most-upstream base in a motif
	# But also check for overlaps with other untenable positions
	if untenablePositions != []:
		if len(untenablePositions) == 1:
			# Case where only one untenable position
			seq = list(seq)
			positionToModify = untenablePositions[0] + motifBase
			oldBase = seq[positionToModify].lower()
			seq[positionToModify] = revComp(oldBase)
			seq = ''.join(seq)
		else:
			# Case where more than one untenable position
			indexIterator = iter(range(0,len(untenablePositions)-1))
			for eachIndex in indexIterator:
				eachPosEnd = untenablePositions[eachIndex] + motifLen
				
				lookAhead = 0
				while eachPosEnd > untenablePositions[eachIndex+lookAhead]:
					if eachIndex+lookAhead >= len(untenablePositions) or (eachIndex+lookAhead)>= len(untenablePositions)-1:
						break
					lookAhead += 1
				
				# Modify the base
				seq = list(seq)
				positionToModify = untenablePositions[eachIndex+lookAhead] + motifBase
				oldBase = seq[positionToModify].lower()
				seq[positionToModify] = revComp(oldBase)
				seq = ''.join(seq)
				
				# Skip positions as necessary				
				for each in range(0,lookAhead):
					next(indexIterator,None)
				
	# Modify the motif if necessary
	currentSite = seq[desiredSuitable:desiredSuitable+motifLen]
	if revComp(currentMotif).lower() == currentMotif.lower():
		orientedMotif = currentMotif
	else:
		# TODO: Have it figure out the right orientation
		orientedMotif = currentMotif
	orientedMotif = list(orientedMotif)
	currentSite = list(currentSite)
	for each in range(0,lastShared-desiredSuitable): # FIXME: think i have a fencepost error here, needs to be lastshared-desiredSuitable+1 I think.
		if orientedMotif[each].lower() in ['g','c','t','a'] and orientedMotif[each].lower() != currentSite[each].lower():
			currentSite[each] = orientedMotif[each]
		elif orientedMotif[each].lower() in ['y','r','w','s','k','m','d','v','h','b'] and hamming(orientedMotif[each].lower(),currentSite[each].lower()) != [0]:
			currentSite[each] = getDegenerateMatch(orientedMotif[each]) # TODO/FIXME: right now it includes the degenerate base in the primer, but should I have it randomly pick a compatible base?
	orientedMotif = ''.join(orientedMotif)
	currentSite = ''.join(currentSite)
	seqLeft = seq[:desiredSuitable]
	seqRight = seq[(desiredSuitable+motifLen):]
	seq = seqLeft + currentSite + seqRight
	
	# Generate primer sub-sequences
	while currentStart >= 0:
		currentPrimer = seq[currentStart:lastShared]
		primerList.append(currentPrimer)
		currentStart -= 1
	for eachPrimer in primerList:
		output.append(estimateTM(eachPrimer))
	
	# Filter for negative values
	filterList = [x for x in range(0,len(output)) if output[x] <= 0]
	primerList = [primerList[x] for x in range(0,len(output)) if x not in filterList]
	output = [output[x] for x in range(0,len(output)) if x not in filterList]
	# TODO: it always checks primer length based on TM right now, have it switch based on user setting
	
	# Find minimum diff
	outputDiff = [abs(x - Settings.TM) for x in output]
	# Choose the best primer sub-sequence based on difference from optimum Tm
	bestPrimer = [primerList[x] for x in range(0,len(outputDiff)) if outputDiff[x] == min(outputDiff)]
	# See if the primer violates any rules
	# Rule: G/T mismatch works at 3' end, but G/A and G/G mismatches don't. (Simsek 2000)
	# Rule: Mismatches to canonical shared sequence less than hamming dist
	lastPrimerBase = bestPrimer[0][-1].lower()
	lastSeqBase = originalSeq[lastShared-1] # wait should this just be lastshared and not -1?
	
	# Check primer against overall hamming distance
	bestPrimerStart = lastShared - len(bestPrimer[0])
	hamTest = hamming(bestPrimer[0],originalSeq[bestPrimerStart:lastShared],False,True)
	if hamTest > Settings.hammingThreshold:
		return(None)
	
	# Check primer against mismatch criteria
	if lastPrimerBase != lastSeqBase and Settings.allowMismatch == True:
		# get the template base, which is the complement of the listed base
		templateBase = revComp(lastSeqBase).lower()
		templateBase = templateBase.lower()
		# compare to the 3' base
		# if template is G and primer 3' is T, ok.
		# Allowable: T/T, T/C, T/G
		if templateBase == 'g' and lastPrimerBase == 't':
			return(bestPrimer)
		elif templateBase == 't' and lastPrimerBase in ['t','c','g']:
			return(bestPrimer)
		else:
			return(None)
	elif lastPrimerBase != lastSeqBase and Settings.allowMismatch == False:
		return(None)
	else:
		return(bestPrimer)

def crisprEdit(seq,position):
	"""
	Makes edits to sequence and returns edited sequences.
	
	Input
	seq = string of valid DNA bases
	position = integer indicating cut site in provided 
		sequence. usually between the 17th and 18th bases 
		in the 20 bp gRNA target sequence.
	organism = string indicating organism being mutated	
	
	Output
	output = list of modified strings
	"""
	output = []
	seqNew = seq
		
	# Get the edit profile for the current organism
	editsToMake = Settings.organismEditProfile
	
	for eachEdit in editsToMake:
		# Check to see if it's an insertion or deletion
		if eachEdit > 0:
			# Add in X bases at the indicated position
			currentSeqLeft = seq[:position]
			currentSeqRight = seq[position:]
			seqNew = currentSeqLeft + 'X'*eachEdit + currentSeqRight
			output.append(seqNew)
		else:
			# Generate all subsequences
			subsequences = [[x, eachEdit-x] for x in range(eachEdit+1)]
			for eachSet in subsequences:
				currentSeqLeft = seq[:(eachEdit-eachSet[0])]
				currentSeqRight = seq[(eachEdit+eachSet[1]):]
				seqNew = currentSeqLeft + currentSeqRight
				output.append(seqNew)
	return(output)


# OLDER DESIGN NOTES - Might be some mistakes in here or abandoned ideas. Mostly keeping in case I need it later.

# Some old cruft from earlier versions but useful to have some sequences around for testing purposes.
# Sample sequences from Leah
#seq1 = "TGTGTGTGCAGGGAGAAGCCAAATGTGGATTTTGACAGGGTGGACTCGTATGTGCATCAG"
#seq2 = "TGTGTGTGCAGGGAGAAGCCAAATGTGGATTTGACAGGGTGGACTCGTATGTGCATCAG"
# other random sequences
#seq1 = 'GCGGAgggcccCTCAAGATCCCGAGTgggTCTTATcccCAGTTTCTTGGCTCTGTTA'
#seq2 = 'GCGGAgggcccCTCAAGATCCCGAGTgggcccCAGTTTCTTGGCTCTGTTA'
#currentMotif = 'gggccc'

# overall region to test:
# [lastShared - (motifLen-2) + iterator] to [lastShared + 1 + iterator]
# iterator is [0 .. motifLen-2]
# shared region to test is [lastShared - (motifLen-2) + iterator] to [lastShared]
# unshared region to test is [lastShared + 1] to [lastShared + 1 + iterator]
# motifLen is the number of characters in the motif
# lastShared is the index of the last shared character


# cases that work
# extant exact match overlapping with unshared region that isn't present in any upstream or downstream shared region. easiest csae
# extant exact match overlapping with unshared region that is present in an upstream or downstream shared region, shared match must be altered or avoided. medium difficulty case
# match created by altering one or more bases in the upstream shared half, with or without upstream/downstream matches. hardest case.
# match created where altered bases are in the unshared region. this is an edge case, and relies on, say, there being this pattern: sssUsU, where s is shared base and U is unshared base, and you need to change the first unshared base so that the second unshared base is discriminatory between the two. THIS WILL BE DIFFICULT TO FIND AND DESIGN

# Example:
# sssCCATT
# sssCTCAT
# where desired motif is GCAT
# primer would be: sssCG, so you have as amplicons:
# CG	ATT
# CG	CAT
# the first C is a shared base, so it's naturally part of the primer, but the final G in the primer is actually in the unshared region, and the resulting motif does not overlap at all with the shared region.