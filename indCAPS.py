#!/usr/bin/env python2
"""@package indCAPS
Examine multiple strings representing DNA bases for cases where 
diagnostic PCR primers can be constructed to distinguish the strings.
"""
# -*- coding: utf-8 -*-
#seq1 = 'GCGGAgggcccCTCAAGATCCCGAGTgggTCTTATcccCAGTTTCTTGGCTCTGTTA' #arguments[1]
#seq2 = 'GCGGAgggcccCTCAAGATCCCGAGTgggcccCAGTTTCTTGGCTCTGTTA' #arguments[2]
#currentMotif = 'gggccc'
#hamNum = 4 #arguments[3]
from math import log

# Sample sequences from Leah
seq1 = "TGTGTGTGCAGGGAGAAGCCAAATGTGGATTTTGACAGGGTGGACTCGTATGTGCATCAG"
seq2 = "TGTGTGTGCAGGGAGAAGCCAAATGTGGATTTGACAGGGTGGACTCGTATGTGCATCAG"

## Enzymes
# Will not include nicking enzymes or intron-encoding enzymes or homing enzymes
# Also does not include two-site cutters: BaeI, BsaXI, CspCI, BcgI
from enzymeList import enzymes


## Custom Classes
class SettingsObject:
	"""
	This object should be in the global scope and hold all the 
	invariant information on things like parameters of primer design.
	
	TM = numerical, desired melting temperature of primer
	ampliconLength = integer, length of dCAPS amplicon
	primerType = string, choice of 'tm' or ''
	primerLength = integer, length of primer to design
	allowMismatch = boolean, whether to allow 3' mismatch in primrs
	hammingThreshold = integer, mismatch threshold
	organism = string, species name
	organismEditProfile = unknown format, for specifying probable crispr edits
	"""
	def __init__(self,TM=55,ampliconLength=70,primerType='tm',primerLength=25,allowMismatch=False,hammingThreshold=1,organism=None,sodiumConc=0.05,primerConc=50*10**(-9)):
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
		
		commonEdits =  {"Athaliana":[],
					"Osativa":[],
					"Hsapiens":[],
					"Mmusculus":[]}
		if self.organism in commonEdits:
			self.organismEditProfile = commonEdits[organism]
		else:
			self.organismEditProfile = [1,-1,2,-2,3,-3]


## Custom Functions
def hammingBool(seq1,seq2,allResults=False):
	# Instead of returning a singular value, returns lists of bools indicating whether there's a match
	return(False)
	
def proportionalDistance(mut,motif):
	"""
	very incomplete function
	"""
	# Mut will have X characters.
	# Returns assessment of most optimistic match and also proportion of all matches that it will match
	possibleMatches = {'a':['a','r','w','m','d','h','v','n'],
				  'g':['g','r','s','k','b','d','v','n'],
				  'c':['c','y','s','m','b','h','v','n'],
				  't':['t','y','w','k','b','d','h','n'],
				  'r':['a','g','n'],
				  'y':['c','t','n'],
				  's':['g','c','n'],
				  'w':['a','t','n'],
				  'k':['g','t','n'],
				  'm':['a','c','n'],
				  'b':['c','g','t','n'],
				  'd':['a','g','t','n'],
				  'h':['a','c','t','n'],
				  'v':['a','c','g','n'],
				  'n':['a','g','c','t','r','y','s','w','k','m','b','d','h','v','n']}
	
	mut = mut.lower()
	motif = motif.lower()
	
	sets = [ [mut,motif] , [revComp(mut),motif] , [mut, revComp(motif)] ]
	
	setOutput = [0,0] * len(sets)
	
	for eachSet in sets:
		currentMut = eachSet[0]
		currentMotif = eachSet[1]
		currentHam = 0
		currentProp = 1
		for eachPos in range(0,len(currentMut)):
			mutBase = currentMut[eachPos]
			motifBase = currentMotif[eachPos]
			
			if mutBase == "x":
				if mutBase in ['g','c','a','t']:
					currentProp = currentProp * 1/4
				
				elif mutBase in ['r','y','s','w','k','m']:
					currentProp = currentProp * 1/2
				
				elif mutBase in ['b','d','h','v']:
					currentProp = currentProp * 3/4
				
				elif mutBase == "n":
					return(None)
				
			else:
				return(None)
	return(None)
	
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
	
	if len(seq1) != len(seq2):
		warnings.warn("Comparing sequences of unequal length.")
		print("Seq 1 is: ["+seq1+"]")
		print("Seq 2 is: ["+seq2+"]")
#	
#	distance1 = len(seq1) - sum([x == y for (x,y) in zip(seq1,seq2)])
#	# These next two are for cases of non-palindromic motifs
#	distance2 = len(seq1) - sum([x == y for (x,y) in zip(seq1,revComp(seq2))])
#	distance3 = len(seq1) - sum([x == y for (x,y) in zip(revComp(seq1),seq2)])
#	return(min(distance1,distance2,distance3))
	possibleMatches = {'a':['a','r','w','m','d','h','v','n'],
				  'g':['g','r','s','k','b','d','v','n'],
				  'c':['c','y','s','m','b','h','v','n'],
				  't':['t','y','w','k','b','d','h','n'],
				  'r':['a','g','n'],
				  'y':['c','t','n'],
				  's':['g','c','n'],
				  'w':['a','t','n'],
				  'k':['g','t','n'],
				  'm':['a','c','n'],
				  'b':['c','g','t','n'],
				  'd':['a','g','t','n'],
				  'h':['a','c','t','n'],
				  'v':['a','c','g','n'],
				  'n':['a','g','c','t','r','y','s','w','k','m','b','d','h','v','n']}
	# NOTE: degenerate bases should only occur in enzyme recognition motifs
	# This should allow motifs with degenerate bases to match given bases in a genomic sequence
	# However, if this were coded for general use, it would produce incorrect results
	# eg, 'r' is returned here as a valid match for 'n', when that is not strictly true
	
	sets = [ [seq1,seq2] , [revComp(seq1).lower(),seq2] , [seq1, revComp(seq2).lower()] ]
	
	setOutput = [0] * len(sets)
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
			 'b':'v'}
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
	# I really don't need this function, I dont know why I included it.
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
	# Units kcal/mole
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
	return(None)

def basicTemp(seq):
	return(None)

def nearestNeighbor(seq):
	# typical salt conc is [50 mM], keep it between [0.01 M] and [1 M] 
	# typical primer conc is 50 nM
	# keep it above 8 bases
	# From http://biotools.nubic.northwestern.edu/OligoCalc.html and cited sources
	R = 1.9872036 # kcal / (mole*kelvin)
	T = (1000*(-deltaH(seq))) / ((-deltaS(seq))+R*log(1/Settings.primerConc))-272.9+16.6*log(Settings.sodiumConc,10)
	return(T)


def estimateTM(seq,func='nearestNeighbor'):
	# Default function is 'nearestNeighbor'
	# see http://biotools.nubic.northwestern.edu/OligoCalc.html
	# Find number of each base
	
	if func is 'nearestNeighbor':
		return(nearestNeighbor(seq))
	
	# primerConc is 50 nM, pH 7.0
	
	# If less than 8 bases in oligo:
	
	# If 13 or fewer bases in oligo:
	
	# If 14 or more bases in oligo:
	return(None)


def lastSharedBase(seq1,seq2,direction):
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
	
	# Reverse strings if necessary
	if direction.lower() == "right":
		seq1 = revComp(seq1)
		seq2 = revComp(seq2)
	
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
		return([[0],[0]])
	
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

def scanSingle(seq,cutSite,currentMotif,direction):
	"""
	Modified version of scanSequence() that examines only one sequence.
	
	seq = string, valid DNA bases
	cutSite = integer, index of cut site
	currentMotif = string, restriction enzyme motif
	direction = string, choice of "left" or "right"
	hammingThreshold = integer, number of tolerable mismatches
	ampliconLength = integer, desired length of dCAPS amplicon
	TM = desired primer TM
	primerType = string, choice of "" or ""
	minLength = integer, length of primer
	"""
	# consider the sequence "abcdefghijklmnop" where a cut will occur between "f" and "g"
	# the last shared base will then be "f", at index 5 (0-based)
	# the reverse last shared base is "g", at index 9 (for reversed string) and index -10 (  0-(seqLen-lastShared)+1, for forward string)
	motifLen = len(currentMotif)
	lastShared = cutSite 
	lastSharedReverse = 0 - (len(seq) - cutSite) + 1
	
	possibleDirections = ["left","right"]
	directionComparison = [x in direction for x in possibleDirections]
	otherDirection = [d for d,s in zip(possibleDirections,directionComparison) if not s][0]
	
	if direction is 'right':
		seq = revComp(seq)
		lastShared = 0 - lastShared
		lastSharedReverse = 0 - lastSharedReverse - 1
	
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

	# TODO: Check all this, I've only done a first pass on this to remove seq1/seq2 refs
	for eachPosition in range(primerCutoff,lastShared - motifLen):
		currentRegion = seq[eachPosition : motifLen + eachPosition]
		
		regionHam = hamming(currentRegion,currentMotif)
		
		if regionHam == 0:
			untenablePositions.append(eachPosition)
		
		if (regionHam <= hammingThreshold):
			suitablePositions.append(eachPosition)
		
	# ===============================
	# Check to see if the current motif is present on the right side of the sequence (defined to be downstream here)
	rightSeq = revComp(seq[lastSharedReverse:])
	rightLastShared = abs(lastSharedReverse)
	rightSharedMotif = False
	
	# Figure out limit of right side
	rightSideCutoff = max(0,rightLastShared - ampliconLength - len(tempPrimer[0]))	
	
	for eachPosition in range(rightSideCutoff,rightLastShared-motifLen):
		currentRightRegion = rightSeq[eachPosition : motifLen + eachPosition]
		
		rightRegionHam = hamming(currentRightRegion,currentMotif)
		
		if 0 in [rightRegionHam]:
			rightSharedMotif = True
	
	# TODO: might be useful to have a way to 
	if rightSharedMotif == True:
		# This motif is unsuable from this direction
		return([[0],[0]])
	
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

def evaluateMutations(seq,targetSeq,enzymeDict,enzymeName):
	"""
	Top-level function that searches a sequence to be 
	edited with CRISPR/Cas9 methods for sites likely to be 
	usable as screening sites.
	"""
	# ALGO: check enzyme for case where it cuts in WT. if it doesn't, skip it. if it does, check the cut position. if the cut position is flanked by Ns, skip checking. if the cut position isn't, make probable edits.

	# This calls scanSequence(), which returns a two-element list of form:
	# [untenable positions, suitable positions]
	
	# Scan targetSeq across seq to find matching position
	# I assume there is only one cut site because a check should have performed for multiple match sites before it ever got this far. Also assuming there is exactly one because of same pre-checks.
	for eachPosition in range(0,len(seq)-(len(target)-1)):
		currentSubset = seq[eachPosition:(eachPosition+len(target))]
		comparisons = hamming(currentSubset,target,True,True)
		if 0 in comparisons:
			break
	
	# Identify cut site
	# Hard assumption that the cut site is at the -3 position from the 3' end of the provided 5'->3' target sequence
	if comparisons[0] == 0:
		cutPosition = eachPosition + len(target) - 3
	elif comparisons[1] == 0 or comparisons[2] == 0:
		cutPosition = eachPosition + 3
	
	
	# Check each enzyme for cutting
	for eachEnzyme in enzymeDict:
		enzymeInfo = enzymeDict[eachEnzyme]
		enzymeName = eachEnzyme
		currentMotif = enzymeInfo[0]
		rightSeqNum = 0
		leftSeqNum = 0
		rightSeqs = []
		leftSeqs = []
		canCutLeft = False
		canCutRight = False
		
		# Check if this enzyme cuts in the WT
		sitesLeft = scanSingle(seq,cutPosition,currentMotif,'left')
		sitesRight = scanSingle(seq,cutPosition,currentMotif,'right')
		
		# See if a cut site was found
		if sitesLeft[1] != []:
			canCutLeft = True
		# see if it cuts on the right side
		if sitesRight[1] != []:
			canCutRight = True
		
		# If either can cut, work with this system
		if canCutLeft or CanCutRight:
			alteredSeq = 
			# Take the sequence that can be cut (extant or modified) and mutate it
			altSeqs = crisprEdit(alteredSeq,cutPosition,organism)
			seqNum = len(altSeqs)
			for eachSequence in altSeqs:
				sitesLeft = scanSequence(alteredSeq,eachSequence,currentMotif,'left')
				sitesRight = scanSequence(alteredSeq,eachSequence,currentMotif,'right')
				# each is form of [[untenable positions],[possible positions]]
				
				if sitesLeft[1] != []:
					leftSeqNum += 1
					leftSeqs.append(sitesLeft)
				if sitesRight[1] != []:
					rightSeqNum += 1
					rightSeqs.append(sitesRight)
			# TODO: need to do some kind of counting of 'eachEnzyme in sitesLeft|sitesRight'
			if (rightSeqNum/seqNum) >= seqThreshold:
				x=1
				# return output
			if (leftSeqNum/seqNum) >= seqThreshold:
				x=1
				# return output
	return(None)

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
	# This calls scanSequence(), which returns a two-element list of form:
	# [untenable positions, suitable positions]
	# Untenable positions are positions in the upstream shared region where there is an exact match
	# The primer will need to edit these positions in order to work.
	# Suitable positions are positions where there is either an exact match for a motif or a match below the threshold
	# Suitable positions mark places that are diagnostic for one sequence or the other
	
	# If [[0],[0]] is returnd, then the motif cannot be used for these sequences because of a downstream exact match

	# TODO: get Settings object
	
	currentMotif = enzymeInfo[0]
	# Search from Left
	sitesLeft = scanSequence(seq1,seq2,currentMotif,'left')
	# Is [ [untenable positions] , [suitable positions] ]
	# Search from Right
	sitesRight = scanSequence(seq1,seq2,currentMotif,'right')
	# Is [ [untenable positions] , [suitable positions] ]
	
	output = []
	motifLen = len(currentMotif)
	for eachSet in [sitesLeft,sitesRight]:
		currentSet=[[],[]]
		if eachSet == [[0],[0]] or eachSet == [[],[]]:
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
			rejectedPrimers = 0
			for eachIndex in currentSet[1]:
				newPrimer = generatePrimer(currentSeq1,currentSet[0],eachIndex,currentLastShared,currentMotif)
				# TODO: check hamdist? or does generatePrimer do that
				if newPrimer is not None:
					lastPrimerBase = newPrimer[0][-1]
					lastSequenceBase = currentSeq1[currentLastShared-1]
					currentOut.append(" "*(12+currentLastShared-len(newPrimer[0]))+newPrimer[0])
					currentOut.append(estimateTM(newPrimer[0]))
				else:
					rejectedPrimers += 1
					continue
			output.append(currentOut)
			if rejectedPrimers == len(currentSet[1]):
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
	originalSeq = seq
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
	if lastPrimerBase != lastSeqBase and Settings.allowMismatch == True:
		# get the template base, which is the complement of the listed base
		templateBase = revComp(lastSeqBase).lower()
		templateBase = templateBase.lower()
		# compare to the 3' base
		# if template is G and primer 3' is T, ok.
		# Allowable: T/T, T/C, T/G
		if templateBase == 'g' and lastPrimerBase == 't':
			print('allowing a mismatching primer')
			return(bestPrimer)
		elif templateBase == 't' and lastPrimerBase in ['t','c','g']:
			print('allowing a mismatching primer')
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


# overall region to test:
# [lastShared - (motifLen-2) + iterator] to [lastShared + 1 + iterator]
# iterator is [0 .. motifLen-2]
# shared region to test is [lastShared - (motifLen-2) + iterator] to [lastShared]
# unshared region to test is [lastShared + 1] to [lastShared + 1 + iterator]
# motifLen is the number of characters in the motif
# lastShared is the index of the last shared character

# At each iteration, check whether: hammingdistance in unshared is 0 in one sequence and not the other, and hamming distance in shared distance in both sequences is less than threshold

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

# mismatches! G/T mismatches are ok. G/A and G/G mismatches are not good! will not amplify!

# so in order to account for indels i really can only butt up against them, i can't do what you do for SNPs which is put them right in the middle