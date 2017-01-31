#!/usr/bin/env python2
from indCAPS import lastSharedBase, revComp

def checkBases(seqToTest):
	if all([x in ['g','c','a','t'] for x in seqToTest.lower()]):
		return(None)
	elif any([x in ['r','y','s','w','k','m','b','d','h','v','n'] for x in seqToTest.lower()]):
		return('Degenerate/nonspecific bases or gaps present in input.')
	else:
		return('Non-bases present in input.')

def removeWhitespace(seqToTest):
	changed = False
	if any([x in ['\t',' '] for x in seqToTest]):
		seqToTest = ''.join([x for x in "gcat tacg" if x not in ['\t', ' ']])
		changed = True
	return([seqToTest,changed])

def nonBasePresent(seqToTest):
	# Returns true if there are non-bases
	if seqToTest == None:
		return(True)
	else:
		return(any([x in ['g','c','a','t','r','y','s','w','k','m','b','d','h','v','n'] for x in seqToTest.lower()]))

def evaluateInput(seq1,seq2=None):
	# Returns [seq1,seq2,notes] where notes is a list of warnings and notes on modifications
	notes = []
	
	# Test seq1 for non-bases
	whitespaceRemoveResults1 = removeWhitespace(seq1)
	if whitespaceRemoveResults1[1] is True:
		seq1 = whitespaceRemoveResults1[0]
		notes.append('Removed whitespace from Sequence 1.')
	notes1 = checkBases(seq1)
	if notes1 is not None:
		notes.append(notes1)
		
	# Test seq2 for non-bases and check last shared bases
	if seq2 is not None:	
		whitespaceRemoveResults2 = removeWhitespace(seq2)
		if whitespaceRemoveResults2[1] is True:
			seq2 = whitespaceRemoveResults2[0]
			notes.append('Removed whitespace from Sequence 2.')
		notes2 = checkBases(seq2)
		if notes2 is not None:
			notes.append(notes2)

		lastSharedLeft = lastSharedBase(seq1,seq2,'left')
		lastSharedRight = lastSharedBase(seq1,seq2,'right')
		if lastSharedLeft <= 25:
			notes.append('Less than 25 identical shared bases on the left side of the inputted sequences. Primer design may fail or may not reach sufficient length or Tm.')
		if lastSharedRight <= 25:
			notes.append('Less than 25 identical shared bases on the right side of the inputted sequences. Primer design may fail or may not reach sufficient length or Tm.')

	return([seq1,seq2,notes])
