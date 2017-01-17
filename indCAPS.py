#!/usr/bin/env python2

#seq1 = 'GCGGAgggcccCTCAAGATCCCGAGTgggTCTTATcccCAGTTTCTTGGCTCTGTTA' #arguments[1]
#seq2 = 'GCGGAgggcccCTCAAGATCCCGAGTgggcccCAGTTTCTTGGCTCTGTTA' #arguments[2]
#currentMotif = 'gggccc'
hamNum = 4 #arguments[3]

# Sample sequences from Leah
seq1 = "TGTGTGTGCAGGGAGAAGCCAAATGTGGATTTTGACAGGGTGGACTCGTATGTGCATCAG"
seq2 = "TGTGTGTGCAGGGAGAAGCCAAATGTGGATTTGACAGGGTGGACTCGTATGTGCATCAG"

## Enzymes
# Will not include nicking enzymes or intron-encoding enzymes or homing enzymes
# Also does not include two-site cutters: BaeI, BsaXI, CspCI, BcgI
enzymes = {"AclI": ["AACGTT",2,4],
"HindIII": ["AAGCTT",1,5],
"SspI": ["AATATT",3,3],
"MluCI": ["AATT",0,4],
"Tsp509I": ["AATT",0,4],
"PciI": ["ACATGT",1,5],
"AgeI": ["ACCGGT",1,5],
"BspMI": ["ACCTGC",10,-8],
"BfuAI": ["ACCTGC",10,-8],
"SexAI": ["ACCWGGT",1,6],
"MluI": ["ACGCGT",1,5],
"BceAI": ["ACGGC",17,-14],
"HpyCH4IV": ["ACGT",1,3],
"HpyCH4III": ["ACNGT",3,2],
"AflIII": ["ACRYGT",1,5],
"SpeI": ["ACTAGT",1,5],
"BsrI": ["ACTGG",6,1],
"BglII": ["AGATCT",1,5],
"AfeI": ["AGCGCT",3,3],
"AluI": ["AGCT",2,2],
"StuI": ["AGGCCT",3,3],
"ScaI": ["AGTACT",3,3],
"ClaI": ["ATCGAT",2,4],
"BspDI": ["ATCGAT",2,4],
"NsiI": ["ATGCAT",5,1],
"AseI": ["ATTAAT",2,4],
"SwaI": ["ATTTAAAT",4,4],
"MfeI": ["CAATTG",1,5],
"PmlI": ["CACGTG",3,3],
"DraIII": ["CACNNNGTG",6,3],
"AleI": ["CACNNNNGTG",5,5],
"PvuII": ["CAGCTG",3,3],
"AlwNI": ["CAGNNNCTG",6,3],
"TspRI": ["NNCASTGNN",9,0],
"NdeI": ["CATATG",2,4],
"NlaIII": ["CATG",4,0],
"CviAII": ["CATG",1,3],
"FatI": ["CATG",0,4],
"MslI": ["CAYNNNNRTG",5,5],
"XcmI": ["CCANNNNNNNNNTGG",8,7],
"BstXI": ["CCANNNNNNTGG",8,4],
"PflMI": ["CCANNNNNTGG",7,4],
"NcoI": ["CCATGG",1,5],
"SmaI": ["CCCGGG",3,3],
"XmaI": ["CCCGGG",1,5],
"TspMI": ["CCCGGG",1,5],
"SacII": ["CCGCGG",4,2],
"MspI": ["CCGG",1,3],
"HpaII": ["CCGG",1,3],
"ScrFI": ["CCNGG",2,3],
"BssKI": ["CCNGG",0,5],
"StyD4I": ["CCNGG",0,5],
"BsaJI": ["CCNNGG",1,5],
"BslI": ["CCNNNNNNNGG",7,4],
"BtgI": ["CCRYGG",1,5],
"NciI": ["CCSGG",2,3],
"AvrII": ["CCTAGG",1,5],
"SbfI": ["CCTGCAGG",6,2],
"Bsu36I": ["CCTNAGG",2,5],
"EcoNI": ["CCTNNNNNAGG",5,6],
"BstNI": ["CCWGG",2,3],
"PspGI": ["CCWGG",0,5],
"StyI": ["CCWWGG",1,5],
"PvuI": ["CGATCG",4,2],
"BstUI": ["CGCG",2,2],
"EagI": ["CGGCCG",1,5],
"RsrII": ["CGGWCCG",2,5],
"BsiEI": ["CGRYCG",4,2],
"BsiWI": ["CGTACG",1,5],
"Hpy99I": ["CGWCG",5,0],
"MspA1I": ["CMGCKG",3,3],
"SgrAI": ["CRCCGGYG",2,6],
"BfaI": ["CTAG",1,3],
"XhoI": ["CTCGAG",1,5],
"PaeR7I": ["CTCGAG",1,5],
"TliI": ["CTCGAG",1,5],
"PstI": ["CTGCAG",5,1],
"DdeI": ["CTNAG",1,4],
"SfcI": ["CTRYAG",1,5],
"AflII": ["CTTAAG",1,5],
"SmlI": ["CTYRAG",1,5],
"AvaI": ["CYCGRG",1,5],
"BsoBI": ["CYCGRG",1,5],
"XmnI": ["GAANNNNTTC",5,5],
"EcoRI": ["GAATTC",1,5],
"AatII": ["GACGTC",5,1],
"ZraI": ["GACGTC",3,3],
"Tth111I": ["GACNNNGTC",4,5],
"PflFI": ["GACNNNGTC",4,5],
"PshAI": ["GACNNNNGTC",5,5],
"AhdI": ["GACNNNNNGTC",6,5],
"DrdI": ["GACNNNNNNGTC",7,5],
"Eco53kI": ["GAGCTC",3,3],
"SacI": ["GAGCTC",5,1],
"HinfI": ["GANTC",1,4],
"EcoRV": ["GATATC",3,3],
"MboI": ["GATC",0,4],
"Sau3AI": ["GATC",0,4],
"DpnII": ["GATC",0,4],
"BfuCI": ["GATC",0,4],
"DpnI": ["GATC",2,2],
"BsaBI": ["GATNNNNATC",5,5],
"TfiI": ["GAWTC",1,4],
"BstAPI": ["GCANNNNNTGC",7,4],
"SphI": ["GCATGC",5,1],
"SrfI": ["GCCCGGGC",4,4],
"NaeI": ["GCCGGC",3,3],
"NgoMIV": ["GCCGGC",1,5],
"BglI": ["GCCNNNNNGGC",7,4],
"AsiSI": ["GCGATCGC",5,3],
"HinP1I": ["GCGC",1,3],
"HhaI": ["GCGC",3,1],
"BssHII": ["GCGCGC",1,5],
"NotI": ["GCGGCCGC",2,6],
"Fnu4HI": ["GCNGC",2,3],
"Cac8I": ["GCNNGC",3,3],
"MwoI": ["GCNNNNNNNGC",7,7],
"NheI": ["GCTAGC",1,5],
"BmtI": ["GCTAGC",5,1],
"BlpI": ["GCTNAGC",2,5],
"TseI": ["GCWGC",1,4],
"ApeKI": ["GCWGC",1,4],
"Bsp1286I": ["GDGCHC",5,1],
"BamHI": ["GGATCC",1,5],
"HaeIII": ["GGCC",2,2],
"PhoI": ["GGCC",2,2],
"FseI": ["GGCCGGCC",6,2],
"SfiI": ["GGCCNNNNNGGCC",8,5],
"NarI": ["GGCGCC",2,4],
"KasI": ["GGCGCC",1,5],
"SfoI": ["GGCGCC",3,3],
"PluTI": ["GGCGCC",5,1],
"AscI": ["GGCGCGCC",2,6],
"ApaI": ["GGGCCC",5,1],
"PspOMI": ["GGGCCC",1,5],
"Sau96I": ["GGNCC",1,4],
"NlaIV": ["GGNNCC",3,3],
"KpnI": ["GGTACC",5,1],
"Acc65I": ["GGTACC",1,5],
"BstEII": ["GGTNACC",1,6],
"AvaII": ["GGWCC",1,4],
"BanI": ["GGYRCC",1,5],
"BaeGI": ["GKGCMC",5,1],
"BsaHI": ["GRCGYC",2,4],
"BanII": ["GRGCYC",5,1],
"RsaI": ["GTAC",2,2],
"CviQI": ["GTAC",1,3],
"BstZ17I": ["GTATAC",3,3],
"SalI": ["GTCGAC",1,5],
"ApaLI": ["GTGCAC",1,5],
"AccI": ["GTMKAC",2,4],
"Hpy166II": ["GTNNAC",3,3],
"Tsp45I": ["GTSAC",0,5],
"HpaI": ["GTTAAC",3,3],
"PmeI": ["GTTTAAAC",4,4],
"HincII": ["GTYRAC",3,3],
"BsiHKAI": ["GWGCWC",5,1],
"ApoI": ["RAATTY",1,5],
"NspI": ["RCATGY",5,1],
"BsrFI": ["RCCGGY",1,5],
"BsrFalphaI": ["RCCGGY",1,5],
"BstYI": ["RGATCY",1,5],
"HaeII": ["RGCGCY",5,1],
"CviKI1": ["RGCY",2,2],
"EcoO109I": ["RGGNCCY",2,5],
"PpuMI": ["RGGWCCY",2,5],
"SnaBI": ["TACGTA",3,3],
"BspHI": ["TCATGA",1,5],
"BspEI": ["TCCGGA",1,5],
"TaqalphaI": ["TCGA",1,3],
"NruI": ["TCGCGA",3,3],
"Hpy188I": ["TCNGA",3,2],
"Hpy188III": ["TCNNGA",2,4],
"XbaI": ["TCTAGA",1,5],
"BclI": ["TGATCA",1,5],
"HpyCH4V": ["TGCA",2,2],
"FspI": ["TGCGCA",3,3],
"MscI": ["TGGCCA",3,3],
"BsrGI": ["TGTACA",1,5],
"MseI": ["TTAA",1,3],
"PacI": ["TTAATTAA",5,3],
"PsiI": ["TTATAA",3,3],
"BstBI": ["TTCGAA",2,4],
"DraI": ["TTTAAA",3,3],
"PspXI": ["VCTCGAGB",2,6],
"BsaWI": ["WCCGGW",1,5],
"BsaAI": ["YACGTR",3,3],
"EaeI": ["YGGCCR",1,5],
"BssSI": ["CACGAG",1,1],
"BssSalphaI": ["CACGAG",1,1],
"BmgBI": ["CACGTC",3,3],
"BmrI": ["ACTGGG",11,-4],
"EciI": ["GGCGGA",17,-9],
"BsmFI": ["GGGAC",15,-14],
"EarI": ["CTCTTC",7,-4],
"AcuI": ["CTGAAG",22,-14],
"EcoP15I": ["CAGCAG",31,-27],
"BtgZI": ["GCGATG",16,-14],
"BsaI": ["GGTCTC",7,-5],
"HphI": ["GGTGA",13,-7],
"AlwI": ["GGATC",9,-5],
"NmeAIII": ["GCCGAG",27,-19],
"BtsI": ["CAGTG",7,0],
"MutI": ["CAGTG",7,0],
"FspEI": ["CC",14,-16],
"BccI": ["CCATC",9,-5],
"MboII": ["GAAGA",13,-7],
"BbsI": ["GAAGAC",8,-6],
"LpnPI": ["CCDG",14,-14],
"BseRI": ["GAGGAG",16,-8],
"PleI": ["GAGTC",9,-5],
"MlyI": ["GAGTC",10,-5],
"BsrDI": ["GCAATG",8,0],
"BciVI": ["GTATCC",12,-5],
"BspCNI": ["CTCAG",14,-7],
"BsgI": ["GTGCAG",22,-14],
"SapI": ["GCTCTTC",8,-4],
"BspQI": ["GCTCTTC",8,-4],
"FokI": ["GGATG",14,-13],
"BtsCI": ["GGATG",7,0],
"MnlI": ["CCTC",11,-6],
"BbvI": ["GCAGC",13,-12],
"BtsI": ["GCAGTG",8,0],
"BtsalphaI": ["GCAGTG",8,0],
"FauI": ["CCCGC",9,-6],
"HpyAV": ["CCTTC",11,-5],
"BsmBI": ["CGTCTC",7,-5],
"MspJI": ["CNNR",13,-13],
"BsmAI": ["GTCTC",6,-5],
"BcoDI": ["GTCTC",6,-5],
"BpmI": ["CTGGAG",22,-14],
"BpuEI": ["CTTGAG",22,-14],
"HgaI": ["GACGC",10,-10],
"SfaNI": ["GCATC",10,-9],
"MmeI": ["TCCRAC",26,-18],
"BseYI": ["CCCAGC",1,1],
"BsmI": ["GAATGC",7,1],
"BbvCI": ["CCTCAGC",2,2],
"Bpu10I": ["CCTNAGC",2,2],
"AciI": ["CCGC",1,1],
"BsrBI": ["CCGCTC",3,3]}

## Errors to Check For
# Non-bases in input
# Cut site outside amplicon
													

## Custom Functions
		
def hammingBool(seq1,seq2,allResults=False):
	# Instead of returning a singular value, returns lists of bools indicating whether there's a match
	return(False)
	
def hamming(seq1,seq2,allResults=False,allComparisons=True):
	# Calculates the hamming distance between two equal-length sequences
	# Expects a character string, will complain if it doesn't get it
	# Returns an integer
	
	# HUGE WARNING. Be careful how you set up sequences with N's and call this.
	# Ex: hamming('gggnnn','grcnnn',allResults=True) returns [4,0,0]
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
	# Takes a string, returns reverse complement
	seq1 = seq1.lower()
	seqNew = []
			
	# seqNew = ""
	# dnaForward = ["g","c","a","t"]
	# dnaReverse = ["c","g","t","a"]
	# for eachBase in seq1:
		# oldPosition = [y for x,y in zip(dnaForward,[0,1,2,3]) if x is eachBase][0]
		# seqNew = seqNew + dnaReverse[oldPosition]
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
	return(seqString.upper())

	
## Melting temperature calculations
def baseNumbers(seq):
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
# from 1. Sugimoto, N., Nakano, S. I., Yoneyama, M. & Honda, K. I. Improved thermodynamic parameters and helix initiation factor to predict stability of DNA duplexes. Nucleic Acids Res. 24, 4501–4505 (1996).

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

def nearestNeighbor(seq,sodiumConc,primerConc):
	# typical salt conc is [50 mM], keep it between [0.01 M] and [1 M] 
	# typical primer conc is 50 nM
	# keep it above 8 bases
	# From http://biotools.nubic.northwestern.edu/OligoCalc.html and cited sources
	R = 1.9872036 # kcal / (mole*kelvin)
	T = (1000*(-deltaH(seq))) / ((-deltaS(seq))+R*log(1/primerConc))-272.9+16.6*log(sodiumConc,10)
	return(T)


def estimateTM(seq,saltConc=0.05,primerConc=50*10**(-9),func='nearestNeighbor'):
	# Default function is 'nearestNeighbor'
	# see http://biotools.nubic.northwestern.edu/OligoCalc.html
	# Find number of each base
	if func is 'nearestNeighbor':
		return(nearestNeighbor(seq,saltConc,primerConc))
	
	# primerConc is 50 nM, pH 7.0
	
	# If less than 8 bases in oligo:
	
	# If 13 or fewer bases in oligo:
	
	# If 14 or more bases in oligo:

## Custom Classes
class comparison:
	# Holds the two sequences for comparison
	
	## Function Description
	# [alignSequences] looks at two unequal length sequences and finds the regions at each side where there is shared sequence
	# [overlapExtent] is called by [alignSequences] and is used to determine how long the shared region is in aligned sequences
	
	def __init__(self):
		return(None)
		# run alignSequences()
		# If no shared regions exist on either side, throw an error
		# Check Tm of shared sequences, warn user if Tm may be too low
	
	def alignSequences(self,seq1,seq2):
		return(None)
	
	def compareEnzyme(self,enzyme):
		return(None)
		# Takes in an enzyme, sees if its motif is similar to the end of the shared region/start of unshared region



# Find the shared and unshared region from both sides
	# See if there's any exact motif matches in the unshared region
# Pick an enzyme's motif
# Look at all substrings that have at least one shared and unshared base
# Find all cases where there's an exact match in the unshared region and downstream shared regions as applicable, and there's a match in the upstream shared region that is less than the hamming distance cutoff


def lastSharedBase(seq1,seq2,direction):
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
	
def scanUnshared(seq,currentMotif,lastShared,lastSharedReverse,hammingThreshold):
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
			elif leftHam <= hammingThreshold and rightHam == 0:
				seqUSpos.append(eachPosition+lastShared - motifLen + 1)

		else:
			currentSeq1 = seqUnshared[eachPosition:eachPosition+motifLen]
			currentHam = hamming(currentSeq1,currentMotif)
			
			if currentHam == 0:
				seqUSpos.append(eachPosition+lastShared - motifLen + 1)
				seqUSexact.append(eachPosition+lastShared - motifLen + 1)
	return([seqUSpos,seqUSexact])


def scanSequence(seq1,seq2,currentMotif,direction,hammingThreshold):
	# This is going to scan the sequence from left to right, with the assumption that 
	# left is the 5' side of the primer and anything right is the downstream direction
	# If 'right' is specified as the direction, the sequence is reversed before processing
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
	for eachPosition in range(0,lastShared - motifLen):
		currentRegion1 = seq1[eachPosition : motifLen + eachPosition]
		currentRegion2 = seq2[eachPosition : motifLen + eachPosition]
		
		region1Ham = hamming(currentRegion1,currentMotif)
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
	
	for eachPosition in range(0,rightLastShared-motifLen):
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
	seq1Results = scanUnshared(seq1,currentMotif,lastShared,lastSharedReverse,hammingThreshold)
	seq2Results = scanUnshared(seq2,currentMotif,lastShared,lastSharedReverse,hammingThreshold)
	
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
#		for eachPosition in newSuitable:
#			print(seq1[eachPosition:eachPosition+motifLen])
#			print(seq2[eachPosition:eachPosition+motifLen])
#			print(currentMotif)
		return([newUntenable,newSuitable])
	else:
#		for eachPosition in suitablePositions:
#			print(seq1[eachPosition:eachPosition+motifLen])
#			print(seq2[eachPosition:eachPosition+motifLen])
#			print(currentMotif)
		return([untenablePositions,suitablePositions])


	
def evaluateSites(seq1,seq2,enzymeInfo,hammingThreshold,enzymeName):
	# This calls scanSequence(), which returns a two-element list of form:
	# [untenable positions, suitable positions]
	# Untenable positions are positions in the upstream shared region where there is an exact match
	# The primer will need to edit these positions in order to work.
	# Suitable positions are positions where there is either an exact match for a motif or a match below the threshold
	# Suitable positions mark places that are diagnostic for one sequence or the other
	
	# If [[0],[0]] is returnd, then the motif cannot be used for these sequences because of a downstream exact match

	currentMotif = enzymeInfo[0]
	# Search from Left
	sitesLeft = scanSequence(seq1,seq2,currentMotif,'left',hammingThreshold)
	# Is [ [untenable positions] , [suitable positions] ]
	# Search from Right
	sitesRight = scanSequence(seq1,seq2,currentMotif,'right',hammingThreshold)
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
			
			if 0 in motifDiff1 or 0 in motifDiff2:
				currentOut.append("CAPS primer possible")
			# TODO: Make sure this is correct
			currentLastShared = lastSharedBase(currentSeq1,currentSeq2,'left')
			for eachIndex in currentSet[1]:
				newPrimer = generatePrimer(currentSeq1,currentSet[0],eachIndex,currentLastShared,58,hammingThreshold,currentMotif)
				currentOut.append(" "*(12+currentLastShared-len(newPrimer[0]))+newPrimer[0])
				currentOut.append(estimateTM(newPrimer[0]))
			output.append(currentOut)
			return(output)


def generatePrimer(seq,untenablePositions,desiredSuitable,lastShared,TM,hammingThreshold,currentMotif):
	# TM is a value for Degree C
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
	for each in range(0,lastShared-desiredSuitable):
		if orientedMotif[each] in ['g','c','t','a'] and orientedMotif[each] != currentSite[each]:
			currentSite[each] = orientedMotif[each]
	orientedMotif = ''.join(orientedMotif)
	currentSite = ''.join(currentSite)
	seqLeft = seq[:desiredSuitable]
	seqRight = seq[(desiredSuitable+motifLen):]
	seq = seqLeft + currentSite + seqRight
	
	# See if the primer violates any rules
	# TODO: rule check
	# Rule: G/T mismatch works at 3' end, but G/A and G/G mismatches don't. (Simsek 2000)
	# Rule: Mismatches to canonical shared sequence less than hamming dist
	
	# Generate primer sub-sequences
	while currentStart > 0:
		currentPrimer = seq[currentStart:lastShared]
		primerList.append(currentPrimer)
		currentStart -= 1
	for eachPrimer in primerList:
		output.append(estimateTM(eachPrimer))
	
	# Filter for negative values
	filterList = [x for x in range(0,len(output)) if output[x] <= 0]
	primerList = [primerList[x] for x in range(0,len(output)) if x not in filterList]
	output = [output[x] for x in range(0,len(output)) if x not in filterList]
	
	
	# Find minimum diff
	outputDiff = [abs(x - TM) for x in output]
	# Choose the best primer sub-sequence based on difference from optimum Tm
	bestPrimer = [primerList[x] for x in range(0,len(outputDiff)) if outputDiff[x] == min(outputDiff)]
	return(bestPrimer)



def compareEnzymes(seq1,seq2,enzymes,hammingThreshold):
	return(None)

def makePrimer(seq1,seq2,enzyme,suitablePositions,untenablePositions,hammingThreshold):
	return(None)

#for eachEnzyme in enzymes:
#	enzymeName = eachEnzyme
#	enzymeValue = enzymes[eachEnzyme]
#	results = evaluateSites(seq1,seq2,enzymeValue,1,enzymeName)
#	if results is not None:
#		for eachEnzyme in results:
#			for eachLine in eachEnzyme:
#				print(eachLine)


#x = generatePrimer(seq1,[5,6],26,29,60,hamNum,"ggnncc")
#print(x)
#print(estimateTM(x[0]))

# possibilities for unshared regions
# both sequences have unshared region of same length
# one sequence has unshared region and one doesn't
# ASSUMPTION: already checked that both sequences aren't identical
# REMEMBER TO: keep track of hamming distance count as we go

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