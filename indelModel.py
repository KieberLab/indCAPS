import numpy


def indelModel(species):
	return(None)
	# choose whether insertions or deletions happen
	# Using negative binomial with a high probability of stopping
	# Calling for 1000 samples even though it will almost certainly be less than 10 indels
	if species is "ArabidopsisThaliana":
		n = 1000
		p = 0.95
	elif species is "OryzaSativa":
		n = 1000
		p = 0.90
	return(numpy.random.binomial(n,p))