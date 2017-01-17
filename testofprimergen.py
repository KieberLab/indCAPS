# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 16:18:14 2017

@author: hodgens
"""



seq1 = 'GCGGAgggcccCTCAAGATCCCGAGTgggTCTTATcccCAGTTTCTTGGCTCTGTTA' #arguments[1]
seq2 = 'GCGGAgggcccCTCAAGATCCCGAGTgggcccCAGTTTCTTGGCTCTGTTA' #arguments[2]
hamNum = 4 #arguments[3]

x = generatePrimer(seq1,[5,6],26,29,60,hamNum,"ggnncc")
print(x)
print(estimateTM(x[0]))