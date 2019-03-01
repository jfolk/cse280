import numpy as np
import math


def patternToNumber(pattern):
	if len(pattern) == 0:
		return 0
	#get the last letter of k-mer
	symbol = pattern[-1]
	#get the k-mer without the last letter
	prefix = pattern[:-1]
	return 4 * patternToNumber(prefix) + symbolToNumber(symbol)


def symbolToNumber(symbol):
	if symbol == 'A' or symbol == 'a':
		return 0
	elif symbol == 'C' or symbol == 'c':
		return 1
	elif symbol== 'G' or symbol == 'g':
		return 2
	else:
		return 3




#given long sequence of base pairs, and k for the length of the k-mer
def sequenceToNumber(sequence, k):
	basePairs = {'A' : 0, 'a' : 0, 'C': 1, 'c': 1, 'G': 2, 'g': 2, 'T': 3, 't': 3}
	referenceArray = [0] * (4**k)

	print(referenceArray)
	kmerToNumber = 0
	kmer = sequence[0: k]
	firstBasePair = 0

	#generate kmertToNumber for first kmer base pairs 
	for idx, letter in enumerate(kmer):
		kmerToNumber += basePairs[letter] * (4 ** (k - idx - 1))


	referenceArray[kmerToNumber] += 1	


	#generate kmerToNumber without going through all k base pairs
	for i in range(k, len(sequence)):

		firstBasePair = basePairs[sequence[i - k]] * 4 ** (k - 1)
		print("firstBasePair: ", firstBasePair)
		#subtract the first base pair's value
		kmerToNumber -= firstBasePair
		kmerToNumber *= 4
		kmerToNumber += basePairs[sequence[i]]

		#save the value of the new first base pair in new k-mer to subtract later
		

		print("kmerToNumber: ", kmerToNumber)
		referenceArray[kmerToNumber] += 1


	return referenceArray

	




print(sequenceToNumber('AGTAGTTTTT', 3))