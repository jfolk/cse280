import numpy as np
import math



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