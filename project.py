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



print(patternToNumber('AGT'))