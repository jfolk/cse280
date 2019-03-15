#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import pysam
import json


# In[2]:


#given long sequence of base pairs, and k for the length of the k-mer
#output a dictionary that will keep track of the counts of each k-mer seen, using their numeric value
def sequenceToNumber(sequence, k):
    basePairs = {'A' : 0, 'a' : 0, 'C': 1, 'c': 1, 'G': 2, 'g': 2, 'T': 3, 't': 3}
    
    referenceArray = {}
    
    
    kmerToNumber = 0
    kmer = sequence[0: k]
    firstBasePair = 0

    for idx, letter in enumerate(kmer):
        kmerToNumber += basePairs[letter] * (4 ** (k - idx - 1))
        
    if kmerToNumber in referenceArray:
        referenceArray[kmerToNumber] += 1
    else:
        referenceArray[kmerToNumber] = 1
        

    for i in range(k, len(sequence)):
        firstBasePair = basePairs[sequence[i - k]] * 4 ** (k - 1)
        # print("firstBasePair: ", firstBasePair)
        #subtract the first base pair's value
        kmerToNumber -= firstBasePair
        kmerToNumber *= 4
        kmerToNumber += basePairs[sequence[i]]
        
        #save the value of the new first base pair in new k-mer to subtract later

        # print("kmerToNumber: ", kmerToNumber)
        if kmerToNumber in referenceArray:
            referenceArray[kmerToNumber] += 1
        else:
            referenceArray[kmerToNumber] = 1
    
    return referenceArray


# In[3]:


def insertIntoGenome(genome, viralSequence, locations):
    locations.sort()
    finalLocs = []
    count = 0
    for l in locations:
        loc = l + count * len(viralSequence)
        finalLocs.append(loc)
        genome = genome[0:loc] + viralSequence + genome[loc:]
        count += 1
    f = open('chr11.fa', 'w+')
    f.write('>chr11\n')
    writeChunks = [genome[i: i+50] for i in range(0, len(genome), 50)]
    for w in writeChunks:
        f.write(w + "\n")
    return genome, finalLocs


# In[4]:


def createRefArray(sequence, k, fileName):
    referenceArray = {}
    
    for i in range(0, len(sequence)-k):
        string = ""
        if type(sequence) is list:
            string = "".join(sequence[i: i+k])
        else:
            string = sequence[i: i+k]
            
        if string in referenceArray:
            referenceArray[string] += 1
        else:
            referenceArray[string] = 1

#     with open(fileName, 'w+') as file:
#         file.write(json.dumps(referenceArray))
    return referenceArray


# In[5]:


def flipSequence(sequence):
    output = ''
    for bp in sequence:
        if bp == 'A':
            output += 'T'
        elif bp == 'T':
            output += 'A'
        elif bp == 'G':
            output += 'C'
        elif bp == 'C':
            output += 'G'
        elif bp == 'N':
            output += 'N'
    return output[::-1]


# In[6]:


def readTestSeq(file, finalLocs, viralString, viralRef, humanRef):
    with open(file) as fp:
        numError = 0
        numErrorHybrid = 0
        numHybrid = 0 
        N = 0.0
        isHybrid = False
        readDirection = 1
        for line in fp:
#             print(line)            
            readNumber = 0
            if line.startswith("@chr"):
                isHybrid = False
                line = line.strip()
                readNumber = int(line[-1])
                readEnd = 0
                readHeader = line.split('_')
                readStart = int(readHeader[int(readNumber)])
                if int(readHeader[1]) > int(readHeader[2]):
                    if readNumber == 1:
                        #print('readDirection = -1')
                        readDirection = -1
                        readEnd = readStart + 150
#                         temp = readStart
#                         readStart = readEnd
#                         readEnd = temp
                    else:
                        #print('readDirection = 1')
                        readDirection = 1
                        readEnd = readStart + 150
                else:
                    if readNumber == 1:
                        #print('readDirection = 1')
                        readDirection = 1
                        readEnd = readStart + 150
                    else:
                        #print('readDirection = -1')
                        readDirection = -1
                        readEnd = readStart + 150
#                         temp = readStart
#                         readStart = readEnd
#                         readEnd = temp
                for loc in finalLocs:
                    if (readStart < loc <= readEnd) or (readStart <= loc + len(viralString) < readEnd):
                        print("hybridRead found")
                        isHybrid = True
                        numHybrid += 1
                        break
                
            if line[0] in validBasePairs:
                N += 1.0
                sequence = line
                #print('read direction ' + str(readDirection))
                if readDirection < 0:
                    #print('i flip')
                    sequence = flipSequence(line)
                
                count = 0.0
                predictHybrid = False
                ref = createRefArray(sequence, 15, "150bpRead")
                hasUniqueViral = False
                hasUniqueHuman = False
                for kmer in list(ref.keys()):
                    inViralRef = kmer.lower() in viralRef
                    inHumanRef = kmer.upper() in humanRef
                    
#                     if isHybrid is True:
#                         print(kmer.lower())
                    
                    #in viral but not in human
                    if inViralRef and not inHumanRef:
                        count += 1.0
                        hasUniqueViral = True
                    #in viral and in human
                    if inViralRef and inHumanRef:
#                         print("in both")
                        count += 0.5
                    #in human but not in viral
                    if inHumanRef and not inViralRef:
                        hasUniqueHuman = True
#                 if 135.0*(1.0/8.0) <= count <= 135 * .9:

                if hasUniqueViral and hasUniqueHuman:
#                     print("has unique viral and unique human")
                    predictHybrid = True
                if 5 <= count <= 140:
                    predictHybrid = True
#                     print("Hybrid Read")
                else:
                    predictHybrid = False
                    
                if predictHybrid != isHybrid:
                    print("count: " + str(count) + "\tpredictHybrid: " + str(predictHybrid) + "\tisHybrid: " + str(isHybrid))
                    if isHybrid is True:
                        print('hasUniqueViral: ' + str(hasUniqueViral) + '\thasUniqueHuman: ' + str(hasUniqueHuman))
                        print('read sequence ' + sequence)
                        print('read direction ' + str(readDirection))
                        numErrorHybrid += 1
                    numError += 1
                    
#         print('errorRate: ' + str(numError/N))
        return numError, numErrorHybrid, numHybrid, N


# In[7]:


def readJson(file):
    dictionary = {}
    with open(file, 'r') as fp:
        dictionary = json.load(fp)
        
    return dictionary


# In[8]:


validBasePairs = ['A', 'a', 'C', 'c', 'G', 'g', 'T', 't']


# In[9]:


filepath = 'hpv68'
viralSequence = None
viralString = ""
with open(filepath) as fp:
    for line in fp:
        line = line.strip()
        viralString += line.upper()
#         viralSequence = sequenceToNumber(line, 15)


# In[10]:


genomeFilepath = "../projectSequences/HumanGenome/chr11.fa"
genome = ""
humanStringSequence = []
humanSequence = None
with open(genomeFilepath) as fp:
    for line in fp:
        line = line.strip()
        if not line.startswith('>chr11'):
            for char in line:
                if char in validBasePairs:
                    humanStringSequence.append(char.upper())
                genome += char.upper()
#     humanSequence = sequenceToNumber(humanStringSequence, 15)

# print(humanSequence)


# In[26]:


g, finalLocs = insertIntoGenome(genome, viralString, [6000000, 600050, 1000000, 1100000, 1200000, 1600000, 1700000, 1800000, 1900000])


# In[21]:


print(str(len(viralString)))


# In[13]:


viralRef = None
with open('hpv68') as fp:
    for line in fp:
        viralRef = createRefArray(line, 15, "viralRefArray")


# In[21]:


# humanRef = createRefArray(humanStringSequence, 15, "humanRefArray1")


# In[22]:


# print(len(humanRef))


# In[14]:


readResult = readJson("humanRefArray")


# In[15]:


print("HumanRef: " + str(len(readResult)) + "\tViralRef: " + str(len(viralRef)))


# In[16]:


humanRef = readResult


# In[24]:


print(finalLocs)


# In[22]:


numError1, numErrorHybrid1, numHybrid1, N1 = readTestSeq("read.bwa.read1.fastq", finalLocs, viralString, viralRef, humanRef)


# In[ ]:


numError2, numErrorHybrid2, numHybrid2, N2 = readTestSeq("read.bwa.read2.fastq", finalLocs, viralString, viralRef, humanRef)


# In[ ]:


numError3, numErrorHybrid3, numHybrid3, N3 = readTestSeq("read.bwa.read1.fastq1", finalLocs, viralString, viralRef, humanRef)


# In[ ]:


k = flipSequence('AACAAACAAAAAAATTGCTCTTTATCAGGTGGGAGATATTCACTTCAATTCCTACTTTGTTGATTAGTATAGAGAACTGCTGTGTTCAGCTTTATATACACCGTTTTCGGTCGTGACCGTTTTCGGTCCCACCCTTTTTTTATATAGAAT')
print(k)


# In[ ]:


print("numError1: " + str(numError1))
print("numErrorHybrid1: " + str(numErrorHybrid1))
print("numHybrid1: " + str(numHybrid1))
print("N1: " + str(N1))


# In[19]:


print("numError2: " + str(numError2))
print("numErrorHybrid2: " + str(numErrorHybrid2))
print("numHybrid2: " + str(numHybrid2))
print("N2: " + str(N2))


# 
