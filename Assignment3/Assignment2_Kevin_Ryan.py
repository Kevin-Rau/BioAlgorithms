######################################################################################
#	File 			: Assignment3_Kevin_Ryan.py
#	Purpose         : Computing Global and Local Alignments.
#					: Create GLOBAL alignments between two sequences using  a dynamic programming alignment algorithm 
#					: Create a distance matrix for a set of sequences using LOCAL alignment.
#	Developers       : Kevin Rau & Ryan O'Connel March 2017
#	Python 			: Version 2.7.10
######################################################################################
#   
#   How to run this application with command line arguments:
#   python Assignment3_Kevin_Ryan.py -F <filename>
#	The -f specifies the fasta file and loction to be read and analyzed.
#
#	***** EXAMPLE RUN *****
#
#	python Assignment2_Kevin_Ryan.py -F <filename>
#
######################################################################################
# 
#   References: Kevin Rau Jan 2017
#               Kevin_Rau_Assignment2.py
#				Formatting and flag handling & File Parsing 
#
######################################################################################
#!/usr/bin/python -O

import sys, getopt, os.path, operator, re

from collections import defaultdict

import itertools


import numpy as np
from numpy import array

A, C, G, T = 0, 1, 2, 3
int_to_char = {0:'A', 1:'C', 2:'G', 3:'T'}



indel = -2 

scoring = array([[3,-2,-1,-2],
                 [-2,3,-2,-1],
                 [-1,-2,3,-2],
                 [-2,-1,-2,3]])


class AlignmentFinder(object):
    def __init__(self, seq1, seq2):
        self.seq1 = seq2
        self.seq2 = seq1
        self.D = None

    def find_gobal_alignment(self):
        self.D = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int16)
        self._compute_array()
        print self.D
        return self._traceback()

    def _compute_array(self):
        for i in xrange(self.seq1.size+1):
            self.D[i,0] = i*indel
        for j in xrange(self.seq2.size+1):
            self.D[0,j] = j*indel
        for i in xrange(1, self.seq1.size+1):
            for j in xrange(1, self.seq2.size+1):
                self.D[i,j] = max(  self.D[i-1, j-1] + self._get_score(i, j),
                                    self.D[i-1, j] + indel,
                                    self.D[i, j-1] + indel)
    def _get_score(self, i, j):
        ''' The indexing is quite tricky because the matrix as one more row & column.
        That causes a shift between the matrix index and the sequence indices.
        Therefore, to obtain the correct nucleotide in the sequence, we must
        substract 1 to the matrix index. '''
        return scoring[self.seq1[i-1], self.seq2[j-1]]
    
    def _get_aligned_pair(self, i, j):
        n1 = int_to_char[self.seq1[i-1]] if i>0 else '-'
        n2 = int_to_char[self.seq2[j-1]] if j>0 else '-'
        return (n1, n2)

    def _traceback(self):
        alignment= []
        i = self.seq1.size
        j = self.seq2.size
        while i >0 and j>0:
            if self.D[i-1, j-1] + self._get_score(i, j) == self.D[i,j]:
                alignment.append(self._get_aligned_pair(i, j))
                i -= 1
                j -= 1
            elif self.D[i-1, j] + indel == self.D[i,j]:
                alignment.append(self._get_aligned_pair(i, 0))
                i -= 1
            else:
                alignment.append(self._get_aligned_pair(0, j))
                j -= 1
        while i > 0:
            alignment.append(self._get_aligned_pair(i, 0))
            i -= 1
        while j > 0:
            alignment.append(self._get_aligned_pair(0, j))
            j -= 1
        alignment.reverse()
        return alignment  

def usage():
    
	# Print the usage of the application
    print 'python overlap.py -F <filename>'
    print "-f specifies the FASTA file containing sequencing reads"

def main(argv):
	
	inputFile = ''   # The input file name
	
	try:
		# go through the input parameters, make sure there is a -f
		opts, args = getopt.getopt(argv,"hf:")
	# If one of the arguments is missing, print an error message and exit the program
	except getopt.GetoptError:
		print "ERROR! An input file (-f). Correct usage:"
		usage()
		sys.exit(2)
	for opt, arg in opts:
		# process each input parameter
		if opt == '-h':
			print "Correct usage is:"
			usage()
			sys.exit()
		if opt == '-f':
			if os.path.isfile(arg):
				inputFile = arg
			# If the input file given does not exist, print an error message and exit the program
			else:
				print "ERROR! Input file must exist. Correct usage:"
				usage()
				sys.exit(2)
	genome = readInput(inputFile)
	completeWork(genome)

def completeWork(genome):
	i = 0
	one = 0
	two = 1
	for i in xrange(len(genome)/2):

		iterateOverFastaFile(one,two, genome)
		one +=2
		two +=2


def iterateOverFastaFile(a,b, g):
	
	seq = list(g[a])
	seq2 = list(g[b])	

	computeMultipleAlignment(seq, seq2)


def computeMultipleAlignment(sequence, sequence2):

	seq_a, seq_b = convertListToArray(sequence, sequence2)
	s1 = array(seq_a, dtype=np.int16)
	s2 = array(seq_b, dtype=np.int16)
	aligner = AlignmentFinder(s1, s2)
	pairs = aligner.find_gobal_alignment()
	top, bottom = formatSequences(pairs)
	printSequences(top, bottom)

def printSequences(t, b):

	total_len = len(t)
	total_len_lines = len(t) * '='
	print t
	print b
	print total_len_lines,'(',total_len,')'


def convertListToArray(sequence, sequence2):

	temp = []


	i = 0
	for i in xrange(len(sequence)):
		if sequence[i] == 'A':
			sequence[i] = A
		if sequence[i] == 'T':
			sequence[i] = T
		if sequence[i] == 'G':
			sequence[i] = G
		if sequence[i] == 'C':
			sequence[i] = C

	x = 0
	for x in xrange(len(sequence2)):
		if sequence2[x] == 'A':
			sequence2[x] = A
		if sequence2[x] == 'T':
			sequence2[x] = T
		if sequence2[x] == 'G':
			sequence2[x] = G
		if sequence2[x] == 'C':
			sequence2[x] = C


	return sequence, sequence2


def formatSequences(pairs):
    top_seq = []
    bottom_seq = []
    formatTop = []
    formatBottom = []
    for (b, t) in pairs:
        bottom_seq.append(b)
        top_seq.append(t)

    formatTop = ''.join(top_seq)
    formatBottom = ''.join(bottom_seq)

    return formatTop, formatBottom

#This will parse the Fasta File and then return a list of sequences only 
def readInput(inFile):
    
	f = open(inFile, 'r')
	sequenceNames = []
	sequences = []
	sequenceSoFar = []
	output = []
 
	for line in f.xreadlines():
		line = line.strip('\n')
		if line != '' and line[0] == '>':
			# Separate each sequence header by spaces and store the name (string before the first space)
			name = line[1:].split()
			sequenceNames.append(name[0])
			# Store the sequence string corresponding to the previous sequence in the file
			# Ensure the corresponding name and sequence stored with the same index
			if sequenceSoFar != []:
				sequences.append(''.join(sequenceSoFar))
				sequenceSoFar = []
		else:
			# Collect all the lines of a sequence in here
			sequenceSoFar.append(line)
   
	# Add the last sequence string into sequences
	sequences.append(''.join(sequenceSoFar))

	return sequences

if __name__ == "__main__":
	main(sys.argv[1:])