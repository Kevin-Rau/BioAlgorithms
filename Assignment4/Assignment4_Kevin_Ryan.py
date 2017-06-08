######################################################################################
#	File 			: Assignment3_Kevin_Ryan.py
#	Purpose         : Computing Global and Local Alignments.
#					: Create GLOBAL alignments between two sequences using a dynamic programming alignment algorithm 
#					: Create a distance matrix for a set of sequences using LOCAL alignment.
#	Developers      : Kevin Rau & Ryan O'Connell March 2017
#	Python 			: Version 2.7.10
######################################################################################
#   
#   How to run this application with command line arguments:
#   python Assignment3_Kevin_Ryan.py -f <filename>
#	The -f specifies the fasta file and loction to be read and analyzed.
#
#	***** EXAMPLE RUN *****
#
#	python Assignment2_Kevin_Ryan.py -f <filename>
#
######################################################################################
# 
#   References: Kevin Rau Jan 2017
#               Kevin_Rau_Assignment2.py
#				Formatting and flag handling & File Parsing 
#
######################################################################################

import sys, getopt, os.path, operator, re, math
import numpy as np
from numpy import array

######################################################################################
#
#	Nucleobase to number
#	This section of code sets up the nucleobase to number for matrix handling.
#	Also sets up the default scoring matrix to be used.
#
######################################################################################

A, C, G, T = 0, 1, 2, 3
int_to_char = {0:'A', 1:'C', 2:'G', 3:'T'}

indel = long(1)
scoring = []

######################################################################################
#
#	AlignmentFinder class
#	This class handles the global sequqnce comparison.
#	The class has six definitions:
#		init: Takes 3 parameters - self, seq1, seq2
#		find_global_alignment: Takes 1 parameter - self
#		_computer_array: Takes 1 parameter - self
#		_get_score: Takes 3 parameters - self, i, j
#		_get_aligned_pair: Takes 3 parameters - self, i, j
#		_traceback: takes 1 paremeter1 - self
#
######################################################################################

class AlignmentFinder(object):

		# Sets up the parameters for the class calls to use
    def __init__(self, seq1, seq2):
        self.seq1 = seq2 
        self.seq2 = seq1
        self.D = None

        # Finds the global alignment as compared to the seq1 and seq2 data
    def find_gobal_alignment(self):
        self.D = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int16)
        self._compute_array()
        #print self.D
        return self._traceback(), self.D # continue the traceback

        # Computes the array. Changes depending on the length of the sequence
    def _compute_array(self):
        for i in xrange(self.seq1.size+1):
            self.D[i,0] = i*long(indel)
        for j in xrange(self.seq2.size+1):
            self.D[0,j] = j*long(indel)
        for i in xrange(1, self.seq1.size+1):
            for j in xrange(1, self.seq2.size+1):
                self.D[i,j] = max(  self.D[i-1, j-1] + self._get_score(i, j),
                                    self.D[i-1, j] + long(indel),
                                    self.D[i, j-1] + long(indel))

    	# Adjustments for matrix
    def _get_score(self, i, j):
        #	The indexing is can be tricky. This is because the matrix as one more row and column.
        #	To obtain the correct nucleotide in the sequence, we must
        #	substract 1 to the matrix index.
        return scoring[self.seq1[i-1], self.seq2[j-1]]
    
    	# Checks for the aligned pair. Then returns the next one
    def _get_aligned_pair(self, i, j):
        n1 = int_to_char[self.seq1[i-1]] if i>0 else '-'
        n2 = int_to_char[self.seq2[j-1]] if j>0 else '-'
        return (n1, n2)

        # The traceback function to find the best alignment
    def _traceback(self):
        alignment= []	# Array to return with best alignment
        i = self.seq1.size
        j = self.seq2.size
        while i >0 and j>0:
            if self.D[i-1, j-1] + self._get_score(i, j) == self.D[i,j]:
                alignment.append(self._get_aligned_pair(i, j))
                i -= 1
                j -= 1
            elif self.D[i-1, j] + long(indel) == self.D[i,j]:
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
        alignment.reverse()	# Need to reverse it because it finds it backwards.
        return alignment  

######################################################################################
#
#	uage()
#	This function is used when the user needs help running the program.
#	The user just needs to pass in the -h flag to see this function.
#
######################################################################################

def usage():
    
	# Print the usage of the application
    print 'python overlap.py -f <filename>'
    print "-f specifies the FASTA file containing sequencing reads"

######################################################################################
#
#	main(argv)
#	This is the main function that handles the flow of the application.
#	It takes in 1 parameter:
#		argv = The arguements that the user passes into the application.
#
######################################################################################

def main(argv):
	
	inputFile = ''   # The input file name
	gap_penalty = 0
	scoring_matrix = ''

	try:
		# go through the input parameters, make sure there is a -f
		opts, args = getopt.getopt(argv,"hf:s:g:t:")
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
		if opt == '-s':
			if os.path.isfile(arg):
				scoring_matrix = arg
			# If the input file given does not exist, print an error message and exit the program
			else:
				print "ERROR! Input scoring matrix must exist. Correct usage:"
				usage()
				sys.exit(2)
		if opt == '-g':
			gap_penalty = arg
		if opt == '-t':
			newick = arg


	#print scoring_matrix
	global indel
	indel = long(gap_penalty)

	global scoring
	scoringNew = array((split((fix_scoring(scoring_matrix)),4)))
	scoring = scoringNew
	(genome, names) = readInput(inputFile)
	distanceMatrix2(genome, names)
	matrix = distanceMatrix(genome)
	M = stripMatrix(matrix)
	format = UPGMA(M, names)
	print("\nTree in Newick Format:\n")
	print format
	print ''
	text_file = open(newick, "w")
	text_file.write(format)
	text_file.close()

def split(arr, size):
     arrs = []
     while len(arr) > size:
         pice = arr[:size]
         arrs.append(pice)
         arr   = arr[size:]
     arrs.append(arr)
     return arrs

def fix_scoring(scoring):
	f = open(scoring, 'r')
	matrix = []
	new_matrix = []
	for line in f.readlines():
		matrix.append(line[2:].replace('\t', ',').replace('\n', '').split(','))
	for list in matrix:
		for number in list:
			new_matrix.append(int(number))
	return new_matrix


def stripMatrix(matrix):
	final_array = []
	for i in xrange(len(matrix)):
		for x in xrange(len(matrix)):
			if matrix[i][x] == 0.0:
				final_array.append(matrix[i][0:x])
	return final_array

def lowest_cell(table):
    # Set default to infinity
    min_cell = float("inf")
    x, y = -1, -1

    # Go through every cell, looking for the lowest
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min_cell:
                min_cell = table[i][j]
                x, y = i, j

    # Return the x, y co-ordinate of cell
    return x, y


# join_labels:
#   Combines two labels in a list of labels
def join_labels(labels, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # Join the labels in the first index
    labels[a] = "(" + labels[a] + "," + labels[b] + ")"

    # Remove the (now redundant) label in the second index
    del labels[b]


# join_table:
#   Joins the entries of a table on the cell (a, b) by averaging their data entries
def join_table(table, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # For the lower index, reconstruct the entire row (A, i), where i < A
    row = []
    for i in range(0, a):
        row.append((table[a][i] + table[b][i])/2)
    table[a] = row
    
    # Then, reconstruct the entire column (i, A), where i > A
    #   Note: Since the matrix is lower triangular, row b only contains values for indices < b
    for i in range(a+1, b):
        table[i][a] = (table[i][a]+table[b][i])/2
        
    #   We get the rest of the values from row i
    for i in range(b+1, len(table)):
        table[i][a] = (table[i][a]+table[i][b])/2
        # Remove the (now redundant) second index column entry
        del table[i][b]

    # Remove the (now redundant) second index row
    del table[b]


# UPGMA:
#   Runs the UPGMA algorithm on a labelled table
def UPGMA(table, labels):
    # Until all labels have been joined...
    while len(labels) > 1:
        # Locate lowest cell in the table
        x, y = lowest_cell(table)

        # Join the table on the cell co-ordinates
        join_table(table, x, y)

        # Update the labels accordingly
        join_labels(labels, x, y)

    # Return the final label
    return labels[0]



## A test using an example calculation from http://www.nmsr.org/upgma.htm

# alpha_labels:
#   Makes labels from a starting letter to an ending letter
def alpha_labels(start, end):
    labels = []
    for i in range(ord(start), ord(end)+1):
        labels.append(chr(i))
    return labels

# Test table data and corresponding labels

######################################################################################
#
#	iterateOverFastaFile2(a,b, g)
#	This function will call set up the seq and seq2 arrays
#	It takes in 3 parameters:
#		a = This is the index number to be used with array g for seq.
#		b = This is the index number to be used with array g for seq2.
#		g = This is the array that will be indexed.
#
######################################################################################

def iterateOverFastaFile2(a,b, g):
	seq = list(g[a])
	seq2 = list(g[b])	
	return seq, seq2

######################################################################################
#
#	computeMultipleAlignment2(sequence, sequence2)
#	This function will compute the mutiple alignment.
#	Same as the one above but only returns the top and bottom.
#	It takes in 2 parameters:
#		sequence = This the first sequence to be computed in the multiple alignment.
#		sequence2 = This the second sequence to be computed in the multiple alignment.
#
######################################################################################

def computeMultipleAlignmentPart2(sequence, sequence2):
	seq_a, seq_b = convertListToArray(sequence, sequence2)
	s1 = array(seq_a, dtype=np.int16)
	s2 = array(seq_b, dtype=np.int16)
	aligner = AlignmentFinder(s1, s2)
	pairs, matrix = aligner.find_gobal_alignment()
	top, bottom = formatSequences(pairs)
	return top, bottom

def distanceMatrix(sequences):
	dist_matrix = [[0 for a in range(len(sequences))]for b in range(len(sequences))]
	total_alignments = []
	for i in range(0, len(sequences)):
		for j in range(0, len(sequences)):
			matches = 0
			set1, set2 = iterateOverFastaFile2(i, j, sequences)
			total_alignments = computeMultipleAlignmentPart2(set1, set2)
			for k in range(0, len(total_alignments[0])):
				if total_alignments[0][k] == total_alignments[1][k]:
					matches += 1
				else:
					pass
			new_float = 1 - (float(matches) /len(total_alignments[0]))

			five_point = round(new_float, 5)

			dist_matrix[i][j] = five_point
	return dist_matrix

######################################################################################
#
#	distanceMatrix(sequences)
#	This function calculates the scoring matrix for all the sequences.
#	It takes in 1 parameter:
#		sequences = This is list that holds each sequqnce in the fasta file.
#
######################################################################################

def distanceMatrix2(sequences, names):
	print("\nDistance Matrix Between All Sequences in Fasta File:\n")
	dist_matrix = [[0 for a in range(len(sequences) + 1)]for b in range(len(sequences) + 1)]
	for i in range(0, len(sequences)):
		if i == 0:
			dist_matrix[0][i] = ''
		dist_matrix[0][i + 1] = names[i]
		for j in range(0, len(sequences)):
			if j == 0:
				dist_matrix[j][0] = ''
			dist_matrix[j + 1][0] = names[j]
	total_alignments = []
	for i in range(0, len(sequences)):
		for j in range(0, len(sequences)):
			matches = 0
			set1, set2 = iterateOverFastaFile2(i, j, sequences)
			total_alignments = computeMultipleAlignmentPart2(set1, set2)
			for k in range(0, len(total_alignments[0])):
				if total_alignments[0][k] == total_alignments[1][k]:
					matches += 1
				else:
					pass
			new_float = 1 - (float(matches) /len(total_alignments[0]))
			five_point = round(new_float, 5)
			five_point_zero = ('%.5f'%five_point)
			dist_matrix[i + 1][j + 1] = five_point_zero

	#print dist_matrix
	s = [[str(e) for e in row] for row in dist_matrix]
	lens = [max(map(len, col)) for col in zip(*s)]
	fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
	table = [fmt.format(*row) for row in s]
	print '\n'.join(table)

######################################################################################
#
#	convertListToArray(sequence, sequence2)
#	This function converts the string array into a char array.
#	It takes in 2 parameters:
#		sequences = This is the sequnce that needs to converted to a char array.
#		sequences2 = This is the second sequnce that needs to converted to a char array.
#
######################################################################################

def convertListToArray(sequence, sequence2):
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

######################################################################################
#
#	formatSequences(pairs)
#	This function will format the top and the bottom sequences for alignment.
#	It takes in 1 parameter1:
#		pairs = This is a list of the top and bottom sequence to be paired together for alignment.
#
######################################################################################

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

######################################################################################
#
#	readInput(inFile)
#	This function reads the inputed FASTA file and parses it.
#	It takes in 1 parameter:
#		inFile = The FASTA file to open and search for.
#
######################################################################################

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
	return sequences, sequenceNames


######################################################################################
#
#	__name__
#	This function calls the main function and passes in the parameters.
#
######################################################################################

if __name__ == "__main__":
	main(sys.argv[1:])