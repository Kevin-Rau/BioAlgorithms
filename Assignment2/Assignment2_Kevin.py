######################################################################################
#	File 			: Assignment2_Kevin.py
#	Purpose         : Read a set of short sequences from a FASTA file.
#					: Find the most common k-mers appearing in all the sequences. 
#					: Find the most common k-mers appearing in the most sequences. 
#					: Find the speed and memory limitations of implementations.
#	Developer       : Kevin Rau Feb 2017
#	Python 			: Version 2.7.10
######################################################################################
#   
#   How to run this application with command line arguments:
#   python Assignment2_Kevin.py -F <filename> -L <size of the k-mer>
#	The -f specifies the fasta file and loction to be read and analyzed.
#	The -k is the length of the k-mer to find in the given fasta file.
#
#	***** EXAMPLE RUN *****
#
#	python Assignment2_Kevin.py -F <filename> -L <size of the k-mer>
#
######################################################################################
# 
#   References: Kevin Rau Jan 2017
#               Kevin_Rau_Assignment1.py
#				Formatting and flag handling & File Parsing 
#
######################################################################################
#
#	Time-Completxity: This will run in worst case O((n*m)^2)					
#					  because the slowest part of this will be when we need to 
#					  find all of the kmers in every sequence
#					  This will move n*m, n being the ength of the list (longest worst case)
#					  and m being the number of sequences
#					  plus this having to do the same work twice.
######################################################################################
######################################################################################

import sys, getopt, os.path, operator, re

from collections import defaultdict

import itertools
######################################################################################
#	usage() 
#	This function helps the user if they don't know how to use the program
######################################################################################
def usage():
    
	# Print the usage of the application
    print 'python overlap.py -F <filename> -L <kmer_size>'
    print "-F specifies the FASTA file containing sequencing reads"
    print "-L specifies the smallest allowed integer matching in the overlap regions (3 =< k <= 8)"

def main(argv):
	
	inputFile = ''   # The input file name
	kmer = 0   # The length of the kmer
	
	try:
		# go through the input parameters, make sure there is a -f, and -k option
		opts, args = getopt.getopt(argv,"hF:L:")
	# If one of the arguments is missing, print an error message and exit the program
	except getopt.GetoptError:
		print "ERROR! An input file (-F) and a length of a kmer (-L) must be specified. Correct usage:"
		usage()
		sys.exit(2)
	for opt, arg in opts:
		# process each input parameter
		if opt == '-h':
			print "Correct usage is:"
			usage()
			sys.exit()
		if opt == '-F':
			if os.path.isfile(arg):
				inputFile = arg
			# If the input file given does not exist, print an error message and exit the program
			else:
				print "ERROR! Input file must exist. Correct usage:"
				usage()
				sys.exit(2)
		if opt == '-L':
			# If the kmer length is not between 3 and 8 inclusive, print an error message and exit the program
			if int(arg)<3 or int(arg)>8:
				print "ERROR! Length of kmer (-L command) must be between 3 and 8 inclusive. Correct usage:"
				usage()
				sys.exit(2)
			else:
				kmer = int(arg)

	# Run the readInput function to format the input file.

	genome = readInput(inputFile)
	listNames = sequenceCount(genome, kmer)
	combined_list = list(itertools.chain.from_iterable(listNames))
	combined_list2 = list(itertools.chain.from_iterable(listNames))

	sequence_final_list = unique(combined_list2)

	totals = defaultdict(int)
	totals2 = defaultdict(int)
	# Here we are taking the built list of lists and making one list of all kmers
	# Then counting the number of occurences of each kmer in the list
	for x in sequence_final_list:
		totals2[x] += 1

	for x in combined_list:
		totals[x] += 1

	#Sort the counts from largest to smallest 
	sorted_kmer = sorted(totals.items(), key=operator.itemgetter(1), reverse=True)
	sorted_kmer2 = sorted(totals2.items(), key=operator.itemgetter(1), reverse=True)

	#Make calls to helper functions 
	top_five_list, top_five_list_names, top_five_list_occurences = findTopFive(sorted_kmer)
	top_sequencesCount, top_sequencesKmers = sequenceFrequency2(sorted_kmer2, genome)
	temp = dict(zip(top_sequencesKmers, top_sequencesCount))
	final_sequence_list = sorted(temp.items(), key=operator.itemgetter(1), reverse=True)
	top_sequences2, top_sequences2_names, top_sequences2_occurences = findTopFive(final_sequence_list)
	formatPrint(top_five_list_names, top_sequences2_names, top_sequences2_occurences, top_five_list_occurences, inputFile, kmer)

def unique(items):
    found = set([])
    keep = []

    for item in items:
        if item not in found:
            found.add(item)
            keep.append(item)

    return keep

def formatPrint(Kmers,Kmers2, Frequences, Occurances, file, length):
	print ""
	print "File:\t", file
	print "k-mer Length:\t", length
	print ""

	print "k-mers", "\t", "Occurances"
	for x in xrange(0, len(Kmers)):
		print Kmers[x], "\t", Occurances[x]


	print "k-mers", "\t", "Seq Count"
	for x in xrange(0, len(Kmers2)):
		print Kmers2[x], "\t", Frequences[x]

#Parameters
#	my_list = the sorted full list of the counts of kmers in all the sequences
def findTopFive(my_list):

	final_list = []
	grabKmer = []
	list_kmers = []
	list_occur = []

	count = 5

	for i in xrange(0, 10):
		if my_list[4+i][1] == my_list[5+i][1]:
			count += 1
		if my_list[4+i][1] != my_list[5+i][1]:
			break

	for x in xrange(0, count):

		grabKmer = my_list[x]
		final_list.append(grabKmer)


	for i in xrange(0, len(final_list)):
		list_kmers.append(final_list[i][0])
		
	for i in xrange(0, len(final_list)):
		list_occur.append(final_list[i][1])
		
	return final_list, list_kmers, list_occur

#Parameters: 
#	top_five: The list of the top  sorted kmers found in the fasta file
#	sequence_list: the list of all the sequences in the fasta file
def sequenceFrequency2(list_of_kmers, sequence_list):


	count = []
	count2 = []
	for x in xrange(0,len(list_of_kmers)):
		acum = 0
		for i in xrange(0, len(sequence_list)):
			if list_of_kmers[x][0] in sequence_list[i]:

				acum += 1
			if list_of_kmers[x][0] not in sequence_list[i]:
				acum +=0
		count.append(acum)	

	for x in xrange(0, len(list_of_kmers)):
		count2.append(list_of_kmers[x][0])

		
	return count, count2

#Parameters: 
#	seqList = the list of all the sequences in the fasta file
#	kmer 	= the full list of kmers found in one particular sequence in the fasta file 
#Using the two parameters we use make calls 'i' calls for how many sequences are in the file
#This will then create a temp list of all the kmers in one sequence
#We take all the kmers in every sequence and join them into one list return the list.
def sequenceCount(seqlist, kmer):

	kmers_extract = []
	full_kmer_list = []
	for i in xrange(0,len(seqlist)):
		kmers_extract = getSubstrings(seqlist[i],kmer) 
		full_kmer_list.append(kmers_extract)

	return full_kmer_list

#Parameters
#	g: the particular one sequence being worked on passed in from the fasta file
#	k: the length of the kmers that we are looking for
def getSubstrings(g, k):

    substrings = list()
    # Start from first character, split into 'k' size chunks
    # Move along one character and repeat. No sense carrying on beyond
    # a starting point of 'k' since that will be the first iteration again.
    for i in range(k):
        line = g[i:]
        substrings += [line[i:i + k]
                       for i in range(0, len(line), k) if i + k <= len(line)]
        
    return substrings # returning the list of every substring kmer in a sequence

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