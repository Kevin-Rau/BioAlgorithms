import os, sys, getopt, math

######################################################################################
#	File 			: Kevin_Rau_Assignment1.py
#	Purpose         : Use command line parameters to control the behavior of an application. 
#				   	: Read a sequence file in FASTA format.
#				   	: Find a specific sequence to be processed within the file.
#				  	: Calculate the GC content across an entire sequence.
#				  	: Search the sequence for exact matches to a given k-mer.
#				  	: Report those matches using an annotation file format (GFF).
#	Developer       : Kevin Rau Jan 2017
#	Python 			: Version 2.7.10
######################################################################################
#   
#   Running the Program:
#   python Kevin_Rau_Assignment1.py -f <filename> -s <seq name> -k <k-mer>
#	The -f is the fasta file to be used
#	The -s is the sequence to be looked at
#	!!! Must be input in the format '>chrIII' or it will not work, string-argument converting issues !!!
#	The -k is the size of kmer (integer between 3 and 8 inclusive) to search.
#
#	***** EXAMPLE RUN *****
#
#	python Kevin_Rau_Assignment1.py -f S288C_R64.fasta -s '>chrI' -k CGGCC
#
######################################################################################
# 
#   References: Alex Okeson Jan 2016
#               alexokeson_hw1.py
#				Formatting and flag handling
#				Stack Overflow for parsing FASTA files and reading/writing these files
#
######################################################################################
#	usage() 
#	This function helps the user if they don't know how to use the program
######################################################################################

def usage():
	print "python Keivn_Rau_Assignment1.py -f <filename> -s <seq name> -k <k-mer>"
	print "The -f specifies the fasta file to be read and analyzed."
	print "The -s is the sequence name to find. Must look something like this: '>chrIII'"
	print "The -k is the size of kmer (integer between 3 and 8 inclusive) to search."

######################################################################################
#	read_fasta
#	This is explained in the seq_find() section following this function
#	Basically, here we are going to read in the file and store in a 2d array
#	Thank you stack overflow for these next two functions
######################################################################################
def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):	# Lines that start with a > indicate a new sequence
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

######################################################################################
#	seq_find
#	This fuction takes two arguments:
#	The first being the FASTA file we are parsing and the sequence we want to work on
#	We are going to parse and store Header Names and Sequence information in a 2D
#	Array through calling the read_fasta funciton
#	Once this is done we compare the squence we want to work on to the saved Headers
#	When there is a match, this is the sequence we want to work and and return it
######################################################################################

def seq_find(fp, seqInsert):
	with open(fp) as fpp:
		for name, seq in read_fasta(fpp):	
			if(seqInsert == name):
				return seq

######################################################################################
#	gff_print()
#	This function prints the output of the file in GFF format
#	The parameters that are passed in are the ones we all calculated in sequence_work()
#	First is the input file we are working on, next is the sequence that was given
#	as a command line arguement, following this is the kmer we have specified
#	seqFound, is the sequence we are actually working on giving the matched Header
#	Lastly, we pass in the calculated GC content.
#	
######################################################################################

def gff_print(inputFile, seqInsert, kmerFind, seqFound, cg):
	print "File:\t", inputFile
	print "Seq:\t", seqInsert
	print "k-mer:\t", kmerFind
	print ""
	print "Seq Length:\t", len(seqFound)
	print "GC Content:\t", round(cg), "%"	
	print ""

######################################################################################
#	reverse_compl()
#	This function takes in our specifed kmer and returns the reverse-compliment of it
#	Huge thanks to Ryan for help on this one
######################################################################################
def reverse_compl(kmerFind):
	reverseComp = kmerFind
	reverse = reverseComp[::-1]
	reverseComp = reverseComp.replace('A','a')	
	reverseComp = reverseComp.replace('T','A')
	reverseComp = reverseComp.replace('a','T')
	reverseComp = reverseComp.replace('C','c')
	reverseComp = reverseComp.replace('G','C')
	reverseComp = reverseComp.replace('c','G')
	return reverseComp

######################################################################################
#	sequence_work()
#
#	Passing in the parameters in order:
#	The sequence specifed, the sequence data, the kmer specified, and the file to used
#	First we are going to calculate the GC content and the length of the sequence
#	This is as simple as running through the array in a for loop and counting 
#	The instances of G's and C's and using len() for arithmatic
#	Pass this in to our print formatter and then seperatly print the kmer information
#	using .find() this made is easy to run through only the array once with the 
#	5-3 sequence and then the reverse compliment of the given kmer.
#
######################################################################################
def sequence_work(seqInsert, seq_data, kmerFind, inputFile):

    		G = float(0.0)
    		C = float(0.0)
    		#k_mer = [G,C,C]
    		length = float(len(seq_data))
    		#where all the work is going to be done
    		for i in xrange(1,len(seq_data)):

    			if seq_data[i] == 'G':
    				G += 1
    			
    			if seq_data[i] == 'C':
    				C += 1

    		#print(seq)

    		GC_content = (G + C)/(length) * float(100.0)
    		#print 'Sequence Length:', len(seq)
    		#print 'GC Content:',int(GC_content),'%'
    		#print 'G:',int(G)
    		#print 'C:',int(C)
    		gff_print(inputFile, seqInsert, kmerFind, seq_data, GC_content)

    		kmerFound = 0
    		kmerLocation = 1
    		reverseComp = reverse_compl(kmerFind)	# finds the reverse compliment of a kmer
    		addDist = len(kmerFind)
    		for i in xrange(1,len(seq_data)):
    			if seq_data.find(kmerFind, kmerFound) != -1:	# -1 means it didn't find anything. Anything positove means it found something
    				kmerFound = (seq_data.find(kmerFind, kmerFound) + 1)
    				print seqInsert, "\t", "Kevin_Rau_Assignment1.py\t", "match\t", kmerFound, "\t", kmerFound + addDist, "\t", "100.\t", "+\t", ".\t", "ID:", kmerLocation
    				kmerLocation += 1
    			if seq_data.find(reverseComp, kmerFound) != -1:
    				kmerFound = (seq_data.find(reverseComp, kmerFound) + 1)
    				print seqInsert, "\t", "Kevin_Rau_Assignment1.py\t", "match\t", kmerFound, "\t", kmerFound + addDist, "\t", "100.\t", "-\t", ".\t", "ID:", kmerLocation
    				kmerLocation += 1
    			else:
    				break
######################################################################################
#	main()
#	Going to grab the arguments and pass them in acoordingly to our functions
#	If there are formatting issues we throw errors
###################################################################################### 

def main(argv):
	inputFile = ''
	kmerFind = ''
	seqInsert = ''
	try:
		opts, args = getopt.getopt(argv, "hf:s:k:")	# This searches for the proper flags in the arguements
	except getopt.GetoptError:
		print "ERROR! Missing and/or incorrect Flags. Correct usage:"
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':	# help option for the user if needed
			print "Help"
			usage()
			sys.exit()
		if opt == '-s':	# User must pass in a valid sequence with single quotes and a > character
			seqInsert = arg
		if opt == '-f':
			inputFile = arg
			dataFile = open(inputFile)
			data = dataFile.read()
		if opt == '-k':
			kmerFind = arg
	seq = seq_find(inputFile, seqInsert)	# The code below handles the flow of the application and printing to the screen
	sequence_work(seqInsert, seq, kmerFind, inputFile)

######################################################################################
#	__name__
#	This function calls the main function and passes in the parameters.
######################################################################################

if __name__ == "__main__":
	main(sys.argv[1:])
