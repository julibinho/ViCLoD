import os
import sys
from collections import Counter
from optparse import OptionParser
import time

#####################################################################

def get_header_index(first_line):
	"""
	Find the required columns in the AIRR file
	input first_line:	str	first line of the AIRR file
	output columns:	Dict()	key=column name; values=index of this column in the AIRR file
	"""
	columns = {'sequence_id':-1, 'productive':-1, 'v_call':-1, 'j_call':-1, 'cdr3_aa':-1, 'v_sequence_alignment':-1, 'j_sequence_alignment':-1, 'v_sequence_alignment_aa':-1, 'j_sequence_alignment_aa':-1, 'occurrence':-1}
	header = (first_line.split('\n')[0]).split('\t')
	for i in range(len(header)) :
		if header[i] in columns :
			columns[header[i]] = i
	return columns

#####################################################################

def is_integer(n) :
	"""
	Returns the value of the variable if it is an int otherwise an empty string
	input n:	str	number to test
	output:	int	1 if the given string can't be an integer 
	"""
	try:
		int(n)
	except ValueError:
		return 1
	return int(n)

#####################################################################
def dico_V_J_CDR3_format(AIRR, alleles):

	Dico={}

	airr_file = open(AIRR, 'r')
	
	first_line = airr_file.readline()
	headers = get_header_index(first_line)
	separator = ' '
	
	missing = [key for key, value in headers.items() if value == -1] 
	
	if len(missing) == 0 :
	
		if alleles == "0":
			separator = '*'
		gene = ['v_call', 'j_call']
		alignment = ['j_sequence_alignment', 'v_sequence_alignment', 'j_sequence_alignment_aa', 'v_sequence_alignment_aa']
	
		for line in airr_file :
		
			sequence = (line.split("\n")[0]).split("\t")
			seq_id = sequence[headers['sequence_id']]
			Dico[seq_id] = []
			
			Dico[seq_id].append(sequence[headers['productive']])
			for name in gene :
				Dico[seq_id].append(sequence[headers[name]].split(separator)[0])

			Dico[seq_id].append((sequence[headers['cdr3_aa']]).replace("#", "."))
			
			for name in alignment :
				Dico[seq_id].append((sequence[headers[name]]).replace('.', ''))
			
			Dico[seq_id].append(is_integer(sequence[headers['occurrence']]))		

	airr_file.close()

	return Dico

#####################################################################				
def write_file(Dico_VJCDR3,output_file):

	#the output file contains the following information for each sequence 
	#sequence_id	productive	v_call	j_call	cdr3_aa	j_sequence_alignment	v_sequence_alignment	j_sequence_alignment_aa	v_sequence_alignment_aa	occurence
	#the gaps in all sequences have been removed 

	outputname = output_file+"_seq_Fo_V_CDR3_Jseq.txt"
	f = open(outputname,"w")
	for key in Dico_VJCDR3.keys():
		line = key + "\t" + Dico_VJCDR3[key][0] + "\t" + Dico_VJCDR3[key][1] + "\t" + Dico_VJCDR3[key][2] + "\t" + Dico_VJCDR3[key][3] + "\t" + Dico_VJCDR3[key][4]+"\t" + Dico_VJCDR3[key][5]+"\t"+ Dico_VJCDR3[key][6]+"\t" + Dico_VJCDR3[key][7]+"\t" + str(Dico_VJCDR3[key][8])+"\n"
		f.write(line)
	f.close()
	return 0


####################################################################
def main():
	start_time = time.time()
	usage = "usage: format_labeling_imgt_airr.py -a AIRR -o output_file -g gene_alleles"
	parser = OptionParser(usage)
	parser.add_option("-a", "--AIRR", dest="AIRR",
	      help="read data from AIRR")	      
	parser.add_option("-o", "--output_file",dest="output_file",
	      help="write data to output_file")
	parser.add_option("-g", "--alleles_taken_into_account",dest="gene_alleles",
	      help="0 if V and J gene alleles are not taken into account for the clustering")
	(options, args) = parser.parse_args()
	
	if len(sys.argv) != 7:
		parser.error("incorrect number of arguments")
	
	AIRR_file = options.AIRR
	output_file = options.output_file
	gene_alleles = options.gene_alleles

	Dico_VJCDR3 = dico_V_J_CDR3_format(AIRR_file, gene_alleles)
	if Dico_VJCDR3 != {} :
		write_file(Dico_VJCDR3,output_file)
	print("  The AIRR file reading step execution time : %s seconds " % (time.time() - start_time))

#####################################################################

if __name__ == "__main__":
	main()
