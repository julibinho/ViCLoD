"""
Author: Lucile Jeusset, Laboratory of Computational and Quantitative Biology, UPMC, Paris.
Email: lucile_jeusset@hotmail.fr

Align sub-clones sequences of the same clone and determine the HAUS sequence

"""

import os
import sys
import time
from optparse import OptionParser
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO


def convert_parameter_to_int(nbr_clone, default):
	"""
	Convert a string passed as a parameter to integer
	input nbr_clone:	str	number of clone or sub-clone to analyze, it can be a int or "all"
	input default:		int	value taken by default
	output:		int	number of clone or sub-clone to analyze
	"""	
	if nbr_clone == "all":
		return -1
	else :
		try:
			int(nbr_clone)
		except ValueError:
			return default
		return int(nbr_clone)
		
#-------------------------------------------------------------------#
	
def convert_to_int(n) :
	"""
	Returns the value of the variable if it is an int otherwise 0
	input n:	str	number to test
	output:	int	1 if the given string can't be an integer 
	"""
	try:
		int(n)
	except ValueError:
		return 1
	return int(n)

#-------------------------------------------------------------------#
	
def convert_to_float(n) :
	"""
	Returns the value of the variable if it's a float otherwise 0
	input n:	str	number to test
	output:	int	0 if the given string can't be a float 
	"""
	try:
		float(n)
	except ValueError:
		return 0
	return float(n)

#-------------------------------------------------------------------#
	
def format_region_pos(list_region):
	"""
	Format the start and end values for a region
	input list_region:	list()	positions (start and end) of the different regions (fwr1, cdr1, ..)
	output list_pos:	list()	positions corrected
	"""
	prev_value = 1
	list_pos = []
	for pos in list_region:
		try:
			int(pos)
			list_pos.append(int(pos))
			prev_value = int(pos)
		except ValueError:
			list_pos.append(prev_value)
	return list_pos
	
#-------------------------------------------------------------------#

def get_repertoire_infos(file_name, nbr_clone, nbr_sub_clone):
	"""
	Read the file (file_name) to get information about a certain number of clones (nbr_clone).
	input file_name:		str		file (repertoire_two_levels_info) containing clones and sub-clones information
	input nbr_clone:		int		number of clone to analyze (-1 if all the clone must be analyzed)
	input nbr_clone:		int		number of sub-clone to analyze (-1 if all the sub-clone must be analyzed)
	output repertoire:		Dict()		key=clone id, value=dict of sub-clone information
	output germline_infos:		Dict()		key=clone id, value=list of tuple (sub-clone ID, V identity, J identity, D identity, number of read, germline sequence)
	"""
	repertoire = {}
	germline_infos = {}
	i=1
	with open (file_name,'r') as f:
		for line in f:
			if i > nbr_clone:
				break
			clone_infos = (line.split("\n")[0]).split("\t")
			sub_clones = clone_infos[4].split(" ")
			
			j=1
			clone_id = clone_infos[0].split(" ")[1]
			repertoire[clone_id] = {}
			germline_infos[clone_id]= []
			for s in sub_clones:
				if j > nbr_sub_clone:
					break
				sub_clone_infos = s.split(",")
				sub_clone_id = sub_clone_infos[0].split(":")[1]
				repertoire[clone_id][sub_clone_id] = {}
				repertoire[clone_id][sub_clone_id]["ref_seq"] = sub_clone_infos[1]
				repertoire[clone_id][sub_clone_id]["occurrence"] = convert_to_int(sub_clone_infos[2])
				repertoire[clone_id][sub_clone_id]["productive"] = sub_clone_infos[3]
				repertoire[clone_id][sub_clone_id]["cdr3_aa"] = sub_clone_infos[4]
				repertoire[clone_id][sub_clone_id]["region"] = format_region_pos(sub_clone_infos[5].split("-"))
				repertoire[clone_id][sub_clone_id]["v_identity"] = convert_to_float(sub_clone_infos[6])
				repertoire[clone_id][sub_clone_id]["d_identity"] = convert_to_float(sub_clone_infos[7])
				repertoire[clone_id][sub_clone_id]["j_identity"] = convert_to_float(sub_clone_infos[8])
				repertoire[clone_id][sub_clone_id]["whole_seq"] = sub_clone_infos[9].upper()
				repertoire[clone_id][sub_clone_id]["germline_seq"] = sub_clone_infos[10].upper()
				germline_infos[clone_id].append((sub_clone_id,repertoire[clone_id][sub_clone_id]["v_identity"],repertoire[clone_id][sub_clone_id]["j_identity"],repertoire[clone_id][sub_clone_id]["d_identity"],repertoire[clone_id][sub_clone_id]["occurrence"],repertoire[clone_id][sub_clone_id]["germline_seq"]))
				j+=1
			i += 1
	return repertoire, germline_infos
	
#-------------------------------------------------------------------#

def select_naive_seq(list_seq):
	"""
	Determine the HAUS sequence for a set of sub-clones
	input list_seq:	list()		list of all the sub-clones in a clone (sub-clone ID, V identity, J identity, D identity, number of read, germline sequence)
	output:		tuple		(sub-clone id, germline sequence of this sub-clone)	
	"""
	# sort the sequences based on 1) V identity 2) J identity 3) D identity 4) abundance 
	sorted_list = sorted(list_seq, key=lambda x: (-x[1], -x[2], -x[3], -x[4]))
	i = 0
	germline_seq = sorted_list[i][5]
	while (germline_seq=="" and (i+1)<len(sorted_list)) :
		i += 1
		germline_seq = sorted_list[i][5]
	return (sorted_list[i][0],germline_seq)

#-------------------------------------------------------------------#

def write_clonaltree_align(germline_seq, sub_clones, file_name):
	"""
	Write sequence of all subclones and their HAUS sequence to fasta file
	input germline_seq:	tuple		(sub-clone id, germline sequence of this sub-clone)
	input sub_clones:	Dict()		key=sub-clone id, value=dict containing sub-clone information
	input file_name:	str		name of the output file
	"""
	output_file = file_name+"_selected_seq.fasta"
	filetowrite=open(output_file,"w")

	germline = ">germline"+ "\n" +germline_seq[1]+ "\n"
	filetowrite.write(germline)
	
	for seq in sub_clones:
		filetowrite.write( ">" + seq + "@" + str(sub_clones[seq]["occurrence"]) + "\n" + sub_clones[seq]["whole_seq"] + "\n")

	filetowrite.close()

#-------------------------------------------------------------------#

def alignment(repertoire_name):
	file_name = repertoire_name+"_selected_seq.fasta"
	muscle_cline = MuscleCommandline(input=file_name,out=os.path.splitext(file_name)[0]+".aln")
	muscle_cline()
	align = AlignIO.read(os.path.splitext(file_name)[0]+".aln", "fasta")
	count = SeqIO.write(align, os.path.splitext(file_name)[0]+"_uniq.aln.fa", "fasta")
	aligned_seq = dict([(seq_record.id.split("@")[0],seq_record.seq) for seq_record in SeqIO.parse(os.path.splitext(file_name)[0]+"_uniq.aln.fa","fasta")])
	return aligned_seq

#-------------------------------------------------------------------#

def column_from_residue_number(sequence, regions):
	"""
	AAAA AA AAAA ===> [3,6]
	A-AAAA-AAAA---A
	A-AAA/A-A /AAA---A ===>[4, 8]
	i=0, 1, 2, 3, 4
	j=0, 0, 1, 2, 3
	"""
	#Region variable : "startV","startCDR1","endCDR1","startCDR2","endCDR2","startCDR3","endV","startD","endD","startJ","endCDR3","endJ"
	nb_regions = len(regions)
	pos_without_gap = 0
	name_index = 0
	regions_alignment = []
	for i in range(len(sequence)) :
		if sequence[i] != '-' :
			pos_without_gap += 1
			while((name_index < nb_regions) and (pos_without_gap == regions[name_index])) :
				regions_alignment.append(i+1)
				name_index += 1
	#The position of one (or more) region is greater than the length of the aligned sequence
	while(name_index < nb_regions):
		regions_alignment.append(i+1)
	
	return regions_alignment
	
#-------------------------------------------------------------------#

def get_aligned_region_position(sequences_aligned, sub_clones, germline_id) :
	"""
	Get position of each region from aligned sequences
	input sequences_aligned:	Dict()		key=sub-clone id, value=aligned sequence
	input sub_clones:		Dict()		key=sub-clone id, value=dict containing sub-clone information
	input germline_id:		str		id of the sub-clone to which the germline corresponds
	output alignment_infos:		Dict()		key=sub-clone id, value=dict containing aligned sequence and region positions
	"""
	alignment_infos = {}
	for seq_id in sequences_aligned:
		alignment_infos[seq_id] = {"sequence": str(sequences_aligned[seq_id])}
		if seq_id == "germline":
			alignment_infos[seq_id]["region"] = column_from_residue_number(alignment_infos[seq_id]["sequence"], sub_clones[germline_id]["region"])
		else:
			alignment_infos[seq_id]["region"] = column_from_residue_number(alignment_infos[seq_id]["sequence"], sub_clones[seq_id]["region"])
	return alignment_infos

#-------------------------------------------------------------------#

def write_sub_clone_alignment(alignment_infos, file_name):
	"""
	Writes the aligned sequences and the position of their regions (fwr1, cdr1, ...) to a file
	input alignment_infos:		Dict()		key=sub-clone id, value=dict containing aligned sequence and region positions
	input file_name:		str		name of the output file
	"""
	output_file = file_name + "_aligned_regions.txt"
	with open (output_file, 'w') as f:
		for seq_id in alignment_infos:
			line = seq_id + "\t" + alignment_infos[seq_id]["sequence"] + "\t" + "-".join(map(str, alignment_infos[seq_id]["region"])) + "\n"
			f.write(line)

#-------------------------------------------------------------------#

def main():
	start_time = time.time()
	usage = "usage: alignment_intraclonal.py -r repertoire_two_levels_info -c number_of_clone -s number_of_sub_clone -o output"
	parser = OptionParser(usage)
	parser.add_option("-r", "--repertoire_two_levels_info", dest="repertoire_two_levels_info", help="read data from repertoire_two_levels_info")
	parser.add_option("-c", "--number_of_clone",dest="number_of_clone", help="the number of clone to analyze, if use all, there will be no selection for the analyzing clonotypes otherwise choose the number of n most abondant clonotype ")
	parser.add_option("-s", "--number_of_sub_clone", dest="number_of_sub_clone", help="the number of sub-clone to analyze, if use all, there will be no selection for the analyzing sub-clones otherwise choose the number of n most abondant sub-clones ")
	parser.add_option("-o", "--output", dest="output", help="path and name of the output file")

	(options, args) = parser.parse_args()
	if len(sys.argv) != 9:
		parser.error("incorrect number of arguments")

	repertoire_file = options.repertoire_two_levels_info
	nbr_clone = options.number_of_clone
	nbr_sub_clone = options.number_of_sub_clone
	output_name = options.output
	num = 0
	
	nbr_clone = convert_parameter_to_int(nbr_clone, 5)
	nbr_sub_clone = convert_parameter_to_int(nbr_sub_clone, 200)
	repertoire, germline_infos = get_repertoire_infos(repertoire_file, nbr_clone, nbr_sub_clone)
	
	for clone in repertoire:
		num += 1
		
		if len(repertoire[clone].keys()) > 1:
			file_name = output_name+"_"+str(num)+"_"+str(nbr_sub_clone)
			germline_seq = select_naive_seq(germline_infos[clone])
		
			if germline_seq[1] != "" :
				write_clonaltree_align(germline_seq, repertoire[clone], file_name)
				aligned_seq = alignment(file_name)
				alignment_infos = get_aligned_region_position(aligned_seq, repertoire[clone], germline_seq[0])
				write_sub_clone_alignment(alignment_infos, file_name)

	print("  Sub-clone sequence alignment step execution time : %s seconds " % (time.time() - start_time))
	
#-------------------------------------------------------------------#

if __name__ == "__main__":
	main()
