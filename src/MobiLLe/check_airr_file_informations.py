#!/usr/bin/env python3
#-*-coding:Utf-8-*-

#Change the id of the sequences in the AIRR file. Two files must be pass in argument : the tsv file in AIRR format and a output file to store the match between old and new id. 

import sys
import os
from optparse import OptionParser

# ------------------------------------------------------------------------------------

def read_file(file_name):
	input_file = open(file_name,"r")
	lines = input_file.readlines()
	input_file.close()
	return lines

# ------------------------------------------------------------------------------------

def save_txt_format(file_name, data):
	output_file = open(file_name, "w")
	for key, value in data.items():
		output_file.write(key+"\t"+value+"\n")
	output_file.close()

# ------------------------------------------------------------------------------------

def save_AIRR_file(file_name,columns,data) :
	output_file = open(file_name, "w")
	output_file.write(columns)
	for values in data.values():
		output_file.write("\t".join(values)+"\n")
	output_file.close()

# ------------------------------------------------------------------------------------

def create_log_file(filename, all_seq, unnanotated_seq):
	with open(filename+"_log.txt", "w") as f:
		f.write("Total sequence count : " + str(all_seq) + "\n")
		f.write("Number of unannotated sequences that have been eliminated from the analysis : " + str(unnanotated_seq) + "\n")

# ------------------------------------------------------------------------------------

def check_columns(file, request_columns):
	columns_index = {}
	with open(file) as f:
		first_line = f.readline().rstrip()
	header = first_line.split("\t")
	for i in range(len(header)) :
		if header[i] in request_columns :
			columns_index[header[i]] = i
			request_columns.remove(header[i])
	return columns_index	

# ------------------------------------------------------------------------------------

def get_missing_columns(remaining_columns, cluster_columns, tree_columns) :
	for column in remaining_columns :
		if column in cluster_columns :
			return False
		if column in tree_columns :
			return False
	return True

# ------------------------------------------------------------------------------------

def check_columns_tree(sequences, required_columns) :

	all_columns = True
	sequences["whole_seq"] = ""
	sequences["germline_seq"] = ""
	
	for column in required_columns :
		if sequences[column] == "" :
			all_columns = False
			break
			
	if (all_columns) :
		sequences["whole_seq"] = sequences["fwr1"] + sequences["cdr1"] + sequences["fwr2"] + sequences["cdr2"] + sequences["fwr3"] + sequences["cdr3"] + sequences["fwr4"]
		sequences["germline_seq"] = sequences["v_germline_alignment"] + sequences["np1"] + sequences["d_germline_alignment"] + sequences["np2"] + sequences["j_germline_alignment"]
		sequences['germline_seq'] = sequences['germline_seq'].replace('.','')

# ------------------------------------------------------------------------------------

def check_columns_cluster(sequences, required_columns, patterns) :
	for column in required_columns :
		if column in patterns :
			sequences[column] = check_value_for_clustering(sequences[column], patterns[column])
		if sequences[column] == "" :
			return False
	return True

# ------------------------------------------------------------------------------------

def get_region_id(region_column, pattern) :
	for element in region_column :
		if element[:4].upper()==pattern :
			return element.replace(',', '')
	return ""

# ------------------------------------------------------------------------------------

def check_value_for_clustering(region, pattern) :
	if (pattern == "AA") :
		elements = set('ACDEFGHIKLMNPQRSTVYWX*-#.')
	elif (pattern == "nt") :
		elements = set('ACTGN-.')
	elif (pattern == "quality") :
		elements = set('N')
	else :
		return get_region_id(region.split(" "), pattern)
	
	if (all(base.upper() in elements for base in region)) :
		return region
	else :
		return ""
		
# ------------------------------------------------------------------------------------

#modify the ID of all the sequences and create a new file
def get_valid_sequences_ids(AIRR_file,file_name,columns,cluster_columns,tree_columns,patterns):

	lines = read_file(AIRR_file)
	sequence_element = {}
	seq_number = 0
	unnanotated_seq = {}
	ids = {} #contain the match between the old id and the new one

	output_file_name = file_name+"_new_id.tsv" 
	output_file = open(output_file_name, "w")
	output_file.write("\t".join(columns.keys())+"\twhole_seq\tgermline_seq\n")
	for l in range(1,len(lines)):
		elements = lines[l].split("\n")[0]
		properties = elements.split("\t")
		for name in columns :
			if columns[name] < len(properties) :
				sequence_element[name] = properties[columns[name]]
			else :
				sequence_element[name] = ""

		valid_column_cluster = check_columns_cluster(sequence_element,cluster_columns,patterns)

		if valid_column_cluster :	
			check_columns_tree(sequence_element, tree_columns)
			seq_number += 1
			if "sequence_id" in sequence_element :
				ids[sequence_element["sequence_id"]]="seq"+str(seq_number)
			sequence_element["sequence_id"]="seq"+str(seq_number)
			output_file.write("\t".join(sequence_element.values())+"\n")
		else :
			unnanotated_seq[l] = list(sequence_element.values())
		
	output_file.close()
	
	return ids, unnanotated_seq, len(lines)-1	

# ------------------------------------------------------------------------------------

def main():
	usage = "usage: check_airr_file_informations.py -a AIRR -o output_file"
	parser = OptionParser(usage)
	parser.add_option("-a", "--AIRR", dest="AIRR", help="read data from AIRR")
	parser.add_option("-o", "--output_file",dest="output_file", help="write data to output_file containing the id")
	(options, args) = parser.parse_args()

	if(len(sys.argv)==5):
		AIRR_file = options.AIRR
		output_file = options.output_file
		request_columns = ['sequence_id', 'sequence', 'productive', 'v_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'np1', 'np2', 'cdr1', 'cdr2', 'cdr3', 'cdr3_aa', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'v_identity', 'd_identity', 'j_identity', 'v_sequence_start', 'v_sequence_end', 'd_sequence_start', 'd_sequence_end', 'j_sequence_start', 'j_sequence_end', 'cdr1_start', 'cdr1_end', 'cdr2_start', 'cdr2_end', 'cdr3_start', 'cdr3_end', 'fwr1_start', 'fwr1_end', 'fwr2_start', 'fwr2_end', 'fwr3_start', 'fwr3_end', 'fwr4_start', 'fwr4_end', 'v_sequence_alignment', 'd_sequence_alignment', 'j_sequence_alignment', 'v_germline_alignment', 'd_germline_alignment', 'j_germline_alignment']
		cluster_essential_columns = ['v_call', 'j_call','cdr3_aa']
		tree_essential_columns = ['cdr1', 'cdr2', 'cdr3', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'v_identity', 'd_identity', 'j_identity', 'v_sequence_start', 'v_sequence_end', 'd_sequence_start', 'd_sequence_end', 'j_sequence_start', 'j_sequence_end', 'cdr1_start', 'cdr1_end', 'cdr2_start', 'cdr2_end', 'cdr3_start', 'cdr3_end', 'd_sequence_alignment', 'v_germline_alignment', 'd_germline_alignment', 'j_germline_alignment']
		patterns = {"v_call":"IGHV","j_call":"IGHJ", "cdr3_aa":"AA"}
		ids = {}
		unnanotated_seq = {} 
		nb_seq = 0
		
		columns_index = check_columns(AIRR_file,request_columns)
		columns_present = get_missing_columns(request_columns,cluster_essential_columns,tree_essential_columns)
		
		if columns_present :
			ids, unnanotated_seq, nb_seq = get_valid_sequences_ids(AIRR_file,output_file, columns_index,cluster_essential_columns,tree_essential_columns,patterns)
			create_log_file(output_file, nb_seq, len(unnanotated_seq))
			save_txt_format(output_file+"_sequence_id.txt",ids)
		else :
			print("The analysis can't be run because in the AIRR file the following columns are missing : "+",".join(request_columns)+".")
			
		columns = "\t".join(columns_index.keys())+"\twhole_seq\tgermline_seq\n"
		save_AIRR_file(output_file+"_unannotated_seq.tsv",columns,unnanotated_seq)
		
	else:
		parser.error("incorrect number of arguments")
	
	print("  # AIRR file informations checked and ID changed")

# ------------------------------------------------------------------------------------

if __name__ == "__main__":
	main()
