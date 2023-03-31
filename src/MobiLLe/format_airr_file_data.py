#!/usr/bin/env python3
#-*-coding:Utf-8-*-

"""
"""

import sys
from optparse import OptionParser
import time

column_headers = ['sequence_id', 'sequence', 'productive', 'v_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'np1', 'np2', 'cdr1', 'cdr2', 'cdr3', 'cdr3_aa', 'fwr1', 'fwr2', 'fwr3', 'fwr4', 'v_identity', 'd_identity', 'j_identity', 'v_sequence_start', 'v_sequence_end', 'd_sequence_start', 'd_sequence_end', 'j_sequence_start', 'j_sequence_end', 'cdr1_start', 'cdr1_end', 'cdr2_start', 'cdr2_end', 'cdr3_start', 'cdr3_end', 'fwr1_start', 'fwr1_end', 'fwr2_start', 'fwr2_end', 'fwr3_start', 'fwr3_end', 'fwr4_start', 'fwr4_end', 'v_sequence_alignment', 'd_sequence_alignment', 'j_sequence_alignment', 'v_germline_alignment', 'd_germline_alignment', 'j_germline_alignment','v_sequence_alignment_aa', 'd_sequence_alignment_aa', 'j_sequence_alignment_aa', 'whole_seq', 'germline_seq', 'occurrence']
cluster_annotations = ['sequence_id','sequence', 'productive', 'v_call', 'j_call', 'cdr3_aa', 'v_sequence_alignment', 'j_sequence_alignment', 'v_sequence_alignment_aa', 'j_sequence_alignment_aa']
tree_annotations = ['v_identity', 'd_identity', 'j_identity', 'v_sequence_start', 'v_sequence_end', 'd_sequence_start', 'd_sequence_end', 'j_sequence_start', 'j_sequence_end', 'cdr1_start', 'cdr1_end', 'cdr2_start', 'cdr2_end', 'cdr3_start', 'cdr3_end', 'v_germline_alignment', 'd_germline_alignment', 'j_germline_alignment','germline_seq']
# global variable for output file names:
output_files = {"matching_seq_file": "_AIRR_seq_name_matching.txt", "unique_seq_name_file": "_unique_seq_id_matching.txt", "unannotated_seq_file": "_unannotated_seq.tsv", "annotated_seq_file": "_airr_new_id.tsv"}

# ------------------------------------------------------------------------------------

def check_header(line, columns) :
	"""
	Find the required columns in the AIRR file
	input line:	str	first line of the AIRR file
	input columns:	Dict()	key=column name; values=index of this column in the AIRR file
	output:	str	name of the part of the processing that can't be performed due to missing columns
	output missing_columns:	list()	list containing the name of the missing columns
	"""
	cluster_header = ['sequence', 'productive', 'v_call', 'j_call', 'cdr3_aa', 'v_sequence_alignment', 'j_sequence_alignment', 'v_sequence_alignment_aa', 'j_sequence_alignment_aa']
	tree_header = ['np1', 'np2', 'v_identity', 'd_identity', 'j_identity', 'v_sequence_start', 'v_sequence_end', 'd_sequence_start', 'd_sequence_end', 'j_sequence_start', 'j_sequence_end', 'cdr1_start', 'cdr1_end', 'cdr2_start', 'cdr2_end', 'cdr3_start', 'cdr3_end', 'v_germline_alignment', 'd_germline_alignment', 'j_germline_alignment']

	header = (line.split('\n')[0]).split('\t')
	for i in range(len(header)) :
		if header[i] in columns :
			columns[header[i]] = i
	missing_columns = [key for key, value in columns.items() if value == -1] 
	for element in missing_columns :
		if element in cluster_header :
			return 'cluster', missing_columns
		elif element in tree_header :
			return 'tree', missing_columns
	return 'none', missing_columns

# ------------------------------------------------------------------------------------

def get_sequence_informations(line, columns, seq_name) :
	"""
	Get the information needed for a line in the AIRR file
	input line:	str	line of the AIRR file
	input columns:	Dict()	key=column name; values=index of this column in the AIRR file
	input seq_name:	str	name of the sequence
	output sequence_info:	Dict()	contains all the information of a sequence, key=column name, value=information in the AIRR file for this column
	"""
	infos_type = {'sequence_id':'str', 'sequence':'nt', 'productive':'bool', 'v_call':'IGHV', 'j_call':'IGHJ', 'sequence_alignment':'nt', 'germline_alignment':'nt', 'junction':'nt', 'np1':'nt', 'np2':'nt', 'cdr1':'nt', 'cdr2':'nt', 'cdr3':'nt', 'cdr3_aa':'AA', 'fwr1':'nt', 'fwr2':'nt', 'fwr3':'nt', 'fwr4':'nt', 'v_identity':'float', 'd_identity':'float', 'j_identity':'float', 'v_sequence_start':'int', 'v_sequence_end':'int', 'd_sequence_start':'int', 'd_sequence_end':'int', 'j_sequence_start':'int', 'j_sequence_end':'int', 'cdr1_start':'int', 'cdr1_end':'int', 'cdr2_start':'int', 'cdr2_end':'int', 'cdr3_start':'int', 'cdr3_end':'int', 'fwr1_start':'int', 'fwr1_end':'int', 'fwr2_start':'int', 'fwr2_end':'int', 'fwr3_start':'int', 'fwr3_end':'int', 'fwr4_start':'int', 'fwr4_end':'int', 'v_sequence_alignment':'nt', 'd_sequence_alignment':'nt', 'j_sequence_alignment':'nt', 'v_germline_alignment':'nt', 'd_germline_alignment':'nt', 'j_germline_alignment':'nt','v_sequence_alignment_aa':'AA', 'd_sequence_alignment_aa':'AA', 'j_sequence_alignment_aa':'AA', 'sequence_old_ids':'str'}
	sequence_info = {}
	line_info = (line.split("\n")[0]).split("\t")
	num_column = len(line_info)
	for element in columns :
		if(columns[element] == -1):
			sequence_info[element] = ""
		elif(columns[element] > num_column):
			sequence_info[element] = ""
		else:
			sequence_info[element] = check_type_of_info(line_info[columns[element]], infos_type[element])
	sequence_info['sequence_id']=seq_name
	sequence_info['whole_seq']=sequence_info['fwr1']+sequence_info['cdr1']+sequence_info['fwr2']+sequence_info['cdr2']+sequence_info['fwr3']+sequence_info['cdr3']+sequence_info['fwr4']
	sequence_info['germline_seq']=sequence_info['v_germline_alignment']+sequence_info['np1']+sequence_info['d_germline_alignment']+sequence_info['np2']+sequence_info['j_germline_alignment']
	sequence_info['germline_seq']=sequence_info['germline_seq'].replace('.','')
	return sequence_info

# ------------------------------------------------------------------------------------

def save_name_matching (seq_info):
	"""
	Writes in a file the correspondence between the new and the old name of the sequence and add the occurrence of this sequence in the dict
	input seq_info:	Dict()	contains all the information of a sequence, key=column name, value=information in the AIRR file for this column
	"""
	seq_info['occurrence'] = 1
	list_old_ids = (seq_info['sequence_old_ids'][:-1]).split(";")
	if len(list_old_ids) > 0 :
		seq_info['occurrence'] = len(list_old_ids)
	
	filename = output_files["matching_seq_file"]
	with open(filename, 'a') as f:
		line = seq_info['sequence_id']+"\t"+seq_info.pop('sequence_old_ids')+"\n"
		f.write(line)

# ------------------------------------------------------------------------------------

def is_integer(text) :
	"""
	Returns the value of the variable if it is an int otherwise an empty string
	input text:	str	
	output text:	str
	"""
	try:
		int(text)
	except ValueError:
		text = ""
	return text

# ------------------------------------------------------------------------------------

def is_float(text) :
	"""
	Returns the value of the variable if it is a float otherwise an empty string
	input text:	str	
	output text:	str
	"""
	try:
		float(text)
	except ValueError:
		text = ""
	return text

# ------------------------------------------------------------------------------------

def format_boolean_values(value) :
	"""
	Formats boolean values : T for true and F for false
	input value:	bool	boolean value
	output:	str	T or F depending on the value
	"""
	value = value.strip()
	value = value.lower()
	if(value == "t" or value == "true"):
		return "T"
	elif(value == "f" or value == "false"):
		return "F"
	else :
		return ""

# ------------------------------------------------------------------------------------

def get_region_id(region_split, pattern) :
	"""
	Returns the first V or J gene annotation without the species annotation
	input region_split:	list()	annotation of the gene contained in the AIR file (separation of this content in list by cutting at the level of the spaces)
	input pattern:		str	gene annotation type (IGHV, IFHJ)	
	output:		str	V or J gene annotation
	"""
	for element in region_split :
		if element[:4].upper()==pattern :
			return element.replace(',', '')
	return ""

# ------------------------------------------------------------------------------------

def check_type_of_info(info, pattern) :
	"""
	Checks that the sequence information matches the right type (AA,nt,int,float,str,IGH)
	input info:	str	information of a sequence for specific columns (v_call, j_call, cdr3_aa, ...)
	input pattern:	str	type of data expected (IGHV, IGHJ, AA, nt, ...)
	output:	str	information or empty string if it's not in the right type		
	"""
	if (pattern == "AA") :
		elements = set('ACDEFGHIKLMNPQRSTVYWX*-#.')
		info = info.upper()
	elif (pattern == "nt") :
		elements = set('ACTGN-.')
		info = info.lower()
	elif (pattern == "int") :
		return is_integer(info)
	elif (pattern == "float") :
		return is_float(info)
	elif (pattern == "str") :
		return info
	elif (pattern == "bool") :
		return format_boolean_values(info)
	else :
		return get_region_id(info.split(" "), pattern)
	
	if (all(base.upper() in elements for base in info)) :
		return info
	else :
		return ""
		
# ------------------------------------------------------------------------------------

def check_contains_all_infos(seq_infos, quality_filter) :
	"""
	Check all the required informations are provide for a given sequence
	input seq_infos:	Dict()	contains all the information of a sequence, key=column name, value=information in the AIRR file for this column
	input quality_filter:	str	if False, sequences contaning N are kept for the analysis
	output data_filter:	Dict()	key=filter, value=boolean False if the filter conditions are not met
	"""
	data_filter = {'cluster annotation':True, 'tree annotation':True, 'quality':True, 'minimum reads': True, 'occurrence': seq_infos['occurrence']}
	
	for h1 in cluster_annotations :
		if seq_infos[h1] == "":
			data_filter['cluster annotation'] = False
			for h2 in tree_annotations :
				if seq_infos[h2] == "":
					data_filter['tree annotation'] = False
					break
			break
	if quality_filter != "False" :
		if (seq_infos['whole_seq'] == "") or ('N' in seq_infos['whole_seq'].upper()) :
			data_filter['quality'] = False
	
	return data_filter
			
# ------------------------------------------------------------------------------------

def get_unique_sequences(sequence_info, unique_sequences) :
	"""
	Regroup the sequence ID that have the same sequence in nt
	input seq_infos:	Dict()	contains all the information of a sequence, key=column name, value=information in the AIRR file for this column
	input unique_sequences:	Dict()	key=whole sequence, value=sequence id
	"""
	if sequence_info['whole_seq'] in unique_sequences :
		unique_sequences[sequence_info['whole_seq']].append(sequence_info['sequence_id'])
	else :
		unique_sequences[sequence_info['whole_seq']] = [sequence_info['sequence_id']]

# ------------------------------------------------------------------------------------

def missing_information_message(filters) :
	"""
	Message indicating the type of missing information depending on the filters information
	input filters:	Dict()	key=filter, value=boolean False if the filter conditions are not met
	output:	str	message indicating which data are missing
	"""
	if (not filters['cluster annotation'] and not filters['quality']) :
		return "Missing annotation and poor quality sequence"
	elif not filters['cluster annotation'] :
		return "Missing annotation"
	elif not filters['quality'] :
		return "Poor quality sequence"
	else :
		return ""
			
# ------------------------------------------------------------------------------------

def save_sequences(sequence_info, filename, **kwargs) :
	"""
	Save sequences and their annotations to a file
	input seq_infos:	Dict()	contains all the information of a sequence, key=column name, value=information in the AIRR file for this column
	input analysis_name:	str	name of the file
	"""
	line = ""
	with open(filename, 'a') as f:
		for element in column_headers :
			line += str(sequence_info[element])+"\t"
		if 'criteria' in kwargs.keys():
			line += kwargs['criteria']+"\t"
		f.write(line[:-1]+"\n")

# ------------------------------------------------------------------------------------
		
def get_minimum_read(threshold_type, threshold, total_reads) :
	"""
	Get the minimum number of reads requested by the user from the total number of reads
	input threshold_type:	str	indicates with which measurement the user gives the threshold (% or number)
	input threshold:	str	the reads occurrence below this threshold must be deleted
	input total_reads:	int	number of reads in all the AIRR file
	output:	int	occurrence number of read from which a read is taken into account
	"""
	if threshold_type == "percentage" :
		if is_float(threshold) != "" and float(threshold)<2 :
			return int(round((float(threshold)*total_reads)/100))
		else :
			return 1
	elif threshold_type == "number" :
		if is_integer(threshold) != "" and int(threshold)<(total_reads/2) :
			return int(threshold)
		else :
			return 1
	else :
		threshold_type = "percentage"
		threshold = 0.005
		return int(round(threshold*total_reads)/100)

# ------------------------------------------------------------------------------------
		
def update_minimum_read_stats(stats, min_read) :
	"""
	Update stats information regarding the minimum occurrence of reads needed
	input stats:	Dict()	key=sequence id, value=statistical information on this sequence
	input min_read:	int	occurrence number of read from which a read is taken into account
	"""
	for seq_id in stats :
		if stats[seq_id]['occurrence'] < min_read :
			stats[seq_id]['minimum reads'] = False

# ------------------------------------------------------------------------------------

def save_unique_seq_ids(ref_seq, list_seq_ids) :
	"""
	Save match between the reference id and all the ids assigned to a same sequence
	input ref_seq:	str	reference id of the sequence
	input list_seq_ids:	list()	all the other id given to this sequence
	"""
	line = ref_seq+"\t"+";".join(list_seq_ids)+"\n"
	filename = output_files["unique_seq_name_file"]
	with open(filename, 'a') as f:
		f.write(line)
	
# ------------------------------------------------------------------------------------

def update_sequence_occurrence(ref_seq, list_seq_ids, annotated_seq) :
	"""
	Changes the number of occurrence of reads depending on the grouping of unique sequences
	input ref_seq:	str	reference id of the sequence
	input list_seq_ids:	list()	all the other id given to this sequence
	input annotated_seq:	Dict()	key=sequence id, value=Dict containing the different annotations of the sequence
	"""
	occurrence = 0
	for seq_id in list_seq_ids :
		if seq_id in annotated_seq :
			occurrence += annotated_seq[seq_id]['occurrence']
	annotated_seq[ref_seq]['occurrence'] += occurrence

# ------------------------------------------------------------------------------------

def check_greater_than_minimum_read(seq_infos, min_read) :
	"""
	Check the tree annotation and that the occurrence of reads is greater than the set minimum
	input seq_infos:	Dict()	contains all the information of a sequence, key=column name, value=information in the AIRR file for this column
	input min_read:	int	occurrence number of read from which a read is taken into account
	output data_filter:	Dict()	key=filter, value=boolean False if the filter conditions are not met
	"""
	data_filter = {'cluster annotation':True, 'tree annotation':True, 'quality':True, 'minimum reads': True, 'occurrence': seq_infos['occurrence']}
	
	for h in tree_annotations :
		if seq_infos[h] == "":
			data_filter['tree annotation'] = False
			break
	
	if seq_infos['occurrence'] <= min_read :
		data_filter['minimum reads'] = False
	
	return data_filter

# ------------------------------------------------------------------------------------

def save_stats_data(stats, analysis_name, total_seqs, min_read, threshold, threshold_type) :
	"""
	Save the statistical information on the sequences contained in the AIRR file
	input stats:	Dict()	key=sequence id, value=statistical information on this sequence
	input analysis_name:	str	name of the analysis used to name the file
	input total_seq:	int	total number of sequences
	input min_read:	int	occurrence number of read from which a read is taken into account
	input threshold:	str	minimum number of reads taken into account
	input threshold_type:	str	measure to which the minimum number of reads is given
	"""
	filtered_seq = 0
	only_unannotated_tree = 0
	all_stats = {'cluster annotation': 0, 'tree annotation': 0, 'quality': 0, 'minimum reads': 0}
	multiple_criteria = [0,0,0,0,0]
	for seq_id in stats :
		accumulation = 0
		for element in all_stats :
			if not stats[seq_id][element] :
				all_stats[element] += stats[seq_id]['occurrence']
				accumulation += 1
		multiple_criteria[accumulation] += stats[seq_id]['occurrence']
		
		if accumulation == 1 and not stats[seq_id]['tree annotation'] :
			only_unannotated_tree += stats[seq_id]['occurrence']
		else :
			filtered_seq += stats[seq_id]['occurrence']
	
	read_analysed = total_seqs - filtered_seq
	minimum_read = str(min_read)+" read(s): "
	
	line = "Total sequences: " + str(total_seqs) + "\n"
	line += "Threshold for eliminating infrequent reads: " + str(threshold) + "\t" + str(threshold_type) + "\n"
	line += "Sequence to be analysed: " + str(read_analysed) + "\t" + str(round((read_analysed*100)/total_seqs,2)) + "% \n"
	line += "Sequence filtered: " + str(filtered_seq) + "\t" + str(round((filtered_seq*100)/total_seqs,2)) + "% \n"
	line += "Sequence unannotated for clustering: " + str(all_stats['cluster annotation']) + "\t" + str(round((all_stats['cluster annotation']*100)/total_seqs,2)) + "% \n"
	line += "Sequence unannotated for tree analysis: " + str(all_stats['tree annotation']) + "\t" + str(round((all_stats['tree annotation']*100)/total_seqs,2)) + "% \n"
	line += "Low quality sequence: " + str(all_stats['quality']) + "\t" + str(round((all_stats['quality']*100)/total_seqs,2)) + "% \n"
	line += "Sequence with an abundance below " + minimum_read + str(all_stats['minimum reads']) + "\t" + str(round((all_stats['minimum reads']*100)/total_seqs,2)) + "% \n"
	line += "Sequence that combines several criteria: One;" + str(multiple_criteria[1]) + "\tTwo;" + str(multiple_criteria[2]) + "\tThree;" + str(multiple_criteria[3]) + "\tFour;" + str(multiple_criteria[4]) + "\n"
	line += "Sequence only unannotated for tree analysis: " + str(only_unannotated_tree) + "\t" + str(round((only_unannotated_tree*100)/total_seqs,2)) + "% \n"
	
	filename = analysis_name+"_log.txt"
	with open(filename, 'w') as f:
		f.write(line)

# ------------------------------------------------------------------------------------

def create_output_files (analysis_name) :
	"""
	Create and name the 4 output files (XXX_seq_name_matching.txt, XXX_unique_seq_id_matching.txt, XXX_unannotated_seq.tsv and XXX_airr_new_id.tsv) and add a header if the type of files need one.
	input analysis_name:	str	name of the analysis
	"""
	global output_files
	headers = { "unannotated_seq_file": "\t".join(column_headers)+"\tcriteria", "annotated_seq_file": "\t".join(column_headers) }
		
	for filename in output_files :
		file_headers = ""
		if filename in headers :
			file_headers = headers[filename]+"\n"
		output_files[filename] = analysis_name + output_files[filename]
		with open(output_files[filename], 'w') as f :
			f.write(file_headers)
	
# ------------------------------------------------------------------------------------

def format_sequences_information(airr_file_name, columns, analysis_name, threshold, threshold_type, quality_filter) :
	"""
	Browse the AIRR file to retrieve sequences matching certain criteria (annotations, quality, number of reads) and group them by unique sequences
	input airr_file_name:	str	path of the AIRR file
	input columns:	Dict()	key=column name; values=index of this column in the AIRR file
	input analysis_name:	str	path of the output files
	input threshold:	str	minimum number of reads taken into account
	input threshold_type:	str	measure to which the minimum number of reads is given
	input quality_filter:	str	if False, sequences contaning N are kept for the analysis
	"""
	
	annotated_seq = {}
	unique_sequences = {}
	stats = {}
	total_seq = 0
	num_line = 1
	
	airr_file = open(airr_file_name, 'r')
	first_line = airr_file.readline()
	
	for line in airr_file :
		seq_name = "seq"+str(num_line)
		num_line += 1
		sequence_info = get_sequence_informations(line, columns, seq_name)
		save_name_matching (sequence_info)
		filters = check_contains_all_infos(sequence_info, quality_filter)
		if (filters['cluster annotation'] and filters['quality']) :
			get_unique_sequences(sequence_info, unique_sequences)
			annotated_seq[sequence_info['sequence_id']] = sequence_info
		else :
			stats[sequence_info['sequence_id']] = filters
			missing = missing_information_message(filters)
			save_sequences(sequence_info, output_files["unannotated_seq_file"], criteria=missing)
		total_seq += sequence_info['occurrence']

	airr_file.close()

	minimum_read = get_minimum_read(threshold_type, threshold, total_seq)
	update_minimum_read_stats(stats, minimum_read)
	
	for sequence in unique_sequences :
		ref_seq = unique_sequences[sequence][0]
		if len(unique_sequences[sequence])>1 :
			save_unique_seq_ids(ref_seq, unique_sequences[sequence][1:])
			update_sequence_occurrence(ref_seq, unique_sequences[sequence][1:], annotated_seq)
		
		sequence_info = annotated_seq[ref_seq]

		filters = check_greater_than_minimum_read(sequence_info, minimum_read)
		if filters['minimum reads'] :
			save_sequences(sequence_info, output_files["annotated_seq_file"])
			if not filters['tree annotation'] :
				stats[sequence_info['sequence_id']] = filters
		else :
			save_sequences(sequence_info, output_files["unannotated_seq_file"], criteria="Number of reads below the threshold")
			stats[sequence_info['sequence_id']] = filters
	save_stats_data(stats, analysis_name, total_seq, minimum_read, threshold, threshold_type)
	
# ------------------------------------------------------------------------------------

def main():

	start_time = time.time()

	usage = "usage: format_airr_file_data.py -a AIRR -o output_file -m min_read -t type_min -q quality_filter"
	parser = OptionParser(usage)
	parser.add_option("-a", "--AIRR", dest="AIRR", help="read data from AIRR")
	parser.add_option("-o", "--output_file", dest="output_file", help="path of the output files")
	parser.add_option("-m", "--min_read", dest="min_read", help="minimum number of reads taken into account")
	parser.add_option("-t", "--type_min", dest="type_min", help="measure to which the minimum number of readings is given")
	parser.add_option("-q", "--quality_filter", dest="quality_filter", help="If True, sequences contaning N are discarded from the analysis")
	(options, args) = parser.parse_args()

	if(len(sys.argv)==11):
		airr_file_name = options.AIRR
		analysis_name = options.output_file
		threshold = options.min_read
		threshold_type = options.type_min
		quality_filter = options.quality_filter

		columns = {'sequence_id':-1, 'sequence':-1, 'productive':-1, 'v_call':-1, 'j_call':-1, 'sequence_alignment':-1, 'germline_alignment':-1, 'junction':-1, 'np1':-1, 'np2':-1, 'cdr1':-1, 'cdr2':-1, 'cdr3':-1, 'cdr3_aa':-1, 'fwr1':-1, 'fwr2':-1, 'fwr3':-1, 'fwr4':-1, 'v_identity':-1, 'd_identity':-1, 'j_identity':-1, 'v_sequence_start':-1, 'v_sequence_end':-1, 'd_sequence_start':-1, 'd_sequence_end':-1, 'j_sequence_start':-1, 'j_sequence_end':-1, 'cdr1_start':-1, 'cdr1_end':-1, 'cdr2_start':-1, 'cdr2_end':-1, 'cdr3_start':-1, 'cdr3_end':-1, 'fwr1_start':-1, 'fwr1_end':-1, 'fwr2_start':-1, 'fwr2_end':-1, 'fwr3_start':-1, 'fwr3_end':-1, 'fwr4_start':-1, 'fwr4_end':-1, 'v_sequence_alignment':-1, 'd_sequence_alignment':-1, 'j_sequence_alignment':-1, 'v_germline_alignment':-1, 'd_germline_alignment':-1, 'j_germline_alignment':-1,'v_sequence_alignment_aa':-1, 'd_sequence_alignment_aa':-1, 'j_sequence_alignment_aa':-1, 'sequence_old_ids':-1}

		airr_file = open(airr_file_name, 'r')
		first_line = airr_file.readline()
		airr_file.close()

		missing_info, missing_columns = check_header(first_line,columns)
		if missing_info == 'cluster' :
			if len(missing_info) > 1 :
				print("Your data can't be analyzed because these columns are missing : "+", ".join(missing_columns)+".")
			else :
				print("Your data can't be analyzed because this column is missing : "+missing_columns[0]+".")
		else :
			create_output_files(analysis_name)
			format_sequences_information(airr_file_name, columns, analysis_name, threshold, threshold_type, quality_filter)
	else:
		parser.error("incorrect number of arguments")
	
	print("  AIRR file information formatted. Execution time %s seconds " % (time.time() - start_time))		
	
# ------------------------------------------------------------------------------------

if __name__ == "__main__":
	main()
		
		
	
