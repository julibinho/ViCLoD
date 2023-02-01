"""
Author: Lucile Jeusset, Laboratory of Computational and Quantitative Biology, UPMC, Paris.
Email: lucile_jeusset@hotmail.fr

Group sub-clones of clonal lineage and retieve their caracteristics (sequence_id, sequence, sequence_aa, ...). 

"""

import sys
import time
from optparse import OptionParser

columns_header = ['productive', 'v_call', 'j_call', 'cdr3_aa', 'v_identity', 'd_identity', 'j_identity', 'v_sequence_start', 'v_sequence_end', 'd_sequence_start', 'd_sequence_end', 'j_sequence_start', 'j_sequence_end', 'cdr1_start', 'cdr1_end', 'cdr2_start', 'cdr2_end', 'cdr3_start', 'cdr3_end', 'fwr1_start', 'fwr1_end', 'fwr2_start', 'fwr2_end', 'fwr3_start', 'fwr3_end', 'fwr4_start', 'fwr4_end', 'whole_seq', 'germline_seq', 'occurrence']

def get_header_index(first_line):
	"""
	Find the required columns in the AIRR file
	input first_line:	str		first line of the AIRR file
	output columns:	Dict()		key=column name; values=index of this column in the AIRR file
	"""
	columns = {'sequence_id':-1, 'productive':-1, 'v_call':-1, 'j_call':-1, 'cdr3_aa':-1, 'v_identity':-1, 'd_identity':-1, 'j_identity':-1, 'v_sequence_start':-1, 'v_sequence_end':-1, 'd_sequence_start':-1, 'd_sequence_end':-1, 'j_sequence_start':-1, 'j_sequence_end':-1, 'cdr1_start':-1, 'cdr1_end':-1, 'cdr2_start':-1, 'cdr2_end':-1, 'cdr3_start':-1, 'cdr3_end':-1, 'fwr1_start':-1, 'fwr1_end':-1, 'fwr2_start':-1, 'fwr2_end':-1, 'fwr3_start':-1, 'fwr3_end':-1, 'fwr4_start':-1, 'fwr4_end':-1, 'whole_seq':-1, 'germline_seq':-1, 'occurrence':-1}
	header = (first_line.split('\n')[0]).split('\t')
	for i in range(len(header)) :
		if header[i] in columns :
			columns[header[i]] = i
	return columns

#-------------------------------------------------------------------#

def read_AIRR_file(AIRR, alleles):

	sequences = {}

	airr_file = open(AIRR, 'r')
	
	first_line = airr_file.readline()
	headers = get_header_index(first_line)
	separator = ' '
	
	missing = [key for key, value in headers.items() if value == -1] 
	
		
	if len(missing) == 0 :
	
		if alleles == "0":
			separator = '*'
	
		for line in airr_file :
		
			annotations = (line.split("\n")[0]).split("\t")
			seq_id = annotations[headers['sequence_id']]
			sequences[seq_id] = []
			
			for name in columns_header :
				if name == 'v_call' or name == 'j_call':
					sequences[seq_id].append(annotations[headers[name]].split(separator)[0])
				else:
					sequences[seq_id].append(annotations[headers[name]])
			
	airr_file.close()

	return sequences

#-------------------------------------------------------------------#

def get_cluster_seq(filename):

	clusters = {}
	
	with open (filename, 'r') as f:
		for line in f:
			info = line.split("\t")
			cluster_id = info[0].rstrip()
			cluster_seq = info[1].rstrip().split(" ")
			if cluster_id in clusters :
				clusters[cluster_id].append(cluster_seq)
			else :
				clusters[cluster_id] = cluster_seq

	return clusters

#-------------------------------------------------------------------#
	
def convert_to_int(n) :
	"""
	Returns the value of the variable if it is an int otherwise 1
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

def get_regions_position(regions):
	"""
	Adjust the values ​​of the beginning and the end of the regions according to the beginning of the V region
	input regions:		Dict()		key=region name, value=position
	output position:	Dict()		key=region name, value=position corrected according to the beginning of V region
	"""
	
	regions_name = ['cdr1_start', 'cdr1_end', 'cdr2_start', 'cdr2_end', 'cdr3_start', 'v_sequence_end', 'd_sequence_start', 'd_sequence_end', 'j_sequence_start', 'cdr3_end', 'j_sequence_end']
	
	position = {}
	position['v_sequence_start'] = 1
	for name in regions_name :
		position[name] = (regions[name] - regions['v_sequence_start'])+1
	
	previous_pos = 'v_sequence_start'
	for pos in position :
		if ( (position[pos] < position[previous_pos]) or (position[pos] < 1) ) :
			position[pos] = position[previous_pos]
		previous_pos = pos

	return position

#-------------------------------------------------------------------#

def format_data(sequence_data):
	"""
	Check the sequence annotations to put them in the right format
	input sequence_data:	list()		list containing the annotation of a sequence
	output sequence_data:	list()		corrected sequence annotations
	"""
	regions = {'v_sequence_start':1, 'cdr1_start':1, 'cdr1_end':1, 'cdr2_start':1, 'cdr2_end':1, 'cdr3_start':1, 'v_sequence_end':1, 'd_sequence_start':1, 'd_sequence_end':1, 'j_sequence_start':1, 'cdr3_end':1, 'j_sequence_end':1}
	identity = ['v_identity', 'd_identity', 'j_identity']
	
	for header in regions:
		column_index = columns_header.index(header)
		regions[header] = convert_to_int(sequence_data[column_index])
	positions = get_regions_position(regions)
	for header in positions:
		column_index = columns_header.index(header)
		sequence_data[column_index] = positions[header]
	
	for header in identity:
		column_index = columns_header.index(header)
		sequence_data[column_index] = convert_to_float(sequence_data[column_index])
	
	occurrence_index = columns_header.index('occurrence')
	sequence_data[occurrence_index] = convert_to_int(sequence_data[occurrence_index])
	
	return sequence_data

#-------------------------------------------------------------------#

def get_clones_and_sub_clone(clusters, sequences):
	"""
	Retrieve clustered sequence information from AIRR file, and group this sequences into sub-clones (same V, J and CDR3 AA). 
	input clusters:	Dict()		key=cluster label, value=list of sequence IDs within the cluster
	input sequences:	Dict()		key=sequence ID, value=list containing sequence information
	output repetoire:	Dict()		key=cluster label, value=dict of sub-clones (key=V,J and cdr3 id, value=dict of sequences)
	"""
	
	index_vcall = columns_header.index('v_call')
	index_jcall = columns_header.index('j_call')
	index_cdr3aa = columns_header.index('cdr3_aa')
	
	seq_unknown = []
	repertoire = {}
	
	for c in clusters:
	
		seq_ids = clusters[c]
		repertoire[c] = {}
	
		for seq in seq_ids :
	
			if seq in sequences:
				infos = format_data(sequences[seq])
				
				v_j_cdr3 = infos[index_vcall]+"_"+infos[index_jcall]+"_"+infos[index_cdr3aa]
	
				if v_j_cdr3 in repertoire[c]:
					repertoire[c][v_j_cdr3][seq] = infos
				else:
					repertoire[c][v_j_cdr3] = {}
					repertoire[c][v_j_cdr3][seq] = infos
				del sequences[seq]
	
			else:
				seq_unknown.append((c,seq))

	return repertoire

#-------------------------------------------------------------------#

def sort_clones_and_sub_clones_by_abundance(repertoire):
	"""
	Sort clones and sub-clones in descending order according to their abundance.
	input repertoire:	Dict()		key=cluster label, value=dict of sub-clones (key=V,J and cdr3 id, value=dict of sequences)
	output repertoire:	Dict()		(key=cluster label, value=dict of sub-clones (key=V,J and cdr3 id, value=dict of sequences)) -> the data are ordered by abundance
	"""
	
	index_occurrence = columns_header.index('occurrence')
	
	for clone in repertoire.keys() :
	
		for sub_clone in repertoire[clone].keys():
	
			if(len(repertoire[clone][sub_clone].keys()) > 1):
				repertoire[clone][sub_clone] = {k: v for k, v in sorted(repertoire[clone][sub_clone].items(), key=lambda item: item[1][index_occurrence], reverse=True)}
			repertoire[clone][sub_clone]["total"] = sum(v[index_occurrence] for v in repertoire[clone][sub_clone].values())
	
		repertoire[clone] = {k: v for k, v in sorted(repertoire[clone].items(), key=lambda item: item[1]["total"], reverse=True)}
		repertoire[clone]["total"] = sum(v["total"] for v in repertoire[clone].values())
	
	repertoire = {k: v for k, v in sorted(repertoire.items(), key=lambda item: item[1]["total"], reverse=True)}
	
	return repertoire	

#-------------------------------------------------------------------#

def write_repertoire_informations(repertoire, file_name):
	"""
	Write clone and sub-clone information to one file and clonotype information to another file.
	input repertoire:	Dict()		key=cluster label, value=dict of sub-clones (key=V,J and cdr3 id, value=dict of sequences)
	input file_name:	str		name of the output file
	"""

	repertoire_file = file_name + '_repertoire_two_levels_info.txt'
	clonotype_file = file_name + '_clonotypes_proportion.txt'
	
	with open (repertoire_file, 'w') as f1, open(clonotype_file, 'w') as f2:
	
		for clone in repertoire:
			sub_clone_id = 0
			line = "Clone "+clone+"\t"+str(repertoire[clone]["total"])+"\t"
			most_ab_sub_clone = next(iter(repertoire[clone]))
			most_ab_clonotype = next(iter(repertoire[clone][most_ab_sub_clone]))
			line += repertoire[clone][most_ab_sub_clone][most_ab_clonotype][columns_header.index("v_call")]+"\t"+repertoire[clone][most_ab_sub_clone][most_ab_clonotype][columns_header.index("j_call")]+"\t"
			
			for sub_clone in repertoire[clone]:
				seq_found = False
				if sub_clone != 'total' :
					clonotypes = clone+"_"+str(sub_clone_id)+"\t"
					
					for clonotype in repertoire[clone][sub_clone]:
						if clonotype != 'total' :
							clonotype_info = repertoire[clone][sub_clone][clonotype]
							
							if (not seq_found) and (clonotype_info[columns_header.index("whole_seq")] != "") and (clonotype_info[columns_header.index("germline_seq")] != ""):
								representative_seq = clonotype
								seq_found = True
								
							clonotypes += clonotype+","+str(clonotype_info[columns_header.index("occurrence")])+" "
					f2.write(clonotypes.rstrip()+"\n")
					
					if not seq_found:
						representative_seq = next(iter(repertoire[clone][sub_clone]))

					infos = repertoire[clone][sub_clone][representative_seq]
					line += "Sub-clone:"+clone+"_"+str(sub_clone_id)+","+representative_seq+","+str(repertoire[clone][sub_clone]["total"])+","+infos[columns_header.index("productive")]+","+infos[columns_header.index("cdr3_aa")]+","+str(infos[columns_header.index("v_sequence_start")])+"-"+str(infos[columns_header.index("cdr1_start")])+"-"+str(infos[columns_header.index("cdr1_end")])+"-"+str(infos[columns_header.index("cdr2_start")])+"-"+str(infos[columns_header.index("cdr2_end")])+"-"+str(infos[columns_header.index("cdr3_start")])+"-"+str(infos[columns_header.index("v_sequence_end")])+"-"+str(infos[columns_header.index("d_sequence_start")])+"-"+str(infos[columns_header.index("d_sequence_end")])+"-"+str(infos[columns_header.index("j_sequence_start")])+"-"+str(infos[columns_header.index("cdr3_end")])+"-"+str(infos[columns_header.index("j_sequence_end")])+","+'{0:g}'.format(infos[columns_header.index("v_identity")])+","+'{0:g}'.format(infos[columns_header.index("d_identity")])+","+'{0:g}'.format(infos[columns_header.index("j_identity")])+","+str(infos[columns_header.index("whole_seq")])+","+str(infos[columns_header.index("germline_seq")])+" "
					
					sub_clone_id += 1

			f1.write(line.rstrip()+"\n")

#-------------------------------------------------------------------#

def main():
	start_time = time.time()
	usage = "usage: two_level_info.py -a AIRR_file -c clustering_file -o output -g gene_alleles"
	parser = OptionParser(usage)
	parser.add_option("-a", "--AIRR_file", dest="airr_file", help="file in AIRR format")
	parser.add_option("-c", "--clustering_file",dest="clustering_file", help="read data from clustering_file")
	parser.add_option("-o", "--output", dest="output", help="path and name of the output file")
	parser.add_option("-g", "--alleles_taken_into_account",dest="gene_alleles", help="0 if V and J gene alleles are not taken into account for the clustering")

	(options, args) = parser.parse_args()
	if len(sys.argv) != 9:
		parser.error("incorrect number of arguments")

	airr_file = options.airr_file
	clustering_file = options.clustering_file
	output_name = options.output
	gene_alleles = options.gene_alleles
	
	sequences = read_AIRR_file(airr_file, gene_alleles)
	clusters = get_cluster_seq(clustering_file)
	
	repertoire = get_clones_and_sub_clone(clusters, sequences)
	repertoire = sort_clones_and_sub_clones_by_abundance(repertoire)
	
	write_repertoire_informations(repertoire, output_name)
#{k: v for k, v in sorted(x.items(), key=lambda item: item[1][index_occurrence])}

	print("  Two level informations step execution time : %s seconds " % (time.time() - start_time))
	
#-------------------------------------------------------------------#

if __name__ == "__main__":
	main()
