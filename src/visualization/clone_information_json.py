"""
Author: Lucile Jeusset, Laboratory of Computational and Quantitative Biology, UPMC, Paris.
Email: lucile_jeusset@hotmail.fr

Translate the txt format of the file containing clone whith their sub-clones in json format. It recover the output file provide from MobiLLe analysis.
Three files must be pass in argument : the file to convert in txt format, the txt file containing the total number of reads and a file containning the html colors to be assigned to the clones and sub-clones.

"""

import sys
import os
import time
from optparse import OptionParser
import json

#Match between the lines of the log file and the information sought
log_file_line = {"total":1, "analyzed":2, "removed":3, "unannotated_for_clustering":4, "unannotated_for_tree":5, "low_quality":6, "low_abundance":7, "cumulated_criteria":8, "only_unannotated_for_tree":9}

# ------------------------------------------------------------------------------------

def convertible_into_int(value):
	"""
	Check if a string can be convert to integer
	input value:	str	string to convert
	output:	bool	True if the string can be convert to int
	"""
	try:
		int(value)
		return True
	except ValueError:
		return False

# ------------------------------------------------------------------------------------

def color_list(file_name):
	"""
	Put the color code contained in the file in a list
	input file_name:		str		path of the file containing the color code
	output colors:		list()		color code
	"""
	colors = []
	with open(file_name) as f:
		for line in f :
			colors.append(line.split("\n")[0])
	return colors

# ------------------------------------------------------------------------------------

def change_index(tab_index,nb_colors,nb_stroke,reverse):
	"""
	Increment the index used to determine the fill color, stroke style and stroke color (employed to represent a clone or subclone)
	input tab_index:	list()		list containing the 3 indices used to respectively determine the fill color, the stroke color and the stroke style
	input nb_colors:	int		maximum number of colors available
	input nb_stroke:	int		maximum number of stroke style available
	input reverse:		bool		determine the order in which the colors should be assigned
	output tab_index:	list()		list containing the 3 indices modified
	"""
	#reading direction of the color table
	if(reverse):
		last_color = 0
		first_color = nb_colors-1
		new_color = tab_index[0] - 1
		new_stroke = tab_index[1] + 1
	else:
		last_color = nb_colors-1
		first_color = 0
		new_color = tab_index[0] + 1
		new_stroke = tab_index[1] - 1
	#modify the index 
	if(tab_index[0]!=last_color):
		tab_index[0] = new_color
	elif(tab_index[1]!=first_color):
		tab_index[0] = first_color
		tab_index[1] = new_stroke
	elif(tab_index[2]!=nb_stroke-1):
		tab_index[0] = first_color
		tab_index[1] = last_color
		tab_index[2] += 1
	else:
		tab_index[0] = first_color
		tab_index[1] = last_color
		tab_index[2] = 0

	return tab_index

# ------------------------------------------------------------------------------------

def get_proportion_uniq_seq (file_name) :
	"""
	Get occurrence of unique sequence in all the subclones
	input file_name:	str		path of the file containing the unique sequences id of all the sub-clones
	output clonotypes:	Dict()		key=subclone id (MobiLLe), value=dict of clonotype sequence ids and their abundance
	"""
	clonotypes = {}
	
	if (os.path.isfile(file_name)):
		with open (file_name, 'r') as f:
			for line in f:
				columns = (line.split("\n")[0]).split("\t")
				if len(columns) >= 2:
					seq = columns[1].split(" ")
					unique_seq = {}
					for s in range(len(seq)):
						seq_info = seq[s].split(",")
						if len(seq_info) >= 2 :
							unique_seq[seq_info[0]] = seq_info[1]
						if unique_seq:
							clonotypes[str(columns[0])] = unique_seq

	return clonotypes

# ------------------------------------------------------------------------------------

def total_number_of_reads(log_file):
	"""
	Retrieve statistical information about the data provide in the AIRR file
	input log_file:	str	path of the file containing sequences information
	output reads_infos:	Dict()	key=criteria, value=number of sequences corresponding to that citeria
	"""
	reads_info = {"total":1, "analyzed":1, "removed":0, "unannotated_for_clustering":0, "unannotated_for_tree":0, "low_quality":0, "low_abundance":0, "cumulated_criteria":[0,0,0,0], "only_unannotated_for_tree":0}
	line_with_two_info = ["analyzed", "removed", "unannotated_for_clustering", "unannotated_for_tree", "low_quality", "low_abundance", "only_unannotated_for_tree"]

	if (os.path.isfile(log_file)) : 
		with open(log_file) as file:
			lines = file.readlines()
		for info in line_with_two_info :
			line_info = lines[log_file_line[info]-1].split("\n")[0]
			line_value = (line_info.split(": ")[1]).split("\t")
			reads_info[info] = int(line_value[0])
		
		reads_info["total"] = int((lines[log_file_line["total"]-1].split("\n")[0]).split(": ")[1])
		
		criteria = (lines[log_file_line["cumulated_criteria"]-1].split("\n")[0]).split("\t")
		for i in range(len(criteria)) :
			value = criteria[i].split(";")
			reads_info["cumulated_criteria"][i] = int(value[1])

	return reads_info

# ------------------------------------------------------------------------------------

def get_percentage_of_clonotypes_in_subclone(subclone_id, clonotypes, nbr_subclone):
	"""
	Get percentage of clonotype for a given subclone
	input subclone_id:	str		id of the subclone (MobiLLe)
	input clonotypes:	Dict()		key=subclone id (MobiLLe), value=dict of clonotype sequence ids and their abundance
	input nbr_subclone:	int		number of sequence in the subclone
	ouptut proportion:	list()		list of clonotypes in this subclone
	"""
	proportion=[]
	if(subclone_id in clonotypes.keys()):
		rank = 0
		for seq in clonotypes[subclone_id] :
			if (convertible_into_int(clonotypes[subclone_id][seq])):
				percentage = (int(clonotypes[subclone_id][seq])*100)/nbr_subclone
				proportion.append({"seq_id":seq,"reads":clonotypes[subclone_id][seq],"percentage":percentage,"rank":rank})
				rank += 1
	if len(proportion)==0:
		proportion.append({"seq_id":subclone_id,"reads":nbr_subclone,"percentage":"100","rank":0})
	return proportion

# ------------------------------------------------------------------------------------

def get_repertoire_informations(file_name, colors, reads, clonotypes):
	"""
	Get clones and sub-clones information from the repertoire_two_levels file
	input file_name:		str		path of the file containing repertoire information 
	input colors:			list()		color code
	input reads_infos:		Dict()		key=criteria, value=number of sequences corresponding to that citeria
	input clonotypes:		Dict()		key=subclone id (MobiLLe), value=dict of clonotype sequence ids and their abundance
	output repertoire:		Dict()		3 level: first -> dict of repertoire information, second -> dict of clone information, third -> dict of sub-clone information
	output repertoire_name:	Dict()		key=cluster id, value=dict containing clone id and a dict of sub-clone ids
	"""
	stroke_style = ["none", "5,5", "1,5"]
	clone_index = [0,len(colors)-1,0]
	clones=[]
	num_clone = 1;
	productivity = {"T": "productive","F" : "unproductive", "" : "/"}
	matching_name = {}

	if (os.path.isfile(file_name)) : 
		with open(file_name) as f:
			for line in f:
				line = line.split("\n")
				[clone_id, occurrence, v_call, j_call, subclones]=line[0].split("	")
				clone_name = "C"+str(num_clone)
				cluster_id = clone_id.split(" ")[1]
				matching_name[cluster_id] = {"viclod_name" : clone_name}
				subclones_list = subclones.split(" ")
				clone_abundance = (int(occurrence)/reads["analyzed"])*100
				#save the abundance of the clone with 3 decimals or in scientific format if it's inferior to 0.001
				if(clone_abundance<0.01):
					ab = "%.2e"%clone_abundance
				else :
					ab = round(clone_abundance,2)
				
				#the clone doesn't have sub-clones
				if(len(subclones_list)==1):
					[subclone_id, seq_id, sc_occurrence, productive, cdr3_aa, regions, v_identity, d_identity, j_identity, whole_seq, germline_seq] = subclones_list[0].split(",")
					uniq_seq = get_percentage_of_clonotypes_in_subclone(subclone_id.split(":")[1], clonotypes, int(sc_occurrence))
					[v_start, cdr1_start, cdr1_end, cdr2_start, cdr2_end, cdr3_start, v_end, d_start, d_end, j_start, cdr3_end, j_end] = regions.split("-")
					clones.append({"name" : clone_name, "value" : clone_abundance, "reads" : occurrence, "ab" : ab, "idV" : v_call, "idJ" : j_call, "cdr3" : cdr3_aa, "functionality" : productivity[productive], "seq" : seq_id, "v_identity" : v_identity, "d_identity" : d_identity, "j_identity" : j_identity, "region":{"V":[v_start,v_end], "D":[d_start,d_end], "J":[j_start,j_end], "CDR1":[cdr1_start,cdr1_end], "CDR2":[cdr2_start,cdr2_end], "CDR3":[cdr3_start,cdr3_end]}, "uniq_seq": uniq_seq, "sequence" : whole_seq, "germline": germline_seq, "color" : colors[clone_index[0]], "stroke" : colors[clone_index[1]], "style" : stroke_style[clone_index[2]]})
					clone_index = change_index(clone_index,len(colors),len(stroke_style),False)	#update the values of the index
				#the clone have sub-clones
				else:
					subclone_info, subclones_id = get_subclones_in_clone(clone_name, colors, [len(colors)-1,clone_index[0],0], subclones_list, int(occurrence), reads["analyzed"], clonotypes)

					matching_name[cluster_id]["subclones_id"] = subclones_id
					clones.append({"name" : clone_name, "value" : clone_abundance, "reads" : occurrence, "ab" : ab, "idV" : v_call, "idJ" : j_call, "cdr3" : subclone_info[0]["cdr3"], "functionality" : subclone_info[0]["functionality"], "color" : colors[clone_index[0]], "stroke" : colors[clone_index[1]], "style" : stroke_style[clone_index[2]], "children" : subclone_info})
					clone_index = change_index(clone_index,len(colors),len(stroke_style),False)	#update the values of the index
				num_clone += 1

	repertoire = {"name" : "repertoire", "value" : 100, "reads": reads, "color":"#808080", "stroke":"#808080", "style":"none", "children": clones}	
	
	return repertoire, matching_name

# ------------------------------------------------------------------------------------

def get_subclones_in_clone(clone_name, colors, subclone_index, subclones_list, nbr_seq_clone, nbr_seq_repertoire, clonotypes):
	"""
	Get sub-clone information from file
	input clone_name:		str		name of the clone to which the sub-clones belong
	input colors:			list()		color code		
	input subclone_index:		list()		list containing the 3 indices used to respectively determine the fill color, the stroke color and the stroke style
	input subclones_list:		list()		list of all sub-clones contained in the clone
	input nbr_seq_clone:		int		total number of sequences in the clone
	input nbr_seq_reperoire:	int		total number of sequence in the repertoire
	input clonotypes:		Dict()		key=subclone id (MobiLLe), value=dict of clonotype sequence ids and their abundance
	output subclones:		list()		list containing a dictionary for each sub-clone with its information
	output matching_name:		Dict()		key=name of sub-clone in MobiLLe, value=name of sub-clone in ViCloD
	"""
	
	matching_name = {}
	stroke_style = ["none", "5,5", "1,5"]
	productivity = {"T": "productive","F" : "unproductive", "" : "/"}
	subclones = []
	
	for i in range(len(subclones_list)):
		[subclone_id, seq_id, sc_occurrence, productive, cdr3_aa, regions, v_identity, d_identity, j_identity, whole_seq, germline_seq] = subclones_list[i].split(",")
		[v_start, cdr1_start, cdr1_end, cdr2_start, cdr2_end, cdr3_start, v_end, d_start, d_end, j_start, cdr3_end, j_end] = regions.split("-")
		
		subclone_name = clone_name+"-"+str(i+1)
		matching_name[subclone_id.split(":")[1]] = subclone_name
		rep_ab = (int(sc_occurrence)/nbr_seq_repertoire)*100
		clone_ab = (int(sc_occurrence)/nbr_seq_clone)*100
		#save the abundance of the clonotype in the repertoire with 3 decimals or in scientific format if it's inferior to 0.001
		if(rep_ab<0.01):
			ab_subclone_rep = "%.2e"%(rep_ab)
		else :
			ab_subclone_rep = round(rep_ab,2)
		#save the abundance of the clonotype in the clone with 3 decimals or in scientific format if it's inferior to 0.001
		if(clone_ab<0.01):
			ab_subclone = "%.2e"%(clone_ab)
		else :
			ab_subclone = round(clone_ab,2)
		
		uniq_seq = get_percentage_of_clonotypes_in_subclone(subclone_id.split(":")[1], clonotypes, int(sc_occurrence))
		
		subclones.append({"name":subclone_name, "value":clone_ab, "reads":sc_occurrence, "ab":ab_subclone, "ab_rep":ab_subclone_rep, "functionality":productivity[productive], "cdr3":cdr3_aa, "seq":seq_id, "v_identity" : v_identity, "d_identity" : d_identity, "j_identity" : j_identity, "region":{"V":[v_start,v_end], "D":[d_start,d_end], "J":[j_start,j_end], "CDR1":[cdr1_start,cdr1_end], "CDR2":[cdr2_start,cdr2_end], "CDR3":[cdr3_start,cdr3_end]}, "uniq_seq": uniq_seq, "sequence" : whole_seq, "germline": germline_seq, "color":colors[subclone_index[0]], "stroke":colors[subclone_index[1]], "style":stroke_style[subclone_index[2]]})
		
		subclone_index = change_index(subclone_index,len(colors),len(stroke_style),True)	#update the values of the index

	return subclones, matching_name

# ------------------------------------------------------------------------------------

def save_json_format(file_name, data):
	"""
	Write data to a file in JSON format
	input file_name:	str		path of file to write
	input data:		Dict()		data to save in JSON format		
	"""
	sorted_output_file = open(file_name+".json", "w")
	json.dump(data, sorted_output_file, indent = 4, sort_keys = False)
	sorted_output_file.close()

# ------------------------------------------------------------------------------------

def save_matching_names(file_name, matching_names) :
	"""
	Write in a file the correspondence between the clone and sub-clone names given by MobiLLe and those used by ViCloD
	input file_name:	str	path of file to write
	input matching_names:	Dict()	correspondance between the name given by MobiLLe and ViCloD		
	"""
	with open ( file_name, 'w') as f :
		for name in matching_names:
			line = name + ":" + matching_names[name]["viclod_name"] + "\t"
			if "subclones_id" in matching_names[name]:
				for subclone in matching_names[name]["subclones_id"] :
					line += subclone + ":" + matching_names[name]["subclones_id"][subclone] + " "
			f.write(line.rstrip() + "\n")
			
# ------------------------------------------------------------------------------------

def main():

	start_time = time.time()
	usage = usage = "python clone_information_json.py -r <repertoire_file> -l <log_file> -p <clonotype_proportion> -a <analysis_id> -c <colors_file> \n"
	parser = OptionParser(usage)
	parser.add_option("-r", "--repertoire_file", dest="repertoire_file",  help="txt file containing repertoire information")
	parser.add_option("-l", "--log_file", dest="log_file",  help="log file containing information on reads")
	parser.add_option("-p", "--clonotype_proportion", dest="clonotype_proportion",  help="file containing the proportion of clonotype of each subclone")
	parser.add_option("-a", "--analysis_id", dest="analysis_id",  help="id of the analysis")
	parser.add_option("-c", "--colors_file", dest="colors_file",  help="txt file containing a list of colors")
		
	(options, args) = parser.parse_args()
	
	if len(sys.argv) < 11:
		parser.error("Incorrect number of arguments")
	
	repertoire_file = options.repertoire_file
	log_file = options.log_file
	clonotype_file = options.clonotype_proportion
	analysis = options.analysis_id
	colors_file = options.colors_file

	reads = total_number_of_reads(log_file)
	colors = color_list(colors_file)
	clonotype_proportion = get_proportion_uniq_seq(clonotype_file)
	repertoire, matching_names = get_repertoire_informations(repertoire_file, colors, reads, clonotype_proportion)	

	path = os.path.dirname(repertoire_file)
	user_path = path + "/" + analysis
	repertoire_file_name = user_path + "_repertoire" 

	save_json_format(repertoire_file_name, repertoire)
	save_matching_names(user_path+"_matching_names_MobiLLe_ViCloD.txt", matching_names)

	print("  JSON file of the repertoire created. Execution time : %s seconds " % (time.time() - start_time))

# ------------------------------------------------------------------------------------

if __name__ == "__main__":
	main()
