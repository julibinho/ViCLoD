"""
Author: Lucile Jeusset, Laboratory of Computational and Quantitative Biology, UPMC, Paris.
Email: lucile_jeusset@hotmail.fr

Adds to the tree in newick format the information of the sub-clones needed to visualize the tree (color, name, abundance value, cdr3, ...)

"""

import sys
import os
import json
import time
from optparse import OptionParser
from newick import load


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

def read_json_file(file_name) :
	"""
	Put information of a JSON file in a dict
	input file_name:	str		path of the JSON file
	output data:		Dict()		information of the file in a dict
	"""
	data = {}
	if (os.path.isfile(file_name)) :
		with open(file_name) as f :
			data = json.load(f)
	return data

# ------------------------------------------------------------------------------------

def get_subclones_id(file_name, nbr_clone) :
	"""
	Get the correspondence between the id of the JSON file and the id given in MobiLLe for the subclones
	input file_name:		str		path of the file with the correspondance
	input nbr_clone:		int		number of clone to analyze
	output selected_subclones:	Dict()		key=ID of clone to analyze (ViCloD), value=dict (key=sub-clone ID in MobiLLe, value=sub-clone ID in ViCloD)
	"""
	selected_subclones = {}
	
	f = open(file_name, 'r')
	for i in range(nbr_clone):
		line = f.readline()
		columns = (line.split("\n")[0]).split("\t")
		clone_id = columns[0].split(":")
		if len(clone_id) >= 2:
			selected_subclones[clone_id[1]] = {}
			if len(columns) >= 2:
				subclone = columns[1].split(" ")
				for s in range(len(subclone)):
					subclone_id = subclone[s].split(":")
					if len(subclone) >= 2:
						selected_subclones[clone_id[1]][subclone_id[0]] = subclone_id[1]
	f.close()
	
	return selected_subclones

# ------------------------------------------------------------------------------------

def get_subclone_info_from_repertoire(repertoire, selected_subclones, complexity) :
	"""
	Retrieve information about clones to be analyzed and subclones from repertoire data (JSON file)
	input repertoire:		Dict()		information on clone and sub-clone (name, ab, idV, idJ, cdr3, ...)
	input selected_subclones:	Dict()		key=ID of clone to analyze (ViCloD), value=dict (key=sub-clone ID in MobiLLe, value=sub-clone ID in ViCloD)
	input complexity:		str		complexity of the tree (all, simplification1 or simplification2)
	output repertoire_tree:	list()		list of dict containing clone
	output clones_data:		Dict()		key=sub-clone ID (ViCloD), value=information on this sub-clone (name, ab, idV, idJ, cdr3, ...)
	"""
	clone_ids = [k for k in selected_subclones.keys()]
	nbr_clone = len(repertoire['children'])
	i = 0
	repertoire_tree = []
	subclones_data = {}
	
	if 'children' in repertoire:	
		while i<nbr_clone and clone_ids:
			name = repertoire['children'][i]['name']
			if name in clone_ids:
				repertoire_tree.append({"name": repertoire["children"][i]["name"], "idV": repertoire["children"][i]["idV"], "idJ": repertoire["children"][i]["idJ"], "ab": repertoire["children"][i]["ab"], "reads": repertoire["children"][i]["reads"], "tree_type": complexity})
				if 'children' in repertoire['children'][i] :
					for subclone in repertoire['children'][i]['children'] :
						subclones_data[subclone["name"]] = subclone
				index = clone_ids.index(name)
				del clone_ids[index]
			i += 1
	
	return repertoire_tree, subclones_data

# ------------------------------------------------------------------------------------

def get_subclones_info_from_repertoire(clone, subclones_data):
	for subclone in clone :
		subclones_data[subclone["name"]] = subclone
# ------------------------------------------------------------------------------------

def read_newick_file(file_name) :
	"""
	Read file in newick format and return the tree
	input file_name:	str		path of the newick file
	output:		Dict()		tree
	"""
	with open(file_name) as f :
		return load(f)

# ------------------------------------------------------------------------------------

def get_tree_branches(file_name) :
	"""
	Get the length of the tree branches
	input file_name:	str		path of the csv file containing the length of the tree branches
	output branches:	Dict()		key=(name of node1, name of node2) (MobiLLe ID), value=distance between the two nodes
	"""
	global branches
	branches = {}
	branch_length = ["","",0]
	if (os.path.isfile(file_name)) :
		with open(file_name) as f :
			for line in f :
				branch_length = line.split("\n")[0].split(",")
				branches[(branch_length[0],branch_length[1])]=branch_length[2]

# ------------------------------------------------------------------------------------

def get_aligned_sequences(file_name):
	"""
	Get the aligned sequences of the sub-clones and the position of the V,D,J regions on these sequences
	input file_name:	str		path of the file containing the aligned sequences of the sub-clones and the position of the different regions on these sequences  
	output sequences:	Dict()		key=sub-clone id (MobiLLe), value=dict with the aligned sequence and a list of the positions of the different regions
	"""
	global sequences
	sequences = {}
	if (os.path.isfile(file_name)) :
		with open(file_name) as f :
			for line in f :
				properties = (line.split("\n")[0]).split("\t")
				if (len(properties) >= 3) :
					regions = properties[2].split("-")
					if(len(regions) >= 12):
						seq_id = properties[0]
						aligned_seq = properties[1]
						sequences[seq_id] = {'aligned_seq': aligned_seq, 'regions': [int(x) for x in regions]}

# ------------------------------------------------------------------------------------

def format_regions_info(regions):
	"""
	Format the list containing information concerning the position of the genes
	input regions:		list()		list containing the position of the region
	ouput positions:	Dict()		key=name of the region, value=position in the nucleotide sequence
	"""
	positions = {}
	if len(regions) >= 12:
		positions["V_region"] = [regions[0],regions[6]]
		positions["CDR1"] = [regions[1],regions[2]]
		positions["CDR2"] = [regions[3],regions[4]]
		positions["CDR3"] = [regions[5],regions[10]]
		positions["D_region"] = [regions[7],regions[8]]
		positions["J_region"] = [regions[9],regions[11]]
	return positions	

# ------------------------------------------------------------------------------------

def branch_length (prev_node, current_node, silent_node) :
	"""
	Return branch length between two nodes in the tree
	input current_node:	str		ID of the node (MobiLLe)
	input prev_node:	str		ID of the previous node in the tree (MobiLLe)
	input silent_node:	int		number of silent node
	output:		str		branch length between current_node and prev_node
	"""
	if((prev_node, current_node) in branches.keys()):
		length = str(float(branches[(prev_node, current_node)])-silent_node)
		del branches[(prev_node, current_node)]
		return length
	else:
		return "0"

# ------------------------------------------------------------------------------------

def add_subclone_info(tree, clone_reads, current_node, prev_node, silent_node) :
	"""
	Add information to a sub-clone (aligned sequence, position of genes, branch length, unique sequences)
	input tree:		Dict()		tree located at node level
	input clone_reads:	int		sequence number of the clone from which these sub-clones are descended
	input current_node:	str		ID of the node (MobiLLe)
	input prev_node:	str		ID of the previous node in the tree (MobiLLe)
	input silent_node:	int		number of silent node
	"""
	if current_node in sequences:
		del tree["sequence"]
		del tree["germline"]
		tree["aligned_seq"] = sequences[current_node]["aligned_seq"]
		tree["region"] = format_regions_info(sequences[current_node]["regions"])
		sequences[current_node]["occurrence"] = tree["reads"]
		sequences[current_node]["name"] = tree["name"]
		
	tree["length"] = branch_length (prev_node, current_node, silent_node) 
	

# ------------------------------------------------------------------------------------

def add_germline_info (tree) :
	"""
	Add information to the germline node in the tree
	input tree:		Dict()		tree located at node level (germline)
	"""
	tree["name"]="HAUS"
	tree["value"]=1
	tree["color"]="#808080"
	tree["stroke"]="#808080"
	tree["style"]="none"
	tree["seq"]="germline"
	tree["uniq_seq"]=[{"seq_id":"germline","percentage":"100","rank":"0"}]
	if "germline" in sequences:
		tree["aligned_seq"] = sequences["germline"]["aligned_seq"]
		tree["region"] = format_regions_info(sequences["germline"]["regions"])
		sequences["germline"]["occurrence"] = 1
		sequences["germline"]["name"] = "HAUS"

# ------------------------------------------------------------------------------------

def add_silent_node_info(tree, node_name):
	"""
	Add information to the silent node in the tree
	input tree:		Dict()		tree located at node level (silent node)
	input node_name:	str		name of the silent node
	"""
	tree["value"]=1
	tree["name"]= node_name
	tree["seq"]= node_name
	tree["color"]="#FFFFFF"
	tree["stroke"]="#808080"
	tree["style"]="2,2"
	tree["length"]=1
	tree["uniq_seq"]=[{"seq_id":node_name,"percentage":"100","rank":"0"}]

# ------------------------------------------------------------------------------------

def add_info_to_tree(newick_tree, clones_data, subclone_ids, clone_reads, silent_node, prev_node):
	"""
	Add information of the sub-clones to newick tree
	input newick_tree:		Dict()		tree
	input clones_data:		Dict()		key=sub-clone ID (ViCloD), value=information on this sub-clone (name, ab, idV, idJ, cdr3, ...)
	input subclones_id:		Dict()		key=sub-clone ID in MobiLLe, value=sub-clone ID in ViCloD
	input clone_reads:		int		sequence number of the clone from which these sub-clones are descended
	input silent_node:		int		number of silent node	
	input prev_node:		str		ID of the previous node in the tree (MobiLLe)
	output tree:			Dict()		tree containing the modifications (key="children", value=list of dict representing sub-clones)
	"""
	tree = {}
	#subclone_ids = [k for v in subclone_ids.values() for k in v]
	
	if newick_tree.name in subclone_ids.keys():
		subclone = clones_data[subclone_ids[newick_tree.name]]
		tree = {**tree, **subclone}
		
		add_subclone_info(tree, clone_reads, newick_tree.name, prev_node, silent_node)
		del subclone_ids[newick_tree.name]
		silent_node = 0
		prev_node = newick_tree.name
		
	elif(newick_tree.name=="germline"):
		add_germline_info(tree)
		prev_node = newick_tree.name

	else :
		add_silent_node_info(tree, newick_tree.name)
		silent_node+=1

	if(newick_tree.descendants!=[]) :
		tree["children"]=newick_tree.descendants	
		for i in range(len(newick_tree.descendants)):	 
			subclone_info = add_info_to_tree(newick_tree.descendants[i], clones_data, subclone_ids, clone_reads, silent_node, prev_node)
			tree["children"][i] = subclone_info	
	return tree
			
# ------------------------------------------------------------------------------------

def save_json_format(file_name, data):
	"""
	Write in a file the data in JSON format
	input data:	Dict()		
	"""
	output_file = open(file_name+".json", "w")
	json.dump(data, output_file, indent = 4, sort_keys = False)
	output_file.close()
	
# ------------------------------------------------------------------------------------

def save_aligned_sequences_in_fasta(file_name1, file_name2):
	"""
	Write aligned sequences of the sub-clones in fasta format (these 2 files will be used to generate the VDJ logo)
	input file_name1:	str	path of the fasta file which will contain the HAUS sequence
	input file_name2:	str	path of the fasta file which will contain all the sequences of the sub-clones
	"""
	with open (file_name1+".fasta", 'w') as f1, open (file_name2+".fasta", 'w') as f2:
		for subclone in sequences:
			if "name" in sequences[subclone]:
				subclone_name = sequences[subclone]["name"]
				if subclone == "germline":
					f1.write(">" + subclone_name + "\n" + sequences[subclone]["aligned_seq"] + "\n")
				else :
					line = ">" + subclone_name + "@"
					if sequences[subclone]["occurrence"] and convertible_into_int(sequences[subclone]["occurrence"]):
						line += sequences[subclone]["occurrence"] + "\n"
					else:
						line += "1\n"
					line += sequences[subclone]["aligned_seq"] + "\n"
					f2.write(line)

# ------------------------------------------------------------------------------------

def main():

	start_time = time.time()
	usage = usage = "python subclone_information_json.py -c <number_clone> -s <number_of_subclone> -o <output> -a <analysis_id> -n <newick_file> -t <tree_complexity> \n"
	parser = OptionParser(usage)
	parser.add_option("-c", "--number_clone", dest="number_clone",  help="number of clone to analyze")
	parser.add_option("-s", "--number_subclone", dest="number_subclone", help="the number of sub-clone to analyze, if use all, there will be no selection for the analyzing sub-clones otherwise choose the number of n most abondant sub-clones ")
	parser.add_option("-o", "--output", dest="output",  help="path of the output folder")
	parser.add_option("-a", "--analysis_id", dest="analysis_id",  help="id of the analysis")
	parser.add_option("-n", "--newick_file", dest="newick_file",  help="end of newick file")
	parser.add_option("-t", "--tree_complexity", dest="tree_complexity",  help="complexity of the tree (all, simplification1 or simplification2)")
	
	(options, args) = parser.parse_args()
	
	if len(sys.argv) != 13:
		parser.error("Incorrect number of arguments")
	
	nbr_clone = convert_parameter_to_int(options.number_clone,5)
	nbr_subclone = convert_parameter_to_int(options.number_subclone,200)
	path = options.output
	analysis_name = options.analysis_id
	end_newick_file = options.newick_file
	complexity = options.tree_complexity
	file_path = path + '/' + analysis_name
	
	repertoire = read_json_file(file_path + "_repertoire.json")
	selected_subclone_ids = get_subclones_id(file_path + "_matching_names_MobiLLe_ViCloD.txt",nbr_clone)
	
	repertoire_tree, clones_data = get_subclone_info_from_repertoire(repertoire, selected_subclone_ids, complexity)
	
	for i in range(len(repertoire_tree)):
		tree_file = file_path + "_" + str(i+1) + "_" + end_newick_file

		if (os.path.isfile(tree_file)) :
			
			clone_name = repertoire_tree[i]["name"]
			newick_tree = read_newick_file(tree_file)
			get_tree_branches(file_path + "_" + str(i+1) + "_" + str(nbr_subclone) + ".nk.csv")
			get_aligned_sequences(file_path + "_" + str(i+1) + "_" + str(nbr_subclone) + "_aligned_regions.txt")

			tree = add_info_to_tree(newick_tree[0].descendants[0], clones_data, selected_subclone_ids[clone_name], int(repertoire_tree[i]["reads"]), 0, "germline")
			repertoire_tree[i]["children"] = tree
			
			save_json_format(file_path + "_tree_" + complexity, repertoire_tree)

			if (complexity=="simplification2") :
				save_aligned_sequences_in_fasta(file_path+"_"+clone_name+"_HAUS_sequence", file_path+"_"+clone_name+"_subclones_sequences")
			
			print("    Tree(" + complexity + ") in json format of clone C"+str(i+1)+" has been generated.")
		
	print("  Generation of sub-clones files for visualization : done. Execution time : %s seconds " % (time.time() - start_time))

#===================================================================================
if __name__ == "__main__":
	main()

