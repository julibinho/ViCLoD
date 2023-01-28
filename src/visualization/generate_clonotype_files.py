#!/usr/bin/env python3
#-*-coding:Utf-8-*-

"""
Author: Lucile Jeusset, Laboratory of Computational and Quantitative Biology, UPMC, Paris.
Email: lucile_jeusset@hotmail.fr

Recover the informations contained for the tree representation (color, name, abundance value, cdr3, ...)

"""

import sys
import os
import json
from optparse import OptionParser
from newick import load

#transform the class format of the tree into dictionary format and add the abundance value
def change_format(newick_tree, num_clone, clonotype_seq_id, silent_node, prev_node, matching_name):
	tree = {}
	if newick_tree.name in clonotype_seq_id.keys():	#add the value of abundance if the name of the clone is find in the dictionary of abundance
		clonotype = repertoire["children"][num_clone]["children"][clonotype_seq_id[newick_tree.name]]
		add_clonotype_info(tree, clonotype, newick_tree.name, prev_node, silent_node, clonotype_seq_id)	
		silent_node = 0
		prev_node = newick_tree.name
		#add the name to the sequences of this clonotype
		if(newick_tree.name in sequences.keys()):
			sequences[newick_tree.name].insert(0,tree["name"])
			matching_name[newick_tree.name] = tree["name"]

	elif(newick_tree.name=="germline"):
		add_germline_info (tree)
		prev_node = newick_tree.name
		#add the name and the length informations to the sequence
		if(newick_tree.name in sequences.keys()):
			sequences[newick_tree.name].insert(0,"naive")
			matching_name[newick_tree.name] = "naive"

	else :
		add_silent_node_info (tree, newick_tree.name)
		silent_node+=1

	if(newick_tree.descendants!=[]) :
		tree["children"]=newick_tree.descendants	
		for i in range(len(newick_tree.descendants)):	 
			clonotype_info = change_format(newick_tree.descendants[i], num_clone, clonotype_seq_id, silent_node, prev_node, matching_name)	#repeat the same actions as its parents (name, value, length)
			tree["children"][i] = clonotype_info	
	return tree

# ------------------------------------------------------------------------------------

def add_clonotype_info (tree, clonotype, current_node, prev_node, silent_node, clonotype_seq_id) :
	tree["name"]= clonotype["name"]
	tree["value"]=str(clonotype["value"])
	tree["reads"]=clonotype["reads"]
	tree["ab"]=clonotype["ab"]
	tree["seq"]=current_node
	tree["ab_rep"]=clonotype["ab_rep"]
	tree["productivity"]=clonotype["productivity"]
	tree["cdr3"]=clonotype["cdr3"]
	tree["color"]=clonotype["color"]
	tree["stroke"]=clonotype["stroke"]
	tree["style"]=clonotype["style"]
	del clonotype_seq_id[current_node]	
	tree["length"] = branch_length (prev_node, current_node, silent_node) #add the length to the node's properties if the length exist in the csv file
	#add the percentage of each uniq seq
	if(current_node in clonotype_uniq_seq.keys()):
		tree["uniq_seq"]=[]
		rank = 0
		for seq in clonotype_uniq_seq[current_node] :
			tree["uniq_seq"].append({"seq_id":seq,"percentage":clonotype_uniq_seq[current_node][seq],"rank":rank})
			rank += 1
			
# ------------------------------------------------------------------------------------

def add_germline_info (tree) :
	tree["name"]="naive"
	tree["value"]=1
	tree["color"]="#808080"
	tree["stroke"]="#808080"
	tree["style"]="none"
	tree["seq"]="germline"
	tree["uniq_seq"]=[{"seq_id":"germline","percentage":"100","rank":"0"}]

# ------------------------------------------------------------------------------------

def add_silent_node_info (tree, node_name) :
	tree["value"]=1
	tree["name"]= node_name
	tree["seq"]= node_name
	tree["color"]="#FFFFFF"
	tree["stroke"]="#808080"
	tree["style"]="2,2"
	tree["length"]=1
	tree["uniq_seq"]=[{"seq_id":node_name,"percentage":"100","rank":"0"}]

# ------------------------------------------------------------------------------------

def branch_length (prev_node, current_node, silent_node) :
	if((prev_node, current_node) in branches.keys()):
		length = str(float(branches[(prev_node, current_node)])-silent_node)	#recover the length of the tree object
		del branches[(prev_node, current_node)]
		return length
	else:
		return "0"

# ------------------------------------------------------------------------------------

#modify the name of the sequences in the fasta file
def changeNameSequences(fasta_file, matching_name, regions):
	seq_name_info = ""
	seq_name = ""
	sequence_info = {}
	sequence = {}
	if (os.path.isfile(fasta_file)) :
		with open(fasta_file) as file:
			for line in file:
				if(line[0]=='>'):
					if seq_name != "":
						sequence[seq_name] += "\n"
					seq_id=line[1:].split('@')[0]
					if(seq_id in matching_name):
						seq_name = matching_name[seq_id]
						if seq_id in regions :
							seq_name_info = matching_name[seq_id] + "|" + "-".join(regions[seq_id])
						else :
							seq_name_info = matching_name[seq_id] + "|undefined-undefined-undefined-undefined-undefined-undefined"
						sequence_info[seq_name_info] = ""
						sequence[seq_name] = ""
					elif(seq_id=="germline\n"):
						seq_name_info = "naive" + "|" + "-".join(regions["germline"])
						seq_name = "naive"
						sequence_info[seq_name_info] = ""
						sequence[seq_name] = ""
					else:
						seq_name_info = ""
						seq_name = ""
				else:
					if(seq_name_info!=""):
						sequence_info[seq_name_info] += line
					if(seq_name!=""):
						sequence[seq_name] += line.split("\n")[0]
	return sequence_info, sequence

# ------------------------------------------------------------------------------------

def save_json_format(file_name, data):
	output_file = open(file_name+".json", "w")
	json.dump(data, output_file, indent = 4, sort_keys = False)
	output_file.close()

# ------------------------------------------------------------------------------------

def save_sequences(file_name, data):
	output_file = open(file_name+".txt", "w")
	for values in data.values():
		if len(values)==6 :
			output_file.write("	".join(values))
	output_file.close()

# ------------------------------------------------------------------------------------

#load repertoire information
def read_json_file (file) :
	global repertoire
	repertoire = {}
	if (os.path.isfile(file)) :
		with open(file) as f :
			repertoire = json.load(f)

# ------------------------------------------------------------------------------------

def read_newick_file (file) :
	with open(file) as f :
		return load(f)

# ------------------------------------------------------------------------------------

def match_sequence_id (sequences_id_file, clone_number) :
	clonotype_seq_id = {}
	matching_id = {}
	if (os.path.isfile(sequences_id_file)) :
		with open(sequences_id_file) as file :
			for line in file :
				ids = line.split("\n")[0].split("	")
				matching_id[ids[0]]=ids[1]
		#the clone has clonotype
		if("children" in repertoire["children"][clone_number-1]):
			for i in range(len(repertoire["children"][clone_number-1]["children"])) :
				seq = repertoire["children"][clone_number-1]["children"][i]["seq"] 
				if seq in matching_id :
					clonotype_seq_id[matching_id[seq]]=i 	#match the name and the sequence ID of a clonotype
	return clonotype_seq_id

# ------------------------------------------------------------------------------------

def find_clone (clone_number, tree_type) :
	clone_tree = []
	clones = []
	if (clone_number > len(repertoire["children"])):
		clone_number = len(repertoire["children"])
	if ("children" in repertoire) :
		for i in range(0,clone_number) :
			if ("children" in repertoire["children"][i]) :
				clone_tree.append({"name": repertoire["children"][i]["name"], "idV": repertoire["children"][i]["idV"], "idJ": repertoire["children"][i]["idJ"], "ab": repertoire["children"][i]["ab"], "tree_type": tree_type, "children": []})
				clones.append(i)
			else :
				clone_tree.append({"name": repertoire["children"][i]["name"], "idV": repertoire["children"][i]["idV"], "idJ": repertoire["children"][i]["idJ"], "ab": repertoire["children"][i]["ab"]})
	return clone_tree, clones
	
# ------------------------------------------------------------------------------------

def get_tree_branches (file) :
	global branches
	branches = {}
	branch_length = ["","",0]
	if (os.path.isfile(file)) :
		with open(file) as f :
			for line in f :
				branch_length = line.split("\n")[0].split(",")
				branches[(branch_length[0],branch_length[1])]=branch_length[2]

# ------------------------------------------------------------------------------------

def get_sequences (file) :
	global sequences
	sequences = {}
	if (os.path.isfile(file)) :
		with open(file) as f :
			for line in f :
				properties = line.split("	")
				if (len(properties)==6) :
					seqId=properties.pop(0)
					sequences[seqId] = properties
					
# ------------------------------------------------------------------------------------

def get_percentage_uniq_seq (file) :
	global clonotype_uniq_seq
	clonotype_uniq_seq = {}
	if (os.path.isfile(file)) :
		#open the file containing the percentage of uniq sequence in each clonotypes
		with open(file) as f :
			for line in f :
				uniq_seqs = line.split("\n")[0].split("	")
				clonotype_uniq_seq[uniq_seqs[0]]={}
				for seq in uniq_seqs[1:]:
					seq_properties = seq.split(";")
					clonotype_uniq_seq[uniq_seqs[0]][seq_properties[0]] = seq_properties[1]

# ------------------------------------------------------------------------------------

def main():

	usage = usage = "python generate_clonotype_files.py -c <number_clone> -o <output> -a <analysis_id> -n <newick_file> -t <tree_complexity> \n"
	parser = OptionParser(usage)
	parser.add_option("-c", "--number_clone", dest="number_clone",  help="number of clone to analyze")
	parser.add_option("-o", "--output", dest="output",  help="path of the output folder")
	parser.add_option("-a", "--analysis_id", dest="analysis_id",  help="id of the analysis or user id")
	parser.add_option("-n", "--newick_file", dest="newick_file",  help="end of newick file")
	parser.add_option("-t", "--tree_complexity", dest="tree_complexity",  help="complexity of the tree (all, simplification1 or simplification2")
	
	(options, args) = parser.parse_args()
	
	if len(sys.argv) != 11:
		parser.error("Incorrect number of arguments")
	
	total_clone_analyzed = int(options.number_clone)
	path = options.output
	analysis_name = options.analysis_id
	read_json_file(path + analysis_name + "_repertoire.json")
	end_newick_file = options.newick_file
	file_name = options.tree_complexity
		
	clone_tree, clone_number = find_clone(total_clone_analyzed, file_name)

	for i in clone_number :
		tree_file = path + analysis_name + "_" + str(i+1) + end_newick_file
		fasta_file = path + analysis_name + "_" + str(i+1) + "_200_selected_seq_uniq.aln.fa"

		if (os.path.isfile(tree_file)) :
			matching_name = {}
			newick_tree = read_newick_file(tree_file)
			clonotype_seq_id = match_sequence_id(path + analysis_name + "_" + str(i+1) + "_200_matching_sequences.txt", i+1)
			get_tree_branches(path + analysis_name + "_" + str(i+1) + "_200.nk.csv")
			get_sequences (path + analysis_name + "_" + str(i+1) + "_200_aligned_regions.txt")
			get_percentage_uniq_seq (path + analysis_name + "_" + str(i+1) + "_200_proportion_uniq_sequences.txt")
				
			tree = change_format(newick_tree[0].descendants[0],i,clonotype_seq_id,0,"germline",matching_name)
			clone_tree[i]["children"] = tree			
			clone_name = clone_tree[i]["name"]
			
			save_json_format(path + analysis_name + "_tree_" + file_name, clone_tree)

			file_name_beggining = analysis_name + "_" + clone_name + "_" + file_name

			save_sequences(path + file_name_beggining +"_sequences", sequences)

	print("  # Generation of clonotype files for visualization : done.")

#===================================================================================
if __name__ == "__main__":
	main()

