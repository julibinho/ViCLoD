#!/usr/bin/env python3
#-*-coding:Utf-8-*-

"""
Author: Lucile Jeusset, Laboratory of Computational and Quantitative Biology, UPMC, Paris.
Email: lucile_jeusset@hotmail.fr

Translate the txt format of the file containing clone whith their clonotypes in json format. It recover the output file provide from MobiLLe analysis.
Three files must be pass in argument : the file to convert in txt format, the txt file containing the total number of reads and a file containning the html colors to be assigned to the clones and clonotypes.

"""

import sys
import os
from optparse import OptionParser
import json

# ------------------------------------------------------------------------------------

def color_list(file):
	colors = []
	with open(file) as f:
		for line in f :
			colors.append(line.split("\n")[0])
	return colors

# ------------------------------------------------------------------------------------

#change the index of the color, stroke and style
def change_index(tab_index,nb_colors,nb_stroke,reverse):
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

#return the number of reads analyzed
def total_number_of_reads(log_file):
	total = 1
	quality = 0
	abundance = 0
	unannotated = 0
	analyzed = 1
	if (os.path.isfile(log_file)) : 
		with open(log_file) as file:
			lines = file.readlines()
		total = int((lines[0].split("\n")[0]).split(": ")[1])
		unannotated = int((lines[1].split("\n")[0]).split(": ")[1])
		quality = int((lines[2].split("\n")[0]).split(": ")[1])
		abundance = int((lines[3].split("\n")[0]).split(": ")[1])
		analyzed = int((lines[4].split("\n")[0]).split(": ")[1])
	reads = { "total" : total, "low_quality" : quality, "low_abundance" : abundance, "unannotated" : unannotated, "analyzed" : analyzed }
	return reads

# ------------------------------------------------------------------------------------

def save_json_format(file_name, data):
	sorted_output_file = open(file_name+".json", "w")
	json.dump(data, sorted_output_file, indent = 4, sort_keys = False)
	sorted_output_file.close()

# ------------------------------------------------------------------------------------

#create a object with all the informations of the repertoire
def recover_all_repertoire_informations(repertoire_file,colors,reads):
	stroke_style = ["none", "5,5", "1,5"]
	clone_index = [0,len(colors)-1,0]
	clones=[]
	num_clone = 1;
	productivity = {"T": "yes","F" : "no", "" : "/"}	#match for the productivity to complete the table of data
	if (os.path.isfile(repertoire_file)) : 
		with open(repertoire_file) as file:
			for line in file :
				line = line.split("\n")
				[clone_name, clone_abundance, nb_reads, clonotype, idV, idJ, cdr3]=line[0].split("	")	#separates the columns
				clone_name = "C"+str(num_clone)
				clonotypes_list = clonotype.split("Clonotype ")
				clonotypes_list = clonotypes_list[1].split(" ")
				cdr3_seq_list = cdr3.split(" ")
				abundance = (int(nb_reads)/reads["analyzed"])*100
				#save the abundance of the clone with 3 decimals or in scientific format if it's inferior to 0.001
				if(abundance<0.001):
					ab = "%.2e"%abundance
				else :
					ab = round(abundance,3)

				#the clone doesn't have clonotypes
				if(len(clonotypes_list)==1):
					cdr3_seq = cdr3_seq_list[0].split(",")
					clones.append({"name" : clone_name, "value" : abundance, "reads" : nb_reads, "ab" : ab, "idV" : idV, "idJ" : idJ, "cdr3" : cdr3_seq[0], "productivity" : productivity[clonotypes_list[0].split(",")[1]], "color" : colors[clone_index[0]], "stroke" : colors[clone_index[1]], "style" : stroke_style[clone_index[2]], "seq" : cdr3_seq[1]})
					clone_index = change_index(clone_index,len(colors),len(stroke_style),False)	#update the values of the index
				#the clone have clonotype
				else:
					clonotypes = all_clonotypes_in_clone(clone_name, clonotypes_list, cdr3_seq_list, colors, stroke_style, abundance, clone_abundance, clone_index, nb_reads)

					clones.append({"name" : clone_name, "value" : abundance, "reads" : nb_reads, "ab" : ab, "idV" : idV, "idJ" : idJ, "cdr3" : cdr3_seq_list[0].split(",")[0], "color" : colors[clone_index[0]], "stroke" : colors[clone_index[1]], "style" : stroke_style[clone_index[2]], "children" : clonotypes})
					clone_index = change_index(clone_index,len(colors),len(stroke_style),False)	#update the values of the index
				num_clone += 1

	repertoire = {"name" : "repertoire", "value" : 100, "reads": {"total": reads["total"], "analyzed": reads["analyzed"], "low_quality": reads["low_quality"], "low_abundance": reads["low_abundance"], "unannotated": reads["unannotated"]}, "color":"#808080", "stroke":"#808080", "style":"none", "children": clones}	#contain the informations of the repertoire (clones, clonotypes)
	return repertoire

# ------------------------------------------------------------------------------------

#recover all the informations of the clonotypes of one clone
def all_clonotypes_in_clone(clone_name,clonotypes_list,cdr3_seq_list,colors,stroke_style,abundance,clone_abundance,clone_index,clone_reads):
	productivity = {"T": "yes","F" : "no", "" : "/"}	#match for the productivity to complete the table of data
	clonotype_index = [len(colors)-1,clone_index[0],0]
	clonotypes = []
	for i in range(len(clonotypes_list)):
		#the cdr3 is the same for all the clonotypes
		if(len(cdr3_seq_list)==1):
			cdr3_seq = cdr3_seq_list[0].split(",")
		else:
			cdr3_seq = cdr3_seq_list[i].split(",")

		clonotype_info = clonotypes_list[i].split(",") #contains the abundance of the clonotypes and the productivity
		nb_reads = round(float(clonotype_info[0]) * int(clone_reads))
		#save the abundance of the clonotype in the repertoire with 3 decimals or in scientific format if it's inferior to 0.001
		if(abundance*float(clonotype_info[0])<0.001):
			ab_clonotype_rep = "%.2e"%(float(clone_abundance)*100*float(clonotype_info[0]))
		else :
			ab_clonotype_rep = round(float(clone_abundance)*100*float(clonotype_info[0]),3)
		#save the abundance of the clonotype in the clone with 3 decimals or in scientific format if it's inferior to 0.001
		if(float(clonotype_info[0])*100<0.001):
			ab_clonotype = "%.2e"%(float(clonotype_info[0])*100)
		else :
			ab_clonotype = round(float(clonotype_info[0])*100,3)
		clonotypes.append({"name":clone_name+"-"+str(i+1), "value":float(clonotype_info[0])*100, "reads":nb_reads, "ab":ab_clonotype, "ab_rep":ab_clonotype_rep, "productivity":productivity[clonotype_info[1]], "cdr3":cdr3_seq[0], "color":colors[clonotype_index[0]], "stroke":colors[clonotype_index[1]], "style":stroke_style[clonotype_index[2]], "seq":cdr3_seq[1]})
		clonotype_index = change_index(clonotype_index,len(colors),len(stroke_style),True)	#update the values of the index

	return clonotypes

# ------------------------------------------------------------------------------------

def main():

	usage = usage = "python clone_informations_json.py -r <repertoire_file> -l <log_file> -a <analysis_id> -c <colors_file> \n"
	parser = OptionParser(usage)
	parser.add_option("-r", "--repertoire_file", dest="repertoire_file",  help="txt file containing repertoire inforamtions")
	parser.add_option("-l", "--log_file", dest="log_file",  help="log file containing information on reads")
	parser.add_option("-a", "--analysis_id", dest="analysis_id",  help="id of the analysis or user id")
	parser.add_option("-c", "--colors_file", dest="colors_file",  default="colors.txt", help="txt file containing a list of colors")
		
	(options, args) = parser.parse_args()
	
	if len(sys.argv) < 7:
		parser.error("Incorrect number of arguments")
	
	repertoire_file = options.repertoire_file
	log_file = options.log_file
	analysis = options.analysis_id
	colors_file = options.colors_file

	reads = total_number_of_reads(log_file)
	colors = color_list(colors_file) 
	repertoire = recover_all_repertoire_informations(repertoire_file,colors,reads)

	path = os.path.dirname(repertoire_file)
	file_name = path + "/" + analysis + "_repertoire" 

	save_json_format(file_name, repertoire)

	print("  # JSON file of the repertoire created")

# ------------------------------------------------------------------------------------

if __name__ == "__main__":
	main()
