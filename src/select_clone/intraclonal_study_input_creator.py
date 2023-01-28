import sys
from optparse import OptionParser
import operator
import collections
import pandas
import time
import numpy as np
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

#####################################################################

def read_file (nomFi):
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	return lines

#-------------------------------------------------------------------#

def read_cluster_in_order(repertoire_two_levels_info,nb_clone):
	list_clone=[]
	lines=read_file (repertoire_two_levels_info)
	if(len(lines)<int(nb_clone)):
		nb_clone = len(lines)
	for l in range(0,int(nb_clone)):
		cluster_name = lines[l].split("\t")[0].split(" ")[2]
		list_clone.append(cluster_name)
	return list_clone


def read_clusters_seq_info(final_clusters_seq_info):
	cluster_seq = {}
	lines=read_file (final_clusters_seq_info)
	for l in lines:
		if l.split("\t")[0].split("_")[0] in cluster_seq.keys():
			cluster_seq[l.split("\t")[0].split("_")[0]].append(l)
		else :
			cluster_seq[l.split("\t")[0].split("_")[0]]= [l]
	return cluster_seq
#-------------------------------------------------------------------#

def extract_seq_from_clone_id(repertoire_name,list_clone,cluster_seq):
	for clone in range(len(list_clone)) :
		write_clonaltree_align(cluster_seq[list_clone[clone]],repertoire_name,clone+1)
	return 0




def write_clonaltree_align(list_seq,repertoire_name,clone_name):
	file_name = repertoire_name+"_"+str(clone_name)+"_seq_info.txt"
	filetowrite=open(file_name,"w")
	for seq in list_seq:
			filetowrite.write(seq)
	filetowrite.close()
	return 0



#####################################################################
def main():
	start_time = time.time()
	usage = "usage: intraclonal_study_input_creator.py -c final_clusters_seq_info -r repertoire_two_levels_info -n repertoire_name -s number_of_clone_to_analyze"
	parser = OptionParser(usage)
	parser.add_option("-c", "--final_clusters_seq_info", dest="final_clusters_seq_info",
	      help="read data from final_clusters_seq_info")
	parser.add_option("-r", "--repertoire_two_levels_info",dest="repertoire_two_levels_info",
	      help="read data from repertoire_two_levels_info")
	parser.add_option("-n", "--repertoire_name",dest="repertoire_name",
	      help="repertoire_name")
	parser.add_option("-s", "--number_of_clonotype_to_analyze",dest="number_of_clone_to_analyze",
	      help="the number of clone to analyze, if use all, there will be no selection for the analyzing clonotypes otherwise choose the number of n most abondant clonotype ")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 9:
		parser.error("incorrect number of arguments")
	
	final_clusters_seq_info = options.final_clusters_seq_info
	repertoire_two_levels_info = options.repertoire_two_levels_info
	repertoire_name = options.repertoire_name
	nb_clone =  options.number_of_clone_to_analyze

	list_clone = read_cluster_in_order(repertoire_two_levels_info,nb_clone)
	cluster_seq = read_clusters_seq_info(final_clusters_seq_info)
	extract_seq_from_clone_id(repertoire_name,list_clone,cluster_seq)
	print("  # The intraclonal study input creator execution time : %s seconds " % (time.time() - start_time))

#####################################################################
if __name__ == "__main__":
	main()
