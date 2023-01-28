import os
import sys
from optparse import OptionParser
import operator
import collections
import time

#####################################################################
def read_file (nomFi):
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	return lines
#####################################################################
def read_vjunction(lines):
	seq_V_CDR3_J = {}
	len_max_CDR3 = 0
	len_max_J = 0
	for l in lines:
		seq = l.split("\t")
		#seq_V_CDR3_J[seq id] = [functionality,V_id,J_id,CDR3_seq_aa,J_seq_nt]
		seq_V_CDR3_J[seq[0]] = [seq[1].rstrip(),seq[2].rstrip(),seq[3].rstrip(),seq[4].rstrip(),seq[5].rstrip()]
		if len(seq[3].rstrip()) > len_max_CDR3 :
			len_max_CDR3 = len(seq[3].rstrip())
		if len(seq[4].rstrip()) > len_max_J :
			len_max_J = len(seq[4].rstrip())
	return seq_V_CDR3_J, len_max_CDR3,len_max_J
#####################################################################
def readClusteringResults(nomFi):
	lines = read_file (nomFi)
	Clustering_lables = {}
	for l in range(len(lines)):
		Seq_nom = lines[l].split("\t")[1].rstrip().split(" ")
		cluster = lines[l].split("\t")[0].rstrip()
		if cluster in Clustering_lables.keys():
			Clustering_lables[cluster].append(Seq_nom)
		else:
			Clustering_lables[cluster] = Seq_nom
	return Clustering_lables
#####################################################################
def calculate_abundance(list_of_occurance):
	sum_all = 0
	abundance = []
	for clonotype in list_of_occurance :
		occurance = clonotype[0]
		functionality = clonotype[1]
		CDR3 = clonotype[2]
		sum_all += occurance
	for clonotype in list_of_occurance :
		occurance = clonotype[0]
		functionality = clonotype[1]
		CDR3 = clonotype[2]
		abundance.append(("%.5f" % (occurance/float(sum_all)),functionality,CDR3))
	return abundance
#####################################################################
def gather_clone_clonotype_info(Dicoresult,dicoSeq,output_file):
	clusters = read_file (output_file)
	clonotypes = {}
	clonal_info = {}
	for c in Dicoresult.keys():
		clonotypes[c] ={}
		seq = Dicoresult[c]
		number_seq_clone = len(seq)
		clone_abundance = number_seq_clone/float(len(dicoSeq))
		#read_count = int(len(dicoSeq)*clone_abundance)
		clonal_info[c] = [("%.5f" % clone_abundance),number_seq_clone,dicoSeq[seq[0]][1],dicoSeq[seq[0]][2]]
		for s in seq :
			if dicoSeq[s][3] in clonotypes[c].keys() :
				clonotypes[c][dicoSeq[s][3]]. append(s)
			else:
				clonotypes[c][dicoSeq[s][3]] =[s]
				clonal_info[c].append(dicoSeq[s][3])
	for clone in clonotypes.keys() :
		list_clonotypes =[]
		for clonotype in clonotypes[clone].keys():
			list_clonotypes.append((len(clonotypes[clone][clonotype]),dicoSeq[clonotypes[clone][clonotype][0]][0],dicoSeq[clonotypes[clone][clonotype][0]][3]))
		clonal_info[clone].append(calculate_abundance(list_clonotypes))
	return clonal_info,clonotypes
#####################################################################
def write_clonal_info(clonal_info_dico_unsorted,output_file,clonotype_dico):
	file_name = output_file+"_repertoire_two_levels_info.txt"
	filetowrite=open(file_name,"w")
	clonal_info_dico = dict(sorted(clonal_info_dico_unsorted.items(), key=lambda t: t[1][0],reverse=True))
	for clone in clonal_info_dico.keys():
		phrase_2 = ''
		phrase = "Clone number " + str(clone) + "\t" + str(clonal_info_dico[clone][0])+ "\t" +str(clonal_info_dico[clone][1])+ "\t" + "Clonotype" 
		sorted_clonotype = sorted(clonal_info_dico[clone][-1],key=lambda tup: tup[0] ,reverse=True)
		for clonotype in  sorted_clonotype:
			phrase += " " + str(clonotype[0])+","+str(clonotype[1])
			phrase_2+=clonotype[2]+","+clonotype_dico[clone][clonotype[2]][0]+" "
		phrase += "\t" + str(clonal_info_dico[clone][2]) + "\t" + str(clonal_info_dico[clone][3]) + "\t"
		phrase_2 = phrase_2[:-1]
		phrase +=phrase_2
		phrase += "\n"
		filetowrite.write(phrase)
	filetowrite.close()
	return 0

#####################################################################
def add_clone_clonotype_to_seq_info(clonotype_dico,dicoSeq,output_file) :
	file_name = os.path.splitext(output_file)[0]+"_final_clusters_seq_info.txt"
	filetowrite=open(file_name,"w")
	for clone in clonotype_dico.keys():
		clonotype_number = 0
		for clonotype in clonotype_dico[clone].keys():
			for seq in clonotype_dico[clone][clonotype] :
				phrase = str(clone) + "_" + str(clonotype_number) + "\t" + str(seq) + "\t" + dicoSeq[seq][0] + "\t" + dicoSeq[seq][1] + "\t" + dicoSeq[seq][2] + "\t" + dicoSeq[seq][3] + "\t" + dicoSeq[seq][4] + "\n"   
				filetowrite.write(phrase)
			clonotype_number += 1
	filetowrite.close()
	return 0
#####################################################################
def main():
	start_time = time.time()
	usage = "usage: two_level_clonal_info.py -f formatted_IMGT_annotation_output -c ClusteringFile -n repertoire_name"
	parser = OptionParser(usage)
	parser.add_option("-f", "--formatted_IMGT_annotation_output", dest="IMGT_seq_info",
	      help="read data from formatted_IMGT_annotation_output")
	parser.add_option("-c", "--ClusteringFile",dest="ClusteringFile",
	      help="read data from ClusteringFile")
	parser.add_option("-n", "--repertoire_name",dest="repertoire_name",
	      help="repertoire_name")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 7:
		parser.error("incorrect number of arguments")
	
	IMGT_seq_info = options.IMGT_seq_info
	ClusteringFile = options.ClusteringFile
	repertoire_name = options.repertoire_name

	lines = read_file (IMGT_seq_info)
	dico_vjunc, len_max_CDR3,len_max_J =read_vjunction(lines)
	dico_result = readClusteringResults(ClusteringFile)
	clonal_info_dico,clonotype_dico = gather_clone_clonotype_info(dico_result,dico_vjunc,IMGT_seq_info)
	write_clonal_info(clonal_info_dico,repertoire_name,clonotype_dico)
	add_clone_clonotype_to_seq_info(clonotype_dico,dico_vjunc,repertoire_name)
	print("  # Creating input file for repertoire visualization execution time : %s seconds " % (time.time() - start_time))
#####################################################################
if __name__ == "__main__":
	main()
