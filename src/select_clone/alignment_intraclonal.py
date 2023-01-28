import os
import sys
from optparse import OptionParser
import operator
import collections
import pandas
import numpy as np
import time
import copy
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
#=============================================================================#

def read_AIRR(nomFi):
	#col_list = ['sequence_id','sequence','sequence_aa','rev_comp','productive','complete_vdj','vj_in_frame','stop_codon','locus','v_call','d_call','j_call','c_call','sequence_alignment','sequence_alignment_aa','germline_alignment','germline_alignment_aa','junction','junction_aa','np1','np1_aa','np2','np2_aa','cdr1','cdr1_aa','cdr2','cdr2_aa','cdr3','cdr3_aa','fwr1','fwr1_aa','fwr2','fwr2_aa','fwr3','fwr3_aa','fwr4','fwr4_aa','v_score','v_identity','v_support','v_cigar','d_score','d_identity','d_support','d_cigar','j_score','j_identity','j_support','j_cigar','c_score','c_identity','c_support','c_cigar','v_sequence_start','v_sequence_end','v_germline_start','v_germline_end','v_alignment_start','v_alignment_end','d_sequence_start','d_sequence_end','d_germline_start','d_germline_end','d_alignment_start','d_alignment_end','j_sequence_start','j_sequence_end','j_germline_start','j_germline_end','j_alignment_start','j_alignment_end','cdr1_start','cdr1_end','cdr2_start','cdr2_end','cdr3_start','cdr3_end','fwr1_start','fwr1_end','fwr2_start','fwr2_end','fwr3_start','fwr3_end','fwr4_start','fwr4_end','v_sequence_alignment','v_sequence_alignment_aa','d_sequence_alignment','d_sequence_alignment_aa','j_sequence_alignment','j_sequence_alignment_aa','c_sequence_alignment','c_sequence_alignment_aa','v_germline_alignment','v_germline_alignment_aa','d_germline_alignment','d_germline_alignment_aa','j_germline_alignment','j_germline_alignment_aa','c_germline_alignment','c_germline_alignment_aa','junction_length','junction_aa_length','np1_length','np2_length','n1_length','n2_length','p3v_length','p5d_length','p3d_length','p5j_length','consensus_count','duplicate_count','cell_id','clone_id','rearrangement_id','repertoire_id','rearrangement_set_id','sequence_analysis_category','d_number','5prime_trimmed_n_nb','3prime_trimmed_n_nb','insertions','deletions','junction_decryption']
	df = pandas.read_csv(nomFi, sep='\t')
	df.replace(np.nan,'', inplace=True)
	#df["whole_seq"] = df["fwr1"].astype(str) + df["cdr1"].astype(str) +df["fwr2"].astype(str)+ df["cdr2"].astype(str)+df["fwr3"].astype(str)+df["cdr3"].astype(str)+df["fwr4"].astype(str)
	#df["germline_seq"] = df["v_germline_alignment"].astype(str) + df["np1"].astype(str) + df["d_germline_alignment"].astype(str) + df["np2"].astype(str) + df["j_germline_alignment"].astype(str)
	#df['germline_seq'] = df['germline_seq'].str.replace('.','')
	df.replace('', np.nan, inplace=True)
	df.dropna(axis=0, how='any', thresh=None, subset=["cdr1_start","cdr2_start","cdr3_start","cdr1_end","cdr2_end","cdr3_end","v_sequence_start","v_identity","j_identity","d_sequence_end","d_sequence_start"], inplace=True)
	#print(df.loc[22,"v_identity"])
	#print(df.loc[1,"sequence"],"     ",df.loc[1,"germline_seq"])
	return df
#=============================================================================#

def read_seq_info(lines,airr_df,nb_clonotype):
#the seq of each clonotypes
#germline : sequence with the highest v_identity 
	clonotype_seq = {}
	uniq_seq_proportion = {}
	list_selected_clonotype = []
	list_representative_seq_clonotype = []
	list_seq_same_v_identity = []
	sequences_matching = {}
	hashtable = {}
	for l in lines:
		seq = l.split("\t")
		#print("seq : ",seq)
		clonotype = seq[0].split("_")[1]
		if clonotype in clonotype_seq.keys():
			clonotype_seq[clonotype].append(seq[1])
		else: 
			clonotype_seq[clonotype] = [seq[1]]
	#print("clonotype_seq",clonotype_seq)
	#major_clonotype = (sorted(clonotype_seq, key=lambda k: len(clonotype_seq[k]), reverse=True)[0])
	#germline = airr_df.loc[airr_df['sequence_id'] == clonotype_seq[major_clonotype][0]]["germline_seq"].values[0]
	if nb_clonotype == 'all':
		list_selected_clonotype = list((sorted(clonotype_seq, key=lambda k: len(clonotype_seq[k]), reverse=True)))
	else :
		list_selected_clonotype = list((sorted(clonotype_seq, key=lambda k: len(clonotype_seq[k]), reverse=True)[0:int(nb_clonotype)]))
	
	
	hashtable = {key: list(value) for key, value in airr_df.groupby('whole_seq')['sequence_id']}
	Dico_seq_abundance = hashtable
	hashtable=collections.OrderedDict(sorted(hashtable.items(), key=lambda k: len(k[1]),reverse=True))
	hashtable_seq_id = {value[0]: len(value) for key, value in hashtable.items()}

	for clonotype in list_selected_clonotype :
		#clonotype_sequences = copy.copy(clonotype_seq[clonotype])
		clonotype_uniq_seq = most_abundant_uniq_sequence_in_clonotype(clonotype_seq[clonotype],hashtable,uniq_seq_proportion)
		if clonotype_uniq_seq!="" :
			uniq_seq_proportion[clonotype_uniq_seq]={"sequences":{},"total":0}
			get_proportion_of_uniq_seq_in_clonotype(clonotype_seq[clonotype],hashtable_seq_id,uniq_seq_proportion,clonotype_uniq_seq)
			list_representative_seq_clonotype.append( (clonotype, clonotype_uniq_seq, float(airr_df.loc[airr_df['sequence_id'] == clonotype_uniq_seq]["v_identity"].values[0])) )
			index_seq = clonotype_seq[clonotype].index(clonotype_uniq_seq)
			#store the old and new id of the representative sequence for the clonotype
			sequences_matching[clonotype_seq[clonotype][0]]=clonotype_seq[clonotype][index_seq]
			clonotype_seq[clonotype].insert(0,clonotype_seq[clonotype].pop(index_seq)) #modify the reference sequence in clonotype_seq
	sorted_based_on_V_identity = sorted(list_representative_seq_clonotype, key=lambda element: element[2])
	max_V_identity = sorted_based_on_V_identity[-1][2]
	list_max_V_identity = [element for element in sorted_based_on_V_identity if element[2] == max_V_identity ]
	for seq in list_max_V_identity:
		seq_to_count = airr_df.loc[airr_df['sequence_id'] == seq[1]]["whole_seq"].values[0]
		list_seq_same_v_identity.append([seq[1],len(Dico_seq_abundance[seq_to_count]),float(airr_df.loc[airr_df['sequence_id'] == seq[1]]["d_identity"].values[0]),float(airr_df.loc[airr_df['sequence_id'] == seq[1]]["j_identity"].values[0])])
	
	#test=[['S7863', 730, 90.0, 10.12], ['S1797', 75, 0.0, 96.08], ['S38571', 23, 0.0, 94.12], ['S4205', 74, 10.0, 96.08], ['S945', 52, 0.0, 96.08], ['S28539', 13, 0.0, 96.08], ['S17660', 13, 0.0, 96.08], ['S9991', 9, 0.0, 96.08]]
	seq_closest_germline_id = select_naive_seq(list_seq_same_v_identity)
	germline = airr_df.loc[airr_df['sequence_id'] == seq_closest_germline_id]["germline_seq"].values[0]
	return clonotype_seq,germline,list_selected_clonotype,seq_closest_germline_id,sequences_matching, uniq_seq_proportion

#=============================================================================#
def select_naive_seq(list_seq_same_v_identity):
	# sort the sequences based on 1) J identity 2) D identity 3) abundance 
	sorted_list = sorted(list_seq_same_v_identity, key=lambda x: (-x[3], -x[2], -x[1]))
	# return the sequence with highest J identity
	return sorted_list[0][0]

#=============================================================================#
def most_abundant_uniq_sequence_in_clonotype(clonotype_seqs,hashtable,uniq_seqs):
	seq_abundant=""
	#uniq_seq_proportion={}
	for sequence, liste_seqs in hashtable.items() :
		#print(sequence)
		if len (clonotype_seqs)>=len(liste_seqs) :
			for seq in liste_seqs :
				if seq!="" and seq in clonotype_seqs :
					seq_abundant=seq
					#uniq_seqs[seq]={"sequences":{seq:len(liste_seqs)},"total":len(liste_seqs)}
					del hashtable[sequence]
					break
		if seq_abundant!="" :
			break
	return seq_abundant

#=============================================================================#
def get_proportion_of_uniq_seq_in_clonotype(clonotype_seqs,hashtable_seq_name,uniq_seqs,clonotype_id):
	uniq_seq_of_clonotype = {}
	for seq_id in clonotype_seqs:
		#sequence = airr_df.loc[airr_df['sequence_id'] == seq_id]["whole_seq"]
		#if(not sequence.empty):
			#if(sequence.values[0] in hashtable):
		if seq_id in hashtable_seq_name :
			uniq_seq_of_clonotype[seq_id] = hashtable_seq_name[seq_id]
			#uniq_seqs[clonotype_id]["sequences"][seq_id] = hashtable_seq_name[seq_id]
			uniq_seqs[clonotype_id]["total"] += hashtable_seq_name[seq_id]
			del hashtable_seq_name[seq_id]
	uniq_seqs[clonotype_id]["sequences"] = {k: v for k, v in sorted(uniq_seq_of_clonotype.items(), key=lambda item: item[1])}

#=============================================================================#
def compare_clonotype_to_germline(region):

	region_compered_to_germ = {}
	region_compered_to_germ['germline'] = region['germline']
	germline_seq = region['germline'][0]
	del region['germline']
	#print("region key",region.keys())
	for key in region.keys():
		aligned_seq = ''
		for nt in range(len(germline_seq)):
				if germline_seq[nt] == region[key][0][nt]:
					aligned_seq+='.'
				else :
					aligned_seq+= region[key][0][nt]
		region_compered_to_germ[key] = [aligned_seq,region[key][1]]

	return region_compered_to_germ

#=============================================================================#

def write_clonotype_align(region,repertoire_name,coressp_dico,total_seq_clone,airr_df):

	file_name = repertoire_name+"_aligned_regions.txt"
	filetowrite=open(file_name,"w")
	#print("region",region)
	cdr_regions = str(region["germline"][1][0])
	seq_regions = str(region["germline"][2][0])
	for i in range(1,len(region["germline"][1])):
		cdr_regions += "-" + str(region["germline"][1][i])
	for i in range(1,len(region["germline"][2])):
		seq_regions += "-" + str(region["germline"][2][i])
	G = "germline" + "\t" + "_" +" \t" + "_" + "\t" + cdr_regions + "\t" + seq_regions + "\t" + region['germline'][0] + "\n"
	filetowrite.write(G)
	for key in region.keys():
		if key != 'germline':
			V_percent = airr_df.loc[airr_df['sequence_id'] == key]["v_identity"].values[0]
			d= airr_df.loc[airr_df['sequence_id'] == key]["d_sequence_alignment"]
			cdr3 = airr_df.loc[airr_df['sequence_id'] == key]["cdr3"]
			percentage = ((len(coressp_dico[key])+1)*100)/float(total_seq_clone)
			if(percentage<0.01):
				percentage = "<0.01"
			else:
				percentage = str("%.2f" %percentage)
			#print(len(coressp_dico[key])+1, float(total_seq_clone),(len(coressp_dico[key])+1)/float(total_seq_clone))
			cdr_regions = str(region[key][1][0])
			seq_regions = str(region[key][2][0])
			for i in range(1,len(region[key][1])):
				cdr_regions += "-" + str(region[key][1][i])
			for i in range(1,len(region[key][2])):
				seq_regions += "-" + str(region[key][2][i])
			clonotype_info = key+"\t"+str(len(coressp_dico[key])+1)+"("+percentage+")"+ "\t"+ str(V_percent) +"\t"+cdr_regions+"\t"+seq_regions+"\t"+region[key][0]+"\n"
			filetowrite.write(clonotype_info)
	filetowrite.close()
	return 0

#=============================================================================#

def write_clonaltree_align(clonotype_seq,list_selected_clonotype,airr_df,repertoire_name,germline):
	coressp_dico = {}
	file_name = repertoire_name+"_selected_seq.fasta"

	filetowrite=open(file_name,"w")

	G = ">germline"+ "\n" +germline+ "\n"
	filetowrite.write(G)
	for clonotype in list_selected_clonotype:
		if len(airr_df.loc[airr_df['sequence_id'] == str(clonotype_seq[clonotype][0])]["whole_seq"].values) != 0 :
			seq = ">"+str(clonotype_seq[clonotype][0])+"@"+str(len(clonotype_seq[clonotype])) + "\n" + airr_df.loc[airr_df['sequence_id'] == str(clonotype_seq[clonotype][0])]["whole_seq"].values[0]+ "\n"
			coressp_dico[clonotype_seq[clonotype][0]] =[]
			for dup in clonotype_seq[clonotype][1:]:
				coressp_dico[clonotype_seq[clonotype][0]].append(dup)
			filetowrite.write(seq)
	filetowrite.close()
	return coressp_dico

#=============================================================================#

def alignment(repertoire_name):
	file_name = repertoire_name+"_selected_seq.fasta"
	#clustalw= "/Users/nikaabdollahi/opt/anaconda3/envs/python3.6/bin/clustalw2"
	#cline = ClustalwCommandline(clustalw, infile=file_name, outfile= "nika.aln")
	muscle_cline = MuscleCommandline(input=file_name,out=os.path.splitext(file_name)[0]+".aln")
	muscle_cline()
	align = AlignIO.read(os.path.splitext(file_name)[0]+".aln", "fasta")
	#print(align)
	count = SeqIO.write(align, os.path.splitext(file_name)[0]+"_uniq.aln.fa", "fasta")
	aligned_seq = [(seq_record.id,seq_record.seq) for seq_record in SeqIO.parse(os.path.splitext(file_name)[0]+"_uniq.aln.fa","fasta")]
	return aligned_seq

#=============================================================================#
def readFastaMul(nomFi):
	"""read the fasta file of input sequences"""	
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()

	seq=""
	nom=""
	lesSeq={}
	for l in lines:
		if l[0] == '>':
			if seq != "":
				lesSeq[nom] = seq
			nom=l[1:-1]
			seq=""
		else:
			seq=seq+l[:-1]
	if seq != "":

		lesSeq[nom.rstrip()] = seq.rstrip()
	return lesSeq

#=============================================================================#


def write_all_aligned(repertoire_name,coressp_dico ):
	file_name = repertoire_name+"_all.aln.fa"
	filetowrite=open(file_name,"w")
	dico_fasta =readFastaMul(repertoire_name+"_selected_seq_uniq.aln.fa")
	for seq in dico_fasta.keys():
		sequence = ">"+str(seq)+ "\n" +dico_fasta[seq]+ "\n"
		filetowrite.write(sequence)
		if seq in coressp_dico.keys():
			for dup in coressp_dico[seq]:
				sequence = ">"+str(dup)+ "\n" +dico_fasta[seq]+ "\n"
				filetowrite.write(sequence)
	filetowrite.close()
	return 0

#=============================================================================#

def write_clonotype_align_seq(airr_df,repertoire_name,aligned_seq,sorted_based_on_V_identity_seq_id):

	region = {}
	for seq in aligned_seq:

		if seq[0].split("@")[0] != 'germline':
			a = airr_df.loc[airr_df['sequence_id'] == seq[0].split("@")[0]]
			if len(a["cdr1_start"].values) != 0 and len(a["d_sequence_start"].values)!= 0:
				region[seq[0].split("@")[0]] = []
				startCDR1 = int(a["cdr1_start"].values[0]) - int(a["v_sequence_start"].values[0])
				endCDR1	= (int(a["cdr1_end"].values[0]) - int(a["v_sequence_start"].values[0])) +1 
				startCDR2 = int(a["cdr2_start"].values[0]) - int(a["v_sequence_start"].values[0])
				endCDR2 = (int(a["cdr2_end"].values[0]) - int(a["v_sequence_start"].values[0]))+1
				startCDR3 = int(float(a["cdr3_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
				startD = int(float(a["d_sequence_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
				endD = (int(float(a["d_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
				endCDR3 = (int(a["cdr3_end"].values[0]) - int(a["v_sequence_start"].values[0]))+1
				startV = 0
				endV = (int(float(a["v_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
				startJ = int(float(a["j_sequence_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
				endJ = (int(float(a["j_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
				list_loc =[startCDR1,endCDR1,startCDR2,endCDR2,startCDR3,startD,endD,endCDR3]
				lim = column_from_residue_number(str(seq[1]),[startCDR1,endCDR1,startCDR2,endCDR2,startCDR3,startD,endD,endCDR3])
				region[seq[0].split("@")[0]] = [str(seq[1]),lim,[startV,endV,startCDR3,endCDR3,startJ,endJ]]
		else :
			a = airr_df.loc[airr_df['sequence_id'] == sorted_based_on_V_identity_seq_id]
			region['germline'] = []
			startCDR1 = int(a["cdr1_start"].values[0]) - int(a["v_sequence_start"].values[0])
			endCDR1	= (int(a["cdr1_end"].values[0]) - int(a["v_sequence_start"].values[0])) +1 
			startCDR2 = int(a["cdr2_start"].values[0]) - int(a["v_sequence_start"].values[0])
			endCDR2 = (int(a["cdr2_end"].values[0]) - int(a["v_sequence_start"].values[0]))+1
			startCDR3 = int(float(a["cdr3_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
			startD = int(float(a["d_sequence_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
			endD = (int(float(a["d_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
			endCDR3 = (int(a["cdr3_end"].values[0]) - int(a["v_sequence_start"].values[0]))+1
			startV = 0
			endV = (int(float(a["v_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
			startJ = int(float(a["j_sequence_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
			endJ = (int(float(a["j_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
			lim2 = column_from_residue_number(str(seq[1]),[startCDR1,endCDR1,startCDR2,endCDR2,startCDR3,startD,endD,endCDR3])
			region[seq[0]] = [str(seq[1]),lim2,[startV,endV,startCDR3,endCDR3,startJ,endJ]]
	return region

#=============================================================================#

def column_from_residue_number(seq,res_no_list):
	#print("seq",seq)
	#print("res_no_list",res_no_list)
	"""
	AAAA AA AAAA ===> [3,6]
	A-AAAA-AAAA---A
	A-AAA/A-A /AAA---A ===>[4, 8]
	i=0, 1, 2, 3, 4
	j=0, 0, 1, 2, 3
	"""
	list_region_seq = []
	pos_without_gap = -1
	gap = 0
	for i in range(len(seq)):
		#print ("i =",i)
		if seq[i] != '-' :
			pos_without_gap +=1
			#print("pos_without_gap =",pos_without_gap)
			if pos_without_gap == res_no_list[0]:
				#print(seq[0:i])
				list_region_seq.append(i)
			elif pos_without_gap ==res_no_list[1]:
				#print(seq[list_region_seq[0]:i])
				list_region_seq.append(i)
			elif pos_without_gap ==res_no_list[2]:
				#print(seq[list_region_seq[0]:i])
				list_region_seq.append(i)
			elif pos_without_gap ==res_no_list[3]:
				#print(seq[list_region_seq[0]:i])
				list_region_seq.append(i)
			elif pos_without_gap ==res_no_list[4]:
				#print(seq[list_region_seq[0]:i])
				#particular case where the start of D and the start of cdr3 are the same
				if(res_no_list[4]==res_no_list[5]):
					list_region_seq.append(i)
				list_region_seq.append(i)
			elif pos_without_gap ==res_no_list[5]:
				#print(seq[list_region_seq[0]:i])
				list_region_seq.append(i)
			elif pos_without_gap ==res_no_list[6]:
				#print(seq[list_region_seq[0]:i])
				#particular case where the end of D and the end of cdr3 are the same
				if(res_no_list[6]==res_no_list[7]):
					list_region_seq.append(i)
				list_region_seq.append(i)
			elif pos_without_gap ==res_no_list[7]:
				#print(seq[list_region_seq[0]:i])
				list_region_seq.append(i)
	#print(seq[list_region_seq[1]:len(seq)])
	#print("list_region_seq",list_region_seq)
	return list_region_seq 

#=============================================================================#

#save the old and new id of the representative seq in a txt file
def write_matching_seq(repertoire_name, matching_seq):
	file_name = repertoire_name+"_matching_sequences.txt"
	filetowrite = open(file_name,"w")
	for seq_ID in matching_seq:
		ids = seq_ID + "	" + matching_seq[seq_ID] + "\n"
		filetowrite.write(ids)
	filetowrite.close()
	return 0

#=============================================================================#

#save the old and new id of the representative seq in a txt file
def write_uniq_seq_proportion(repertoire_name, all_seq_clonotype):
	file_name = repertoire_name+"_proportion_uniq_sequences.txt"
	filetowrite = open(file_name,"w")
	for representative_seq in all_seq_clonotype:
		proportion = representative_seq
		for uniq_seq in all_seq_clonotype[representative_seq]["sequences"]:
			percentage = (all_seq_clonotype[representative_seq]["sequences"][uniq_seq]/all_seq_clonotype[representative_seq]["total"])*100
			proportion += "	" + uniq_seq + ";" + str(percentage)
		proportion += "\n"
		filetowrite.write(proportion)
	filetowrite.close()
	return 0

#=============================================================================#

def write_region(repertoire_name, regions):
	file_name = repertoire_name+"_regions.txt"
	filetowrite = open(file_name,"w")
	seq_regions = "germline" +  "\t"
	for i in range(len(regions["germline"][2])):
		seq_regions += str(regions["germline"][2][i]) + " "
	seq_regions += "\n"
	filetowrite.write(seq_regions)
	for seq in regions:
		seq_regions = seq +  "\t"
		for i in range(len(regions[seq][2])):
			seq_regions += str(regions[seq][2][i]) + " "
		seq_regions += "\n"
		filetowrite.write(seq_regions)
	filetowrite.close()
	return 0


#####################################################################
def main():
	start_time = time.time()
	usage = "usage: alignment_intraconal.py -a AIRR_IMGT_annotation_output -f final_seq_info -n repertoire_name -s number_of_clonotype_to_analyze"
	parser = OptionParser(usage)
	parser.add_option("-a", "--AIRR_IMGT_annotation_output", dest="IMGT_seq_info",
	      help="read data from AIRR_IMGT_annotation_output")
	parser.add_option("-f", "--final_seq_info",dest="final_seq_info",
	      help="read data from final_seq_info")
	parser.add_option("-n", "--repertoire_name",dest="repertoire_name",
	      help="repertoire_name")
	parser.add_option("-s", "--number_of_clonotype_to_analyze",dest="number_of_clonotype_to_analyze",
	      help="the number of clonotype to analyze, if use all, there will be no selection for the analyzing clonotypes otherwise choose the number of n most abondant clonotype ")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 9:
		parser.error("incorrect number of arguments")
	
	IMGT_seq_info = options.IMGT_seq_info
	final_seq_info = options.final_seq_info
	repertoire = options.repertoire_name
	nb_clonotype =  options.number_of_clonotype_to_analyze
	repertoire_name = repertoire+"_"+str(nb_clonotype)
	airr_df = read_AIRR(IMGT_seq_info)
	lines_clonotype_seq = read_file (final_seq_info)

	clonotype_seq,germline,list_selected_clonotype,seq_closest_germline_id, matching_seq, proportion_of_uniq_seq = read_seq_info(lines_clonotype_seq,airr_df,nb_clonotype)
	write_matching_seq(repertoire_name, matching_seq)
	write_uniq_seq_proportion(repertoire_name, proportion_of_uniq_seq)
	coressp_dico = write_clonaltree_align(clonotype_seq,list_selected_clonotype,airr_df,repertoire_name,germline)
	aligned_seq = alignment(repertoire_name)
	write_all_aligned(repertoire_name,coressp_dico)

	#column_from_residue_number("A-AAAA-AAAA---A",[3,6])
	region = write_clonotype_align_seq(airr_df,repertoire_name,aligned_seq,seq_closest_germline_id)
	#write_region(repertoire_name, region)
	#region_compered_to_germ=compare_clonotype_to_germline(region)
	total_seq_clone =len(lines_clonotype_seq)
	write_clonotype_align(region,repertoire_name,coressp_dico,total_seq_clone,airr_df)
	clone=repertoire_name.split("_"+str(nb_clonotype))[0].split("/")[-1]
	print("    The intraclonal alignment step execution time for the ",clone," ", nb_clonotype, "clonotypes : %s seconds " % (time.time() - start_time))
#####################################################################
if __name__ == "__main__":
	main()
