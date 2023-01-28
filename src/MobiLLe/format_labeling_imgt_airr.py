import os
import sys
import pandas
import numpy as np
import math
from collections import Counter
from optparse import OptionParser
import time

#####################################################################
def read_AIRR(nomFi):
	df = pandas.read_csv(nomFi, sep='\t')
	return df
	
#####################################################################
def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines

#####################################################################
def filter_dico(Dico,remain_seq) :
	list_all_seq = Dico.keys()
	to_delet = set(list_all_seq) - set(remain_seq)
	for seq in to_delet:
		del Dico[seq]
	return Dico, len(to_delet)
	
#####################################################################
def create_log_file(filename, threshold, all_seq, low_quality, deleted_seq):
	with open(filename, "w") as f:
		f.write("Total sequence count : "+str(all_seq)+"\n")
		f.write("Number of unannotated sequences that have been eliminated from the analysis : 0\n")
		f.write("Number of low quality sequences that have been eliminated from the analysis : "+str(low_quality)+"\n")
		f.write("Number of low abundant sequences that have been eliminated from the analysis, using the threshold of "+str(threshold)+": "+str(deleted_seq)+"\n")
		f.write("Number of sequences to be analysed : "+str(all_seq-low_quality-deleted_seq)+"\n")

#####################################################################

def read_log_file(filename) :
	if os.path.isfile(filename) :
		with open(filename,'r') as f:
			lines = f.readlines()
		if (len(lines)==2) :
			line1 = lines[0].split(":")
			line2 = lines[1].split(":")
			if (len(line1)==2 and len(line2)==2) :
				total = line1[1].strip()
				unnanotated = line2[1].strip()
				if(total.isdigit() and unnanotated.isdigit()) :
					return {"total": int(total), "unnanotated" : int(unnanotated)}
	return None
						
#####################################################################

def modify_log_file(filename, all_seq, threshold, low_quality, deleted_seq) :
	
	log_info = read_log_file(filename)
	if log_info == None :
		create_log_file(filename, threshold, all_seq, low_quality, deleted_seq)
	else :
		seq_analysed = log_info["total"]-(log_info["unnanotated"]+low_quality+deleted_seq)
		with open(filename,'a') as f:
			f.write("Number of low quality sequences that have been eliminated from the analysis : "+str(low_quality)+"\n")
			f.write("Number of low abundant sequences that have been eliminated from the analysis, using the threshold of "+str(threshold)+": "+str(deleted_seq)+"\n")
			f.write("Number of sequences to be analysed : "+str(seq_analysed )+"\n")

#####################################################################
def dico_V_J_CDR3_format(AIRR,threshold_percentage,output_file,pattern):
	df = read_AIRR(AIRR)
	sequences_number = len(df)
	abundance_threshold = round((sequences_number*threshold_percentage)/100)
	if(abundance_threshold<1):
		abundance_threshold = 1
	low_quality_seq, remain_seq =[], []
	Dico={}
	Dico_uniq_adundance = {}
	for i in df.index:
		seq_alignment = df["sequence_alignment"][i]
		if pandas.isnull(seq_alignment):
			seq_alignment = ""
		#the sequence alignment, is better than the sequence because it starts at the beginig of the V and ends at the end of J 
		if 'n' in seq_alignment.lower():
			low_quality_seq.append(i)
		else :
			functionality,V,J,CDR3,Jseq = "_","_","_","_","_"
			sequence_id = df["sequence_id"][i]
			if not pandas.isnull(df["productive"][i]):
				functionality = df["productive"][i]
			if not pandas.isnull(df["v_call"][i]):
				V = df["v_call"][i]
			if not pandas.isnull(df["j_call"][i]):
				J = df["j_call"][i]
			if not pandas.isnull(df["cdr3_aa"][i]):
				CDR3 = df["cdr3_aa"][i].replace("#", ".")
			if not pandas.isnull(df["j_sequence_alignment"][i]):
				Jseq = df["j_sequence_alignment"][i]
			Dico[sequence_id] = [functionality,V,J,CDR3,Jseq]
			
			# keep the sequences based on the read abundace
			
			if seq_alignment in Dico_uniq_adundance.keys() :
				Dico_uniq_adundance[seq_alignment].append(sequence_id)
			else : 
				Dico_uniq_adundance[seq_alignment] = [sequence_id]	
		
	for key in Dico_uniq_adundance.keys():
		if len(Dico_uniq_adundance[key]) > abundance_threshold : 
			remain_seq.append(Dico_uniq_adundance[key])
	flat_list_remain_seq = [item for sublist in remain_seq for item in sublist]
	#delete low abundant unique sequences 
	final_dico, deleted = filter_dico(Dico,flat_list_remain_seq)
	number_of_analyzed_seq = len(Dico.keys()) - deleted
	modify_log_file(output_file+"_log.txt", sequences_number, abundance_threshold, len(low_quality_seq), deleted)

	return final_dico
	'''
	lines = read_output_file(AIRR)
	abundance_threshold = round(((len(lines)-1)*threshold_percentage)/100)
	if(abundance_threshold<1):
		abundance_threshold = 1
	low_quality_seq, remain_seq =[], []
	Dico={}
	Dico_uniq_adundance = {}
	for l in range(1,len(lines)):
		split=lines[l].split("\t")
		#print(l)
		#split[13] is the sequence alignment, is better than the sequence because it starts at the beginig of the V and ends at the end of J 
		if 'n' in split[13].lower():
			low_quality_seq.append(l)
			#print (split[13].lower())
		else :
			functionality,V,J,CDR3,Jseq = "_","_","_","_","_"
			sequence_id = split[0]
			if split[4] != "":
				functionality = split[4]
			if split[9] != "":
				V_columns = split[9].split(" ")
				V = getRegionId(V_columns,pattern['V'])
				#print ("V",V)
			if split[11] != "":
				J_column = split[11].split(" ")
				J = getRegionId(J_column, pattern['J'])
				#print ("J",J)
			if split[28] != "":
				CDR3 = split[28].replace("#", ".")
			if split[89] != "":
				Jseq = split[89]
				#print (CDR3)
			Dico[sequence_id] = [functionality,V,J,CDR3,Jseq]
			print(Dico[sequence_id])
			# keep the sequences based on the read abundace
			if split[13] in Dico_uniq_adundance.keys() :
				Dico_uniq_adundance[split[13]].append(sequence_id)
			else : 
				Dico_uniq_adundance[split[13]] = [sequence_id]
	for key in Dico_uniq_adundance.keys():
		if len(Dico_uniq_adundance[key]) > abundance_threshold : 
			remain_seq.append(Dico_uniq_adundance[key])
	flat_list_remain_seq = [item for sublist in remain_seq for item in sublist]
	#delete low abundant unique sequences 
	final_dico, deleted = filter_dico(Dico,flat_list_remain_seq)
	number_of_analyzed_seq = len(Dico.keys()) - deleted
	create_log_file(output_file+"_log.txt", abundance_threshold, len(lines)-1, len(low_quality_seq), deleted)

	return final_dico
	'''

#####################################################################				
def write_file(Dico_VJCDR3,output_file):
	outputname = os.path.splitext(output_file)[0]+"_V_CDR3_Jseq.txt"
	f = open(outputname,"w")
	for key in Dico_VJCDR3.keys():
		line = key + "\t" + Dico_VJCDR3[key][0] + "\t" + Dico_VJCDR3[key][1] + "\t" + Dico_VJCDR3[key][2] + "\t" + Dico_VJCDR3[key][3] + "\t" + Dico_VJCDR3[key][4] + "\n"
		f.write(line)
	f.close()
	return 0


####################################################################
def main():
	start_time = time.time()
	usage = "usage: format_labeling_imgt.py -a AIRR -o output_file -t threshold_percentage"
	parser = OptionParser(usage)
	parser.add_option("-a", "--AIRR", dest="AIRR",
	      help="read data from AIRR")
	parser.add_option("-o", "--output_file",dest="output_file",
	      help="write data to output_file")
	parser.add_option("-t", "--threshold_percentage",dest="threshold_percentage",
	      help="the threshold in percentage of unique sequences to be analyzed, default = 0.005")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 7:
		parser.error("incorrect number of arguments")
	
	AIRR_file = options.AIRR
	output_file = options.output_file
	threshold_percentage = options.threshold_percentage
	output_file_name = output_file+"_seq_Fo.txt"
	patterns = {"V":"IGHV","J":"IGHJ"}
	time_start = time.perf_counter()
	Dico_VJCDR3 = dico_V_J_CDR3_format(AIRR_file,float(threshold_percentage),output_file, patterns)
	write_file(Dico_VJCDR3,output_file_name)
	print("  # The AIRR file reading step execution time : %s seconds " % (time.time() - start_time))


#####################################################################
if __name__ == "__main__":
	main()
