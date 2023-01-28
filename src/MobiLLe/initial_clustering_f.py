"""
The program gather all sequences with the same V gene and allele, J gene and
Junction AA with an identity percentage value higher than Z into one clone.
"""

import sys
from optparse import OptionParser
import operator
from Levenshtein import distance as levenshtein_distance
import skbio 
from skbio import TabularMSA
from ClusterCDR3Only import FaIR_CDR3Only
import time

# =============================================================================
#               	Read IMGT high/vquest output 
# =============================================================================	

def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines

def filter_seq(lines,output_file):
	empty_cell = "\t "
	file_name = output_file + "_unannotated_seq.tsv"
	filetowrite=open(file_name,"a")
	filtered_lines = []
	number_unannotated = 0
	for l in range(0,len(lines)):
		line = lines[l].split("\n")
		seq= line[0].split("\t")
		if seq[2].rstrip() == "_" or seq[3].rstrip() == "_" or seq[4].rstrip() == "_" :
			seq_unannotated =  seq[0]+  empty_cell + "\t" + seq[1] + "\t" + seq[2] + "\t" + seq[3] + empty_cell*8 + "\t" + seq[4] + empty_cell*29 + "\t" + seq[5] + empty_cell*5 +"\n"
			number_unannotated += 1
			filetowrite.write(seq_unannotated)
		else :
			filtered_lines.append(lines[l])
	filetowrite.close()
	return filtered_lines, number_unannotated

def delete_duplicate(lines):
	uniq_seq_dico = {}
	filtered_lines =[]
	dup_corresp = {}

	for l in range(0,len(lines)):

		seq = lines[l].split("\t")

		dup_id = seq[2].rstrip()+"_"+seq[3].rstrip()+"_"+seq[4].rstrip()

		if dup_id in  dup_corresp.keys() :

			dup_corresp[dup_id].append(seq[0])
		else :
			dup_corresp[dup_id] = [seq[0]]
			filtered_lines.append(lines[l])

	for key in dup_corresp.keys():
		uniq_seq_dico[dup_corresp[key][0]] = dup_corresp[key][1:]
	return uniq_seq_dico,filtered_lines


#=============================================================================#

def group_same_VJ(lines,formatted_file):

	dico_same_VJ = {} # same Vgene, same J gene
	dicoSeq = {}

	for l in range(0,len(lines)):
		NumClone= lines[l].split("\t")
		dicoSeq[NumClone[0]] = [NumClone[1].rstrip(),NumClone[2].rstrip(),NumClone[3].rstrip(),NumClone[4].rstrip()]
		Clone_identity = ""
		Clone_identity = str(NumClone[2].split("*")[0]+"_"+NumClone[3].split("*")[0]) # same V gene  + same J gene
		if Clone_identity in dico_same_VJ.keys():
			dico_same_VJ[Clone_identity].append(NumClone[0])
		else:
			dico_same_VJ[Clone_identity] = [NumClone[0]]
	return dico_same_VJ,dicoSeq

# =============================================================================

def hamming_distance(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))
    
# =============================================================================	
def group_clone_VJ_cdr3(dico_same_VJ,dicoSeq,Clone_threshold):
	VJ_ID_diff_CDR3 = {}
	dicoclone_vj_cdr3 = {}
	loclist =[]
	for VJ_id in dico_same_VJ.keys():
		if VJ_id.split("_")[1] == "IGHJ6" : 
			#for J6 we toleart 1 AA of length difference
			loclist =[]
			for i in range(len(dico_same_VJ[VJ_id])):
				VJ_ID_diff_CDR3[VJ_id]={}
				loclist.append({'Id':dico_same_VJ[VJ_id][i] ,'CDR3' : dicoSeq[dico_same_VJ[VJ_id][i]][3] , 'Length' : len(dicoSeq[dico_same_VJ[VJ_id][0]][3]) })
			clusters=FaIR_CDR3Only(loclist, th=Clone_threshold*100,tolerance=1, mth=Clone_threshold*100,split_size=250)
			for c in clusters :
				CDR3_seq=c[0]["CDR3"]
				VJ_ID_diff_CDR3[VJ_id][CDR3_seq] =[]
				for seq in c:
					VJ_ID_diff_CDR3[VJ_id][CDR3_seq].append(seq["Id"])

		else : 
			loclist =[]
			for i in range(len(dico_same_VJ[VJ_id])):
				VJ_ID_diff_CDR3[VJ_id]={}
				loclist.append({'Id':dico_same_VJ[VJ_id][i] ,'CDR3' : dicoSeq[dico_same_VJ[VJ_id][i]][3] , 'Length' : len(dicoSeq[dico_same_VJ[VJ_id][0]][3]) })
			clusters=FaIR_CDR3Only(loclist, th=Clone_threshold*100,tolerance=0, mth=Clone_threshold*100,split_size=250)
			#print(clusters,"nikaika")
			for c in clusters :
				#print(c)
				CDR3_seq=c[0]["CDR3"]
				VJ_ID_diff_CDR3[VJ_id][CDR3_seq] =[]
				for seq in c:
					VJ_ID_diff_CDR3[VJ_id][CDR3_seq].append(seq["Id"])
	#print(VJ_ID_diff_CDR3)
	return VJ_ID_diff_CDR3



def write_clone_VJ_cdr3_all(VJ_ID_diff_CDR3,dicoSeq,output_file,Clone_threshold):
	file_name = output_file+"_sameVJ_noallele_CDR3_"+ Clone_threshold +".txt"
	filetowrite=open(file_name,"w")
	clone_number = 0
	for VJ in VJ_ID_diff_CDR3.keys():
		#print (VJ,"VJ")
		for cdr3 in VJ_ID_diff_CDR3[VJ]:
			#print (cdr3,"cdr3")
			for seq in VJ_ID_diff_CDR3[VJ][cdr3] :
				#print (dicoSeq[seq.rstrip()],"dicoSeq")
				sequence = str(clone_number) + "\t" + str(seq) + "\t" + dicoSeq[seq.rstrip()][0] + "\t" +dicoSeq[seq.rstrip()][1] +"\t" +dicoSeq[seq.rstrip()][2]+ '\t' + dicoSeq[seq.rstrip()][3]+"\n"
				filetowrite.write(sequence)
			clone_number += 1
	filetowrite.close()
	return 0

def find_max_list(list):
    list_len = [len(i) for i in list]
    return(max(list_len))

#=============================================================================#

def modify_log_file(filename, unnanotated_seq):
	with open(filename,'r+') as file:
		lines = file.readlines()
		file.seek(0)
		for line in lines:
			elements = line.split(':')
			if (elements[0].rstrip() == "Number of sequences to be analysed"):
				file.write(elements[0].rstrip() + " : " + str(int(elements[1])-unnanotated_seq) + '\n')
			elif (elements[0].rstrip() == "Number of unannotated sequences that have been eliminated from the analysis"):
				file.write(elements[0].rstrip() + " : " + str(unnanotated_seq) + '\n')
			else : 
				file.write(line)
		file.truncate()


#=============================================================================#

def main():
    start_time = time.time()
    usage = "python  initial_clustering.py  -i <formated IMGT highvquest statistics output> -o <output file name> -s <Clone identity between zero and one> \n "
    parser = OptionParser(usage)
    parser.add_option("-i", "--hv_stat_output", dest="hv_stat_output",
          help="formated IMGT highvquest statistics output")
    parser.add_option("-o", "--output_file_name", dest="output_file_name",
          help="the name for the file to write")
    parser.add_option("-s", "--Clone_threshold", dest="Clone_threshold",
          help="Clone identity percentage between 0 and 1")
    (options, args) = parser.parse_args()
    if len(sys.argv) != 7:
        parser.error("incorrect number of arguments")

    IMGT_file = options.hv_stat_output
    output_file_name = options.output_file_name
    Clone_threshold = options.Clone_threshold


    CloneID = read_output_file(IMGT_file)
    annotated_CloneID, number_unannotated = filter_seq(CloneID,output_file_name)
    uniq_seq_dico,filtered_lines = delete_duplicate(annotated_CloneID)


    dico_same_VJ,dicoSeq = group_same_VJ(annotated_CloneID,output_file_name)

    
    VJ_ID_diff_CDR3 = group_clone_VJ_cdr3(dico_same_VJ,dicoSeq,float(Clone_threshold))

    write_clone_VJ_cdr3_all(VJ_ID_diff_CDR3,dicoSeq,output_file_name ,Clone_threshold)
    if number_unannotated > 0 :
    	modify_log_file(output_file_name+"_log.txt", number_unannotated)
    	
    print("  # The initial clustering step execution time : %s seconds " % (time.time() - start_time))


#=============================================================================#

if __name__ == "__main__":
    main()
