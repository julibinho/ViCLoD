
#Have to install optparse and tqdm

import sys
import time
import copy

from collections import Counter
from optparse import OptionParser
from tkinter import N
from distances import *

#####################################################################
indS = {'Productivity': 0, 'V_label' : 1,  'J_label':  2, 'CDR3_AA': 3, 'J_seq_nuc' : 4, 'V_seq_nuc': 5, 'J_seq_AA': 6, 'V_seq_AA': 7}


#####################################################################
def read_file (nomFi):
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	return lines
#####################################################################
def delete_duplicate(nomFi):
	"""
	Remove duplicates (identical lignes) from sequence information file
	input nomFi: str				sequence information file name	
	output uniq_seq_dico: Dict()	key=sequence ID; values=list of sequence Ids having sequence identical information to the sequence key
	output filtered_lines: list		list of unique sequence information
	"""
	lines = read_file (nomFi)
	uniq_seq_dico = {}
	filtered_lines =[]
	dup_corresp = {}
	for l in range(0,len(lines)):
		seq= lines[l].split("\t")
		dup_id = seq[2].rstrip()+"_"+seq[3].rstrip()+"_"+seq[4].rstrip()
		if dup_id in  dup_corresp.keys() :
			dup_corresp[dup_id].append(seq[0])
		else :
			dup_corresp[dup_id] = [seq[0]]
			filtered_lines.append(lines[l])
	for key in dup_corresp.keys():
		uniq_seq_dico[dup_corresp[key][0]] = dup_corresp[key][1:]
	return uniq_seq_dico, filtered_lines

#####################################################################
def read_vjunction(lines):
	"""
	Load sequence information 
	input lines: list 			List of unique sequence information
	output seq_V_CDR3_J: Dict()	key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 aa sequence, J nuc sequence, V nuc sequence, J aa sequence, V aa sequence)
	output len_max: Tuple(int, int, int)  max length of CDR3, J, and V sequences
	"""
	seq_V_CDR3_J = {}
	len_max_CDR3 = 0
	len_max_J = 0
	len_max_V = 0
	for l in lines:
		seq = l.split("\t")
		#seq_V_CDR3_J[seq id] = [Functionality,V_id,J_id,CDR3_seq_aa,J_seq_nt,V_seq_nt,J_seq_aa,V_seq_aa]
		seq_V_CDR3_J[seq[0]] = [seq[1].rstrip(),seq[2].rstrip(),seq[3].rstrip(),seq[4].rstrip(),seq[5].rstrip(),seq[6].rstrip(),seq[7].rstrip(),seq[8].rstrip()]
		if len(seq[4].rstrip()) > len_max_CDR3 :
			len_max_CDR3 = len(seq[4].rstrip())
		if len(seq[5].rstrip()) > len_max_J :
			len_max_J = len(seq[5].rstrip())
		if len(seq[6].rstrip()) > len_max_V :
			len_max_V = len(seq[6].rstrip())
	return seq_V_CDR3_J, (len_max_CDR3, len_max_J, len_max_V)
	
#####################################################################
# read the clustering result
def readClusteringResults(nomFi):
	"""
	Load intial clustering information 
	input nomFi: list 					intial clustering file name
	output Clustering_labels: Dict()	key=cluster ID, value= list of sequence Ids within the cluster
	"""

	lines = read_file (nomFi)
	Clustering_labels = {}
	for l in range(len(lines)):
		Seq_nom = lines[l].split("\t")[1].rstrip().split(" ")
		cluster = lines[l].split("\t")[0].rstrip()
		if cluster in Clustering_labels.keys():
			Clustering_labels[cluster].append(Seq_nom)
		else:
			Clustering_labels[cluster] = Seq_nom
	
	return Clustering_labels
	
#####################################################################
def filter_clustering_output (Clustering_labels, uniq_seq_dico):
	"""
	Filter duplicates (identical sequences) from clustering output
	input Clustering_labels: 	Dict()		key=cluster label; values=list of sequence Ids within the cluster
	input uniq_seq_dico:     	Dict()		key=sequence ID; values=list of sequence Ids having sequence identical to the sequence key
	output filtered_clus_label: Dict()		key=cluster label; values=list of unique sequence Ids within the cluster
	"""
		
	filtered_clus_label = {}
	
	for cluster in Clustering_labels.keys():
		
		filtered_clus_label[cluster] = []
		for seq in Clustering_labels[cluster]:
			if seq in uniq_seq_dico.keys():
				filtered_clus_label[cluster].append(seq)
	
	return filtered_clus_label


#####################################################################
def CalculateMedoid(dico_vjunc, Dicoresult):
	"""
	Compute de medoid CDR3 sequence of each cluster
	input dico_vjunc: Dict()	key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicoresult: Dict()	key=cluster label, value= list of unique sequence IDs within the cluster
	output centroid:  Dict()	key=cluster label, value= CDR3 aa sequence representing the medoid
	"""
	
	centroid={}
	for key in Dicoresult.keys() :
		listloc=[]
		for seq in Dicoresult[key]:
			
			if seq.rstrip() in dico_vjunc.keys():
				listloc.append(dico_vjunc[seq.rstrip()][indS['CDR3_AA']])
		if len(listloc) != 0 :
			centroid[key]=Levenshtein.median(listloc)
	return centroid
#####################################################################
def Creat_dico_neighbour(Dicocentroid):
	"""
	Compute the closeness between pair of clusters, finding the nearest neighbor 
	input Dicocentroid:   Dict()	 key=cluster label, value= CDR3 aa sequence representing the medoid
	output dicoNeighbour: Dict()	 key=cluster label, value= cluster label of its closest neighbor
	"""
	Cluster_list=[]
	Centroid_list=[]
	dicoNeighbour={}
		
	for cluster in Dicocentroid.keys():
		Cluster_list.append(cluster)
		Centroid_list.append(Dicocentroid[cluster])
	
	for i in range(len(Centroid_list)) :	
		minDist=1000; closer = ""
		for j in range(len(Centroid_list)):
			if i != j:
				distN = Levenshtein.distance(Centroid_list[i], Centroid_list[j])
				if distN < minDist:
					minDist = distN
					closer = Cluster_list[j]
				
		dicoNeighbour[Cluster_list[i]]= str(closer)
		
	return dicoNeighbour
#####################################################################
def mergeCluster(Dicoresult, inconsistence, Dicofasta, max_lens):
	"""
	Merge clusters when there is inconsistences
	input Dicoresult: 		Dict() 		key=cluster label, value= list of unique sequence IDs within the cluster
	input Dicofasta: 		Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input inconsistence:   	Dict()	 	key=cluster labels of clusters to be mergerd, value= number of inconsistences
	input max_lens:			Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences		
	output DicoCluster: 	Dict()	 	key=cluster label, value= list of unique sequence IDs within the cluster
	"""

	merged = {}
	#CDRBefore = computeCDR(Dicofasta, Dicoresult, len_max_CDR3, len_max_J); print ("CDRBefore : ", CDRBefore )
	DicoCluster = copy.deepcopy(Dicoresult)
	
	for k,v in inconsistence.items():
		cluster1 = k.split('-')[0]; cluster2 = k.split('-')[1]
		#verify if merge will improve CDR
		if cluster1 not in merged.keys() and cluster2 not in merged.keys():
			for seq in DicoCluster[cluster1]:
				DicoCluster[cluster2].append(seq)
			del DicoCluster[cluster1]
			merged[cluster1] = True; merged[cluster2] = True; 
	#CDRAfter= computeCDR(Dicofasta, DicoCluster, len_max_CDR3, len_max_J); print ("CDRAfter : ", CDRAfter )
	#if CDRAfter <= CDRBefore:
	#	return DicoCluster

	return DicoCluster
	

#####################################################################
def run_refinement(Dicofasta, Dicocentroid, Dicoresult, DicoNeighbour, max_lens, typeDistanceT, IT):
	"""
	Run refinement until stop condition
	input Dicofasta: 		Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicocentroid: 	Dict() 		key=cluster label, value= CDR3 aa sequence representing the medoid	
	input Dicoresult: 		Dict() 		key=cluster label, value= list of unique sequence IDs within the cluster
	input max_lens:			Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences		
	input typeDistance: 	Tuple(int, int, int, int) type distance of each componnent
	input IT: int			number of max iterations
	output new_dico_res		key=cluster label, value= list of unique sequence IDs within the cluster
	"""
	
	i=0; j=1

	new_dico_res = Dicoresult
	new_dico_neig = DicoNeighbour
	new_dico_cent = Dicocentroid

	while (j !=0 and i <=IT and new_dico_neig != {}) :
		#before = new_dico_res
		incons = findInconsistence(Dicofasta, new_dico_res, new_dico_neig, max_lens, typeDistanceT)
		if incons:
			new_dico_res = mergeCluster(new_dico_res, incons, Dicofasta, max_lens)
			new_dico_cent = CalculateMedoid(Dicofasta, new_dico_res)
			new_dico_neig = Creat_dico_neighbour(new_dico_cent)
			j +=1
		else :
			j = 0
		i +=1
	
	return new_dico_res, new_dico_neig, new_dico_cent

#####################################################################
def mergeByUniformity(Dicofasta, Dicoresult, DicoNeighbour, max_lens, typeDistanceT):
	"""
	Try to merge singletons if it improve cluster uniformity
	input Dicofasta: 		Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicoresult: 		Dict() 		key=cluster label, value= list of unique sequence IDs within the cluster
	input DicoNeighbour: 	Dict()	 	key=cluster label, value= cluster label of its closest neighbor
	input max_lens:			Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences
	input typeDistance: 	Tuple(int, int, int, int) type distance of each componnent
	"""
	
	DicoCluster = copy.deepcopy(Dicoresult)
	DicoDel = {}
	
	for cluster in DicoCluster.keys() :
		seqs = DicoCluster[cluster]
		
		#verify if it will be merged with its closest neigbour
		clusterNeig = DicoNeighbour[cluster]

		if clusterNeig not in DicoDel.keys():
			seqNeig = DicoCluster[clusterNeig]
			uniformityBef = computeUniformity(Dicofasta, seqNeig, max_lens, typeDistanceT)
			
			newElems = seqs + seqNeig
			
			uniformityAft = computeUniformity(Dicofasta, newElems, max_lens, typeDistanceT)

			diff = uniformityAft - 	uniformityBef		
			if uniformityAft <= uniformityBef and uniformityAft < 0.1 and checkConditions(Dicofasta, seqNeig[0], seqs[0]):
								
				Dicoresult[clusterNeig] = newElems
				del Dicoresult[cluster]
				DicoDel[cluster] = 'del'
	
	return Dicoresult

def checkConditions(Dicofasta, seq1, seq2):
	"""
	Check if cluster sequences has the same CDR3 lenght
	input Dicofasta: 		Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input seq1: 			str 		sequence ID
	input seq2: 			str 		sequence ID
	output boolean
	"""
	if len(Dicofasta[seq1][indS['CDR3_AA']].split(" ")[0]) == len(Dicofasta[seq2][indS['CDR3_AA']].split(" ")[0]):
		return True
	return False


#####################################################################
def mergeSinglentons(Dicofasta, Dicoresult, DicoNeighbour, max_lens, typeDistanceT):
	"""
	Try to merge singletons if it improve cluster uniformity
	input Dicofasta: 		Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicoresult: 		Dict() 		key=cluster label, value= list of unique sequence IDs within the cluster
	input dicoNeighbour: 	Dict()	 	key=cluster label, value= cluster label of its closest neighbor
	input max_lens:			Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences
	input typeDistance: 	Tuple(int, int, int, int) type distance of each componnent
	"""
	
	DicoCluster = copy.deepcopy(Dicoresult)
	DicoMerged = {}
	
	for cluster in DicoCluster.keys() :
		seqs = DicoCluster[cluster]
		if cluster in DicoMerged.keys():
			seqs = Dicoresult[cluster]

		if len(seqs) == 1:
			bi, to_move, cluster_move = computeInterClonesDistance(Dicofasta, DicoCluster, DicoNeighbour, cluster, seqs[0], max_lens, typeDistanceT)
			
			seqsMove =  DicoCluster[cluster_move]
			if cluster_move in DicoMerged.keys():
				seqsMove = Dicoresult[cluster_move]

			uniformityBef = computeUniformity(Dicofasta, seqsMove, max_lens, typeDistanceT)
			
			newElems = copy.deepcopy(seqsMove)
			newElems.append(seqs[0])

			uniformityAft = computeUniformity(Dicofasta, newElems, max_lens, typeDistanceT)
			diff = uniformityAft - 	uniformityBef
			
			if diff <= 0.05 and uniformityAft < 0.1 and checkConditions(Dicofasta, seqs[0], to_move):

				Dicoresult[cluster_move].append(seqs[0])
				del Dicoresult[cluster]

				new_dico_cent = CalculateMedoid(Dicofasta, Dicoresult)
				DicoNeighbour = Creat_dico_neighbour(new_dico_cent)

				DicoMerged[cluster_move] = ""

	return Dicoresult


#####################################################################
def findInconsistence(Dicofasta, Dicoresult, DicoNeighbour, max_lens, typeDistanceT):
	"""
	Detect inconsistences between pair of clusters, when inter clonal distance is smaller than intra clonal distances
	input Dicofasta: 		Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicoresult: 		Dict() 		key=cluster label, value= list of unique sequence IDs within the cluster
	input dicoNeighbour: 	Dict()	 	key=cluster label, value= cluster label of its closest neighbor
	input max_lens:			Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences
	output inconsistence	Dict()	 	key=pair of cluster labels, value= number of inconsistences bi <=ai
	
	"""
	summe=0
	

	inconsistence = {}
	for cluster in Dicoresult.keys() :
	
		for seq in Dicoresult[cluster]:
			ai = computeIntraClonalDistance(Dicofasta, Dicoresult, cluster, seq, max_lens, typeDistanceT)
			bi, to_move, cluster_move = computeInterClonesDistance(Dicofasta, Dicoresult, DicoNeighbour, cluster, seq, max_lens, typeDistanceT)
			
			if bi<ai :
				#verify if CDR3 has the same length
				if len(Dicofasta[seq][indS['CDR3_AA']].split(" ")[0]) == len(Dicofasta[to_move][indS['CDR3_AA']].split(" ")[0]):
					#listI = [cluster, cluster_move]; listI.sort(); keyI =  '-'.join(listI)
					keyI = cluster + "-" + cluster_move
					if keyI in inconsistence.keys():
						inconsistence[keyI] = inconsistence[keyI] + 1
					else:
						inconsistence[keyI] = 1
					
	return inconsistence

	


#####################################################################
def computeUniformity(Dicofasta, seqsCluster, max_lens, typeDistanceT):
	"""
	Compute the uniformity of a cluster
	input Dicofasta: 	Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input seqsCluster: 	list		list of sequence IDs of a cluster
	input max_lens:		Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences
	output avgDensity: 	float		the average of local density of cluster sequences
	
	"""	
	
	seq_number = len(seqsCluster)
	if (seq_number>1) :
		local_density = localDensity(Dicofasta, seqsCluster, max_lens, typeDistanceT)
		average_density = sum(local_density)/seq_number
		return average_density
	
	return 0
	
	
#####################################################################
def computeCDR(Dicofasta, Dicoresult, max_lens, typeDistanceT):
	"""
	Contiguous Density Region (CDR) index calculation steps
	input Dicofasta: 	Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicoresult: 	Dict() 		key=cluster label, value= list of unique sequence IDs within the cluster
	input max_lens:		Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences
	output CDR: float				CDR index identify the largest difference between the clusters regarding their respective densities
	
	"""	
	
	total_seq = 0
	uniform_seq = 0
	for cluster in Dicoresult.keys() :
		seq_number = len(Dicoresult[cluster])
		if (seq_number>1) :
			local_density = localDensity(Dicofasta, Dicoresult[cluster], max_lens, typeDistanceT)
			average_density = sum(local_density)/seq_number
			if (average_density>0) :
				uniformity = sum([x-average_density for x in local_density])/average_density
				uniform_seq += (seq_number*uniformity)
		total_seq+=seq_number
	CDR = uniform_seq/total_seq
	
	return CDR
	
#####################################################################
def localDensity(Dicofasta, seqsCluster, max_lens, typeDistanceT):
	"""
	compute the local density of each sequences of a cluster 
	input Dicofasta: 	Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input seqsCluster: 	list		list of sequence IDs of a cluster
	input max_lens:		Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences
	input typeDistance: Tuple(int, int, int, int) type distance of each componnent
	output local_distances: list()	list of the nearest neighbor of each sequence for a cluster
	"""
	local_distances = []
	for seq1 in seqsCluster :
		distance_neigbor_seq = []
		for seq2 in seqsCluster :
			if seq1 != seq2 :
				distance = computeDistance(Dicofasta, seq1, seq2, max_lens, typeDistanceT)
				distance_neigbor_seq.append(distance)
		local_distances.append(min(distance_neigbor_seq))
	return local_distances

#####################################################################
def write_clone_V_cdr3_(Dicoresult,dicoSeq,uniq_seq_dico,output_file):
	file_name = output_file+"_clone_V_CDR3_J.txt"
	filetowrite=open(file_name,"w")
	clone_number = 0
	for cluster in Dicoresult.keys():
		for seq in Dicoresult[cluster]:
				sequence = str(clone_number) + "\t" + str(seq) + "\t" +dicoSeq[seq][0] + "\t" + dicoSeq[seq][1] + "\t" +dicoSeq[seq][2].split(" ")[0] +"\t" +dicoSeq[seq][3]+"\t" +dicoSeq[seq][4]+ "\n"
				filetowrite.write(sequence)
				if len(uniq_seq_dico[seq]) != 0 :
					for dup in uniq_seq_dico[seq] :
						sequence = str(clone_number) + "\t" + str(dup) + "\t" +dicoSeq[seq][0] + "\t" + dicoSeq[seq][1] + "\t" +dicoSeq[seq][2].split(" ")[0] +"\t" +dicoSeq[seq][3]+"\t" +dicoSeq[seq][4]+ "\n"
						filetowrite.write(sequence)
		clone_number += 1

	filetowrite.close()
	return 0
####################################################################

def main():
	start_time = time.time()
	usage = "usage: refinement.py -i InputFile -l ClusteringFile -v VDistance -j JDistance -c CDR3Distance -d CombineDistance -m MergeSingleton -o outputname"
	parser = OptionParser(usage)
	parser.add_option("-i", "--InputFile", dest="InputFile", help="file containing V, J and CDR3 informations for each sequence")
	parser.add_option("-l", "--ClusteringFile",dest="ClusteringFile", help="read data from ClusteringFile")
	parser.add_option("-v", "--VDistance",dest="VDistance", help="method for calculating the V gene distance, 1- binaire, 2-levenstein,3- GIANA, 4-K-mers")
	parser.add_option("-j", "--JDistance",dest="JDistance", help="method for calculating the J gene distance, 1- binaire, 2-levenstein,3- GIANA, 4-K-mers")
	parser.add_option("-c", "--CDR3Distance",dest="CDR3Distance", help="method for calculating the CDR3 region distance, 1- binaire, 2-levenstein,3- GIANA, 4-K-mers")
	parser.add_option("-d", "--CombineDistance",dest="CombineDistance", help="method for combining three distances (1-mean, 2-harmonic mean,  ponderate mean(requires a Tuple with three weights), )")
	parser.add_option("-m", "--MergeSingleton", dest="MergeSingleton", help="Merging the singletons")
	parser.add_option("-o", "--output", dest="outputname", help="output name")


	IT = 10
	
	
	(options, args) = parser.parse_args()
	if len(sys.argv) != 17:
		parser.error("incorrect number of arguments")

	InputFile = options.InputFile
	outputname = options.outputname
	ClusteringFile = options.ClusteringFile
	VDistance = options.VDistance
	JDistance = options.JDistance
	CDR3Distance = options.CDR3Distance
	CombineDistance = options.CombineDistance
	MergeSingleton = options.MergeSingleton

	

	typeDistanceT = (VDistance, CDR3Distance, JDistance, CombineDistance)
	#typeDistanceT = (1, 2, 2, 1)

	uniq_seq_dico, filtered_lines = delete_duplicate(InputFile)
	dico_vjunc, max_lens =read_vjunction(filtered_lines)
	
	(len_max_CDR3, len_max_J, len_max_V)= max_lens
	dico_result = readClusteringResults(ClusteringFile)
	filtered_clustering_label = filter_clustering_output (dico_result, uniq_seq_dico)

	#There is only one cluster
	if(len(filtered_clustering_label)==1):
		write_clone_V_cdr3_(filtered_clustering_label, dico_vjunc, uniq_seq_dico, outputname)
	else:
	
		######
		#Modify the following functions to integrate the type of distance, for V, J and CDR3 (using VDistance,JDistance,CDR3Distance), as well as the way to calculate the combined distance
		######

		dico_centroid = CalculateMedoid(dico_vjunc, filtered_clustering_label)
		dico_neighbour = Creat_dico_neighbour(dico_centroid)
		
		Dicoresult, dico_neighbour, dico_centroid = run_refinement(dico_vjunc, dico_centroid, filtered_clustering_label, dico_neighbour, max_lens, typeDistanceT, IT)
	
		######
		#merge singletons
		######
	
		if MergeSingleton == '1':
			Dicoresult = mergeSinglentons(dico_vjunc, Dicoresult, dico_neighbour, max_lens, typeDistanceT)

		write_clone_V_cdr3_(Dicoresult, dico_vjunc, uniq_seq_dico, outputname)
	print("  The refinement step execution time : %s seconds " % (time.time() - start_time))
#####################################################################
if __name__ == "__main__":
	main()
  


