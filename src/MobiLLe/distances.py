import Levenshtein 
import sys
import re
from GIANA_distance import * 

#####################################################################
indS = {'Productivity': 0, 'V_label' : 1,  'J_label':  2, 'CDR3_AA': 3, 'J_seq_nuc' : 4, 'V_seq_nuc': 5, 'J_seq_AA': 6, 'V_seq_AA': 7}

#####################################################################
def hamming_similarity(s1,s2,l):
	# computes the hamming similarity between two sequences
	# s1 and s2 are sequencs and l is the length of the sequences
	delta = 0.0
	for a1, a2 in zip(s1, s2):
		if a1 == a2:
			delta = delta + 1
	return (delta / l) * 100


#####################################################################
def generateKmers(seq, k=5):
	"""
	Generate all K-mers of a sequence
	input seq   	: str		sequence
	input k		: int		length of K-mers
	output dico 	: Dict()	key=K-mer, value= frenquence
	"""
	dico = {}
	for i in range(len(seq) -k + 1):
		s = seq[i:i+k]
		if s in dico.keys():
			dico[s]= dico[s] + 1
		else:
			dico[s]= 1
	return dico
	
#####################################################################
def fractionalCommonkMerDistance(dicoS1, dicoS2, minConst):
	"""
	Compute the fractional common k-mers distance, the distance is computed based on the minimum count of every k-mer in the two sequences, 
	thus if two sequences are very different, the minimums will all be small
	input dicoS1   	: Dict()	key=K-mer, value= frenquence
	input dicoS2	: Dict()	key=K-mer, value= frenquence
	output minConst	: int 		for normalisation => min(n,m) - k + 1
	REF: Muscle: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics, 5:113.	
	"""

	s = 0
	for kmer in dicoS1.keys():
		if kmer in dicoS2.keys():
			s += 1
	print (s, minConst)
	return 1 - s/float(minConst)
	
#####################################################################
"""def CalculeClusterDistances(Dicofasta, Dicocentroid, clusters, DicoNeighbour, len_max_CDR3, len_max_J, dico_vjunc):
	
	pre-compute all sequence distances
	input Dicofasta:	Dict()	key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicocentroid: Dict() 	key=cluster label, value= CDR3 aa sequence representing the medoid	
	input Dicoresult:   Dict() 	key=cluster label, value= list of unique sequence IDs within the cluster
	input len_max_CDR3: int		Max length of CDR3 sequence 
	input len_max_J:	int		Max length of IGHJ sequence 
	input dico_vjunc		idem to Dicofasta WHY IS IT TWICE 
	output new_dico_res:Dict()	key=cluster label, value= list of unique sequence IDs within the cluster
	
	#print (clusters)
	DicoDistances = {}
	typeDistanceT = (1, 2, 2, 1)
	for idCluster, idSeqs in clusters.items():
		for seq in idSeqs:
			ai = computeIntraClonalDistance(Dicofasta, clusters, idCluster, seq, len_max_CDR3, len_max_J, typeDistanceT)
			bi, to_move, clusterMove = computeInterClonesDistance(Dicofasta, clusters, DicoNeighbour, idCluster, seq, len_max_CDR3, len_max_J, typeDistanceT)
			DicoDistances[seq] = (ai, idCluster, bi, to_move, clusterMove)
			#print (seq, ai, idCluster, bi, to_move)
			
		
	return DicoDistances
"""
#####################################################################
def IGHVDistance(Dicofasta, ID1, ID2, len_max_V,  typeDistance=1):
	"""
	Compute de IGHV distance between two sequences
	input Dicofasta: Dict()	key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input ID1: str		sequence identifier
	input ID2: str		sequence identifier
	input typeDistance: int	1- binaire, 2-levenstein, 3-K-mers, 4-GIANA (IF WE HAVE AA SEQS
	output V_component: float	IGHV distance 
	"""
	k = 5
	typeDistance = int(typeDistance)
	if typeDistance == 1: #Binary distance
		V_component = 0
		# same V, without considering allel
		if Dicofasta[ID1][indS['V_label']].split("*")[0] != Dicofasta[ID2][indS['V_label']].split("*")[0] :
			V_component = 1
	elif typeDistance == 2:
		V_component = float(Levenshtein.distance(Dicofasta[ID1][indS['V_seq_nuc']], Dicofasta[ID2][indS['V_seq_nuc']])) / len_max_V
	elif typeDistance == 3:
		minConst = min(len(Dicofasta[ID1][indS['V_seq_nuc']]), len(Dicofasta[ID2][indS['V_seq_nuc']])) - k + 1
		dicS1 = generateKmers(Dicofasta[ID1][indS['V_seq_nuc']], k)
		dicS2 = generateKmers(Dicofasta[ID2][indS['V_seq_nuc']], k)
		V_component = fractionalCommonkMerDistance(dicS1, dicS2, minConst)
		#print (Dicofasta[ID1][indS['V_seq_nuc']]); print (Dicofasta[ID2][indS['V_seq_nuc']]); print (V_component)
	elif typeDistance == 4 : 
		V_component = GIANA_distance(Dicofasta[ID1][indS['V_seq_AA']], Dicofasta[ID2][indS['V_seq_AA']])
	return V_component


#####################################################################
def CDR3Distance(Dicofasta, ID1, ID2, len_max_CDR3, typeDistance=2):
	"""
	Compute de CDR3 distance between two sequences
	input Dicofasta: Dict()	key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input ID1: str		sequence identifier
	input ID2: str		sequence identifier
	input len_max_CDR3: int	Max length of CDR3 sequence 
	input typeDistance: int	1- binaire, 2-levenstein, 3-K-mers, 4-GIANA 
	output CDR_component: float	CDR3 distance 
	"""
	k = 3
	CDR3Seq1 = Dicofasta[ID1][3].split(" ")[0]; CDR3Seq2 = Dicofasta[ID2][3].split(" ")[0]; 
	typeDistance = int(typeDistance)
	if typeDistance == 1: #Binary distance
		if CDR3Seq1 == CDR3Seq2:
			CDR_component = 0
		else:
			CDR_component = 1
	elif typeDistance == 2: #Normalized Levenshtein distance
		CDR_component = float(Levenshtein.distance(CDR3Seq1, CDR3Seq2)) / len_max_CDR3
	elif typeDistance == 3:
		minConst = min(len(Dicofasta[ID1][indS['CDR3_AA']]), len(Dicofasta[ID2][indS['CDR3_AA']])) - k + 1
		dicS1 = generateKmers(Dicofasta[ID1][indS['CDR3_AA']], k)
		dicS2 = generateKmers(Dicofasta[ID2][indS['CDR3_AA']], k)
		CDR_component = fractionalCommonkMerDistance(dicS1, dicS2, minConst)
	elif typeDistance == 4 : 
		CDR_component = GIANA_distance(Dicofasta[ID1][indS['CDR3_AA']], Dicofasta[ID2][indS['CDR3_AA']])
		
	return CDR_component

#####################################################################
def IGHJDistance(Dicofasta, ID1, ID2, len_max_J, typeDistance=2):
	"""
	Compute de IGHJ distance between two sequences
	input Dicofasta: Dict()	key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input ID1: str		sequence identifier
	input ID2: str		sequence identifier
	input len_max_J: int	Max length of IGHJ sequence 
	input typeDistance: int	1- binaire, 2-levenstein, 3-K-mers, 4-GIANA ???
	output J_component: float	IGHJ distance 
	"""
	k = 5
	typeDistance = int(typeDistance)
	if typeDistance == 1: #Binary distance
		if Dicofasta[ID1][indS['J_seq_nuc']] == Dicofasta[ID2][indS['J_seq_nuc']]:
			J_component = 0
		else:
			J_component = 1
	elif typeDistance == 2: #Normalized Levenshtein distance
		J_component = 1- get_similarity_score(Dicofasta[ID1][indS['J_seq_nuc']], Dicofasta[ID2][indS['J_seq_nuc']], len(Dicofasta[ID1][indS['J_seq_nuc']]), len(Dicofasta[ID2][indS['J_seq_nuc']])) 
	elif typeDistance == 3:
		minConst = min(len(Dicofasta[ID1][indS['J_seq_nuc']]), len(Dicofasta[ID2][indS['J_seq_nuc']])) - k + 1
		dicS1 = generateKmers(Dicofasta[ID1][indS['J_seq_nuc']], k)
		dicS2 = generateKmers(Dicofasta[ID2][indS['J_seq_nuc']], k)
		J_component = fractionalCommonkMerDistance(dicS1, dicS2, minConst)
	elif typeDistance == 4 : 
		J_component = GIANA_distance(Dicofasta[ID1][indS['J_seq_AA']], Dicofasta[ID2][indS['J_seq_AA']])
		
	return J_component

#####################################################################
def computeDistance(Dicofasta, ID1, ID2, max_lens, typeDistanceT):
	"""
	Compute de distance between two sequences
	input Dicofasta: 	Dict()	key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input ID1: 			str		sequence identifier
	input ID2: 			str		sequence identifier
	input max_lens:		Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences 
	input typeDistanceT:Tuple(int, int, int, int) type distance of each componnent
	output dist: 		float	composed distance 
	"""
	dist = 0
	(len_max_CDR3, len_max_J, len_max_V)= max_lens

	tDV, tDCDR3, tDJ, tMean = typeDistanceT 

	V_component = IGHVDistance(Dicofasta, ID1, ID2, len_max_V, tDV)
	CDR_component = CDR3Distance(Dicofasta, ID1, ID2, len_max_CDR3, tDCDR3)
	J_component = IGHJDistance(Dicofasta, ID1, ID2, len_max_J, tDJ)
	
		
	if tMean == '1': 
		dist = (CDR_component + V_component + J_component)/3.0
	else:
		new_tMean = re.sub(r"[^a-zA-Z0-9]","", tMean)
		t = list(new_tMean); #print (t)
		dist = (V_component*float(t[0]) + CDR_component*float(t[1]) + J_component*float(t[2]))/(float(t[0]) + float(t[1]) + float(t[2]))
		
	#elif tMean == '2': 
	#	dist = 3/(pow(CDR_component,-1) + pow(V_component,-1) + pow(J_component,-1))

	return dist
	
#####################################################################
def computeIntraClonalDistance(Dicofasta, Dicoresult, cluster, seqi, max_lens, typeDistanceT):
	"""
	compute intraclonal distance for a sequence i within its cluster
	input Dicofasta: 	Dict()	key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicoresult: 	Dict() key=cluster label, value= list of unique sequence IDs within the cluster
	input cluster: 		str	cluster identifier Ci
	input seqi: 		str		sequence identifier
	input max_lens:		Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences 
	input typeDistanceT:Tuple(int, int, int, int) type distance of each componnent
	output ai: float	the average distance of seqi and all other sequences within the same cluster Ci
	"""
	dist_intra =  0; ai = 0
	
	for seq_same_clust in Dicoresult[cluster] :
		if seq_same_clust != seqi :
			dist_intra += computeDistance(Dicofasta, seqi, seq_same_clust, max_lens, typeDistanceT)
	if len(Dicoresult[cluster]) != 1 :
		ai = float(dist_intra) / (len(Dicoresult[cluster]) - 1)
	
	return ai
	
	
#####################################################################
def computeInterClonesDistance(Dicofasta, Dicoresult, DicoNeighbour, cluster, seqi, max_lens, typeDistanceT):
	"""
	compute interclone distance for a sequence i out its cluster, to simplyfy we consider only the nearest cluster
	input Dicofasta: 		Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicoresult: 		Dict() 		key=cluster label, value= list of unique sequence IDs within the cluster
	input DicoNeighbour: 	Dict() 		key=cluster label, value= cluster label of its closest neighbour
	input cluster: 			str			cluster identifier
	input seqi: 			str			sequence identifier
	input max_lens:			Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences 
	input typeDistance: 	Tuple(int, int, int, int) type distance of each componnent
	output bi:				float		the min distance of seqi and all other sequences of the closest neighbour cluster
	output neighbour_seq:	str			ID of the sequence with minimum distance bi
	output neighbour_clus	str			iD of the cluster with minimum distance bi
	"""

	dist_neighb = {}
	for seq_neighb in Dicoresult[DicoNeighbour[cluster]] :
		dist_neighb[seq_neighb] = computeDistance(Dicofasta, seqi, seq_neighb, max_lens, typeDistanceT)
	bi = min(dist_neighb.values())
	neighbour_seq = (list(dist_neighb.keys())[list(dist_neighb.values()).index(bi)])
	return bi, neighbour_seq, DicoNeighbour[cluster]


#####################################################################
#computes the sequence similairy between s1 with length l1 and s2
# with length l2. Hamiing similarity if l1==l2 and edit distance/levenstein 
# distance if l1!=l2
def get_similarity_score(s1, s2, l1, l2):
	
	l_avg=(l1+l2)/2.0
	if l1==l2:
		sim=hamming_similarity(s1,s2,l1)
	else:
		if l1>l2 :
			long_seq = s1
			short_seq = s2
		else : 
			long_seq = s2
			short_seq = s1
		five_prime = long_seq[ : len(short_seq)]
		three_prime = long_seq[-(len(short_seq)) :]
		sim_five=hamming_similarity(five_prime,short_seq,len(short_seq))
		sim_three=hamming_similarity(three_prime,short_seq,len(short_seq))
		sim_lev=100.0*(1.0-Levenshtein.distance(s1,s2)/l_avg)
		sim = max(round(sim_five),round(sim_three),round(sim_lev))

	#if sim > 100 or sim <0:
	#	print ("WARNING ", s1, s2, l1, l2)
	return sim/100
