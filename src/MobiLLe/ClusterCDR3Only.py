"""
Copyright (c) 2019 Bishnu Sarker (bishnukuet@gmail.com), Nika Abdollahi, Juliana Silva Bernardes

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""




from clusterIO import save_pkl
#from distance_measure import delta_similarity, hamming_similarity
from Levenshtein import distance
from concensus_machine import get_con_seq
import multiprocessing as MP
from clusterIO import print_unmerged


#========Similarity Function======================

def hamming_similarity(s1,s2,l):
    # computes the hamming similarity between two sequences
    # s1 and s2 are sequencs and l is the length of the sequences
    delta = 0.0
    for a1, a2 in zip(s1, s2):
        if a1 == a2:
            delta = delta + 1
    return (delta / l) * 100

def delta_similarity(s1, s2):
    #computes the hamming similarity between the two sequences
    # s1 and s2 are the sequences 
    l1 = len(s1)
    l2 = len(s2)
    delta = 0.0
    if l1 != l2:
        return 0
    for a1, a2 in zip(s1, s2):
        if a1 == a2:
            delta = delta + 1
    return (delta / l1) * 100
def get_similarity_score(s1,s2, l1, l2):
    #computes the sequence similairy between s1 with length l1 and s2
    # with length l2. Hamiing similarity if l1==l2 and edit distance/levenstein 
    # distance if l1!=l2
    #l1=len(s1)
    #l2=len(s2)
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
       
        #d=100.0*(1.0-editdistance.eval(s1,s2)/float(max(l1,l2)))
        #sim=100.0*(1.0-nlevenshtein(s1,s2,method=1))
        sim_lev=100.0*(1.0-distance(s1,s2)/l_avg)
        sim = max(round(sim_five),round(sim_three),round(sim_lev))

    return sim


#=========================================================================
# Given a  list of sequence, it filters out the duplicate and group the 
# same lenghth sequences together.
#=========================================================================
def duplicate_filter(sequences,directory):
    # the variable sequences is a dictionary with an  key "CDR3" that holds 
    # the peptide sequence/protein sequence
    #print("Eliminating Duplicate Sequence\n")
    c=0
    Clusters = {}
    unique_sequences = {}
    duplicate_dico = {}

    for s in sequences:
        aa = s["CDR3"]
        l = len(aa)
        if l not  in unique_sequences:

            #Clusters.__setitem__(l, [s])
            unique_sequences[l]={aa:[s]}

        else:

            length_dic = unique_sequences[l]
            if aa not in length_dic:
                #print ("Unique Sequence Added\n")
                length_dic.update({aa:[s]})
            else:
                #print ("Duplicate Sequence\n")
                #print (s)
                s_list=length_dic[aa]
                s_list.append(s)
                length_dic[aa]=s_list

            unique_sequences[l]=length_dic
        c=c+1
    # it saves the sequence clusters to retrieve the duplicate sequences later 
    # when printing the result.
    adress = directory+"Sorted_clusters.p"
    save_pkl(unique_sequences, adress)

    for key in unique_sequences.keys():
        #print "preparing for Pre-cluster....\n"
        for k in unique_sequences[key].keys():
            duplicate_dico[unique_sequences[key][k][0]['Id']] = unique_sequences[key][k][0:]
        Clusters[key]=[unique_sequences[key][x][0] for x in unique_sequences[key].keys()]
    del unique_sequences
    return Clusters,duplicate_dico
def precluster_pool(data):
    #print (data,"kjsgcqgsmjcjkqsbmk")
    # This function is the worker function for Multiprocessing Pool
    # the primary task is  performing a greedy clustering over the small
    # chunks of data send  to it.
    S,key,th=data
    clusters = []
    nClust = 0
    le=len(S)
    cc=0
    for p in S:

        if clusters == []:
            clusters.append([p])
            nClust += 1
        else:
            l = len(clusters)
            i = 0
            maxSim = -1
            maxClust = -1
            while i < l:
                t = clusters[i]
                con_seq = get_con_seq(t)  # compute the representative sequence of the cluster
                #f, r = average_similarity2(p, t, threshold)
                f= hamming_similarity(p["CDR3"],con_seq, p["Length"])
                if f >= maxSim:
                    maxSim = f
                    maxClust = i
                i = i + 1
            if maxSim <th:
                clusters.append([p])
                nClust += 1
            else:
                g = clusters[maxClust]
                clusters.remove(g)
                g.append(p)
                clusters.append(g)
        cc=cc+1

    return (key,clusters)
   # print "finished", key, len(S)

#===============================================
# this function does the splitting task
# L list of items, key of the final dictionary
# limit is the maximum  block size
#===============================================
def crunch(L, key, limit):
    # Splits the dataL into several chunks of size limit
    l=len(L)
    final_l={}
    key=str(key)+"_"
    splits=int(l/limit)
    left=l-splits*limit
    i=0
    while i<splits:
        j=i*limit
        p=[]
        while j<i*limit+limit:
            p.append(L[j])
            j+=1
        final_l.update({key+str(i):p})
        i+=1
    if  left>0:
        final_l.update( {key+str(i): [L[k] for k in range(i*limit,l)]})
    
    return final_l

#=========================================================================
# this function splits the sequences into smaller block of split_limit
# It helps to process the faster by individual cores
#=========================================================================
def spliting_filter(Clusters, split_limit=3000):
    # Clusters={l1:[s1,s2,s3,.....,sn],l2:[s1,s2,s3,.....,sn]}
    pre_clusters={}
    for key in Clusters:
        if len(Clusters[key])>split_limit:
            pre_clusters.update(crunch(Clusters[key],key,split_limit))
        else:
            pre_clusters.update({str(key)+"_0":Clusters[key]})
    return pre_clusters

def fast_preclustering_pool(pre_groups, th, split_size=250):
    # Performs a multicore  pre_clustering
    cpc=MP.cpu_count()
    process_pool=MP.Pool(cpc)
    splited_group=spliting_filter(pre_groups, split_limit=split_size)
    keys=splited_group.keys()

    pre_clusters=process_pool.map(precluster_pool, [(splited_group[key], key,th) for key in keys ])
    process_pool.close()
    process_pool.join()
    return dict(pre_clusters)

def  super_merge (pre_clusters, duplicate_dico, t=0, mth=90.0):
    # Perform the merging of the pre_clusters
    clusters=[]
    finalclusters = []
    #print("Merging.....")
    for key in pre_clusters.keys():
        #print (key)
        leftout=[]
        if clusters==[]:
            clusters=clusters+pre_clusters[key]
            #print clusters
        else:
            nC=pre_clusters[key]
            for c1 in nC:
                s1=c1[0]['CDR3']
                l1=c1[0]['Length']
                a=[l1+x for x in range(0,t+1)]+[l1-x for x in range(0,t+1)]
                candidates=list(set(a))
                smax=-1
                index=0
                i=0
                L=len(clusters)
                while i<L:
                    #print "i",i
                    c2=clusters[i]
                    s2=c2[0]['CDR3']
                    l2=c2[0]['Length']
                    if not l2 in candidates:
                        i=i+1
                        continue
                    s=get_similarity_score(s1,s2,l1, l2)
                    if s>=smax:
                        smax=s
                        index=i
                    i=i+1
                if smax>=mth:
                    
                    clusters[index]=clusters[index]+c1
                    
                else:
                    leftout.append(c1)
            clusters=clusters+leftout
    for c in clusters:
        listloc=[]
        for s in c :
            for i in duplicate_dico[s['Id']]:
                listloc.append(i)
        finalclusters.append(listloc)
    return finalclusters



#Final Procedure accumulating all the functions for clustering the sequences based on CDR3 only
def FaIR_CDR3Only(sequences, directory,th=60.0,tolerance=1, mth=60.0,split_size=250):
    #print ("Eliminating Duplicate Sequences...") 
    pre_group,duplicate_dico=duplicate_filter(sequences,directory) 
    del sequences 
    keys=pre_group.keys() 
    #print ("Pre Clustering....\n")
    pre_clusters=fast_preclustering_pool(pre_group,th,split_size)
    pre_keys=pre_clusters.keys()
    #print(" Merging is started..\n")
    clusters=super_merge(pre_clusters,duplicate_dico, t=tolerance,mth=mth )
    
    return clusters







