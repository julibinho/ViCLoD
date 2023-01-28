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


# this file contains some utilty function helpful for the  processing data
# Mainly for reading fasta file, 

import pickle
from concensus_machine import get_con_seq



def save_pkl(object, filename):
    pickle.dump(object,open(filename,'wb'))

def load_pkl(filename):
    return pickle.load(open(filename,'rb'))

def csv2fasta(S,outfile):
    with open(outfile,'w') as W:
        for s in S:
            W.write('>'+s["Id"]+'\n'+s['CDR3']+'\n')

def fasta2dict(infile):
    Data=[]

    R=open(infile).readlines()
    i=0
    id=0
    while i<len(R):
        s={"Id":None,"CDR3":None,"Label":None,"Length":None}
        s["Label"]=R[i].strip('>\n')
        cdr=R[i+1].strip()
        s["CDR3"]=cdr
        s["Length"]=len(cdr)
        s['Id']='ANT-'+str(id+1000)
        Data.append(s)
        id+=1
        i+=2
    return Data

def dic2csv(data):
    with open("testdata.txt",'w') as W:
        for d in data:
            W.write(d["Id"]+'\t'+d['CDR3']+'\t'+str(len(d['CDR3']))+'\t'+d['IGHV']+'\t'+d['IGHJ']+'\t'+d['IGHD']+'\t'+str(d['Identity'])+'\n')

def readCSV(filename):
    Point={"Id":None,"IGHV":None,"IGHD":None,"IGHJ":None,"Length":None,"Identity":None,"CDR3":None}
    data=[]
    with open(filename) as F:
        for l in F:
            line=l.strip().split()
            #print (line)
            Point["Id"]=line[0].strip()
            Point["CDR3"] = line[1].strip()
            Point["Length"] = int(line[2].strip())
            Point["IGHV"]=line[3].strip()
            Point["IGHJ"]=line[4].strip()
            Point["IGHD"] = line[5].strip()
            Point["Identity"]=float(line[6].stripprint())
            data.append(Point)
            Point={"Id":None,"IGHV":None,"IGHD":None,"IGHJ":None,"Length":None,"Identity":None,"CDR3":None}
    return data

def readFasta(filename):
    Point={"Id":None,"IGHV":None,"IGHD":None,"IGHJ":None,"Length":None,"Identity":None,"CDR3":None}
    data=[]; info = ""; seq=""
    with open(filename) as F:
        for l in F:
        	l = l.strip()
        	if l.startswith(">"):
        		
        		if (seq != ""):
        			Point["Id"] = ID
        			Point["CDR3"] = seq.upper()
        			Point["Length"] = len(seq)
        			data.append(Point)    		
        		ID  = l.replace(">", "")
        		seq = ""
        		Point={"Id":None,"IGHV":None,"IGHD":None,"IGHJ":None,"Length":None,"Identity":None,"CDR3":None}
        	else:
        		seq += l
        		
        if (seq != ""):
        	Point["Id"] = ID
        	Point["CDR3"] = seq
        	Point["Length"] = len(seq)
        	data.append(Point)
    #print data
    return data

def printList(L):
    for l in L:
        print (l)

def printClusters(S,Cluster_Labels):
    labels=list(set(Cluster_Labels))
    Clusters={}
    #print len(Cluster_Labels),len(S)
    #a=[]
    for i in range(len(S)):
        cls=Cluster_Labels[i]
        #print i,S[i]
        #print cls
        #a=[]
        if Clusters.__contains__(cls):
            #pass
            a=list(Clusters[cls])
            a.append(S[i])
            #print a
            Clusters.__setitem__(cls,a)
        else:
            Clusters.__setitem__(cls,[S[i]])

    return Clusters

def print_unmerged(unmerged_clusters, filename="fast_unmerged_Tcell2.txt"):

    with open(filename, 'w') as W:
        i = 0
        j = 0
        for c in unmerged_clusters.keys():
            c_clusts = unmerged_clusters[c]

            for cc in c_clusts:
                cons_seq = get_con_seq(cc)
                W.write("cluster_" + str(j)+":"+cons_seq + '\n')
                for l in cc:
                    #l = get_real_seq(unmerged_clusters, ccc['Id'])
                    #W.write(str(j) + '\t' + l["Id"] + "\t" + l["CDR3"] + '\t' + l["IGHV"] + '\t' + l["IGHJ"] + '\t' + l[
                     #   "IGHD"] + '\n')
                    W.write(str(j) + '\t' + l["Id"] + "\t" + l["CDR3"]+'\n')
                    j += 1
def write_cluster(Clusters,filename,dbfile="Sorted_clusters.p"):
    original=load_pkl(dbfile)
    L=0
    with open(filename, "w") as F:
        i=1
        for c in Clusters:
            #if len(c)==1:
            #    continue
            #F.write("Cluster"+str(i)+":"+str(len(c))+'\n')
            for l in c:
                #F.write(str(i)+'\t'+l["Id"]+"\t"+l["CDR3"]+'\t'+ l["IGHV"]+'\t'+l["IGHJ"]+'\t'+l["IGHD"]+'\n')

                seq=l["CDR3"]
                l1=len(seq)

                items=original[l1][seq]
                L=L+len(items)
                for item in items:
                    F.write(str(i) + '\t' + item["Id"] + "\t" + item["CDR3"] + '\n')
                #F.write(str(i)+'\t'+l["Id"] + '\n')
            i+=1
    #print ("#Sequence written in the file:", L,'\n')
def prepare_for_kld(clusters):
    final_cluster=[]

    for c in clusters:
        if len(c)>1:
            final_cluster.append(c)
    return final_cluster

