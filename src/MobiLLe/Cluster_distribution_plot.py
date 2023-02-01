
import sys
import math
import time
import resource
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from collections import Counter
from optparse import OptionParser
from palettable.cartocolors.diverging import TealRose_7
#####################################################################
# read the clustering result
def readClusteringResults(nomFi):
	# read the clustering result
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	nom=""
	Clustering_lables={}
	for l in range(len(lines)-1):
		cluster=lines[l].split("\t")[0].rstrip()
		Seq_nom=lines[l].split("\t")[1].split(" ")
		Seq_nom[-1]= Seq_nom[-1][:-1]
		Clustering_lables[cluster] = Seq_nom
	return Clustering_lables

#####################################################################
def gini(arr):
	## first sort
	sorted_arr = arr.copy()
	sorted_arr.sort()
	n = arr.size
	coef_ = 2. / n
	const_ = (n + 1.) / n
	weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
	return coef_*weighted_sum/(sorted_arr.sum()) - const_

#####################################################################
def read_uniformity_file (file_name) :
	data = {}
	
	f=open(file_name,"r")
	lines=f.readlines()
	f.close()
	
	for l in lines :
		line = l.split("\n")[0]
		columns = line.split("\t")
		data[columns[0]] = columns[1]
	return data

#####################################################################
def get_match_uniformity_abundance(dico_result, uniformity_file) :
	uniformity_data = read_uniformity_file(uniformity_file)
	uniform_abundance = []
	for cluster in dico_result :
		if cluster in uniformity_data :
			uniform_abundance.append((cluster,len(dico_result[cluster]),float(uniformity_data[cluster])))
		else :
			uniform_abundance.append((cluster,len(dico_result[cluster]),0))
	uniform_abundance.sort(key=lambda x: x[1])
	return uniform_abundance

#####################################################################
def Plot(Clustering_lables,Dataset_name,adress_tempo_file,uniformity_data):
	plot_letter = [ '(A)' ,'(B)', '(C)', '(D)']
	file = open(adress_tempo_file+Dataset_name.split(".")[0] +"_cluster_distribution.txt","w")
	fig = plt.gcf()
	fig.subplots_adjust(hspace=0.4)
	fig.set_size_inches((20,12))
	fig.suptitle(Dataset_name.split(".")[0]+"'s repertoire", fontsize=20)
	x=[]

	for key in Clustering_lables.keys():
		x.append(len(Clustering_lables[key]))
		if len(Clustering_lables[key]) == 0:
			print (key)
	X=np.array(sorted(x))
	abundance = [i[1] for i in uniformity_data]
	uniformity = [i[2] for i in uniformity_data]
	num_var = len(abundance)

	ax = fig.add_subplot(221)
	ax.annotate(str(plot_letter[0]), xy=(-0.08,1.1 ), xycoords='axes fraction', fontsize=14,horizontalalignment='right', verticalalignment='top')
	N = len(uniformity_data)	
	if num_var >= 20 :
		begin = num_var-20
	else :
		begin = 0
	
	s = [1]*num_var
	for i in range(begin,num_var) :
		s[i] = (abundance[i]/min(abundance))
		
	colors = np.random.rand(N)
	plt.scatter(abundance,uniformity,s=s,c=colors,cmap=TealRose_7.mpl_colormap,alpha=0.5)
	plt.xlabel('Number of sequence per cluster')
	plt.ylabel('Uniformity')
	plt.xlim(xmax = abundance[-1]+abundance[-1]/3, xmin = -(abundance[-1]/10))

	ax = fig.add_subplot(222)
	ax.annotate(str(plot_letter[1]), xy=(-0.08,1.1 ), xycoords='axes fraction', fontsize=14,horizontalalignment='right', verticalalignment='top')
	l=plt.plot(X,'k.')
	plt.setp(l, markersize=3)
	ax.set_yscale("log", nonpositive='clip')
	plt.xlabel('Clusters')
	plt.ylabel('Number of sequence per cluster')

	ax = fig.add_subplot(223, aspect='equal')
	ax.annotate(str(plot_letter[2]), xy=(-0.08,1.1 ), xycoords='axes fraction', fontsize=14,horizontalalignment='right', verticalalignment='top')
	X_lorenz = X.cumsum() /float( X.sum())
	X_lorenz = np.insert(X_lorenz, 0, 0)
	X_lorenz[0], X_lorenz[-1]
	ax.scatter(np.linspace(0.0, 1.0, num = X_lorenz.size, endpoint=True),X_lorenz, marker='.', color='#009392', s=100)
	ax.plot([0,1], [0,1], color='k')
	ax.text(0.0, 0.9, "Gini coefficient = %.2f" % gini(X),fontsize=10)
	plt.xlabel('Cumulative share of clusters')
	plt.ylabel('Cumulative share of sequences')

	ax = fig.add_subplot(224)
	ax.annotate(str(plot_letter[3]), xy=(-0.08,1.1 ), xycoords='axes fraction', fontsize=14,horizontalalignment='right', verticalalignment='top')
	rev=np.fliplr([X])[0]
	clone_list = list(rev)
	seq_number = sum(clone_list)
	j = 1
	for i in clone_list :
		line_to_write = "Cluster " + str(j) + "\t"+ str(i/float(seq_number)) + "\n"
		file.write(line_to_write)
		j += 1
	file.close()

	if len(X) > 100 :
		axe=list(range(1,101))
		s = [(float(n)/(sum(X))) for n in rev[0:100]]
		plt.bar(axe ,s,color='#2b89c8')
		plt.axhline(y=0.05 ,linewidth=1, color='r')
		plt.xlim(0,100)
	else :
		axe=list(range(1,len(X)+1))
		s = [(float(n)/(sum(X))) for n in rev[0:len(X)]]
		plt.bar(axe ,s,color='#2b89c8')
		plt.axhline(y=0.05 ,linewidth=1, color='r')
		plt.xlim(0,len(X))

	plt.xlabel('Clusters')
	plt.ylabel('Abundance')
	fig.savefig(adress_tempo_file+Dataset_name.split(".")[0]+'_repertoire.png', dpi=300)

	#plt.show()
#####################################################################
def main():
	usage = "usage: %prog -n Dataset_name -c ClusteringFile -a adress_tempo_file -u uniformity_file"
	parser = OptionParser(usage)
	parser.add_option("-n", "--Dataset_name", dest="Dataset_name",
	      help="name of the dataset to analyze")
	parser.add_option("-c", "--ClusteringFile",dest="ClusteringFile",
	      help="read data from ClusteringFile")
	parser.add_option("-a", "--adress_tempo_file", dest="adress_tempo_file",
          help="adress for temporary files")
	parser.add_option("-u", "--uniformity_file",dest="uniformity_file",
	      help="read data from uniformity_file")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 9:
		parser.error("incorrect number of arguments")

	Dataset_name = options.Dataset_name
	ClusteringFile = options.ClusteringFile
	adress_tempo_file = options.adress_tempo_file
	uniformity_file = options.uniformity_file

	Dicoresult=readClusteringResults(ClusteringFile)
	uniformity_data = get_match_uniformity_abundance(Dicoresult, uniformity_file)
	
	Plot(Dicoresult,Dataset_name,adress_tempo_file,uniformity_data)
	#print ("Done!")

#####################################################################
if __name__ == "__main__":
	main()


