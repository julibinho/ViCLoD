
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
def Plot(Clustering_lables,Dataset_name):
	plot_letter = [ '(A)' ,'(B)', '(C)', '(D)']
	file = open(Dataset_name.split(".")[0] +"_cluster_distribution.txt","w")
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


	ax = fig.add_subplot(221)
	ax.annotate(str(plot_letter[0]), xy=(-0.08,1.1 ), xycoords='axes fraction', fontsize=14,horizontalalignment='right', verticalalignment='top')
	N = len(x)
	y = [0]*len(x)
	s = [(n/min(X)) for n in X]
	colors = np.random.rand(N)
	plt.scatter(X,y,s=s,c=colors,cmap=TealRose_7.mpl_colormap,alpha=0.5)
	ax.get_yaxis().set_visible(False) 
	plt.axis('equal')
	plt.xlabel('Number of sequence per cluster')
	plt.xlim(xmax = X[-1]+X[-1]/3, xmin = -(X[-1]/10))

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
	ax.scatter(np.arange(0,1,1/X_lorenz.size),X_lorenz, marker='.', color='#009392', s=100)
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
	fig.savefig(Dataset_name.split(".")[0]+'_repertoire.png', dpi=300)

	#plt.show()
#####################################################################
def main():
	usage = "usage: %prog -n Dataset_name -c ClusteringFile"
	parser = OptionParser(usage)
	parser.add_option("-n", "--Dataset_name", dest="Dataset_name",
	      help="name of the dataset to analyze")
	parser.add_option("-c", "--ClusteringFile",dest="ClusteringFile",
	      help="read data from ClusteringFile")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 5:
		parser.error("incorrect number of arguments")

	Dataset_name = options.Dataset_name
	ClusteringFile = options.ClusteringFile
	Dicoresult=readClusteringResults(ClusteringFile)

	Plot(Dicoresult,Dataset_name)
	#print ("Done!")

#####################################################################
if __name__ == "__main__":
	main()


