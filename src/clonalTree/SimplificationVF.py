from ete3 import Tree
import BasicSeq
import BasicTree
import numpy as np
import sys
import os
from optparse import OptionParser

#Fonction auxiliaire de détermination des 10 clonotypes les plus abondants
def dix_plus_abondants(Abondance): 
    return(sorted(Abondance.items(),key=lambda x: x[1])[-10:])

#Première Simplification

def arbre_simplifie(NKfile,FastaFile,Seuil_Abondance, nombre_min_noeuds) : #Fonction renvoyant la liste des noms des clonotypes à conserver après simplification

    Abondance=BasicSeq.readFastaAbundance(FastaFile)[3]
    labels = BasicSeq.readFastaAbundance(FastaFile)[0]
    keptSequencesNames=[]
    tot=0 ## On calcule le nombre total de séquences observées pour calculer le pourcentage de représentation de chaque séquence unique
    for i in Abondance :
        tot+=Abondance[i]
    t=BasicTree.readNKTree(NKfile) ##arbre que l'on souhaite simplifier
    for node in t.traverse() :
        if node.name not in Abondance :
            Abondance[node.name]=0
        if node.name not in labels :
            labels.append(node.name)
    if len(Abondance)<nombre_min_noeuds : #si l'arbre de base est plus petit que la taille minimale de l'arbre que l'on veut obtenir
        nombre_min_noeuds=len(Abondance)-1
    L=t.get_leaves()
    nodes_done=[]   #noeuds dont la descendance a déjà été traitée (permet de diminuer la complexité temporelle)
    for i in L :
        p=i.up
        if p.name not in nodes_done :
            nodes_done.append(p.name)
            l=p.get_children()
            if (len(l)>1) : #on a ici plus d'un descendant, donc un effet peigne à considérer
                abund={} #sous dictionnaire contenant seulement les clonotypes à considérer dans cette partie de la boucle
                descendance={} #sous dictionnaire contenant les noms des clonotypes considérés auxquels est associé true si le clonotype a une descendance, false sinon
                for j in l :
                    name=j.name
                    if name !="" :
                        abund[name]=Abondance[name]
                    if j.is_leaf():
                        descendance[name]=False
                    else :
                        descendance[name]=True
                k="" ##nom de l'abondance maximale
                max=0 ##abondance maximale
                ajout=False         #ajout permet de savoir si on a déjà ajouté un représentant de cette lignée
                for j in abund :
                    if abund[j]>max : #calcul du max de l'abondance pour savoir quel clonotype conserver
                        k=j
                        max = abund[j]
                    
                    if (descendance[j]==True and j not in keptSequencesNames ) : 
                        keptSequencesNames.append(j)
                        ajout=True
                        
                    if (abund[j]/tot > Seuil_Abondance and j not in keptSequencesNames) : 
                        ajout=True
                        keptSequencesNames.append(j)
                        
                if (k not in keptSequencesNames and  ajout==False) :
                    keptSequencesNames.append(k)
                    

            else :
                if l[0].name not in keptSequencesNames and l[0].name !='':
                    keptSequencesNames.append(l[0].name)
               
    
    for k in dix_plus_abondants(Abondance) : #on s'assure de conserver les 10 clonotypes les plus abondants dans le résultat final
        if k[0] not in keptSequencesNames :
            keptSequencesNames.append(k[0])
        
    Seq_a_ajouter =[] # on ajoute les ancêtres de tous les clonotypes conservés
    for i in keptSequencesNames : 
        p=t.search_nodes(name=i)[0]
        while p.is_root()==False :
            p=p.up
            if (p.name!="" and p.name not in keptSequencesNames and p.name not in Seq_a_ajouter) :
                Seq_a_ajouter.append(p.name)
    
    keptSequencesNames=keptSequencesNames+Seq_a_ajouter

    i=0
    abondance_triee=sorted(Abondance.items(),key=lambda p: p[1], reverse=True)
    while len(keptSequencesNames)<nombre_min_noeuds : #on complète l'arbre avec les clonotypes les plus abondants
        if i >=len(Abondance) :
            break
        seq_a_ajouter=abondance_triee[i][0]
        node = t.search_nodes(name=seq_a_ajouter)[0]
        while node.is_root()==False :
            if node.name not in keptSequencesNames :
                keptSequencesNames.append(node.name)
            node=node.up
        i+=1
    
    return keptSequencesNames


def simplification_abondance(NKfile,FastaFile,Seuil_Abondance,nombre_min_noeuds): #permet de récupérer une chaîne de caractère correspondant à l'arbre simplifié, au format nk
    t=Tree(NKfile,format=1)
    L=arbre_simplifie(NKfile,FastaFile,Seuil_Abondance,nombre_min_noeuds)
    t.prune(L)
    return (t.write(format=1))


def generation_fasta_simplifie(NKfile,FastaFile,Seuil_Abondance,nombre_min_noeuds ) : #fonction permettant de générer un fichier fasta ne contenant que l'information des séquences que l'on souhaite conserver
    keptSequencesNames=arbre_simplifie(NKfile,FastaFile,Seuil_Abondance,nombre_min_noeuds)
    labels = BasicSeq.readFastaAbundance(FastaFile)[0]
    arraySeqs = BasicSeq.readFastaAbundance(FastaFile)[2]
    Abondance = BasicSeq.readFastaAbundance(FastaFile)[3]
    string =""
    for i,j in enumerate(labels) :
        if j in keptSequencesNames :
            string = string + ">"+str(j)+"@"+str(Abondance[j])+"\n"+str(arraySeqs[i])+"\n"
    return string.rstrip(string[-1])


def save_fasta_format (file_name,data) :
    output_file = open(file_name+".fa", "w")
    output_file.write(data)
    output_file.close()

def save_nk_format (file_name,data) :
    output_file = open(file_name+".nk", "w")
    output_file.write(data)
    output_file.close()

#Deuxième simplification

def simplification_v2_arret_nombre_de_noeuds(NkFile,FastaFile,Seuil_Abondance,nb_noeuds) : #simplification forçant la réduction de l'arbre à nb_noeuds clonotypes, attention perte probable d'information, retourne directement une chaine de caractère au format nk
    Abondance=BasicSeq.readFastaAbundance(FastaFile)[3]
    arraySeqs=BasicSeq.readFastaAbundance(FastaFile)[2]
    labels=BasicSeq.readFastaAbundance(FastaFile)[0]
    dix_abondants=[]
    for k in dix_plus_abondants(Abondance):
        dix_abondants.append(k[0])
    tot=0 ## On calcule le nombre total de séquences observées pour calculer le pourcentage de représentation de chaque séquence unique
    for i in Abondance :
        tot+=Abondance[i]
    M_distance = BasicSeq.createAdjMatrix(arraySeqs)    
    vmax=np.max(M_distance)
    t=Tree(NkFile,format=1)
    kept_sequences=[]
    for node in t.traverse() :
        kept_sequences.append(node.name)
    while len(kept_sequences)>nb_noeuds :
        imin=0
        jmin=0
        vmin = vmax+1
        L=t.get_leaves()
        deleted = False #variable permettant de savoir si on a déjà enlevé une feuille de l'arbre
        for i in L :
            if i.name not in labels : #on aura ici un noeud silencieux qui se retrouverait en position de feuille : pas pertinent donc on la taille
                del kept_sequences[kept_sequences.index(i.name)]
                deleted = True
                
        if deleted == False :
            que_abondant=True
            que_abondant_max=True
            for i in L :
                if i.name not in dix_abondants : que_abondant_max=False
                if Abondance[i.name]/tot<Seuil_Abondance : que_abondant=False
            if que_abondant :
                break
            if que_abondant_max :
                break
            for k in range(len(L)-1):
                for l in range (k+1,len(L)):
                    leavei=L[k]
                    leavej=L[l]
                    i=labels.index(leavei.name)
                    j=labels.index(leavej.name)
                    if M_distance[i][j]==vmin :
                        if Abondance[leavei.name]+Abondance[leavej.name]<Abondance[labels[imin]]+Abondance[labels[jmin]] :
                            vmin = M_distance[i][j]
                            imin=i
                            jmin=j
                    elif M_distance[i][j]<vmin :
                        vmin = M_distance[i][j]
                        imin=i
                        jmin=j

            M_distance[imin][jmin]=vmax+10
        
            if Abondance[labels[imin]]>Abondance[labels[jmin]] :
                if Abondance[labels[jmin]]/tot<Seuil_Abondance and labels[jmin] not in dix_abondants :
                    del kept_sequences[kept_sequences.index(labels[jmin])]
                    
                
                    for i in range(len(M_distance)):
                        M_distance[i][jmin]=vmax+10
                        M_distance[jmin][i]=vmax+10
                
            else :
                if Abondance [labels[imin]]/tot<Seuil_Abondance and labels[imin] not in dix_abondants:
                    del kept_sequences[kept_sequences.index(labels[imin])]
                    
                    for j in range(len(M_distance)):
                        M_distance[imin][j]=vmax+10
                        M_distance[j][imin]=vmax+10
        
        
        t.prune(kept_sequences)
        

    return t.write(format=1)

def simplification_v2_arret_nombre_de_noeuds_retour_fasta(NkFile,FastaFile,Seuil_Abondance,nb_noeuds) : #même algorithme mais renvoyant cette fois juste la liste des séquences conservées, dans l'objectif de réécrire un fichier fasta
    Abondance=BasicSeq.readFastaAbundance(FastaFile)[3]
    arraySeqs=BasicSeq.readFastaAbundance(FastaFile)[2]
    labels=BasicSeq.readFastaAbundance(FastaFile)[0]
    dix_abondants=[]
    for k in dix_plus_abondants(Abondance):
        dix_abondants.append(k[0])
    tot=0 ## On calcule le nombre total de séquences observées pour calculer le pourcentage de représentation de chaque séquence unique
    for i in Abondance :
        tot+=Abondance[i]
    M_distance = BasicSeq.createAdjMatrix(arraySeqs)    
    vmax=np.max(M_distance)
    t=Tree(NkFile,format=1)
    kept_sequences=[]
    for node in t.traverse() :
        kept_sequences.append(node.name)
    while len(kept_sequences)>nb_noeuds :
        imin=0
        jmin=0
        vmin = vmax+1
        L=t.get_leaves()
        deleted = False #variable permettant de savoir si on a déjà enlevé une feuille de l'arbre
        for i in L :
            if i.name not in labels : #on aura ici un noeud silencieux qui se retrouverait en position de feuille : pas pertinent donc on la taille
                del kept_sequences[kept_sequences.index(i.name)]
                deleted = True
        if deleted == False :
            que_abondant=True
            que_abondant_max=True
            for i in L :
                if i.name not in dix_abondants : que_abondant_max=False
                if Abondance[i.name]/tot<Seuil_Abondance : que_abondant=False
            if que_abondant :
                break
            if que_abondant_max :
                break
            for k in range(len(L)-1):
                for l in range (k+1,len(L)):
                    leavei=L[k]
                    leavej=L[l]
                    i=labels.index(leavei.name)
                    j=labels.index(leavej.name)
                    if M_distance[i][j]==vmin :
                        if Abondance[leavei.name]+Abondance[leavej.name]<Abondance[labels[imin]]+Abondance[labels[jmin]] :
                            vmin = M_distance[i][j]
                            imin=i
                            jmin=j
                    elif M_distance[i][j]<vmin :
                        vmin = M_distance[i][j]
                        imin=i
                        jmin=j

            M_distance[imin][jmin]=vmax+10
        
            if Abondance[labels[imin]]>Abondance[labels[jmin]] :
                if Abondance[labels[jmin]]/tot<Seuil_Abondance and labels[jmin] not in dix_abondants :
                    del kept_sequences[kept_sequences.index(labels[jmin])]
                
                    for i in range(len(M_distance)):
                        M_distance[i][jmin]=vmax+10
                        M_distance[jmin][i]=vmax+10
                
            else :
                if Abondance [labels[imin]]/tot<Seuil_Abondance and labels[imin] not in dix_abondants:
                    del kept_sequences[kept_sequences.index(labels[imin])]
                
                    for j in range(len(M_distance)):
                        M_distance[imin][j]=vmax+10
                        M_distance[j][imin]=vmax+10
        
        
        t.prune(kept_sequences)
    return kept_sequences

def generation_fasta_simplifie_advanced(NKfile,FastaFile,Seuil_Abondance,nombre_min_noeuds ) : #génère le fichier fasta correspondant à la simplification avancée
    keptSequencesNames=simplification_v2_arret_nombre_de_noeuds_retour_fasta(NKfile,FastaFile,Seuil_Abondance,nombre_min_noeuds)
    labels = BasicSeq.readFastaAbundance(FastaFile)[0]
    arraySeqs = BasicSeq.readFastaAbundance(FastaFile)[2]
    Abondance = BasicSeq.readFastaAbundance(FastaFile)[3]
    string =""
    for i,j in enumerate(labels) :
        if j in keptSequencesNames :
            string = string + ">"+str(j)+"@"+str(Abondance[j])+"\n"+str(arraySeqs[i])+"\n"
    return string.rstrip(string[-1])


def enchainement_simplification(NkFile,FastaFile,file_name,Seuil_Abondance=0.01,nombre_min_noeuds=30,nombre_noeuds=30): #enchainement des deux simplifications, mettre en entrée l'arbre initial à simplifier
    level1_data_nk = simplification_abondance(NkFile,FastaFile,Seuil_Abondance,nombre_min_noeuds)
    level1_data_fa = generation_fasta_simplifie(NkFile,FastaFile,Seuil_Abondance,nombre_min_noeuds)
    save_fasta_format(file_name+"_simplification1",level1_data_fa)
    level2_data_nk = simplification_v2_arret_nombre_de_noeuds(level1_data_nk,file_name+"_simplification1.fa",Seuil_Abondance,nombre_noeuds)
    level2_data_fa = generation_fasta_simplifie_advanced(level1_data_nk,file_name+"_simplification1.fa",Seuil_Abondance,nombre_min_noeuds)
    save_fasta_format(file_name+"_simplification2",level2_data_fa)
    return level1_data_nk, level2_data_nk

def nettoyage_arbre_simplifie(tree, file_name) :
	t=Tree(tree,format=1)
	for node in t.traverse():
		if node.is_root()==False and node.name !="germline": #on parcourt tous les noeuds qui ne sont pas la racine
			parent=node.up
			if str(parent.name)[0]=="n": #le parent de ce noeud est un noeud silencieux
				children = parent.get_children()
				prev_node = parent.up #noeud parent du noeud silencieux
				if len(children)==1 : #ce noeud silencieux n'a qu'un enfant dans la version simplifiée
					parent.detach()
					for child in children:
						prev_node.add_child(child)
			if(str(node.name)[0]=="n" and node.is_leaf()) : #la feuille de l'arbre correspond à un noeud silencieux
				node.delete()								
	save_nk_format(file_name,t.write(format=1))

def main():
	usage = usage = "python Simplification.py -n <newick_file> -a <alignment_file> \n"
	parser = OptionParser(usage)
	parser.add_option("-n", "--newick_file", dest="newick_file",  help="file containing the tree in newick format")
	parser.add_option("-a", "--alignment_file", dest="alignment_file",  help="file containing the alignment of the sequences")
			
	(options, args) = parser.parse_args()
	
	if len(sys.argv) != 5:
		parser.error("Incorrect number of arguments")
	
	newick_file = options.newick_file
	alignment_file = options.alignment_file

	file_name = os.path.splitext(newick_file)[0]
	level1_simplification, level2_simplification = enchainement_simplification(newick_file,alignment_file,file_name)
	nettoyage_arbre_simplifie(level1_simplification, file_name+"_simplification1")
	nettoyage_arbre_simplifie(level2_simplification, file_name+"_simplification2")
	print("    # Simplification : done")

if __name__ == "__main__":
	main()
