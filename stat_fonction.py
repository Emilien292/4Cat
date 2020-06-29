"""
4Cat
this software analyse 4C data from Bam 
Copyright (C) 2020 N'guyen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
"""
# -*- coding: utf-8 -*-
import numpy as np
from sklearn.isotonic import IsotonicRegression
from statistics import *
import scipy.stats as scipy
from math import *
#from moduleLectureFiltreBam import *

#Adaptation des donnees pour le package scikit learn
def formatage_donnee(data):
	#Creation des labels une ligne par replicat et une colonne par position
	data_comptage_y = [fenetre[2:] for fenetre in data]
	indexFenetre = [fenetre[:2] for fenetre in data]
	#Transformer les données de comptage en tableau numpy pour utiliser scikit learn
	data_comptage_y = np.asarray(data_comptage_y)
	data_comptage_y = np.transpose(data_comptage_y)
	#print(data_comptage_y)
	return indexFenetre,data_comptage_y


#Adaptation des donnees pour le package scikit learn
def formatage_donnee_init(data):
	#Creation des labels une ligne par replicat et une colonne par position
	data_comptage_y = []
	for i in range(len(data[1])):
		replicat = []
		for j in range(len(data)):
			replicat.append(data[j,i])
		data_comptage_y.append(replicat)
	#Création de la table dentraînement -> valeur moyenne dans tout les réplicats à la position considerer
	data_comptage_x = []
	for i in range(len(data)):
		moyenne = mean(data[i,1:len(data[i])])
		data_comptage_x.append((data[i,0],moyenne))
	#Transformer les données de comptage en tableau numpy pour utiliser scikit learn
	data_comptage_y = np.asarray(data_comptage_y)
	data_comptage_x = np.asarray(data_comptage_x)
	return data_comptage_x,data_comptage_y

def regression_monotone(data_target,data):
	regression_data = []
	
	#Realisation de la regression monotone qui nous donnes un bruit de fond théorique
	#On fit et transform les données pour chacun des replicats indépendamment
	ir = IsotonicRegression()
	#data_target_transpose = np.transpose(data_target)
	for replicat in data_target:
		regression = ir.fit_transform(np.arange(0,len(replicat),1),replicat)
		regression_data.append(regression)
	regression_data = np.asarray(regression_data)
	#print(len(regression_data))
	return regression_data
	
def regression_monotone_initial(data_training,data_target,data):
	regression_data = []
	#Realisation de la regression monotone qui nous donnes un bruit de fond théorique
	#On fit et transform les données pour chacun des replicats indépendamment
	ir = IsotonicRegression()
	for i in range(len(data_target)-1):
		regression = ir.fit_transform(data_training[0:len(data_training),1],data_target[i+1,0:len(data)])
		regression_data.append(regression)
	regression_data = np.asarray(regression_data)
	return regression_data

#Calcul du ratio et de la différence entre données observé et bruit de fond théorique
def ratio_delta(regression_data,data_comptage_y,data):
	ratio = []
	delta = []
	for j in range(len(regression_data)):
		replicatRatio = []
		replicatDelta = []
		for i in range(len(regression_data[0])):
			replicatRatio.append((data_comptage_y[j,i]+1)/(regression_data[j,i]+1))
			replicatDelta.append(data_comptage_y[j,i]-regression_data[j,i])
		delta.append(replicatDelta)
		ratio.append(replicatRatio)
	return ratio, delta



#Calcul le produit de rang pour chaque positions a partir des differences calculee dans la fonction précédente
def rang(delta):
	list_sort = []
	list_score_rang = []
	list_score_rang.append(delta[0])
	# On parcourt le tableau des difference et on donne dans chaque replicat un rang pour chaque position le rang 1 etant la valeur la plus forte
	for i in range(len(delta)):
		rank=(len(delta[i])-scipy.rankdata(delta[i]))+1 # on soustrait a la taille de la liste pour que la plus forte difference soit au rang 1 la 2eme au rang etc
		list_sort.append(rank)
	list_sort = np.asarray(list_sort)
	# On cacul le produit de rang pour chaque position, axis=0 pour calculer par colonne
	list_score_rang.append(np.prod(list_sort,axis=0))
	list_score_rang = np.asarray(list_score_rang)
	return(list_score_rang)

#Calcul de la p-value a partir de l'approximation de la loi de permutation.
def p_value(rank,K):
	pvalue = []
	#print(rank[0])
	#Calcul a partir des produits de rangs du nombre de position et du nombre de replicat des quantiles de la densité cumule de la loi gamma pour calculer les p-values
	p = (rank[1,:]/(len(rank[1])+1)**K)
	for i in p:
		ln = -log(i)
		#Calcul de la p-value
		pvalue.append(1 - scipy.gamma.cdf(ln,K))
	return(pvalue)

#calcul des pvalue ajusté (FDR)
def correct_pvalues_for_multiple_testing(pvalues):
	from numpy import array, empty		
	pvalues = array(pvalues) 
	n = pvalues.shape[0]
	new_pvalues = empty(n)	 
	values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]	
	#met les p-value dans l'odre de leur rang		  
	values.sort()
	values.reverse()		  
	new_values = []
	#procedure de Benjamini_Hochberg
	for i, vals in enumerate(values):		 
		rank = n - i
		pvalue, index = vals		  
		new_values.append((n/rank) * pvalue)	
	#enlever les valeurs superieurs a 1	  
	for i in range(0, int(n)-1):  
		if new_values[i] < new_values[i+1]:
			new_values[i+1] = new_values[i]
	#remet les p_values dans l'ordre de leur positions
	for i, vals in enumerate(values):
		pvalue, index = vals
		new_pvalues[index] = new_values[i]		  
	return new_pvalues

#Recherche des positions significatives
def fragment_significatif(fdr_list,pos,FDR):
	#Recherche les index des positions significatif
	index_fdr_significatif = np.where(fdr_list < FDR)
	pos = np.asarray(pos)
	#On garde uniquement les positions correspondante au index
	pos_significative = pos[index_fdr_significatif]
	return(pos_significative,index_fdr_significatif[0])

#Filtre selon l'effet taille a partir de la moyennnes des ratios pour chaque position
def filtre_effect_size(ratio,QWR):
	ratio = np.asarray(ratio)
	mean = np.sum(ratio[1:,],axis=0)/len(ratio[1:,])
	Q3 = np.quantile(mean, 0.75)
	IQR = np.quantile(mean, 0.75)-np.quantile(ratio, 0.25)
	#Calul du seuil	
	threshold = Q3 + QWR*IQR
	#Recherche les index des positions superieur au seuil
	index_mean = np.where(mean > threshold)
	#On garde uniquement les positions correspondante au index
	#pos_significative = position[index_mean]
	return(index_mean[0])



#On garde les n meilleurs p-value si on ne trouve aucunes positions significative
def best_pvalue(pvalue,nb_meilleur_candidat,pos):
	pval = np.asarray(pvalue)
	meilleur_candidat = pval[0:nb_meilleur_candidat]
	meilleur_candidat.sort()
	meilleur_candidat = np.asarray(meilleur_candidat)
	pos_candidat = []
	#On cherche les n meilleurs p-value
	for i in pval[nb_meilleur_candidat:,]:
		if i < max(meilleur_candidat) and i not in meilleur_candidat:
			index = np.where(i < meilleur_candidat)
			meilleur_candidat = np.insert(meilleur_candidat,index[0][0],i)
			meilleur_candidat = np.delete(meilleur_candidat,nb_meilleur_candidat)
	#On cherche les postions correspondantes au n meilleur p-value	
	for i in meilleur_candidat:
		index = np.where(i == pvalue)
		#print(index[0])
		pos_candidat.append(pos[index[0][0]])
	return(pos_candidat,meilleur_candidat)


#Fonction de lancement de l'analyse
def main_stat(liste_pose_Candidat,chromosome,alphaFDR=0.01,QWR = 1,nb_meilleur_candidat = 5,basename = "default_peak", repertoire = "resultats", code = 1000):
	data = liste_pose_Candidat
	#data = np.asarray(data)
	indexFenetre,Y = formatage_donnee(data)
	regression_data = regression_monotone(Y,data)
	#print(len(indexFenetre))
	print("calcul du bruit de fond terminé")

	ratio, delta = ratio_delta(regression_data,Y,data)
	rank = rang(delta)
	

	pvalue = p_value(rank,len(data[1])-2)
	print("calcul des P-values terminé")

	fdr = correct_pvalues_for_multiple_testing(pvalue)
	pos_significative,fdr_index_sinificatif = fragment_significatif(fdr,indexFenetre,alphaFDR)
	resultats = []
	#Si aucunes positions significatives on ecrit dans le fichier de sortie les n meilleurs p-value 
	if len(pos_significative) == 0:
		pos_candidat, meilleur_pvalue = best_pvalue(pvalue,nb_meilleur_candidat,indexFenetre)
	
		pos = np.asarray(indexFenetre)
		index_significatif = np.where(pos_candidat == pos)
		produit_rang = []

		for i in pos_candidat:
			index_significatif = np.where(i == pos)
			produit_rang.append(rank[1,index_significatif[0]])
	
		with open(repertoire + "/" + basename + "_peak_non_significatif.txt", "w") as fichier:
			fichier.write("chromosome"+ "\t"+ "start"+ "\t" + "end" + "\t" + "pvalue" + "\t" +  "\t" +  "produit de rang" + "\n")
			for i, pos  in enumerate(pos_candidat):
				fichier.write(str(pos[0]) + "\t" +str(pos[1]) + "\t" + str(meilleur_pvalue[i]) + "\t" + str(produit_rang[i]) + "\n")
			print("Fin de la recherche des pics" + "\n" +"Aucun pic significatif trouvé les "+str(nb_meilleur_candidat)+" meilleurs candidats non significatif seront donnés dans le fichier de sortie "+repertoire + "/" + basename + "_peak_non_significatif.txt dans le dossier resultats")
	else:
		taille_effet = filtre_effect_size(ratio,QWR)
		significatif = list(set(taille_effet) & set(fdr_index_sinificatif))
		if(len(significatif)==0):
			significatif = list(set(fdr_index_sinificatif))
		#print(fdr_index_sinificatif)
		#significatif = fdr_index_sinificatif
		with open (repertoire + "/" + basename + "_peak_significatif.txt","w") as fichier:
			fichier.write("chromosome"+ "\t"+ "start"+ "\t" + "end" + "\t" + "pvalue" + "\t" +  "\t" +  "produit de rang" + "\n")
			for i, pos  in enumerate(significatif):
				resultats.append([chromosome,indexFenetre[pos][0],indexFenetre[pos][1],fdr[pos],code])
				fichier.write(chromosome+'\t'+str(indexFenetre[pos][0])+'\t'+str(indexFenetre[pos][1]) + "\t" + str(fdr[pos]) + "\t" + str(rank[1,pos]) + "\n")
		print("Fin de la recherche des pics." +  "\n" +"les pics significatifs se trouvent dans le fichier " + repertoire + "/" + basename + "_peak_significatif.txt dans le dossier resultats")
	return resultats

#main_stat()
