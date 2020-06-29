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
import re
import pysam 
import matplotlib.pylab as plt
import ruptures as rpt
import numpy as np
from copy import copy, deepcopy
import math
import functools
"""@package docstring
Documentation for this module.
 
More details.
"""

def extractSiteToGenome(genome,enzyme1,enzyme2):
    """! Documentation pour la fonction extractSiteToGenome
        cette fonction prends en parametre la sequence d'un chromosome (un string),les 2 enzyme E1 et E2 (string)
        Elle renvoie une liste de 2 liste, la premiere renvoie le debut de chaque occurence de la premiere enzyme
        la deuxieme renvoie l'ocurence de la deuxieme Enzyme

        @param genome String
        @param enzyme1 String
        @param enzyme2 String

        @return List of 2 list 
    """
    l2 = [m.start(0) for m in re.finditer(enzyme1,genome)]
    l1 = [m.end(0) for m in re.finditer(enzyme2,genome)]
    return [l1,l2]

def create2Fenetre(l1,l2,distFrag = 20):
    """! Documentation pour la fonction create2fenetre
        cette fonction prends en parametre 2 liste composer de postion des sites de restriction trié respectivement E1 et E2 par ordre croissant
        et un int qui corespond a la distance minamle requise pour crée un "fragment" valide.
        La fonction renvoie une liste qui corespond au fragment au fragment valide E1 à E2 ou E2 à E1 
        et une seconde liste possedant les fragments valide E1 à E1.

        @param l1 liste de int correspondant au position des sites d'enzymes de la premiere Enzyme E1
        @param l2 liste de int correspondant au position des sites d'enzymes de la deuxieme Enzyme E2
        @param distFrag int corespondant a la distance minimale requise entre 2 position pour crée un fragment valide

        @return Une liste des 2 liste des 2 type de fragment valide on parlera de fenetre: [positionStart,positionEnd]
    """
    pl1 = 0 #pointeur de la fenetre l1
    pl2 = 0 #pointeur de la fenetre l2
    maxL1 = len(l1)
    maxL2 = len(l2)
    lF1 = []
    lF2 = []
    while(l1[pl1] > l2[pl2]):
        #on initialise les compteurs
        pl2+=1
    while(pl1 < maxL1):
        #print(l1[pl1])
        #print(l2[pl2])
        #print(l1[pl1])
        if (pl2 >= maxL2):
            return [lF1,lF2]
        if (l1[pl1] < l2[pl2]-distFrag):
            # on a un fragment valide possible  E1 à E2
            if(pl1 + 1 >= maxL1):#extremiter du genome
                return [lF1,lF2]
            elif(l1[pl1+1] > l2[pl2]):
                lF1.append([l1[pl1],l2[pl2]]) # fenetre E1 à E2
                #on va crée la fenetre E1 à E1 avec un fragment valide
                lF2.append([l1[pl1],l1[pl1+1]])# ici on ne s'interresse pas a savoir si ce fragment
                #est composer de 1 ou 2 fragment valide.
            pl1 += 1 
        elif(l1[pl1] > l2[pl2]):
            # on cherche a savoir maintenant si on a un fragment E2 à E1
            while (l2[pl2] < l1[pl1]):
                if (l2[pl2]+distFrag <= l1[pl1]):
                    if (pl2+1 >= maxL2):
                        lF1.append([l2[pl2],l1[pl1]]) # fenetre E2 à E1
                        #alors c'est un fragment E1 à E1 avec un fragment valide
                        if(l1[pl1-1]!=lF2[-1][0]):
                            lF2.append([l1[pl1-1],l1[pl1]])
                        return [lF1,lF2]
                    if (l2[pl2+1] > l1[pl1]):
                        lF1.append([l2[pl2],l1[pl1]]) # fenetre E2 à E1
                        #alors c'est un fragment E1 à E1 avec un fragment valide
                        if(l1[pl1-1]!=lF2[-1][0]):
                            lF2.append([l1[pl1-1],l1[pl1]])
                        pl2 += 1
                    elif (pl2+1 >= maxL2):
                        return [lF1,lF2]
                    else:
                        pl2 += 1
                elif (pl2+1 >= maxL2):
                    return [lF1,lF2]
                else:
                    pl2 += 1
        else:
            pl1 += 1
    return [LF1,LF2]

def calculFacteurNormalisation(listeBam):
    """! Documentation pour la fonction calculFacteurNormalisation
        cette fonction prends en parametre la liste des bam et calcule le rapport 
        entre le nombre de read du fichier et le nombre de read augmenter 

        @param listeBam list de string correspondant au chemin des fichier bam

        @return une liste avec pour chaque fichier le facteur de normalisation nécéssaire  
    """
    
    listCorrecteur = []
    for filename in listeBam:
        nbrRead = 0
        tmp = pysam.idxstats(filename)
        for  l in tmp.split(('\n'))[:-1]:
            #obtient le nombre total de read mapper
            nbrRead += int(l.split('\t')[2])
        listCorrecteur.append(nbrRead)
    maxNumberRead = max(listCorrecteur)
    #print(listCorrecteur)
    for i in range(len(listCorrecteur)):
        listCorrecteur[i] = maxNumberRead/listCorrecteur[i]       
    return(listCorrecteur)

def filterRead(fileBam,listeFenetreE1aE2,listeFenetreE1aE1,chrFilter,mapQ = 20,dist = 10):
    """! Documentation pour la fonction filterRead
        cette fonction prends en parametre un fichier bam trier et indexer les listes des 2 types de fragments valide
        le chromosome sur lequel on filtre, mapQ est la qualiter limite minamale requis pour accepter un read,
        dist est la distance maximale entre les extremiter des reads et des fragments de type E1 a E2 ou E2 à E1

        @param fileBam le fichier bam (son chemin)
        @param listeFenetreE1aE2 une liste de meme type de fenetre [start,end,rep1,rep2,...]
        @param listeFenetreE1aE1 une liste de meme type de fenetre [start,end,rep1,rep2,...]
        @param chrFilter le nom du chromosome sur lequel on quantifie le nombre de read pour un replicat sur chque fenetre
        @param mapQ int seuil poru lequel un read est accepter dans le comptage
        @param dist int qui donne le maximun de distance entre l'une des extremiter du read et de la fenetre

        @return Une liste des deux liste mis en parametre en rajoutant le comptage du bam dans chaque fenetre a la derniere position ainsi que le nom du chromosome du bam
    """
    bamFile = pysam.AlignmentFile(fileBam,"rb")
    plF1 = 0
    plF2 = 0
    count = 0
    #on initialise le comptage
    for i in listeFenetreE1aE1:
        i.append(0)
    for i in listeFenetreE1aE2:
        i.append(0)
    headerStr = str(bamFile.header)
    #print(headerStr)
    present = False
    for i in headerStr.split("\n"):
        j = i.split("\t")
        if (j[0] == "@SQ"):
            if (j[1].split(":")[1].strip().replace("chr","") == chrFilter.replace("chr","")):
                chrFilter = j[1].split(":")[1].strip()
                tailleChr = j[2].split(":")[1].strip()
                present = True
                break
    if(present != True):
        return [listeFenetreE1aE2,listeFenetreE1aE1,chrFilter]
    for read in bamFile.fetch(chrFilter):
        if (read.mapping_quality >= mapQ):
            #count +=1
            start = read.reference_start
            end = read.reference_end
            while(start > listeFenetreE1aE2[plF1][1]):
                #passe a la fenetre suivante
                plF1 += 1
            #print([start,end])
            #print(listeFenetreE1aE2[plF1])
            if (end > listeFenetreE1aE2[plF1][0]):
                #possible chevauchement
                if ((start < listeFenetreE1aE2[plF1][0] + dist)and(start >listeFenetreE1aE2[plF1][0] - dist)):
                    #on chevauche bien au niveau des sites d'enzymes +- une distance 
                    listeFenetreE1aE2[plF1][-1] += 1
                    count += 1
                    while(listeFenetreE1aE1[plF2][1] <= listeFenetreE1aE2[plF1][0]):
                        plF2 += 1
                    listeFenetreE1aE1[plF2][-1] += 1
                elif ((end < listeFenetreE1aE2[plF1][1] + dist)and(end >listeFenetreE1aE2[plF1][1] - dist)):
                    listeFenetreE1aE2[plF1][-1] += 1
                    count += 1
                    while(listeFenetreE1aE1[plF2][1] <= listeFenetreE1aE2[plF1][0]):
                        plF2 += 1
                    listeFenetreE1aE1[plF2][-1] += 1
                #print(listeFenetreE1aE2[plF1])   
    return [listeFenetreE1aE2,listeFenetreE1aE1,chrFilter]

def separeList(listeFenetre,cordSeparate):
    """! Documentation pour la fonction separeList
        cette fonction prends en parametre une liste corespondant a un type de fragment et une liste,tuple de deux int qui corespond au viewpoint.
        Elle permettera de separer la liste en deux selon les cordonnées du viewpoint.
        @param cordSeparate une liste de 2 int debut et fin (ordre croissant)

        @return list de deux sousliste provenant de la liste passez en parametre()
    """
    # cette fonction permet de séparer une fenetre en 2 selon les cordonnées du viewpoints
    for vp,fenetre in enumerate(listeFenetre):
        if (fenetre[1] >= cordSeparate[0]):
            i = vp
            if(fenetre[1] >= cordSeparate[1]):
                j = vp+1
                break
            else:
                j = vp+2
                break
    return([listeFenetre[:i],listeFenetre[j:]])




def SepareInfluenceVp(listeFenetre,methods = "ruptures",n_bkps = 2,nbrJump = 100):
    """! Documentation pour la fonction SepareInfluenceVp
        cette fonction prends en parametre une liste de fenetre et cherche a calculer des points de cassure par une methode. 
        La methode utiliser pour le moment est uniquements ruptures. 

        @param listeFenetre une liste de fenetre
        @param methods string contenant le nom de la methode. Pour le moment seul ruptures est disponible
        @param n_bkps int corespondant au nombre de points de cassure chercher par default 2 on cherche a separer 3 zone en fonction de l'influence du viewpoint
        @param nbrjump un int pour des raison de cout en memoire et en temps rupture ne calcule pas un coeficient entre chaque element de la liste, il saute de position en position par defaut ça valeur est 1000

        @return list de deux 3 sous liste de fenetre correspont au 3 zone peu influencer, influencer fortement influencer.
    """
    # creation of data
    signal = []
    for i in listeFenetre:
        tmp = []
        for j in i[2:]:
            tmp.append(math.log(float(j)+1.0))
        signal.append(tmp)
    signal = np.asarray(signal)
    # change point detection
    model = "l1"  #"l1"  # "l2", "rbf"
    jump = int(len(listeFenetre)/nbrJump)
    algo = rpt.Dynp(model=model, min_size=3, jump=jump).fit(signal)
    my_bkps = algo.predict(n_bkps=n_bkps)
    # show results
    rpt.show.display(signal, my_bkps[:-1], my_bkps[:-1], figsize=(10, 6))
    #plt.show()
    j = 0
    resultatsSeparer = []
    for i in my_bkps:
        resultatsSeparer.append(listeFenetre[j:i-1])
        j = i

    return resultatsSeparer

def createFenetreConsecutif(listeFenetreE1aE1,nbConsecutif = 10,distmax = 1):
    """! Documentation pour la fonction createFenetreConsecutif
        cette fonction cherche a regrouper un nombre de fenetre en une seule avec pour critere le nombre de fenetre souhaiter et la distance maximale entre 2 fragment pour etre regrouper. les fenetre seront chevauchante


        @param listeFenetre une liste de fenetre que l'on cherche a regrouper
        @param nbConsecutif int corespondant au nombre de fenetre que l'on cherche à regrouper
        @param distmax float corespendant a la distance maximal en kb entre deux fenetre regrouper

        @return list de fenetre chevauchante prevenant d'un regroupement entre fenetre de la liste en parametre.
    """
    distmax = distmax * pow(10,3)
    #print(distmax)
    listeFenetreConsecutif = []
    pLE1  = 0
    fin = len(listeFenetreE1aE1)-nbConsecutif
    start = listeFenetreE1aE1[0][0]
    end = listeFenetreE1aE1[0][1]
    pnC = 0
    #print(len(listeFenetreE1aE1))

    while (pLE1 < fin) :
        pnC = 1
        start = listeFenetreE1aE1[pLE1][0]
        end = listeFenetreE1aE1[pLE1][1]
        tmpEnd = end
        comptage = copy(listeFenetreE1aE1[pLE1][2:])
        tmpComptage = copy(listeFenetreE1aE1[pLE1][2:])
        while ((end - start < distmax) and (pnC < nbConsecutif)):
            tmpEnd = end
            tmpComptage = copy(comptage)
            for j,nbrRead in enumerate(comptage):
                comptage[j] += listeFenetreE1aE1[pLE1+pnC][2+j]
            pnC += 1
            if(pnC < nbConsecutif):
                end = listeFenetreE1aE1[pLE1+pnC][1]
        if(end - start < distmax):
            fenetre = [start,end]
            for i in comptage:
                fenetre.append(i)
            if (len(listeFenetreConsecutif) != 0):
                if(listeFenetreConsecutif[-1][1]!=end):
                    listeFenetreConsecutif.append(fenetre)   
            else:
                listeFenetreConsecutif.append(fenetre)
        else:
            fenetre = [start,tmpEnd]
            for i in tmpComptage:
                fenetre.append(i)
            if (len(listeFenetreConsecutif) != 0):
                if(listeFenetreConsecutif[-1][1]!=tmpEnd):
                    listeFenetreConsecutif.append(fenetre)
            else:
                listeFenetreConsecutif.append(fenetre)
        pLE1 += 1
    return listeFenetreConsecutif
