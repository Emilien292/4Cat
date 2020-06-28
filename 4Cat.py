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
from moduleBam import *
from Zscore import *
from utils import *
import argparse
from stat_fonction import main_stat

parser = argparse.ArgumentParser(description='Analyse de détection de pics a partir de données 4cSeq')
optional = parser._action_groups.pop() # Edited this line

# remove this line: optional = parser...
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional) # added this line
required.add_argument("-e1",'--enzyme1', type=str,required=True, 
                    help="Sequence de l'enzyme de restriction 1")
required.add_argument("-e2",'--enzyme2', type=str,required=True, 
                    help="Sequence de l'enzyme de restriction 2")
optional.add_argument('--alphaFdr',"-FDR", type=float, default=0.05,
                    help="Seuil de fdr en dessous duquel les pics seront considérés comme significatif (default 0.05) ")
optional.add_argument('--nb_meilleur_candidat', type=float, default=5,
                    help="Si aucun pic n'est trouvé avec un fdr inferieur au seuil choisis lors d'une analyses statistique, le programme donnera les n pics ayant les plus basses P-values (default 5)")
optional.add_argument('--QWR', type=float, default=0.1,
                    help="Permet de filtrer les résultats en fonctions de la taille de l’effet, plus il est grand plus on est strict (default 1)")
optional.add_argument('--mappingQuality',"-mapQ", type=float, default=20,
                    help="définit la qualiter minimale requise pour un alignements, tous les read ayant une valeur inférieur seront retiré")
optional.add_argument('--nbrConsecutif',"-nbrC", type=int, default=10,
                    help="définit le nombre de fragment de type E1 a E1 regrouper, ce parametre a une influence sur les regions rechercher")
optional.add_argument('--maxDist',"-mD", type=int, default=5,
                    help="définit la distance maximun de regroupement entre fragment de type E1 a E1, cela evite le regroupement de ffragment trop éloigner")                    
optional.add_argument('--createBedgraph',"-cBg", type=bool, default=False,
                    help="Crée pour chaque bam un bedgraph avec uniquement les read filtré")                    
optional.add_argument('--repertoireBedgraph',"-rBg", type=str, default="resultats",
                    help="donne le repertoire de la création des bedgraph")  
optional.add_argument('--cordinateViewpoint',"-cordVp", type=str,
                    help="indique les cordonnées présuposer du viewpoint au format start-end ou uniquement une position")                 
optional.add_argument('--allChrom',"-allChrom", type=bool, default = False,
                    help="indique si on souhaite analyser tous les chromosomes ou uniquement le chromosome du viewpoint")                 
required.add_argument('--listeBam', nargs="+",required=True,
                    help="Chemin(s) vers le(s) fichier(s) BAMs utilisé(s)")
required.add_argument("-g",'--genome', type=str,required=True,
                    help="Chemin vers le génome au format fasta utilisé pour les alignements")
required.add_argument("-chr",'--chrViewPoint', type=str,required=True,
                    help="Nom du chromosome sur lequel se trouve le viewPoint")
args = parser.parse_args()

#print(args.accumulate(args.enzyme))
listeBam = args.listeBam
print(listeBam)
genome = args.genome
enzyme1 = args.enzyme1
enzyme2 = args.enzyme2
chrViewPoint = args.chrViewPoint
mapQ = args.mappingQuality
nbrConsecutif = args.nbrConsecutif
maxDist = args.maxDist
bigWig = args.createBedgraph
repBigwig = args.repertoireBedgraph
allChrom = args.allChrom
nbpks = 2
cordVp = args.cordinateViewpoint
cordVp = "23044000"
#######################################
###
alphaFDR = args.alphaFdr
nb_meilleur_candidat = args.nb_meilleur_candidat
QWR = args.QWR
fileGenome = open(genome,"r")
nomChroLocal = fileGenome.readline().strip()
sequenceLocal = ""
correcteur = [1 for i in listeBam]
correcteur = calculFacteurNormalisation(listeBam)
#Analyse 
for line in fileGenome:
    if(">" in line):
        nomChroLocal = nomChroLocal.split(" ")[0].replace("chr","").replace("Chr","").replace(">","")
        if(nomChroLocal == chrViewPoint.replace("chr","").replace("Chr","")):
            print("analyse du viewpoint")
            listeSite1,listeSite2 = extractSiteToGenome(sequenceLocal.upper(),enzyme1,enzyme2)
            ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1 = create2Fenetre(listeSite1,listeSite2,25)
            
            ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1,chrName = filterRead(listeBam[0],ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1,nomChroLocal,mapQ,correcteur[0])
            ###Calcul du vp si l'utilisateur n'a pas mis ses cordonnées:
            if("-" in cordVp):
                tmp = cordVp.split("-")
                cordVp = [int(tmp[0]),int(tmp[1])]
            elif(cordVp != ""):
                tmp = int(cordVp)
                cordVp = [tmp,tmp]
            else:
                l = [i[2] for i in ListeFenetreE1aE1] 
                j = l.index(max(l))
                cordVp = ListeFenetreE1aE1[j][0:2]
            
            for i,fileBam in enumerate(listeBam[1:]):
                ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1,chrName = filterRead(fileBam,ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1,nomChroLocal,mapQ,correcteur[i+1])                
            
            if (bigWig):
                i = 0
                for Bam in listeBam:
                    createBedgraph(ListeFenetreE1aE2ouE2aE1,1+i,Bam.split("/")[-1].split(".")[0]+".bedgraph",chrName)
                    i+=1
            ListeFenetreE1aE2ouE2aE1 = separeList(ListeFenetreE1aE2ouE2aE1,cordVp)
            ListeFenetreE1aE1 = separeList(ListeFenetreE1aE1,cordVp)
            ListeFenetreConsecutif = [createFenetreConsecutif(ListeFenetreE1aE1[0],nbrConsecutif,maxDist),createFenetreConsecutif(ListeFenetreE1aE1[1],nbrConsecutif,maxDist)]
            ListeFenetreE1aE2ouE2aE1[1].reverse()
            ListeFenetreConsecutif[1].reverse()
            ListeFenetreE1aE1[1].reverse()
            ListeFenetreSeparer = [SepareInfluenceVp(ListeFenetreConsecutif[0]),SepareInfluenceVp(ListeFenetreConsecutif[1])]
            regroupListe = []
            regroupListe += main_stat(ListeFenetreE1aE2ouE2aE1[0],nomChroLocal,alphaFDR = alphaFDR,basename = "downE1aE2",code = 1)
            regroupListe += main_stat(ListeFenetreE1aE2ouE2aE1[1],nomChroLocal,alphaFDR = alphaFDR,basename = "upE1aE2",code = 1)
            regroupListe += main_stat(ListeFenetreE1aE1[0],nomChroLocal,alphaFDR = alphaFDR,basename = "downE1aE1",code = 10)
            regroupListe += main_stat(ListeFenetreE1aE1[1],nomChroLocal,alphaFDR = alphaFDR,basename = "upE1aE1",code = 10)
            regroupListe += main_stat(ListeFenetreConsecutif[0],nomChroLocal,alphaFDR = alphaFDR,basename = "downConsecutif",code = 100)
            regroupListe += main_stat(ListeFenetreConsecutif[1],nomChroLocal,alphaFDR = alphaFDR,basename = "upConsecutif",code = 100)
            regroupListe += main_stat(ListeFenetreSeparer[0][-1],nomChroLocal,alphaFDR = alphaFDR,basename = "downSeparer",code = 1000)
            regroupListe += main_stat(ListeFenetreSeparer[1][-1],nomChroLocal,alphaFDR = alphaFDR,basename = "upSeparer",code = 1000)
            regroupListe2 = sorted(regroupListe, key=lambda x: x[1])
            regroupFenetre(regroupListe2)
            
            for i in range(nbpks):
                detectionPeakZscore(ListeFenetreSeparer[1][i],nomChroLocal)
                detectionPeakZscore(ListeFenetreSeparer[0][i],nomChroLocal)
        elif(allChrom):
            try:
                listeSite1,listeSite2 = extractSiteToGenome(sequenceLocal.upper(),enzyme1,enzyme2)
                print(nomChroLocal)
                ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1 = create2Fenetre(listeSite1,listeSite2,25)
                
                ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1,chrName = filterRead(listeBam[0],ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1,nomChroLocal,mapQ)
                for fileBam in listeBam[1:]:
                    ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1,chrName = filterRead(fileBam,ListeFenetreE1aE2ouE2aE1,ListeFenetreE1aE1,nomChroLocal,mapQ)                
                ListeFenetreConsecutif = createFenetreConsecutif(ListeFenetreE1aE1,nbrConsecutif,maxDist)
                if (bigWig):
                    for i,bam in enumerate(listeBam):
                        createBedgraph(ListeFenetreE1aE2ouE2aE1,1+i,bam.split("/")[-1].split(".")[0]+".bedgraph",chrName)
                if (len(ListeFenetreConsecutif)>=1):
                    detectionPeakZscore(ListeFenetreConsecutif,nomChroLocal)
            except:
                print(nomChroLocal)
            #"analyse classique"
            #print(nomChroLocal)
            #print(chrViewPoint.replace("chr","").replace("Chr",""))
            #print("le chromosome "+nomChroLocal+ " ne possedent pas assez de site de restriction pour former un fragment valide et ou n'est pas présent dans les bam")
        nomChroLocal = line.strip()
        sequenceLocal = ""
    else:
        sequenceLocal += line.strip()