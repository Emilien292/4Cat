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
from scipy.stats import zscore,norm
#p_values = scipy.stats.norm.sf(abs(z_scores))
def detectionPeakZscore(listePeakPotentielle,chromosome, FDR = 0.05, basename = "Zscore", repertoire = "resultats", methods = "intersection_majority",corrector = 1):
    ListeReplicats = []
    nbrReplicat = len(listePeakPotentielle[0]) - 2
    for i in range(2,len(listePeakPotentielle[0])):
        replicats = []
        for j in listePeakPotentielle:
            replicats.append(j[i])
        ListeReplicats.append(zscore(replicats))
    for i in range(len(ListeReplicats[0])):
        nbrPvalueSinificative = 0
        for j in range(nbrReplicat):
            #print(norm.sf(abs(ListeReplicats[j][i])))
            if(norm.sf(abs(ListeReplicats[j][i])) <= FDR):
                #print(ListeReplicats[j][i])
                #print(norm.sf(abs(ListeReplicats[j][i])))
                nbrPvalueSinificative +=1
        if (methods == "intersection_majority"):
            if(nbrPvalueSinificative > (nbrReplicat/2)):
                a = 1+1
                #print(listePeakPotentielle[i])
        elif (methods == "union"):
            if(nbrPvalueSinificative >= 1):
                print(listePeakPotentielle[i])
        elif (methods == "intersection_total"):
            if (nbrPvalueSinificative == nbrReplicat):
                print(listePeakPotentielle[i])