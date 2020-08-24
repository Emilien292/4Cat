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
"""@package docstring
Documentation for this module.
 
More details.
"""
def createBedgraph(listeFenetre,numeroBam,basename,chrName, repertoire = "resultats"):
    """! Documentation pour la fonction createBedgraph
        cette fonction construit un bedgraph provenant d'une liste de fenetre. la taille des fragments varies. Chaque fenetre est un fragment.
        Dans le Cas d'un fichier existant. on ecrit a la suite du fichier. Cela permet d'eviter de mettre en memoire chaque liste de fenetre d'un chromosome.

        @param listeFenetre une liste de fenetre
        @param numeroBam int corespondant au numero du bam pour connaitre la valeur a utiliser dans la fenetre.
        @param basename string choisissant le nom du fichier que l'on va crée
        @param chrName string correspondant au nom du chromosome correspondant a la fenetre
        @param repertoire string correspondant au repertoire ou le fichier va s'ecrire par defaut: resultats

        @return pas de retour juste la création du fichier en mode ecriture "a"
    """
    #print(repertoire + "/" + basename)
    bg = open(repertoire + "/" + basename,"a")   
    for fenetre in listeFenetre[1:]:
        bg.write(chrName + '\t' + str(fenetre[0]) + '\t' + str(fenetre[1]) + '\t' +str(fenetre[1+numeroBam])+'\n')
    bg.close()

def regroupFenetre(listeListeFenetre,pourcentage = 25, repertoire = "resultats",basename = "regroupsPeak.bed"):
    """! Documentation pour la fonction regroupFenetre
        cette foanction permet de créer a partir de toute les liste de fenetre un fichier bed regroupant les resultats en ajoutant le fdr moyen
        pour chque peak trouver et un code permettant de retrouver la type de fenetre d'ou provient le peak. le code est un int de 4 chiffre inspirer d'un code bianire correspondant a chaque type de fenetre.
        1 présent 0 absent 

        @param listeListeFenetre une liste de liste de fenetre.
        @param pourcentage int corespondant au pourcentage de chevauchement entre 2 peak pour interpreter cela comme une intersection
        @param basename string choisissant le nom du fichier que l'on va crée
        @param repertoire string correspondant au repertoire ou le fichier va s'ecrire par defaut: resultats

        @return pas de retour juste la création du fichier en mode ecriture "w"
    """
    nouvelleFenetre = []
    chr = listeListeFenetre[0][0]
    
    for i in range(len(listeListeFenetre)):
        start = listeListeFenetre[i][1]
        end = listeListeFenetre[i][2]
        fdr = listeListeFenetre[i][3]
        code = listeListeFenetre[i][4]
        if(i != len(listeListeFenetre)-1):
            j = i+1
            if(end > listeListeFenetre[j][1]):
                end =  listeListeFenetre[j][1]
        else:
            j = i
        if(len(nouvelleFenetre) >= 1):
            if(nouvelleFenetre[-1][2] < start):
                nouvelleFenetre.append([chr,start,end,fdr,code])
            elif (nouvelleFenetre[-1][2] < end):
                nouvelleFenetre.append([chr,nouvelleFenetre[-1][2],end,fdr,code])
        else:
            nouvelleFenetre.append([chr,start,end,fdr,code])
    f = open(repertoire+'/'+basename,"w")
    for fenetre in nouvelleFenetre:
        #print(fenetre)
        start = fenetre[1]
        end = fenetre [2]
        code = fenetre[4]
        fdr = fenetre[3]
        nbrValidation = 1
        for tmpFenetre in listeListeFenetre:
            if(tmpFenetre[1] >= fenetre[2]):
               break
            elif(tmpFenetre[2]>fenetre[1]):
                #il y a intersection on verifie le pourcentage d'intersection
                if(tmpFenetre[2]<fenetre[2]):
                    longueurIntersection = tmpFenetre[2]-fenetre[1]
                else:
                    longueurIntersection = fenetre[2]-tmpFenetre[1]
                F1 = ((fenetre[2]-fenetre[1])/longueurIntersection) * 100
                F2 = ((tmpFenetre[2]-tmpFenetre[1])/longueurIntersection) * 100
                if((F1 >= pourcentage) or (F2 >= pourcentage)):
                    fdr += tmpFenetre[3]
                    nbrValidation += 1
                    if(tmpFenetre[4]==1):
                        if (code%10 < 1):
                            code += tmpFenetre[4]
                    elif((tmpFenetre[4]==10)):
                        if (code % 100 < 10):
                            code += tmpFenetre[4]
                    elif((tmpFenetre[4]==100)):
                        if (code % 1000 < 100):
                            code += tmpFenetre[4]
                    elif((tmpFenetre[4]==1000)):
                        if (code  < 1000):
                            code += tmpFenetre[4]
        bed = str(chr)+'\t'+str(start)+'\t'+str(end)+'\t'+str(fdr/nbrValidation)+'\t'+str(code)+" "+str(nbrValidation)+'\n'
        f.write(bed)
    #print(listeListeFenetre[-1])
    f.close()