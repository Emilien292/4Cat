# 4Cat
outils d'analyse de données 4C

## Description
L'outils permet actuellement uniquement de faire de la detection de peak à partir de bam triée avec son index sur tous les chromosomes (trans analyse) et spécifiquement le chromosome que l'on appelle le viewpoint. L'outils à comme but de laissez chaque personne insérer une méthode d'analyse statistique  

### L'architecture de l'outils. 
L'outils est séparer en plusieurs script. Chaque script aura un but particulier.

- moduleBam: C'est le script qui s'occupe de la modification, utilisation des reads 
- utils: C'est le scripts qui permmettent de crée des fichiers suplémentaire pour le controle qualiter, et simplifier la lecture des résultats
- Zscore: C'est le script qui s'occupent de calculer le zscore pour la detection de peak dans les liste ou l'influence du viewpoint est faible
- stat_function: c'est le script qui s'ocuppent de la detection de peak sur les listes ou l'influence du viewpoint est présent
- 4Cat: C'est le script qui regroupent toute l'analyse et le lance. Il récupere les options choisis par les utilisateur et lance en fonction des options choisi.

## Mise en place de l'outils

### Besoin
l'outils utilise python >= 3.5 en utilisant les packages suivants:
- joblib 0.12.0
- matplotlib 3.0.3
- rupture 1.4.5
- sklearn >= 0.0
- numpy 1.18.5 >= 1.11.5
- spicy >= 0.17.0 

il suffira par la suite d'avoir tous les scripts au meme endroits.

## Utilisation
l'outils utilise des options obligatoire comme une liste des fichier bam, le chemin du genome utiliser, les enzymes, le chromosome du viewpoint. Le lancement de l'outils peut etre fait par cette commande.
python3 4Cat.py -e1 CATG -e2 GATC --listeBam data/*A.bam -g dm6.fa -chr chr2R 
python3 4Cat.py --help permet d'avoir une aide plus détailler de l'outils
