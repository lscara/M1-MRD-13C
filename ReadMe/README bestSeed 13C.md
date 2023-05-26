# BestSeed

Ces scripts utilisent R pour effectuer des analyses sur des données de cytométrie en flux normalisé.

## Configuration requise

- R version 4.X.X ou supérieure
- Packages R requis : `flowCore`, `sp`, `Biobase`, `flowViz`,` flowSOM`,`permute`,`igraph`,`rgl`,`GEOmap`,`Mass`

## Installation

1. Assurez-vous d'avoir R installé sur votre système.
2. Installez les packages R requis en exécutant la commande suivante dans R :
         install.packages(c("flowCore", "sp", "Biobase", "flowViz","flowSOM","permute","igraph","rgl","GEOmap","Mass"))
         
## Utilisation

1. Téléchargez ou clonez le script sur votre machine.
2. Ouvrez le script R dans RStudio ou tout autre environnement R de votre choix.
3. Modifiez les paramètres des scripts selon vos besoins, cette option peut être intéractive.
4. Exécutez les scripts et sélectionnez vos fichiers pour déterminer les meilleurs parmètre du flowSOM.

## Sortie

Le script génère des fichiers FCS (Flow Cytometry Standard) contenant les résultats des analyses ainsi que des fichier RData. 
Les fichiers de sortie seront enregistrés dans un nouveau répertoire :BSD result.
Assurez-vous de spécifier le chemin de sortie approprié dans les scripts si vous souhaitez changer l'emplacement du fichier de sortie.

## Avertissement

Ces scripts sont fournis à titre d'exemple et nécessitent des données spécifiques pour fonctionner correctement. 
Veuillez vous assurer d'adapter le code à vos propres données et de comprendre les implications des analyses effectuées.

