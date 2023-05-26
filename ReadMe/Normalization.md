# Normalization

Ces scripts utilisent R pour effectuer des analyses sur des données de cytométrie en flux.

## Configuration requise

- R version 4.X.X ou supérieure
- Packages R requis : `flowCore`, `sp`, `Biobase`, `flowViz`

## Installation

1. Assurez-vous d'avoir R installé sur votre système.
2. Installez les packages R requis en exécutant la commande suivante dans R :  
         install.packages(c("flowCore", "sp", "Biobase", "flowViz"))
         
## Utilisation

1. Téléchargez ou clonez les scripts sur votre machine.
2. Ouvrez les scripts R dans RStudio ou tout autre environnement R de votre choix.
3. Modifiez les paramètres des scripts selon vos besoins.
4. Exécutez les scripts et sélectionnez vos fichiers pour effectuer les analyses sur les données de cytométrie en flux.
5. Les scripts sont interactifs :  
   5.1. Assurez-vous de bien sélectionner à l'aide de la souris les cellules indiquées dans le titre du graphe biparamétrique.  
         5.1.1. Il s'agit des lymphocytes dans les scripts "Normalization NBM 13C.R" et "Normalization FU 13C.R".  
         5.1.2. Il s'agit des lymphocytes puis des cellules blastiques dans le script "Normalization DIAG 13C.R".  
   5.2 Assurez-vous de bien sélectionner à l'aide de la souris les pics négatifs dans les histogrammes de densité.  
       Il s'agit des pics ayant la valeur en abscisse la moins élevée, peu importe leur densité.

## Sortie

Les scripts génèrent des fichiers FCS (Flow Cytometry Standard) contenant les résultats des analyses.  
Les fichiers de sortie seront enregistrés dans le répertoire courant dans un nom de format prédéfini.  
Assurez-vous de spécifier le chemin de sortie approprié dans les scripts si vous souhaitez changer l'emplacement du fichier de sortie.

## Avertissement

Ces scripts sont fournis à titre d'exemple et nécessitent des données spécifiques pour fonctionner correctement.  
Veuillez vous assurer d'adapter le code à vos propres données et de comprendre les implications des analyses effectuées.
