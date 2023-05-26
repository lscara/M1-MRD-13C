# Free and 

Ces scripts utilisent R pour effectuer des analyses sur des données de cytométrie en flux normalisé.

## Configuration requise

- R version 4.X.X ou supérieure
- Packages R requis : `flowCore`, `sp`, `Biobase`, `flowViz`,` flowSOM`,`permute`,`igraph`,`rgl`,`GEOmap`,`Mass`,`ade4`

## Installation

1. Assurez-vous d'avoir R installé sur votre système.
2. Installez les packages R requis en exécutant la commande suivante dans R :
         install.packages(c("flowCore", "sp", "Biobase", "flowViz","flowSOM","permute","igraph","rgl","GEOmap","Mass","ade4))
         
## Utilisation

1. Téléchargez ou clonez le script sur votre machine.
2. Ouvrez le scripts R dans RStudio ou tout autre environnement R de votre choix.
3. Modifiez les paramètres du script selon vos besoins, cette option peut être intéractive.
4. Exécutez le scripts et sélectionnez vos fichiers pour effectuer le FlowSOM free et le FlowSOM Frozen.
5. Lors de la selection du dossier assurez vous que cela soit le dossier qui contient comme sous dossier BSD result.

## Sortie

Le script génère un fichiers FCS (Flow Cytometry Standard) contenant le résultat des flowSOM et de l'ACP. 
Le fichier de sortie sera enregistré dans le répertoire courant dans un nom de format prédéfini.
Assurez-vous de spécifier le chemin de sortie approprié dans les scripts si vous souhaitez changer l'emplacement du fichier de sortie.

## Avertissement

Ces scripts sont fournis à titre d'exemple et nécessitent des données spécifiques pour fonctionner correctement. 
Veuillez vous assurer d'adapter le code à vos propres données et de comprendre les implications des analyses effectuées.