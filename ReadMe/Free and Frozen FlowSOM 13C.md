# Documentation Free and Frozen flowSOM

## Selection du fichier DIAG, FU, NBM

On selectionne dans l'ordre les fichier fcs Diag, FU, NBM, puis on les lit. on choisit aussi le dossier ou se trouve le résultat du bestSEED

## Selection des paramètre

Un tableau de données est créer ou l'on rentre tous les paramètre que l'on doit spécifier.On entre aussi le résultat du best Seed que l'on souhaite.

## Create merge

On fusionne les 3 fichier FCS.

## Output folder

On définie le dossier de sortie et on se place dedans

## Identification

On identifie chaque data set par leur TAG

## Frozen integration

On intègre le flowSOM créer par le bestSeed et on injecte individuellement chaque FCS dans le flowSOM. A partir de cela on construit un MST par flowSOM et on récupère dans un dataFrame les donné du MST.

## Free flowSOM

On crée un flowSOM free sur les données fusionner et on determine un MST.

## selection des noeud d'interêt

On compte le nombre de cellule dans chaque noeud ou des cellule du diagnostique sont présente et on détermine quelle noeuds à selectionner à partir des paramètre définis dans la selection des paramètre. Ensuite pour chaque cellule on détermine dans un dataframe si elle fait partie d'un noeuds d'intéret ou non.

## ACP

On détermine une ACP sur les noeud pour selectionner les noeud d'inetérets

## Construction du fichier FCS

Dans un seul fichier FCS on spécifie on intègre au data frame contennant les données fusionner, les coordonné du MST fussioné, les noeuds d'intéret, les noeuds, leur meta cluster et les donneé de l'ACP. Puis on integre les données du MST Frozen. On sauve le tout dans un fichier FCS.