# Documentation bestSeed

## Selection du fichier NBM
On selection le fichier fcs de NBM, puis on le lit

## Selection des paramètre

Un tableau de données est créer ou l'on rentre tous les paramètre que l'on doit spécifier.On entre aussi les cannaux de permutation.

## Réduction du nombre de cellue

On réduis le nombre de cellule.

## Selection du folder

On choisis le dossier de stockage et on cree les ficchier et dossier de stockage. On cree le data frame contenant les données du fichier FCS

## Calcul flowSOM

Pour chaque canal de permutation on détermine un flowSOM, on construit un MST avec des coordonné x et y remise à échelle. On crée ensuite le fichier FCS de sortie en fusionnant le MST avec le data frame du fichier NBM.

