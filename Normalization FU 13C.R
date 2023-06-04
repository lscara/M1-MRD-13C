# Normalization for Flow Cytometry Data
# Loredan Scarabello
# 26-05-2023
# https://github.com/lscara
# Script adaped from https://github.com/FlowFun/M1-MRD-en-13C/blob/main/Legacy%2010C/Normalization%20FU%2010C%20navios.R by Benoit DUPONT and https://gist.github.com/yannabraham/c1f9de9b23fb94105ca5


#############################

#Télécharger les bibliothèques
library(flowCore)
library(sp)
library(Biobase)
library(flowViz)

#Demander la date
date <- "JJ-MMM-YYYY <- Exemple : 02-OCT-2023"
date <- edit(date)

#Demander le nom du patient
nom_patient <- "Prenom NOM"
nom_patient <- edit(nom_patient)

#Lecture du fichier
logicle <- readRDS(file.choose())

#Récupération des colonnes d'intérêt
col_names <- colnames(logicle)[grep("-A$", colnames(logicle))]
new_logicle <- logicle[, col_names]

TAG.channel <- 500000 #DIAG=250000 // FU=500000 // NBM=750000
final.cell.number <- 500000 #DIAG=30000 // FU=500000 // NBM=50000
CD45.parameter.nb <- 12
SSCINT.parameter.nb <-  2

#Mise à l'échelle des paramètres "FSC" et "SSC"
max_range_FS <- max(new_logicle[,"FSC-A"])
max_range_SS <- max(new_logicle[,"SSC-A"])
new_logicle[,'FSC-A'] <- new_logicle[,'FSC-A'] / max_range_FS
new_logicle[,'SSC-A'] <- new_logicle[,'SSC-A'] / max_range_SS

#Nommer les paramètres
parameter.name <- c("FSC INT", 
                    "SSC INT",
                    "CD7 FITC B-525",
                    "CD13 PE B-585",
                    "HLA-DR ECD B-610",
                    "CD33 PC5.5 B-690",
                    "CD38 HB7 B-780",
                    "CD34 APC R-660",
                    "CD56 APC700 R-712",
                    "CD19 APC750 R-780",
                    "CD117 SNv428 V-450",
                    "CD45 KO V-525",
                    "CD133 BV605 V-610",
                    "CD54 BV650 V-660",
                    "CD36 BV786 V-780") 

colnames(new_logicle) <- parameter.name

#Définir des valeurs expérimentales pour chaque paramètre 
reference.neg.pop.mode <- c(0.2,
                            0.2,
                            0.2, 
                            0.2, 
                            0.2,
                            0.2,
                            0.2,
                            0.2,
                            0.2,
                            0.2,
                            0.2,
                            0.7,
                            0.2, 
                            0.2, 
                            0.2)

target.column <- 1:15

#Affichage du graphe biparamétrique ou dot plot
windows()
plot(new_logicle [,CD45.parameter.nb ],
     new_logicle [,SSCINT.parameter.nb ],
     col = densCols(new_logicle [,CD45.parameter.nb ],
                    new_logicle [,SSCINT.parameter.nb ],
                    colramp = colorRampPalette(c("blue2", "green2","red2", "yellow"))),
     pch = 19, cex=0.1, xlim=c(0,1),
     main = "Draw a region on lymphocytes population")

options(locatorBell = FALSE)
coord <- locator(500, type="o", col="red", pch=3)
lympho.pop <- point.in.polygon (point.x = new_logicle[,CD45.parameter.nb], 
                                point.y = new_logicle[,SSCINT.parameter.nb],
                                pol.x = coord$x, 
                                pol.y = coord$y, 
                                mode.checked = FALSE)
dev.off()

#Récupération des événements 
ff.logicle.gated.CD45 <- new_logicle [lympho.pop==1,]
negative.population.mode <- c()

#Affichage des histogrammes de densité 
for(i in target.column) {
  windows()
  x1 <- ff.logicle.gated.CD45[,i]
  d1 <- density(x1)
  plot(d1, xlab= colnames(new_logicle)[i], main="Click on negative mode peak", cex.main=0.7, xlim=c(0,1))
  coord.1 <- locator(1, type = "l", col="red")
  negative.population.mode[i] <- coord.1$x
  abline(v = negative.population.mode[i],col="red",lwd=1,lty=5)
  dev.off()
}

#Normalisation
csv.logicle.corrected <- lapply(1:15, function(i)
  new_logicle[,i] + 
    (reference.neg.pop.mode[i] - negative.population.mode[i]))

csv.logicle.corrected <- do.call(cbind, csv.logicle.corrected)

#Récupération d'un échantillon de taille définie
if(nrow(csv.logicle.corrected) >final.cell.number){
  sample_indices <- sample(1:nrow(csv.logicle.corrected), final.cell.number)
  csv.logicle.corrected <- csv.logicle.corrected[sample_indices,]
}

#Récupération de la colonne "TIME"
col_names2 <- colnames(logicle)[grep("TIME", colnames(logicle))]
TIME <- logicle[, col_names2]
TIME <- head(TIME, nrow(csv.logicle.corrected))

#Renommage des colonnes 
tempo <- colnames(new_logicle)
colnames(csv.logicle.corrected) <- c(tempo)

#Création d'une colonne "TAG" 
TAG.fcs <- rep(TAG.channel, nrow(csv.logicle.corrected))

#Fusion des différents objets
csv.DF.numeric <- cbind(csv.logicle.corrected, TAG.fcs, TIME)

#Création du FCS
dta <- csv.DF.numeric

meta <- data.frame(name=colnames(dta),
                   desc=paste(colnames(dta))
)
meta$Range <- apply(apply(dta, 2, FUN='range'), 2, FUN='diff')
meta$minRange <- apply(dta, 2, min)
meta$maxRange <- apply(dta, 2, max)

ff <- new("flowFrame",
          exprs=dta,
          parameters=AnnotatedDataFrame(meta)
)

#Nom du fichier de sortie
file_name <- paste(nom_patient, date, "FU", ".fcs", sep = "_")
file_name <- gsub("_\\.", ".", file_name)
write.FCS(ff, file = file_name)
