

library("GEOmap")
library(rgl)
library(MASS)
library(flowCore)
library(FlowSOM)
library(igraph)
library(permute)
library(tcltk)
library(Biobase)
library(flowViz)

###selec NBM file
my.nNBM.files <- choose.files(default = "", 
                              caption = "Select the nNBM fcs",
                              multi = TRUE, filters = NULL, index = 1)
###determine th parameter
DF.Field <- data.frame(INPUT=rep("   ",32), 
                       var2=rep("  ",32), 
                       Parameters.nb=rep("   ",32), 
                       Parameter.name=rep("  ",32), 
                       stringsAsFactors = FALSE)

DF.Field [1,1] <- "Enter below the nmNBM Merge file name: "
DF.Field [2,1] <- "FSOM64_nmNBM LAMA"

DF.Field [4,1] <- "Enter below the final number of the nm NBM file: "
DF.Field [5,1] <- 4000000


DF.Field [7,1] <- "Enter below the set seed you want to use: "
DF.Field [8,1] <- "1,2,3,4"

DF.Field [10,1] <- "Enter below the FlowSOM parameter number you want to use: "
DF.Field [11,1] <- "2,3,4,5,6,7,8,9,10,11,12,13,14,15"

DF.Field [13,1] <- "Enter below the number of FlowSOM nodes you want to set: "
DF.Field [14,1] <- "100"

DF.Field [16,1] <- "Enter below the final number of clusters you want to get: "
DF.Field [17,1] <- "5"

DF.Field [19,1] <- "Enter below the dispersion coefficient you want to set: "
DF.Field [20,1] <- "60"

DF.Field [22,1] <- "Enter below the TIME parameter number: "
DF.Field [23,1] <- "17"

DF.Field [25,1] <- "Enter below the CD45 parameter number: "
DF.Field [26,1] <- "12"

DF.Field [28,1] <- "Enter below the SSC parameter number: "
DF.Field [29,1] <- "2"

DF.Field [31,1] <- "Enter below the nmNBM Best Set seed file name: "
DF.Field [32,1] <- "FSOM64 nmNBM LAMA Best set seed FROZEN Screening FlowSOM"

flowSet.nNBM <- read.flowSet(my.nNBM.files)
ben <- (flowSet.nNBM[[1]])
ben.1 <- ben@parameters@data$name
ben.2 <- names(ben.1)
DF.Field [1:length(ben.2),3] <- ben.2
DF.Field [1:length(ben.2),4] <- ben.1

DF.Field <- edit(name = DF.Field)##changement manuel de parametre

sample.name <- (DF.Field [2,1])
nmNBM.cells.final.nb <- as.numeric(DF.Field [5,1])

FlowSOM.parameter.1 <- DF.Field [11,1]
FlowSOM.parameter <- as.numeric(unlist(strsplit(FlowSOM.parameter.1,",")))

permutation.1 <- DF.Field [8,1]
permutation <- as.numeric(unlist(strsplit(permutation.1,",")))
nb.nodes.fsom <- as.numeric(DF.Field [14,1])
nb.clusters.fsom <- as.numeric(DF.Field [17,1])
dispersion.coef <- as.numeric(DF.Field [20,1])

CD45.parameter <- as.numeric(DF.Field [26,1])
SSC.parameter <- as.numeric(DF.Field [29,1])
TIME.parameter <- as.numeric(DF.Field [23,1])

sample.name <- DF.Field [32,1]

nmNBM.csv <- fsApply(flowSet.nNBM, exprs)

####rÃ©duire alÃ©atoirement le nombre de cellules de nmNBM

if(nrow(nmNBM.csv) > nmNBM.cells.final.nb){
  nmNBM.cells.final.nb.index <- sample(nrow(nmNBM.csv), nmNBM.cells.final.nb)
  nmNBM.csv <- nmNBM.csv[nmNBM.cells.final.nb.index,]
}

## construction du fichier fcs
nb.parameter <- ncol(nmNBM.csv)

csv.DF <- as.data.frame(nmNBM.csv)

csv.DF.numeric <- as.matrix(csv.DF)


##selection du folder ou ranger les données
output.folder.final <- paste0(choose.dir(), "/", "BSD Results" ,"/")
dir.create(output.folder.final)

setwd(output.folder.final)



##creation du folder 
tempo.1 <- unlist(strsplit(output.folder.final, "/"))
tempo.2 <- length(tempo.1)
tempo.3 <- paste(tempo.1[1:(tempo.2-2)], collapse = "/")
folder.nmNBM <- paste0(tempo.3,"/", "4-Files NBM Merge Norm/")
dir.create(path= folder.nmNBM )

tempo.4 <- paste0(tempo.1[1:(tempo.2-2)], collapse = "/")
tempo.5 <- paste0(tempo.4, "/","nmNBM.folder.name.RDS")
saveRDS(folder.nmNBM, tempo.5)

##creation du fichier rds pour aller rechercher les infos
saveRDS(output.folder.final, paste0(tempo.4,"/", "Best_set_seed_folder_location.RDS"))

list.nmNBM <- list.files(folder.nmNBM, include.dirs = F, full.names = T, recursive = T)

file.remove(list.nmNBM)

##preparation du fichier fcs
cParams<-ncol(csv.DF.numeric)
nb.cells <- cEvents <-nrow(csv.DF.numeric)
csv.File2be.Converted <- csv.DF.numeric

fcs.description.name <- colnames(csv.DF.numeric)
fcs.description.desc <- colnames(csv.DF.numeric)


ff <- flowFrame(csv.File2be.Converted)
linear <- ff

table.to.be.computed <- linear@exprs
table.to.be.used.as.fcs.file <- linear@exprs

ff <- linear

fcs.description.name <- colnames(table.to.be.used.as.fcs.file)
fcs.description.desc <-  colnames(table.to.be.used.as.fcs.file)


######FlowSOM calculation

colsToCluster <- FlowSOM.parameter 
ff <- flowFrame(table.to.be.computed)

all.permutation <- permute::allPerms(permutation)
all.permutation <- as.data.frame(all.permutation)

nb.permutation <- nrow(all.permutation)
all.permutation[nb.permutation+1,] <- permutation

x.dim <- as.numeric(round(sqrt(nb.nodes.fsom)))
y.dim <- as.numeric(round(sqrt(nb.nodes.fsom)))

ben=1
for(ben in 1 : (nb.permutation+1)){
  permutation.1 <- paste0(all.permutation[ben,], collapse = "")
  permutation.2 <- as.numeric(permutation.1)
  
  name.file.flowSOM.fixed.res <- paste0(output.folder.final, "/", "essai",permutation.2,"_fSOM_fixed_res_NBM_SSC.Rdata")
  
  set.seed(permutation.2)
  fSOM.frozen_res <- FlowSOM(ff, compensate=FALSE, transform=FALSE, 
                             scale=F,  
                             colsToUse = colsToCluster, 
                             xdim= x.dim, 
                             ydim= y.dim, 
                             nClus = nb.clusters.fsom, 
                             rlen= 10, 
                             mst = 1, 
                             distf = 2)
  print(name.file.flowSOM.fixed.res)
  save(fSOM.frozen_res, file= name.file.flowSOM.fixed.res)
 
  
  fSOM <- fSOM.frozen_res


  nodes.mapping <- fSOM[["map"]][["mapping"]][,1]
  
  layout <- fSOM[["MST"]][["l"]]
  ###construction MST
  nb.cells <- nrow(table.to.be.computed)
  node.fSOM <- vector()
  node.fSOM.x <- NULL
  node.fSOM.y <- NULL
  
  layout.fcs <- fSOM[["MST"]][["l"]]
  mapping.fcs <- fSOM[["map"]][["mapping"]][,1]
  
  
  vertex.size.fcs <-fSOM[["map"]][["pctgs"]]*500
  
  nb.cluster <- (fSOM[["map"]][["xdim"]]) * (fSOM[["map"]][["ydim"]])
  
  nb.nodes <- nb.cluster
  nb.points.per.line <- nb.cells /nb.nodes
  
  for (i in 1:nb.cluster) {
    tempo <-i
    tempo_indices <- which(mapping.fcs==i)
  
    random.nb <- runif(n=tempo,min = 0, max = 1 )
  
    random.radius <-  random.nb * vertex.size.fcs[i]/2
    random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
    
    node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
    node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
  }
  
   ###mise à l'echelle
  xmin.graph <-  par("usr")[1]
  xmax.graph <- par("usr")[2]
  
  ymin.graph <- par("usr")[3]
  ymax.graph <- par("usr")[4]
  
  node.fSOM.x.min <- xmin.graph
  node.fSOM.x.max <- xmax.graph
  
  node.fSOM.y.min <- ymin.graph
  node.fSOM.y.max <- ymax.graph
  
  node.fSOM.x.coef <- ((1024*1024)-100000)/(node.fSOM.x.max-node.fSOM.x.min)
  node.fSOM.y.coef <- ((1024*1204)-100000)/(node.fSOM.y.max-node.fSOM.y.min)
  
  node.fSOM.x.ind.term <- 100000 - (node.fSOM.x.coef * node.fSOM.x.min)
  node.fSOM.y.ind.term <- 100000- (node.fSOM.y.coef * node.fSOM.y.min)
  
  node.fSOM.x.scaled <- (node.fSOM.x.coef * node.fSOM.x) + node.fSOM.x.ind.term
  node.fSOM.y.scaled <- (node.fSOM.y.coef * node.fSOM.y) + node.fSOM.y.ind.term
  
  
  matrix.node.fSOM <- cbind(node.fSOM.x, node.fSOM.y)
  matrix.node.fSOM.scaled <-igraph::norm_coords(matrix.node.fSOM, xmin=100000, xmax=2^20-100000, ymin=100000, ymax=2^20-100000)
  colnames(matrix.node.fSOM.scaled) <- c("FlowSOM.1","FlowSOM.2")
  
  node.fSOM.x.scaled <- matrix.node.fSOM.scaled[,1]
  node.fSOM.y.scaled <- matrix.node.fSOM.scaled[,2]
  ####preparation FCS
  col.of.interest <- 1:ncol(table.to.be.used.as.fcs.file)
  col.names.fcs <- colnames(table.to.be.used.as.fcs.file)
  
  fcs.Gated.CD45.comp.csv  <- table.to.be.used.as.fcs.file
  
  csv.1 <- cbind(fcs.Gated.CD45.comp.csv, node.fSOM.x.scaled, node.fSOM.y.scaled)
  
  
  colnames(csv.1) <- c(colnames(table.to.be.used.as.fcs.file),
                       "FlowSOM Nodes.x", 
                       "FlowSOM Nodes.y")
  
  
  
  csv.merge <- csv.1
  
  fcs.description.name.1 <- colnames(table.to.be.used.as.fcs.file) #name
  
  
  
  ####Define the complete output fcs file name
  
  
  colnames(csv.merge) <- c(fcs.description.name.1,"FlowSOM Nodes.x", "FlowSOM Nodes.y") 

  dta <- csv.merge
  
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
  
  write.FCS(ff,paste0(output.folder.final, "/", "essai",permutation.2,"_fSOM_fixed_res_NBM_SSC.FCS"))
}
