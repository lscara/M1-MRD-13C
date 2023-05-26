library("GEOmap")
library(rgl)
library(MASS)
library(flowCore)
library(FlowSOM)
library(igraph)

library(tcltk)
library(sp)
library(grDevices)
library(ade4)
library(Biobase)
library(flowViz)

###
output.folder=paste0(choose.dir("",caption="Select the folder where BSD "),"/")

my.nDIAG <- choose.files(default = "", 
                         caption = "Select the DIAG fcs",
                         multi = TRUE, filters = NULL, index = 1)
my.nFU <-  choose.files(default = "", 
                        caption = "Select the FU fcs",
                        multi = TRUE, filters = NULL, index = 1)
my.nmNBM <-  choose.files(default = "", 
                          caption = "Select the nNBM fcs",
                          multi = TRUE, filters = NULL, index = 1)

name.FU<-basename(my.nFU)
patient.name.date<- sub("_FU\\.fcs$", "", name.FU)

###read fcs
ff.nDIAG <- read.FCS(my.nDIAG)
ff.nFU <- read.FCS(my.nFU)
ff.nmNBM <- read.FCS(my.nmNBM)



ben.1 <- ff.nDIAG@parameters@data$name
ben.2 <- names(ben.1)

###determine the parameter

DF.Field <- data.frame(INPUT=rep("   ",44), 
                       var2=rep("  ",44), 
                       Parameters.nb=rep("   ",44),
                       Parameters.name=rep("   ",44),
                       var5=rep("  ",44),
                       best.set.seed.Rdata_list=rep("  ",44), 
                       stringsAsFactors = FALSE)

DF.Field [1,1] <- "Enter below the nDIAG nFU nmNBM Merge file name: "
DF.Field [2,1] <- patient.name.date

DF.Field [4,1] <- "Enter below the best set seed file you want to use: "
DF.Field [5,1] <- "essai1234_fSOM_fixed_res_NBM_SSC.Rdata"


DF.Field [7,1] <- "Enter below the set seed you want to use: "
DF.Field [8,1] <- "1234"

DF.Field [10,1] <- "Enter below the FlowSOM parameter number you want to use: "
DF.Field [11,1] <- "2,3,4,5,6,7,8,9,10,11,12,13,14,15"

DF.Field [13,1] <- "Enter below the number of FlowSOM nodes you want to set: "
DF.Field [14,1] <- "64" #default 64

DF.Field [16,1] <- "Enter below the final number of Metaclusters you want to get: "
DF.Field [17,1] <- "30"

DF.Field [19,1] <- "Enter below the dispersion coefficient you want to set: "
DF.Field [20,1] <- "70"

DF.Field [22,1] <- "Enter below the TIME parameter number: "
DF.Field [23,1] <- "17"

DF.Field [25,1] <- "Enter below the CD45 parameter number: "
DF.Field [26,1] <- "12"

DF.Field [28,1] <- "Enter below the SSC parameter number: "
DF.Field [29,1] <- "2"

DF.Field [31,1] <- "Enter below the TAG parameter number: "
DF.Field [32,1] <- "16"

DF.Field [34,1] <- "Enter below the final number of the nDIAG nFU nmNBM merge file: "
DF.Field [35,1] <- "4000000"

DF.Field [37,1] <- "Enter below the NOI DIAG pct threshold: "
DF.Field [38,1] <- "1"

DF.Field [40,1] <- "Enter below the NOI Ratio threshold: "
DF.Field [41,1] <- "20"

DF.Field [43,1] <- "Enter below the NOI DIAG cell nb threshold: "
DF.Field [44,1] <- "10"


##name the RDATa containning the NBM flowSOM
tempo.1 <- unlist(strsplit(output.folder, "/"))
tempo.2 <- length(tempo.1)
tempo.3 <- paste(tempo.1[1:(tempo.2-1)], collapse = "/")
tempo.4 <- paste0(tempo.3,"/", "Best_set_seed_folder_location.RDS")
folder.NBM <- readRDS(tempo.4)

tempo.5 <- list.files(path = folder.NBM, pattern = "\\.Rdata$",full.names = FALSE)
DF.Field [1:length(tempo.5),6] <- tempo.5

DF.Field [1:length(ben.2),3] <- ben.2
DF.Field [1:length(ben.2),4] <- ben.1

#DF.Field <- edit(name = DF.Field)

#####name the paramater tu use
sample.name <- (DF.Field [2,1])
essai_fSOM_fixed_res_NBM_SSC_Rdata <- DF.Field [5,1]
final.cell.number <- as.numeric(DF.Field [5,1])
CD45.parameter <- as.numeric(DF.Field [26,1])
SSC.parameter <- as.numeric(DF.Field [29,1])
TAG.fcs.parameter <- as.numeric(DF.Field [32,1])
TIME.parameter <- as.numeric(DF.Field [23,1])

FlowSOM.parameter.1 <- DF.Field [11,1]
FlowSOM.parameter <- as.numeric(unlist(strsplit(FlowSOM.parameter.1,",")))

###NOI parameter
nb.nodes.freefsom <- as.numeric(DF.Field [14,1]) 
nb.clusters.freefsom <- as.numeric(DF.Field [17,1])
dispersion.coef <- as.numeric(DF.Field [20,1])  #default=30 how big you want the nodes

set.seed(as.numeric(DF.Field [8,1]))
final.cell.number <- as.numeric(DF.Field [35,1])

DIAG.pct.threshold <- as.numeric(DF.Field [38,1])
Ratio.threshold <- as.numeric(DF.Field [41,1])
DIAG.cell.nb.threshold <- as.numeric(DF.Field [44,1])


###create merge
csv.merge <-rbind(ff.nDIAG@exprs, ff.nFU@exprs, ff.nmNBM@exprs)
if(nrow(csv.merge) >final.cell.number){
  sample_indices <- sample(nrow(csv.merge),final.cell.number)	
  sample_indices <- sample(1:nrow(csv.merge), final.cell.number)
  csv.merge <- csv.merge[sample_indices,]
}

nb.parameter <- ncol(csv.merge)
csv.DF <- as.data.frame(csv.merge)
csv.DF.numeric <- as.matrix(csv.DF)


file.name.nDIAGnFUnmNBM <- sample.name

#####calculate number of events and number of parameters

cParams<-ncol(csv.DF.numeric)
nb.cells <- cEvents <-nrow(csv.DF.numeric)
csv.File2be.Converted <- csv.DF.numeric
ff <- flowFrame(csv.File2be.Converted)

linear <- csv.File2be.Converted
cEvents <- nrow(csv.File2be.Converted)
cParams <- ncol(csv.File2be.Converted)
table.to.be.computed <- linear
table.to.be.used.as.fcs.file <- linear
fcs.description.name <- colnames(csv.DF.numeric)
fcs.description.desc <- colnames(csv.DF.numeric)
##outputfolder
title.frozenfsom <- c("DIAG", "FU1","NBM merge")
sample.name <- paste(sample.name, "_", "FSOM64_fSOM free_fSOM frozen_NOI" ,nb.nodes.freefsom,sep="")

tempo.1 <- unlist(strsplit(output.folder, "/"))
tempo.2 <- length(tempo.1)
tempo.3 <- paste(tempo.1[1:(tempo.2-1)], collapse = "/")
tempo.4 <- paste0(tempo.3,"/", "Best_set_seed_folder_location.RDS")
folder.NBM <- readRDS(tempo.4)
name.file.flowSOM.fixed.res <- paste0(folder.NBM, essai_fSOM_fixed_res_NBM_SSC_Rdata)
table.to.be.computed <- linear

output.folder.final <-  paste0(output.folder,sample.name,"_results/")
print(output.folder.final)
dir.create(path= output.folder.final )
setwd(output.folder.final)

sFilename <-paste(output.folder.final,file.name.nDIAGnFUnmNBM,".fcs", sep="")
sFilename.2 <-paste(file.name.nDIAGnFUnmNBM,".fcs", sep="")





###############################################################################################
#                    identification of each dataset
#For Navios DIAG @250000 / FU @ 500000 / NBM @ 750000
###############################################################################################
##### START TAG Files
DIAG_indices <- which(linear[,TAG.fcs.parameter]==250000)
FU_indices <- which(linear[,TAG.fcs.parameter]==500000)
NBM_indices <- which(linear[,TAG.fcs.parameter]==750000)

csv <- linear
csv.logicle <- linear


par.default <- par() 

ff <- flowFrame(as.matrix(table.to.be.computed))

fcs.description.name <- colnames(table.to.be.used.as.fcs.file)
fcs.description.desc <-  colnames(table.to.be.used.as.fcs.file)


###FlowSOM FROZEN integration
load(file = name.file.flowSOM.fixed.res)

####reference MST plot
fSOM.fixed.NBM <- fSOM.frozen_res
#### Building MST tree nodes from FlowSOM and scaled them on 2^20 
nb.cells <- nrow(table.to.be.computed)
node.fSOM <- vector()
node.fSOM.x <- NULL
node.fSOM.y <- NULL

layout.fcs <- fSOM.fixed.NBM[["MST"]][["l"]]
mapping.fcs <- fSOM.fixed.NBM[["map"]][["mapping"]][,1]


vertex.size.fcs <-fSOM.fixed.NBM[["map"]][["pctgs"]]*500

nb.cluster <- (fSOM.fixed.NBM[["map"]][["xdim"]]) * (fSOM.fixed.NBM[["map"]][["ydim"]])

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


matrix.node.fSOM <- cbind(node.fSOM.x, node.fSOM.y)
matrix.node.fSOM.scaled <- igraph::norm_coords(matrix.node.fSOM, xmin=100000, xmax=2^20-100000, ymin=100000, ymax=2^20-100000)
colnames(matrix.node.fSOM.scaled) <- c("FlowSOM.x","FlowSOM.y")

node.fSOM.x.scaled <- matrix.node.fSOM.scaled[,1]
node.fSOM.y.scaled <- matrix.node.fSOM.scaled[,2]




####inject the data of each dataset in the FlowSOM FROZEN
csv1.0 <- linear[DIAG_indices,] # csv.DIAG.0

csv2.0 <- linear[FU_indices,] # csv.FU.0

csv3.0 <- linear[NBM_indices,] #csv.NBM.0

ff1 <- flowFrame(csv1.0) #ff.DIAG
ff2 <- flowFrame(csv2.0) #ff.FU
ff3 <- flowFrame(csv3.0) #ff.NBM

fSOM1.FROZEN <- NewData(fSOM.fixed.NBM, ff1) #fSOM.FROZEN.DIAG
fSOM2.FROZEN <- NewData(fSOM.fixed.NBM, ff2)#fSOM.FROZEN.FU
fSOM3.FROZEN <- NewData(fSOM.fixed.NBM, ff3)#fSOM.FROZEN.NBM  
csv.nodes <- list()

csv.list <- list(csv1.0, csv2.0, csv3.0)

fSOM.list <- list(fSOM1.FROZEN, fSOM2.FROZEN ,fSOM3.FROZEN) 

for (m in 1:3){
  
  csv.fsom <- csv.list[[m]]
  nb.cells <- nrow(csv.fsom)
  
  
  # **FlowSOM: build the MST layout and Vertex size nodes for the fcs file**
  #############################################################################################
  #                     Build fcs file with embedded MST plot setting                         #
  #############################################################################################
  
  #### Building MST tree nodes from FlowSOM and scaled them on 2^20 
  node.fSOM <- vector()
  node.fSOM.x <- NULL
  node.fSOM.y <- NULL
  
  layout.fcs <- fSOM.list[[m]][["MST"]][["l"]]
  xmin.frozen.NBM <- par("usr")[1]
  xmax.frozen.NBM <- par("usr")[2]
  ymin.frozen.NBM <- par("usr")[3]
  ymax.frozen.NBM <- par("usr")[4]
  
  mapping.fcs <- fSOM.list[[m]][["map"]][["mapping"]]
  vertex.size.fcs <- fSOM.list[[m]][["map"]][["pctgs"]]*500
  nb.cluster <- (fSOM.list[[m]][["map"]][["xdim"]]) * (fSOM.list[[m]][["map"]][["ydim"]])
  nb.nodes <- nb.cluster
  nb.points.per.line <- nb.cells /nb.nodes
  
  for (i in 1:nb.cluster) {
    tempo <- mapping.fcs[mapping.fcs==i]
    tempo_indices <- which(mapping.fcs==i)
    random.nb <- runif(n=tempo,min = 0, max = 1 )
    random.radius <-  random.nb * vertex.size.fcs[i]/2
    random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
    
    node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
    node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
    
  }
  
  
  #############################################################################################
  #                  Building fcs file: MST layout + Vertex Size nodes                        #
  #############################################################################################
  
  matrix.node.fSOM <- cbind(node.fSOM.x, node.fSOM.y)
  matrix.node.fSOM.scaled <- igraph::norm_coords(matrix.node.fSOM, xmin=100000, xmax=2^20-100000, ymin=100000, ymax=2^20-100000)
  colnames(matrix.node.fSOM.scaled) <- c("FlowSOM.x","FlowSOM.y")
  
  node.fSOM.x.scaled <- matrix.node.fSOM.scaled[,1]
  node.fSOM.y.scaled <- matrix.node.fSOM.scaled[,2]
  
  #############################################################################################
  #                     START building the fcs file with new parameters                       #
  #############################################################################################
  
  csv.nodes[[m]] <- cbind(node.fSOM.x.scaled, node.fSOM.y.scaled)
  
  colnames(csv.nodes[[m]]) <- c("FROZEN FlowSOM Nodes.x", 
                                "FROZEN FlowSOM Nodes.y") 
  
}
csv.FROZEN <- do.call(rbind,csv.nodes)


###############################################################################
#FlowSOM Free calculation
###############################################################################
colsToCluster <- FlowSOM.parameter 

x.dim <- as.numeric(round(sqrt(nb.nodes.freefsom)))
y.dim <- as.numeric(round(sqrt(nb.nodes.freefsom)))

fSOM_res <- FlowSOM(ff, compensate=FALSE, transform=FALSE, 
                    scale=T,  
                    colsToUse = colsToCluster, 
                    xdim= x.dim, 
                    #xdim =8,
                    ydim= y.dim, 
                    #ydim = 8,
                    nClus = nb.clusters.freefsom, 
                    rlen= 10, 
                    mst = 1, 
                    distf = 2)

fSOM <- fSOM_res

nodes.mapping <- fSOM[["map"]][["mapping"]]

layout <- fSOM[["MST"]][["l"]]
layout.norm <- norm_coords(layout, 
                           xmin = 100000, xmax = (2^20-100000), 
                           ymin = 100000, ymax = (2^20-100000)-10)
####metaclustering
metacluster <- fSOM 
metacluster <- metaClustering_consensus(fSOM[["map"]][["codes"]], k = nb.clusters.freefsom)
data.metacluster <- metacluster[fSOM[["map"]][["mapping"]][,1]]
data.metacluster <- as.matrix(data.metacluster)
colnames(data.metacluster) <- "Metaclustering Consensus"


######determination MST

target.channel <- FlowSOM.parameter 
marker.name.desc <- fcs.description.desc[target.channel]
marker.name.name <- fcs.description.name[target.channel]


#####START building the fcs file with new parameters - nodes
nb.cells <- nrow(table.to.be.computed)

lDiag<-length(ff.nDIAG@exprs[,1])
lFU<-length(ff.nFU@exprs[,1])+lDiag
lNBM<-length(ff.nmNBM@exprs[,1])+lFU
nb.cells <- nrow(table.to.be.computed)

####parameter node diag
node.fSOM <- vector()
node.fSOM.x <- NULL
node.fSOM.y <- NULL

layout.fcs <- fSOM[["MST"]][["l"]]
mapping.fcs <- fSOM[["map"]][["mapping"]][1:lDiag]
vertex.size.fcs <- fSOM[["map"]][["pctgs"]]*500

nb.cluster <- (fSOM[["map"]][["xdim"]]) * (fSOM[["map"]][["ydim"]])

nb.nodes <- nb.cluster
nb.points.per.line <- nb.cells /nb.nodes


for (i in 1:nb.cluster) {
  tempo <- mapping.fcs[mapping.fcs==i]
  tempo_indices <- which(mapping.fcs==i)
  
  random.nb <- runif(n=tempo,min = 0, max = 1 )
  
  random.radius <-  random.nb * vertex.size.fcs[i]/2
  random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
  
  node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
  node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
  
}

matrix.node.fSOM <- cbind(node.fSOM.x, node.fSOM.y)
matrix.node.fSOM.scaled <- igraph::norm_coords(matrix.node.fSOM, xmin=100000, xmax=2^20-100000, ymin=100000, ymax=2^20-100000)
colnames(matrix.node.fSOM.scaled) <- c("FlowSOM.x","FlowSOM.y")

node.fSOM.x.scaled.Diag <- matrix.node.fSOM.scaled[,1]
node.fSOM.y.scaled.Diag <- matrix.node.fSOM.scaled[,2]

###parameter nodeFU
node.fSOM <- vector()
node.fSOM.x <- NULL
node.fSOM.y <- NULL
layout.fcs <- fSOM[["MST"]][["l"]]
mapping.fcs <- fSOM[["map"]][["mapping"]][(lDiag+1):lFU]
vertex.size.fcs <- fSOM[["map"]][["pctgs"]]*500

nb.cluster <- (fSOM[["map"]][["xdim"]]) * (fSOM[["map"]][["ydim"]])

nb.nodes <- nb.cluster
nb.points.per.line <- nb.cells /nb.nodes


for (i in 1:nb.cluster) {
  tempo <- mapping.fcs[mapping.fcs==i]
  tempo_indices <- which(mapping.fcs==i)
  
  random.nb <- runif(n=tempo,min = 0, max = 1 )
  
  random.radius <-  random.nb * vertex.size.fcs[i]/2
  random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
  
  node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
  node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
  
}


matrix.node.fSOM <- cbind(node.fSOM.x, node.fSOM.y)
matrix.node.fSOM.scaled <- igraph::norm_coords(matrix.node.fSOM, xmin=100000, xmax=2^20-100000, ymin=100000, ymax=2^20-100000)
colnames(matrix.node.fSOM.scaled) <- c("FlowSOM.x","FlowSOM.y")

node.fSOM.x.scaled.FU <- matrix.node.fSOM.scaled[,1]
node.fSOM.y.scaled.FU<- matrix.node.fSOM.scaled[,2]

#####parameter node NBM
node.fSOM <- vector()
node.fSOM.x <- NULL
node.fSOM.y <- NULL
layout.fcs <- fSOM[["MST"]][["l"]]
mapping.fcs <- fSOM[["map"]][["mapping"]][(lFU+1):lNBM]
vertex.size.fcs <- fSOM[["map"]][["pctgs"]]*500

nb.cluster <- (fSOM[["map"]][["xdim"]]) * (fSOM[["map"]][["ydim"]])

nb.nodes <- nb.cluster
nb.points.per.line <- nb.cells /nb.nodes


for (i in 1:nb.cluster) {
  tempo <- mapping.fcs[mapping.fcs==i]
  tempo_indices <- which(mapping.fcs==i)
  
  random.nb <- runif(n=tempo,min = 0, max = 1 )
  
  random.radius <-  random.nb * vertex.size.fcs[i]/2
  random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
  
  node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
  node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
  
}


matrix.node.fSOM <- cbind(node.fSOM.x, node.fSOM.y)
matrix.node.fSOM.scaled <- igraph::norm_coords(matrix.node.fSOM, xmin=100000, xmax=2^20-100000, ymin=100000, ymax=2^20-100000)
colnames(matrix.node.fSOM.scaled) <- c("FlowSOM.x","FlowSOM.y")

node.fSOM.x.scaled.NBM <- matrix.node.fSOM.scaled[,1]
node.fSOM.y.scaled.NBM<- matrix.node.fSOM.scaled[,2]
par(mfrow=c(1,1))

node.fSOM.x.scaled<-c(node.fSOM.x.scaled.Diag,node.fSOM.x.scaled.FU,node.fSOM.x.scaled.NBM)
node.fSOM.y.scaled<-c(node.fSOM.y.scaled.Diag,node.fSOM.y.scaled.FU,node.fSOM.y.scaled.NBM)


#Blasts selection


ben.1 <- table.to.be.computed[DIAG_indices,]
ben.2 <- linear[DIAG_indices,]


blasts.nodes <- cell.node[DIAG_indices] 
blasts.nodes.sorted <- rle(sort(blasts.nodes))
blasts.nodes.of.interrest <- blasts.nodes.sorted$values


############### table count DIAG
fSOM1 <- FlowSOMSubset(fSOM, DIAG_indices)
fSOM2 <- FlowSOMSubset(fSOM, FU_indices)

table.count.DIAG <- rle(sort(fSOM1[["map"]][["mapping"]][,1]))

table.count.DIAG_indices <- rep(0,  nb.cluster)

table.count.DIAG.1 <- data.frame("nodes"=table.count.DIAG$values, 
                                 "DIAG.nb.cells"=table.count.DIAG$lengths,
                                 row.names = NULL)
cells.SUM.DIAG <- apply(table.count.DIAG.1,2,sum)

table.count.DIAG.1$DIAG.pct.cells <- table.count.DIAG.1 [,2]/cells.SUM.DIAG[2] *1000

table.count.DIAG.1 <- as.matrix(table.count.DIAG.1)

nodes <- 1: fSOM[["map"]][["nNodes"]]
nb.nodes <- fSOM[["map"]][["nNodes"]]
new.DF <- NULL
new.DF <- matrix(ncol=3, nrow=nb.nodes)
new.DF[,1] <- 1:nb.nodes

loop.2 <- nrow(table.count.DIAG.1)

for(k in 1:nb.nodes){
  tempo1 <- nodes[k]
  for(l in 1:loop.2){
    tempo2 <- table.count.DIAG.1[l,1]
    if (tempo2==tempo1) {
      tempo3 <- unname(table.count.DIAG.1[l,])
      new.DF[k,] <- tempo3
    }
  }
}

table.count.DIAG.2 <- new.DF 
colnames(table.count.DIAG.2) <- colnames(table.count.DIAG.1)

table.count.DIAG.2[is.na(table.count.DIAG.2)] <- 0

######################## table NBM
fSOM3 <- FlowSOMSubset(fSOM,NBM_indices)

table.count.NBM <- rle(sort(fSOM3[["map"]][["mapping"]][,1]))
table.count.NBM.1 <- data.frame("nodes"=table.count.NBM$values, 
                                "NBM.nb.cells"=table.count.NBM$lengths,
                                row.names = NULL)
cells.SUM.NBM <- apply(table.count.NBM.1, 2,sum)
table.count.NBM.1$NBM.pct.cells <- table.count.NBM.1 [,2]/cells.SUM.NBM[2] *100
table.count.NBM.1 <- as.matrix(table.count.NBM.1)

new.DF <- matrix(ncol=3, nrow=nb.nodes)
new.DF[,1] <- 1:nb.nodes

for(k in 1:nb.nodes){
  tempo1 <- nodes[k]
  for(l in 1:nrow(table.count.NBM.1)){
    tempo2 <- table.count.NBM.1[l,1]
    if (tempo2==tempo1) {
      tempo3 <- unname(table.count.NBM.1[l,])
      new.DF[k,] <- tempo3
    }
  }
}
table.count.NBM.2 <- new.DF 
colnames(table.count.NBM.2) <- colnames(table.count.NBM.1)

table.count.NBM.2[is.na(table.count.NBM.2)] <- 0


###table summary DIAG // NBM , creation d'une table de comptage
table.3 <- cbind("nodes"=table.count.DIAG.2[,1],
                 "DIAG % cells"=table.count.DIAG.2[,3], 
                 "NBM % cells" =table.count.NBM.2[,3],
                 "Ratio % DIAG NBM"= table.count.DIAG.2[,3]/table.count.NBM.2[,3],
                 "DIAG nb cells"=table.count.DIAG.2[,2],
                 "NBM nb cells"=table.count.NBM.2[,2])


table.3 <- table.3 [blasts.nodes.of.interrest,]

write.csv(table.3, file= "table.csv", row.names=FALSE)


#####extraction of node of interrest based on boolean equation
nodes.interest <- c()


nodes.interest <- (table.3 [,2] > DIAG.pct.threshold) & (table.3 [,4] >= Ratio.threshold) & (table.3[,5]>DIAG.cell.nb.threshold)

node.interest <- table.3[nodes.interest,1]


####Perfrom a PCA to be able to select NOI

PCA.nodes.ade4_result <- dudi.pca(table.3[,2:3], scannf = FALSE, nf= 2,
                                  scale=T, center = T)

(PCA.nodes.ade4_result$eig)

par.default; par(mar= c(5.1, 4.1, 4.1, 2.1));par(mfrow=c(1,1))
windows()
plot(PCA.nodes.ade4_result$li[,1], PCA.nodes.ade4_result$li[,2],
      pch=20, cex=1.5,
      col= table.3[,1], 
      main="Dudi.PCA center=T, scale=F on FlowSOM nodes results")
      text(PCA.nodes.ade4_result$li[,1], PCA.nodes.ade4_result$li[,2], 
      labels = table.3[,1], cex=0.7, pos = 3)
box()


matrix.PCA <- cbind(PCA.nodes.ade4_result$li[,1], PCA.nodes.ade4_result$li[,2])
matrix.PCA.scaled <- igraph::norm_coords(matrix.PCA, xmin=100000, xmax=2^20-100000, ymin=100000, ymax=2^20-100000)
colnames(matrix.PCA.scaled) <- c("PCA.1","PCA.2")

PCA1.scaled <- matrix.PCA.scaled[,1]
PCA2.scaled <- matrix.PCA.scaled[,2]

PCA.table <- cbind(table.3[,1],PCA1.scaled,PCA2.scaled)
colnames(PCA.table) <- c("nodes","PCA1.scaled", "PCA2.scaled")


##PCA integration in fcs file
PCA.table.1 <- matrix(ncol=4, nrow=length(nodes.mapping))
PCA.table.1[,1] <- 1:length(nodes.mapping)
PCA.table.1[,2] <- nodes.mapping
PCA.table.1[,3] <- 0 #PCA1
PCA.table.1[,4] <- 0 #PCA2
colnames(PCA.table.1) <- c("indices", 
                           "nodes in selected targeted region",
                           "PCA1 on nodes in selected region",
                           "PCA2 on nodes in selected region")

for(i in 1:nrow(PCA.table)){
  tempo.PCA <- PCA.table.1[,2]%in%PCA.table[i,1]
  tempo.PCA.indices <- which(tempo.PCA==TRUE)
  
  random.nb <- runif(n=length(tempo.PCA.indices), min = 0, max = 1)
  random.radius <- random.nb * 800000
  random.angle <- runif(n=length(tempo.PCA.indices),min = 0, max = 2*pi)
  
  PCA.table.1[tempo.PCA.indices,3] <- 
    ((random.radius*cos(random.angle))/dispersion.coef) + PCA.table[i,2] 
  PCA.table.1[tempo.PCA.indices,4] <- 
    ((random.radius*sin(random.angle))/dispersion.coef) + PCA.table[i,3]
}


###NOI integration in FCS file

DF.NOI <- data.frame(Nodes=rep("   ",length(node.interest)), 
                     your.choice=rep("   ",length(node.interest)),
                     stringsAsFactors = FALSE)

DF.NOI [1:length(node.interest),1] <- node.interest
DF.NOI [1:length(node.interest),2] <- rep(1,length(node.interest) )

DF.NOI<-edit(DF.NOI)
which(DF.NOI[,2]==1)




if(!is.null(length(node.interest))){
  
  BM.node.per.cell <- fSOM[["map"]][["mapping"]][,1]
  BM.node.per.cell.1 <- which(BM.node.per.cell%in%node.interest)
  BM.node.per.cell.2 <- sort( BM.node.per.cell.1)
}else{
  BM.node.per.cell.2<- 0
  nodes.interest <- 0
}
#remove variable from vector matrix list dataframe from RAM memory
#remove(ff, fSOM_res, fSOM, fSOM1, fSOM3)

col.of.interest <- 1:ncol(table.to.be.used.as.fcs.file)
col.names.fcs <- colnames(table.to.be.used.as.fcs.file)

fcs.Gated.CD45.comp.csv  <- table.to.be.used.as.fcs.file

###attribution d'une valeur binaire
nb.cells <- nrow(table.to.be.used.as.fcs.file)
nodes.of.interrest <- rep(100,nb.cells)
nodes.of.interrest[BM.node.per.cell.2] <- 200

#Free flowSOM integration
csv.free <- cbind(fcs.Gated.CD45.comp.csv, 
                  #"tag, 
                  node.fSOM.x.scaled, node.fSOM.y.scaled, 
                  nodes.of.interrest,fSOM[["map"]][["mapping"]][,1],data.metacluster,
                  PCA.table.1[,3], PCA.table.1[,4])




colnames(csv.free) <- c(colnames(table.to.be.used.as.fcs.file),
                        #"Tag", 
                        "FlowSOM Nodes.x", 
                        "FlowSOM Nodes.y",
                        "Nodes of Interrest","Nodes", "MetaClusters",
                        "PCA.1 on fSOM nodes in targeted region",
                        "PCA.2 on fSOM nodes in targeted region")


csv.merge <- cbind(csv.free, csv.FROZEN)

fcs.description.name.1 <- colnames(table.to.be.used.as.fcs.file) #name
fcs.description.name <- c(fcs.description.name.1, 
                          #"Tag cluster MST", 
                          "FlowSOM FREE.x", "FlowSOM FREE.y", "Nodes of Interrest","Nodes", "MetaClusters",
                          "PCA.1 on fSOM nodes in targeted region",
                          "PCA.2 on fSOM nodes in targeted region",
                          "FlowSOM FROZEN.x", "FlowSOM FROZEN.y")

fcs.description.desc.1 <- colnames(table.to.be.used.as.fcs.file)
fcs.description.desc <-  c(fcs.description.desc.1, 
                           #"Tag cluster MST",
                           "FlowSOM FREE.x", "FlowSOM FREE.y",  "Nodes of Interrest","Nodes","MetaClusters",
                           "PCA.1 on fSOM nodes in targeted region",
                           "PCA.2 on fSOM nodes in targeted region",
                           "FlowSOM FROZEN.x", "FlowSOM FROZEN.y")


#Define the complete output fcs file name


sFilename <- paste(output.folder.final, sample.name, ".fcs", sep="")
print(sFilename)


colnames(csv.merge) <- c(fcs.description.name.1,
                         "FlowSOM FREE.x", "FlowSOM FREE.y", "Nodes of Interrest","Nodes", "MetaClusters",
                         "PCA.1 on fSOM nodes in targeted region",
                         "PCA.2 on fSOM nodes in targeted region",
                         "FlowSOM FROZEN.x", "FlowSOM FROZEN.y") 

csv.File2be.Converted <- csv.merge

#calculate number of events and number of parameters
cParams <- ncol(csv.File2be.Converted)
nb.cells <- nrow(csv.File2be.Converted)
cEvents <- nrow(csv.File2be.Converted)


ff <- flowFrame(csv.File2be.Converted)

#write.FCS(ff, sFilename)

fcs.description.name <- colnames(csv.File2be.Converted)
fcs.description.desc <- colnames(csv.File2be.Converted)

dta <-csv.File2be.Converted

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

write.FCS(ff,sFilename)
