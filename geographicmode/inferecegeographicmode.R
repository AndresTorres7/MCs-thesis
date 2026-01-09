
library(gstat)
library(abc)
library(caret)
library(DDD)
library(gridExtra)
library(ggplot2)
library(MASS)
library(phytools)
library(RPANDA)
library(RRphylo)
library(stringi)
library(stringr)
library(tidyverse)
library(rgbif)
library(readxl)
library(terra)
library(RANN)

scripts <- c("rangeDispersal.R", "nicheEvolution.R","speciateAllopatric.R","speciateSympatric.R",
             "speciateParapatric.R","speciateDispersal.R","seedSpecies.R","environmentalChange.R",
             "nicheRecenter.R","DREaD.R","generateSummaryStatistics.R", "helperFunctions.R","findSisters.R","summary_statsitics_functions.R")



for (i in 1:length(scripts)){
  scripts[[i]] <- paste("~/Biologia/SkeelsShareddiversificationTest/Code/AuxiliaryScripts/DREAD/",scripts[[i]],sep="")
}
lapply(scripts, source)

rangeDispersal <- disperseRange
required.packages <- (c("raster","gstat", "SpaDES", "ape","phytools","geiger",
                        "phyloclim","ggplot2","gridExtra","moments",
                        "apTreeshape","parallel", "doSNOW", "rgeos","knitr",
                        "data.table", "fossil", "ENMTools"))

lapply(required.packages, require, character.only=T)

setwd("~/Biologia/Maestria/Tesis/official analysis/MCs-thesis/geographicmode")
load("MicraRasters.rda")

sp.rasters <- k
spat.rasters <- k
names(sp.rasters) <- gsub("Micrathena", "M", names(sp.rasters))
names(spat.rasters) <- gsub("Micrathena", "M", names(spat.rasters))


for( f in 1:length(sp.rasters)) {
  
  origin(sp.rasters[[f]])<-origin(sp.rasters[[1]])
  
  
}

clade.trees <- readNexus("k.BEAST-100-posterior-trees.nex")
newclade.trees <- list()

for (i in seq_along(clade.trees)) {
  tr <- clade.trees[[i]]
  tr$tip.label <-   sub("^(([^_]+_[^_]+)).*", "\\1", tr$tip.label)
  remove<-tr$tip.label[!tr$tip.label %in% names(sp.rasters)] 
  tr<-drop.tip(tr, tip = remove)
  newclade.trees[[i]] <- tr

}

names(newclade.trees) <- names(clade.trees)


df <- data.frame(clade = NA, ntips = NA, crown.age=NA,
                 RO0=NA, RO50=NA, RO75=NA, RO90=NA, RO100=NA, ROmean=NA, ROslope=NA, ROintercept=NA, ROskew=NA, ROkurtosis=NA,
                 ASYMmean=NA, ASYMslope=NA, ASYMintercept=NA,
                 BIMOD50=NA, BIMOD75=NA, BIMOD90=NA, BIMOD100=NA,
                 RSskew=NA, RSmean=NA, RSsd=NA, CI=NA, Beta=NA, Gamma=NA, SI=NA,
                 lambda.x=NA, lambda.y=NA, TD=NA )

results_df <- data.frame()

for (w in 1:length(newclade.trees)){
  phy = newclade.trees[[w]]
  Clade= names(newclade.trees)[w]
  
  root.age.true <- max(findMRCA(phy, type=c("height")))
  phy$edge.length<-phy$edge.length/root.age.true
  root.age<-1
  
  # get sister species
  sisters <- findSisters(phy, solve.polytomies=T)
  
  # find divergence times
  sisters.split.time <- sisterSpeciesSplitTimes(phy, sisters, root.age=1)
  
  
  #get sister species range overlap
  RO_ss <- sisterSpeciesOverlap(sp.rasters, sisters)
  
  #get sister species range asymmetry
  Asym_ss <- sisterSpeciesRangeAsymmetry(sp.rasters, sisters)
  
  # get sister species range distances
  print(paste("running range distance ... takes a min"))
  #trace("sisterSpeciesRangeDistances", edit = T)  # 3 times min(pointDistance(p1, p2, longlat=F) to min(nn2(data = p2, query = p1, k = 1)$nn.dists)
  RD_ss <- sisterSpeciesRangeDistances(sp.rasters, sisters, aggregate_size = 2)
  
  # get sister-species-outgroup overlap metrics
  #TO <- outgroupOverlap(sp.rasters, sisters, sisters.overlap = RO_ss)
  
  #! first batch of summary statistics - mean and SD
  range.overlap.av <- mean(RO_ss, na.rm=T)
  range.overlap.sd <- sd(RO_ss, na.rm=T)
  
  range.asymmetry.av <- mean(Asym_ss, na.rm=T)
  range.asymmetry.sd <- sd(Asym_ss, na.rm=T)
  
  range.distance.av <- mean(RD_ss, na.rm=T)
  range.distance.sd <- sd(RD_ss, na.rm=T)
  
  #TO.av <- mean(TO, na.rm=T)
  #TO.sd <- sd(TO, na.rm=T)
  
  #! second batch of summary statistics - linear models
  
  lm.asym    <-     lm(Asym_ss ~ sisters.split.time)
  lm.overlap <-  lm(RO_ss ~ sisters.split.time)
  lm.distance <- lm(RD_ss ~ sisters.split.time)
  
  range.overlap.slope <- lm.overlap$coefficients[[2]]
  range.overlap.intercept <- lm.overlap$coefficients[[1]]
  
  range.distance.slope <- lm.distance$coefficients[[2]]
  range.distance.intercept <- lm.distance$coefficients[[1]]
  
  range.asymmetry.slope<-lm.asym$coefficients[[2]]
  range.asymmetry.intercept<-lm.asym$coefficients[[1]]
  
  #! third batch of summary statistics - species range sizes
  
  species.range.sizes <- rangeSize(sp.rasters)
  species.range.sizes.stand <- species.range.sizes/max(species.range.sizes)
  
  range.size.skew <- skewness(species.range.sizes.stand)
  range.size.av <- mean(species.range.sizes.stand, na.rm=T)
  range.size.sd <- sd(species.range.sizes.stand, na.rm=T)
  
  #! fourth batch of summary statistics - Phylogenetic balance and treeshape metrics
  
  phy_treeshape <- as.treeshape(phy)
  beta <- maxlik.betasplit(phy_treeshape, up = 10, remove.outgroup = FALSE, confidence.interval = "none", conf.level = 0.95, size.bootstrap = 100)$max_lik
  collessI <- colless(phy_treeshape)
  sackinI <- sackin(phy_treeshape)
  ltt <- ltt(phy, plot=FALSE)
  gamma <-ltt$gamma
  
  #! fifth batch of summary statistics - further overlap summary metrics
  
  symp0 <-   length(which(RO_ss == 0))/   length(RO_ss)
  symp50 <-  length(which(RO_ss >= 0.5))/ length(RO_ss)
  symp75 <-  length(which(RO_ss >= 0.75))/length(RO_ss)
  symp90 <-  length(which(RO_ss >= 0.9))/ length(RO_ss)
  symp100 <- length(which(RO_ss == 1))/   length(RO_ss)
  
  range.overlap.kurt<- kurtosis(RO_ss)
  range.overlap.skew <- skewness(RO_ss)
  
  #! sixth batch of summary statistics - biomodality
  sisterpairs <- length(sisters)
  
  bimodality100 <-((symp0*sisterpairs) *  (symp100*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  bimodality90  <-((symp0*sisterpairs) *  (symp90*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  bimodality75  <-((symp0*sisterpairs) *  (symp75*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  bimodality50  <-((symp0*sisterpairs) *  (symp50*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  
  #! seventh batch of summary statistics - ARC
  print(paste("running ARC ... takes a min"))
  
 #trace("geog.range.overlap", edit = T)# Change SpatRaster to RasterLayer to fix the version problem of raster and Terra
  
  ARC <- calculateARC(phy, sp.rasters)
  
  
  
  ARCint <- ARC$coefficients[[1,1]]#What's this?
  ARCslope <- ARC$coefficients[[1,2]]#What's this?
  
  
  
  # put it together
  
  tmp_df <- data.frame(clade         = Clade,
                       ntips         = length(phy$tip.label), 
                       ROmean        = range.overlap.av,
                       ROsd          = range.overlap.sd,
                       ROslope       = range.overlap.slope, 
                       ROintercept   = range.overlap.intercept, 
                       ROskew        = range.overlap.skew, 
                       ROkurtosis    = range.overlap.kurt,
                       RO0           = symp0, 
                       RO50          = symp50, 
                       RO75          = symp75, 
                       RO90          = symp90, 
                       RO100         = symp100,
                       ASYMmean      = range.asymmetry.av, 
                       ASYMsd        = range.asymmetry.sd,
                       ASYMslope     = range.asymmetry.slope, 
                       ASYMintercept = range.asymmetry.intercept,
                       RDmean        = range.distance.av, 
                       RDsd          = range.distance.sd, 
                       RDintercept   = range.distance.intercept, 
                       RDslope       = range.distance.slope,
                       BIMOD50       = bimodality50, 
                       BIMODE75      = bimodality75, 
                       BIMOD90       = bimodality90, 
                       BIMOD100      = bimodality100,
                       RSskew        = range.size.skew, 
                       RSmean        = range.size.av, 
                       RSsd          = range.size.sd, 
                       CI            = collessI, 
                       Beta          = beta, 
                       Gamma         = gamma, 
                       SI            = sackinI,
                       ARCslope      = ARCslope, 
                       ARCint        = ARCint) 
  
  results_df <- rbind(results_df, tmp_df)
  
 
}

setwd("~/Biologia/Maestria/Tesis/official analysis/MCs-thesis/geographicmode")
write.csv(results_df,file="Empirical_data_Micra.csv")



### Results and figures

setwd("~/Biologia/Maestria/Tesis/official analysis/MCs-thesis/geographicmode")

dat <- read.csv("simulation_results.csv",stringsAsFactors = F)
variabs <- c("asymintercept","ROintercept","RSmean",#These are the 14 statistics used in model selection
             "ROmean", "SI", "RDintercept", "ARCslope",
             "RO90", "RSsd","RO100", "RDsd", "ARCint" ,"asymslope" )


dat <- dat[!grepl("EXTINCT", dat[,14]),]#Delete simulations that resulted in extinction
for (i in 1:ncol(dat)){
  dat[,i] <- gsub("NaN",NA,dat[,i])
}
dat <- na.omit(dat)
dat[,14:45] <- as.numeric(as.character(unlist(dat[,14:45])))
df <- data.frame(dat$geo.mode,dat[,which(colnames(dat) %in% variabs)])
colnames(df)[[1]] <- "geo.mode"
models <- as.character(df[,1])
df <- df[,2:length(df)]

clados <- names(clade.trees)

for (l in 1:length(cladoss)) {
  


#Model selection to infer the predominant geographic mode of speciation: ABC

clades <- names(clados[[l]])
all.res.list <- list()

setwd("~/Biologia/Maestria/Tesis/official analysis/MCs-thesis/geographicmode")

for (i in 1:length(clades)){
  clade <- clades[i]
  
  emp.dat <- read.csv(paste("./Empirical_data_",Clade,".csv",sep=""),stringsAsFactors = F)
  
  colnames(emp.dat)[[which(colnames(emp.dat)=="ASYMslope")]] <- "asymslope"
  colnames(emp.dat)[[which(colnames(emp.dat)=="ASYMintercept")]] <- "asymintercept"
  
  emp.dat <- emp.dat[,which(colnames(emp.dat) %in% variabs)]
  
  emp.dat <- emp.dat[,(colnames(df))]
  
  ###
  
  emp.log <- postpr(emp.dat, models, df, tol=.05, method="mnlogistic")
  emp.neu <- postpr(emp.dat, models, df, tol=.05, method="neuralnet")
  
  table.tmp <- cbind(emp.log$pred,emp.neu$pred)
  all.res.list[[i]] <- table.tmp
  names(all.res.list)[i] <- names(clades)[i]
}

all.res.table <- cbind(all.res.list[[1]])#,all.res.list[[2]],all.res.list[[3]])
colnames(all.res.table) <- rep(clades,each=2)
colnames(all.res.table) <- paste(colnames(all.res.table),rep(c("log","neu"),1),sep="_")
all.res.table <- cbind(rownames(all.res.table),as.data.frame(all.res.table))
colnames(all.res.table)[1] <- "Mode"
all.res.table[,1] <- as.character(all.res.table[,1])


setwd("~/Biologia/Maestria/Tesis/official analysis/MCs-thesis/geographicmode/resultspersp")
dir.create(names(clades))
setwd(paste("~/Biologia/Maestria/Tesis/official analysis/MCs-thesis/geographicmode/resultspersp",
            names(clades), sep = "/"))
write.csv(all.res.table, "MicraABCresults.csv",row.names = F)

all.res.table 



#Model selection to infer the predominant geographic mode of speciation: LDA

df.lda <- cbind(df,models)
lda.res <- lda(models ~ asymintercept + ROintercept + RSmean +
                 ROmean + SI + RDintercept + ARCslope + ARCint +
                 RO90 + RSsd + RO100 + RDsd + asymslope, df.lda)
emp.table <- data.frame(matrix(ncol=ncol(emp.dat),nrow=length(clades))) # ARCslope + ARCint +
for (i in 1:length(clades)){
  
  
  clade <- clades[i]
  setwd("~/Biologia/Maestria/Tesis/official analysis/MCs-thesis/geographicmode")
  emp.dat <- read.csv(paste("./Empirical_data_",clade,".csv",sep=""),stringsAsFactors = F)
  colnames(emp.dat)[[which(colnames(emp.dat)=="ASYMslope")]] <- "asymslope"
  colnames(emp.dat)[[which(colnames(emp.dat)=="ASYMintercept")]] <- "asymintercept"
  
  emp.dat <- emp.dat[,which(colnames(emp.dat) %in% variabs)]
  
  emp.dat <- emp.dat[,(colnames(df))]
  emp.table[i,] <- emp.dat
}
colnames(emp.table) <- colnames(emp.dat)
colnames(emp.table)[[which(colnames(emp.table)=="ASYMslope")]] <- "asymslope" #
colnames(emp.table)[[which(colnames(emp.table)=="ASYMintercept")]] <- "asymintercept" #
emp.table <- emp.table[,which(colnames(emp.table) %in% variabs)]
emp.table <- emp.table[,(colnames(df))]
emp.table$models <- rep("allopatric",nrow(emp.table))
predictions <- lda.res %>% predict(emp.table)
lda.plot <- cbind(df.lda, predict(lda.res)$x)

setwd(paste("~/Biologia/Maestria/Tesis/official analysis/MCs-thesis/geographicmode/resultspersp",
            names(clades), sep = "/"))
write.table(predictions,  "MicraLDA_results.txt")
write.csv(lda.plot, "backgroundpoints.csv", row.names = F)

a<-data.table(Reg = 1:15, 
              value = c(all.res.table[,2], all.res.table[,3], c(predictions$posterior)), ########
              Mode = c(all.res.table$Mode, all.res.table$Mode, all.res.table$Mode), 
              model = c(rep( "mnL", times = 5), rep( "NN", times = 5), rep( "LDA", times = 5) )
)

all.res.table
a
write.csv(a, paste(Clade,"Allvalues.csv",sep = "_"), row.names = F)


tryCatch({
  
  pdf(paste(clade, "ABC.pdf"))
  
  plo <- ggplot(a, aes(fill=Mode, y=value, x=model)) + 
    geom_bar(position="stack", stat="identity") +
    ggtitle("Micrathena")
  print(plo)
  dev.off()
  
}, error = function(e) {
  message("Error creating ABC.pdf for ", clade, ": ", e$message)
}, finally = {
  if (dev.cur() != 1) dev.off()
})

Sys.sleep(0.2)



li<-predictions$x[,1:2]
u<-as.data.frame(li)
p<-as.data.frame(t(u))

tryCatch({
  
  pdf( "MicraLDAgraphic.pdf")
  plo <- ggplot(lda.plot, aes(LD1, LD2)) +
    geom_point(aes(color = models)) +
    stat_ellipse(aes(color = models),level=0.95) +   
    theme_bw() +
    geom_point(data= p, color = "black") +
    ggtitle("Micrathena")
  print(plo)
  dev.off()
  
}, error = function(e) {
  message("Error creating LDAgraphic.pdf for ", clade, ": ", e$message)
}, finally = {
  if (dev.cur() != 1) dev.off()
})

} 