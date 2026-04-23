#### Morphological evolution

#Functions 

mayorentredoscol <- function(x,y) {
  finalvector <- numeric()
  
  for (i in 1:length(x)) {
    value <- max(x[i],y[i])
    finalvector<-c(finalvector,value)
    
  }
  return(finalvector)
}


  
  
  


library(phytools)
library(diversitree)
library(geiger)
library(dplyr)
library(RPANDA)



## Reading tree

micratree <- readNexus("k.BEAST-100-posterior-trees.nex")
newclade.trees <- list()

for (i in seq_along(micratree)) {
  tr <- micratree[[i]]
  
  tr<-drop.tip(tr,"M_horrida_122") # M horrida is repeated
  tr$tip.label <- sub("_", ".", tr$tip.label)
  tr$tip.label <- sub("\\_.*", "", tr$tip.label)
  
  newclade.trees[[i]] <- tr
  
}

names(newclade.trees) <- names(micratree)

micratree <- newclade.trees


## Reading morphological data

morphdata <- read.csv("datosespinas.csv", sep= ";")
morphdata$Species <- paste("M",morphdata$Species,sep = ".") 

## Morph data preparation; average of spine's measurements, ratio of spine length 
## and carapace length and average of the ratio per species

spinedata <- data.frame(species = morphdata[,1], carapace = morphdata[,11])

spinedata$firstspine <- rowMeans(morphdata[,c(6,7)])
spinedata$secondspine <- rowMeans(morphdata[,c(8,9)])
spinedata$longestspine <- mayorentredoscol(spinedata$firstspine,spinedata$secondspine) 
spinedata$spinelength <- spinedata$longestspine / spinedata$carapace
spinedata$fspinelength <- spinedata$firstspine / spinedata$carapace
spinedata$sspinelength <- spinedata$secondspine / spinedata$carapace

#the longest spine
spinedataready<-aggregate(spinelength ~ species, data = spinedata, FUN = mean)
row.names(spinedataready)<- spinedataready$species
spinedataready <- setNames(spinedataready$spinelength, rownames(spinedataready))
# 
#the first spine
fspinedataready<-aggregate(fspinelength ~ species, data = spinedata, FUN = mean)
row.names(fspinedataready)<- fspinedataready$species
fspinedataready <- setNames(fspinedataready$fspinelength, rownames(fspinedataready))
# 
# Check names in both data sets

# nameschecked<-name.check(phy = micratree[[1]],data = fspinedataready)
# for(g in seq_along(newclade.trees)) {
#   newclade.trees[[g]] <- drop.tip(newclade.trees[[g]],nameschecked$tree_not_data)
# }
# micratree <- newclade.trees
# fspinedataready <- fspinedataready[!names(fspinedataready) == nameschecked$data_not_tree]
# name.check(phy = micratree[[1]],data = fspinedataready)
# # 
# #the second spine
# sspinedataready<-aggregate(sspinelength ~ species, data = spinedata, FUN = mean)
# row.names(sspinedataready)<- sspinedataready$species
# sspinedataready <- setNames(sspinedataready$sspinelength, rownames(sspinedataready))

# Check names in both data sets

# nameschecked<-name.check(phy = micratree,data = sspinedataready)
# treeready <- drop.tip(micratree,nameschecked$tree_not_data)
# sspinedataready <- sspinedataready[!names(sspinedataready) == nameschecked$data_not_tree]
# name.check(phy = treeready,data = sspinedataready)
# 
# spinas <- as.data.frame(spinedataready)
#write.csv(spinas, "spinedatacleaned.csv", row.names = T)

# Check names in both data sets

nameschecked<-name.check(phy = micratree[[1]],data = spinedataready)

for(g in seq_along(newclade.trees)) {
  newclade.trees[[g]] <- drop.tip(newclade.trees[[g]],nameschecked$tree_not_data)
}
micratree <- newclade.trees
spinedataready <- spinedataready[!names(spinedataready) == nameschecked$data_not_tree]
name.check(phy = micratree[[66]],data = spinedataready)



# data <- read.csv("spinedatacleaned.csv", row.names = 1)
# data$fspine <- fspinedataready
# data$sspine <- sspinedataready
# head(data,3)


options(warn = 1)


#################################################### Here start results ##########----

################   max spine length


aics <- data.frame(B = numeric() ,EB= numeric(),OU=numeric())
modelspertree <- list()

for(i in 1:length(micratree)) {
  
  cat("Fitting tree number",i,"\n")
  
  ## fit Brownian motion model using fitContinuous
  fitBM_gs<-fitContinuous(micratree[[i]],spinedataready, control=list(niter=50), ncores=2)
  cat("B okay")
  ## fit EB model to genome size
  fitEB_gs<-fitContinuous(micratree[[i]],spinedataready,
                          model="EB", bounds = list(a = c(-0.5,10^-20)), control=list(niter=50), ncores=2)
  cat("EB okay")
  ## fit OU model to genome size
  
  fitOU_gs<-fitContinuous(micratree[[i]],spinedataready,
                          model="OU",bounds=list(alpha=c(0,10)), control=list(niter=100), ncores=2)
  cat("OU okay")
  
  
  models<-list(B = fitBM_gs, EB = fitEB_gs, OU = fitOU_gs)
  #saveRDS(models, paste("models",i,".rds", sep = "")) 
  modelspertree[[i]] <- models
  
  ## accumulate AIC scores from our three models into
  ## a vector
  
  aics[i,1] <- AIC(fitBM_gs) 
  aics[i,2] <- AIC(fitEB_gs)
  aics[i,3] <- AIC(fitOU_gs)
  
}


saveRDS(modelspertree, "newmodelscorrected.rds")



####################### first spine

modelspertree <- list()

for(i in 1:length(micratree)) {
  
  cat("Fitting tree number",i,"\n")
  
  ## fit Brownian motion model using fitContinuous
  fitBM_gs<-fitContinuous(micratree[[i]],fspinedataready , control=list(niter=50), ncores=2)
  cat("B okay")
  ## fit EB model to spine length
  fitEB_gs<-fitContinuous(micratree[[i]],fspinedataready ,
                          model="EB", bounds = list(a = c(-0.5,10^-20)), control=list(niter=50), ncores=2)
  cat("EB okay")
  ## fit OU model to genome size
  
  fitOU_gs<-fitContinuous(micratree[[i]],fspinedataready ,
                          model="OU",bounds=list(alpha=c(0,10)), control=list(niter=50), ncores=2)
  cat("OU okay")
  
  
  models<-list(B = fitBM_gs, EB = fitEB_gs, OU = fitOU_gs)
  #saveRDS(models, paste("models",i,".rds", sep = "")) 
  modelspertree[[i]] <- models
  
  
}



phylosig(micratree[[i]], spinedataready, method = "lambda", test = T)

micratree_scaled <- lapply(micratree, function(tr) {
  tr$edge.length <- tr$edge.length / max(node.depth.edgelength(tr))
  tr
})

fitOU_gs<-fitContinuous(micratree_scaled[[i]],spinedataready,
                        model="OU",bounds=list(alpha=c(0,100)), control=list(niter=100), ncores=2)






modelos <- modelspertree

modelos <- readRDS("modelsmorphevolution.rds")

aicsmodelos <- data.frame(Tree = numeric(), BMAIC = numeric(), EBAIC = numeric(), OUAIC = numeric())

for(i in seq_along(modelos)) {
  
  aicsmodelos[i,1] <- as.numeric(i)
  aicsmodelos[i,2] <- modelos[[i]]$B$opt$aic
  aicsmodelos[i,3] <- modelos[[i]]$EB$opt$aic
  aicsmodelos[i,4] <- modelos[[i]]$OU$opt$aic
}


write.csv(aicsmodelos, "aicvalues.csv", row.names = F)

bestmodel <- character()

for(i in 1:nrow(aicsmodelos)) {
  aic_values <- aicsmodelos[,-1]
  
  aic <- min(aic_values[i,])
   aicp <- which(aic == aic_values[i,])
   
   if(aicp == 3) {
     
     bestmodel <- c(bestmodel, "OU")
   } else if (aicp == 2) {
     
     bestmodel <- c(bestmodel, "EB")
   } else if (aicp == 1) {
     
     bestmodel <- c(bestmodel, "BM")
   } else {
     
     bestmodel <- c(bestmodel, "problem")
   }
  
}


alphas <- numeric()
sigs <- numeric()

for(i in which(bestmodel == "OU")) {
  
alphas <-  c(alphas, modelspertree[[i]]$OU$opt$alpha)
sigs <- c(sigs,modelspertree[[i]]$OU$opt$sigsq)
  
  
}

mean(alphas)
sd(alphas)
min(alphas)
max(alphas)

mean(sigs)
sd(sigs)
min(sigs)
max(sigs)


