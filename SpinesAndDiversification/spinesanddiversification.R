####  EVOLUTION OF SPINES #############-----------------------------------------


#Functions 

mayorentredoscol <- function(x,y) {
  finalvector <- numeric()
  
  for (i in 1:length(x)) {
    value <- max(x[i],y[i])
    finalvector<-c(finalvector,value)
    
  }
  return(finalvector)
}

sampFit <- function(x, frac, spinedata) {
  
  modelspertree <- list()
  
  for(i in 1:length(x)) {
    
    cat("Fitting models for tree number ",i,"\n")
    tot_time<-max(node.age(x[[i]])$ages)
    
    ## run fit.bd and fitContinuous to get starting
    ## values for our QuaSSE optimization
    
    bd<-fit.bd(x[[i]],rho=frac)
    bm<-fitContinuous(x[[i]],spinedata)
    p<-setNames(c(bd$b,bd$d,bm$opt$sigsq),
                c("lambda","mu","diffusion"))
    
    ### linear function to fit the QuaSSE model
    ## define range of x
    
    xr<-range(spinedata)+c(-1,1)*2*
      p["diffusion"]
    
    ## make linear model for QuaSSE
    
    linear.x<-make.linear.x(xr[1],xr[2])
    
    ## make QuaSSE likelihood function forvariable lambda and constrain
    
    lik.lambda<-make.quasse(x[[i]],spinedata,lambda=linear.x,mu=constant.x,
                            sampling.f=frac,states.sd=0.1)
    lik.lambda<-constrain(lik.lambda,drift~0) # We constrained one parameter, drift~0.The drift parameter inmake.quasse
                                               #allows the mean of our character evolution process to change throughtime.We’renotgoingto
                                                 # worryaboutthathere.
    ## subsample starting parameter values to match themodel we’re fitting
    
    pp<-setNames(c(p["lambda"],0,p["mu"],
                   p["diffusion"]),argnames(lik.lambda))

    
    ## fitour first QuaSSE model Mu constant and labda changes
    cat("Mu constant and labda changes")
    lambda.mle<-find.mle(lik.lambda,x.init=pp,
                         control=list(parscale=0.1),lower=rep(0,4))
    
    print(lambda.mle)
    
    ## fixed lambda model
    ## make QuaSSE likelihood function for variable mu and constrain
    
    lik.mu<-make.quasse(x[[i]],spinedata,lambda=constant.x,mu=linear.x,
                        sampling.f=frac,states.sd=0.1)
    lik.mu<-constrain(lik.mu,drift~0)
    
    ## fit variable mu model
    pp<-setNames(c(p[c("lambda","mu")],0,
                   p["diffusion"]),argnames(lik.mu))
    
    cat("Mu changes and labda constant")
    
    mu.mle<-find.mle(lik.mu,x.init=pp,
                     control=list(parscale=0.1),lower=rep(0,4))
    print(mu.mle)
    
    ## create full likelihood function
    cat("full model")
    lik.full<-make.quasse(x[[i]],spinedata,lambda=linear.x,mu=linear.x,
                          sampling.f=frac,states.sd=0.1)
    lik.full<-constrain(lik.full,drift~0)
    
    ## fitfull QuaSSE model
    pp<-setNames(c(lambda.mle$par[1:2],mu.mle$par[2:3],
                   p["diffusion"]),argnames(lik.full))
    
    full.mle<-find.mle(lik.full,x.init=pp,
                       control=list(parscale=0.1),lower=rep(0,5))
  print(full.mle)
    
  ## likelihood function for character
  ## independent model
  lik.cid<-make.quasse(x[[i]],spinedata,lambda=constant.x,mu=constant.x,
                       sampling.f=frac,states.sd=0.1)
  lik.cid<-constrain(lik.cid,drift~0)

  
  ## fitCID QuaSSE model and print coefficients
  cat("null model")
  cid.mle<-find.mle(lik.cid,x.init=p,
                    control=list(parscale=0.1),lower=rep(0,3))
  print(cid.mle)
   
    models<-list(null = cid.mle, b = lambda.mle, d = mu.mle, full = full.mle)
    saveRDS(models, paste("models",i,".rds", sep = "")) 
    modelspertree[[i]] <- models
    gc()
  }
  
  return(modelspertree)
}

## Required libraries

library(phytools)
library(diversitree)
library(geiger)
library(phytools)
library(dplyr)
library(RPANDA)
library(writexl)
library(readxl)
library(ggplot2)

## Reading morphological data

morphdata <- read.csv("datosespinas.csv", sep= ";")
morphdata$Species <- paste("M",morphdata$Species,sep = ".") 
## Morph data preparation; average of spine's measurements, ratio of spine length 
## and carapace length and average of the ratio per species

spinedata <- data.frame(species = morphdata[,1], carapace = morphdata[,11])

spinedata$firstspine <- rowMeans(morphdata[,c(6,7)])
spinedata$secondspine <- rowMeans(morphdata[,c(8,9)])
spinedata$longestspine <- mayorentredoscol(spinedata$firstspine,spinedata$secondspine) 
spinedata$spinelength <- spinedata$longestspine / spinedata$carapace # longest spine
spinedata$fspinelength <- spinedata$firstspine / spinedata$carapace  # first spine
spinedata$sspinelength <- spinedata$secondspine / spinedata$carapace  # second spine

# longest spine
spinedataready<-aggregate(spinelength ~ species, data = spinedata, FUN = mean)
row.names(spinedataready)<- spinedataready$species
spinedataready <- setNames(spinedataready$spinelength, rownames(spinedataready))
# first spine
fspinedataready<-aggregate(fspinelength ~ species, data = spinedata, FUN = mean)
row.names(fspinedataready)<- fspinedataready$species
fspinedataready <- setNames(fspinedataready$fspinelength, rownames(fspinedataready))
# second spine
sspinedataready<-aggregate(sspinelength ~ species, data = spinedata, FUN = mean)
row.names(sspinedataready)<- sspinedataready$species
sspinedataready <- setNames(sspinedataready$sspinelength, rownames(sspinedataready))


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


# Check names in both data sets

nameschecked<-name.check(phy = newclade.trees[[1]] ,data = spinedataready)

for(g in seq_along(newclade.trees)) {
  newclade.trees[[g]] <- drop.tip(newclade.trees[[g]],nameschecked$tree_not_data)
  
}

sspinedataready <- sspinedataready[!names(sspinedataready) == nameschecked$data_not_tree]
fspinedataready <- fspinedataready[!names(fspinedataready) == nameschecked$data_not_tree]
spinedataready <- spinedataready[!names(spinedataready) == nameschecked$data_not_tree]
name.check(phy = newclade.trees[[1]],data = spinedataready)

## Visualize the trait in the tree

# spine.map <-contMap(treeready, spinedataready, plot = F)
# spine.map <- setMap(spine.map, c("yellow", "blue"))
# plot(spine.map, lwd=c(2,5), outline = FALSE, ftype="off",
#      leg.txt="LongestSpine/carapace ratio", legend = 60)


## sampling fraction p for fitting QuaSSE model

samplingFraction <- (100/150) #fraction of extant species included in the phylogeny Pedro says are 150

spineresults <- sampFit(newclade.trees, samplingFraction,spinedataready)




saveRDS(spineresults, "spineresultscorrected.rds")

    spineresults <- readRDS("spineresultscorrected.rds") # DATA


# for(h in 1:47) {
#   
#   recover <- readRDS(paste("models", h, ".rds", sep = ""))
#   spineresults[[h]] <- recover
# }

    
    
    
results<- data.frame()
### anova among models

for (i in seq_along(spineresults)) {

  anovas <- anova(spineresults[[i]]$null,variable.lambda=spineresults[[i]]$b,
      variable.mu=spineresults[[i]]$d,full.model=spineresults[[i]]$full)

  results <- rbind(results,anovas)

}
x <- rep(1:100, each = 4)
results$treenum <- x
write.csv(results, "modelsinfo.csv")
results<- read.csv("modelsinfo.csv", row.names = 1)
colnames(results)[5] <- "Pvalue"



# calculating AIC among models

resultsaic<- data.frame()
for (i in seq_along(spineresults)) {
  
  aics <- AIC(spineresults[[i]]$null,spineresults[[i]]$b,spineresults[[i]]$d,spineresults[[i]]$full)
  
  resultsaic <- rbind(resultsaic,aics)
  
}
resultsaic$treenum <- x
rownames(resultsaic) <- row.names(results)
resultsaic$modelname <- rep(c("minimal", "variable.lambda","variable.mu","full.model"), 100) 



  
#### Selectiong the best model in each treee -----------------------------



library(dplyr)

best_models <- resultsaic %>%
   group_by(treenum) %>%
  mutate(
    minAIC = min(AIC),
    has_close_model = sum((AIC - minAIC) < 2 & (AIC - minAIC) > 0) > 0
  ) %>%
  slice_min(AIC, n = 1, with_ties = FALSE) %>%
  ungroup()

write.csv(best_models, "bestmodelsaic.csv", row.names = F )
write_xlsx(best_models, "bestmodelsaic.xlsx")

best_models <- read.csv("bestmodelsaic.csv") #DATA
bestmodelsaiccorrected <- read_xlsx("newbestmodelsaic.xlsx") #DATA


# frequency figure
setwd("~/Biologia/Maestria/Tesis/official analysis/MCs-thesis/spines/figures")





#p1 <- 
  ggplot(best_models, aes(x = reorder(modelname, modelname,
                               function(x)-length(x)))) +
  geom_bar(width = 0.7) +
  xlab("Mejor modelo de diversificación") +
  ylab("Frecuencia absoluta") +
   scale_x_discrete(labels = c("minimal" = "Nulo",
                               "variable.lambda" = "lambda variable",
                               "variable.mu" = "mu variable",
                               "full.model" = "Completo")) +
  theme_classic() +
  scale_fill_grey(start = 0.8, end = 0.3) +
  scale_y_continuous(
    limits = c(0, 80),
    breaks = seq(0, 120, by = 10)
  ) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 10)
  )
ggsave("modelfrequency.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("modelfrequency.pdf", plot = p1, width = 8, height = 6)


#p2 <- 
  ggplot(bestmodelsaiccorrected, aes(x = reorder(modelname, modelname,
                                    function(x)-length(x)))) +
  geom_bar(width = 0.7) +
  xlab("Mejor modelo de diversificación") +
  ylab("Frecuencia absoluta") +
  scale_x_discrete(labels = c("minimal" = "Nulo",
                              "variable.lambda" = "lambda variable",
                              "variable.mu" = "mu variable",
                              "full.model" = "Completo")) +
  theme_classic() +
  scale_fill_grey(start = 0.8, end = 0.3)+
  scale_y_continuous(
    limits = c(0, 80),
    breaks = seq(0, 120, by = 10)
  ) + theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 10)
  )

ggsave("modelfrequencyaic2.png", plot = p2, width = 8, height = 6, dpi = 300)
ggsave("modelfrequencyaic2.pdf", plot = p2, width = 8, height = 6)

## Let´s make just one dual column bar figure for these two.

dualcol <- rbind(best_models,bestmodelsaiccorrected)
dualcol$selectionaic <- c(rep("AIC mín.", 100), rep("Delta AIC > 2", 100))

p2 <-
ggplot(dualcol, aes(x = reorder(modelname, modelname,
                                    function(x)-length(x)), fill = selectionaic)) +
  geom_bar(width = 0.7, position = position_dodge(preserve = "single")) +
  xlab("Mejor modelo de diversificación") +
  ylab("Frecuencia absoluta") +
  scale_x_discrete(labels = c("minimal" = "Nulo",
                              "variable.lambda" = "lambda variable",
                              "variable.mu" = "mu variable",
                              "full.model" = "Completo")) +
  theme_classic() +
  scale_fill_grey(start = 0.8, end = 0.3, name = NULL) +
  scale_y_continuous(
    limits = c(0, 80),
    breaks = seq(0, 120, by = 10)
  ) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 10)
  )

ggsave("newmodelfrequencydual.png", plot = p2, width = 8, height = 6, dpi = 300)
ggsave("newmodelfrequencydual.pdf", plot = p2, width = 8, height = 6)


## subdivide plotting area ------------------------------------

minimalno2 <- best_models[best_models$modelname=="minimal",]
lambdano2 <- best_models[best_models$modelname=="variable.lambda",]
muno2 <- best_models[best_models$modelname=="variable.mu",]
fullno2 <- best_models[best_models$modelname=="full.model",]

minimalwt2 <- bestmodelsaiccorrected[bestmodelsaiccorrected$modelname=="minimal",]
lambdawt2 <- bestmodelsaiccorrected[bestmodelsaiccorrected$modelname=="variable.lambda",]
muwt2 <- bestmodelsaiccorrected[bestmodelsaiccorrected$modelname=="variable.mu",]
fullwt2 <- bestmodelsaiccorrected[bestmodelsaiccorrected$modelname=="full.model",]

png("parameters.png", width = 2000, height = 1500, res = 300)
pdf("parameters.pdf", width = 8, height = 6)

par(mfrow=c(2,4), xpd = T)



### Ploting best models without the <2 correction

## a) plotconstant rate (CID) QuaSSE model
plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",#"Longitud de espina",
     ylab=expression(paste(lambda," o ",mu)))
mtext("(a)",line=1,adj=-0.7)
clip(min(spinedataready),max(spinedataready),0,1.5)
for(i in minimalno2$treenum) {
  
abline(h=spineresults[[i]]$null$par["lambda"],lwd=0.05)
abline(h=spineresults[[i]]$null$par["mu"],lwd=0.05,col="gray")
}
legend(x=-1.5, y= 0.4, legend = "n = 52", bty = "n" )

## b) plotvariable lambda QuaSSE model

plot(NULL,xlim=range(spinedataready),ylim=c(0,0.6),bty="n",
     xlab="",#"Longitud de espina",
     ylab="") #expression(paste(lambda," o ",mu))
mtext("(b)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,1.5)

for(i in lambdano2$treenum) {
abline(a=coef(spineresults[[i]]$b)["l.c"],
       b=coef(spineresults[[i]]$b)["l.m"],lwd=0.2)
abline(h=coef(spineresults[[i]]$b)["m.c"],lwd=0.2,
       col="gray")
}

legend(x=-1.5, y= 0.8, legend = "n = 28", bty = "n" )


## c) plot variable mu QuaSSE model
plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",
     ylab="")
mtext("(c)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,1.5)

for(i in muno2$treenum) {
abline(h=coef(spineresults[[i]]$d)["l.c"],lwd=0.05)
abline(a=coef(spineresults[[i]]$d)["m.c"],
       b=coef(spineresults[[i]]$d)["m.m"],lwd=0.05,col="gray")
}
legend(x=--1.5, y= 0.4, legend = "n = 19", bty = "n" )


## d) plot variable lambda and mu QuaSSE model

plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",
     ylab="")
mtext("(d)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,3)

for(i in fullno2$treenum) {
abline(a=coef(spineresults[[i]]$full)["l.c"],
       b=coef(spineresults[[i]]$full)["l.m"],lwd=2)
abline(a=coef(spineresults[[i]]$full)["m.c"],
       b=coef(spineresults[[i]]$full)["m.m"],lwd=2,col="gray")
}

legend(x=-1.5, y= 0.4, legend = "n = 1", bty = "n" )

### Ploting best models with the <2 correction

#par(mfrow=c(2,2))

## a) plotconstant rate (CID) QuaSSE model
plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",
     ylab=expression(paste(lambda," o ",mu)))
mtext("(e)",line=1,adj=-0.7)
clip(min(spinedataready),max(spinedataready),0,1.5)
for(i in minimalwt2$treenum) {
  
  abline(h=spineresults[[i]]$null$par["lambda"],lwd=0.05)
  abline(h=spineresults[[i]]$null$par["mu"],lwd=0.05,col="gray")
}

legend(x=-1.5, y= 0.4, legend = "n = 79", bty = "n" )

## b) plotvariable lambda QuaSSE model

plot(NULL,xlim=range(spinedataready),ylim=c(0,0.6),bty="n",
     xlab="",#"Longitud de espina",
     ylab="")
mtext("(f)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,1.5)

for(i in lambdawt2$treenum) {
  abline(a=coef(spineresults[[i]]$b)["l.c"],
         b=coef(spineresults[[i]]$b)["l.m"],lwd=0.2)
  abline(h=coef(spineresults[[i]]$b)["m.c"],lwd=0.2,
         col="gray")
  
}
legend(x=-1.5, y= 0.8, legend = "n = 9", bty = "n" )


## c) plot variable mu QuaSSE model
plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="Longitud de espina",
     ylab="")
mtext("(g)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,1.5)

for(i in muwt2$treenum) {
  abline(h=coef(spineresults[[i]]$d)["l.c"],lwd=0.05)
  abline(a=coef(spineresults[[i]]$d)["m.c"],
         b=coef(spineresults[[i]]$d)["m.m"],lwd=0.05,col="gray")
}
legend(x=-1.5, y= 0.4, legend = "n = 12", bty = "n" )

## d) plot variable lambda and mu QuaSSE model

plot(NULL,xlim=range(spinedataready),ylim=c(0,0.5), type = "n", axes = FALSE, ann = FALSE)
clip(min(spinedataready),max(spinedataready),0,1.5)

legend("topright",
       c(expression(lambda),expression(mu)),
       lwd=5,col=c("black","gray"),bty="n")

# plot(NULL,xlim=range(spinedataready),ylim=c(0,1.5),bty="n",
#      xlab="",
#      ylab="")
# mtext("(h)",line=1,adj=-1)
# clip(min(spinedataready),max(spinedataready),0,3)
# 
# for(i in fullwt2$treenum) {
#   abline(a=coef(spineresults[[i]]$full)["l.c"],
#          b=coef(spineresults[[i]]$full)["l.m"],lwd=2)
#   abline(a=coef(spineresults[[i]]$full)["m.c"],
#          b=coef(spineresults[[i]]$full)["m.m"],lwd=2,col="gray")
# }
# 
# legend(x=-1.5, y= 2, legend = "n = 4", bty = "n" )


dev.off()


### diversification and turnover rates


transform <- function(m,b,otherRate,inputrate,outputrate){
 
   # extracting 2 points from lambda or mu
  
  pointy1 <- m*2+b
  pointy2 <- m*5+b

  # transform the points into diversification or turnover

  if(outputrate == "d") {# diversification
    if(inputrate == "l") { # inputlambda
      newpointy1 <- pointy1 - otherRate
      newpointy2 <- pointy2 - otherRate
      
    } else if(inputrate == "m") {  # inpútmu 
      newpointy1 <- otherRate - pointy1  
      newpointy2 <- otherRate - pointy2
    }
  } else if(outputrate == "t") { # turnover
      newpointy1 <- pointy1 + otherRate
      newpointy2 <- pointy2 + otherRate
  }
  
  # calculating slope of diversification or turnover
  
  slope <- (newpointy2-newpointy1)/(5-2)
  return(slope)
  
}




png("diversiparameters.png", width = 2000, height = 1500, res = 300)
pdf("diversiparameters.pdf", width = 8, height = 6)

par(mfrow=c(2,4))

## a) plotconstant div rate (CID) QuaSSE model
plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",#"Longitud de espina",
     ylab= "Tasa de diversificación neta", cex.lab=1.5)
mtext("(a)",line=1,adj=-0.7)
clip(min(spinedataready),max(spinedataready),0,1.5)
for(i in minimalno2$treenum) {
  
  abline(h= (spineresults[[i]]$null$par["lambda"]- spineresults[[i]]$null$par["mu"]),lwd=0.05,col="red")
  
}
#legend(x=-1.5, y= 0.4, legend = "n = 52", bty = "n" )

## b) plotvariable lambda QuaSSE model

plot(NULL,xlim=range(spinedataready),ylim=c(0,0.6),bty="n",
     xlab="",#"Longitud de espina",
     ylab="") #expression(paste(lambda," o ",mu))
mtext("(b)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,1.5)

for(i in lambdano2$treenum) {
  abline(a=(coef(spineresults[[i]]$b)["l.c"] - coef(spineresults[[i]]$b)["m.c"]),
         b=transform(coef(spineresults[[i]]$b)["l.m"],coef(spineresults[[i]]$b)["l.c"]
                     ,coef(spineresults[[i]]$b)["m.c"],"l","d"),lwd=0.2, col = "red")
}

#legend(x=-1.5, y= 0.8, legend = "n = 28", bty = "n" )



## c) plot variable mu QuaSSE model
plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",
     ylab="")
mtext("(c)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,1.5)

for(i in muno2$treenum) {
  abline(a=(coef(spineresults[[i]]$d)["l.c"] - coef(spineresults[[i]]$d)["m.c"]),
         b=transform(coef(spineresults[[i]]$d)["m.m"],coef(spineresults[[i]]$b)["m.c"]
                     ,coef(spineresults[[i]]$d)["l.c"],"m","d"),lwd=0.2, col = "red")
  
  
}
#legend(x=--1.5, y= 0.4, legend = "n = 19", bty = "n" )

# d) plot variable lambda and mu QuaSSE model

plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",
     ylab="")
mtext("(d)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,3)

for(i in fullno2$treenum) {
  abline(a= (coef(spineresults[[i]]$full)["l.c"] - coef(spineresults[[i]]$full)["m.c"]),
         b= (coef(spineresults[[i]]$full)["l.m"] - coef(spineresults[[i]]$full)["m.m"]),lwd=2, col="red")
  
}

#legend(x=-1.5, y= 0.4, legend = "n = 1", bty = "n" )


## a) plotconstant div rate (CID) QuaSSE model
plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",#"Longitud de espina",
     ylab="Tasa de recambio", cex.lab=1.5)
mtext("(e)",line=1,adj=-0.7)
clip(min(spinedataready),max(spinedataready),0,1.5)
for(i in minimalno2$treenum) {
  
  
  abline(h= (spineresults[[i]]$null$par["lambda"] + spineresults[[i]]$null$par["mu"]),lwd=0.05,col="blue")
}
#legend(x=-1.5, y= 0.4, legend = "n = 52", bty = "n" )

## b) plotvariable lambda QuaSSE model

plot(NULL,xlim=range(spinedataready),ylim=c(0,0.6),bty="n",
     xlab="",#"Longitud de espina",
     ylab="") #expression(paste(lambda," o ",mu))
mtext("(f)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,1.5)

for(i in lambdano2$treenum) {
  
  abline(a=(coef(spineresults[[i]]$b)["l.c"] - coef(spineresults[[i]]$b)["m.c"]),
         b=transform(coef(spineresults[[i]]$b)["l.m"],coef(spineresults[[i]]$b)["l.c"]
                     ,coef(spineresults[[i]]$b)["m.c"],"l","t"),lwd=0.2, col = "blue")
}

#legend(x=-1.5, y= 0.8, legend = "n = 28", bty = "n" )



## c) plot variable mu QuaSSE model
plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",
     ylab="")
mtext("(g)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,1.5)

for(i in muno2$treenum) {
  
  
  abline(a=(coef(spineresults[[i]]$d)["l.c"] - coef(spineresults[[i]]$d)["m.c"]),
         b=transform(coef(spineresults[[i]]$d)["m.m"],coef(spineresults[[i]]$b)["m.c"]
                     ,coef(spineresults[[i]]$d)["l.c"],"m","t"),lwd=0.2, col = "blue")
}
#legend(x=--1.5, y= 0.4, legend = "n = 19", bty = "n" )

# d) plot variable lambda and mu QuaSSE model

plot(NULL,xlim=range(spinedataready),ylim=c(0,0.3),bty="n",
     xlab="",
     ylab="")
mtext("(h)",line=1,adj=-1)
clip(min(spinedataready),max(spinedataready),0,3)

for(i in fullno2$treenum) {

  
  abline(a= (coef(spineresults[[i]]$full)["l.c"] + coef(spineresults[[i]]$full)["m.c"]),
         b= (coef(spineresults[[i]]$full)["l.m"] + coef(spineresults[[i]]$full)["m.m"]),lwd=2, col="blue")
}

#legend(x=-1.5, y= 0.4, legend = "n = 1", bty = "n" )



dev.off()

