

library(TreeTools)
library(phytools)
library(dplyr)

scripts <- c("rangeDispersal.R", "nicheEvolution.R","speciateAllopatric.R","speciateSympatric.R",
             "speciateParapatric.R","speciateDispersal.R","seedSpecies.R","environmentalChange.R",
             "nicheRecenter.R","DREaD.R","generateSummaryStatistics.R", "helperFunctions.R","findSisters.R","summary_statsitics_functions.R")



for (i in 1:length(scripts)){
  scripts[[i]] <- paste("~/Biologia/diversificationTest/Code/AuxiliaryScripts/DREAD",scripts[[i]],sep="/")
}
lapply(scripts, source)

rangeDispersal <- disperseRange
required.packages <- (c("raster","gstat", "SpaDES", "ape","phytools","geiger",
                        "phyloclim","ggplot2","gridExtra","moments",
                        "apTreeshape","parallel", "doSNOW", "rgeos","knitr",
                        "data.table", "fossil", "ENMTools"))

lapply(required.packages, require, character.only=T)

setwd("~/Biologia/Maestria/Tesis/Final/disparidad morfológica")

# Read TNT matrix
morpho_matrix <-  ReadTNTCharacters("S17_morphological_matrix.tnt")

# Read tree

setwd("~/Biologia/Maestria/Tesis/Final/speciationmode")

tree <- read.tree("micra.tre")
tree$tip.label <- sub("^(([^_]+_[^_]+)).*", "\\1", tree$tip.label)
tree$tip.label <- gsub("_",".",tree$tip.label)

clade.trees <- list(genero = tree) 
names(clade.trees) <- "Micrathena"

# Sisters pairs

findSisters <- function(phy,solve.polytomies=FALSE,include.polytomies=TRUE){
  require(phytools)
  n <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (n < 4)
    stop("tree has fewer than 4 tips")
  if(solve.polytomies){
    phy=multi2di(phy)
    include.polytomies=FALSE
  }
  if(include.polytomies){
    sister.ancestors <- which(tabulate(phy$edge[, 1][phy$edge[, 2] <= n])>=2)
    cherries <- vector("list",length(sister.ancestors))
    for(i in 1:length(cherries)){
      x <- getDescendants(phy, sister.ancestors[i])
      cherries[[i]] <- phy$tip.label[x]
      if(length(which(is.na(cherries[[i]])))>0) cherries[[i]]<-NULL  #remove internal nodes
    }
  }
  if(!include.polytomies){
    sister.ancestors <- which(tabulate(phy$edge[, 1][phy$edge[, 2] <= n])==2)
    cherries <- vector("list",length(sister.ancestors))
    for(i in 1:length(cherries)){
      x <- getDescendants(phy, sister.ancestors[i])
      cherries[[i]] <- c(phy$tip.label[x[1]], phy$tip.label[x[2]])
    }
  }
  cherries[!sapply(cherries,is.null)]
  return(cherries)
}


sisters <- findSisters(tree)

# alopatric and nonalopatric

#make the first part of "inference of speciation / inference of geographic modes 
#of speciation script till RO_ss

RO_ss <- sisterSpeciesOverlap(sp.rasters, sisters)
tabla <- as.data.frame(RO_ss)
tabla <- cbind(SpeciesPairs = rownames(tabla), tabla)
rownames(tabla) <- NULL
tabla$distribution <- ifelse(tabla$RO_ss == 0, "allopatric", "nonallopatric")
tabla$realdistribution <- case_when(tabla$RO_ss == 0 ~ "allopatric", tabla$RO_ss == 1 ~ "sympatric",
                                    0 < tabla$RO_ss & tabla$RO_ss< 1  ~ "parapatric")


# Extract all matches
matches <- regmatches(tabla$SpeciesPairs, gregexpr("\"([A-Za-z0-9_]+)\"", tabla$SpeciesPairs, perl = TRUE))
matches <- unlist(lapply(matches, function(m) gsub("\"", "", m)))
matches <- gsub("_",".",matches)

tabla$pair1 <-  matches[seq(1,74, 2)]
tabla$pair2 <-  matches[seq(2,74, 2)]
tabla$totaldisimi <- NA
tabla$somaticdisim <- NA
tabla$sexdisimi<-NA
tabla$age <- NA



# Sexual vs somatic characters
setwd("~/Biologia/Maestria/Tesis/Final/disparidad morfológica")

sexinfo <- read.csv("somaticvssexual.csv")
sexinfo$Character <- 1:192

sexcharacter <- morpho_matrix[,sexinfo$Type=="se"]
somacharacter<- morpho_matrix[,sexinfo$Type=="so"]



# branching.times returns ages of nodes (distance from tips)
bt <- branching.times(tree)

for(i in 1:nrow(tabla) ) {

v1 <- morpho_matrix[tabla$pair1[i],]
v2 <-morpho_matrix[tabla$pair2[i],]
v1[v1 %in% c("-", "?")] <- NA
v2[v2 %in% c("-", "?")] <- NA

keep <- !is.na(v1) & !is.na(v2) 
v1 <- v1[keep]
v2 <- v2[keep]

tabla$totaldisimi[i] <- sum(v1 != v2)/length(v1)

# sexual characters

vse1 <- sexcharacter[tabla$pair1[i],]
vse2 <- sexcharacter[tabla$pair2[i],]

vse1[vse1 %in% c("-", "?")] <- NA
vse2[vse2 %in% c("-", "?")] <- NA

keep <- !is.na(vse1) & !is.na(vse2) 
vse1 <- vse1[keep]
vse2 <- vse2[keep]

tabla$sexdisimi[i] <- sum(vse1 != vse2)/length(vse1)

# somatic characters

vso1 <- somacharacter[tabla$pair1[i],]
vso2 <- somacharacter[tabla$pair2[i],]

vso1[vso1 %in% c("-", "?")] <- NA
vso2[vso2 %in% c("-", "?")] <- NA

keep <- !is.na(vso1) & !is.na(vso2) 
vso1 <- vso1[keep]
vso2 <- vso2[keep]

tabla$somaticdisim[i] <- sum(vso1 != vso2)/length(vso1)

# Extracting species pair age.

pairs <- c(tabla$pair1[i], tabla$pair2[i])

# 1. Get the MRCA node number
mrca_node <- getMRCA(tree, pairs)



# Age of MRCA
tabla$age[i] <- bt[as.character(mrca_node)]


}

write.csv(tabla, "datossimilitudages.csv")

tabla$diff <- tabla$somaticdisim - tabla$sexdisimi

# Organizing data for lmer analysis

sexdata <- tabla[,c(1,3,8,9,11)]
somadata <- tabla[,c(1,3,7,9)]
totaldata <- tabla[,c(1,3,7,8,9)]

#boxplots
library(cowplot)

p1 <- ggplot(tabla, aes(x=distribution, y=totaldisimi)) +
  geom_boxplot() +
  ylab("Total disparity")

p2<-ggplot(tabla, aes(x=distribution, y=somaticdisim)) +
  geom_boxplot() +
  ylab("Somatic disparity")

p3<-ggplot(tabla, aes(x=distribution, y=sexdisimi)) +
  geom_boxplot() +
  ylab("Sexual disparity")

p4<-ggplot(tabla, aes(x=distribution, y=diff)) +
  geom_boxplot() +
  ylab("Disparity difference somatic - sexual")

plot_grid(p1, p2, p3,p4, labels = c('A', 'B','C','D'), label_size = 12)

# scatterplots


p1<- ggplot(tabla, aes(x=RO_ss, y=totaldisimi)) +
  geom_point()+
  ylab("Total disparity")


p2 <-ggplot(tabla, aes(x=RO_ss, y=somaticdisim)) +
  geom_point()+
  ylab("Somatic disparity")

p3<-ggplot(tabla, aes(x=RO_ss, y=sexdisimi)) +
  geom_point()+
  ylab("Sexual disparity")

p4<-ggplot(tabla, aes(x=RO_ss, y=diff)) +
  geom_point()+
  ylab("Disparity difference somatic-sexual") +
  theme(axis.title.y = element_text(size = 10))

plot_grid(p1, p2, p3,p4, labels = c('A', 'B','C','D'), label_size = 12)


#scatterplots with time

ggplot(tabla, aes(x=age, y=totaldisimi)) +
  geom_point()

ggplot(tabla, aes(x=age, y=somaticdisim)) +
  geom_point()

ggplot(tabla, aes(x=age, y=sexdisimi)) +
  geom_point()

ggplot(tabla, aes(x=age, y=diff)) +
  geom_point()



# prueba de supuestos de modelo lineal


library(car) # para el test de durbin-watson 
library(psych) # para la función pair panels
library(lmtest) # para el test bptest

cor.test(tabla$RO_ss, tabla$totaldisimi)
cor.test(tabla$RO_ss, tabla$sexdisimi)
cor.test(tabla$RO_ss, tabla$somaticdisim)

cor.test(tabla$age, tabla$totaldisimi)
cor.test(tabla$age, tabla$sexdisimi)
cor.test(tabla$age, tabla$somaticdisim)

cor.test(tabla$RO_ss, tabla$diff)


# Without age as a factor, the model dont pass the normality test.

# Ro and age
modelo <- lm(sexdisimi ~ RO_ss , data = tabla)
modelo <- lm(sexdisimi ~ RO_ss + age, data = tabla)
modelo <- lm(sexdisimi ~ RO_ss * age, data = tabla)
summary(modelo)


# discrete distribution and age
modelo <- lm(totaldisimi ~ distribution, data = tabla)

modelo <- lm(sexdisimi ~ distribution + age, data = tabla)

modelo <- lm(sexdisimi ~ distribution * age, data = tabla)

summary(modelo)




#  diff - discrete distribution and age
modelo <- lm(diff ~ distribution, data = tabla)

modelo <- lm(diff ~ distribution + age, data = tabla)

modelo <- lm(diff ~ distribution * age, data = tabla)

summary(modelo)


mean(modelo$residuals)
shapiro.test(modelo$residuals)

qqPlot(modelo$residuals, 
       distribution = "norm",
       main = "Q-Q PLOT de residuos en mpg ~ wt",
       xlab = "cuantiles teóricos",
       ylab = "cuantiles de la muestra",
       id = FALSE, grid = TRUE,
       envelope = 0.95, col = carPalette()[1], col.lines = carPalette()[3],
       pch = 20,
       cex = 1,
       lwd = 2)
qqline(modelo$residuals,
       col = "blue",
       lty = 1,
       lwd = 2)  

bptest(modelo)


# generalized model, beta distribution is better for values among 0 and 1

# values cant be exactly 0 nor 1

N <- nrow(tabla)
tabla$sexdisimitrans <- (tabla$sexdisimi * (N - 1) + 0.5) / N


library(betareg)

# Basic model: effect of RO_ss
modelo1 <- betareg(sexdisimitrans ~ RO_ss, data = tabla)

# Model with an additional factor (e.g. distribution)
modelo2 <- betareg(sexdisimitrans ~ RO_ss + age, data = tabla)

# Model with interaction
modelo3 <- betareg(sexdisimitrans ~ RO_ss * age, data = tabla)


# Summaries
summary(modelo1)
summary(modelo2)
summary(modelo3)


## GLM


# Install required packages if not installed
# install.packages(c("DHARMa", "statmod", "tweedie", "cplm", "car"))

library(DHARMa)
library(statmod)
library(tweedie)
library(cplm)
library(car)

# ---- 1. Estimate optimal Tweedie variance power ----
# This profiles the likelihood to find the best var.power (p)
profile <- tweedie.profile(
  sexdisimi ~ distribution + age,
  data = tabla,
  p.vec = seq(1.1, 2, by = 0.1),   # search range for p
  do.plot = TRUE
)

best_p <- profile$p.max
cat("Estimated Tweedie variance power:", best_p, "\n")
library(statmod)

modelo1 <- glm(totaldisimi ~ distribution ,
               data = tabla,
               family = tweedie(var.power = 1.1, link.power = 0))

modelo2 <- glm(totaldisimi ~ distribution + age,
              data = tabla,
              family = tweedie(var.power = 1.1, link.power = 0))

modelo3 <- glm(sexdisimi ~ distribution * age,
               data = tabla,
               family = tweedie(var.power = 1.1, link.power = 0))

modelo <- glm(sexdisimi ~  age,
               data = tabla,
               family = tweedie(var.power = 1.1, link.power = 0))

summary(modelo)
summary(modelo2)
summary(modelo3)


# 1. Residuals vs fitted
plot(fitted(modelo), residuals(modelo), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

# 2. QQ-plot of residuals
qqnorm(residuals(modelo))
qqline(residuals(modelo))

# 3. Scale-location (check variance)
sqrt_abs_resid <- sqrt(abs(residuals(modelo)))
plot(fitted(modelo), sqrt_abs_resid, 
     xlab = "Fitted values", ylab = "Sqrt(|Residuals|)")

summary(modelo2)$dispersion
