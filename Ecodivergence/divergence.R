## FUnctions ---------------------------------------

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

scripts <- c("rangeDispersal.R", "nicheEvolution.R","speciateAllopatric.R","speciateSympatric.R",
             "speciateParapatric.R","speciateDispersal.R","seedSpecies.R","environmentalChange.R",
             "nicheRecenter.R","DREaD.R","generateSummaryStatistics.R", "helperFunctions.R","findSisters.R","summary_statsitics_functions.R")



for (i in 1:length(scripts)){
  scripts[[i]] <- paste("~/Biologia/diversificationTest/Code/AuxiliaryScripts/DREAD",scripts[[i]],sep="/")
}
lapply(scripts, source)

rangeDispersal <- disperseRange

## libraries-------------------------------------------

library(ggplot2)
library(car)
library(MASS)
library(psych)
library(cowplot)
require(ape)
require(diverge)
require(PBSmapping)
require(parallel)
library(TreeTools)
library(phytools)
library(dplyr)
required.packages <- (c("raster","gstat", "SpaDES", "ape","phytools","geiger",
                        "phyloclim","ggplot2","gridExtra","moments",
                        "apTreeshape","parallel", "doSNOW", "rgeos","knitr",
                        "data.table", "fossil", "ENMTools"))

lapply(required.packages, require, character.only=T)


## Datos------------------------------------------------

#tree  I think the best will be to use only one tree in this case

micratree <- read.tree("micra.tre")
micratree<-drop.tip(micratree,"M_horrida_122") # M horrida is repeated
micratree$tip.label <- sub("_", ".", micratree$tip.label)
micratree$tip.label <- sub("\\_.*", "", micratree$tip.label)
micratree$tip.label <- sub("M.", "Micrathena_", micratree$tip.label)

# newclade.trees <- list()
# 
# for (i in seq_along(micratree)) {
#   tr <- micratree[[i]]
#   
#   tr<-drop.tip(tr,"M_horrida_122") # M horrida is repeated
#   tr$tip.label <- sub("_", ".", tr$tip.label)
#   tr$tip.label <- sub("\\_.*", "", tr$tip.label)
#   
#   newclade.trees[[i]] <- tr
#   
# }
# 
# names(newclade.trees) <- names(micratree)
# 
# micratree <- newclade.trees



# range rasters
load("MicraRasters.rda")

sp.rasters <- k
rm(k)
for( f in 1:length(sp.rasters)) {
  
  origin(sp.rasters[[f]])<-origin(sp.rasters[[1]])
  
  
}
## Selections sister species  ### Not necesary to doit again, chechk the sisterSpecies.rds archive >50 -------------------------

# 
# spairs <- list()
# 
# for(i in seq_along(micratree)) {
#   
#   sisters <- findSisters(micratree[[i]])
#   a <- character()
#   
#   for(g in seq_along(sisters)) {
#     
#     pair <- paste(sisters[[g]][1],sisters[[g]][2], sep = "_")
#     
#     a <- c(a,pair)
#     
#   }
#   
#   spairs[[i]] <- a
#   
# }
# 
# #ALL in one
# all <- character()
# for(i in seq_along(micratree)) {
#   
#   sisters <- findSisters(micratree[[i]])
#  
#   
#   for(g in seq_along(sisters)) {
#     
#     pair <- paste(sisters[[g]][1],sisters[[g]][2], sep = "_")
#     
#     all <- c(all,pair)
#     
#   }
#   
# }
# 
# speciespairs <- table(all)
# 
# betstsuportsp <- speciespairs[speciespairs > 50]
# 
# sisters <- strsplit(names(betstsuportsp), "_")
# 
# saveRDS(sisters, "sisterSpecies.rds") 
# 
# sisters <- readRDS("sisterSpecies.rds")

#making the basic data set ------------------------------------------------------

# Find sister in only one tree ----------------------------------------------



basicdata <- data.frame(sp1=character(), sp2=character(), pair_age = numeric(), overlap = numeric(), patry = character())


pairs <- extract_sisters(micratree, sis_age = T)

# alopatric and nonalopatric ----------------------------------

#make the first part of "inference of speciation / inference of geographic modes 
#of speciation script till RO_ss

phy <- micratree
root.age.true <- max(findMRCA(micratree, type=c("height")))
micratree$edge.length<-micratree$edge.length/root.age.true
root.age<-1
sisters <- findSisters(micratree, solve.polytomies=T)
RO_ss <- sisterSpeciesOverlap(sp.rasters, sisters)

tabla <- as.data.frame(RO_ss)
tabla <- cbind(SpeciesPairs = rownames(tabla), tabla)
rownames(tabla) <- NULL
tabla$distribution <- ifelse(tabla$RO_ss <=  0.25, "Alopatrico", "No alopatrico")

pairs$overlap <- tabla$RO_ss
pairs$patry <- tabla$distribution

write.csv(pairs, "paresespecies.csv",row.names = F)


### Disparidad morfológica discreta

# Read TNT matrix
morpho_matrix <-  ReadTNTCharacters("S17_morphological_matrix.tnt")
morpho_matrix <- morpho_matrix[-c(1:13),]
rownames(morpho_matrix)<-gsub("M.","Micrathena_", rownames(morpho_matrix), fixed = T)  

sexinfo <- read.csv("somaticvssexual.csv")
sexinfo$Character <- 1:192

sexcharacter <- morpho_matrix[,sexinfo$Type=="se"]
somacharacter<- morpho_matrix[,sexinfo$Type=="so"]

for(i in 1:nrow(pairs) ) {

v1 <- morpho_matrix[pairs$sp1[i],]
v2 <-morpho_matrix[pairs$sp2[i],]

#Keeping only non NA, missing or no aplica values
v1[v1 %in% c("-", "?")] <- NA
v2[v2 %in% c("-", "?")] <- NA
keep <- !is.na(v1) & !is.na(v2) 
v1 <- v1[keep]
v2 <- v2[keep]

# Calculating disimilitud total

pairs$totaldisimi[i] <- sum(v1 != v2)/length(v1)
 
# only sexual characters

vse1 <- sexcharacter[pairs$sp1[i],]
vse2 <- sexcharacter[pairs$sp2[i],]

vse1[vse1 %in% c("-", "?")] <- NA
vse2[vse2 %in% c("-", "?")] <- NA

keep <- !is.na(vse1) & !is.na(vse2) 
vse1 <- vse1[keep]
vse2 <- vse2[keep]

pairs$sexdisimi[i] <- sum(vse1 != vse2)/length(vse1)

# somatic characters

vso1 <- somacharacter[pairs$sp1[i],]
vso2 <- somacharacter[pairs$sp2[i],]

vso1[vso1 %in% c("-", "?")] <- NA
vso2[vso2 %in% c("-", "?")] <- NA

keep <- !is.na(vso1) & !is.na(vso2) 
vso1 <- vso1[keep]
vso2 <- vso2[keep]

pairs$somaticdisim[i] <- sum(vso1 != vso2)/length(vso1)




}


write.csv(pairs, "datosdisimilitudes.csv", row.names = F)
pairs <- read.csv("datosdisimilitudes.csv")

pairs$patry[pairs$patry== "Alopatrico"] <- "Alopátrico"
pairs$patry[pairs$patry== "No alopatrico"] <- "No alopátrico"



#boxplots

ylims <- range(
  pairs$totaldisimi,
  pairs$somaticdisim,
  pairs$sexdisimi,
  na.rm = TRUE
)



p1 <- ggplot(pairs, aes(x=patry, y=totaldisimi)) +
  geom_boxplot() +
  scale_y_continuous(limits = ylims) +
  ylab("Disparidad") +
  xlab("") +
  theme_classic() 

p2<-ggplot(pairs, aes(x=patry, y=somaticdisim)) +
  geom_boxplot() +
  scale_y_continuous(limits = ylims) +
  ylab("") +
  xlab("Solape de rangos de distribución") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
  )

p3<-ggplot(pairs, aes(x=patry, y=sexdisimi)) +
  geom_boxplot() +
  scale_y_continuous(limits = ylims) +
  ylab("") +
  xlab("") +
   theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
  )

setwd("figures")

png("disparidadmorfo.png", width = 7, height = 4, units = "in", res = 600)

 plot_grid(p1, p2, p3, labels = c('A', 'B','C'), label_size = 12, nrow = 1)

dev.off()





  library(DHARMa)
  library(statmod)
  library(tweedie)
  library(cplm)
  library(car)

tabla <- read.csv("datosdisimilitudes.csv")


# ---- 1. Estimate optimal Tweedie variance power ----
# This profiles the likelihood to find the best var.power (p)
profile <- tweedie.profile(
  totaldisimi ~ overlap + pair_age,
  data = tabla,
  p.vec = seq(1.1, 2, by = 0.1),   # search range for p
  do.plot = TRUE
)

best_p <- profile$p.max
cat("Estimated Tweedie variance power:", best_p, "\n")
library(statmod)

modelo1t <- glm(totaldisimi ~ overlap ,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))

modelo2t <- glm(totaldisimi ~ overlap + pair_age,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))

modelo2t <- glm(sexdisimi ~ overlap + pair_age,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))

modelo2t <- glm(somaticdisim ~ overlap + pair_age,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))

modelo3t <- glm(totaldisimi ~ overlap * pair_age,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))

modelot <- glm(totaldisimi ~  pair_age,
               data = tabla,
               family = tweedie(var.power = 1.1, link.power = 0))

summary(modelo1t)
summary(modelo2t)
summary(modelo3t)
summary(modelot)


profiles <- tweedie.profile(
  sexdisimi ~ overlap + pair_age,
  data = tabla,
  p.vec = seq(1.1, 2, by = 0.1),   # search range for p
  do.plot = TRUE
)


best_p <- profiles$p.max
cat("Estimated Tweedie variance power:", best_p, "\n")
library(statmod)

modelo1s <- glm(sexdisimi ~ overlap ,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))

modelo2s <- glm(sexdisimi ~ overlap + pair_age,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))

modelo3s <- glm(sexdisimi ~ overlap * pair_age,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))

modelos <- glm(sexdisimi ~  pair_age,
               data = tabla,
               family = tweedie(var.power = 1.1, link.power = 0))

summary(modelo1s)
summary(modelo2s)
summary(modelo3s)
summary(modelos)



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



#### Modelos en tesis
library(statmod)

profile1 <- tweedie.profile(
  totaldisimi ~ overlap + pair_age,
  data = tabla,
  p.vec = seq(1.1, 2, by = 0.1),   # search range for p
  do.plot = TRUE
)

best_p1 <- profile1$p.max
cat("Estimated Tweedie variance power:", best_p1, "\n")

modelo1 <- glm(totaldisimi ~ overlap * pair_age,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))






profile2 <- tweedie.profile(
  sexdisimi ~ overlap + pair_age,
  data = tabla,
  p.vec = seq(1.1, 2, by = 0.1),   # search range for p
  do.plot = TRUE
)


best_p2 <- profile2$p.max
cat("Estimated Tweedie variance power:", best_p2, "\n")



modelo2 <- glm(sexdisimi ~ overlap * pair_age,
                data = tabla,
                family = tweedie(var.power = 1.1, link.power = 0))







profile3 <- tweedie.profile(
  somaticdisim ~ overlap + pair_age,
  data = tabla,
  p.vec = seq(1.1, 2, by = 0.1),   # search range for p
  do.plot = TRUE
)


best_p3 <- profile3$p.max
cat("Estimated Tweedie variance power:", best_p3, "\n")



modelo3 <- glm(somaticdisim ~ overlap * pair_age,
               data = tabla,
               family = tweedie(var.power = 1.1, link.power = 0))

summary(modelo1)
summary(modelo3)
summary(modelo2)































## Nicho climático

datos <- read.csv("climaticdivergencedata.csv", row.names = 1)

png("nichegeooverlap.png", width = 5, height = 4, units = "in", res = 600)

ggplot(datos, aes(x = RO_ss, y = nicheOverlap, color = pair_age)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", name = "Edad del Par") +
  labs(title = "",
       x = "Superposición Geográfica",
       y = "Superposición de Nicho") +
  theme_minimal() +
  theme(legend.position = "right")

dev.off()

png("boxplot.png", width = 5, height = 5, units = "in", res = 600)

boxplot(datos$nicheOverlap ~ datos$distribution, 
        main = "",
        xlab = "", ylab = "Superposición de nicho", 
        names = c("Alopátrico", "No Alopátrico")) 

dev.off()

# Matriz de correlaciones
cor_matrix <- cor(datos[, c("RO_ss", "nicheOverlap", "pair_age")], 
                  use = "complete.obs")
print(cor_matrix)

# Pruebas de correlación
cor.test(datos$RO_ss, datos$nicheOverlap)
cor.test(datos$pair_age, datos$nicheOverlap)
cor.test(datos$RO_ss, datos$pair_age)

## Supuestos pearson

# Normalidad

# Pruebas de normalidad para cada variable
cat("=== PRUEBAS DE NORMALIDAD (Shapiro-Wilk) ===\n")
cat("RO_ss: W =", shapiro.test(datos$RO_ss)$statistic, 
    "p =", shapiro.test(datos$RO_ss)$p.value, "\n")
cat("nicheOverlap: W =", shapiro.test(datos$nicheOverlap)$statistic, 
    "p =", shapiro.test(datos$nicheOverlap)$p.value, "\n")
cat("pair_age: W =", shapiro.test(datos$pair_age)$statistic, 
    "p =", shapiro.test(datos$pair_age)$p.value, "\n")

# Gráficos QQ-plot
par(mfrow = c(1, 3))
qqPlot(datos$RO_ss, main = "QQ-plot: RO_ss", ylab = "RO_ss")
qqPlot(datos$nicheOverlap, main = "QQ-plot: nicheOverlap", ylab = "nicheOverlap")
qqPlot(datos$pair_age, main = "QQ-plot: pair_age", ylab = "pair_age")

# Histogramas con curva normal
par(mfrow = c(1, 3))
hist(datos$RO_ss, main = "Distribución de RO_ss", freq = FALSE, 
     xlab = "RO_ss", col = "lightblue")
curve(dnorm(x, mean = mean(datos$RO_ss), sd = sd(datos$RO_ss)), 
      add = TRUE, col = "red", lwd = 2)

hist(datos$nicheOverlap, main = "Distribución de nicheOverlap", freq = FALSE, 
     xlab = "nicheOverlap", col = "lightgreen")
curve(dnorm(x, mean = mean(datos$nicheOverlap), sd = sd(datos$nicheOverlap)), 
      add = TRUE, col = "red", lwd = 2)

hist(datos$pair_age, main = "Distribución de pair_age", freq = FALSE, 
     xlab = "pair_age", col = "lightcoral")
curve(dnorm(x, mean = mean(datos$pair_age), sd = sd(datos$pair_age)), 
      add = TRUE, col = "red", lwd = 2)



datos$log_niche <- log1p(datos$nicheOverlap)  # log(1 + x)
shapiro.test(datos$log_niche)

# Cargar datos
datos <- read.csv("climaticdivergencedata.csv")

##### pearson correlation que está en la tesis de la maestria



# Pearson con Bootstrapping
library(boot)

# Función para bootstrapping
boot_cor <- function(data, indices) {
  muestra <- data[indices, ]
  return(cor(muestra$RO_ss, muestra$nicheOverlap, method = "pearson"))
}

# Aplicar bootstrapping
set.seed(123)
resultado_boot <- boot(datos, boot_cor, R = 2000)

# Calcular intervalo de confianza 95%
ic <- boot.ci(resultado_boot, type = "all")

# Mostrar resultados
cat("=== RESULTADOS PEARSON CON BOOTSTRAP ===\n")
cat("Correlación (r):", round(resultado_boot$t0, 3), "\n")
cat("IC 95%: [", round(ic$percent[4], 3), ", ", round(ic$percent[5], 3), "]\n\n")

# Decisión: ¿Hay correlación?
if (ic$percent[4] > 0 && ic$percent[5] > 0) {
  cat("✅ SÍ hay correlación positiva significativa\n")
  cat("   (El intervalo de confianza no incluye 0)\n")
} else if (ic$percent[4] < 0 && ic$percent[5] < 0) {
  cat("✅ SÍ hay correlación negativa significativa\n")
  cat("   (El intervalo de confianza no incluye 0)\n")
} else {
  cat("❌ NO hay correlación significativa\n")
  cat("   (El intervalo de confianza incluye 0)\n")
}















######################## Artículo ------------------------------------------

# FPS = First posterior spine, SPS = Second posterior spine, BL = Body length, FL = Femur length

micra.dat <- data.frame(sp1 = character(), sp2 = character(), pair_age = numeric(), sp1FPS = numeric(), 
                        sp1SPS = nuemric(), sp1BL = numeric(), sp1FL = numeric(), sp2FPS = numeric(), 
                        sp2SPS = nuemric(), sp2BL = numeric(), sp2FL = numeric(), overlap = numeric() )


# Load the dataset
adams.dat= read.csv("adams_salampairs_formatted.csv")

# vars to be tested
colnames(adams.dat)
# first the values for which the difference needs to be calculated
testvars = c("SVL", "TL", "HL", "S.E", "BW", "FLL", "HLL", "cor_logpc1", "cor_logpc2", "cov_logpc1", "cov_logpc2")
# second the values that are already differences
vars2d = c("cor_log_2d", "cov_log_2d")

# Trim dataset to just the allopatric pairs 
adams.allos = adams.dat[which(adams.dat$overlap < 20),] #27

## PART A. MODEL FITTING WITHOUT LATITUDE

# MODEL FITTING 1: PCs
fits = mclapply(testvars[-c(1:7)], function(x) modelfitter(datset=adams.allos, varb=x, mods=models))
names(fits) = testvars[-c(1:7)]

# MODEL FITTING 2: log-transformed linear variables
logfits = mclapply(testvars[1:7], function(x) modelfitter(datset=adams.allos, varb=x, logs=TRUE, me.levels=meas.err, mods=models))
names(logfits) = testvars[1:7]

# MODEL FITTING 3: 2D EUCLIDEAN DISTANCES
fits2d = mclapply(vars2d, function(x) modelfitter(datset=adams.allos, varb=x, is.difference=TRUE, mods=models))
names(fits2d) = vars2d

# SUMMARIZE MODEL FITS
adamsfits = summarizer(taxa="salamanders", dataset="adams2009", reslist=fits, transf=NA)
adamslog_fits = summarizer(taxa="salamanders", dataset="adams2009", reslist=logfits, transf="log")
adams2d_fits = summarizer(taxa="salamanders", dataset="adams2009", reslist=fits2d, transf=NA)
adams_modelfit_summary = rbind(adamsfits, adamslog_fits, adams2d_fits)
# Save long-form results
write.csv(adams_modelfit_summary, file="adams_modelfit_summary.csv", row.names=FALSE)
any(adams_modelfit_summary$convergence>0)

# SHORTFORM SUMMARY
adams_modelfit_short = res.reporter(summ.mat=adams_modelfit_summary, traits=c(testvars,vars2d))
adams_modelfit_short
# Save short-form results
write.csv(adams_modelfit_short, file="adams_modelfit_short.csv", row.names=FALSE)

# VISUALIZE MODEL FITS
# raw traits
# cov = "lat.diff"
ylim=NULL
var.index = 8 # only use 8 to 11; will fail otherwise
fit.comparer(res=adams_modelfit_summary, dataset=adams.allos, trait=testvars[var.index], ylim=ylim)

# log-transformed traits
ylim=NULL
var.index = 5 # only use 1 to 7; will fail otherwise
fit.comparer(res=adams_modelfit_summary, dataset=adams.allos, transf="log", trait=testvars[var.index], ylim=ylim)

# 2d traits
ylim=NULL
var.index = 1 # only use 1 to 2; will fail otherwise
fit.comparer(res=adams_modelfit_summary, dataset=adams.allos, is.diff=TRUE, trait=vars2d[var.index], ylim=ylim)




