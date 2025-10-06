
# sister species
setwd("~/Biologia/Maestria/Tesis/Final/disparidad morfológica")
sisters <- read.csv("datossimilitud.csv") 
sisters <- sisters[, c(3,4,5,6)]
sisters$pair1 <- gsub("M.", "M_", sisters$pair1)
sisters$pair2 <- gsub("M.", "M_", sisters$pair2)
sisters <- sisters[-1,]

# data points 
setwd("D:/climatic niche micrathena")
poins <- read.csv("alloccurrences.csv", row.names = 1)
names(poins)[1] <- "species"


# only species with more thatn 4 points

newsisters <- c(sisters$pair1,sisters$pair2)
newpoins <- poins[poins$species %in% newsisters,]
cant <- table(newpoins$species)
cant <- cant[cant>=4]
sisters <- sisters[sisters$pair1 %in% names(cant) & sisters$pair2 %in% names(cant),]

for(i in 1:nrow(sisters) ) {

  setwd("D:/climatic niche micrathena/occurrences")
  dir.create(paste(i,sisters$pair1[i], sep = "_"))
  setwd(paste(i,sisters$pair1[i], sep = "_"))
  sppoints <- poins[poins$species == sisters$pair1[i],]
  write.csv(x=sppoints, paste(sisters$pair1[i], "csv", sep="."),row.names = F)
  
  setwd("D:/climatic niche micrathena/occurrences")
  dir.create(paste(i,sisters$pair2[i],sep = "_"))
  setwd(paste(i,sisters$pair2[i],sep = "_"))
  sppoints <- poins[poins$species == sisters$pair2[i],]
  write.csv(x=sppoints, paste(sisters$pair2[i], "csv", sep="."), row.names = F)
  
}

# data visualization
climaticdata <- sisters[-c(24, 18, 13, 4),]
climaticdata$nicheOverlap <- NA

for(i in 1:length(climaticdata$pair2)) {
  
  setwd("D:/climatic niche micrathena/occurrences")
  setwd(list.files(pattern = climaticdata$pair2[i]))

  overdata <- read.csv("overlap_result.csv")
  climaticdata$nicheOverlap[i] <- overdata$Overlap_MVE    
  
  
}

setwd("D:/climatic niche micrathena")
write.csv(climaticdata, "climaticdivergencedata.csv")
climaticdata <- read.csv("climaticdivergencedata.csv") 


ggplot(climaticdata, aes(x = distribution, y = nicheOverlap)) +
  geom_boxplot()

ggplot(climaticdata, aes(x = RO_ss, y = nicheOverlap)) +
  geom_point()


correlation <- cor.test(climaticdata$RO_ss, climaticdata$nicheOverlap, method = "p", exact = F)

# statistical analysis

setwd("~/Biologia/Maestria/Tesis/Final/disparidad morfológica")
tabla <- read.csv("datossimilitud.csv")
ages <- tabla[,c(6,10)]
ages$pair2 <- gsub("M.","M_", ages$pair2)
merged <- merge(climaticdata, ages, by = "pair2", all.x = TRUE)

modelo <- lm(nicheOverlap ~ RO_ss , data = climaticdata)
summary(modelo)

modelo <- lm(nicheOverlap ~ RO_ss * age , data = merged)
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






