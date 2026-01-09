# Figure and results of time dependent diversification


library(ggplot2)


# Data loading

datos <- readRDS("resultscorrected.rds")

##### Selecting best model based only in Aic values in every tree #########-------------------------------

a<-list()
for(i in 1:length(datos)) {

  aic_values <- c(datos[[i]]$y$aicc, datos[[i]]$ct$aicc, datos[[i]]$eb$aicc, datos[[i]]$lb$aicc, datos[[i]]$ed$aicc,
                  datos[[i]]$ld$aicc, datos[[i]]$ebd$aicc, datos[[i]]$lbd$aicc)
  
  a[[i]] <- datos[[i]][aic_values == min(aic_values)]

}

modelos <- character()

for(i in 1:length(a)) { 

  modelos <- c(modelos, names(a[[i]])) 
  }

models <- data.frame(Tree = 1:100, Best_Model = modelos)

# frequency table

ggplot(models, aes(x = reorder(Best_Model, Best_Model,
                               function(x)-length(x)))) +
  geom_bar(width = 0.7) +
  xlab("Mejor modelo de diversificación") +
  ylab("Frecuencia absoluta") +
  scale_x_discrete(labels = c("lbd" = "Completo(lin)",
                               "y" = "Yule",
                               "lb" = "Lambda/Mu_var(lin)",
                              "ed" = "Mu_var(expo)")) +
  theme_classic() +
  scale_fill_grey(start = 0.8, end = 0.3)




## Selection best model based on both aic and eexclusion of negative rates ######_-----------------------------

kth_smallest <- function(x, k) {
  sort(x, na.last = NA)[k]
}

b<-list()

datos <- readRDS("resultscorrected.rds")

for(i in 1:length(datos)) {
aic_values <- c(datos[[i]]$y$aicc, datos[[i]]$ct$aicc, datos[[i]]$eb$aicc, datos[[i]]$lb$aicc, datos[[i]]$ed$aicc,
                                 datos[[i]]$ld$aicc, datos[[i]]$ebd$aicc, datos[[i]]$lbd$aicc)

aic <- min(aic_values)
name <- names(datos[[i]][aic_values == aic])
if (!(name %in% c("y","ct"))) {
  
  l <- datos[[i]][aic_values == aic]
  r <- min(l[[name]]$lamb_par[1],l[[name]]$mu_par[1])
  z<-1
  while (r<0) {
    
    z<-z+1
  aic <- kth_smallest(aic_values,z)
  name <- names(datos[[i]][aic_values == aic])
  
  if (!(name %in% c("y","ct"))) {
    
    l <- datos[[i]][aic_values == aic]
    r <- min(l[[name]]$lamb_par[1],l[[name]]$mu_par[1])
  } else {
    r<-1
  }
  
  }
  
  
  b[[i]] <- datos[[i]][aic_values == kth_smallest(aic_values,z)]
} else {
b[[i]] <- datos[[i]][aic_values == min(aic_values)]
}
}

modelos <- character()

for(i in 1:length(b)) { 
  
  modelos <- c(modelos, names(b[[i]])) 
}

models <- data.frame(Tree = 1:100, Best_Model = modelos)

# frequency table

ggplot(models, aes(x = reorder(Best_Model, Best_Model,
                               function(x)-length(x)))) +
  geom_bar(width = 0.7) +
  xlab("Mejor modelo de diversificación") +
  ylab("Frecuencia absoluta") +
  # scale_x_discrete(labels = c("lbd" = "Completo(lin)",
  #                             "y" = "Yule",
  #                             "lb" = "Lambda/Mu_var(lin)",
  #                             "ed" = "Mu_var(expo)")) +
  theme_classic() +
  scale_fill_grey(start = 0.8, end = 0.3)


