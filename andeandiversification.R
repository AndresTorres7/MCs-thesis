## Testing andean dependet diversification in a sample of trees

library(RPANDA)
library(phytools)
library(pspline)

## Andean data 
#(Lagomarsino et al., 2016


env <- read.table("AndeanElevation.txt", header = T)

env_mean <- aggregate(
  Altitude ~ Age,
  data = env,
  FUN = mean
)



# Fit smoothing spline
ps <- smooth.Pspline(
  x = env_mean$Age,
  y = env_mean$Altitude,
  norder = 2
)

paleo_fun <- function(t) {
  predict(ps, x = t)
}

time_grid <- seq(
  from = 0,
  to   = max(env_mean$Age),
  length.out = 1000
)

env_RPANDA <- data.frame(
  Time = time_grid,
  Env  = predict(ps, x = time_grid)
)

write.csv(env_RPANDA, "andeselevaciondatos.csv", col.names = F)
env_RPANDA <- read.csv("andeselevaciondatos.csv")

#Visualize the smooth curve (important sanity check)

plot(
  env_mean$Age, env_mean$Altitude,
  pch = 16, col = "gray",
  xlab = "Time (Ma)", ylab = "Andean elevation (m)",
  xlim = rev(range(env_mean$Age))
)

curve(
  paleo_fun(x),
  from = min(env_mean$Age),
  to   = max(env_mean$Age),
  add = TRUE, col = "black", lwd = 2
)




# Function for fitting multiple models to multiple trees

sampFit <- function(x, frac) {
  
  modelspertree <- list()
  
  
  for(i in 1:length(x)) {
    
    cat("Fitting models for tree number ",i,"\n")
    tot_time<-max(node.age(x[[i]])$ages)
    
    #### fitting lambda dependent on andes and mu fixed
    # expo 
    
    f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
    f.mu<-function(t,x,y){0}
    lamb_par_init<-c(0.1,0.001)
    mu_par_init<- c()
    cat("Starting mu expo fixed model")
    a_expo_mufixed <-fit_env(x[[i]],env_RPANDA,tot_time,f.lamb,f.mu,
                           lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)
    
    print(a_expo_mufixed)
    
    #linear
    
    f.lamb<-function(t,x,y){y[1]+y[2]*x}
    f.mu<-function(t,x,y){0}
    lamb_par_init<-c(0.1,0.001)
    mu_par_init<- c()
    cat("Starting mu linear fixed model")
    a_lin_mufixed <-fit_env(x[[i]],env_RPANDA,tot_time,f.lamb,f.mu,
                             lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)
    
    print(a_lin_mufixed)
    
    #### fitting lambda dependent on andes and mu constant
    #expo
    
    f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
    f.mu<-function(t,x,y){y[1]}
    lamb_par_init<-c(0.10,0.001)
    mu_par_init<- c(0.01)
    cat("Starting expoBAnd model")
    expoBAnd <-fit_env(x[[i]],env_RPANDA,tot_time,f.lamb,f.mu,
                       lamb_par_init,mu_par_init,f=1,cst.mu=TRUE,dt=1e-3)
    
    print(expoBAnd)
    
    #linear
    
    f.lamb<-function(t,x,y){y[1]+y[2]*x}
    f.mu<-function(t,x,y){y[1]}
    lamb_par_init<-c(0.10,0.001)
    mu_par_init<- c(0.01)
    cat("Starting linearBAnd model")
    linearBAnd <-fit_env(x[[i]],env_RPANDA,tot_time,f.lamb,f.mu,
                         lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)
    
    print(linearBAnd)
    
    #### fitting lambda constant  and mu dependent on andes
    
    # expo
    
    f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
    f.lamb<-function(t,x,y){y[1]}
    lamb_par_init<-c(0.1)
    mu_par_init<- c(0.01,0.001)
    cat("Starting expoDAnd model")
    expoDAnd <-fit_env(x[[i]],env_RPANDA,tot_time,f.lamb,f.mu,
                       lamb_par_init,mu_par_init,f=1,cst.lamb=TRUE,dt=1e-3)
   
    print(expoDAnd)
    
    #linear
    
    f.lamb<-function(t,x,y){y[1]}
    f.mu<-function(t,x,y){y[1]+y[2]*x}
    lamb_par_init<-c(0.1)
    mu_par_init<- c(0.01,0.0001)
    cat("Starting linearDAnd model")
    linearDAnd <-fit_env(x[[i]],env_RPANDA,tot_time,f.lamb,f.mu,
                         lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)
   
    print(linearDAnd)
    
    #### fitting lambda and mu dependent on temperature
    
    #expo 
    
    f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
    f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
    lamb_par_init<-c(0.10,0.001)
    mu_par_init<- c(0.10,0.001)
    cat("Starting expoBDAnd model")
    expoBDAnd <-fit_env(x[[i]],env_RPANDA,tot_time,f.lamb,f.mu,
                        lamb_par_init,mu_par_init,f=1,dt=1e-3)
    
    print(expoBDAnd)
    
    #Linear
    
    f.lamb<-function(t,x,y){y[1]+y[2]*x}
    f.mu<-function(t,x,y){y[1]+y[2]*x}
    lamb_par_init<-c(0.10,0.001)
    mu_par_init<- c(0.10,0.001)
    cat("Starting linearBDAnd model")
    linearBDAnd <-fit_env(x[[i]],env_RPANDA,tot_time,f.lamb,f.mu,
                          lamb_par_init,mu_par_init,f=1,dt=1e-3)
    
    
    print(linearBDAnd)
    
    
    
    models<-list(Aem0 = a_expo_mufixed, Alm0 = a_lin_mufixed, aeb = expoBAnd, alb = linearBAnd,
                 aed =  expoDAnd, ald = linearDAnd , aebd = expoBDAnd, albd = linearBDAnd)
    
    modelspertree[[i]] <- models
    gc()
  }
  
  return(modelspertree)
}


# Data entry

micratree <- readNexus("k.BEAST-100-posterior-trees.nex")
samplingFraction <- (120/150) #fraction of extant species included in the phylogeny Pedro says are 150
#115

results <- sampFit(micratree, samplingFraction)

saveRDS(results, "andesresults.rds")

results <- readRDS("andesresults.rds")