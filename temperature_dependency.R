## Testing time dependet diversification in a sample of trees

library(RPANDA)
library(phytools)

# Function for fitting multiple models to multiple trees

sampFit <- function(x, frac) {
  
  modelspertree <- list()
  
  for(i in 1:length(x)) {
    
    ##Fitting constant lambda and fixed mu YUle
    cat("Fitting models for tree number ",i,"\n")
    tot_time<-max(node.age(x[[i]])$ages)
    
    f.lamb<-function(t,y){y[1]}
    f.mu<-function(t,y){0}
    lamb_par_init<-c(0.05)
    mu_par_init<-c()
    cat("Starting yule model")
    yule <-fit_bd(x[[1]],tot_time,f.lamb,f.mu,lamb_par_init,
                  mu_par_init,f=frac,cst.lamb=TRUE,fix.mu=TRUE)
    print(yule)
    ##Fitting constant lambda and mu
    
    f.lamb<-function(t,y){y[1]}
    f.mu<-function(t,y){y[1]}
    lamb_par_init<-c(0.05)
    mu_par_init<-c(0.005)
    cat("Starting constant model")
    constant <-fit_bd(x[[i]],tot_time,f.lamb,f.mu,lamb_par_init,
                      mu_par_init,f=frac,cst.lamb=TRUE,cst.mu=TRUE)
    print(constant)
    #### fittimg variable lambda and constant mu
    
    ### exponential
    
    f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
    f.mu<-function(t,y){y[1]}
    lamb_par_init<-c(0.05,0.01)
    mu_par_init<-c(0.005)
    cat("Starting expoBvariable model")
    expoBvariable <-fit_bd(x[[i]],tot_time,f.lamb,f.mu,lamb_par_init,  
                           mu_par_init,f=frac,expo.lamb=TRUE,cst.mu=TRUE)
    print(expoBvariable)
    #### lineal 
    
    f.lamb.lin<-function(t,y){y[1]+y[2]*t}
    f.mu<-function(t,y){y[1]}
    lamb_par_init<-c(0.05,0.01)
    mu_par_init<-c(0.005)
    cat("Starting linearBvariable model")
    linearBvariable <-fit_bd(x[[i]],tot_time,f.lamb.lin,f.mu,lamb_par_init,    
                             mu_par_init,f=frac,expo.lamb=FALSE,cst.mu=TRUE)
    print(linearBvariable)
    ###### #### fittimg constant lambda and variable mu
    
    ### exponential
    
    f.lamb<-function(t,y){y[1]}
    f.mu<- function(t,y){y[1]*exp(y[2]*t)}
    lamb_par_init<-c(0.05)
    mu_par_init<-c(0.005, 0.01)
    cat("Starting expoDvariable model")
    expoDvariable <-fit_bd(x[[i]],tot_time,f.lamb,f.mu,lamb_par_init,
                           mu_par_init,f=frac,cst.lamb=TRUE,expo.mu=TRUE)
    print(expoDvariable)
    #### lineal 
    
    f.lamb<-function(t,y){y[1]}
    f.mu.lin<-function(t,y){y[1]+y[2]*t}
    lamb_par_init<-c(0.05)
    mu_par_init<-c(0.005,0.01)
    cat("Starting linearDvariable model")
    linearDvariable <-fit_bd(x[[i]],tot_time,f.lamb,f.mu.lin,lamb_par_init,
                             mu_par_init,f=frac,cst.lamb=TRUE, expo.mu = FALSE)
    print(linearDvariable)
    #### fitting both variables
    
    ### Exponential
    
    f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
    f.mu<- function(t,y){y[1]*exp(y[2]*t)}
    lamb_par_init<-c(0.05, 0.01)
    mu_par_init<-c(0.005, 0.01)
    cat("Starting expoBDvariable model")
    expoBDvariable <-fit_bd(x[[i]],tot_time,f.lamb,f.mu,lamb_par_init,
                            mu_par_init,f=frac,expo.lamb=TRUE,expo.mu=TRUE)
    print(expoBDvariable)
    ### linear  
    
    f.lamb<-function(t,y){y[1]+y[2]*t}
    f.mu<- function(t,y){y[1]+y[2]*t}
    lamb_par_init<-c(0.05, 0.01)
    mu_par_init<-c(0.005, 0.01)
    cat("Starting linearBDvariable model")
    linearBDvariable <-fit_bd(x[[i]],tot_time,f.lamb,f.mu,lamb_par_init,
                              mu_par_init,f=frac)
    print(linearBDvariable)
    models<-list(y = yule, ct = constant, eb = expoBvariable, lb = linearBvariable, ed = expoDvariable, 
                 ld = linearDvariable , ebd = expoBDvariable, lbd = linearBDvariable)
    
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

saveRDS(results, "results.rds")

results <- readRDS("results.rds")
