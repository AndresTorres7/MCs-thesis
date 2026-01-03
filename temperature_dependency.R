## Testing time dependet diversification in a sample of trees

library(RPANDA)
library(phytools)

## Temperature data 
#Zachos, Dickens & Zeebe 2008
data(InfTemp)

# Function for fitting multiple models to multiple trees

sampFit <- function(x, frac) {
  
  modelspertree <- list()
  data(InfTemp)
  
  for(i in 1:length(x)) {
    
    cat("Fitting models for tree number ",i,"\n")
    tot_time<-max(node.age(x[[i]])$ages)
    
    #### fitting lambda dependent on temperature and mu fixed
    # expo 
    
    f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
    f.mu<-function(t,x,y){0}
    lamb_par_init<-c(0.10,0.01)
    mu_par_init<- c()
    cat("Starting mu expo fixed model")
    t_exp_mufixed <-fit_env(x[[i]],InfTemp,tot_time,f.lamb,f.mu,
                            lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)
    print(t_exp_mufixed)
    
    #linear
    f.lamb<-function(t,x,y){y[1]+y[2]*x}
    f.mu<-function(t,x,y){0}
    lamb_par_init<-c(0.10,0.01)
    mu_par_init<- c()
    cat("Starting mu linear fixed model")
    t_lin_mufixed <-fit_env(x[[i]],InfTemp,tot_time,f.lamb,f.mu,
                              
                              lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)
    
    print(t_lin_mufixed)
    
    #### fitting lambda dependent on temperature and mu constant
    
    f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
    f.mu<-function(t,x,y){y[1]}
    lamb_par_init<-c(0.10,0.01)
    mu_par_init<- c(0.01)
    cat("Starting expoBTemp model")
    expoBTemp <-fit_env(x[[i]],InfTemp,tot_time,f.lamb,f.mu,
                        
                        lamb_par_init,mu_par_init,f=1,cst.mu=TRUE,dt=1e-3)
    
    print(expoBTemp)
    #linear
    
    f.lamb<-function(t,x,y){y[1]+y[2]*x}
    f.mu<-function(t,x,y){y[1]}
    lamb_par_init<-c(0.10,0.01)
    mu_par_init<- c(0.01)
    cat("Starting linearBTemp model")
    linearBTemp <-fit_env(x[[i]],InfTemp,tot_time,f.lamb,f.mu,
                          
                          lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)
    
    print(linearBTemp)
    #### fitting lambda constant  and mu dependent on temperature
    
    # expo
    
    f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
    f.lamb<-function(t,x,y){y[1]}
    lamb_par_init<-c(0.01)
    mu_par_init<- c(0.10,0.01)
    cat("Starting expoDTemp model")
    expoDTemp <-fit_env(x[[i]],InfTemp,tot_time,f.lamb,f.mu,
                        
                        lamb_par_init,mu_par_init,f=1,cst.lamb=TRUE,dt=1e-3)
    print(expoDTemp)
    #linear
    
    f.lamb<-function(t,x,y){y[1]}
    f.mu<-function(t,x,y){y[1]+y[2]*x}
    lamb_par_init<-c(0.01)
    mu_par_init<- c(0.10,0.01)
    cat("Starting linearDTemp model")
    linearDTemp <-fit_env(x[[i]],InfTemp,tot_time,f.lamb,f.mu,
                          
                          lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)
    
    print(linearDTemp)
    #### fitting lambda and mu dependent on temperature
    
    #expo 
    
    f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
    f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
    lamb_par_init<-c(0.10,0.01)
    mu_par_init<- c(0.10,0.01)
    cat("Starting expoBDTemp model")
    expoBDTemp <-fit_env(x[[i]],InfTemp,tot_time,f.lamb,f.mu,
                         
                         lamb_par_init,mu_par_init,f=1,dt=1e-3)
    print(expoBDTemp)
    
    #Linear
    
    f.lamb<-function(t,x,y){y[1]+y[2]*x}
    f.mu<-function(t,x,y){y[1]+y[2]*x}
    lamb_par_init<-c(0.10,0.01)
    mu_par_init<- c(0.10,0.01)
    cat("Starting linearBTemp model")
    linearBDTemp <-fit_env(x[[i]],InfTemp,tot_time,f.lamb,f.mu,
                           
                           lamb_par_init,mu_par_init,f=1,dt=1e-3)
    
    print(linearBDTemp)
    
    
    
    models<-list(em0 = t_exp_mufixed, lm0 = t_lin_mufixed, teb = expoBTemp, tlb = linearBTemp,
                 ted =  expoDTemp, tld = linearDTemp , tebd = expoBDTemp, tlbd = linearBDTemp)
    
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

saveRDS(results, "temperatureresults.rds")

results <- readRDS("temperatureresults.rds")
