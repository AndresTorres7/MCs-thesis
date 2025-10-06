
library(RPANDA)

micrathena<- read.tree("micra.tre")
tot_time<-max(node.age(micrathena)$ages)


##Fitting constant lambda and fixed mu YUle

f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){0}

lamb_par_init<-c(0.05)
mu_par_init<-c()

yule <-fit_bd(micrathena,tot_time,f.lamb,f.mu,lamb_par_init,
                  
                  mu_par_init,f=1,cst.lamb=TRUE,fix.mu=TRUE)

##Fitting constant lambda and mu

f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y[1]}

lamb_par_init<-c(0.05)
mu_par_init<-c(0.005)

constant <-fit_bd(micrathena,tot_time,f.lamb,f.mu,lamb_par_init,
            
            mu_par_init,f=1,cst.lamb=TRUE,cst.mu=TRUE)


res1<-fit_bd(micrathena,tot_time,f.lamb,f.mu,lamb_par_init,
            
            mu_par_init,f=1)


#### fittimg variable lambda and constant mu

### exponential

 f.lamb<-function(t,y){y[1]*exp(y[2]*t)}

 f.mu<-function(t,y){y[1]}
 
 lamb_par_init<-c(0.05,0.01)
 mu_par_init<-c(0.005)

 expoBvariable <-fit_bd(micrathena,tot_time,f.lamb,f.mu,lamb_par_init,
                           
                           mu_par_init,f=1,expo.lamb=TRUE,cst.mu=TRUE)



 #### lineal 
 
 f.lamb.lin<-function(t,y){y[1]+y[2]*t}
 
 f.mu<-function(t,y){y[1]}
 
 lamb_par_init<-c(0.05,0.01)
 mu_par_init<-c(0.005)
 
 linearBvariable <-fit_bd(micrathena,tot_time,f.lamb.lin,f.mu,lamb_par_init,
                        
                        mu_par_init,f=1,expo.lamb=FALSE,cst.mu=TRUE)
 
 
 

###### #### fittimg constant lambda and variable mu
 
 ### exponential
 
 f.lamb<-function(t,y){y[1]}
 
 f.mu<- function(t,y){y[1]*exp(y[2]*t)}
 
 lamb_par_init<-c(0.05)
 mu_par_init<-c(0.005, 0.01)
 
 expoDvariable <-fit_bd(micrathena,tot_time,f.lamb,f.mu,lamb_par_init,
                        
                        mu_par_init,f=1,cst.lamb=TRUE,expo.mu=TRUE)
 
 
 
 
 #### lineal 
 
 f.lamb<-function(t,y){y[1]}
 
 f.mu.lin<-function(t,y){y[1]+y[2]*t}
 
 lamb_par_init<-c(0.05)
 mu_par_init<-c(0.005,0.01)
 
 linearDvariable <-fit_bd(micrathena,tot_time,f.lamb,f.mu.lin,lamb_par_init,
                          
                          mu_par_init,f=1,cst.lamb=TRUE, expo.mu = FALSE)
 
 
 
 #### fitting both variables
 
 ### Exponential
 
 f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
 
 f.mu<- function(t,y){y[1]*exp(y[2]*t)}
 
 lamb_par_init<-c(0.05, 0.01)
 mu_par_init<-c(0.005, 0.01)
 
 expoBDvariable <-fit_bd(micrathena,tot_time,f.lamb,f.mu,lamb_par_init,
                        
                        mu_par_init,f=1,expo.lamb=TRUE,expo.mu=TRUE)
 
 ### linear   ###falta
 
 f.lamb<-function(t,y){y[1]+y[2]*t}
 
 f.mu<- function(t,y){y[1]+y[2]*t}
 
 lamb_par_init<-c(0.05, 0.01)
 mu_par_init<-c(0.005, 0.01)
 
 linearBDvariable <-fit_bd(micrathena,tot_time,f.lamb,f.mu,lamb_par_init,
                         
                         mu_par_init,f=1)
 

 
 
 
 ##### Fitting temperature 
 
 
 data(InfTemp)


#### fitting lambda dependent on temperature and mu fixed

 ## exponential
 
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}


lamb_par_init<-c(0.10,0.01)

mu_par_init<- c()

expoyuleBTemp <-fit_env(micrathena,InfTemp,tot_time,f.lamb,f.mu,
             
             lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)


expoBTemp$lamb_par[1]

expoyuleBTemp$lamb_par[2]

plot_fit_env(expoyuleBTemp,InfTemp,tot_time)


#### fitting lambda dependent on temperature and mu fixed

#linear
 
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){0}


lamb_par_init<-c(0.10,0.01)

mu_par_init<- c()

linearyuleBTemp <-fit_env(micrathena,InfTemp,tot_time,f.lamb,f.mu,
                        
                        lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)



#### fitting lambda dependent on temperature and mu constant

f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}


lamb_par_init<-c(0.10,0.01)

mu_par_init<- c(0.01)

expoBTemp <-fit_env(micrathena,InfTemp,tot_time,f.lamb,f.mu,
                        
                        lamb_par_init,mu_par_init,f=1,cst.mu=TRUE,dt=1e-3)


#### fitting lambda dependent on temperature and mu constant

#linear

f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]}


lamb_par_init<-c(0.10,0.01)

mu_par_init<- c(0.01)

linearBTemp <-fit_env(micrathena,InfTemp,tot_time,f.lamb,f.mu,
                          
                          lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)


#### fitting lambda constant  and mu dependent on temperature

# expo
  
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}

f.lamb<-function(t,x,y){y[1]}

lamb_par_init<-c(0.01)

mu_par_init<- c(0.10,0.01)

expoDTemp <-fit_env(micrathena,InfTemp,tot_time,f.lamb,f.mu,
                    
                    lamb_par_init,mu_par_init,f=1,cst.lamb=TRUE,dt=1e-3)


res$lamb_par[1]

res$lamb_par[2]

#### fitting lambda constant  and mu dependent on temperature
#linear

f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]+y[2]*x}


lamb_par_init<-c(0.01)

mu_par_init<- c(0.10,0.01)

linearDTemp <-fit_env(micrathena,InfTemp,tot_time,f.lamb,f.mu,
                      
                      lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)

#### fitting lambda and mu dependent on temperature

#expo 

f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}

lamb_par_init<-c(0.10,0.01)

mu_par_init<- c(0.10,0.01)

expoBDTemp <-fit_env(micrathena,InfTemp,tot_time,f.lamb,f.mu,
                    
                    lamb_par_init,mu_par_init,f=1,dt=1e-3)


#### fitting lambda and mu dependent on temperature

#Linear

f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]+y[2]*x}

lamb_par_init<-c(0.10,0.01)

mu_par_init<- c(0.10,0.01)

linearBDTemp <-fit_env(micrathena,InfTemp,tot_time,f.lamb,f.mu,
                     
                     lamb_par_init,mu_par_init,f=1,dt=1e-3)





#### fitting Andes
install.packages("pspline")
library(pspline)

paleodata <- read.table("Lagomarsino_et_al_PastAndeanElevation.txt", header  = T)

papi <- smooth.Pspline(paleodata$Age, paleodata$Altitude) 




#### fitting lambda dependent on andes and mu fixed

# expo 

f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}


lamb_par_init<-c(0.1,0.001)

mu_par_init<- c()

expoyuleBAnd <-fit_env(micrathena,paleodata,tot_time,f.lamb,f.mu,
                        
                        lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)


#### fitting lambda dependent on andes and mu fixed

#linear

f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){0}


lamb_par_init<-c(0.1,0.001)

mu_par_init<- c()

linearyuleBAnd <-fit_env(micrathena,paleodata,tot_time,f.lamb,f.mu,
                          
                          lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)


#### fitting lambda dependent on andes and mu constant
# expo

f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}


lamb_par_init<-c(0.10,0.001)

mu_par_init<- c(0.01)

expoBAnd <-fit_env(micrathena,paleodata,tot_time,f.lamb,f.mu,
                    
                    lamb_par_init,mu_par_init,f=1,cst.mu=TRUE,dt=1e-3)



#### fitting lambda dependent on andes and mu constant

#linear

f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]}


lamb_par_init<-c(0.10,0.001)

mu_par_init<- c(0.01)

linearBAnd <-fit_env(micrathena,paleodata,tot_time,f.lamb,f.mu,
                      
                      lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)


#### fitting lambda constant  and mu dependent on andes

#expo

f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}

f.lamb<-function(t,x,y){y[1]}

lamb_par_init<-c(0.1)

mu_par_init<- c(0.01,0.001)

expoDAnd <-fit_env(micrathena,paleodata,tot_time,f.lamb,f.mu,
                    
                    lamb_par_init,mu_par_init,f=1,cst.lamb=TRUE,dt=1e-3)



#### fitting lambda constant  and mu dependent on andes
#linear

f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]+y[2]*x}



lamb_par_init<-c(0.1)

mu_par_init<- c(0.01,0.0001)

linearDAnd <-fit_env(micrathena,paleodata,tot_time,f.lamb,f.mu,
                      
                      lamb_par_init,mu_par_init,f=1,fix.mu=TRUE,dt=1e-3)


#### fitting lambda and mu dependent on andes

#expo

f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}

lamb_par_init<-c(0.10,0.001)

mu_par_init<- c(0.10,0.001)

expoBDAnd <-fit_env(micrathena,paleodata,tot_time,f.lamb,f.mu,
                     
                     lamb_par_init,mu_par_init,f=1,dt=1e-3)



#### fitting lambda and mu dependent on andes

#Linear

f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]+y[2]*x}

lamb_par_init<-c(0.10,0.001)

mu_par_init<- c(0.10,0.001)

linearBDAnd <-fit_env(micrathena,paleodata,tot_time,f.lamb,f.mu,
                       
                       lamb_par_init,mu_par_init,f=1,dt=1e-3)





model <- c("Constant", "BvarDcst", "BvarDcst", "BcstDvar","BcstDvar" , "BvarDvar",
           "BtempDfixed", "BtempDcst", "BcstMtemp", "BtempDtemp", "BandDfixed", "BandDcst",
           "BcstMtemp", "BandDand")

mode <- c("-", "Exponential", "linear", "Exponential", "Linear", "Exponential" ,
          "Exponential", "Exponential", "Exponential", "Exponential", "Exponential", "Exponential",
          "Exponential", "Exponential")

lh <- c(constant$LH, expoBvariable$LH, linearBvariable$LH, expoDvariable$LH,
        linearDvariable$LH, expoBDvariable$LH, expoyuleBTemp$LH, linearyuleBTemp$LH,
        expoBTemp$LH, expoDTemp$LH, expoBDTemp$LH, expoyuleBAnd$LH,  expoBAnd$LH, expoDAnd$LH,
        expoBDAnd$LH)

Aic <- c(constant$aicc, expoBvariable$aicc, linearBvariable$aicc, expoDvariable$aicc,
         linearDvariable$aicc, expoBDvariable$aicc, expoyuleBTemp$aicc, linearyuleBTemp$aicc,
         expoBTemp$aicc, expoDTemp$aicc, expoBDTemp$aicc, expoyuleBAnd$aicc, expoBAnd$aicc,
         expoDAnd$aicc, expoBDAnd$aicc)



datos <-data.frame(Model = model, ModeofDependency = mode , logL = lh , AICc = Aic)

install.packages("writexl")
library(writexl)

write_xlsx(datos, "modelcomparison.xlsx")

