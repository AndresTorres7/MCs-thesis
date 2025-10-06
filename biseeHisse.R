


## Chacracter dependent diversification ####
library(readxl)  
library(geiger) 
library(phytools)
library(diversitree)
library(geiger)
library(phytools)

setwd("~/Biologia/Maestria/Tesis/Avances/Diversification/characterDependentDiversification")

micratree <- read.tree("micra.tre")
micratree<-drop.tip(micratree,"M_horrida_122") # M horrida is repeated

#apdatos<- read_xlsx("absentpresentdata.xlsx")
#colnames(apdatos)<-gsub(",", ".",colnames(apdatos)) 
#write.csv(apdatos, "apdatos.csv", quote = F)


apdatos<- read.csv("apdatos.csv", row.names = 2)
apdatos <- apdatos[,-1]

apindex <- read.table("absent_present characters.txt", sep = ";")
apindex <- apindex[,-2]
apindex <- gsub(" .*$", "", apindex)
apindex <- substr(apindex, 2,nchar(apindex))
apindex <- as.numeric(apindex) + 1

apdatos <- apdatos[-c(1:13),c(1,apindex)]  # remove outgroups and non binary characters

micratree$tip.label<-sub("_", ".", micratree$tip.label)
micratree$tip.label<-sub("\\_.*", "", micratree$tip.label)


name.check(phy = micratree,data = apdatos) 


apdatos <-apdatos[micratree$tip.label,] 





## extract habitat data
hab<-as.numeric(apdatospure[,])
## set names
names(hab)<-rownames(apdatospure)
hab[49] <- 2
## plot ourtree
plotTree(micratree,ftype="i",fsize=0.7,
         offset=0.5)
## add tiplabels
tiplabels(pie=to.matrix(hab,0:2)[micratree$tip.label,],
          piecol=c("red","blue", "white"),cex=0.4)
## create legend
legend("bottomleft",c("absent","present"),
       pch=21,pt.cex=1.6,
       cex=0.8,bty="n",
       pt.bg=c("red","blue"))

text(7,122, "Carapace..thoracic.region..lateral.dimples")

## make BiSSE likelihood function
bisse.model<-make.bisse(micratree,hab)
## find reasonable parameter values for
## optimization
p<-starting.point.bisse(micratree)
p  
# optimize BiSSE model
bisse.mle<-find.mle(bisse.model,p)
bisse.mle   

## create constrained null model
bissenull.model<-constrain(bisse.model,
                           lambda1~lambda0,mu1~mu0)
## optimize nullmodel
bissenull.mle<-find.mle(bissenull.model,
                        p[c(-2,-4)])
coef(bissenull.mle)  

# run likelihood-ratio test
bisseAnova<-anova(bisse.mle,
                  null=bissenull.mle)
bisseAnova

aicw(setNames(bisseAnova$AIC,
              rownames(bisseAnova)))

#


########## Checking what character has a signal

p<-starting.point.bisse(micratree)

datosmodelos <- data.frame(namecharacter = character(), character = integer(), AICnull = double(), AICfull = double(),
                           pr_chi = double())

##removing the next characters because they online have either 0 or 1 data

a<- apply(apdatos, 2, table)

conteo <- data.frame(name = character(), oin = logical(), iin = logical())

for (z in 1:length(a)) {
  cat(z)
  conteo[z,1] <- names(a[z])
  conteo[z,2] <- "0" %in% names(a[[z]])
  conteo[z,3] <- "1" %in% names(a[[z]])
  
}

conteo$junto <- conteo$oin & conteo$iin
apdatospure <- apdatos[,conteo$junto]
newa <- a[conteo$junto]




colnames(apdatos)[58]
apdatos[,12]



for (i in 1:ncol(apdatospure)) {
  
  character<-as.numeric(apdatospure[,i])
  names(character)<-rownames(apdatospure)
  datosmodelos[i,1] <- colnames(apdatospure)[i]
  datosmodelos[i,2] <- i
  cat(i)
  bisse.model<-make.bisse(micratree,character) 
  
  
  # optimize BiSSE model
  bisse.mle<-find.mle(bisse.model,p)
  ## create constrained null model
  bissenull.model<-constrain(bisse.model,
                             lambda1~lambda0,mu1~mu0)
  ## optimize nullmodel
  bissenull.mle<-find.mle(bissenull.model,
                          p[c(-2,-4)])
  
  # run likelihood-ratio test
  bisseAnova<-anova(bisse.mle,
                    null=bissenull.mle)
  
  datosmodelos[i,3] <- bisseAnova[2,3]
  datosmodelos[i,4] <- bisseAnova[1,3]
  datosmodelos[i,5] <- bisseAnova[2,5]
  
}

write.csv(datosmodelos, "valoresp_modelos_bisseVsNull_para_todos_los_caracteres.csv")

datosmodelos <-read.csv("valoresp_modelos_bisseVsNull_para_todos_los_caracteres.csv", row.names = 1)

significative <- datosmodelos[datosmodelos$pr_chi<0.05,]




#### Graficar modelos ##########################

library(ggplot2)

sig <- subset(datosmodelos, pr_chi < 0.05)
ggplot(datosmodelos, aes(character, pr_chi)) +
  geom_point() +
  geom_hline(yintercept = 0.05, color = "red", size = 1, linetype = "dashed")+
  geom_point(data = sig, colour = "blue", size = 2) +
  geom_text(data = sig[1,], colour = "red", label = "0.05", nudge_y = 0.05, nudge_x = -8)


##### mApear caracter 11 ###########



## extract habitat data
hab<-as.numeric(apdatospure[,11])
## set names
names(hab)<-rownames(apdatospure)
hab[49] <- 2
## plot ourtree
plotTree(micratree,ftype="i",fsize=0.6,
         offset=0.5)
## add tiplabels
tiplabels(pie=to.matrix(hab,0:2)[micratree$tip.label,],
          piecol=c("red","blue", "white"),cex=0.2, offset = 0.2)
## create legend
legend("bottomleft",c("absent","present"),
       pch=21,pt.cex=1.6,
       cex=0.8,bty="n",
       pt.bg=c("red","blue"))

text(17,122, datosmodelos$namecharacter[11])

##### graficar caracter 25



## extract habitat data
hab<-as.numeric(apdatospure[,25])
## set names
names(hab)<-rownames(apdatospure)
hab[!hab %in% c(0,1)] <- 2 
## plot ourtree
plotTree(micratree,ftype="i",fsize=0.6,
         offset=0.5)
## add tiplabels
tiplabels(pie=to.matrix(hab,0:2)[micratree$tip.label,],
          piecol=c("red","blue", "white"),cex=0.2, offset = 0.2)
## create legend
legend("bottomleft",c("absent","present"),
       pch=21,pt.cex=1.6,
       cex=0.8,bty="n",
       pt.bg=c("red","blue"))

text(17,122, datosmodelos$namecharacter[25])


table(apdatospure[,25])
table(apdatospure[,11])




##### Hisse ##################################

library(hisse)

## create input data frame for hisse
md<-data.frame(Genus.species=rownames(apdatospure),
               x=apdatospure[,11])
head(md)

## create HiSSE design matrix
rates.hisse<-TransMatMakerHiSSE(hidden.traits=1)
rates.hisse

## create hisse design matrix for BiSSE model
rates.bisse<-TransMatMakerHiSSE(hidden.traits=0)
## fit BiSSE model using hisse


md[49,2] <- 1 ## The hisse function doesn't allow 1&0 state for this species so I had to change it
bisse.hmle<-hisse(micratree,md,turnover=c(1,2),
                  eps=c(1,2),hidden.states=FALSE,
                  trans.rate=rates.bisse)

bisse.hmle

## custom function to back-transform turnover and
## extinction-fraction to lambda &mu
repar.bd<-function(object,k=2){
  pars<-object$solution
  tt<-pars[grep("turnover",names(pars))][1:k]
  ee<-pars[grep("eps",names(pars))][1:k]
  lambda<-tt/(1+ee)
  mu<-tt-lambda
  nn<-sapply(strsplit(names(tt),"turnover"),
             function(x) x[2])
  matrix(c(lambda,mu),k,2,dimnames=list(nn,
                                        c("lambda","mu")))
}
repar.bd(bisse.hmle)


## fitCID model using hisse  character independent diversification  = constant rate birth death model
cid.mle<-hisse(micratree,md,turnover=c(1,1),
               eps=c(1,1),hidden.states=FALSE,
               trans.rate=rates.bisse)
cid.mle  

## reparameterize CID model in terms oflambda
## andmu
repar.bd(cid.mle,1)


## create CID-2 design matrix
rates.cid2<-rates.hisse
rates.cid2[!is.na(rates.cid2)]<-1
rates.cid2


## fit CID-2 model using hisse
cid2.mle<-hisse(micratree,md,f=c(1,1),
                turnover=c(1,1,2,2),eps=c(1,1,2,2),
                hidden.states=TRUE,trans.rate=rates.cid2)

## reparameterize to lambda &mu
repar.bd(cid2.mle,k=4)

## obtain marginal reconstructions under CID-2 model
cid2.recon<-MarginReconHiSSE(phy=micratree,data=md,
                             f=cid2.mle$f,pars=cid2.mle$solution,
                             hidden.states=2)


## create a plotof the rates onthe tree
cid2.map<-plot.hisse.states(cid2.recon,
                            rate.param="speciation",
                            show.tip.label=TRUE,type="phylogram",
                            fsize=0.6,legend.position=c(0,0.3,0,0.3))

## createCID-2plotusingphytools
plot(setMap(cid2.map$rate.tree,c("white","red")),
     fsize=c(0.5,0.8),leg.txt="prob.(low/highspeciation)",
     dig=2)

## create design matrix for CID-4 model
rates.cid4<-TransMatMakerHiSSE(hidden.traits=3)
rates.cid4[!is.na(rates.cid4)]<-1
rates.cid4

## fitCID-4 model
cid4.mle<-hisse(micratree,md,f=c(1,1),
                turnover=c(1,1,2,2,3,3,4,4),
                eps=c(1,1,2,2,3,3,4,4),
                hidden.states=TRUE,
                trans.rate=rates.cid4)

## reparameterize model to lambda &mu
repar.bd(cid4.mle,8)

## fitfull HiSSE model
hisse.mle<-hisse(micratree,md,f=c(1,1),
                 hidden.states=TRUE,
                 turnover=c(1,2,3,4,5,6,7,8),
                 eps=c(1,2,3,4,5,6,7,8),
                 trans.rate=rates.cid4)

repar.bd(hisse.mle,8)

## our logLik methods
logLik.hisse.fit<-function(x,...){
  lik<-x$loglik
  attr(lik,"df")<-(x$AIC+2*lik)/2
  lik
}
## print a table of results
data.frame(
  model=c("CID","BiSSE","HiSSE CID-2",
          "HiSSE CID-4","HiSSE full"),
  logL=sapply(list(cid.mle,bisse.hmle,
                   cid2.mle,cid4.mle,hisse.mle),
              logLik),
  k=sapply(list(cid.mle,bisse.hmle,
                cid2.mle,cid4.mle,hisse.mle),
           function(x) attr(logLik(x),"df")),
  AIC=aic<-sapply(list(cid.mle,bisse.hmle,
                       cid2.mle,cid4.mle,hisse.mle),
                  AIC),
  Akaike.weight=unclass(aic.w(aic))
)


### misse ##########################

## setsampling fraction
rho<-1
## fitMiSSE model
misse2.mle<-MiSSE(micratree,f=rho,
                  turnover=c(1,2),eps=c(1,2))


repar.bd(misse2.mle,2)

## conduct marginal reconstruction
misse2.recon<-MarginReconMiSSE(phy=micratree,
                               f=1,
                               pars=misse2.mle$solution,
                               hidden.states=2)

## graph MiSSE model on thetree
misse2.map<-plot.misse.states(misse2.recon,
                              rate.param="speciation",edge.width=3,
                              show.tip.label=TRUE,type="phylogram",
                              fsize=0.5,legend.position=c(0,0.3,0.85,1),
                              rate.colors=c("yellow","darkblue"))

## fitbirth-death model using MiSSE
misse1.mle<-MiSSE(micratree,
                  f=1,
                  turnover=1,eps=1)

repar.bd(misse1.mle,1)

## copylogLik method for different object
## class: "misse.fit"
logLik.misse.fit<-logLik.hisse.fit

## compile our results
data.frame(
  model=c("MiSSE-1","MiSSE-2"),
  logL=sapply(list(misse1.mle,misse2.mle),
              logLik),
  k=sapply(list(misse1.mle,misse2.mle),
           function(x) attr(logLik(x),"df")),
  AIC=aic<-sapply(list(misse1.mle,misse2.mle),
                  AIC),
  Akaike.weight=unclass(aic.w(aic))
)


#
##
###
####
#####
######
#######
########
#########
##########
###########







##readtreeadn data fromfile
contdata<- read.table("contdata.txt")
nameschara <- read.table("contcahracternames.txt")
colnames(contdata)<- nameschara$V2

micratree <- read.tree("micra.tre")
micratree<-drop.tip(micratree,"M_horrida_122") # M horrida is repeated
micratree$tip.label<-sub("_", ".", micratree$tip.label)
micratree$tip.label<-sub("\\_.*", "", micratree$tip.label)


micraspine <- contdata[-c(1:14),c(1,2,11)]
rownames(micraspine) <- micraspine$Species
names(micraspine)<- c("species","Carapacelength","spinelength")

## check tree and data to ensure matching

check <-name.check(phy = micratree,data = micraspine) 
check
summary(check)

## prune mismatched taxa from thetree
micratree<-drop.tip(micratree, check$tree_not_data)
micraspine <- micraspine[micratree$tip.label,]
check <-name.check(phy = micratree,data = micraspine) 
check

micralenspine<-setNames(micraspine$spinelength ,micraspine$species)
micralencarapace <- setNames(micraspine$Carapacelength ,micraspine$species)

minmicraspine <- setNames(as.numeric(sub("-.*", "", micralenspine)), names(micralenspine))
maxmicraspine <- setNames(as.numeric(sub(".*-", "", micralenspine)), names(micralenspine))
promspine <- (minmicraspine + maxmicraspine)/2

mincarapace <- setNames(as.numeric(sub("-.*", "", micralencarapace)), names(micralencarapace))
maxcarapace <- setNames(as.numeric(sub(".*-", "", micralencarapace)), names(micralencarapace))
promcarapace <- (mincarapace + maxcarapace)/2



spinetocara <- log(promspine/promcarapace)

rangemicraspine <- maxmicraspine - minmicraspine
plot(micratree)

##ln.hosts<-setNames(log(scale_insect.data$host.families), ###Should I log transform my data?
  #                 rownames(scale_insect.data))


## visualize a continuous character map of host plant
## number on scale insect tree
host.map<-contMap(micratree,spinetocara,
                  plot=FALSE)
host.map<-setMap(host.map,c("yellow","darkblue"))
plot(host.map,lwd=c(2,5),outline=FALSE,ftype="off",
     leg.txt="max spine length",legend=60)


rho<-Ntip(micratree)/117

#starting.point.quasse(scale_insect.pruned, ln.hosts)

## run fit.bd and fitContinuous to get starting
## values for our QuaSSE optimization
bd<-fit.bd(micratree,rho=rho)
bm<-fitContinuous(micratree,spinetocara)
p<-setNames(c(bd$b,bd$d,bm$opt$sigsq),
            c("lambda","mu","diffusion"))
p


## define range ofx
xr<-range(spinetocara)+c(-1,1)*20*
  p["diffusion"]
## make linear model for QuaSSE
linear.x<-make.linear.x(xr[1],xr[2])
linear.x

## make QuaSSE likelihood function for variable
## lambda and constrain
lik.lambda<-make.quasse(micratree,
                        spinetocara,lambda=linear.x,mu=constant.x,
                        sampling.f=rho,states.sd=0.1)
lik.lambda<-constrain(lik.lambda,drift~0)

## subsample starting parameter values tomatch
## themodel we’re fitting
pp<-setNames(c(p["lambda"],0,p["mu"],
               p["diffusion"]),argnames(lik.lambda))
pp

## fitour first QuaSSE model
lambda.mle<-find.mle(lik.lambda,x.init=pp,
                     control=list(parscale=0.1),lower=rep(0,4))
coef(lambda.mle)


## make QuaSSE likelihood function forvariable
## mu andconstrain
lik.mu<-make.quasse(micratree,
                    spinetocara,lambda=constant.x,mu=linear.x,
                    sampling.f=rho,states.sd=0.1)
lik.mu<-constrain(lik.mu,drift~0)
## fit variable mu model
pp<-setNames(c(p[c("lambda","mu")],0,
               p["diffusion"]),argnames(lik.mu))
mu.mle<-find.mle(lik.mu,x.init=pp,
                 control=list(parscale=0.1),lower=rep(0,4))

coef(mu.mle)


## create full likelihood function
lik.full<-make.quasse(micratree,
                      spinetocara,lambda=linear.x,mu=linear.x,
                      sampling.f=rho,states.sd=0.1)
lik.full<-constrain(lik.full,drift~0)
## fitfull QuaSSE model
pp<-setNames(c(lambda.mle$par[1:2],mu.mle$par[2:3],
               p["diffusion"]),argnames(lik.full))
pp


full.mle<-find.mle(lik.full,x.init=pp,
                   control=list(parscale=0.1),lower=rep(0,5))
## print model coefficients and log-likelihood
coef(full.mle)

## likelihood function for character
## independent model
lik.cid<-make.quasse(micratree,
                     spinetocara,lambda=constant.x,mu=constant.x,
                     sampling.f=rho,states.sd=0.1)
lik.cid<-constrain(lik.cid,drift~0)
argnames(lik.cid)
## fitCID QuaSSE model and print coefficients
cid.mle<-find.mle(lik.cid,x.init=p,
                  control=list(parscale=0.1),lower=rep(0,3))
coef(cid.mle)

anova(cid.mle,variable.lambda=lambda.mle,
      variable.mu=mu.mle,full.model=full.mle)

















