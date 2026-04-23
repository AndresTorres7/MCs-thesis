library(castor)
library(ape)
library(phytools)



pulldivfit <- function(x) {
  
  pullmodelspertree <- list()
  Ngrid <- 2
  
  recent_frac <- 0.5   # last 30% of time
  recent_pts  <- 1     # points near present
  old_pts     <- Ngrid - recent_pts
  
  ####
  
  ltt <- ltt.plot.coords(x[[1]])
  
  
  min_tips <- 10
  
  ltt[,1] <- ltt[,1]*-1
  
  valid_ages <- ltt[,1][ltt[,2] > min_tips]
  max_age <- max(valid_ages)
  
  
  
  #####
  
  for(i in 1:length(x)) {
    cat("Fitting models for tree number ",i,"\n")

    # edad máxima del árbol
#root_age <- max(branching.times(x[[i]]))
root_age <- max_age



age_grid <- c(
  seq(0, root_age * recent_frac, length.out = recent_pts),
  seq(root_age * recent_frac, root_age, length.out = old_pts)[-1]
)
age_grid[2] <- root_age
age_grid <- round(age_grid,4)
age_grid[2] <- age_grid[2] + 0.0001


#null

fitnull = fit_hbd_pdr_on_grid(x[[i]], 
                          age_grid    = root_age,
                          min_PDR     = -100,
                          max_PDR     = +100,
                          condition   = "crown",
                          Ntrials     = 20,	# perform 10 fitting trials
                          Nthreads    = 2,	# use two CPUs
                          max_model_runtime = 1) 	# limit model evaluation to 1 second




fit = fit_hbd_pdr_on_grid(x[[i]], 
                          oldest_age = max_age,
                          age_grid    = age_grid,
                          min_PDR     = -5,
                          max_PDR     = +5,
                          condition   = "crown",
                          Ntrials     = 20,	# perform 10 fitting trials
                          Nthreads    = 2,	# use two CPUs
                          max_model_runtime = 2,
                          verbose = T,
                          ) 	# limit model evaluation to 1 second   fit_control = list(iter.max = 5000)

models<-list(nulo = fitnull, variable = fit)
pullmodelspertree[[i]] <- models
gc()
  }
 return(pullmodelspertree) 
  
}

micratree <- readNexus("k.BEAST-100-posterior-trees.nex")

pullresults <- pulldivfit(micratree)

saveRDS(pullresults, "pullresultscorrected.rds")

  pullresults <- readRDS("pullresultscorrected.rds")



### Checking what models failed

lista <- c()

for (i in seq_along(pullresults)) {
  
lista <- c(lista, pullresults[[i]]$variable$success  )
  
  
}

onesfailed <- pullresults[!lista]
treesfailed <- micratree[!lista]
fixedpullresults <- pulldivfit(treesfailed)




tablaaic <- data.frame(nuloaic = numeric(), varaic = numeric(), delta = numeric(), pdr1 = numeric(), pdr2 = numeric())

for (i in seq_along(pullresults)) {
  
  tablaaic[i,1] <- pullresults[[i]]$nulo$AIC
  tablaaic[i,2] <- pullresults[[i]]$variable$AIC
  tablaaic[i,3] <- tablaaic[i,2] - tablaaic[i,1]
  tablaaic[i,4] <- pullresults[[i]]$variable$fitted_PDR[1]
  tablaaic[i,5] <- pullresults[[i]]$variable$fitted_PDR[2]
  
}
mean(tablaaic$delta)  


mean(tablaaic$pdr1)
sd(tablaaic$pdr1)
mean(tablaaic$pdr2)
sd(tablaaic$pdr2)

pullresults[[1]]$variable$fitted_PDR
evaluate_spline(pullresults[[1]]$variable$age_grid, pullresults[[1]]$variable$fitted_PDR, splines_degree = 3, Xtarget = 4)
