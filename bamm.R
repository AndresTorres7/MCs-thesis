library(BAMMtools) # Assuming you have installed BAMMtools!

setBAMMpriors(read.tree("micra.tre"))
tree <- read.tree("micra.tre")

setwd("millionResults")

edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.1)

mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
post_probs['X'] / post_probs['Y']


shift_probs <- summary(edata)

postfile <- "mcmc_out.txt"
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)
computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)


mysample <- 25
nrow(edata$eventData[[ mysample ]])

shiftnodes <- getShiftNodesFromIndex(edata, index = mysample)

plot.phylo(tree, show.tip.label = F)
nodelabels(node = shiftnodes, pch=21, col="red", cex=1.5)
