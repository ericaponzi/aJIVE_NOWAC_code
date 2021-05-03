####################
# Angle based JIVE #
####################

# load libraries
# devtools::install_github("idc9/r_jive")
library(ajive)
library(SpatioTemporal)
library(ggplot2)

source('functions/ajive.dataprep.R')
data.ajive <- ajive.dataprep(dataNOWAC)
colnames(data.ajive[[1]]) <- rownames(dataNOWAC[[1]])
colnames(data.ajive[[2]]) <- rownames(dataNOWAC[[2]])
colnames(data.ajive[[3]]) <- rownames(dataNOWAC[[3]])

# need to determine initial ranks
# look at screeplots
source('functions/Screeplots.R')
screeplot(data.ajive)

# use profile likelihood
# SVD on each source
singular.values1 <- svd(data.ajive[[1]])[['d']]
singular.values2 <- svd(data.ajive[[2]])[['d']]
singular.values3 <- svd(data.ajive[[3]])[['d']]
# check first if we want to eliminate first 1 or 2 
# if too extremes
singular.values1 <- singular.values1[-c(1,2)]
singular.values2 <- singular.values2[-c(1)]
singular.values3 <- singular.values3[-1]

singular.val <- list(singular.values1, 
                        singular.values2, 
                        singular.values3)
r <- c()
for (s in 1:3){
  singular.values <- singular.val[[s]]
# calculate profile loglikelihood utility function
l <- length(singular.values)
proflik <- c()
for (i in 1: l){
  mu1 <- mean(singular.values[1:i])
  s1 <- sum((singular.values[1:i]- mu1)^2)
  mu2 <- mean(singular.values[(i+1):l])
  s2 <- sum((singular.values[(i+1):l]-mu2)^2)
  if (i == l) s2 <- 0
  proflik[i] <- s1+s2
}
# visualize results
plot(-proflik)
# rank will be arg min 
# +1 because we left one out
r[s] <- which.min(proflik)
}

# add eliminated ones
r <- r + c(2,2,1)

# run ajive
ajiveResults <- ajive(data.ajive, 
                      initial_signal_ranks = r) 

load('ajiveAssocLastproflik.RData')
load('data.ajive.RData')
decomposition_heatmaps(data.ajive, ajiveResults)

# pairwise ajive comparisons
data12 <- list(t(methylation), t(mRNA))
data.ajive12 <- ajive.dataprep(data12)
ajive12 <- ajive(data.ajive12, 
                 initial_signal_ranks = r[c(1,2)])

data23 <- list(t(mRNA), t(miRNA))
data.ajive23 <- ajive.dataprep(data23)
ajive23 <- ajive(data.ajive23, 
                 initial_signal_ranks = r[c(2,3)])

data13 <- list(t(methylation), t(miRNA))
data.ajive13 <- ajive.dataprep(data13)
ajive13 <- ajive(data.ajive13, 
                 initial_signal_ranks = r[c(1,3)])
