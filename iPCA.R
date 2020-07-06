# load libraries 
# to install iPCA you need to download the folder
# or load the functions
# library(iPCA)
source('functions_iPCA.R')
source('ajive.dataprep.R')
library(SpatioTemporal)
library(ggplot2)
library(sparseEigen)

############################################
### iPCA
# implement integrated PCA
# set initial lambdas
lams <- c(1e-4, 1e-2, 1, 100, 10000)

# it needs the transpose wrt jive 
# data.ipca <- lapply(dataNOWAC, function(x) t(x))
# SVD miss for missing data
# center and scale 
# same as ajive 
data.ipca <- ajive.dataprep(dataNOWAC)
# select lambdas
best_lambdas <- choose_lambdas(dat = data.ipca, q = "multfrob", 
                               lams = lams, greedy.search = T)$best_lambdas
best_lambdas
# run Flip FLop Algorithm 
iPCAresults <- FFmleMultFrob(dat = data.ipca, 
                             lamDs = best_lambdas)


# visualize results 
covs$y <- 'control'
covs$y[covs$case_ctrl == 'case'] <- 'case' 
covs$y <- as.factor(covs$y)
plot_ipca(Sig = iPCAresults$Sig, pcs = 1:2)$ipca_plot
plot_ipca(Sig = iPCAresults$Sig, y = covs$y, pcs = 1:2, show_legend = TRUE)

# visualize scores
scores_iPCA <- plot_ipca(Sig = iPCAresults$Sig, pcs = 1:2)$ipca_scores
# variation proportions as in jive
plot_ipca_varexplained(Xs = as.matrix(data.ipca), 
                       Sig = iPCAresults$Sig, 
                       Delts = iPCAresults$Delts)



# compare to ajive

ajive.ipca <- ajive(data.ipca, initial_signal_ranks = c(23, 18, 9))






# run iPCA on pairwise 
data12 <- list(dataNOWAC[[1]], dataNOWAC[[2]])
data.ipca12 <- ajive.dataprep(data12)
# select lambdas
best_lambdas <- choose_lambdas(dat = data.ipca12, q = "multfrob", 
                               lams = lams, greedy.search = T)$best_lambdas
best_lambdas
# run Flip FLop Algorithm 
iPCAresults12 <- FFmleMultFrob(dat = data.ipca12, 
                             lamDs = best_lambdas)


data13 <- list(dataNOWAC[[1]], dataNOWAC[[3]])
data.ipca13 <- ajive.dataprep(data13)
# select lambdas
best_lambdas <- choose_lambdas(dat = data.ipca13, q = "multfrob", 
                               lams = lams, greedy.search = T)$best_lambdas
best_lambdas
# run Flip FLop Algorithm 
iPCAresults13 <- FFmleMultFrob(dat = data.ipca13, 
                               lamDs = best_lambdas)


data23 <- list(dataNOWAC[[2]], dataNOWAC[[3]])
data.ipca23 <- ajive.dataprep(data23)
# select lambdas
best_lambdas <- choose_lambdas(dat = data.ipca23, q = "multfrob", 
                               lams = lams, greedy.search = T)$best_lambdas
best_lambdas
# run Flip FLop Algorithm 
iPCAresults23 <- FFmleMultFrob(dat = data.ipca23, 
                               lamDs = best_lambdas)

