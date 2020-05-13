# load libraries 
# to install iPCA you need to download the folder
# or load the functions

source('functions_iPCA.R')
source('ajive.dataprep.R')
library(SpatioTemporal)
library(ggplot2)
library(sparseEigen)
library(pROC)

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




# predictions caco
load("data/covs.ajive.RData")


# iPCA
scores_iPCA <- plot_ipca(Sig = iPCAresults$Sig, pcs = 1:2)$ipca_scores
jointscores <- scores_iPCA[ , 1:5]


data.logistic <- as.data.frame(cbind(y = covs$case_ctrl, jointscores, covs$age.sample, covs$BMI, 
                                     covs$smoking_status_cat))


data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)


xvars <- names(data.logistic[2:ncol(data.logistic)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ))

xvars0 <- names(data.logistic[12:14])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )


model.fit <- glm(formula  , family = 'binomial', 
                 data = data.logistic)

model.fit0 <- glm(formula0 ,
                  family = 'binomial', 
                  data = data.logistic)


data.logistic$pred <- predict(model.fit, data.logistic, type = 'response')
roc <- roc(data.logistic$y, data.logistic$pred)

data.logistic$pred0 <- predict(model.fit0, data.logistic, type = 'response')
roc0 <- roc(data.logistic$y, data.logistic$pred0)



# ajive
jointscores <- unlist(ajive.ipca$joint_scores)

# prediction on outcome

# case vs control

data.logistic.j <- as.data.frame(cbind(y = covs$case_ctrl, jointscores, covs$age.sample, covs$BMI, 
                                     covs$smoking_status_cat))


data.logistic.j$y[data.logistic.j$y == '2'] <- '0'
data.logistic.j$y <- as.factor(data.logistic.j$y)


xvars <- names(data.logistic.j[2:ncol(data.logistic.j)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )


model.fit.j <- glm(formula ,
                  family = 'binomial', 
                  data = data.logistic.j)


data.logistic.j$pred <- predict(model.fit.j, data.logistic.j, type = 'response')
roc.j <- roc(data.logistic.j$y, data.logistic.j$pred)


par(oma = c(4, 1, 1, 1))
plot(roc, col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.8, identity = FALSE, 
     main = 'Case vs control')
plot(roc0, col = 'darkgrey', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.8)
plot(roc10, col = 'darkgreen', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.8)
plot(roc.j, col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.1, print.auc.cex = 0.8)



# predictions metastasis
clin <- readRDS('NOWAC_covariates_database_laboratory_160320.rds')

# merge covs and clin
clin$patient.id <- paste(clin$CaseStatus, clin$Pairs, sep = ' ')
rownames(clin) <- clin$patient.id

covs.all <- merge(covs, clin)
# response is  metastasis or not
covs.all$metastasis.gr <- factor(covs.all$Metastasis > '0')



# iPCA
scores_iPCA <- plot_ipca(Sig = iPCAresults$Sig, pcs = 1:2)$ipca_scores
jointscores <- scores_iPCA[ , 1:5]
jointscores10 <- scores_iPCA[, 1:10]


data.logistic <- as.data.frame(cbind(y = covs.all$metastasis.gr, jointscores, covs.all$age.sample, covs.all$BMI, 
                                     covs.all$smoking_status_cat))


data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)


xvars <- names(data.logistic[2:ncol(data.logistic)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ))

xvars0 <- names(data.logistic[(ncol(data.logistic)-2):ncol(data.logistic)])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )


model.fit <- glm(formula  , family = 'binomial', 
                 data = data.logistic)

model.fit0 <- glm(formula0 ,
                  family = 'binomial', 
                  data = data.logistic)


data.logistic$pred <- predict(model.fit, data.logistic, type = 'response')
roc10 <- roc(data.logistic$y, data.logistic$pred)

data.logistic$pred0 <- predict(model.fit0, data.logistic, type = 'response')
roc0 <- roc(data.logistic$y, data.logistic$pred0)



# ajive
jointscores <- unlist(ajive.ipca$joint_scores)

# prediction on outcome

# case vs control
data.logistic.j <- as.data.frame(cbind(y = covs.all$metastasis.gr, jointscores, covs.all$age.sample, covs.all$BMI, 
                                     covs.all$smoking_status_cat))




data.logistic.j$y[data.logistic.j$y == '2'] <- '0'
data.logistic.j$y <- as.factor(data.logistic.j$y)


xvars <- names(data.logistic.j[2:ncol(data.logistic.j)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )


model.fit.j <- glm(formula ,
                   family = 'binomial', 
                   data = data.logistic.j)


data.logistic.j$pred <- predict(model.fit.j, data.logistic.j, type = 'response')
roc.j <- roc(data.logistic.j$y, data.logistic.j$pred)


par(oma = c(4, 1, 1, 1))
plot(roc, col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.8, identity = FALSE, 
     main = 'Metastasis (yes vs no)')
plot(roc0, col = 'darkgrey', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.8)
plot(roc10, col = 'darkgreen', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.8)
plot(roc.j, col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.1, print.auc.cex = 0.8)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",bty = "n", inset=c(0,0.02), xpd = TRUE, ncol = 2,  legend=c("iPCA 5 comp",
                                                                            "Patient Covariates", 
                                                                            "iPCA 10 comp",
                                                                            "aJIVE"),
       lty=c(5,5,5, 5, 3), col = c("black", "darkgrey", "darkgreen", "red", "darkgrey"), cex = 0.8)







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

