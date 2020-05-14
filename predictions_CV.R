# prediction and plots
# 10 fold Cross validation

set.seed(052020)

library(pROC)
library(nnet)

load("results/ajiveAssocLastproflik.RData")
load("data/covs.ajive.RData")
load("data/data.ajive.RData")

rj <- ajiveResults$joint_rank
r1 <- ajiveResults$block_decomps[[1]]$individual$rank
r2 <- ajiveResults$block_decomps[[2]]$individual$rank
r3 <- ajiveResults$block_decomps[[3]]$individual$rank

# We only want 5 individual each
r1 <- r2 <- r3 <- 5

jointscores <- unlist(ajiveResults$joint_scores)
indivscores1 <- ajiveResults$block_decomps[[1]]$individual$u[,1:r1]
indivscores2 <- ajiveResults$block_decomps[[2]]$individual$u[,1:r2]
indivscores3 <- ajiveResults$block_decomps[[3]]$individual$u[,1:r3]


# prediction on outcome

# case vs control

data.logistic <- as.data.frame(cbind(y = covs$case_ctrl, jointscores, covs$age.sample, covs$BMI, 
                                     covs$smoking_status_cat,
                                     indivscores1, indivscores2, 
                                     indivscores3))


data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)


xvars <- names(data.logistic[2:(rj+r1+r2+r3+4)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )

xvars0 <- names(data.logistic[(rj+2):(rj+4)])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )

xvars.j <- names(data.logistic[2:(rj+4)])
formula.j <- paste( 'y', '~', paste( xvars.j, collapse=' + ' ) )



# pre set folds
# shuffle data
data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)
# CV
auc <- c()
auc0 <- c()
auc.j <- c()
roc <- list()
roc0 <- list()
roc.j <- list()


for (i in 1:10){
  testIndexes <- which(folds == i, arr.ind = TRUE)
  data.test <- data.logistic[testIndexes, ]
  data.train <- data.logistic[-testIndexes, ]
  model.fit <- glm(formula ,
                   family = 'binomial', 
                   data = data.train)
  model.fit0 <- glm(formula0 ,
                    family = 'binomial', 
                    data = data.train)
  model.fit.j <- glm(formula.j ,
                     family = 'binomial', 
                     data = data.train)
  
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  data.test$pred.j <- predict(model.fit.j, data.test, type = 'response')
  
  roc[[i]] <- roc(data.test$y, data.test$pred)
  roc0[[i]] <- roc(data.test$y, data.test$pred0)
  roc.j[[i]] <- roc(data.test$y, data.test$pred.j)
  
  
  
  auc[i] <- roc(data.test$y, data.test$pred)$auc
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  auc.j[i] <- roc(data.test$y, data.test$pred.j)$auc
  
  lty = c(1, rep(3, 9))
  plot(roc(data.test$y, data.test$pred), col = 'black', add = TRUE, lty = lty[i])
  plot(roc(data.test$y, data.test$pred0), col = 'blue', add = TRUE, lty = lty[i])
  plot(roc(data.test$y, data.test$pred.j), col = 'red', add = TRUE, lty = lty[i])
  
  
}
mean(auc)
mean(auc0)
mean(auc.j)
plot(auc.j-auc0)
abline(h = 0)

which(auc == max(auc))
which(auc0 == max(auc0))
which(auc.j == max(auc.j))


which(auc == min(auc))
which(auc0 == min(auc0))
which(auc.j == min(auc.j))

par(oma = c(4, 1, 1, 1))
plot(roc[[1]], col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.8, identity = FALSE, 
     main = 'Case vs control')
plot(roc0[[1]], col = 'blue', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.8)
plot(roc.j[[1]], col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.8)

plot(roc[[9]], col = 'black', add = TRUE, lty = 3, print.auc = TRUE, print.auc.y = 0.4, 
     print.auc.x = 0.1, print.auc.col = 'darkgrey', print.auc.cex = 0.8)
plot(roc0[[9]], col = 'blue', add = TRUE, lty = 3, print.auc = TRUE, print.auc.y = 0.3, 
     print.auc.x = 0.1, print.auc.col = 'lightblue', print.auc.cex = 0.8)
plot(roc.j[[9]], col = 'red', add = TRUE, lty = 3, print.auc = TRUE, print.auc.y = 0.2, 
    print.auc.x = 0.1, print.auc.col = 'lightcoral', print.auc.cex = 0.8)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",bty = "n", inset=c(0,0.02), xpd = TRUE, ncol = 2,  legend=c("Joint and Individual Components","Joint Components", 
                                                    "Patient Covariates", 
                                                    "(Max AUC)","(Min AUC)"),
       lty=c(5,5,5, 5, 3), col = c("black", "red", "blue", "darkgrey", "darkgrey"), cex = 0.8)


# load clinical covariates
clin <- readRDS('NOWAC_covariates_database_laboratory_160320.rds')

# merge covs and clin
clin$patient.id <- paste(clin$CaseStatus, clin$Pairs, sep = ' ')
rownames(clin) <- clin$patient.id

covs.all <- merge(covs, clin)
# response is  metastasis or not
covs.all$metastasis.gr <- factor(covs.all$Metastasis > '0')


data.logistic <- as.data.frame(cbind(y = covs.all$metastasis.gr,
                                     jointscores, covs.all$age.sample, covs.all$BMI, 
                                     covs.all$smoking_status_cat,
                                     indivscores1, indivscores2, 
                                     indivscores3
))

data.logistic <- data.logistic[which(!is.na(data.logistic$y)), ]
data.logistic$y <- factor(data.logistic$y)



xvars <- names(data.logistic[2:(rj+r1+r2+r3+4)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )

xvars0 <- names(data.logistic[(rj+2):(rj+4)])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )

xvars.j <- names(data.logistic[2:(rj+4)])
formula.j <- paste( 'y', '~', paste( xvars.j, collapse=' + ' ) )



# pre set folds
# shuffle data
data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)
# CV
auc <- c()
auc0 <- c()
auc.j <- c()
roc <- list()
roc0 <- list()
roc.j <- list()


for (i in 1:10){
  testIndexes <- which(folds == i, arr.ind = TRUE)
  data.test <- data.logistic[testIndexes, ]
  data.train <- data.logistic[-testIndexes, ]
  model.fit <- glm(formula ,
                   family = 'binomial', 
                   data = data.train)
  model.fit0 <- glm(formula0 ,
                    family = 'binomial', 
                    data = data.train)
  model.fit.j <- glm(formula.j ,
                     family = 'binomial', 
                     data = data.train)
  
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  data.test$pred.j <- predict(model.fit.j, data.test, type = 'response')
  
  roc[[i]] <- roc(data.test$y, data.test$pred)
  roc0[[i]] <- roc(data.test$y, data.test$pred0)
  roc.j[[i]] <- roc(data.test$y, data.test$pred.j)
  
  
  
  auc[i] <- roc(data.test$y, data.test$pred)$auc
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  auc.j[i] <- roc(data.test$y, data.test$pred.j)$auc
  
  lty = c(1, rep(3, 9))
  plot(roc(data.test$y, data.test$pred), col = 'black', add = TRUE, lty = lty[i])
  plot(roc(data.test$y, data.test$pred0), col = 'blue', add = TRUE, lty = lty[i])
  plot(roc(data.test$y, data.test$pred.j), col = 'red', add = TRUE, lty = lty[i])
  
  
}
mean(auc)
mean(auc0)
mean(auc.j)
plot(auc.j-auc0)
abline(h = 0)

which(auc == max(auc))
which(auc0 == max(auc0))
which(auc.j == max(auc.j))


which(auc == min(auc))
which(auc0 == min(auc0))
which(auc.j == min(auc.j))

par(oma = c(4, 1, 1, 1))
plot(roc[[3]], col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.x = 0.4, print.auc.cex = 0.8, identity = FALSE, 
     main = 'Metastasis (yes vs no)')
plot(roc0[[3]], col = 'blue', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.x = 0.4,print.auc.cex = 0.8)
plot(roc.j[[3]], col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.x = 0.4, print.auc.cex = 0.8)

plot(roc[[2]], col = 'black', add = TRUE, lty = 3, print.auc = TRUE, print.auc.y = 0.4, 
     print.auc.x = 0.1, print.auc.col = 'darkgrey', print.auc.cex = 0.8)
plot(roc0[[2]], col = 'blue', add = TRUE, lty = 3, print.auc = TRUE, print.auc.y = 0.3, 
     print.auc.x = 0.1, print.auc.col = 'lightblue', print.auc.cex = 0.8)
plot(roc.j[[2]], col = 'red', add = TRUE, lty = 3, print.auc = TRUE, print.auc.y = 0.2, 
     print.auc.x = 0.1, print.auc.col = 'lightcoral', print.auc.cex = 0.8)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",bty = "n", inset=c(0,0.02), xpd = TRUE, ncol = 2,  legend=c("Joint and Individual Components","Joint Components", 
                                                                            "Patient Covariates", 
                                                                            "(Max AUC)","(Min AUC)"),
       lty=c(5,5,5, 5, 3), col = c("black", "red", "blue", "darkgrey", "darkgrey"), cex = 0.8)





