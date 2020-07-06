# predictions caco
load("covs.ajive.RData")


# iPCA
load("iPCAresults.RData")
source("functions_iPCA.R")
scores_iPCA <- plot_ipca(Sig = iPCAresults$Sig, pcs = 1:2)$ipca_scores
jointscores <- scores_iPCA[ , 1:10]


# ajive
library(pROC)
load("ajiveipca.RData")




# all data
jointscores <- scores_iPCA[ , 1:10]
jointscores.a <- unlist(ajive.ipca$joint_scores)

data.logistic <- as.data.frame(cbind(y = covs$case_ctrl, jointscores, jointscores.a, covs$age.sample, covs$BMI, 
                                     covs$smoking_status_cat))
data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)


xvars10 <- names(data.logistic[c(2:11, 18:ncol(data.logistic))])
formula10 <- paste( 'y', '~', paste( xvars10, collapse=' + ' ) )
xvars5 <- names(data.logistic[c(2:6, 18:ncol(data.logistic))])
formula5 <- paste( 'y', '~', paste( xvars5, collapse=' + ' ) )
xvarsa <- names(data.logistic[c(12:ncol(data.logistic))])
formula.a <- paste( 'y', '~', paste( xvarsa, collapse=' + ' ) )
xvars0 <- names(data.logistic[c(18:ncol(data.logistic))])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )

# in sample predictions

model.fit10 <- glm(formula10 ,
                   family = 'binomial', 
                   data = data.logistic)
model.fit5 <- glm(formula5 ,
                  family = 'binomial', 
                  data = data.logistic)
model.fit.a <- glm(formula.a ,
                   family = 'binomial', 
                   data = data.logistic)
model.fit0 <- glm(formula0 ,
                  family = 'binomial', 
                  data = data.logistic)


data.logistic$pred10 <- predict(model.fit10, data.logistic, type = 'response')
data.logistic$pred5 <- predict(model.fit5, data.logistic, type = 'response')
data.logistic$pred.a <- predict(model.fit.a, data.logistic, type = 'response')
data.logistic$pred0 <- predict(model.fit0, data.logistic, type = 'response')

roc10 <- roc(data.logistic$y, data.logistic$pred10)
roc5 <- roc(data.logistic$y, data.logistic$pred5)
roc.a <- roc(data.logistic$y, data.logistic$pred.a)
roc0 <- roc(data.logistic$y, data.logistic$pred0)


par(oma = c(4, 1, 1, 1))
plot(roc5, col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.8, identity = FALSE, 
     main = 'Case vs control')
plot(roc0, col = 'darkgrey', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.8)
plot(roc10, col = 'darkgreen', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.8)
plot(roc.j, col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.1, print.auc.cex = 0.8)



# cross validation

auc10 <- c()
auc5 <- c()
auc.a <- c()
auc0 <- c()
roc10 <- list()
roc5 <- list()
roc.a <- list()
roc0 <- list()

data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)


for (i in 1:10){
  testIndexes <- which(folds == i, arr.ind = TRUE)
  data.test <- data.logistic[testIndexes, ]
  data.train <- data.logistic[-testIndexes, ]
  model.fit10 <- glm(formula10 ,
                   family = 'binomial', 
                   data = data.train)
  model.fit5 <- glm(formula5 ,
                    family = 'binomial', 
                    data = data.train)
  model.fit.a <- glm(formula.a ,
                     family = 'binomial', 
                     data = data.train)
  model.fit0 <- glm(formula0 ,
                     family = 'binomial', 
                     data = data.train)
  
  
  
  data.test$pred10 <- predict(model.fit10, data.test, type = 'response')
  data.test$pred5 <- predict(model.fit5, data.test, type = 'response')
  data.test$pred.a <- predict(model.fit.a, data.test, type = 'response')
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  
  
  roc10[[i]] <- roc(data.test$y, data.test$pred10)
  roc5[[i]] <- roc(data.test$y, data.test$pred5)
  roc.a[[i]] <- roc(data.test$y, data.test$pred.a)
  roc0[[i]] <- roc(data.test$y, data.test$pred0)
  
  
  
  auc10[i] <- roc(data.test$y, data.test$pred10)$auc
  auc5[i] <- roc(data.test$y, data.test$pred5)$auc
  auc.a[i] <- roc(data.test$y, data.test$pred.a)$auc
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  
  
  
  
}
mean(auc10)
mean(auc5)
mean(auc.a)
mean(auc0)



# predictions metastasis
clin <- readRDS('NOWAC_covariates_database_laboratory_160320.rds')

# merge covs and clin
clin$patient.id <- paste(clin$CaseStatus, clin$Pairs, sep = ' ')
rownames(clin) <- clin$patient.id

covs.all <- merge(covs, clin, sort = FALSE)
# response is  metastasis or not
covs.all$metastasis.gr <- factor(covs.all$Metastasis > '0')





# all data
jointscores <- scores_iPCA[ , 1:10]
jointscores.a <- unlist(ajive.ipca$joint_scores)

data.logistic <- as.data.frame(cbind(y = covs.all$metastasis.gr, jointscores, jointscores.a, covs$age.sample, covs$BMI, 
                                     covs$smoking_status_cat))

data.logistic <- data.logistic[which(!is.na(data.logistic$y)), ]
data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)




xvars10 <- names(data.logistic[c(2:11, 18:ncol(data.logistic))])
formula10 <- paste( 'y', '~', paste( xvars10, collapse=' + ' ) )
xvars5 <- names(data.logistic[c(2:6, 18:ncol(data.logistic))])
formula5 <- paste( 'y', '~', paste( xvars5, collapse=' + ' ) )
xvarsa <- names(data.logistic[c(12:ncol(data.logistic))])
formula.a <- paste( 'y', '~', paste( xvarsa, collapse=' + ' ) )
xvars0 <- names(data.logistic[c(18:ncol(data.logistic))])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )

# in sample predictions

model.fit10 <- glm(formula10 ,
                   family = 'binomial', 
                   data = data.logistic)
model.fit5 <- glm(formula5 ,
                  family = 'binomial', 
                  data = data.logistic)
model.fit.a <- glm(formula.a ,
                   family = 'binomial', 
                   data = data.logistic)
model.fit0 <- glm(formula0 ,
                  family = 'binomial', 
                  data = data.logistic)


data.logistic$pred10 <- predict(model.fit10, data.logistic, type = 'response')
data.logistic$pred5 <- predict(model.fit5, data.logistic, type = 'response')
data.logistic$pred.a <- predict(model.fit.a, data.logistic, type = 'response')
data.logistic$pred0 <- predict(model.fit0, data.logistic, type = 'response')

roc10 <- roc(data.logistic$y, data.logistic$pred10)
roc5 <- roc(data.logistic$y, data.logistic$pred5)
roc.a <- roc(data.logistic$y, data.logistic$pred.a)
roc0 <- roc(data.logistic$y, data.logistic$pred0)


par(oma = c(4, 1, 1, 1))
plot(roc5, col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.8, identity = FALSE, 
     main = 'Metastasis (yes vs no)')
plot(roc0, col = 'darkgrey', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.8)
plot(roc10, col = 'darkgreen', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.8)
plot(roc.a, col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.1, print.auc.cex = 0.8)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",bty = "n", inset=c(0,0.02), xpd = TRUE, ncol = 2,  legend=c("iPCA 5 comp",
                                                                            "Patient Covariates", 
                                                                            "iPCA 10 comp",
                                                                            "aJIVE"),
       lty=c(5,5,5, 5, 3), col = c("black", "darkgrey", "darkgreen", "red", "darkgrey"), cex = 0.8)



# cross validation

auc10 <- c()
auc5 <- c()
auc.a <- c()
auc0 <- c()
roc10 <- list()
roc5 <- list()
roc.a <- list()
roc0 <- list()

data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)


for (i in 1:10){
  testIndexes <- which(folds == i, arr.ind = TRUE)
  data.test <- data.logistic[testIndexes, ]
  data.train <- data.logistic[-testIndexes, ]
  model.fit10 <- glm(formula10 ,
                     family = 'binomial', 
                     data = data.train)
  model.fit5 <- glm(formula5 ,
                    family = 'binomial', 
                    data = data.train)
  model.fit.a <- glm(formula.a ,
                     family = 'binomial', 
                     data = data.train)
  model.fit0 <- glm(formula0 ,
                    family = 'binomial', 
                    data = data.train)
  
  
  
  data.test$pred10 <- predict(model.fit10, data.test, type = 'response')
  data.test$pred5 <- predict(model.fit5, data.test, type = 'response')
  data.test$pred.a <- predict(model.fit.a, data.test, type = 'response')
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  
  
  roc10[[i]] <- roc(data.test$y, data.test$pred10)
  roc5[[i]] <- roc(data.test$y, data.test$pred5)
  roc.a[[i]] <- roc(data.test$y, data.test$pred.a)
  roc0[[i]] <- roc(data.test$y, data.test$pred0)
  
  
  
  auc10[i] <- roc(data.test$y, data.test$pred10)$auc
  auc5[i] <- roc(data.test$y, data.test$pred5)$auc
  auc.a[i] <- roc(data.test$y, data.test$pred.a)$auc
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  
 
  
}
mean(auc10)
mean(auc5)
mean(auc.a)
mean(auc0)









