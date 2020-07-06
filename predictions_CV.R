# prediction and plots
# 10 fold Cross validation

set.seed(052020)

library(pROC)
library(nnet)

load("ajiveAssocLastproflik.RData")
load("covs.ajive.RData")
load("data.ajive.RData")

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

xvars.i <- names(data.logistic[(rj+2):((rj+r1+r2+r3+4))])
formula.i <- paste( 'y', '~', paste( xvars.i, collapse=' + ' ) )

# pre set folds
# shuffle data
data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)


# CV
auc <- c()
auc0 <- c()
auc.j <- c()
auc.i <- c()
roc <- list()
roc0 <- list()
roc.j <- list()
roc.i <- list()


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
  model.fit.i <- glm(formula.i ,
                     family = 'binomial', 
                     data = data.train)
  
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  data.test$pred.j <- predict(model.fit.j, data.test, type = 'response')
  data.test$pred.i <- predict(model.fit.i, data.test, type = 'response')
  
  roc[[i]] <- roc(data.test$y, data.test$pred)
  roc0[[i]] <- roc(data.test$y, data.test$pred0)
  roc.j[[i]] <- roc(data.test$y, data.test$pred.j)
  roc.i[[i]] <- roc(data.test$y, data.test$pred.i)
  
  
  
  auc[i] <- roc(data.test$y, data.test$pred)$auc
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  auc.j[i] <- roc(data.test$y, data.test$pred.j)$auc
  auc.i[i] <- roc(data.test$y, data.test$pred.i)$auc
  
 
  
}
mean(auc)
mean(auc0)
mean(auc.j)
mean(auc.i)



# load clinical covariates
clin <- readRDS('NOWAC_covariates_database_laboratory_160320.rds')

# merge covs and clin
clin$patient.id <- paste(clin$CaseStatus, clin$Pairs, sep = ' ')
rownames(clin) <- clin$patient.id

covs.all <- merge(covs, clin, sort = FALSE)
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

xvars.i <- names(data.logistic[(rj+5):((rj+r1+r2+r3+4))])
formula.i <- paste( 'y', '~', paste( xvars.i, collapse=' + ' ) )


# pre set folds
# shuffle data
data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)

# CV
auc <- c()
auc0 <- c()
auc.j <- c()
auc.i <- c()
roc <- list()
roc0 <- list()
roc.j <- list()
roc.i <- list()


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
  model.fit.i <- glm(formula.i ,
                     family = 'binomial', 
                     data = data.train)
  
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  data.test$pred.j <- predict(model.fit.j, data.test, type = 'response')
  data.test$pred.i <- predict(model.fit.i, data.test, type = 'response')
  
  roc[[i]] <- roc(data.test$y, data.test$pred)
  roc0[[i]] <- roc(data.test$y, data.test$pred0)
  roc.j[[i]] <- roc(data.test$y, data.test$pred.j)
  roc.i[[i]] <- roc(data.test$y, data.test$pred.i)
  
  
  auc[i] <- roc(data.test$y, data.test$pred)$auc
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  auc.j[i] <- roc(data.test$y, data.test$pred.j)$auc
  auc.i[i] <- roc(data.test$y, data.test$pred.i)$auc
  
  
  
}
mean(auc)
mean(auc0)
mean(auc.j)
mean(auc.i)





