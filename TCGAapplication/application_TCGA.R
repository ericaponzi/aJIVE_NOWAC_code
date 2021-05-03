

# application on TCGA data (from mixomics/rJIVE R package)



showVarExplained.ajive <- function(ajiveResults, blocks,
                                   col = c("grey20", "grey43", "grey65")){
  l <- length(blocks)
  # joint variance
  VarJoint = rep(0, l)
  for (i in 1:l) VarJoint[i] = norm(as.matrix(ajiveResults$block_decomps[[i]]$joint[[1]]), 
                                    type = "F")^2/norm(blocks[[i]], type = "F")^2
  
  # individual variances
  VarIndiv = rep(0, l)
  for (i in 1:l) VarIndiv[i] = norm(as.matrix(ajiveResults$block_decomps[[i]]$individual[[1]]), 
                                    type = "F")^2/norm(blocks[[i]], type = "F")^2
  
  # residual variance
  VarSubtr = 1 - VarJoint - VarIndiv
  # plot
  
  datavar <- data.frame(var = c(VarJoint, VarIndiv, VarSubtr), 
                        source = rep(c('mRNA', 'miRNA', 'Methylation'), 3), 
                        type = c(rep('Joint',3), rep('Individual', 3), rep('Residual', 3)))
  p <- ggplot(datavar, aes(x = source, y = var, fill = type)) +
    geom_bar(stat = 'identity', position = position_dodge())+
    labs(title="Variation Explained")+ scale_fill_brewer(palette="Blues")+ xlab('')
  
  p
}




library(mixOmics)
library(ajive)
library(SpatioTemporal)
library(ggplot2)
library(pROC)
library(nnet)
load(file = 'TCGA.normalised.mixDIABLO')



# extract training and test data
data = list(mRNA = data.train$mrna, 
           miRNA = data.train$mirna,
           methylation = data.train$methylation)



# check dimension
lapply(data, dim)

# outcome
Y = c(as.character(data.train$subtype))
Y <- as.factor(Y)
summary(Y)


# prepare data for ajive
# modify this function (do not need transpose here!)

data2 <- lapply(data, function(l) t(l))

ajive.dataprep <- function(data){
  
  
  data.ajive <- list()
  for (l in 1:length(data)){
    X <- as.data.frame(data[[l]])
    
    # svdmiss
    Xtemp <- SVDmiss(X)[[1]]
    Ximp0 <- Xtemp$u %*% diag(x = Xtemp$d) %*% t(Xtemp$v)
    
    
    # center values
    centerValues <- apply(Ximp0, 1, mean, na.rm = T)
    Ximp <- Ximp0 - matrix(rep(centerValues, 
                               ncol(Ximp0)), nrow = nrow(Ximp0))
    data.ajive[[l]] <- t(Ximp)
    
    
  }
  
  return(data.ajive)
}

data.ajive  <- ajive.dataprep(data2)
data.ajive2 <- data.ajive


# determine ranks via profile likelihood
singular.values1 <- svd(data.ajive2[[1]])[['d']]
singular.values2 <- svd(data.ajive2[[2]])[['d']]
singular.values3 <- svd(data.ajive2[[3]])[['d']]

# remove the first cause it needs a category for itself
# 
singular.values1 <- singular.values1[-c(1)]
singular.values2 <- singular.values2[-c(1)]
#singular.values3 <- singular.values3[-c(1,2)]

singular.val <- list(singular.values1, 
                     singular.values2, 
                     singular.values3)

# calculate profile loglikelihood utility function
r <- c()
for (s in 1:3){
  singular.values <- singular.val[[s]]
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
r <- r + c(1,1, 0)


# run ajive
ajiveResults <- ajive(data.ajive2, initial_signal_ranks = r)
ajiveResults$joint_rank
save(ajiveResults, file = 'C:/Users/ericapo/Desktop/NOWAC/0912/ajive.RData')

load('C:/Users/ericapo/Desktop/NOWAC/Review/ajiveResultsTCGA.RData')

source('C:/Users/ericapo/Desktop/NOWAC/ResultsFromHunt/data_heatmap_modified.R')
decomposition_heatmaps(data.ajive2, ajiveResults)

showVarExplained.ajive(ajiveResults, data.ajive2)



## use joint and individual ajive components to predict Y
rj <- ajiveResults$joint_rank
r1 <- ajiveResults$block_decomps[[1]]$individual$rank
r2 <- ajiveResults$block_decomps[[2]]$individual$rank
r3 <- ajiveResults$block_decomps[[3]]$individual$rank


# pick first five
r1 <- r2 <-  5
jointscores <- unlist(ajiveResults$joint_scores)
indivscores1 <- ajiveResults$block_decomps[[1]]$individual$u[,1:r1]
indivscores2 <- ajiveResults$block_decomps[[2]]$individual$u[,1:r2]
indivscores3 <- ajiveResults$block_decomps[[3]]$individual$u[,1:r3]

# non integrated analysis
mRNA <- data.ajive2[[1]]
miRNA <- data.ajive2[[2]]
met <- data.ajive2[[3]]

dim(mRNA)
dim(miRNA)
dim(met)



# pca on mRNA
mrna.pca <- prcomp(mRNA, center = TRUE, scale = TRUE)

std_dev <- mrna.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex[1:20], xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

mrna.pcs <- mrna.pca$x[,1:5] #try with 7

# pca on miRNA
mirna.pca <- prcomp(miRNA, center = TRUE, scale = TRUE)

std_dev <- mirna.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex[1:20], xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

mirna.pcs <- mirna.pca$x[,1:5] #try with 8

# pca on mRNA
met.pca <- prcomp(met, center = TRUE, scale = TRUE)

std_dev <- met.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex[1:20], xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

met.pcs <- met.pca$x[,1:2] #try with 7


# prediction on outcome
#Y <- clusts
data.logistic <- as.data.frame(cbind(Y,  jointscores, 
                                     indivscores1, indivscores2, indivscores3,
                                     mrna.pcs, 
                                     mirna.pcs, 
                                     met.pcs
                                     ))

names(data.logistic)[(rj+r1+r2+r3+2):ncol(data.logistic)] <- c('PC1mr', 'PC2mr', 'PC3mr','PC4mr', 'PC5mr',
                                                               'PC1mi', 'PC2mi', 'PC3mi', 'PC4mi', 'PC5mi',
                                                            'PC1met', 'PC2met')


# test LumA vs not LumA
#data.logistic[data.logistic$Y <2, ]$Y <- '0'
#data.logistic$Y[data.logistic$Y >1 ] <- '1'
data.logistic$Y <- as.factor(data.logistic$Y)


xvars <- names(data.logistic[2:(rj+r1+r2+r3+1)])
formula <- paste( 'Y', '~', paste( xvars, collapse=' + ' ) )

xvars.j <- names(data.logistic[2:(rj+1)])
formula.j <- paste( 'Y', '~', paste( xvars.j, collapse=' + ' ) )

xvars.i <- names(data.logistic[(rj+2):((rj+r1+r2+r3+1))])
formula.i <- paste( 'Y', '~', paste( xvars.i, collapse=' + ' ) )

xvars.ni <- names(data.logistic[(rj+r1+r2+r3+2):ncol(data.logistic)])
formula.ni <- paste( 'Y', '~', paste( xvars.ni, collapse=' + ' ) )

xvars.ni.m <- names(data.logistic[(rj+r1+r2+r3+2): (rj+r1+r2+r3+6)])
formula.ni.m <- paste( 'Y', '~', paste( xvars.ni.m, collapse=' + ' ) )

#index <- sample(150, 100)
data.train <- data.logistic
model.fit <- multinom(formula , 
                 data = data.train)
model.fit.j <- multinom(formula.j, 
                   data = data.train)
model.fit.i <- multinom(formula.i, 
                   data = data.train)
model.fit.ni <- multinom(formula.ni, 
                    data = data.train)
model.fit.ni.m <- multinom(formula.ni.m, 
                         data = data.train)


data.test <- data.logistic
data.test$pred <- predict(model.fit, data.test, type = 'class')
data.test$pred.j <- predict(model.fit.j, data.test, type = 'probs')
data.test$pred.i <- predict(model.fit.i, data.test, type = 'probs')
data.test$pred.ni <- predict(model.fit.ni, data.test, type = 'class')
data.test$pred.ni.m <- predict(model.fit.ni.m, data.test, type = 'class')

multiclass.roc(data.test$Y, as.numeric(data.test$pred))$auc
multiclass.roc(data.test$Y, as.numeric(data.test$pred.j))$auc
multiclass.roc(data.test$Y, as.numeric(data.test$pred.i))$auc
multiclass.roc(data.test$Y, as.numeric(data.test$pred.ni))$auc




# CV

# pre set folds
# shuffle data
set.seed(202011)
data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)


auc <- c()
auc.j <- c()
auc.i <- c()
auc.ni <- c()
auc.ni.m <- c()
roc <- list()
roc.ni <- list()
roc.j <- list()
roc.i <- list()
roc.ni.m <- list()

for (i in 1:10){
  testIndexes <- which(folds == i, arr.ind = TRUE)
  data.test <- data.logistic[testIndexes, ]
  data.train <- data.logistic[-testIndexes, ]
  model.fit <- multinom(formula, 
                   data = data.train)
  model.fit.ni <- multinom(formula.ni, 
                    data = data.train)
  model.fit.j <- multinom(formula.j, 
                     data = data.train)
  model.fit.i <- multinom(formula.i, 
                     data = data.train)
  model.fit.ni.m <- multinom(formula.ni.m, 
                          data = data.train)
  
  
  data.test$pred <- predict(model.fit, data.test, type = 'class')
  data.test$pred.ni <- predict(model.fit.ni, data.test, type = 'class')
  data.test$pred.j <- predict(model.fit.j, data.test, type = 'class')
  data.test$pred.i <- predict(model.fit.i, data.test, type = 'class')
  data.test$pred.ni.m <- predict(model.fit.ni.m, data.test, type = 'class')
  
  roc[[i]] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred))
  roc.ni[[i]] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.ni))
  roc.j[[i]] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.j))
  roc.i[[i]] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.i))
  roc.ni.m[[i]] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.ni.m))
  
  
  auc[i] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred))$auc
  auc.ni[i] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.ni))$auc
  auc.ni.m[i] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.ni.m))$auc
  auc.j[i] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.j))$auc
  auc.i[i] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.i))$auc
 
  
}
mean(auc)
mean(auc.ni)
mean(auc.ni.m)
mean(auc.j)
mean(auc.i)

library(randomForest)
data.rf <- data.logistic[, 1:(rj+r1+r2+r3+1)]
ind = sample(2, nrow(data.rf), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.rf[ind == 1, ]
test.data = data.rf[ind == 2, ]

rf <- randomForest(Y ~., data = train.data, 
                   ntree = 1000, proximity = TRUE, 
                   na.action = na.omit)
predictions <- predict(rf, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf), train.data$Y)
cm <- table(predictions, test.data$Y)
accuracy = sum(diag(cm))/sum(cm)
print(rf)
varImpPlot(rf, sort = TRUE, n.var = 10)



imp <- importance(rf)
ord <- order(imp, decreasing = TRUE)[1:10]
# 1 = caco
# 2, 3, 4, 5 = joint
# 6, 7, 8, 9, 10 = Ind mRNA
# 11, 12, 13, 14, 15 = Ind miRNA
# 16 = Ind met

shape <- rep('darkgrey', 10)
shape[c(1, 2, 3, 4)] <- 'blue'

varimp <- data.frame(imp = imp[ord], var = c('Joint2', 'Joint1', 'Joint3', 'Joint4', 
                                             'mRNAInd4', 'mRNAInd1', 'MetInd1', 
                                              'miRNAInd1', 'miRNAInd3',
                                             'mRNAInd5'), shape = shape)

varimp <- varimp[order(varimp$imp), ]



qplot(x = imp, 
      y = reorder(var, imp), color = factor(shape), data = varimp, geom = "point") +  geom_point(size = 2) +
  xlab("MeanDecreaseGini") + ylab('') + theme(legend.position = "none") +
  ggtitle('Variable Importance') + geom_segment(aes(x = 0, y = var, xend = imp, yend = var ), lty = 2)


# save all results 
save(ajiveResults, file = 'ajiveResultsTCGA.RData')
l <- list(data.logistic, folds)
save(l, file = 'CV_TCGA.RData')
save(rf, file = 'RF_TCGA.RData')



data.rf <- data.logistic
ind = sample(2, nrow(data.rf), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.rf[ind == 1, ]
test.data = data.rf[ind == 2, ]

rf <- randomForest(Y ~., data = train.data[, c('Y', xvars)], 
                   ntree = 1000, proximity = TRUE, 
                   na.action = na.omit)
predictions <- predict(rf, newdata = test.data[, c('Y', xvars)])
train.data <- na.omit(train.data)
table(predict(rf), train.data$Y)
cm <- table(predictions, test.data$Y)
accuracy = sum(diag(cm))/sum(cm)
print(rf)
varImpPlot(rf, sort = TRUE, n.var = 10)

multiclass.roc(test.data$Y, as.numeric(predictions))


rf.i <- randomForest(Y~., data = train.data[, c('Y', xvars.i)], 
                     ntree = 1000, proximity = TRUE, 
                     na.action = na.omit)
predictions.i <- predict(rf.i, newdata = test.data[, c('Y', xvars.i)])
train.data <- na.omit(train.data)
table(predict(rf.i), train.data$Y)
cm.i <- table(predictions.i, test.data$Y)
accuracy.i = sum(diag(cm.i))/sum(cm.i)
multiclass.roc(test.data$Y, as.numeric(predictions.i))
print(rf.i)

rf.j <- randomForest(Y~., data = train.data[, c('Y', xvars.j)], 
                     ntree = 1000, proximity = TRUE, 
                     na.action = na.omit)
predictions.j <- predict(rf.j, newdata = test.data[, c('Y', xvars.j)])
train.data <- na.omit(train.data)
table(predict(rf.j), train.data$Y)
cm.j <- table(predictions.j, test.data$Y)
accuracy.j = sum(diag(cm.j))/sum(cm.j)
multiclass.roc(test.data$Y, as.numeric(predictions.j))
print(rf.j)

rf.ni <- randomForest(Y~., data = train.data[, c('Y', xvars.ni)], 
                     ntree = 1000, proximity = TRUE, 
                     na.action = na.omit)

predictions.ni <- predict(rf.ni, newdata = test.data[, c('Y', xvars.ni)])
train.data <- na.omit(train.data)
table(predict(rf.ni), train.data$Y)
cm.ni <- table(predictions.ni, test.data$Y)
accuracy.ni = sum(diag(cm.ni))/sum(cm.ni)
multiclass.roc(test.data$Y, as.numeric(predictions.ni))
print(rf.ni)

rf.list <- list(rf, rf.j, rf.i, rf.ni)
names(rf.list) <- c('rf', 'rf.j', 'rf.i', 'rf.ni')
save(rf.list, file = 'RFlistTCGA.RData')

# check impact of mRNA on non-integrative models
data.ni <- as.data.frame(cbind(Y, mrna.pcs, mirna.pcs, met.pcs))
names(data.ni) <- c('Y', 'PC1mr', 'PC2mr', 'PC3mr', 'PC4mr', 'PC5mr',
                    'PC1mi', 'PC2mi', 'PC3mi', 'PC4mi', 'PC5mi',
                    'PC1met', 'PC2met')

xvars.ni <- names(data.ni[2:ncol(data.ni)])
formula.ni <- paste( 'Y', '~', paste( xvars.ni, collapse=' + ' ) )

model.fit.ni <- multinom(formula.ni, 
                         data = data.ni)

summary(model.fit.ni)

xvars.ni.m <- names(data.ni[2:6])
formula.ni.m <- paste( 'Y', '~', paste( xvars.ni.m, collapse=' + ' ) )

model.fit.ni.m <- multinom(formula.ni.m, 
                         data = data.ni)

summary(model.fit.ni.m)

xvars.ni.mi <- names(data.ni[7:11])
formula.ni.mi <- paste( 'Y', '~', paste( xvars.ni.mi, collapse=' + ' ) )

model.fit.ni.mi <- multinom(formula.ni.mi, 
                           data = data.ni)

summary(model.fit.ni.mi)

xvars.ni.me <- names(data.ni[12:13])
formula.ni.me <- paste( 'Y', '~', paste( xvars.ni.me, collapse=' + ' ) )

model.fit.ni.me <- multinom(formula.ni.me, 
                            data = data.ni)

summary(model.fit.ni.me)


data.ni$pred.ni <- predict(model.fit.ni, data.ni, type = 'prob')
multiclass.roc(data.ni$Y, as.numeric(data.ni$pred.ni))$auc
data.ni$pred.ni.m <- predict(model.fit.ni.m, data.ni, type = 'prob')
multiclass.roc(data.ni$Y, as.numeric(data.ni$pred.ni.m))$auc
data.ni$pred.ni.mi <- predict(model.fit.ni.mi, data.ni, type = 'prob')
multiclass.roc(data.ni$Y, as.numeric(data.ni$pred.ni.mi))$auc
data.ni$pred.ni.me <- predict(model.fit.ni.me, data.ni, type = 'prob')
multiclass.roc(data.ni$Y, as.numeric(data.ni$pred.ni.me))$auc

# lasso 
X.full <- cbind(met, mRNA, miRNA)
p <- ncol(X.full)
n <- nrow(X.full)
B=50

# needs p+1 for intercept
#beta.lasso <- matrix(NA, ncol=2*B, nrow=p+1)
auc.lasso <- c()

for(j in 1:B){
  
  #set.seed(j*10)
  ss=sample.int(n,n/2)
  Y1=Y[ss]
  X1=X.full[ss,]
  Y2=Y[-ss]
  X2=X.full[-ss,]
  
  lambda_seq <- 10^seq(2, -2, length.out = 100)
  
  
  cv <- cv.glmnet(X1, Y1, alpha = 1, family = 'multinomial', lambda = lambda_seq)
  lambda.best <- cv$lambda.min
  fit.lasso <- glmnet(X1, Y1, alpha = 1,   family = 'multinomial', 
                      lambda = lambda.best, standardize = TRUE)
  
  
  predictions <- predict(fit.lasso, X2, type = 'class')
  predictions <- as.factor(predictions)
  
  auc.lasso[j] <- as.numeric(multiclass.roc(Y2, as.numeric(predictions))$auc)
  
  
  
  
  
}
mean(auc.lasso) #0.68 with 2/3 of the data, 0. with 1/2 0.647
