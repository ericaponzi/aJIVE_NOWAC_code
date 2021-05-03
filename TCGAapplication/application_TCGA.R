

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


xvars.ni <- names(data.logistic[(rj+r1+r2+r3+2):ncol(data.logistic)])
formula.ni <- paste( 'Y', '~', paste( xvars.ni, collapse=' + ' ) )



#index <- sample(150, 100)
data.train <- data.logistic
model.fit <- multinom(formula , 
                 data = data.train)

model.fit.ni <- multinom(formula.ni, 
                    data = data.train)


data.test <- data.logistic
data.test$pred <- predict(model.fit, data.test, type = 'class')
data.test$pred.ni <- predict(model.fit.ni, data.test, type = 'class')

multiclass.roc(data.test$Y, as.numeric(data.test$pred))$auc
multiclass.roc(data.test$Y, as.numeric(data.test$pred.ni))$auc




# CV

# pre set folds
# shuffle data
set.seed(202011)
data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)


auc <- c()
auc.ni <- c()
roc <- list()
roc.ni <- list()


for (i in 1:10){
  testIndexes <- which(folds == i, arr.ind = TRUE)
  data.test <- data.logistic[testIndexes, ]
  data.train <- data.logistic[-testIndexes, ]
  model.fit <- multinom(formula, 
                   data = data.train)
  model.fit.ni <- multinom(formula.ni, 
                    data = data.train)
  
  
  
  data.test$pred <- predict(model.fit, data.test, type = 'class')
  data.test$pred.ni <- predict(model.fit.ni, data.test, type = 'class')
  
  
  roc[[i]] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred))
  roc.ni[[i]] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.ni))
  
  
  
  auc[i] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred))$auc
  auc.ni[i] <- multiclass.roc(data.test$Y, as.numeric(data.test$pred.ni))$auc
  
 
  
}
mean(auc)
mean(auc.ni)


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

