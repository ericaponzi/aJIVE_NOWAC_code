# check classification from joint components 
# and individual (?)

# load data
load('work/multiomics/erica/data.ajive.RData')
load('work/multiomics/erica/covs.ajive.RData')

# load ajive results
load('work/multiomics/erica/results/ajiveAssocLast.RData')

rj <- ajiveResults$joint_rank
r1 <- ajiveResults$block_decomps[[1]]$individual$rank
r2 <- ajiveResults$block_decomps[[2]]$individual$rank
r3 <- ajiveResults$block_decomps[[3]]$individual$rank

jointscores <- unlist(ajiveResults$joint_scores)
indivscores1 <- ajiveResults$block_decomps[[1]]$individual$u[,1:r1]
indivscores2 <- ajiveResults$block_decomps[[2]]$individual$u[,1:r2]
indivscores3 <- ajiveResults$block_decomps[[3]]$individual$u[,1:r3]


# prediction on outcome

# case vs control

data.logistic <- as.data.frame(cbind(y = covs$case_ctrl, jointscores, covs$age.sample, covs$BMI, 
                                     covs$smoking_status_cat,
                                     indivscores1, indivscores2, 
                                     indivscores3
))

data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)
# repeat until here for rf below


xvars <- names(data.logistic[2:(rj+r1+r2+r3+4)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )
#formula <- y ~ V7 + V8 +V9
model.fit <- glm(formula, family = 'binomial', 
                 data = data.logistic)
summary(model.fit)
data.logistic$pred <- predict(model.fit, data.logistic, type = 'response')
roc(data.logistic$y, data.logistic$pred)$auc

auc <- c()
K <- 10
for (k in 1:K){
  sample <- sample(nrow(data.logistic), floor(nrow(data.logistic)*0.9))
  data.train <- data.logistic[sample, ]
  data.test <- data.logistic[-sample, ]
  model.fit <- glm(formula ,
                   family = 'binomial', 
                   data = data.train)
  #summary(model.fit)
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  library(pROC)
  auc[k] <- roc(data.test$y, data.test$pred)$auc
}
mean(auc) 
# 0.76
# when adding covariates (smoking, age and bmi) 0.77

# smoking status
# three categories 
library(nnet)
data.logistic <- as.data.frame(cbind(y = covs$smoking_status_cat, covs$age.sample, covs$BMI, 
                                     jointscores, 
                                     indivscores1, indivscores2, 
                                     indivscores3
))

data.logistic$y <- as.factor(data.logistic$y)
xvars <- names(data.logistic[2:(rj+r1+r2+r3+3)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )

model.fit <- multinom(formula, 
                      data = data.logistic)
summary(model.fit)
data.logistic$pred <- predict(model.fit, data.logistic, 'class')
multiclass.roc(data.logistic$y~as.numeric(data.logistic$pred))

auc <- c()
K <- 10
for (k in 1:K){
  sample <- sample(nrow(data.logistic), floor(nrow(data.logistic)*0.9))
  data.train <- data.logistic[sample, ]
  data.test <- data.logistic[-sample, ]
  model.fit <- multinom(formula ,
                        data = data.train)
  #summary(model.fit)
  
  data.test$pred <- predict(model.fit, data.test, type = 'class')
  library(pROC)
  auc[k] <- multiclass.roc(data.test$y, as.numeric(data.test$pred))$auc
}
mean(auc)  # 0.73
# same with age and bmi

# age 
data.lm <- as.data.frame(cbind(y = covs$age.diag, jointscores, 
                                     indivscores1, indivscores2, 
                                    indivscores3
))

data.lm$logy <- log(data.lm$y)
xvars <- names(data.lm[2:(rj+r1+r2+r3+1)])
formula <- paste( 'logy', '~', paste( xvars, collapse=' + ' ) )

model.fit <- lm(formula, 
                      data = data.lm)
summary(model.fit)
data.lm$pred <- predict(model.fit, data.lm)
plot(data.lm$y~data.lm$pred)
# age not a good idea
# only between 50 and 60 approx
# same for BMI 
# only R squared between 0.1 and 0.2


## random forest prediction
library(randomForest)
ind = sample(2, nrow(data.logistic), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.logistic[ind == 1, ]
test.data = data.logistic[ind == 2, ]

rf <- randomForest(y ~., data = train.data, 
                   ntree = 100, proximity = TRUE, 
                   na.action = na.omit)
predictions <- predict(rf, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf), train.data$y)
cm <- table(predictions, test.data$y)
accuracy = sum(diag(cm))/sum(cm)
print(rf)
varImpPlot(rf, sort = TRUE, n.var = 10)
# first joint component is the first in variable importance
# second is smoking
# third is the 5th joint component
# all the others are individual


# check prediction in time from diagnosis
# split data in two groups
# cases based on median time to diag
# matched controls
covs$years.to.diag <- covs$age.diag - covs$age.sample
covs$time.to.diag <- as.numeric((covs$diagnosedato1 -covs$inndato)/365)
group1 <- covs[covs$case_ctrl == 'case' & 
                 covs$time.to.diag < median(covs$time.to.diag, na.rm = TRUE), ]
id.group1 <- group1$pairs
group1.caco <- covs[is.element(covs$pairs, id.group1), ]
group2.caco <- covs[!is.element(covs$pairs, id.group1), ]
ids.1 <- which(is.element(covs$pairs, id.group1))
# logistic model auc on group 1
jointscores1 <- jointscores[ids.1, ]
indivscores1.1 <- indivscores1[ids.1, ]
indivscores2.1 <- indivscores2[ids.1, ]
indivscores3.1 <- indivscores3[ids.1, ]
data.logistic <- as.data.frame(cbind(y = group1.caco$case_ctrl, jointscores1, group1.caco$age.sample, 
                                     group1.caco$BMI, 
                                     group1.caco$smoking_status_cat,
                                     indivscores1.1, indivscores2.1, 
                                     indivscores3.1
))

data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)

xvars <- names(data.logistic[2:(rj+r1+r2+r3+4)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )
#formula <- y ~ V7 + V8 +V9
model.fit <- glm(formula, family = 'binomial', 
                 data = data.logistic)
summary(model.fit)
data.logistic$pred <- predict(model.fit, data.logistic, type = 'response')
roc(data.logistic$y, data.logistic$pred)$auc

auc <- c()
K <- 10
for (k in 1:K){
  sample <- sample(nrow(data.logistic), floor(nrow(data.logistic)*0.9))
  data.train <- data.logistic[sample, ]
  data.test <- data.logistic[-sample, ]
  model.fit <- glm(formula ,
                   family = 'binomial', 
                   data = data.train)
  #summary(model.fit)
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  library(pROC)
  auc[k] <- roc(data.test$y, data.test$pred)$auc
}
mean(auc) 
# 0.64

library(randomForest)
data.logistic <- data.logistic[, -43]
ind = sample(2, nrow(data.logistic), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.logistic[ind == 1, ]
test.data = data.logistic[ind == 2, ]

rf <- randomForest(y ~., data = train.data, 
                   ntree = 100, proximity = TRUE, 
                   na.action = na.omit)
predictions <- predict(rf, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf), train.data$y)
cm <- table(predictions, test.data$y)
accuracy = sum(diag(cm))/sum(cm)
# accuracy 0.69
print(rf)
varImpPlot(rf, sort = TRUE, n.var = 10)
# first is smoking, 2nd and 3rd (and 10th) are joint components


# logistic model auc on group 2
jointscores2 <- jointscores[-ids.1, ]
indivscores1.2 <- indivscores1[-ids.1, ]
indivscores2.2 <- indivscores2[-ids.1, ]
indivscores3.2 <- indivscores3[-ids.1, ]
data.logistic <- as.data.frame(cbind(y = group2.caco$case_ctrl, jointscores2, group2.caco$age.sample, 
                                     group2.caco$BMI, 
                                     group2.caco$smoking_status_cat,
                                     indivscores1.2, indivscores2.2, 
                                     indivscores3.2
))

data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)

xvars <- names(data.logistic[2:(rj+r1+r2+r3+4)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )
#formula <- y ~ V7 + V8 +V9
model.fit <- glm(formula, family = 'binomial', 
                 data = data.logistic)
summary(model.fit)
data.logistic$pred <- predict(model.fit, data.logistic, type = 'response')
roc(data.logistic$y, data.logistic$pred)$auc

auc <- c()
K <- 10
for (k in 1:K){
  sample <- sample(nrow(data.logistic), floor(nrow(data.logistic)*0.9))
  data.train <- data.logistic[sample, ]
  data.test <- data.logistic[-sample, ]
  model.fit <- glm(formula ,
                   family = 'binomial', 
                   data = data.train)
  #summary(model.fit)
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  library(pROC)
  auc[k] <- roc(data.test$y, data.test$pred)$auc
}
mean(auc)
# 0.65

library(randomForest)
data.logistic <- data.logistic[, -43]

ind = sample(2, nrow(data.logistic), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.logistic[ind == 1, ]
test.data = data.logistic[ind == 2, ]

rf <- randomForest(y ~., data = train.data, 
                   ntree = 100, proximity = TRUE, 
                   na.action = na.omit)
predictions <- predict(rf, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf), train.data$y)
cm <- table(predictions, test.data$y)
accuracy = sum(diag(cm))/sum(cm)
# 0.76
print(rf)
varImpPlot(rf, sort = TRUE, n.var = 10)
# no smoking here, 2nd and 5th are joint