# check classification of clinical covariates
# from joint components 
# and individual 
set.seed(230420)
library(pROC)
library(randomForest)
library(nnet)

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



# load clinical covariates
clin <- readRDS('work/phenotype/NOWAC/data/NOWAC_covariates_database_laboratory_160320.rds')

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

model.fit <- glm(formula, family = 'binomial', 
                 data = data.logistic)
summary(model.fit)
data.logistic$pred <- predict(model.fit, data.logistic, type = 'response')
roc(data.logistic$y, data.logistic$pred)$auc

data.logistic <- data.logistic[,-ncol(data.logistic)]
auc <- c()
K <- 1000
for (k in 1:K){
  sample <- sample(nrow(data.logistic), floor(nrow(data.logistic)*0.9))
  data.train <- data.logistic[sample, ]
  data.test <- data.logistic[-sample, ]
  model.fit <- glm(formula ,
                   family = 'binomial', 
                   data = data.train)
  #summary(model.fit)
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  
  auc[k] <- roc(data.test$y, data.test$pred)$auc
}
mean(auc) 
#0.59

## random forest prediction

ind = sample(2, nrow(data.logistic), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.logistic[ind == 1, ]
test.data = data.logistic[ind == 2, ]

rf <- randomForest(y ~., data = train.data, 
                   ntree = 500, proximity = TRUE, 
                   na.action = na.omit)
predictions <- predict(rf, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf), train.data$y)
cm <- table(predictions, test.data$y)
accuracy = sum(diag(cm))/sum(cm) #0.86
print(rf)
varImpPlot(rf, sort = TRUE, n.var = 10)
# in the first ten the second joint and the third joint and bmi
# all the others are individual

# same on histological subtypes
# three categories 
smallcellICD <- c(804133, 804134, 804139, 804234, 804239, 804633, 804634, 804639)
squamousICD <- c(807031, 807032, 807033, 807039, 807134, 807233, 807132)
adenoICD <- c(814031, 814032, 814033, 814039, 820039, 824039, 824633, 825031, 825032, 825039, 826039, 848131, 848139, 855032, 856039)
noICD <- c(800034, 800039, 801033, 802034)
covs.all$subtype <- ifelse(covs.all$CaseStatus == 'case' & covs.all$Morfologi_ICDo31 %in% smallcellICD, 'Smallcell', NA)
covs.all$subtype <- ifelse(covs.all$CaseStatus == 'case' & covs.all$Morfologi_ICDo31 %in% squamousICD, 'Squamous', covs.all$subtype)
covs.all$subtype <- ifelse(covs.all$CaseStatus == 'case' & covs.all$Morfologi_ICDo31 %in% adenoICD, 'Adeno', covs.all$subtype)
covs.all$subtype <- ifelse(covs.all$CaseStatus == 'case' & covs.all$Morfologi_ICDo31 %in% noICD, 'Others', covs.all$subtype)



data.logistic <- data.frame(y = covs.all$subtype,
                                     jointscores, covs.all$age.sample, covs.all$BMI, 
                                     covs.all$smoking_status_cat,
                                     indivscores1, indivscores2, 
                                     indivscores3
)
names(data.logistic)
data.logistic <- data.logistic[which(!is.na(data.logistic$y)), ]
data.logistic$y <- as.factor(data.logistic$y)
xvars <- names(data.logistic[2:(rj+r1+r2+r3+4)])
#(
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )

model.fit <- multinom(formula, 
                      data = data.logistic)
summary(model.fit)
data.logistic$pred <- predict(model.fit, data.logistic, 'class')
multiclass.roc(data.logistic$y~as.numeric(data.logistic$pred))

data.logistic <- data.logistic[, -ncol(data.logistic)]
auc <- c()
K <- 1000
for (k in 1:K){
  sample <- sample(nrow(data.logistic), floor(nrow(data.logistic)*0.9))
  data.train <- data.logistic[sample, ]
  data.test <- data.logistic[-sample, ]
  model.fit <- multinom(formula ,
                        data = data.train)
  #summary(model.fit)
  
  data.test$pred <- predict(model.fit, data.test, type = 'class')
 
  auc[k] <- multiclass.roc(data.test$y, as.numeric(data.test$pred))$auc
}
mean(auc)  # 0.67

## random forest prediction

ind = sample(2, nrow(data.logistic), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.logistic[ind == 1, ]
test.data = data.logistic[ind == 2, ]

rf <- randomForest(y ~., data = train.data, 
                   ntree = 500, proximity = TRUE, 
                   na.action = na.omit)
predictions <- predict(rf, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf), train.data$y)
cm <- table(predictions, test.data$y)
accuracy = sum(diag(cm))/sum(cm) #0.5
print(rf)
varImpPlot(rf, sort = TRUE, n.var = 20)

# first  is the first joint, second joint also in, age, smoking and bmi in all others are individual
