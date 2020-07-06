# single PCA on each data source
library(pROC)
library(randomForest)

load('data.ajive.RData')
load('covs.ajive.RData')

methylation <- data.ajive[[1]]
mRNA <- data.ajive[[2]]
miRNA <- data.ajive[[3]]

dim(methylation)
dim(mRNA)
dim(miRNA)

# pca on methylation
met.pca <- prcomp(methylation, center = TRUE, scale = TRUE)
str(met.pca)
met.pcs <- met.pca$x[, 1:5]

# pca on mRNA
mrna.pca <- prcomp(mRNA, center = TRUE, scale = TRUE)
mrna.pcs <- mrna.pca$x[,1:5]

# pca on miRNA
mirna.pca <- prcomp(miRNA, center = TRUE, scale = TRUE)
mirna.pcs <- mirna.pca$x[,1:5]


# prediction models with single pca 
data.logistic <- as.data.frame(cbind(y = covs$case_ctrl, covs$age.sample, covs$BMI, 
                                     covs$smoking_status_cat,
                                     met.pcs, mrna.pcs, 
                                     mirna.pcs))

data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)

names(data.logistic) <- c('y', 'age', 'bmi', 'smoke', 
                          'PC1met', 'PC2met', 'PC3met', 'PC4met', 'PC5met', 
                          'PC1mr', 'PC2mr', 'PC3mr', 'PC4mr', 'PC5mr',
                          'PC1mi', 'PC2mi', 'PC3mi', 'PC4mi', 'PC5mi')

xvars <- names(data.logistic[2:ncol(data.logistic)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )


xvars0 <- names(data.logistic[5:ncol(data.logistic)])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )




# pre set folds
# shuffle data
data.logistic <- data.logistic[sample(nrow(data.logistic)),]

# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)

# CV

auc <- c()
roc <- list()
auc0 <- c()
roc0 <- list()
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
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  
  roc[[i]] <- roc(data.test$y, data.test$pred)
  
  auc[i] <- roc(data.test$y, data.test$pred)$auc
  
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  
  roc0[[i]] <- roc(data.test$y, data.test$pred0)
  
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  
}
mean(auc)
mean(auc0)
sd(auc0)



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
accuracy = sum(diag(cm))/sum(cm)


print(rf)
varImpPlot(rf, sort = TRUE, n.var = 10, lty = 1)
r <- roc(train.data$y, rf$votes[,2])
r

train.data <- train.data[,-c(2,3,4)]
test.data <- test.data[,-c(2,3,4)]
rf <- randomForest(y ~., data = train.data, 
                   ntree = 500, proximity = TRUE, 
                   na.action = na.omit)
predictions <- predict(rf, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf), train.data$y)
cm <- table(predictions, test.data$y)
accuracy0 = sum(diag(cm))/sum(cm)


print(rf)
varImpPlot(rf, sort = TRUE, n.var = 10, lty = 1)
r <- roc(train.data$y, rf$votes[,2])
r


# load clinical covariates
clin <- readRDS('NOWAC_covariates_database_laboratory_160320.rds')
# merge covs and clin
clin$patient.id <- paste(clin$CaseStatus, clin$Pairs, sep = ' ')
rownames(clin) <- clin$patient.id

covs.all <- merge(covs, clin, sort = F)
# response is  metastasis or not
covs.all$metastasis.gr <- factor(covs.all$Metastasis > '0')


data.logistic <- as.data.frame(cbind(y = covs.all$metastasis.gr,covs$age.sample, covs$BMI, 
                                     covs$smoking_status_cat,
                                     met.pcs, mrna.pcs, 
                                     mirna.pcs))

data.logistic <- data.logistic[which(!is.na(data.logistic$y)), ]
data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)

names(data.logistic) <- c('y', 'age', 'bmi', 'smoke', 
                          'PC1met', 'PC2met', 'PC3met', 'PC4met', 'PC5met', 
                          'PC1mr', 'PC2mr', 'PC3mr', 'PC4mr', 'PC5mr',
                          'PC1mi', 'PC2mi', 'PC3mi', 'PC4mi', 'PC5mi')

xvars <- names(data.logistic[2:ncol(data.logistic)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )

xvars0 <- names(data.logistic[5:ncol(data.logistic)])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )


# pre set folds
# shuffle data
data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)



# CV

auc <- c()
roc <- list()

auc0 <- c()
roc0 <- list()

for (i in 1:10){
  testIndexes <- which(folds == i, arr.ind = TRUE)
  data.test <- data.logistic[testIndexes, ]
  data.train <- data.logistic[-testIndexes, ]
  model.fit <- glm(formula ,
                   family = 'binomial', 
                   data = data.train)
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  
  roc[[i]] <- roc(data.test$y, data.test$pred)
  
  auc[i] <- roc(data.test$y, data.test$pred)$auc
  
  model.fit0 <- glm(formula0 ,
                   family = 'binomial', 
                   data = data.train)
  
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  
  roc0[[i]] <- roc(data.test$y, data.test$pred0)
  
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  
}
mean(auc)
mean(auc0)




