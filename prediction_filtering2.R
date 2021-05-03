
set.seed(122020)

library(pROC)
library(nnet)

load(file = 'C:/Users/ericapo/Desktop/NOWAC/0912/201209.ajive.case3rank71.14.10.RData')
load("C:/Users/ericapo/Desktop/NOWAC/ResultsFromHunt/covs.ajive.RData")
load("C:/Users/ericapo/Desktop/NOWAC/0912/data.ajive.case3_0912.RData")

ajiveResults <- ajiveResults2
rj <- ajiveResults$joint_rank
r1 <- ajiveResults$block_decomps[[1]]$individual$rank
r2 <- ajiveResults$block_decomps[[2]]$individual$rank
r3 <- ajiveResults$block_decomps[[3]]$individual$rank

# We only want 5 individual each
r1 <- 5
r2 <- 5
r3 <- 5

jointscores <- unlist(ajiveResults$joint_scores)
indivscores1 <- ajiveResults$block_decomps[[1]]$individual$u[,1:r1]
indivscores2 <- ajiveResults$block_decomps[[2]]$individual$u[,1:r2]
indivscores3 <- ajiveResults$block_decomps[[3]]$individual$u[,1:r3]

# non integrative analysis
methylation <- data.ajive[[1]]
mRNA <- data.ajive[[2]]
miRNA <- data.ajive[[3]]


# pca on methylation
met.pca <- prcomp(methylation, center = TRUE, scale = TRUE)
str(met.pca)
var <- (met.pca$sdev)*2
prop_varex <- var/sum(var)
plot(prop_varex[1:20])
met.pcs <- met.pca$x[, 1:5]

# pca on mRNA
mRNA.pca <- prcomp(mRNA, center = TRUE, scale = TRUE)

var <- (mRNA.pca$sdev)*2
prop_varex <- var/sum(var)
plot(prop_varex[1:20])
mRNA.pcs <- mRNA.pca$x[, 1:5]

# pca on miRNA
miRNA.pca <- prcomp(miRNA, center = TRUE, scale = TRUE)

var <- (miRNA.pca$sdev)*2
prop_varex <- var/sum(var)
plot(prop_varex[1:20])
miRNA.pcs <- miRNA.pca$x[, 1:5]

# prediction on outcome

# case vs control

data.logistic <- as.data.frame(cbind(y = covs$case_ctrl, jointscores, covs$age.sample, covs$BMI, 
                                     as.factor(covs$smoking_status_cat),
                                     indivscores1, indivscores2, 
                                     indivscores3, 
                                     met.pcs, mRNA.pcs, miRNA.pcs))


data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)

names(data.logistic)[(rj+r1+r2+r3+5):ncol(data.logistic)] <- c('PCmet1', 'PCmet2', 'PCmet3', 
                                                               'PCmet4', 'PCmet5', #'PCmet6', 'PCmet7',
                                                               'PCmR1', 'PCmR2', 'PCmR3', 
                                                               'PCmR4', 'PCmR5', #'PCmr6', 'PCmr7',
                                                               'PCmiR1', 'PCmiR2', 'PCmiR3', 
                                                               'PCmiR4', 'PCmiR5')

xvars <- names(data.logistic[2:(rj+r1+r2+r3+4)])
#xvars <- names(data.logistic[c(2:(rj+11), (rj+r1+5):(rj+r1+11), (rj+r1+r2+5):(rj+r1+r2+r3+4))])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )

xvars0 <- names(data.logistic[(rj+2):(rj+4)])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )

xvars.j <- names(data.logistic[2:(rj+4)])
formula.j <- paste( 'y', '~', paste( xvars.j, collapse=' + ' ) )

xvars.i <- names(data.logistic[(rj+2):(rj+r1+r2+r3+4)])
formula.i <- paste( 'y', '~', paste( xvars.i, collapse=' + ' ) )

xvars.ni.cov <- names(data.logistic[c((rj+2):(rj+4), (rj+r1+r2+r3+5):ncol(data.logistic))])
formula.ni.cov <- paste( 'y', '~', paste( xvars.ni.cov, collapse=' + ' ) )

xvars.ni <- names(data.logistic[(rj+r1+r2+r3+5):ncol(data.logistic)])
formula.ni <- paste( 'y', '~', paste( xvars.ni, collapse=' + ' ) )



# in sample ROC
library(pROC)
data.train <- data.logistic
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

model.fit.ni <- glm(formula.ni ,
                    family = 'binomial', 
                    data = data.train)

model.fit.ni.cov <- glm(formula.ni.cov,
                        family = 'binomial', 
                        data = data.train)

data.test <- data.logistic
data.test$pred <- predict(model.fit, data.test, type = 'response')
data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
data.test$pred.j <- predict(model.fit.j, data.test, type = 'response')
data.test$pred.i <- predict(model.fit.i, data.test, type = 'response')
data.test$pred.ni <- predict(model.fit.ni, data.test, type = 'response')
data.test$pred.ni.cov <- predict(model.fit.ni.cov, data.test, type = 'response')


par(oma = c(4,0,0,0), pty = 's')
plot.roc(roc(data.test$y, data.test$pred), 
         col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.6, print.auc.cex = 0.9, identity = FALSE, 
         main = 'Case vs control', legacy.axes = TRUE)
plot(roc(data.test$y, data.test$pred0), col = 'blue', 
     add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.9, legacy.axes = TRUE)
plot(roc(data.test$y, data.test$pred.j), 
     col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.5, print.auc.cex = 0.9, legacy.axes = TRUE)
plot(roc(data.test$y, data.test$pred.i), 
     col = 'darkgreen', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.9, legacy.axes = TRUE)
#plot(roc(data.test$y, data.test$pred.ni), 
 #    col = 'darkorchid', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.ni.cov), 
     col = 'darkorange', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.9, legacy.axes = TRUE)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",bty = "n", inset=c(0,0.02), xpd = TRUE, ncol = 2,  
       
       legend=c("Joint and Individual Components",
                "Patient Covariates", 
                "Joint Components",
                "Individual Components", 
                #"Non integrative", 
                "Non integrative "),
       lty = c(5,5,5,5,5),
       
       col = c("black", 
               "blue", 
               "red", 
               "darkgreen", 
               
               "darkorange")
       
       
       , cex = 0.8)

# 10 fold cross validation
# pre set folds
# shuffle data
data.logistic <- data.logistic[sample(nrow(data.logistic)),]
# 10 equally sized folds
folds <- cut(seq(1, nrow(data.logistic)), breaks = 10, labels = FALSE)
#load('work/multiomics/erica/results/20201118_10CVcaco.RData')
#folds <- l[[2]]
#data.logistic <- l[[1]]
# CV
auc <- c()
auc0 <- c()
auc.j <- c()
auc.i <- c()
auc.ni <- c()
auc.ni.cov <- c()


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
  
  model.fit.ni <- glm(formula.ni ,
                      family = 'binomial', 
                      data = data.train)
  
  model.fit.ni.cov <- glm(formula.ni.cov,
                          family = 'binomial', 
                          data = data.train)
  
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  data.test$pred.j <- predict(model.fit.j, data.test, type = 'response')
  data.test$pred.i <- predict(model.fit.i, data.test, type = 'response')
  data.test$pred.ni <- predict(model.fit.ni, data.test, type = 'response')
  data.test$pred.ni.cov <- predict(model.fit.ni.cov, data.test, type = 'response')
  
  auc[i] <- roc(data.test$y, data.test$pred)$auc
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  auc.j[i] <- roc(data.test$y, data.test$pred.j)$auc
  auc.i[i] <- roc(data.test$y, data.test$pred.i)$auc
  auc.ni[i] <- roc(data.test$y, data.test$pred.ni)$auc
  auc.ni.cov[i] <- roc(data.test$y, data.test$pred.ni.cov)$auc
  
  
  
}
mean(auc) #0.69
mean(auc0) #0.67
mean(auc.j) #0.64
mean(auc.i) #0.69
#mean(auc.ni) #0.63
mean(auc.ni.cov) #0.67

l <- list(data.logistic, folds)
save(l, file = '10CV.case.casegeneloc_last.RData')


# random forest
library(randomForest)
data.rf <- data.logistic[, 1:(rj+r1+r2+r3+4)]
ind = sample(2, nrow(data.rf), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.rf[ind == 1, ]
test.data = data.rf[ind == 2, ]

rf <- randomForest(y ~., data = train.data, 
                   ntree = 1000, proximity = TRUE, 
                   na.action = na.omit)
varImpPlot(rf, sort = TRUE, n.var = 10)
predictions <- predict(rf, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf), train.data$y)
cm <- table(predictions, test.data$y)
accuracy = sum(diag(cm))/sum(cm)
print(rf)




imp <- importance(rf)
ord <- order(imp, decreasing = TRUE)[1:10]
# 1 = caco
# 2, 3, 4, 5, 6 = joint
#  7, 8, 9 =  age, BMi, smoking
#10, 11, 12, 13, 14 = Ind met
# 15, 16, 17, 18, 19= Ind mRNA
# 20, 21, 22, 23, 24= Ind mirna

shape <- rep('darkgrey', 10)
shape[2] <- 'blue'
shape[c(3,6)] <- 'green'

varimp <- data.frame(imp = imp[ord], var = c('MetInd2', 'Smoking',  'Joint4',
                                             'MetInd3', 'miRNAInd3', 'Joint3',
                                             'mRNAInd2', 'mRNAInd4', 'mRNAInd3',
                                              'miRNAInd1'
                                             ), shape = shape)

varimp <- varimp[order(varimp$imp), ]


library(ggplot2)
qplot(x = imp, 
      y = reorder(var, imp), color = factor(shape), data = varimp, geom = "point") +  geom_point(size = 2) +
  xlab("MeanDecreaseGini") + ylab('') + theme(legend.position = "none") +
  ggtitle('Variable Importance') + geom_segment(aes(x = 0, y = var, xend = imp, yend = var ), lty = 2)

# save all results 

save(rf, file = 'RF_cacocase3.RData')

data.rf <- data.logistic
ind = sample(2, nrow(data.rf), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.rf[ind == 1, ]
test.data = data.rf[ind == 2, ]

rf <- randomForest(y ~., data = train.data[, c('y', xvars)], 
                   ntree = 1000, proximity = TRUE, 
                   na.action = na.omit)
predictions <- predict(rf, newdata = test.data[, c('y', xvars)])
train.data <- na.omit(train.data)
table(predict(rf), train.data$y)
cm <- table(predictions, test.data$y)
accuracy = sum(diag(cm))/sum(cm)
print(rf)
roc(test.data$y, as.numeric(predictions))
#roc(test.data$y, rf$votes[,2])


rf.i <- randomForest(y~., data = train.data[, c('y', xvars.i)], 
                     ntree = 1000, proximity = TRUE, 
                     na.action = na.omit)
predictions.i <- predict(rf.i, newdata = test.data[, c('y', xvars.i)])
train.data <- na.omit(train.data)
table(predict(rf.i), train.data$y)
cm.i <- table(predictions.i, test.data$y)
accuracy.i = sum(diag(cm.i))/sum(cm.i)
roc(test.data$y, as.numeric(predictions.i))
print(rf.i)
roc(train.data$y, rf.i$votes[,2])

rf.j <- randomForest(y~., data = train.data[, c('y', xvars.j)], 
                     ntree = 1000, proximity = TRUE, 
                     na.action = na.omit)
predictions.j <- predict(rf.j, newdata = test.data[, c('y', xvars.j)])
train.data <- na.omit(train.data)
table(predict(rf.j), train.data$y)
cm.j <- table(predictions.j, test.data$y)
accuracy.j = sum(diag(cm.j))/sum(cm.j)
roc(test.data$y, as.numeric(predictions.j))
print(rf.j)
roc(train.data$y, rf.j$votes[,2])

rf0 <- randomForest(y~., data = train.data[, c('y', xvars0)], 
                     ntree = 1000, proximity = TRUE, 
                     na.action = na.omit)
predictions0 <- predict(rf0, newdata = test.data[, c('y', xvars0)])
train.data <- na.omit(train.data)
table(predict(rf0), train.data$y)
cm0 <- table(predictions0, test.data$y)
accuracy0 = sum(diag(cm0))/sum(cm0)
roc(test.data$y, as.numeric(predictions0))
print(rf0)
roc(train.data$y, rf0$votes[,2])

rf.ni <- randomForest(y~., data = train.data[, c('y', xvars.ni.cov)], 
                      ntree = 1000, proximity = TRUE, 
                      na.action = na.omit)

predictions.ni <- predict(rf.ni, newdata = test.data[, c('y', xvars.ni.cov)])
train.data <- na.omit(train.data)
table(predict(rf.ni), train.data$y)
cm.ni <- table(predictions.ni, test.data$y)
accuracy.ni = sum(diag(cm.ni))/sum(cm.ni)
roc(test.data$y, as.numeric(predictions.ni))
print(rf.ni)
roc(train.data$y, rf.ni$votes[,2])

roc(train.data$y, rf$votes[,2])
roc(train.data$y, rf.j$votes[,2])
roc(train.data$y, rf.i$votes[,2])
roc(train.data$y, rf0$votes[,2])
roc(train.data$y, rf.ni$votes[,2])


accuracy
accuracy.i
accuracy.j
accuracy0
accuracy.ni


roc(test.data$y, as.numeric(predictions.ni))
roc(test.data$y, as.numeric(predictions.i))
roc(test.data$y, as.numeric(predictions))
roc(test.data$y, as.numeric(predictions.j))
roc(test.data$y, as.numeric(predictions0))



rf.list <- list(rf, rf.j, rf.i, rf.ni, rf0)
names(rf.list) <- c('rf', 'rf.j', 'rf.i', 'rf.ni', 'rf0')
save(rf.list, file = 'RFlistcase3_last.RData')



# metastasis
# load clinical covariates
clin <- readRDS('NOWAC_covariates_database_laboratory_160320.rds')

# merge covs and clin
clin$patient.id <- paste(clin$CaseStatus, clin$Pairs, sep = ' ')
rownames(clin) <- clin$patient.id

covs.all <- merge(covs, clin, sort = F)
# response is  metastasis or not
covs.all$metastasis.gr <- factor(covs.all$Metastasis > '0')

data.logistic <- as.data.frame(cbind(y = covs.all$metastasis.gr, jointscores, covs.all$age.sample, covs.all$BMI, 
                                     covs.all$smoking_status_cat,
                                     indivscores1, indivscores2, 
                                     indivscores3, 
                                     met.pcs, mRNA.pcs, miRNA.pcs))



#r2 <- 5
names(data.logistic)[(rj+r1+r2+r3+5):ncol(data.logistic)] <- c('PCmet1', 'PCmet2', 'PCmet3', 
                                                               'PCmet4', 'PCmet5',
                                                               'PCmR1', 'PCmR2', 'PCmR3', 
                                                               'PCmR4', 'PCmR5', #'PCmR6', 'PCmR7', 
                                                               #'PCmr8', 'PCmr9',
                                                               'PCmiR1', 'PCmiR2', 'PCmiR3', 
                                                               'PCmiR4', 'PCmiR5')





data.logistic <- data.logistic[which(!is.na(data.logistic$y)), ]
data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- factor(data.logistic$y)



xvars <- names(data.logistic[2:(rj+r1+r2+r3+4)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )

xvars0 <- names(data.logistic[(rj+2):(rj+4)])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )

xvars.j <- names(data.logistic[2:(rj+4)])
formula.j <- paste( 'y', '~', paste( xvars.j, collapse=' + ' ) )

xvars.i <- names(data.logistic[(rj+2):(rj+r1+r2+r3+4)])
formula.i <- paste( 'y', '~', paste( xvars.i, collapse=' + ' ) )

xvars.ni.cov <- names(data.logistic[c((rj+2):(rj+4), (rj+r1+r2+r3+5):ncol(data.logistic))])
formula.ni.cov <- paste( 'y', '~', paste( xvars.ni.cov, collapse=' + ' ) )

xvars.ni <- names(data.logistic[(rj+r1+r2+r3+5):ncol(data.logistic)])
formula.ni <- paste( 'y', '~', paste( xvars.ni, collapse=' + ' ) )



# in sample ROC
library(pROC)
data.train <- data.logistic
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

model.fit.ni <- glm(formula.ni ,
                    family = 'binomial', 
                    data = data.train)

model.fit.ni.cov <- glm(formula.ni.cov,
                        family = 'binomial', 
                        data = data.train)

data.test <- data.logistic
data.test$pred <- predict(model.fit, data.test, type = 'response')
data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
data.test$pred.j <- predict(model.fit.j, data.test, type = 'response')
data.test$pred.i <- predict(model.fit.i, data.test, type = 'response')
data.test$pred.ni <- predict(model.fit.ni, data.test, type = 'response')
data.test$pred.ni.cov <- predict(model.fit.ni.cov, data.test, type = 'response')


par(oma = c(4,0,0,0), pty = 's')
plot.roc(roc(data.test$y, data.test$pred), 
         col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.6, print.auc.cex = 0.9, identity = FALSE, 
         main = 'Metastasis (yes vs no)', legacy.axes = TRUE)
plot(roc(data.test$y, data.test$pred0), col = 'blue', 
     add = TRUE, lty = 5, print.auc = TRUE, legacy.axes = TRUE,print.auc.y = 0.4, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.j), 
     col = 'red', add = TRUE, lty = 5, legacy.axes = TRUE,print.auc = TRUE, print.auc.y = 0.5, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.i), 
     col = 'darkgreen', add = TRUE, lty = 5,legacy.axes = TRUE, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.9)
#plot(roc(data.test$y, data.test$pred.ni), 
#    col = 'darkorchid', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.ni.cov), 
     col = 'darkorange', add = TRUE, lty = 5, legacy.axes = TRUE,print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.9)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",bty = "n", inset=c(0,0.02), xpd = TRUE, ncol = 2,  
       
       legend=c("Joint and Individual Components",
                "Patient Covariates", 
                "Joint Components",
                "Individual Components", 
                #"Non integrative", 
                "Non integrative "),
       lty = c(5,5,5,5,5),
       
       col = c("black", 
               "blue", 
               "red", 
               "darkgreen", 
               
               "darkorange")
       
       
       , cex = 0.8)

# 10 fold cross validation
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
auc.ni <- c()
auc.ni.cov <- c()


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
  
  model.fit.ni <- glm(formula.ni ,
                      family = 'binomial', 
                      data = data.train)
  
  model.fit.ni.cov <- glm(formula.ni.cov,
                          family = 'binomial', 
                          data = data.train)
  
  
  data.test$pred <- predict(model.fit, data.test, type = 'response')
  data.test$pred0 <- predict(model.fit0, data.test, type = 'response')
  data.test$pred.j <- predict(model.fit.j, data.test, type = 'response')
  data.test$pred.i <- predict(model.fit.i, data.test, type = 'response')
  data.test$pred.ni <- predict(model.fit.ni, data.test, type = 'response')
  data.test$pred.ni.cov <- predict(model.fit.ni.cov, data.test, type = 'response')
  
  auc[i] <- roc(data.test$y, data.test$pred)$auc
  auc0[i] <- roc(data.test$y, data.test$pred0)$auc
  auc.j[i] <- roc(data.test$y, data.test$pred.j)$auc
  auc.i[i] <- roc(data.test$y, data.test$pred.i)$auc
  auc.ni[i] <- roc(data.test$y, data.test$pred.ni)$auc
  auc.ni.cov[i] <- roc(data.test$y, data.test$pred.ni.cov)$auc
  
  
  
}
mean(auc) #0.66
mean(auc0) #0.73
mean(auc.j) #0.71
mean(auc.i) #0.65
#mean(auc.ni) 
mean(auc.ni.cov) #0.65
l <- list(data.logistic, folds)
save(l, file = '10CVmet_casegeneloc.RData')
