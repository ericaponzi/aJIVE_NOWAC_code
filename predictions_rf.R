# predictions via random forest

set.seed(052020)

library(pROC)
library(randomForest)

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

# 1 = caco
# 2, 3, 4, 5, 6 = joint
# 7 = age
# 8 = bmi
# 9 = smoke
# 10, 11, 12, 13, 14 = Ind met
# 15, 16, 17, 18, 19 = Ind mRNA
# 20, 21, 22, 23, 24 = Ind miRNA


data.logistic$y[data.logistic$y == '2'] <- '0'
data.logistic$y <- as.factor(data.logistic$y)


xvars <- names(data.logistic[2:(rj+r1+r2+r3+4)])
formula <- paste( 'y', '~', paste( xvars, collapse=' + ' ) )

xvars0 <- names(data.logistic[(rj+2):(rj+4)])
formula0 <- paste( 'y', '~', paste( xvars0, collapse=' + ' ) )

xvars.j <- names(data.logistic[2:(rj+4)])
formula.j <- paste( 'y', '~', paste( xvars.j, collapse=' + ' ) )


ind = sample(2, nrow(data.logistic), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.logistic[ind == 1, ]
test.data = data.logistic[ind == 2, ]

rf <- randomForest(y ~., data = train.data, 
                   ntree = 500, proximity = TRUE, 
                   na.action = na.omit)


rf0 <- randomForest(y~., data = train.data[, c('y', xvars0)], 
                    ntree = 500, proximity = TRUE, 
                    na.action = na.omit)

rf.j <- randomForest(y~., data = train.data[, c('y', xvars.j)], 
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


predictions0 <- predict(rf0, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf0), train.data$y)
cm0 <- table(predictions0, test.data$y)
accuracy0 = sum(diag(cm0))/sum(cm0)
print(rf0)
varImpPlot(rf0, sort = TRUE, n.var = 10)
r0 <- roc(train.data$y, rf0$votes[,2])

predictions.j <- predict(rf.j, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf.j), train.data$y)
cm.j <- table(predictions.j, test.data$y)
accuracy.j = sum(diag(cm.j))/sum(cm.j)
print(rf.j)
varImpPlot(rf.j, sort = TRUE, n.var = 10)
r.j <- roc(train.data$y, rf.j$votes[,2])

par(oma = c(4, 1, 1, 1))
plot(r, col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.8, identity = FALSE, 
     main = 'Case vs control')
plot(r0, col = 'blue', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.8)
plot(r.j, col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.8)


imp <- importance(rf)
ord <- order(imp, decreasing = TRUE)[1:10]
# 1 = caco
# 2, 3, 4, 5, 6 = joint
# 7 = age
# 8 = bmi
# 9 = smoke
# 10, 11, 12, 13, 14 = Ind met
# 15, 16, 17, 18, 19 = Ind mRNA
# 20, 21, 22, 23, 24 = Ind miRNA

shape <- rep('darkgrey', 10)
shape[c(1, 3, 5)] <- 'blue'

varimp <- data.frame(imp = imp[ord], var = c('Joint1', 'miRNAInd5', 'Joint4', 
                                             'Smoke', 'Joint5', 'methylInd3', 
                                             'methylInd2', 'mRNAInd4', 'mRNAInd3',
                                             'miRNAInd2'), shape = shape)

varimp <- varimp[order(varimp$imp), ]


  
qplot(x = imp, 
      y = reorder(var, imp), color = factor(shape), data = varimp, geom = "point") +  geom_point(size = 2) +
  xlab("MeanDecreaseGini") + ylab('') + theme(legend.position = "none") +
  ggtitle('Variable Importance: case vs control') + geom_segment(aes(x = 0, y = var, xend = imp, yend = var ), lty = 2)
  
 

# metastasis
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


ind = sample(2, nrow(data.logistic), replace = TRUE, prob = c(0.9, 0.1))
train.data = data.logistic[ind == 1, ]
test.data = data.logistic[ind == 2, ]
l = list(ind, train.data, test.data, rf, rf0, rf.j)
save(l, file= 'rf.met.RData')

rf <- randomForest(y ~., data = train.data, 
                   ntree = 500, proximity = TRUE, 
                   na.action = na.omit)


rf0 <- randomForest(y~., data = train.data[, c('y', xvars0)], 
                    ntree = 500, proximity = TRUE, 
                    na.action = na.omit)

rf.j <- randomForest(y~., data = train.data[, c('y', xvars.j)], 
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


predictions0 <- predict(rf0, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf0), train.data$y)
cm0 <- table(predictions0, test.data$y)
accuracy0 = sum(diag(cm0))/sum(cm0)
print(rf0)
varImpPlot(rf0, sort = TRUE, n.var = 10)
r0 <- roc(train.data$y, rf0$votes[,2])

predictions.j <- predict(rf.j, newdata = test.data)
train.data <- na.omit(train.data)
table(predict(rf.j), train.data$y)
cm.j <- table(predictions.j, test.data$y)
accuracy.j = sum(diag(cm.j))/sum(cm.j)
print(rf.j)
varImpPlot(rf.j, sort = TRUE, n.var = 10)
r.j <- roc(train.data$y, rf.j$votes[,2])

par(oma = c(4, 1, 1, 1))
plot(r, col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.8, identity = FALSE, 
     main = 'Metastasis (yes vs no)')
plot(r0, col = 'blue', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.8)
plot(r.j, col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.8)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",bty = "n", inset=c(0,0.02), xpd = TRUE, ncol = 2,  legend=c("Joint and Individual Components","Joint Components", 
                                                                            "Patient Covariates"),
       lty=c(5,5,5), col = c("black", "red", "blue"), cex = 0.8)


 

imp <- importance(rf)
ord <- order(imp, decreasing = TRUE)[1:10]
# 1 = caco
# 2, 3, 4, 5, 6 = joint
# 7 = age
# 8 = bmi
# 9 = smoke
# 10, 11, 12, 13, 14 = Ind met
# 15, 16, 17, 18, 19 = Ind mRNA
# 20, 21, 22, 23, 24 = Ind miRNA

shape <- rep('darkgrey', 10)
shape[c(5, 8, 10)] <- 'blue'

varimp <- data.frame(imp = imp[ord], var = c('mRNAInd1', 'mRNAInd3', 'methylInd1', 
                                             'miRNAInd2', 'Joint4', 'Age', 
                                             'mRNAInd2', 'Joint5', 'methylInd4',
                                             'Joint3'), shape = shape)

varimp <- varimp[order(varimp$imp), ]



qplot(x = imp, 
      y = reorder(var, imp), color = factor(shape), data = varimp, geom = "point") +  geom_point(size = 2) +
  xlab("MeanDecreaseGini") + ylab('') + theme(legend.position = "none") +
  ggtitle('Variable Importance: metastasis') + geom_segment(aes(x = 0, y = var, xend = imp, yend = var ), lty = 2)


