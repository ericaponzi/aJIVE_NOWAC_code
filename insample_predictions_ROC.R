# roc curves for paper



# case vs control
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

# non integrated analysis
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
std_dev <- met.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex[1:20], xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

# pca on mRNA
mrna.pca <- prcomp(mRNA, center = TRUE, scale = TRUE)
mrna.pcs <- mrna.pca$x[,1:5]
std_dev <- mrna.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex[1:20], xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")


# pca on miRNA
mirna.pca <- prcomp(miRNA, center = TRUE, scale = TRUE)
mirna.pcs <- mirna.pca$x[,1:5]
std_dev <- mirna.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex[1:20], xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")


# prediction on outcome

# case vs control
data.logistic <- as.data.frame(cbind(y = covs$case_ctrl,  jointscores, covs$age.sample, covs$BMI, 
                                     covs$smoking_status_cat,
                                     indivscores1, indivscores2, 
                                     indivscores3, 
                                     met.pcs, mrna.pcs, 
                                     mirna.pcs))

names(data.logistic)[(rj+r1+r2+r3+5):ncol(data.logistic)] <- c('PC1met', 'PC2met', 'PC3met', 'PC4met', 'PC5met', 
                                                               'PC1mr', 'PC2mr', 'PC3mr', 'PC4mr', 'PC5mr',
                                                               'PC1mi', 'PC2mi', 'PC3mi', 'PC4mi', 'PC5mi')

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


xvars.ni.cov <- names(data.logistic[c((rj+2):(rj+4), (rj+r1+r2+r3+5):ncol(data.logistic))])
formula.ni.cov <- paste( 'y', '~', paste( xvars.ni.cov, collapse=' + ' ) )

xvars.ni <- names(data.logistic[(rj+r1+r2+r3+5):ncol(data.logistic)])
formula.ni <- paste( 'y', '~', paste( xvars.ni, collapse=' + ' ) )


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


par(oma = c(4, 1, 1, 1))
plot.roc(roc(data.test$y, data.test$pred), 
         col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.6, print.auc.cex = 0.9, identity = FALSE, 
         main = 'Case vs control', xlim = c(1, 0), xlab = "", ylab = "")
plot(roc(data.test$y, data.test$pred0), col = 'blue', 
     add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.j), 
     col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.5, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.i), 
     col = 'darkgreen', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.9)
#plot(roc(data.test$y, data.test$pred.ni), 
 #    col = 'darkorchid', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.ni.cov), 
     col = 'darkorange', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.9)


# load clinical covariates
clin <- readRDS('NOWAC_covariates_database_laboratory_160320.rds')
# merge covs and clin
clin$patient.id <- paste(clin$CaseStatus, clin$Pairs, sep = ' ')
rownames(clin) <- clin$patient.id

covs.all <- merge(covs, clin, sort = F)

# response is  metastasis or not
covs.all$metastasis.gr <- factor(covs.all$Metastasis > '0')


data.logistic <- as.data.frame(cbind(y = covs.all$metastasis.gr,  jointscores, covs.all$age.sample, covs.all$BMI, 
                                     covs.all$smoking_status_cat,
                                     indivscores1, indivscores2, 
                                     indivscores3, 
                                     met.pcs, mrna.pcs, 
                                     mirna.pcs))

names(data.logistic)[(rj+r1+r2+r3+5):ncol(data.logistic)] <- c('PC1met', 'PC2met', 'PC3met', 'PC4met', 'PC5met', 
                          'PC1mr', 'PC2mr', 'PC3mr', 'PC4mr', 'PC5mr',
                          'PC1mi', 'PC2mi', 'PC3mi', 'PC4mi', 'PC5mi')

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


xvars.ni.cov <- names(data.logistic[c((rj+2):(rj+4), (rj+r1+r2+r3+5):ncol(data.logistic))])
formula.ni.cov <- paste( 'y', '~', paste( xvars.ni.cov, collapse=' + ' ) )

xvars.ni <- names(data.logistic[(rj+r1+r2+r3+5):ncol(data.logistic)])
formula.ni <- paste( 'y', '~', paste( xvars.ni, collapse=' + ' ) )


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


par(oma = c(4, 1, 1, 1))
plot.roc(roc(data.test$y, data.test$pred), 
         col = 'black', lty = 5, print.auc = TRUE, print.auc.y = 0.6, print.auc.cex = 0.9, identity = FALSE, 
         main = 'Metastasis (yes vs no)', xlim = c(1, 0), xlab = "", ylab = "")
plot(roc(data.test$y, data.test$pred0), col = 'blue', 
     add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.4, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.j), 
     col = 'red', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.5, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.i), 
     col = 'darkgreen', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.3, print.auc.cex = 0.9)
#plot(roc(data.test$y, data.test$pred.ni), 
 #    col = 'darkorchid', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.9)
plot(roc(data.test$y, data.test$pred.ni.cov), 
     col = 'darkorange', add = TRUE, lty = 5, print.auc = TRUE, print.auc.y = 0.2, print.auc.cex = 0.9)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",bty = "n", inset=c(0,0.02), xpd = TRUE, ncol = 2,  
       
       legend=c("Joint and Individual Components",
                "Patient Covariates", 
                "Joint Components",
                "Individual Components", 
                "Non integrative"),
       lty = c(5,5,5,5,5),
       
       col = c("black", 
               "blue", 
               "red", 
               "darkgreen", 
               "darkorange")
       
       
       , cex = 0.7)


