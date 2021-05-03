

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




