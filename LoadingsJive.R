# this function extracts loadings from jive
# joint and individual
# it also returns the variables with higher loadings
# k is the number of variables we want to extract


extract_loadings_jive <- function(jiveResults, k = 10){
  
  # extract joint loadings
  n_joint <- jiveResults$rankJ
  SVD = svd(do.call(rbind,jiveResults$joint), nu=n_joint, nv=n_joint)
  joint.load <- SVD$u
  
  # extract individual loadings for each source
  indiv.load <- list()
  n_indiv <- list()
  for (i in 1:length(jiveResults$rankA)){
    n_indiv[[i]] <- jiveResults$rankA[i] 
    SVDI = svd(jiveResults$individual[[i]],nu=n_indiv[[i]],nv=n_indiv[[i]])
    indiv.load[[i]] <- SVDI$u
  }
  
  loadings <- list(joint.load, indiv.load)
  names(loadings) <- c('Joint', 'Indiv')
  
  # visualize variables with higher importance in joint
  abs.loadings.joint <- abs(joint.load)
  topload.joint <- order(abs.loadings.joint, decreasing = TRUE)
  varnames <- c(rownames(dataNOWAC[[1]]), rownames(dataNOWAC[[2]]), 
                rownames(dataNOWAC[[3]]))
  top.joint.var <- varnames[topload.joint[1:k]]
  
  # visualize variables with higher importance in individual
  abs.loadings.indiv <- lapply(indiv.load, abs)
  topload.indiv <- list()
  top.indiv.var <- list()
  for (i in 1:length(dataNOWAC)){
    topload.indiv[[i]] <- apply(abs.loadings.indiv[[i]],2, 
                                function(x) order(x, decreasing = TRUE))
    
    top.indiv.var[[i]] <- rownames(dataNOWAC[[i]])[topload.indiv[[i]][1:k]]
  }
  
  top.var <- list(top.joint.var, top.indiv.var)
  names(top.var) <- c('Joint', 'Indiv')
  
  res <- list(loadings, top.var)
  names(res) <- c('Loadings', 'TopVariables')
  return(res)
  
}
