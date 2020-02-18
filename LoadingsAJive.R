# this function extracts loadings from ajive
# joint and individual
# it also returns the variables with higher loadings
# k is the number of variables we want to extract


extract_loadings_ajive <- function(data, ajiveResults, k = 10){
  joint.load <- list()
  indiv.load <- list()
  
  for (i in 1:length(data)){
    joint.load[[i]] <- ajiveResults$block_decomps[[i]][['joint']][['v']]
    indiv.load[[i]] <- ajiveResults$block_decomps[[i]][['individual']][['v']]
  }
  
  loadings <- list(joint.load, indiv.load)
  names(loadings) <- c('Joint', 'Indiv')
  
  # visualize variables with higher importance 
  abs.loadings.indiv <- lapply(indiv.load, abs)
  abs.loadings.joint <- lapply(joint.load, abs)
  topload.indiv <- list()
  top.indiv.var <- list()
  topload.joint <- list()
  top.joint.var <- list()
  for (i in 1:length(length(data))){
    topload.indiv[[i]] <- apply(abs.loadings.indiv[[i]],2, 
                                function(x) order(x, decreasing = TRUE))
    
    top.indiv.var[[i]] <- colnames(data[[i]])[topload.indiv[[i]][1:k]]
    topload.joint[[i]] <- apply(abs.loadings.joint[[i]],2, 
                                function(x) order(x, decreasing = TRUE))
    
    top.joint.var[[i]] <- colnames(data[[i]])[topload.joint[[i]][1:k]]
  }
  top.var <- list(top.joint.var, top.indiv.var)
  names(top.var) <- c('Joint', 'Indiv')
  res <- list(loadings, top.var)
  names(res) <- c('Loadings', 'TopVariables')
  return(res)
  
}
