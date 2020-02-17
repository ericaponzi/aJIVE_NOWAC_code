# ajive wrapper
# prepares the data for ajive main function and runs i
# extracts results 

# we need to perform SVDmiss oon the data first
# need to center and scale the data

ajive.dataprep <- function(data){
  
  # needs transpose of the standard data for jive
  data.ajive <- list()
  for (l in 1:length(data)){
    X <- t(as.data.frame(data[[l]]))
    # svdmiss
    Xtemp <- SVDmiss(X)[[1]]
    Ximp0 <- Xtemp$u %*% diag(x = Xtemp$d) %*% t(Xtemp$v)
    # center values
    centerValues <- apply(Ximp0, 1, mean, na.rm = T)
    Ximp <- Ximp0 - matrix(rep(centerValues, 
                                 ncol(Ximp0)), nrow = nrow(Ximp0))
    # scale values
    n <- nrow(Ximp) * ncol(Ximp) 
    Ximp <- Ximp/norm(Ximp, type = "f") * sqrt(sum(n))
    data.ajive[[l]] <- Ximp
  }
  
    return(data.ajive)
}
