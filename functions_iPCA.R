# Ipca functions 
# all these functions are taken from: DataSlingers/iPCA
library(parallel)
library(ggplot2)
library(gridExtra)
library(randomForest)
library(viridis) 


# TRCMAimpute_iPCA: TRCMAimpute modified for iPCA
# 
# Inputs:
#  -x = list of K data matrices with missing values
#  -lamS = penalty for row covariance matrix (only needed if q == "multfrob")
#  -lamDs = penalty for column covariance matrix
#  -q = type of penalty for (column) covariance matrix, either 1, "addfrob", or "multfrob"
#  -rmvn.ed1s = list of K saved eigenvalue decomposition for Eig(t(xxc) %*% xxc) to speed up computations (optional)
#  -seed = seed
#  -maxit = maximum number of iterations
#  -trace = T/F (default F) whether or not to show loglikelihood stuff
#  -thr.em = threshold for EM algorithm
#  -thr.glasso = threshold for GLASSO
#  -maxit.glasso = max number of iterations for GLASSO
#
# Outputs:
#  -xhat = list of K imputed data matrices
#  -xh1 = initial imputation assuming Delta = I (not computed)
#  -xh.init = initial imputation assuming Sig = I
#  -Sigma = estimated Sigma row covariance
#  -Delta = estimated Delta column covariance
#  -Sigma.init = initial estimate of Sigma
#  -Delta.init = initial estimates for Deltas_ks
#  -M = estimated mean matrix
#  -loglike = loglikelihood values after initial imputation and one-step EM algorithm
#  -rmvn.ed1s = list of K saved eigenvalue decomposition for Eig(t(xxc) %*% xxc)

TRCMAimpute_iPCA <- function(x, lamS = NULL, lamDs, q = "multfrob", rmvn.ed1s,
                             
                             
                             seed=1, maxit=50, trace=FALSE, 
                             thr.em=1e-6, thr.glasso=1e-2, maxit.glasso=1) {
  
  # record dimensions
  n <- nrow(x[[1]])
  pks <- sapply(x, FUN = ncol)
  p <- sum(pks)
  K <- length(x)
  n_cores <- min(K,detectCores())
  
  set.seed(seed)
  
  # initial imputation: impute missing values assuming Sig = I
  cat("Initial Imputation... \n")
  rmvn.ans <- xinit <- list()
  Mhatinit <- lapply(X = x, 
                     FUN = function(X) {
                       mu <- colMeans(X, na.rm = T)
                       return(rep(1,n) %*% t(mu))
                     })
  Siginit <- diag(n)
  if (missing(rmvn.ed1s)) {
    #     for (k in 1:K) {
    #       cat(paste0("... for k = ",k, "\n"))
    #       rmvn.ans[[k]] <- RMVNimpute(x=x[[k]], rho=lamDs[k], q=q, trace=F,
    #                                   maxit=maxit, thr.em=thr.em, maxit.glasso=maxit.glasso, thr.glasso=thr.glasso)
    #     }
    rmvn.ans <- mcmapply(X = x, D = lamDs,
                         FUN = function(X,D) {return(RMVNimpute(X,D,q=q,trace=F,maxit=maxit,thr.em=thr.em,maxit.glasso=25,thr.glasso=thr.glasso))},
                         SIMPLIFY = FALSE,
                         mc.cores = n_cores)
    xinit <- lapply(X = rmvn.ans, FUN = function(X) {return(X$xhat)})
    Delthats <- Deltinits <- lapply(X = rmvn.ans, FUN = function(X) {return(X$Sig)})
    rmvn.ed1s <- lapply(X = rmvn.ans, FUN = function(X) {return(X$rmvn.ed1)})
  }else {
    xinit <- mcmapply(X = x, D = lamDs, E = rmvn.ed1s,
                      FUN = function(X,D,E) {return(RMVNimpute(X,D,q=q,ed1=E,trace=F,maxit=maxit,thr.em=thr.em,maxit.glasso=25,thr.glasso=thr.glasso)$xhat)},
                      SIMPLIFY = FALSE,
                      mc.cores = n_cores)
    Delthats <- Deltinits <- lapply(X = rmvn.ans, FUN = function(X) {return(X$Sig)})
  }
  
  xhat <- xinit
  like <- NA
  like_idx <- 1
  if (trace) {
    like[like_idx] <- MNloglike_iPCA(Xs = x, Ms = Mhatinit, Sig = diag(n), Delts = Delthats, lamS = lamS, lamDs = lamDs, q = q, print = T)
    like_idx <- like_idx + 1
  }
  
  # compute MLEs given the initial imputed xhat
  cat("Computing MLEs... \n")
  xc <- lapply(xhat, FUN = function(X) {return(scale(X,center = T, scale = F))}) # centered data
  muhat <- lapply(xc, FUN = function(X) {return(attr(X, "scaled:center"))}) # estimated column means
  Mhat <- lapply(muhat, FUN = function(X) {return(matrix(rep(X, times = n), nrow = n, byrow = T))})
  
  if (q == 1) {
    ans <- FFmleGlasso(dat = xc, lamDs = lamDs, lamS = lamS, maxit = maxit.glasso, thr = thr.glasso, parallel = T)
  }else if (q == "1_off") {
    ans <- FFmleGlasso(dat = xc, lamDs = lamDs, lamS = lamS, maxit = maxit.glasso, thr = thr.glasso, parallel = T, pen_diag = F)
  }else if (q == "addfrob") {
    ans <- FFmleAddFrob(dat = xc, lamDs = lamDs, lamS = lamS)
  }else if (q == "multfrob") {
    ans <- FFmleMultFrob(dat = xc, lamDs = lamDs)
  }else {
    stop("q must be either 1, 1_off, addfrob, or multfrob.")
  }
  
  Sighat <- ans$Sig
  Delthats <- ans$Delts
  
  # given MLEs for covariances, impute X using ECM algorithm from GA's paper
  cat("Running ECM Algorithm... \n")
  
  xhat <- mcmapply(X = x, pk = as.list(pks), xhatk = xhat, Delthatk = Delthats, Mhatk = Mhat, SIMPLIFY = F, mc.cores = n_cores,
                   FUN = function(X, pk, xhatk, Delthatk, Mhatk) {
                     ind <- 1; iter <- 1
                     rmi <- (1:n)[apply(is.na(X),1,sum)>0] # indices of rows with missing values
                     cmj <- (1:pk)[apply(is.na(X),2,sum)>0] # indices of columns with missing values
                     na_idx <- is.na(X)
                     
                     while (ind > thr.em & iter < maxit) { # same algorithm as GA but specialized to the case where nu = 0
                       oldx <- xhatk
                       oldxc <- scale(oldx, center = T, scale = F)
                       
                       #by rows
                       for(i in rmi) {
                         mi <- na_idx[i,]
                         sinvs <- t(solve(Sighat[-i,-i], Sighat[-i,i]))
                         psi <- Mhatk[i,] + sinvs %*% (xhatk[-i,] - Mhatk[-i,])
                         Gam <- (Sighat[i,i] - sinvs %*% Sighat[-i,i]) %x% Delthatk
                         ginvg <- t(solve(Gam[!mi,!mi], Gam[!mi,mi]))
                         xhatk[i,mi] <- psi[mi] + ginvg %*% (xhatk[i,!mi] - psi[!mi])
                       }
                       
                       #by cols
                       for(j in cmj) {
                         mj <- na_idx[,j] # T/F: missing or not
                         invdd <- solve(Delthatk[-j,-j], Delthatk[-j,j])
                         nu <- Mhatk[,j] + (xhatk[,-j] - Mhatk[,-j]) %*% invdd
                         Phi <- Sighat %x% (Delthatk[j,j] - Delthatk[j,-j] %*% invdd)
                         pinvp <- t(solve(Phi[!mj,!mj], Phi[!mj,mj]))
                         xhatk[mj,j] <- nu[mj] + pinvp %*% (xhatk[!mj,j] - nu[!mj])
                       }
                       
                       ind <- sum((oldx[na_idx] - xhatk[na_idx])^2)/sum(oldxc[na_idx]^2)
                       iter <- iter + 1
                     }
                     return(xhatk)
                   })
  
  if (trace) {
    like[like_idx] <- MNloglike_iPCA(Xs = x, Ms = Mhat, Sig = Sighat, Delts = Delthats, lamS = lamS, lamDs = lamDs, q = q, print = T)
    like_idx <- like_idx + 1
  }
  
  return(list(xhat = xhat, xh.init = xinit, M.init = Mhatinit,
              Sigma.init = Siginit, Delta.init = Deltinits,
              Sigma = Sighat, Delta = Delthats,M = Mhat,
              loglike = like, rmvn.ed1s = rmvn.ed1s))
}

# choose_lambdas: function to choose penalty parameters via ECM-type algorithm
#
# Inputs:
#  -dat = list of K data matrices (centered or uncentered)
#  -q = type of penalty, either 1, addfrob, or multfrob
#  -lams = (optional) vector of possible lams to use for construction of lam_grid
#  -lam_grid = grid of penalty parameters to try; each row is a set of penalty parameters to try;
#              lamS (if needed) must be in the first column, followed by lamDs
#  -prop_miss = proportion of data missing from each Xk (default = 0.05)
#  -trcma = T/F; whether or not to use TRCMAimpute_iPCA method (instead of RMNimpute_iPCA)
#  -maixt = maximum number of iterations for the imputation method
#  -greedy.search = logical T/F; whether or not to do a greedy search to choose penalty parameters
#  -maxit.search = number of iterations for the greedy search
#  -save.filepath = path where to save files (optional)
#  -save.fileext = extension name of files (optional)
#  -seed = numeric for set.seed (optional)
#  
# Outputs:
#  -best_lambdas = optimal set of penalty parameters
#  -errs = frobenius error from the scattered missing data
#  -errs_tab = data.frame with errors per dataset
#  -Xmiss = X with scattered missing elements
choose_lambdas <- function(dat, q, lams, lam_grid, prop_miss = 0.05, trcma = T, maxit = 20,
                           greedy.search = T, maxit.search = 1, save.filepath, save.fileext = "", seed = 100) {
  
  
  # if (!missing(seed)) {
  set.seed(seed)
  # }
  
  Xs <- dat
  
  # record dimensions
  n <- nrow(Xs[[1]])
  pks <- sapply(Xs, FUN = ncol)
  p <- sum(pks)
  K <- length(Xs)
  
  # leave out scattered missing elements in each of the Xks
  Xmiss <- Xs
  for (k in 1:K) {
    Xk <- Xs[[k]]
    miss_idx <- sample(1:length(Xk), prop_miss*length(Xk), replace = F)
    Xmiss[[k]][miss_idx] <- NA
  }
  
  if (!missing(save.filepath)) {
    saveRDS(Xmiss, file = paste0(save.filepath, "Xmiss", save.fileext, ".rds"))
    saveRDS(lam_grid, file = paste0(save.filepath, "lam_grid", save.fileext, ".rds"))
  }
  
  if (!missing(lam_grid)) {
    err <- rep(NA,nrow(lam_grid))
    errs <- matrix(NA, nrow = nrow(lam_grid), ncol = K)
    Xdiff <- list()
  }else {
    err <- NA
    errs <- matrix(NA, nrow = 1, ncol = K)
    Xdiff <- list()
  }
  
  if (!greedy.search) {
    if (!missing(lams)) {
      nreps <- ifelse(q == "multfrob", K, K+1)
      lam_grid <- expand.grid(rep(list(lams),nreps))
    }
    for (l in 1:nrow(lam_grid)) {
      print(l)
      lambdas <- as.matrix(lam_grid[l,])
      if (q == "multfrob" & ncol(lam_grid) == K) {
        lamS <- NULL
        lamDs <- lambdas
      }else if ((q == 1 | q == "1_off" | q == "addfrob") & ncol(lam_grid) == K+1) {
        lamS <- lambdas[1]
        lamDs <- lambdas[-1]
      }else {
        stop("q must be either 1, 1_off, addfrob, or multfrob, and check the number of penalty parameters.")
      }
      
      if (trcma) {
        if (l == 1) {
          ans <- TRCMAimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, seed = seed, maxit = maxit)
          rmvn.ed1s <- ans$rmvn.ed1s
        }else {
          ans <- TRCMAimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, rmvn.ed1s = rmvn.ed1s, seed = seed, maxit = maxit)
        }
      }else {
        ans <- RMNimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, seed = seed, maxit = maxit)
      }
      
      Xhat <- ans$xhat
      
      Xdiff[[l]] <- mapply(X = Xs, Y = Xhat, FUN = function(X,Y) {return(X-Y)},
                           SIMPLIFY = F)
      if (!missing(save.filepath)) {
        saveRDS(Xdiff, file = paste0(save.filepath, "Xdiff", save.fileext, ".rds"))
      }
      
      Xnum <- mapply(X = Xdiff[[l]], I = Xmiss, FUN = function(X, I) { return(X[is.na(I)]) }, SIMPLIFY = F) # numerator
      Xden <- mapply(X = Xs, I = Xmiss, 
                     FUN = function(X, I) {
                       Xc <- scale(X, center = T, scale = F)
                       return(Xc[is.na(I)]) 
                     }, SIMPLIFY = F) # denominator
      
      errs[l,] <- mapply(X = Xnum, Y = Xden,
                         FUN = function(X, Y) {
                           return(as.numeric((sum(X^2)^2) / (sum(Y^2)^2)))
                         })
      
      err[l] <- sum(errs[l,])
      
      if (!missing(save.filepath)) {
        saveRDS(errs, file = paste0(save.filepath, "lambda_errors_table", save.fileext, ".rds"))
        saveRDS(err, file = paste0(save.filepath, "lambda_errors", save.fileext, ".rds"))
      }
    }
    
    best_lambdas <- lam_grid[which.min(err),]
    errs_tab <- errs
    errs <- err
    
  }else { # do greedy search
    if (!missing(lam_grid)) {
      lams <- list()
      best_lambdas <- rep(NA,ncol(lam_grid))
      for (col in 1:ncol(lam_grid)) {
        lams[[col]] <- unique(lam_grid[,col]) # get unique penalty parameters
        best_lambdas[col] <- lams[[col]][1] # initialize
      }
    }else {
      nreps <- ifelse(q == "multfrob", K, K+1)
      lams <- rep(list(lams),nreps)
      best_lambdas <- rep(NA, nreps)
      for (col in 1:nreps) {
        best_lambdas[col] <- lams[[col]][1] # initialize
      }
    }
    
    # fix all but one penalty parameter
    search.iter <- 1; ptr <- 1
    already_searched <- NULL
    while (search.iter <= maxit.search) {
      old_best_lambdas <- best_lambdas
      for (idx in 1:length(best_lambdas)) {
        err <- rep(NA,length(lams[[idx]]))
        for (l in 1:length(lams[[idx]])) {
          print(paste0("Search.iter = ", search.iter, " // idx = ", idx, " // l = ",l))
          lambdas <- best_lambdas; lambdas[idx] <- lams[[idx]][l]
          
          if (q == "multfrob" & length(best_lambdas) == K) {
            lamS <- NULL
            lamDs <- lambdas
          }else if ((q == 1 | q == "1_off" | q == "addfrob") & length(best_lambdas) == K+1) {
            lamS <- lambdas[1]
            lamDs <- lambdas[-1]
          }else {
            stop("q must be either 1, addfrob, or multfrob, and check the number of penalty parameters.")
          }
          
          if (is.null(already_searched)) {
            marker <- T
          }else {
            marker <- sum(apply(X = as.matrix(already_searched), 1, FUN = function(X) {return(identical(X[-1], lambdas))}) == 0)
          }
          
          if (marker) { # if havent already looked at this choice of lambda
            if (trcma) {
              if (search.iter == 1 & idx == 1 & l == 1) {
                ans <- TRCMAimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, seed = seed, maxit = maxit)
                rmvn.ed1s <- ans$rmvn.ed1s
              }else {
                ans <- TRCMAimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, rmvn.ed1s = rmvn.ed1s, seed = seed, maxit = maxit)
              }
            }else {
              ans <- RMNimpute_iPCA(x = Xmiss, lamS = lamS, lamDs = lamDs, q = q, seed = seed, maxit = maxit)
            }
            
            Xhat <- ans$xhat
            
            Xdiff[[ptr]] <- mapply(X = Xs, Y = Xhat, FUN = function(X,Y) {return(X-Y)}, 
                                   SIMPLIFY = F)
            
            #             if (!missing(save.filepath)) {
            #               saveRDS(Xdiff, file = paste0(save.filepath, "Xdiff", save.fileext, ".rds"))
            #             }
            
            Xnum <- mapply(X = Xdiff[[ptr]], I = Xmiss, FUN = function(X, I) { return(X[is.na(I)]) }, SIMPLIFY = F) # numerator
            Xden <- mapply(X = Xs, I = Xmiss, 
                           FUN = function(X, I) { 
                             Xc <- scale(X, center = T, scale = F)
                             return(Xc[is.na(I)]) 
                           }, SIMPLIFY = F) # denominator
            current_errs <- mapply(X = Xnum, Y = Xden,
                                   FUN = function(X, Y) {
                                     return(as.numeric((sum(X^2)^2) / (sum(Y^2)^2)))
                                   })
            errs <- rbind(errs, current_errs)
            err[l] <- sum(current_errs)
            
            already_searched <- rbind(already_searched, c(err[l], lambdas))
            ptr <- ptr + 1
          }else {
            err_idx <- which(apply(X = as.matrix(already_searched), 1, FUN = function(X) {return(identical(X[-1], lambdas))}))
            err[l] <- already_searched[err_idx, 1]
          }
          
          #           if (!missing(save.filepath)) {
          #             saveRDS(errs, file = paste0(save.filepath, "lambda_errors_table", save.fileext, ".rds"))
          #           }
          
          if (l > 1) {
            if (err[l] > err[l-1]) {
              break
            }
          }
        }
        
        # update best lambdas
        best_lambdas[idx] <- lams[[idx]][which.min(err)]
        print(best_lambdas)
      }
      
      search.iter <- search.iter + 1
      
      if (identical(old_best_lambdas, best_lambdas)) { # no change
        search.iter <- maxit.search + 1 # exit while loop
      }
    }
    
    lam_grid <- already_searched[,-1]
    errs_tab <- errs
    errs <- rowSums(errs)
    
  }
  
  return(list(best_lambdas = best_lambdas, errs = errs, errs_tab = errs_tab, lam_grid = lam_grid, Xmiss = Xmiss))
}

# RMVNimpute: function imputes missing values for multivariate normal case
# 
# Inputs:
#  -x = n x p data matrix which arises from multivariate normal
#  -rho = penalty parameter
#  -cov.corr = T/F logical; whether or not to compute corrected covariance matrix
#  -q = penalty type; either 1, "addfrob", "multfrob"
#  -ed1 = saved eigenvalue decomposition for Eig(t(xc) %*% xc) to speed up computations (optional)
#  -maxit = max number of iterations
#  -thr.em = threshold for EM algorithm
#  -trace = T/F; whether or not to show loglikilihood stuff
#  -maxit.glasso = max number of iterations for GLASSO
#  -thr.glasso = threshold for GLASSO
#
# Outputs:
#  -xhat = imputed data matrix
#  -Sig = estimated Sigma (p x p matrix)
#  -mu = estimated mean
#  -loglike = loglikelihood value
#  -rmvn.ed1 = saved eigenvalue decomposition for Eig(t(xc) %*% xc)

RMVNimpute <- function(x,rho,cov.corr=TRUE,q,ed1,
                       maxit=10,thr.em=1e-4,trace=FALSE,maxit.glasso=50,thr.glasso=1e-2) {
  
  n <- nrow(x)
  p <- ncol(x)
  rmi <- (1:n)[apply(is.na(x),1,sum)>0]  # indices of rows with missing values
  
  
  if (n > 1) {muhat <- colMeans(x,na.rm=TRUE)} else{muhat <- rep(0,p)}
  
  xc <- t(t(x) - muhat) # center by column means
  xc[is.na(x)] <- 0 # initialize missing data to 0 (equivalent to imputing column means)
  
  # Initialize Sighat by solving regularization optimization problem
  iter <- 1
  Sighat <- t(xc)%*%xc
  if (q==1){ 
    cat(paste0('Glasso Iteration: ', iter, '\n'))
    quic_res <- QUIC(S = Sighat/n, rho = rho/n, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                     X.init = diag(p), W.init = diag(p))
    Sigihat <- quic_res$X
    Sighat <- quic_res$W
    ed1 <- NULL
  }else if (q == "1_off"){
    cat(paste0('Glasso Iteration: ', iter, '\n'))
    tmp_p <- nrow(Sighat)
    rho_off <- rho/n * (matrix(1, nrow = tmp_p, ncol = tmp_p) - diag(tmp_p))
    quic_res <- QUIC(S = Sighat/n, rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                     X.init = diag(p), W.init = diag(p))
    Sigihat <- quic_res$X
    Sighat <- quic_res$W
    ed1 <- NULL
  }else if (q == "addfrob"){
    if (missing(ed1)) {ed1 <- eigen(Sighat)}
    theta <- (ed1$values + sqrt(ed1$values^2 + 8*n*rho))/(2*n)
    Sighat <- ed1$vectors%*%diag(theta)%*%t(ed1$vectors)
  }else if (q == "multfrob") {
    if (missing(ed1)) {ed1 <- eigen(Sighat)}
    theta <- (ed1$values + sqrt(ed1$values^2 + 8*n*rho*n))/(2*n) # here, n = ||Sigi||_F^2 since we assume Sig = I
    Sighat <- ed1$vectors%*%diag(theta)%*%t(ed1$vectors)
  }else {
    stop("q must be 1, 1_off addfrob, or multfrob.")
  }
  
  # to ensure symmetry
  Sighat <- 1/2*(Sighat + t(Sighat))
  
  if(is.nan(sum(Sighat)) || is.na(sum(Sighat))) {maxit <- 0} # error checking
  
  like <- NA
  if(trace){
    like[iter] <- MVNloglike(x = x, mu = muhat, Sig = Sighat, rho = rho, q = q, print = T)
  }
  
  # repeat until EM converges
  ind <- 1
  xhat <- t(t(xc) + muhat) # undo centering
  na_idx <- is.na(x)
  while (ind > thr.em & iter < maxit) {
    oldxh <- xhat
    oldxc <- scale(oldxh, center = T, scale = F)
    
    iter <- iter + 1
    if (cov.corr) {covc <- matrix(0,p,p)}
    
    #E step
    for(i in rmi) {
      mi <- na_idx[i,]
      sinvs <- t(solve(Sighat[!mi,!mi], Sighat[!mi,mi]))
      xhat[i,mi] <- muhat[mi] + sinvs %*% (xhat[i,!mi] - muhat[!mi])
      if(cov.corr){ # compute corrected covariance matrix (see GA supplemental materials pg2)
        covc[mi,mi] <- covc[mi,mi] + Sighat[mi,mi] - sinvs %*% Sighat[!mi,mi]
      }          
    }
    
    #M Step
    if (n > 1) {muhat <- colMeans(xhat)} # update column means
    xhc <- t(t(xhat) - muhat) # center data
    if (cov.corr){ Sighat <- t(xhc)%*%xhc + covc}  else { Sighat <- t(xhc)%*%xhc }
    if(q == 1) {
      cat(paste0('Glasso Iteration: ', iter, '\n'))
      quic_res <- QUIC(S = Sighat/n, rho = rho/n, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigihat, W.init = Sighat)
      Sigihat <- quic_res$X
      Sighat <- quic_res$W
    }else if (q == "1_off"){
      cat(paste0('Glasso Iteration: ', iter, '\n'))
      tmp_p <- nrow(Sighat)
      rho_off <- rho/n * (matrix(1, nrow = tmp_p, ncol = tmp_p) - diag(tmp_p))
      quic_res <- QUIC(S = Sighat/n, rho = rho_off, tol = thr.glasso, msg = 0, maxIter = maxit.glasso, 
                       X.init = Sigihat, W.init = Sighat)
      Sigihat <- quic_res$X
      Sighat <- quic_res$W
    }else if (q == "addfrob"){
      ed <- eigen(Sighat)
      theta <- (ed$values + sqrt(ed$values^2 + 8*n*rho))/(2*n)
      Sighat <- ed$vectors%*%diag(theta)%*%t(ed$vectors)
    }else if (q == "multfrob") {
      ed <- eigen(Sighat)
      theta <- (ed$values + sqrt(ed$values^2 + 8*n*rho*n))/(2*n) # here, n = ||Sigi||_F^2 since we assume Sig = I
      Sighat <- ed$vectors%*%diag(theta)%*%t(ed$vectors)
    }else {
      stop("q must be 1, addfrob, or multfrob.")
    }
    
    # to ensure symmetry
    Sighat <- 1/2*(Sighat + t(Sighat))
    
    if(trace) {
      like[iter] <- MVNloglike(x = x, mu = muhat, Sig = Sighat, rho = rho, q = q, print = T) 
    }
    
    if (sum((oldxc[na_idx])^2, na.rm = T) == 0) {
      scaled_denom <- 1e-16 # to avoid dividing by 0
    }else {
      scaled_denom <- sum((oldxc[na_idx])^2,na.rm=TRUE)
    }
    
    ind <- sum((xhat[na_idx] - oldxh[na_idx])^2)/scaled_denom
    if(is.nan(sum(Sighat)) || is.na(sum(Sighat))){ind <- 0}
  }
  
  if(!trace){like=NULL}
  return(list(xhat=xhat,Sig=Sighat,mu=muhat,loglike=like,rmvn.ed1=ed1))
}


# Flip Flop algorithm 
# with multiplicative Frobenius estimators
FFmleMultFrob <- function(dat, lamDs,
                          maxit = 1e2, thr = 1e-6, init) {
  
  Xs <- dat
  
  if (missing(init)) {
    init <- c(list(diag(nrow(Xs[[1]]))), lapply(X = Xs, FUN = function(X) {return(diag(ncol(X)))}))
  }
  
  K <- length(Xs)
  n_cores <- min(K, detectCores()) # number of cores to use
  
  # record dimensions
  n <- nrow(Xs[[1]])
  pks <- sapply(Xs, FUN = ncol)
  p <- sum(pks)
  
  # center data
  for (k in 1:K) {
    Xs[[k]] <- scale(Xs[[k]], center = T, scale = F)
  }
  
  # initialize
  Sigi <- init[[1]]
  Deltis <- init[2:(K+1)]
  
  ind <- 1; iter <- 1; inds <- rep(0, each = maxit)
  while (ind > thr & iter < maxit) {
    oldS <- Sigi
    
    # update Sig
    inS <- Reduce("+", mcmapply(X = Xs, D = Deltis,
                                FUN = function(X, D) {return(X %*% D %*% t(X))},
                                SIMPLIFY = FALSE,
                                mc.cores = n_cores))
    inS <- 1/2 * (inS + t(inS)) # to ensure symmetry
    
    Sig_eigs <- eigen(inS)
    Sig_V <- Sig_eigs$vectors; Sig_d <- Sig_eigs$values
    
    sumDs <- sum(lamDs * sapply(X = Deltis, FUN = function(X) {return(norm(X, "F")^2)}))
    gams <- 1/(2*p) * (Sig_d + sqrt(Sig_d^2 + 8*p*sumDs))
    
    Sigi <- Sig_V %*% diag(1/gams) %*% t(Sig_V)
    Sig <- Sig_V %*% diag(gams) %*% t(Sig_V)
    
    # update Deltks
    inD <- mclapply(X = Xs, 
                    FUN = function(X) {return(t(X) %*% Sigi %*% X)},
                    mc.cores = n_cores)
    inD <- mclapply(X = inD,
                    FUN = function(X) {return(1/2 * (X + t(X)))},
                    mc.cores = n_cores) # to ensure symmetry
    
    Delt_eigs <- mclapply(X = inD, FUN = eigen, mc.cores = n_cores)
    Delt_Vs <- mclapply(X = Delt_eigs, FUN = function(X) {return(X$vectors)}, mc.cores = n_cores)
    Delt_ds <- mclapply(X = Delt_eigs, FUN = function(X) {return(X$values)}, mc.cores = n_cores)
    
    for (k in 1:K) {
      d <- Delt_ds[[k]]
      gams <- (1/(2*n)) * (d + sqrt(d^2 + 8*n*lamDs[k]*norm(Sigi, "F")^2))
      Deltis[[k]] <- Delt_Vs[[k]] %*% diag(1/gams) %*% t(Delt_Vs[[k]])
    }
    Deltis <- mclapply(X = Deltis, FUN = function(X) {return(1/2 * (X + t(X)))}, mc.cores = n_cores) # to ensure symmetry
    
    iter <- iter + 1
    ind <- sqrt(mean(lamDs))*norm(oldS - Sigi, type = "F") / norm(oldS, type = "F") # to account for differences in penalties
    inds[iter] <- ind
  }
  
  if (ind > thr) {
    message('Warning: Multiplicative Frobenius Sigma estimate did not converge.')
  }else {
    cat(paste0('Multiplicative Frobenius Estimate converged after ', iter-1, ' iterations. \n'))
  }
  
  # compute Delts
  Delts <- mclapply(X = Deltis, FUN = function(X) {chol2inv(chol(X))}, mc.cores = n_cores) # invert Deltis
  
  # to ensure symmetry
  Sig <- 1/2 * (Sig + t(Sig))
  Delts <- mclapply(X = Delts, FUN = function(X) {return(1/2 * (X + t(X)))}, mc.cores = n_cores)
  
  return(list(Sig = Sig, Delts = Delts, lamDs = lamDs))
}


plot_ipca_varexplained <- function(Xs, Sig, Delts,
                                   plot_title = "iPCA Proportion of Variance Explained", 
                                   show_legend = T, point_size = 1.5, line_size = 1, 
                                   show_plot = T, save_plot = T, mypath, file_ext = "") {
  
  mytheme <- theme(legend.text = element_text(family = "Helvetica", size = rel(1.25)),
                   strip.text = element_text(family = "Helvetica", size = rel(1.25), face = "bold"),
                   axis.title = element_text(family = "Helvetica", size = rel(1.25)), 
                   axis.text = element_text(family = "Helvetica", size = rel(1.25)), 
                   axis.line = element_line(size = 1,colour = "black"), 
                   axis.ticks = element_line(colour="black",size = rel(1)),
                   panel.grid.major = element_line(colour="grey90",size = rel(0.5)), 
                   panel.grid.minor = element_blank(), 
                   panel.background = element_rect(fill = "grey98"), 
                   legend.key = element_rect(fill = "grey98"), 
                   legend.title = element_text(family = "Helvetica", size = rel(1.25)), 
                   plot.title = element_text(face = "bold", size = rel(1.75),family = "Helvetica"))
  
  K <- length(Xs); n <- nrow(Xs[[1]]); pks <- sapply(Xs, FUN = ncol)
  
  # center the data
  Xs <- lapply(X = Xs, FUN = function(X) {return(scale(X, center = T, scale = F))})
  
  
  Sig_eig <- eigen(Sig)
  U <- Sig_eig$vectors
  d <- Sig_eig$values
  Uinv <- solve(U)
  
  props_U <- props_Uk <- matrix(0, nrow = n, ncol = K)
  for (k in 1:K) {
    # regular PCA
    X_svd <- svd(Xs[[k]])
    Uk <- X_svd$u; Dk <- diag(X_svd$d); Vk <- X_svd$v
    
    Deltk_V <- eigen(Delts[[k]])$vectors
    
    tv <- (norm(Xs[[k]], 'F'))^2 # total variance
    
    # plot proportion of variance of X explained per matrix
    for (l in 1:min(n,pks[k])) {
      
      proj_U <- t(U[,1:l]) %*% Xs[[k]] %*% Deltk_V[,1:l]
      var_U <- norm(proj_U, 'F')^2
      
      prop_var_U <- var_U / tv
      props_U[l,k] <- prop_var_U
      
      proj_Uk <- Uk[,1:l] %*% t(Uk[,1:l]) # usual PCA
      var_Uk <- norm(proj_Uk %*% Xs[[k]], 'F')^2
      prop_var_Uk <- var_Uk / tv
      props_Uk[l,k] <- prop_var_Uk
      
    }
  }
  
  props_U_long <- cbind(reshape2::melt(data.frame(m = 1:n, props_U), id = "m"), var_ex = "iPCA (Joint PCs)")
  props_Uk_long <- cbind(reshape2::melt(data.frame(m = 1:n, props_Uk), id = "m"), var_ex = "PCA (Individual PCs)")
  
  df_long <- rbind(props_U_long, props_Uk_long)
  df_long[df_long == 0] <- NA
  
  p1 <- ggplot(df_long) +
    aes(x = m, y = value, group = var_ex) +
    facet_grid(~variable) +
    geom_point(size = point_size, aes(color = var_ex)) +
    geom_line(size = line_size, aes(color = var_ex)) +
    labs(x = "Number of Components", y = "Proportion of Variance Explained", color = "") +
    labs(title = plot_title) +
    ylim(0-1e-5,1+1e-5) +
    mytheme
  
  if (show_legend == F) {
    p1 <- p1 + guides(color = F, shape = F)
  }
  
  if (show_plot) {
    #   annotate_figure(p1,
    #                   top = text_grob(plot_title, color = "black", face = "bold", size = 16))
    print(p1)
  }
  
  if (save_plot) {
    ggsave(paste0(mypath,"ipca_prop_var_explained", file_ext, ".pdf"), p1, width = 11, height = 6)
    saveRDS(p1, file = paste0(mypath,"ipca_prop_var_explained", file_ext, ".rds"))
  }
  
  return(list(ipca_varexplained_plot = p1, df = df_long))
}
plot_ipca <- function(Sig, U, y, plot_title, show_legend = F, point_size = 1, text_size = 12, show_plot = T, pcs) {
  
  mytheme <- theme(legend.text = element_text(family = "Helvetica", size = rel(1)), 
                   axis.title = element_text(family = "Helvetica", size = rel(1)), 
                   axis.text = element_text(family = "Helvetica", size = rel(1)), 
                   axis.line = element_line(size = 1,colour = "black"), 
                   axis.ticks = element_line(colour="black",size = rel(1)),
                   panel.grid.major = element_line(colour="grey90",size = rel(0.5)), 
                   panel.grid.minor = element_blank(), 
                   panel.background = element_rect(fill = "grey98"), 
                   legend.key = element_rect(fill = "grey98"), 
                   legend.title = element_text(family = "Helvetica", size = rel(1)), 
                   plot.title = element_text(face = "bold", size = rel(1.25),family = "Helvetica"))
  
  # myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  # sc <- scale_colour_gradientn(colours = myPalette(100)) # for continuous y labels
  # sc <- scale_colour_viridis(discrete = F, option = "plasma")
  sc <- scale_colour_viridis(discrete = F, option = "plasma", begin = 0, end = 0.95)
  
  if (!(missing(U))) {
    Sig_V <- U
  }else {
    Sig_eig <- eigen(Sig)
    Sig_V <- Sig_eig$vectors; Sig_d <- Sig_eig$values
  }
  
  if (!missing(pcs)) { # rearrange depending on desired pcs
    Sig_V <- cbind(Sig_V[,pcs], Sig_V[,-pcs]) 
  }
  
  if (!missing(y)) {
    p1 <- ggplot(data = data.frame(Sig_V, y = y)) +
      aes(x = X1, y = X2, color = y) +
      geom_point(size = point_size) +
      labs(x = "PC1", y = "PC2", color = "Class") +
      mytheme 
    if (!is.factor(y)) {
      p1 <- p1 + sc
    }
  }else {
    p1 <- ggplot(data = data.frame(Sig_V)) +
      aes(x = X1, y = X2) +
      geom_point(size = point_size) +
      labs(x = "PC1", y = "PC2", color = "Class") +
      mytheme 
  }
  
  if (!(missing(plot_title))) {
    p1 <- p1 + labs(title = plot_title)
  }
  
  if (show_legend == F) {
    p1 <- p1 + guides(color = F, shape = F)
  }
  
  if (show_plot) {
    print(p1)
  }
  
  return(list(ipca_plot = p1, ipca_scores = Sig_V))
}

