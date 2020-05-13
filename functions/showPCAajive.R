showPCA.ajive <- function (results, n, n_joint = 0, 
          Colors = "black", pch = 1) 
{
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  nPCs = n_joint 
  PCs = matrix(nrow = nPCs, ncol = n)
  PC_names = rep("", nPCs)
  if (n_joint > 0) {
    
    PCs[1:n_joint, ] = results$joint_scores
    PC_names[1:n_joint] = paste("Joint ", 1:n_joint)
  }
  
  nplots = (nPCs - 1)^2
  par(mar = c(4, 4, 2, 2))
  layout(matrix(c(1:nplots), nrow = nPCs - 1, ncol = nPCs - 
                  1))
  for (i in 2:nPCs) {
    for (j in 1:(nPCs - 1)) {
      if (j >= i) 
        plot.new()
      else plot(PCs[i, ], PCs[j, ], xlab = PC_names[i], 
                ylab = PC_names[j], col = Colors, pch = pch)
    }
  }
}

