# show variance explained from ajive results 
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
par(mar = c(5.1, 4.1, 4.1, 0))
layout(matrix(c(1, 2), 1, 2), heights = c(5, 5), widths = c(5, 
                                                            2))
barplot(rbind(VarJoint, VarIndiv, VarSubtr), col = col, main = "Variation Explained", 
        names.arg = names(blocks))
par(mar = c(0, 0, 0, 0))
plot.new()
legend(x = 0.05, y = 0.8, legend = c("Joint", "Individual", 
                                     "Residual"), bty = "n", fill = col)

VarProp <- list(VarJoint, VarIndiv, VarSubtr)
names(VarProp) <- c('Joint', 'Indiv', 'Resid')
return(VarProp)
}
