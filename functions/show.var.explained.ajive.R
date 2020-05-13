# show variance explained from ajive results 
showVarExplained.ajive <- function(ajiveResults, blocks,
                                   col = c("grey20", "grey43", "grey65")){
l <- length(blocks)
# joint variance
VarJoint = rep(0, l)
for (i in 1:l) VarJoint[i] = norm(as.matrix(ajiveResults$block_decomps[[i]]$joint[[4]]), 
                                  type = "F")^2/norm(blocks[[i]], type = "F")^2

# individual variances
VarIndiv = rep(0, l)
for (i in 1:l) VarIndiv[i] = norm(as.matrix(ajiveResults$block_decomps[[i]]$individual[[4]]), 
                                  type = "F")^2/norm(blocks[[i]], type = "F")^2

# residual variance
VarSubtr = 1 - VarJoint - VarIndiv
# plot
par(mar = c(5.1, 4.1, 4.1, 0))
layout(matrix(c(1, 2), 1, 2), heights = c(5, 5), widths = c(5, 
                                                            2))

datavar <- data.frame(var = c(VarJoint, VarIndiv, VarSubtr), 
                      source = c(rep(c('methylation', 'mRNA','miRNA'),3)), 
                      type = c(rep('Joint',3), rep('Individual', 3), rep('Residual', 3)))
ggplot(datavar, aes(x = source, y = var, fill = type)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  labs(title="Variation Explained")+ scale_fill_brewer(palette="Blues")


VarProp <- list(VarJoint, VarIndiv, VarSubtr)
names(VarProp) <- c('Joint', 'Indiv', 'Resid')
return(VarProp)
}
