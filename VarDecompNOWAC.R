# load libraries 
# install.packages('r.jive')
library(r.jive)
# devtools::install_github("idc9/r_jive")
library(ajive)
# to install iPCA you need to download the folder
# or load the functions
# library(iPCA)
source('work/multiomics/erica/functions/functions_iPCA.R')

library(SpatioTemporal)
library(ggplot2)
library(sparseEigen)


# Load NOWAC data
# three omic layers: methylation, mRNA and miRNA

# this file contains:
# beta.LC: CpGs methylation
# exprs: mRNA expression
# LC samples: technical info
# covs: covariates
# keeps: vector of id numbers of patients to be included in the analysis
load(file = 'work/multiomics/erica/DNAmethylation_Geneexpr_Covariates_Lung.RData')

# this file contains miRNA expression data
miRNA <- readRDS(file="work/multiomics/erica/020320_Filtered_199miRNAs_NOWAC.rds")

# only patients with all omics collected
covs <- covs[rownames(covs) %in% rownames(miRNA),]
table(covs$case_ctrl, covs$pairs)
# needs complete pairs so remove singles manually
covs <- covs[!rownames(covs) %in% c("case 12", "case 16", "case 20",
                                    "case 24", "case 52", "case 66", 
                                    "ctrl 68", "case 77",
                                    "case 86",  "case 105"),]
table(covs$case_ctrl, covs$pairs)

# only patients with all omics collected
miRNA <- miRNA[rownames(miRNA) %in% rownames(covs),]
exprs <- exprs[,colnames(exprs) %in% rownames(miRNA)]
beta.LC <- beta.LC[,colnames(beta.LC) %in% rownames(miRNA)]
LCsamples <- LCsamples[rownames(LCsamples) %in% rownames(miRNA),]

beta.LC <- beta.LC[, rownames(miRNA)]
exprs <- exprs[, rownames(miRNA)]
covs <- covs[rownames(miRNA),]
LCsamples <- LCsamples[rownames(miRNA),]

# optional: transform beta values into M values 
# filtering can then be done on M instead of beta
Mvalues <- log2(beta.LC/(1-beta.LC))


source('work/multiomics/erica/functions/filtering.R')

# Filtering to the top 50000 most variable CpG probes based on variance
# alternative could be IQR
DNAm.var50000 <- filtering(Mvalues, var, 50000)
# same procedure for mRNA expression (select first 5000)
mRNA.var5000 <- filtering(exprs, var, 5000)

# keep all miRNA
# but make it numeric
miRNA <- lapply(miRNA, function(l) as.numeric(l))


methylation <- as.data.frame(t(DNAm.var50000))
mRNA <- as.data.frame(t(mRNA.var5000))
miRNA <- as.data.frame(miRNA)
rownames(methylation) <- rownames(mRNA) <- rownames(miRNA) <- 1:230
dataNOWAC <- list(t(methylation), t(mRNA), t(miRNA))



############################################
### JIVE
# implement jive on the three data sets
# also store computing time 
# we only need the jive function from the package
# it already performs the SVDmiss for missing data imputation
# it centers and scales the data
t1 <- system.time(jiveResults <- jive(dataNOWAC))
# visualize results
# very simple plot of proportions of variance
showVarExplained(jiveResults)
# this gives an overview
summary(jiveResults)
# PCA on first two components
# joint and individual
showPCA(jiveResults, n_joint = 2, n_indiv = c(2,2,2))

# extract loadings estimated from JIVE
# and top variables for higher loadings
# in terms of both joint and individual
loadings <- extract_loadings_jive(data, jiveResults)



# repeat jive implementation on possible combinations of two sets
# also store computing time 
data12 <- list(t(methylation), t(mRNA))
t12 <- system.time(jiveResults12 <- jive(data12))
# visualize results
# very simple plot of proportions of variance
showVarExplained(jiveResults12)
# this gives an overview
summary(jiveResults12)
# PCA on first two components
# joint and individual
showPCA(jiveResults12, n_joint = 2, n_indiv = c(2,2,2))

# extract loadings estimated from JIVE
# and top variables for higher loadings
# in terms of both joint and individual
loadings12 <- extract_loadings_jive(data12, jiveResults12)


data13 <- list(t(methylation), t(miRNA))
t13 <- system.time(jiveResults13 <- jive(data13))
# visualize results
# very simple plot of proportions of variance
showVarExplained(jiveResults13)
# this gives an overview
summary(jiveResults13)
# PCA on first two components
# joint and individual
showPCA(jiveResults13, n_joint = 2, n_indiv = c(2,2,2))

# extract loadings estimated from JIVE
# and top variables for higher loadings
# in terms of both joint and individual
loadings13 <- extract_loadings_jive(data13, jiveResults13)

 
data23 <- list(t(mRNA), t(miRNA))
t23 <- system.time(jiveResults23 <- jive(data23))
# visualize results
# very simple plot of proportions of variance
showVarExplained(jiveResults23)
# this gives an overview
summary(jiveResults23)
# PCA on first two components
# joint and individual
showPCA(jiveResults23, n_joint = 2, n_indiv = c(2,2,2))

# extract loadings estimated from JIVE
# and top variables for higher loadings
# in terms of both joint and individual
loadings23 <- extract_loadings_jive(data23, jiveResults23)


############################################
### AJIVE
# implement angle based JIVE

# use function to prepare data for ajive

data.ajive <- ajive.dataprep(dataNOWAC)
colnames(data.ajive[[1]]) <- rownames(dataNOWAC[[1]])
colnames(data.ajive[[2]]) <- rownames(dataNOWAC[[2]])
colnames(data.ajive[[3]]) <- rownames(dataNOWAC[[3]])


# need to determine initial ranks
# look at screeplots
screeplot(data.ajive)

# run ajive
ajiveResults <- ajive(data.ajive, 
                      initial_signal_ranks = c(39, 18,11)) 
# variation proportions as in jive
showVarExplained.ajive(ajiveResults, data.ajive)

# extract loadings 
# and top variables for higher loadings
# in terms of both joint and individual
loadings.ajive <- extract_loadings_ajive(data.ajive, ajiveResults)


############################################
### iPCA
# implement integrated PCA
# set initial lambdas
lams <- c(1e-4, 1e-2, 1, 100, 10000)

# it needs the transpose wrt jive 
data.ipca <- lapply(dataNOWAC, function(x) t(x))

# select lambdas
best_lambdas <- choose_lambdas(dat = data.ipca, q = "multfrob", 
                               lams = lams, greedy.search = T)$best_lambdas

# run Flip FLop Algorithm 
iPCAresults <- FFmleMultFrob(dat = data.ipca, 
                          lamDs = best_lambdas)


# visualize results 
plot_ipca(Sig = iPCAresults$Sig, pcs = 1:2)$ipca_plot
plot_ipca(Sig = iPCAresults$Sig, y = covs$case_ctrl_factor, pcs = 1:2)$ipca_plot

# visualize scores
scores_iPCA <- plot_ipca(Sig = iPCAresults$Sig, pcs = 1:2)$ipca_scores
# variation proportions as in jive
plot_ipca_varexplained(Xs = as.matrix(data.ipca), 
                       Sig = iPCAresults$Sig, 
                       Delts = iPCAresults$Delts)

# to extract top variables apply sparse PCA to each 
# covariance matrix (delta) estimated by iPCA
# use sparse pcan on each Delts
top.eigen <- lapply(iPCAresults$Delts, function(x) spEigen(x, 10, 0.6))
pply(iPCAresults$Delts, function(x) spEigen(x, 10, 0.6))
.6))
sults$Delts, function(x) spEigen(x, 10, 0.6))
