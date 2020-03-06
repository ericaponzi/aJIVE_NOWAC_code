# load libraries 
# install.packages('r.jive')
library(r.jive)
# devtools::install_github("idc9/r_jive")
library(ajive)
# to install iPCA you need to download the folder
# or load the functions
# library(iPCA)
source('functions/functions_iPCA.R')
library(SpatioTemporal)
library(ggplot2)
library(sparseEigen)
#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
require("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#BiocManager::install("lumiHumanIDMapping")
library("lumi")

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
miRNA <- readRDS(file="work/multiomics/020320_Filtered_199miRNAs_NOWAC.rds")

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

# filter gene expression first
source('work/multiomics/erica/functions/filtering.R')
mRNA.var5000 <- filtering(exprs, var, 500)
exprs.var <- rownames(mRNA.var5000)


# pick CpGs that are on the same genes as mRNA 

# need CpGs-genes association
## Gene names is listed in the UCSC_RefGene_Name column of the file
## CpGs are listed in Name
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annot450k = as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

# information of mRNA genes
# gene names in 3rd column
mRNAmapping <- nuID2RefSeqID(exprs.var, lib.mapping='lumiHumanIDMapping', returnAllInfo=TRUE)

# extract CpGs on these genes
genes <- intersect(annot450k$UCSC_RefGene_Name, mRNAmapping[,3])
CpGs <- annot450k[is.element(annot450k$UCSC_RefGene_Name, genes),  ]$Name

# filter beta.LC based on these CpGs
selected.CpGs <- intersect(CpGs, rownames(beta.LC))
DNAm.selected <- beta.LC[selected.CpGs, ]
naprop <- apply(is.na(DNAm.selected), 1, sum)
naprop <- naprop/230
DNAm.selected <- DNAm.selected[naprop < 0.4, ]
# transform into M values 
DNAm.selected <- log2(DNAm.selected/(1-DNAm.selected))

miRNA <- lapply(miRNA, function(l) as.numeric(l))


methylation <- as.data.frame(t(DNAm.selected))
mRNA <- as.data.frame(t(mRNA.var5000))
mRNA <- as.data.frame(t(mRNA.var5000))
miRNA <- as.data.frame(miRNA)
miRNA <- miRNA[, -1]
miRNA <- log(miRNA+1)
rownames(methylation) <- rownames(mRNA) <- rownames(miRNA) <- 1:230
dataNOWAC <- list(t(methylation), t(mRNA), t(miRNA))

t1 <- system.time(jiveResults <- jive(dataNOWAC))
# visualize results
# very simple plot of proportions of variance
showVarExplained(jiveResults)
# this gives an overview
summary(jiveResults)
save(jiveResults, file = 'jiveResultsCpGassocMlogmirna.RData')


data12 <- list(t(methylation), t(mRNA))

t12 <- system.time(jiveResults12 <- jive(data12))
# visualize results
# very simple plot of proportions of variance
showVarExplained(jiveResults12)
# this gives an overview
summary(jiveResults12)
# PCA on first two components
# joint and individual

showPCA(jiveResults12,  n_indiv = c(2,2))
save(jiveResults12, file = 'work/multiomics/erica/results/jiveResults12assocM.RData')

methylation <- as.data.frame(t(DNAm.selected))
mRNA <- as.data.frame(t(mRNA.var5000))
mRNA <- as.data.frame(t(mRNA.var5000))
miRNA <- as.data.frame(miRNA)
miRNA <- miRNA[, -1]
miRNA <- log2(miRNA+1)
rownames(methylation) <- rownames(mRNA) <- rownames(miRNA) <- 1:230
dataNOWAC <- list(t(methylation), t(mRNA), t(miRNA))

t1 <- system.time(jiveResults <- jive(dataNOWAC))
# visualize results
# very simple plot of proportions of variance
showVarExplained(jiveResults)
# this gives an overview
summary(jiveResults)
save(jiveResults, file = 'jiveResultsCpGassocMlogmirna.RData')


