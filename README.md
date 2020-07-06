# VarianceDecompositionNOWAC

A set of scripts and auxiliary functions to perform aJIVE and iPCA on different omic layers. 

Main R scripts:

- datafiltering.R : filters the omic data. 
mRNA is filtered based on variance (5000), 
CpGs are selected from the same genes as mRNA, missing values >40% eliminated, extreme M values (>3) eliminated
miRNA all included
needs filtering.R

-aJIVE.R: runs ajive on the filtered data. 
ranks chosen by profile likelihood
needs ajive.dataprep.R, show.var.explained.ajive and aJIVE package
screeplots also plotted (needs Screeplots.R)

-iPCA.R: runs iPCA on a smaller subset of data 
(same as data.filtering but 500 mRNA to start with)
runs aJIVE on the same subset
needs function_iPCA

- insample_predictions_ROC: predictions and ROC curve in sample for 
case vs control and metastasis on the aJIVE results

- predictions_ipca.R: predicts case vs control and metastasis
 on the same subset as iPCA. both in sample and 10 fold CV
needs showPCAajive

- predictions_CV.R: cross validation for predictions of case vs control 
and metastasis on the aJIVE results

- predictions_rf.R: random forests for predictions of case vs control 
 on the aJIVE results

- singlePCA.R: performs non integrative analysis on single data sets
