# ajive_NOWAC

A set of scripts and auxiliary functions to perform aJIVE and iPCA on different omic layers. 

Main R scripts:

- datafiltering.R : filters the omic data. 
mRNA is filtered based on variance (5000), 
CpGs are selected from the same genes as mRNA, missing values >40% eliminated, extreme M values (>3) eliminated
miRNA all included
needs filtering.R

- aJIVE.R: runs ajive on the filtered data. 
ranks chosen by profile likelihood
needs ajive.dataprep.R, show.var.explained.ajive and aJIVE package
screeplots also plotted (needs Screeplots.R)

- predictions.R: in-sample, CV and random forests for predictions of case vs control and metastasis 
 on the aJIVE results
