# VarianceDecompositionNOWAC

A set of functions to perform aJIVE and iPCA on different omic layers. 

Main R scripts: 
- datafiltering: prepares the data in the format required by aJIVE and iPCA. Filtering is based on variance for mRNAs (first 5000 for aJIVE, first 500 for iPCA), CpGs on the same genes of the selected mRNAs are then selected and filtered. miRNAs are all included.
- aJIVE
- iPCA
- predictions_cv: prediction models using joint scores from aJIVE on case vs control and metastasis. 10 fold cross validation
- predictions_rf: prediction models using joint scores from aJIVE on case vs control and metastasis. Random forests with 500 trees. Variable importance plots and random forest diagnostics.
