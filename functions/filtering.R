# function to filter data 
# data is the initial complete dataset to be filtered
# crit can be variance, IQR or other quantities
# keep is the number of variables to keep in the filtered data
# returns filtered data with keep n of variables


filtering <- function(data, crit, keep){
  filtering.crit <- apply(data, 1, crit, na.rm = TRUE)
  topvar <- order(filtering.crit, decreasing = TRUE)
  data.filtered <- data[topvar[1:keep],] 
  rownames(data.filtered) <- rownames(data[topvar[1:keep],])
  return(data.filtered)
}