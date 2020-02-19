# screeplot function

screeplot <- function(data, title=''){
  for(l in 1:length(data)){
  singular.values <- svd(data[[l]])[['d']]
  index <- 1:length(singular.values)
  print(ggplot(as.data.frame(x=cbind(index, singular.values)),
         aes(x=index, y=singular.values)) +
    geom_point() +
    geom_line() +
    labs(y='singular value', x='index', title=title))
  }
  
}
