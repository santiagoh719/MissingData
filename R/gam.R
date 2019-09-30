gam_missing <- function(Est_target, Est_reference,i, n = 4, cord = NULL, dist = NULL){
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' n the degrees of freedom used as a smoothing parameter, default 4 (cubic spline)

  Est_reference <- as.matrix(Est_reference)
  Est_target <- as.numeric(Est_target)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(nrow(Est_reference) != length(Est_target)){stop('Not coerent dimentions of Est_reference and Est_target')}

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  data <- data.frame(Est_target, Est_reference)
  colnames(data) <- c('Target',paste('c',1:ncol(Est_reference),sep=''))
  nom <- paste('c',1:ncol(Est_reference),sep='')
  formulae <- paste(nom,collapse = ' + ')
  formulae <- paste('Target ~',formulae,collapse = ' ')
  gg <- gam(formula(formulae),data = data, na.action =na.gam.replace)
  val <- predict(gg,data[i,])
  return(val)
}
