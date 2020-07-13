gam_missing <- function(Est_target, Est_reference,i){
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' n the degrees of freedom used as a smoothing parameter, default 4 (cubic spline)
  Est_target <- as.numeric(Est_target)
  Est_reference <- as.matrix(Est_reference)
  checks <- makeAssertCollection()
  assert_matrix(x = Est_reference,
                mode = 'numeric',
                all.missing = FALSE, add = checks)
  assert_integerish(x = i, lower = 1,
                    upper = dim(Est_reference)[1], add = checks)
  assert_numeric(x = Est_target,
                 len = dim(Est_reference)[1],
                 finite = TRUE, all.missing = FALSE, add = checks)
  reportAssertions(checks)


  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  data <- data.frame(Est_target, Est_reference)
  colnames(data) <- c('Target',paste('c',1:(ncol(data)-1),sep=''))
  nom <- paste('c',1:(ncol(data)-1),sep='')
  formulae <- paste(nom,collapse = ' + ')
  formulae <- paste('Target ~',formulae,collapse = ' ')

  gg <- gam(formula(formulae),
            data = data,
            na.action =na.gam.replace)

  val <- predict(gg,data[i,])
  return(val)
}
