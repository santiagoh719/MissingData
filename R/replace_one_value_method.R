
AA <- function(Est_reference, i){
  #' Simple aritmetric average
  #' given for the same time
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by stemated  Est_reference <- as.matrix(Est_reference)
  Est_reference <- as.matrix(Est_reference)
  checks <- makeAssertCollection()
  assert_matrix(x = Est_reference,
                mode = 'numeric',
                all.missing = FALSE, add = checks)
  assert_integerish(x = i, lower = 1,
                    upper = dim(Est_reference)[1], add = checks)
  reportAssertions(checks)

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop(paste0('Not enougth stations with data for the ', i ,'time'))
  }
  Est_reference <- as.matrix(Est_reference[,ind])
  return(mean(Est_reference[i,],na.rm = T))
}

AM <- function(Est_reference, i){
  #' median of reference stations
  #' given for the same time
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by stemated
  Est_reference <- as.matrix(Est_reference)
  checks <- makeAssertCollection()
  assert_matrix(x = Est_reference,
                mode = 'numeric',
                all.missing = FALSE, add = checks)
  assert_integerish(x = i, lower = 1,
                    upper = dim(Est_reference)[1], add = checks)
  reportAssertions(checks)

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop(paste0('Not enougth stations with data for the ', i ,'time'))
  }
  Est_reference <- as.matrix(Est_reference[,ind])
  return(quantile(Est_reference[i,],0.5,na.rm = T))
}

NN <- function(Est_reference, i, cord = T, distans = F){
  #' Nearest Nightbore, value of nearest station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by stemated
  #' Either cord is a matrix of two colums, first Longitud, second Latitud
  #' The first row is the cordinate of target station, then by order the reference stations
  #' Or dist, a vector with the distance to Est_target, in order by reference stations

  Est_reference <- as.matrix(Est_reference)
  checks <- makeAssertCollection()
  assert_matrix(x = Est_reference,
                mode = 'numeric',
                all.missing = FALSE, add = checks)
  assert_integerish(x = i, lower = 1,
                    upper = dim(Est_reference)[1], add = checks)
  reportAssertions(checks)

  if(!is.logical(cord)){
    if(ncol(cord)!= 2){
      stop('not acceptable dimentions of cord, must have 2 columns (Long, Lat)')
    }
    if(nrow(cord) != ncol(Est_reference)+1){
      stop('not coherent dimentions of cord width Est_reference')
    }
    distans <-  spDistsN1(as.matrix(cord[2:nrow(cord),]),
                          as.matrix(cord[1,]),
                          longlat = T)
  }else{
    if(length(distans) != ncol(Est_reference)){
      stop('not coerent dimentions of distans width Est_reference')
    }
  }

  Est_reference <- as.matrix(Est_reference)
  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop(paste0('Not enougth stations with data for the ', i ,'time'))
  }

  Est_reference <- as.matrix(Est_reference[,ind])
  distans <- distans[ind]
  return(Est_reference[i,which.min(distans[])])
}

SBE <-  function(Est_target, Est_reference, i){
  #' Single best estimator method, value of most correlated station
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Both Est_target and Est_reference must have the same nrows
  #' Also the same start date.
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

  if(sum(ind) == 1){
    return(Est_reference[i])
  }
  if(all(!(ind))){
    stop(paste0('Not enougth stations with data for the ', i ,'time'))
  }
  Est_reference <- Est_reference[,ind]

  if(sd(Est_target,na.rm = T)<0.0001){
    stop('Target station with almost zero standard deviation, computation fail')
  }

  cols <- apply(Est_reference,2,sd,na.rm=T) > 0.0001

  if(all(!(cols))){
    stop(paste0('Not enougth stations with data for the ', i ,'time'))
  }

  Est_reference <- Est_reference[,cols]
  if(sum(cols) == 1){
    return(Est_reference[i])
  }

  est <- which.max(cor(Est_target,Est_reference,use = 'pairwise.complete.obs'))
  return(Est_reference[i,est])
}

