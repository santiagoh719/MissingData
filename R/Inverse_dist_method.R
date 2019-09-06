IDW <- function( Est_reference, i, k, cord = T , distans = F){
  #' Inverse distance weigth method
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Est_reference must have the same nrows
  #' Also the same start date.
  #' Either cord is a matrix of two colums, first Longitud, second Latitud
  #' The first row is the cordinate of target station, then by order the reference stations
  #' Or dist, a vector with the distance to Est_target, in order by reference stations
  #' k the potencial factor
  Est_reference <- as.matrix(Est_reference)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}

  if(!is.logical(cord)){
    if(ncol(cord)!= 2){stop('not acceptable dimentions of cord, must have 2 columns (Long, Lat)')}
    if(nrow(cord) != ncol(Est_reference)+1){stop('not coerent dimentions of distans width Est_reference')}
    distans <-  spDistsN1(as.matrix(cord[2:nrow(cord),]),as.matrix(cord[1,]),longlat = T)
  }else{
    if(length(distans) != ncol(Est_reference)){stop('not coerent dimentions of distans width Est_reference')}
  }

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  reference <- data.frame(Est_reference[,ind])
  distans <- distans[ind]
  distans <- distans**(-k)
  w <- distans / sum(distans)
  return( sum(reference[i,]*w,na.rm = T) )
}

IDW_simple <- function( Est_reference, i, cord = T , distans = F){
  #' Inverse distance weigth method simple
  #' Wei y McGuinnes 1973
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Est_reference must have the same nrows
  #' Also the same start date.
  #' Either cord is a matrix of two colums, first Longitud, second Latitud
  #' The first row is the cordinate of target station, then by order the reference stations
  #' Or dist, a vector with the distance to Est_target, in order by reference stations
  k = 1
  return(IDW( Est_reference, i, k, cord , distans))
}

RNNWM <- function(Est_target, Est_reference, i, k, reg= 1){
  #' Revised nearest neighbor weighting method
  #' select the station acording similarity paterns and create a new distans
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' k the potencial factor
  #' reg is the amount of values to use to create the paterns, 1 and 2 are possible
  #' make_exp, should IDW be exponential, default FALSE
  Est_reference <- as.matrix(Est_reference)
  Est_target <- as.vector(Est_target)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(nrow(Est_reference) != length(Est_target)){stop('Not coerent dimentions of Est_reference and Est_target')}

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  distans <- patern_distans_change(Est_target, Est_reference,reg)
  return(IDW(Est_reference, i, k, distans = distans))
}
