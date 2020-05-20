UK <-  function(Est_target, Est_reference, i){
  #' UK traditional method
  #' Hasanpur Kashani & Dinpashoh 2012
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Both Est_target and Est_reference must have the same nrows
  #' Also the same start date.

  Est_reference <- as.matrix(Est_reference)
  Est_target <- as.numeric(Est_target)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(nrow(Est_reference) != length(Est_target)){stop('Not coerent dimentions of Est_reference and Est_target')}

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  ind <- apply(matrix(as.numeric(!is.na(Est_reference)),ncol=ncol(Est_reference)),2,sum) > 4
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  if(sum(as.numeric(!is.na(Est_target))) < 4){
    stop('Not enougth stations')
  }

  Est_reference <- Est_reference[,ind]
  mat_cor <- cor(Est_target,Est_reference,use = 'pairwise.complete.obs')
  if(all(is.na(mat_cor))){
    stop('use another method or use more reference stations')
  }
  est <- which.max(mat_cor)

  return(Est_reference[i,est] * mean(Est_target,na.rm = T)/mean(Est_reference[,est],na.rm = T) )
}

NR <- function(Est_target, Est_reference, i){
  #' Normal Ratio method
  #' Proposed by Paulhus y Kohler 1952
  #' Modified by Young 1992
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Both Est_target and Est_reference must have the same nrows
  #' Also the same start date.
  Est_reference <- as.matrix(Est_reference)
  Est_target <- as.numeric(Est_target)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(nrow(Est_reference) != length(Est_target)){stop('Not coerent dimentions of Est_reference and Est_target')}

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  r <- cor(Est_target,Est_reference,use = 'pairwise.complete.obs')
  Est_reference <- Est_reference[,r>0]
  r <- r[r>0]
  n <- apply(matrix(as.numeric(!is.na(Est_reference)),nrow = nrow(Est_reference)), 2, sum)
  w <- r**2 * (n-2) / (1-r**2)
  return(sum(w*Est_reference[i,],na.rm = T)/sum(w))
}

NR_1952 <- function(Est_target, Est_reference, i, tt='month'){
  #' Normal Ratio method
  #' Proposed by Paulhus y Kohler 1952
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Both Est_target and Est_reference must have the same nrows
  #' Also the same start date.
  #' tt is periodicity of sample data, if tt=day is dayly data
  #' if tt = month, is monthly data (default)
  #' if tt = hour, is hourly data
  #' if tt = 6hour, is data sample every 6 hours

  Est_reference <- as.matrix(Est_reference)
  Est_target <- as.numeric(Est_target)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(nrow(Est_reference) != length(Est_target)){stop('Not coerent dimentions of Est_reference and Est_target')}

  n <- length(Est_target)

  if(tt=='day'){
    naux <- 365
  }else if(tt=='month'){
    naux <- 12
  }else if(tt=='hour'){
    naux <- 365*24
  }else if(tt=='6hour'){
    naux <- 365*4
  } else{ stop('No aceptable value of tt given, use either month, day, hour or 6hour')}

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  iters <- ncol(Est_reference)
  pp_reference <- vector(mode='double',length = iters)
  for(j in 1:iters){
    pp_reference[j] <- mean(rollapply(Est_reference[,j],width=naux,by=naux,FUN = sum ),na.rm = T)
  }
  pp_target <- mean(rollapply(Est_target,width=naux,by=naux,FUN = sum ),na.rm = T)

  w <- pp_target / pp_reference

  return(mean(w*Est_reference[i,]))

}

CWM <- function(Est_target, Est_reference, i){
  #' Correlation weighting method
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Both Est_target and Est_reference must have the same nrows
  #' Also the same start date.
  Est_reference <- as.matrix(Est_reference)
  Est_target <- as.numeric(Est_target)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(nrow(Est_reference) != length(Est_target)){stop('Not coerent dimentions of Est_reference and Est_target')}

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  w <- cor(Est_target,Est_reference,use = 'pairwise.complete.obs')
  Est_reference <- Est_reference[,w>0]
  w <- w[w>0]

  return(sum(w*Est_reference[i,])/sum(w))
}

