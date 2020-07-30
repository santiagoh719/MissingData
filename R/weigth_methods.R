UK <-  function(Est_target, Est_reference, i, W = NULL){
  #' UK traditional method
  #' Hasanpur Kashani & Dinpashoh 2012
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
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  if(is.null(W)){
    if(sum(ind) == 1){
      warning('only one station used')
      return(Est_reference[i] * mean(Est_target,na.rm = T) /
               mean(Est_reference,na.rm = T) )
    }

    mat_cor <- cor(Est_target,Est_reference,use = 'pairwise.complete.obs')
    if(all(is.na(mat_cor))){
      stop('use another method or use more reference stations')
    }
    est <- which.max(mat_cor)
    return(Est_reference[i,est] * mean(Est_target,na.rm = T) /
             mean(Est_reference[,est],na.rm = T) )
  }else{
    if(length(W) != length(ind) | !is.numeric(W)){
      stop('If W is given, it must be numeric and of the same length to the amount of reference stations')
    }

    W <- W[ind]
    if(sum(ind) == 1){
      warning('only one station used')
      return(Est_reference[i] * W )
    }

    mat_cor <- cor(Est_target,Est_reference,use = 'pairwise.complete.obs')
    if(all(is.na(mat_cor))){
      stop('use another method or use more reference stations')
    }

    est <- which.max(mat_cor)
    return(Est_reference[i,est] * W[est])
  }


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

  if(sum(ind) == 1){
    warning('only one station used')
    return(Est_reference[i])
  }

  r <- cor(Est_target,Est_reference,use = 'pairwise.complete.obs')
  Est_reference <- Est_reference[,r>0]
  r <- r[r>0]
  n <- apply(matrix(as.numeric(!is.na(Est_reference)),
                    nrow = nrow(Est_reference)), 2, sum)
  w <- r**2 * (n-2) / (1-r**2)
  return(sum(w*Est_reference[i,])/sum(w))
}

NR_1952 <- function(Est_target, Est_reference, i, tt='month', W = NULL){
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
  assert_character(x = tt, len = 1)
  reportAssertions(checks)
  n <- length(Est_target)

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  if(is.null(W)){
    if(tt=='day'){
      naux <- 365
    }else if(tt=='month'){
      naux <- 12
    }else if(tt=='hour'){
      naux <- 365*24
    }else if(tt=='6hour'){
      naux <- 365*4
    } else{
      stop('No aceptable value of tt given, use either month, day, hour or 6hour')
    }
    pp_target <- mean(rollapply(Est_target,
                                width=naux,by=naux,FUN = sum ),na.rm = T)

    if(sum(ind) == 1){
      warning('only one station used')
      pp_reference <- mean(rollapply(Est_reference,
                                     width=naux,by=naux,FUN = sum ),na.rm = T)
      return(Est_reference[i] * pp_target / pp_reference)
    }

    iters <- ncol(Est_reference)
    pp_reference <- vector(mode='double',length = iters)
    for(j in 1:iters){
      pp_reference[j] <- mean(rollapply(Est_reference[,j],
                                        width=naux,by=naux,FUN = sum ),na.rm = T)
    }

    W <- pp_target / pp_reference
  }else{
    if(length(W) != length(ind) | !is.numeric(W)){
      stop('If W is given, it must be numeric and of the same length to the amount of reference stations')
    }
    W <- W[ind]
  }

  return(mean(W*Est_reference[i,]))

}

CWM <- function(Est_target, Est_reference, i){
  #' Correlation weighting method
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
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  if(sum(ind) == 1){
    warning('only one station used')
    return(Est_reference[i])
  }

  w <- cor(Est_target,Est_reference,use = 'pairwise.complete.obs')
  Est_reference <- Est_reference[,w>0]
  w <- w[w>0]

  return(sum(w*Est_reference[i,])/sum(w))
}

