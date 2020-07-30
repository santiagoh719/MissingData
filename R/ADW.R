ADW <- function(Est_reference, i, cord, CDD, m = 4){
  Est_reference <- as.matrix(Est_reference)
  cord <- as.matrix(cord)
  checks <- makeAssertCollection()
  assert_matrix(x = Est_reference,
                mode = 'numeric',
                all.missing = FALSE, add = checks)
  assert_integerish(x = i, lower = 1,
                    upper = dim(Est_reference)[1], add = checks)
  assert_matrix(x = cord, ncols = 2,
                nrows = dim(Est_reference)[2]+1,
                mode = 'numeric', any.missing = FALSE, add = checks)
  assert_numeric(x = CDD, finite = TRUE, len = 1,
                 add = checks)
  assert_numeric(x = m, finite = TRUE, len = 1,
                 add = checks)
  reportAssertions(checks)

  ind <- !is.na(Est_reference[i,])
  if(sum(as.numeric(ind))<= 4){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  cord_ref <- cord[which(ind)+1,]
  cord_target <- cord[1,]

  n_est <- ncol(Est_reference)
  distance <- spDistsN1(as.matrix(cord_ref),
                        as.matrix(cord_target),
                        longlat = T)

  W_first <- (exp(-distance/CDD))^m

  A_k <- vector(mode = 'double', length = n_est )
  for(kk in 1:n_est){
    A_k[kk] <- sum(W_first[-kk] *(1 - cos(.Ang(cord_ref,
                                               cord_target,
                                               kk)
                                          )
                                  )
                   ) /
      sum(W_first[-kk])
  }
  WW <- W_first * (1+A_k)
  WW <- WW / sum(WW)
  return(sum(Est_reference[i,] * (WW)))
}
