
patern_distans_change <- function(Est_target, Est_reference, reg = 1){
  #' select the station acording similarity paterns and create a new distans
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' reg is the amount of values to use to create the paterns
  #' The output is the new "distans" estimated with the paterns
  #' Teegavarapu y chandramouli 2005
  #'
  Est_target <- as.numeric(Est_target)
  Est_reference <- as.matrix(Est_reference)
  checks <- makeAssertCollection()
  assert_matrix(x = Est_reference,
                mode = 'numeric',
                all.missing = FALSE, add = checks)
  assert_integerish(x = reg, lower = 1, len = 1,
                    upper = 2, add = checks)
  assert_numeric(x = Est_target,
                 len = dim(Est_reference)[1],
                 finite = TRUE, all.missing = FALSE, add = checks)
  reportAssertions(checks)


  n <- length(Est_target)
  paterns <- permutations(n= 3,r = reg,v = c(2,0,1),repeats.allowed = T)
  npat <- nrow(paterns)
  targ_pat <- Est_target[1:(n-1)] - Est_target[2:n]
  targ_pat[targ_pat==0] <- 0
  targ_pat[targ_pat > 0] <- 1
  targ_pat[targ_pat < 0] <- 2

  pat_cant <- matrix(0,nrow = npat,ncol = ncol(Est_reference)+1)
  row <- 0
  for(i in 1:npat){
    row <- row+1
    pat_cant[row,1] <- length(.patern_search(patern = paterns[i,],
                                            vector = targ_pat))
  }

  ref_pat <- Est_reference[1:(n-1),] - Est_reference[2:n,]

  for(i in 1:ncol(ref_pat)){
    ref_pat[ref_pat[,i] == 0,i] <- 0
    ref_pat[ref_pat[,i] > 0,i] <- 1
    ref_pat[ref_pat[,i] < 0,i] <- 2
  }
  for(cols in 2:(ncol(Est_reference)+1)){
    row <- 0
    for(i in 1:npat){
      row <- row+1
      pat_cant[row,cols] <- length(.patern_search(patern = paterns[i,],
                                                 vector = ref_pat[,cols-1]))
    }
  }

  pat_cant <- t(t(pat_cant) * (1/colSums(pat_cant)) ) * 100
  distans_mat <- abs(pat_cant[,2:ncol(pat_cant)] - matrix(pat_cant[,1],nrow= nrow(pat_cant), ncol = ncol(pat_cant)-1))
  return(colSums(distans_mat))
}


.patern_search <- function(patern, vector){
  #' patern is the patern to look for in the vector,
  #' character is not support
  #' complex numbers are not supported
  #' the output is a a vector with the index of the pattern

  npat <- length(patern)
  ind <- which(vector == patern[1])
  ind <- ind[ ind < (length(vector)-npat+1)]

  if(npat == 1){return(ind)}
  for(i in 2:npat){
    ind <- ind[vector[ind+i-1] == patern[i]]
  }
  return(ind)
}

.Ang <- function(cord_ref, cord_j, kk){

  lat_ref <- pi/2 - as.numeric(unlist(cord_ref[,2]))*pi/180
  lon_ref <-  as.numeric(unlist(cord_ref[,1]))*pi/180

  lat_j <-  pi/2 - as.numeric(unlist(cord_j[2]))*pi/180
  lon_j <-  as.numeric(unlist(cord_j[1]))*pi/180

  lat_k <- lat_ref[kk]
  lon_k <- lon_ref[kk]

  lat_l <- lat_ref[-kk]
  lon_l <- lon_ref[-kk]

  delta_jk <- lon_j - lon_k

  cos_akj <- cos(lat_k) * cos(lat_j) + sin(lat_j) * sin(lat_k) * cos(delta_jk)

  delta_lj <- lon_l - lon_j
  delta_lk <- lon_l - lon_k

  cos_alj <- cos(lat_j) * cos(lat_l) + sin(lat_j) * sin(lat_l) * cos(delta_lj)
  cos_alk <- cos(lat_k) * cos(lat_l) + sin(lat_k) * sin(lat_l) * cos(delta_lk)

  akj <- acos(cos_akj)
  alj <- acos(cos_alj)

  Ang <- acos(
    (cos_alk - cos_alj * cos_akj) /
      (sin(alj) * sin(akj))
    )

  return(Ang)
}



.rank_val <- function(target,rk){
  #' Get the values of the ranked rk acording the target serie
  #' rk must be numerical
  #' output: a vector of length rk, with the values for that rank
  target <- as.numeric(target)
  rk <- as.numeric(rk)
  est_ord <- target[order(target,na.last = NA)]
  esta <- rank(est_ord,na.last = NA)
  j <- length(rk)
  val <- vector(mode = 'double', length = j)
  for(i in 1:j){
    if(rk[i] <= min(esta)){
      val[i] <- est_ord[which.min(esta)[1]]
    }else if(rk[i] >= max(esta)){
      val[i] <- est_ord[which.max(esta)[1]]
    }else if(any(esta==rk[i])){
      val[i] <- est_ord[which(esta == rk[i])[1]]
    }else{
      y0 <- which.max(esta[(esta-rk[i])<0])
      y1 <- which((esta-rk[i])>0)[1]
      val[i] <- est_ord[y0] + (rk[i]-esta[y0]) * (est_ord[y1]-est_ord[y0])
    }
  }
  return(val)
}
