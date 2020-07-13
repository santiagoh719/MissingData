OI <- function(Est_target, Est_reference, i, cord= T, criterion = 'SBE' ){
  #' Optimal interpolation, Gandin and Kagan 1974 <- E ISCHEID 1999
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Both Est_target and Est_reference must have the same nrows
  #' Also the same start date.
  #' cord, if criterion=NN must be given, is a matrix of two colums, first Longitud, second Latitud
  #' The first row is the cordinate of target station, then by order the reference stations
  #' criterion is wich criterion must be use for first guest stimation,
  #' either Cor for correlation or NN for nearest neigthbor or None for not use firstguest

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
  assert_character(x = criterion, len = 1)
  reportAssertions(checks)

  cols <- apply(Est_reference,2,sd,na.rm=T) > 0.001
  Est_reference <- Est_reference[,cols]
  ind <- !is.na(Est_reference[i,])
  if(sum(ind) < 2){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  if(length(which(!is.na(Est_target)))<5){
    stop('Not possible to calculate correlation with target, to short length of series')
  }
  if(sd(Est_target,na.rm = T)<0.001){
    stop('Not possible to calculate correlation with target, to low standard desviation')
  }

  Cant <- ncol(Est_reference)
  first_guest <- vector(mode = 'double', length = Cant)
  c1 <- cor(Est_target,Est_reference,use = 'pairwise.complete.obs')
  c2 <- cor(Est_reference,use = 'pairwise.complete.obs')
  w <- c1 %*% solve(c2)
  c2 <- c2 - diag(1,nrow = nrow(c2))

  if(criterion == 'SBE'){
    for( j in 1:Cant){
      first_guest[j] <- Est_reference[i,which.max(c2[,j])]
    }
    first_guest_target <- Est_reference[i,which.max(c1)]
  }else if(criterion == 'NN'){

    if(ncol(cord)!= 2){
      stop('not acceptable dimentions of cord, must have 2 columns (Long, Lat)')
    }
    cord <- cord[c(1,1+which(ind)),]
    if(nrow(cord) != ncol(Est_reference)+1){
      stop('not coerent dimentions of distans width Est_reference')
    }

    distans <-  spDistsN1(as.matrix(cord[2:nrow(cord),]),as.matrix(cord[1,]),longlat = T)
    jaux <- which.min(distans)
    first_guest_target <- Est_reference[i,jaux]
    cord <- cord[2:nrow(cord),]
    for( j in 1:Cant){
      aux1 <- rbind(cord[j,],cord[-j,])
      distans <-  spDistsN1(as.matrix(aux1[2:nrow(aux1),]),as.matrix(aux1[1,]),longlat = T)
      jaux <- which.min(distans)
      if(jaux >= j){jaux <- jaux + 1}
      first_guest[j] <- Est_reference[i,jaux]
    }
  }else if(criterion == 'None'){
    first_guest_target <- 0
  }else{
    stop('not useable criterion given')
  }

  val <- first_guest_target + sum(w * (Est_reference[i,] - first_guest))
  return(val)

}
