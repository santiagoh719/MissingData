

TSLR_ranked <- function(Est_target, Est_reference,i, criterion='Max_Cor',cord=T,distans=F){
  #' fit linear models based on Theil-Sen single median, or Siegel repeated medians. It use mblm package
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' criterion is wich criterion is used to select reference estation
  #' Criterion options:
  #' Max_Cor
  #' Nearest; if nearest cordenates or distans must be given
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
  assert_character(x = criterion, len = 1, add = checks)
  reportAssertions(checks)

  ind <- which(!is.na(Est_reference[i,]))
  if(length(ind)<1){
    stop('Not enougth stations')
  }

  Est_reference <- as.matrix(Est_reference[,ind])

  if(criterion=='Max_Cor'){
    est <- Est_reference[,which.max(cor(Est_target,Est_reference,use = 'na.or.complete'))]
  }else if(criterion=='Nearest'){

    if(!is.logical(cord)){
      if(ncol(cord)!= 2){
        stop('not acceptable dimentions of cord')
      }
      cord <- cord[c(1,1+ind),]
      if(nrow(cord) != ncol(Est_reference)+1){
        stop('not coerent dimentions of cord and Est_reference')
      }
      distans <-  spDistsN1(as.matrix(cord[2:nrow(cord),]),as.matrix(cord[1,]),longlat = T)
    }else{

      if(length(distans)!=nrow(Est_reference)){
        stop('not coerent dimentions of distans width Est_reference')
      }
      distans <- distans[ind]
    }

    est <- Est_reference[,which.min(distans[ind])]
  }

  ind <- which(!is.na(est) & !is.na(Est_target))
  est_a <- rank(Est_target,na.last = 'keep')
  est_b <- rank(est,na.last = 'keep')
  a <- est_a[ind]
  b <- est_b[ind]
  mod <- theil_sen_regression(a~b)
  rk <- mod$coefficients[1] + mod$coefficients[2] * est_b[i]
  val <- .rank_val(target=Est_target,rk=rk)
  return(val)
}

TSLR <- function(Est_target, Est_reference,i, criterion='Max_Cor',cord=T,distans=F){
  #' fit linear models based on Theil-Sen single median, or Siegel repeated medians. It use mblm package
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' criterion is wich criterion is used to select reference estation
  #' Criterion options:
  #' Max_Cor
  #' Nearest; if nearest cordenates or distans must be given
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
  assert_character(x = criterion, len = 1, add = checks)
  reportAssertions(checks)

  ind <- which(!is.na(Est_reference[i,]))
  if(length(ind)<1){
    stop('Not enougth stations')
  }

  Est_reference <- as.matrix(Est_reference[,ind])

  if(criterion=='Max_Cor'){
    est <- Est_reference[,which.max(cor(Est_target,Est_reference,use = 'na.or.complete'))]
  }else if(criterion=='Nearest'){

    if(!is.logical(cord)){
      if(ncol(cord)!= 2){
        stop('not acceptable dimentions of cord')
      }
      cord <- cord[c(1,1+ind),]
      if(nrow(cord) != ncol(Est_reference)+1){
        stop('not coerent dimentions of cord and Est_reference')
      }
      distans <-  spDistsN1(as.matrix(cord[2:nrow(cord),]),as.matrix(cord[1,]),longlat = T)
    }else{

      if(length(distans)!=nrow(Est_reference)){
        stop('not coerent dimentions of distans width Est_reference')
      }
      distans <- distans[ind]
    }

    est <- Est_reference[,which.min(distans[ind])]
  }

  ind <- which(!is.na(est) & !is.na(Est_target))
  a <- Est_target[ind]
  b <- est[ind]
  mod <- theil_sen_regression(a~b)
  return(mod$coefficients[1] + mod$coefficients[2] * est[i] )
}
