

TSLR_ranked <- function(Est_target, Est_reference,i, criterion='Max_Cor',cord=T,distans=F){
  #' fit linear models based on Theil-Sen single median, or Siegel repeated medians. It use mblm package
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' criterion is wich criterion is used to select reference estation
  #' Criterion options:
  #' Max_Cor
  #' Nearest; if nearest cordenates or distans must be given
  Est_reference <- as.matrix(Est_reference)
  as.numeric(Est_target)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(nrow(Est_reference) != length(Est_target)){stop('Not coerent dimentions of Est_reference and Est_target')}

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  if(criterion=='Max_Cor'){
    est <- Est_reference[,which.max(cor(Est_target,Est_reference,use = 'na.or.complete'))]
  }else if(criterion=='Nearest'){

    if(!is.logical(cord)){
      if(ncol(cord)!= 2){stop('not acceptable dimentions of cord')}
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
  esta <- rank(Est_target,na.last = 'keep')
  est <- rank(est,na.last = 'keep')
  a <- esta[ind]
  b <- est[ind]
  mod <- theil_sen_regression(a~b)
  rk <- mod$coefficients[1] + mod$coefficients[2] * est[i]
  val <- rank_val(target=Est_target,rk=rk)
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
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i')}
  if(nrow(Est_reference) != length(as.vector(Est_target))){stop('Not coerent dimentions of Est_reference and Est_target')}
  if(!is.logical(distans)){ if(length(distans)!=nrow(Est_reference)){stop('not coerent dimentions of distans width Est_reference')}}
  Est_reference <- as.matrix(Est_reference)
  ind <- which(!is.na(Est_reference[i,]))
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- as.matrix(Est_reference[,ind])

  if(criterion=='Max_Cor'){
    est <- Est_reference[,which.max(cor(Est_target,Est_reference,use = 'na.or.complete'))]
  }else if(criterion=='Nearest'){

    if(!is.logical(cord)){
      if(ncol(cord)!= 2){stop('not acceptable dimentions of cord')}
      cord <- cord[c(1,1+ind),]
      if(nrow(cord) != ncol(Est_reference)+1){stop('not coerent dimentions of cord and Est_reference')}
      distans <-  spDistsN1(as.matrix(cord[2:nrow(cord),]),as.matrix(cord[1,]),longlat = T)
    }else{distans <- distans[ind]}

    est <- Est_reference[,which.min(distans[ind])]
  }

  ind <- which(!is.na(est) & !is.na(Est_target))
  a <- Est_target[ind]
  b <- est[ind]
  mod <- theil_sen_regression(a~b)
  return(mod$coefficients[1] + mod$coefficients[2] * est[i] )
}
