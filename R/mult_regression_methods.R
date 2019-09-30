MLR <- function(Est_target, Est_reference,i, criterion='AIC',alpha=0.05, cord = NULL, dist = NULL){
  #' Use package SignifReg for multiple lineal regretion
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Criterion is wich criterion must be used for model selection
  #' posible uses are AIC, BIC and Cp (mallows Cp)
  if(nrow(Est_reference) != length(as.vector(Est_target))){stop('Not coerent dimentions of Est_reference and Est_target')}
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i')}

  Est_reference <- as.matrix(Est_reference)
  ind <- which(!is.na(Est_reference[i,]))
  if(sum(as.numeric(ind))< 1){
    stop('Not enougth stations')
  } else if(sum(as.numeric(ind)) == 1){
    Est_reference <- Est_reference[,ind]
    datos <- cbind(target = Est_target,predic = Est_reference)
    datos <- as.data.frame(datos)
    model <- lm(target ~ predic, data = datos)

    est <- t(t(datos[,2]) * model$coefficients[2])
  }else{
    Est_reference <- as.data.frame(Est_reference[,ind])

    colnames(Est_reference) <- paste0('V',colnames(Est_reference))
    datos <- cbind(target = Est_target,Est_reference)
    fit_model <- lm(target~.,data = datos)
    model <- SignifReg(fit = fit_model,scope = target~.,alpha=alpha,direction = 'both',criterion = criterion,trace = FALSE)

    est <- t(t(datos[,2:model$rank]) * model$coefficients[2:model$rank])
  }
  return(sum(est[i,])+ model$coefficients[1])
}

RMLR <- function(Est_target, Est_reference,i, M_Method= 'MM', cord = NULL, dist = NULL){
  #' Use package MASS for robust multiple lineal regretion, using M-Estimator propouse by Huber
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' M_Method, either M estimator or MM estimator
  #' First must select reference stations, all will be used
  if(nrow(Est_reference) != length(as.vector(Est_target))){stop('Not coerent dimentions of Est_reference and Est_target')}
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i')}

  Est_reference <- as.matrix(Est_reference)
  ind <- which(!is.na(Est_reference[i,]))
  if(sum(as.numeric(ind))< 1){
    stop('Not enougth stations')
  }else if(sum(as.numeric(ind)) == 1){
    Est_reference <- Est_reference[,ind]
    datos <- cbind(target = Est_target,predic = Est_reference)
    datos <- as.data.frame(datos)

    model <- rlm(target ~ predic, data = datos,method= M_Method)

    est <- t(t(datos[,2]) * model$coefficients[2])
  }else{
    Est_reference <- as.data.frame(Est_reference[,ind])
    colnames(Est_reference) <- paste0('V',colnames(Est_reference))
    datos <- cbind(target = Est_target,Est_reference)
    model <- rlm(target~.,datos, method= M_Method)

    est <- t(t(datos[,2:model$rank]) * model$coefficients[2:model$rank])
  }
  return(sum(est[i,])+ model$coefficients[1])
}

MLR_ranked <- function(Est_target, Est_reference,i, criterion='AIC',alpha=0.05, cord = NULL, dist = NULL){
  #' Use package SignifReg for multiple lineal regretion
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Criterion is wich criterion must be used for model selection
  #' posible uses are AIC, BIC and Cp (mallows Cp)
  Est_reference <- as.matrix(Est_reference)
  Est_target <- as.numeric(Est_target)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(nrow(Est_reference) != length(Est_target)){stop('Not coerent dimentions of Est_reference and Est_target')}

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  if(ncol(Est_reference) == 1){
    warning('only 1 possible predictor, common linear regression done')
    datos <- cbind(target = rank(Est_target,na.last = 'keep'),refe = rank(Est_reference,na.last = 'keep'))
    datos <- as.data.frame(datos)
    model <- lm(target ~ refe,data = datos )
    est <- t(t(datos[,2]) * model$coefficients[2])
    rk <- sum(est[i,])+model$coefficients[1]
  }else{
    Est_reference <- as.data.frame(Est_reference)
    colnames(Est_reference) <- paste0('V',colnames(Est_reference))
    datos <- cbind(target = rank(Est_target,na.last = 'keep'),apply(Est_reference,2,rank,na.last = 'keep'))
    datos <- as.data.frame(datos)
    fit_model <- lm(target~.,data = datos)
    model <- SignifReg(fit = fit_model,scope = target~.,alpha=alpha,direction = 'both',criterion = criterion,trace = FALSE)

    est <- t(t(datos[,2:model$rank]) * model$coefficients[2:model$rank])
    rk <- sum(est[i,])+model$coefficients[1]
  }
  val <- rank_val(target=Est_target,rk=rk)

  return(val)
}

RMLR_ranked <- function(Est_target, Est_reference,i, M_Method= 'MM', cord = NULL, dist = NULL){
  #' Use package MASS for robust multiple lineal regretion, using M-Estimator propouse by Huber
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' M_Method, either M estimator or MM estimator
  #' First must select reference stations, all will be used
  Est_reference <- as.matrix(Est_reference)
  Est_target <- as.numeric(Est_target)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(nrow(Est_reference) != length(Est_target)){stop('Not coerent dimentions of Est_reference and Est_target')}

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  Est_reference <- as.data.frame(Est_reference)
  colnames(Est_reference) <- paste0('V',colnames(Est_reference))

  datos <- cbind(target = rank(Est_target,na.last = 'keep'),apply(Est_reference,2,rank,na.last = 'keep'))
  datos <- as.data.frame(datos)
  model <- rlm(target~.,datos, method= M_Method)

  est <- t(t(datos[,2:model$rank]) * model$coefficients[2:model$rank])
  rk <- sum(est[i,])+model$coefficients[1]

  val <- rank_val(target=Est_target,rk=rk)
  return(val)
}
