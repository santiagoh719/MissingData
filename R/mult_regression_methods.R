MLR <- function(Est_target, Est_reference,i, criterion='AIC',alpha=0.05){
  #' Use package SignifReg for multiple lineal regretion
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Criterion is wich criterion must be used for model selection
  #' posible uses are AIC, BIC and Cp (mallows Cp)
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
  assert_numeric(x = alpha, lower = 0, upper = 1,
                 len = 1, any.missing = FALSE, add = checks)
  assert_character(x = criterion, len = 1, add = checks)
  reportAssertions(checks)

  ind <- !is.na(Est_reference[i,])
  if(sum(as.numeric(ind))< 1){

    stop('Not enougth stations')

  } else if(sum(as.numeric(ind)) == 1){

    Est_reference <- Est_reference[,ind]
    datos <- cbind(target = Est_target,predic = Est_reference)
    datos <- as.data.frame(datos)
    model <- lm(target ~ predic, data = datos)

    est <- t(t(datos[,2]) * model$coefficients[2])
    return(sum(est[i,])+ model$coefficients[1])

  }else{
    Est_reference <- as.data.frame(Est_reference[,ind])

    datos <- cbind(target = Est_target,Est_reference)

    fit_model <- lm(target~.,data = datos)

    model <- SignifReg(fit = fit_model,
                       #scope = datos,
                       alpha=alpha,
                       direction = 'both',
                       criterion = criterion,
                       trace = FALSE)

    if(model$rank != 1){

      est <- t(t(datos[,gsub(x = names(model$coefficients)[-1],
                             pattern = '`',
                             replacement = '')]) *
                 model$coefficients[2:model$rank])

      return(sum(est[i,])+ model$coefficients[1])

    }else{

      return(model$coefficients[1])

    }
  }
}

RMLR <- function(Est_target, Est_reference,i, M_Method= 'MM'){
  #' Use package MASS for robust multiple lineal regretion, using M-Estimator propouse by Huber
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' M_Method, either M estimator or MM estimator
  #' First must select reference stations, all will be used
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
  assert_character(x = M_Method, len = 1, add = checks)
  reportAssertions(checks)

  ind <- (!is.na(Est_reference[i,]))
  if(sum(as.numeric(ind))< 1){

    stop('Not enougth stations')

  }else if(sum(as.numeric(ind)) == 1){

    Est_reference <- Est_reference[,ind]
    datos <- cbind(target = Est_target,predic = Est_reference)
    datos <- as.data.frame(datos)

    model <- rlm(target ~ predic, data = datos,method= M_Method)

    est <- t(t(datos[,2]) * model$coefficients[2])
    return(sum(est[i,])+ model$coefficients[1])

  }else{
    Est_reference <- as.data.frame(Est_reference[,ind])
    colnames(Est_reference) <- paste0('V',colnames(Est_reference))
    datos <- cbind(target = Est_target,Est_reference)
    model <- rlm(target~.,datos, method= M_Method)

    if(model$rank != 1){
      est <- t(t(datos[,names(model$coefficients)[-1]]) * model$coefficients[2:model$rank])
      return(sum(est[i,])+ model$coefficients[1])
    }else{
      return(model$coefficients[1])
    }
  }
}

MLR_ranked <- function(Est_target, Est_reference,i, criterion='AIC',alpha=0.05){
  #' Use package SignifReg for multiple lineal regretion
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' Criterion is wich criterion must be used for model selection
  #' posible uses are AIC, BIC and Cp (mallows Cp)
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
  assert_numeric(x = alpha, lower = 0, upper = 1,
                 len = 1, any.missing = FALSE, add = checks)
  assert_character(x = criterion, len = 1, add = checks)
  reportAssertions(checks)

  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  if(ncol(Est_reference) == 1){

    warning('only 1 possible predictor, common linear regression done')
    datos <- cbind(target = rank(Est_target,na.last = 'keep'),
                   refe = rank(Est_reference,na.last = 'keep'))
    datos <- as.data.frame(datos)
    model <- lm(target ~ refe,data = datos )
    est <- t(t(datos[,2]) * model$coefficients[2])
    return(.rank_val(target = Est_target, rk = sum(est[i,])+ model$coefficients[1]))

  }else{
    Est_reference <- as.data.frame(Est_reference)
    colnames(Est_reference) <- paste0('V',colnames(Est_reference))
    datos <- cbind(target = rank(Est_target,na.last = 'keep'),
                   apply(Est_reference,2,rank,na.last = 'keep'))
    datos <- as.data.frame(datos)
    fit_model <- lm(target~.,data = datos)
    model <- SignifReg(fit = fit_model,
                       #scope = target~.,
                       alpha=alpha,
                       direction = 'both',
                       criterion = criterion,
                       trace = FALSE)

    est <- t(t(datos[,2:model$rank]) * model$coefficients[2:model$rank])
    rk <- sum(est[i,])+model$coefficients[1]

    if(model$rank != 1){
      est <- t(t(datos[,names(model$coefficients)[-1]]) * model$coefficients[2:model$rank])
      return(.rank_val(target = Est_target, rk = sum(est[i,])+ model$coefficients[1]))
    }else{
      return(.rank_val(target = Est_target, rk = model$coefficients[1]))
    }
  }

}

RMLR_ranked <- function(Est_target, Est_reference,i, M_Method= 'MM'){
  #' Use package MASS for robust multiple lineal regretion, using M-Estimator propouse by Huber
  #' Est_target a vector with the data of the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' M_Method, either M estimator or MM estimator
  #' First must select reference stations, all will be used
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
  assert_character(x = M_Method, len = 1, add = checks)
  reportAssertions(checks)
  ind <- !is.na(Est_reference[i,])
  if(all(!(ind))){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  Est_reference <- as.data.frame(Est_reference)
  colnames(Est_reference) <- paste0('V',colnames(Est_reference))

  datos <- cbind(target = rank(Est_target,na.last = 'keep'),
                 apply(Est_reference,2,rank,na.last = 'keep'))
  datos <- as.data.frame(datos)
  model <- rlm(target~.,datos, method= M_Method)

  if(model$rank != 1){
    est <- t(t(datos[,names(model$coefficients)[-1]]) * model$coefficients[2:model$rank])
    return(.rank_val(target = Est_target, rk = sum(est[i,])+ model$coefficients[1]))
  }else{
    return(.rank_val(target = Est_target, rk = model$coefficients[1]))
  }
}
