
interpolate_all <- function(Est_target,CF = FALSE, CF_fun='median',fun,...){
  Est_target <- as.numeric(Est_target)
  checks <- makeAssertCollection()
  assert_logical(x = CF, len = 1, add = checks)
  assert_character(x = CF_fun, len = 1, add = checks)
  assert_function(x = fun, add = checks)
  assert_numeric(x = Est_target,
                 finite = TRUE,
                 all.missing = FALSE,
                 add = checks)
  reportAssertions(checks)

  my_fun <- possibly(.f = fun,otherwise = NA)
  na <- is.na(Est_target)
  if(any(na)){
    CC <- 0
    if(CF){
      new_target <- rep(NA,length(Est_target))
      na_index <- which(!na)
      for(i in na_index){
        new_target[i] <- my_fun(Est_target = Est_target,i=i,...) - Est_target[i]
      }
      if(CF_fun == 'median'){
        CC <- median(new_target,na.rm = T)
      }else if(CF_fun=='mean'){
        CC <- mean(new_target,na.rm = T)
      }else{
        stop('not suported CF_fun, only support median or mean')
      }
    }


    new_target <- Est_target
    na_index <- which(na)
    for(i in na_index){
      new_target[i] <- my_fun(Est_target = Est_target,i=i,...) - CC
    }
    if(any(is.na(new_target))){
      print('Not all NAs could be interpolated for those remaining NA try another method')
      print('This could be because lack of data of other stations for that time (row)')
    }
    return(new_target)

  }else{
    return(Est_target)
  }

}
