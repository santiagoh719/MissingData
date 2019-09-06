
interpolate_all <- function(Est_target,fun,...){
  my_fun <- possibly(.f = fun,otherwise = NA)
  na <- is.na(Est_target)
  if(any(na)){
    new_target <- Est_target
    na_index <- which(na)
    for(i in na_index){
      new_target[i] <- my_fun(Est_target = Est_target,i=i,...)
    }
    if(any(is.na(new_target))){
      print('Not al NAs could be interpolated for those remaining NA try another method')
    }
    return(new_target)
  }else{
    return(Est_target)
  }

}
