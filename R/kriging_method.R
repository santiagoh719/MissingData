Ord_kriging <- function( Est_reference, i, cord, method = NULL){
  #' Ordinary kriging using sp and gstat packge; only interpolate 1 value to the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' cord is a matrix of two colums, first Longitud, second Latitud
  #' The first row is the cordinate of target station, then by order the reference stations
  #' method one of c("Exp", "Mat", "Gau","Sph"), 'Gau' default

  Est_reference <- as.matrix(Est_reference)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(ncol(cord)!= 2){stop('not acceptable dimentions of cord, must have 2 columns (Long, Lat)')}
  if(nrow(cord) != ncol(Est_reference)+1){stop('not coerent dimentions of distans width Est_reference')}



  ind <- !is.na(Est_reference[i,])
  if(sum(as.numeric(ind))<= 5){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  # Values to perform kriging
  predictors <- SpatialPointsDataFrame(coords = cord[which(ind)+1,],
                                       data = data.frame(val=Est_reference[i,]),
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))


  # Cordinates where kriging will be computed
  new_cord <-SpatialPoints(coords = cord[1,],
                           proj4string = CRS("+proj=longlat +datum=WGS84"))

  vario <- variogram(val~1,predictors)
  if(is.null(method)){
    fit_vario <- fit.variogram(vario,vgm(c("Exp", "Mat", "Gau","Sph")),fit.kappa = TRUE)

  }else{
    fit_vario <- fit.variogram(vario,vgm(method),fit.kappa = TRUE)

  }
  if(!exists('fit_vario')){
    warning(paste0('No convergence, to many iterations, prove with other method for time ', i))
    val <- NA
  } else if(fit_vario$range[2] < 0){
    warning(paste0('variogram range can never be negative, prove with other method for time', i))
    val <- NA
  }else{
    ok <- krige(val~1,predictors,new_cord,fit_vario)
    val <- ok$var1.pred[1]
  }

  return(val)
}

Uni_kriging_Alt <- function(Est_reference, i, cord,method = NULL ){
  #' Universal kriging using sp and gstat packge; only interpolate 1 value to the target station
  #' Use as a explicative variable the altitude
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' cord is a matrix of three colums, first Longitud, second Latitud, altitude
  #' The first row is the cordinate of target station, then by order the reference stations
  #' method one of c("Exp", "Mat", "Gau","Sph"), 'Gau' default


  Est_reference <- as.matrix(Est_reference)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(ncol(cord)!= 3){stop('not acceptable dimentions of cord, must have 3 columns (Long, Lat, Z)')}
  if(nrow(cord) != ncol(Est_reference)+1){stop('not coerent dimentions of distans width Est_reference')}



  ind <- !is.na(Est_reference[i,])
  if(sum(as.numeric(ind))<= 5){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  # Values to perform kriging
  predictors <- SpatialPointsDataFrame(coords = cord[which(ind)+1,1:2],
                                       data = data.frame(val=Est_reference[i,], Altitude = cord[which(ind)+1,3]),
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))


  # Cordinates where kriging will be computed
  new_cord <-SpatialPointsDataFrame(coords = cord[1,1:2],
                                    data = data.frame(Altitude=cord[1,3]),
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))


  vario <- variogram(val~Altitude,predictors)
  if(is.null(method)){
    fit_vario <- fit.variogram(vario,vgm(c("Exp", "Mat", "Gau","Sph")),fit.kappa = TRUE)

  }else{
    fit_vario <- fit.variogram(vario,vgm(method),fit.kappa = TRUE)

  }
  if(!exists('fit_vario')){
    warning(paste0('No convergence, to many iterations, prove with other method for time ', i))
    val <- NA
  } else if(fit_vario$range[2] < 0){
    warning(paste0('variogram range can never be negative, prove with other method for time', i))
    val <- NA
  }else{
    ok <- krige(val~Altitude,predictors,new_cord,fit_vario)
    val <- ok$var1.pred[1]
  }
  return(val)
}

Uni_kriging <- function(Est_reference,Est_reference_var, var_target, i, cord, method = NULL){
  #' Universal kriging using sp and gstat packge; only interpolate 1 value to the target station
  #' Use as a explicative variable the altitude
  #' Est_reference a matrix or data.frame with the reference stations by columms width the variable to fit
  #' Est_reference_var a matrix or data.frame with the reference stations by columms width the predictor variable
  #' var_target, single numeric value, the predictor variable at the target station in the ieth index
  #' i the index to by estemated
  #' cord is a matrix of three colums, first Longitud, second Latitud
  #' The first row is the cordinate of target station, then by order the reference stations
  #' method one of c("Exp", "Mat", "Gau","Sph"), 'Gau' default

  Est_reference <- as.matrix(Est_reference)
  Est_reference_var <- as.numeric(Est_reference_var)
  if( i > nrow(Est_reference) | i<0 ){stop('Not aceptable value of i or Est_reference is not matrix nor data.frame')}
  if(ncol(cord)!= 2){stop('not acceptable dimentions of cord, must have 3 columns (Long, Lat, Z)')}
  if(nrow(cord) != ncol(Est_reference)+1){stop('not coerent dimentions of distans width Est_reference')}
  if(dim(Est_reference)[2] != length(Est_reference_var)){stop('not coerent dimentions of Est_reference_var and Est_reference')}

  ind <- !is.na(Est_reference[i,])
  if(sum(as.numeric(ind))<= 5){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]


  # Values to perform kriging
  predictors <- SpatialPointsDataFrame(coords = cord[which(ind)+1,],
                                       data = data.frame(val=Est_reference[i,], var = Est_reference_var[which(ind)]),
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))


  # Cordinates where kriging will be computed
  new_cord <-SpatialPointsDataFrame(coords = cord[1,],
                                    data = data.frame(var=var_target),
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))


  vario <- variogram(val~var,predictors)
  if(is.null(method)){
    fit_vario <- fit.variogram(vario,vgm(c("Exp", "Mat", "Gau","Sph")),fit.kappa = TRUE)

  }else{
    fit_vario <- fit.variogram(vario,vgm(method),fit.kappa = TRUE)

  }
  if(!exists('fit_vario')){
    warning(paste0('No convergence, to many iterations, prove with other method for time ', i))
    val <- NA
  } else if(fit_vario$range[2] < 0){
    warning(paste0('variogram range can never be negative, prove with other method for time', i))
    val <- NA
  }else{
    ok <- krige(val~var,predictors,new_cord,fit_vario)
    val <- ok$var1.pred[1]
  }
  return(val)
}
