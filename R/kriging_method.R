Ord_kriging <- function( Est_reference, i, cord, method = NULL){
  #' Ordinary kriging using sp and gstat packge; only interpolate 1 value to the target station
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' cord is a matrix of two colums, first Longitud, second Latitud
  #' The first row is the cordinate of target station, then by order the reference stations
  #' method one of c("Exp", "Mat", "Gau","Sph"), 'Gau' default
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
  reportAssertions(checks)

  ind <- !is.na(Est_reference[i,])
  if(sum(as.numeric(ind))<= 6){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  # Values to perform kriging
  predictors <- SpatialPointsDataFrame(coords = cord[which(ind)+1,],
                                       data = data.frame(val=Est_reference[i,]),
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))


  # Cordinates where kriging will be computed
  new_cord <-SpatialPoints(coords = matrix(cord[1,],1,2),
                           proj4string = CRS("+proj=longlat +datum=WGS84"))

  vario <- variogram(val~1,predictors)
  if(is.null(method)){
    fit_vario <- fit.variogram(vario,vgm(c("Exp", "Mat", "Gau","Sph")),
                               fit.kappa = TRUE)

  }else{
    fit_vario <- fit.variogram(vario,vgm(method),
                               fit.kappa = TRUE)
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

KED_Alt <- function(Est_reference, i, cord,method = NULL){
  #' kriging with external drift using sp and gstat packge; only interpolate 1 value to the target station
  #' Use as a explicative variable the altitude
  #' Est_reference a matrix or data.frame with the reference stations by columms
  #' i the index to by estemated
  #' cord is a matrix of three colums, first Longitud, second Latitud, altitude
  #' The first row is the cordinate of target station, then by order the reference stations
  #' method one of c("Exp", "Mat", "Gau","Sph"), 'Gau' default

  Est_reference <- as.matrix(Est_reference)
  cord <- as.matrix(cord)
  checks <- makeAssertCollection()
  assert_matrix(x = Est_reference,
                mode = 'numeric',
                all.missing = FALSE, add = checks)
  assert_integerish(x = i, lower = 1,
                    upper = dim(Est_reference)[1], add = checks)
  assert_matrix(x = cord, ncols = 3,
                nrows = dim(Est_reference)[2]+1,
                mode = 'numeric', any.missing = FALSE, add = checks)
  reportAssertions(checks)

  ind <- !is.na(Est_reference[i,])
  if(sum(as.numeric(ind))<= 6){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]

  # Values to perform kriging
  predictors <- SpatialPointsDataFrame(coords = cord[which(ind)+1,1:2],
                                       data = data.frame(val=Est_reference[i,],
                                                         Altitude = cord[which(ind)+1,3]),
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))


  # Cordinates where kriging will be computed
  new_cord <-SpatialPointsDataFrame(coords = matrix(cord[1,1:2],1,2),
                                    data = data.frame(Altitude=cord[1,3]),
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))


  vario <- variogram(val~Altitude, predictors)
  if(is.null(method)){
    fit_vario <- fit.variogram(vario,vgm(c("Exp", "Mat", "Gau","Sph")),
                               fit.kappa = TRUE)

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

KED <- function(Est_reference,
                Est_reference_var,
                var_target, i, cord, method = NULL){
  #' kriging with external drift using sp and gstat packge; only interpolate 1 value to the target station
  #' Use as a explicative variable the altitude
  #' Est_reference a matrix or data.frame with the reference stations by columms width the variable to fit
  #' Est_reference_var a matrix or data.frame with the reference stations by columms width the predictor variable
  #' var_target, single numeric value, the predictor variable at the target station in the ieth index
  #' i the index to by estemated
  #' cord is a matrix of three colums, first Longitud, second Latitud
  #' The first row is the cordinate of target station, then by order the reference stations
  #' method one of c("Exp", "Mat", "Gau","Sph"), 'Gau' default

  Est_reference_var <- as.matrix(Est_reference_var)
  Est_reference <- as.matrix(Est_reference)
  cord <- as.matrix(cord)
  checks <- makeAssertCollection()
  assert_matrix(x = Est_reference,
                mode = 'numeric',
                all.missing = FALSE, add = checks)
  assert_numeric(x = var_target, finite = TRUE,
                 len = 1, add = checks)
  assert_matrix(x = Est_reference_var,
                ncols = ncol(Est_reference),
                nrows = nrow(Est_reference),
                mode = 'numeric',
                all.missing = FALSE, add = checks)
  assert_integerish(x = i, lower = 1,
                    upper = dim(Est_reference)[1], add = checks)
  assert_matrix(x = cord, ncols = 2,
                nrows = dim(Est_reference)[2]+1,
                mode = 'numeric', any.missing = FALSE, add = checks)
  reportAssertions(checks)

  ind <- !is.na(Est_reference[i,]) & !is.na(Est_reference_var[i,])
  if(sum(as.numeric(ind))<= 6){
    stop('Not enougth stations')
  }
  Est_reference <- Est_reference[,ind]
  Est_reference_var <- Est_reference_var[,ind]

  # Values to perform kriging
  predictors <- SpatialPointsDataFrame(coords = cord[which(ind)+1,],
                                       data = data.frame(val = Est_reference[i,],
                                                         var = Est_reference_var[i,]),
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))


  # Cordinates where kriging will be computed
  new_cord <-SpatialPointsDataFrame(coords = matrix(cord[1,],1,2),
                                    data = data.frame(var=var_target),
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))


  vario <- variogram(val~var,predictors)
  if(is.null(method)){
    fit_vario <- fit.variogram(vario,vgm(c("Exp", "Mat", "Gau","Sph")),
                               fit.kappa = TRUE)

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
