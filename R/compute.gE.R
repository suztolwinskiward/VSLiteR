#' Compute just gE from latitude.
#' 
#' \code{compute.gE} computes an estimate of the seasonal growth response to insolation at
#' each month in the year at latitude phi. 
#'
#' The monthly values returned in gE are used to scale the growth responses to temperature
#' and soil moisture in VS-Lite.
#'
#' @param phi Latitude in degrees North.
#'
#' @export

compute.gE <- function(phi){
  
  gE <- matrix(NA,12,1);
  tmp <- daylength.factor.from.lat(phi,TRUE);
  L <- tmp$L;
  ndl <- tmp$ndl
  cdays <- tmp$cdays
  #
  for (t in 1:12){
    gE[t] <-  mean(ndl[(cdays[t]+1):cdays[t+1]]);
  }
  return(gE)
}
