#' Compute normalized daylength.
#' 
#' \code{daylength.factor.from.lat} computes a normalized daylength factor for each month, 
#' used to modulate growth responses in VS-Lite to account for seasonal variations in 
#' insolation.
#'  
#' The small difference in calculation for leap-years is negleced.
#' 
#' @param phi Latitude in degrees.
#' @param return.ndl.and.cdays Logical; return auxiliary variabiles ndl and cdays? Default
#' is FALSE.
#' 
#' @export

daylength.factor.from.lat <- function(phi,return.ndl.and.cdays=FALSE){
  latr <- phi*pi/180;  # change to radians
  ndays <- cbind(0,31,28,31,30,31,30,31,31,30,31,30,31);
  cdays <- cumsum(ndays);
  sd <- t(asin(sin(pi*23.5/180) * sin(pi * (((1:365) - 80)/180))));   # solar declination
  y <- -tan(matrix(latr,365,1)) * t(tan(sd));
  # bound y within (-1,1):
  y[y >= 1] <- 1;
  y[y <= -1] <- -1;
  
  hdl <- acos(y);
  dtsi <- hdl * sin(matrix(latr,365,1)) * t(sin(sd)) +
    (matrix(cos(latr),365,1)) * t(cos(sd)) * sin(hdl);
  ndl <- dtsi/max(dtsi); # normalized day length
  
  # calculate mean monthly daylength (used for evapotranspiration in soil moisture calcs)
  jday <- cdays[1:12] +.5*ndays[2:13];
  m.star <- 1-tan(phi*pi/180)*tan(23.439*pi/180*cos(jday*pi/182.625));
  # bound m.star between 0 and 2:
  m.star[m.star < 0] <- 0
  m.star[m.star > 2] <- 2;
  
  nhrs <- 24*acos(1-m.star)/pi; # the number of hours in the day in the middle of the month
  # mean normalized daylength factor:
  L <- (ndays[2:13]/30) * (nhrs/12)
  if(return.ndl.and.cdays){
    out <- list(L,ndl,cdays)
    names(out) <- c("L","ndl","cdays")
    return(out)
  }else{
    return(L); 
  }
}
