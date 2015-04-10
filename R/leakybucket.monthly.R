#' Run CPC Leaky Bucket model (monthly version).
#' 
#' \code{leakybucket.monthly} simulates soil moisture with coarse monthly time step
#' 
#' Monthly time step version implemented by Suz Tolwinski-Ward in 2011; ported to R
#' by SETW in 2015. Implementation of CPC Leaky Bucket model as described in
#' Huang et al., 'Analysis of Model-Calculated Soil Moisture over the United States 
#' (1931-1993) and Applications to Long-Range Temperature Forecasts,' J. Clim. (1995).
#' Fixes to model consistent with Matlab VSLite version 2.5. 
#' 
#' @param syear Start year of simulation.
#' @param eyear End year of simulation.
#' @param phi Latitude of site (in degrees N).
#' @param T (12 x Nyrs) Matrix of ordered mean monthly temperatures (in degEes C).
#' @param P (12 x Nyrs) Matrix of ordered accumulated monthly precipitation (in mm).
#' @param Mmax Scalar maximum soil moisture held by the soil (in v/v).
#' @param Mmin Scalar minimum soil moisture (for error-catching) (in v/v).
#' @param alph Scalar runoff parameter 1 (in inverse months).
#' @param m.th Scalar runoff parameter 3 (unitless).
#' @param mu.th Scalar runoff parameter 2 (unitless).
#' @param rootd Scalar root/"bucket" depth (in mm).
#' @param M0 Initial value for previous month's soil moisture at t = 1 (in v/v).
#'
#' @return M Soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs).
#' @return potEv Potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm).
#'
#' @seealso \code{\link{leakybucket.submonthly}}
#' 
#' @export
####################################################################################################

## LEAKY BUCKET WITHOUT SUBSTEPPING ##
leakybucket.monthly <- function(syear,eyear,phi,T,P,Mmax = 0.76,Mmin = 0.01,alph = 0.093,
                                m.th = 4.886,mu.th = 5.8,rootd = 1000,M0 = .2){
  
  iyear <- syear:eyear;
  nyrs <- length(iyear);
  
  # Storage for growth response output variables (size [12 x Nyears]):
  M <- potEv <- matrix(NA,12,nyrs);
  
  if(M0 < 0.){M0 <- 200/rootd;}
  
  # Compute normalized daylength (neglecting small difference in calculation for leap-years)
  L <- daylength.factor.from.lat(phi); 
  
  # Pre-calculation of istar and I, using input T to compute the climatology:
  Tm <- rowMeans(T);
  istar <- (Tm/5)^1.514; 
  istar[Tm < 0] <- 0;
  I <- sum(istar);
  
  # precalculation of the exponent alpha in the Thornwaite (1948) equation:
  a < (6.75e-7)*I^3 - (7.71e-5)*I^2 + (1.79e-2)*I + 0.49;
  
  #########################################################################################
  #### -- year cycle -- ####
  # syear = start (first) year of simulation
  # eyear = end (last) year of simulation
  # cyear = year the model is currently working on
  # iyear = index of simulation year
  
  for (cyear in 1:nyrs){     # begin cycling over years
    #########################################################################################
    for (t in 1:12){  # begin cycling over months in a year
      
      ##### Compute potential evapotranspiration for current month after Thornthwaite:
      if ( T[t,cyear] < 0 ){Ep = 0;}
      if ( T[t,cyear] >= 0 && T[t,cyear] < 26.5 ){Ep <- 16*L[t]*(10*T[t,cyear]/I)^a;}
      if ( T[t,cyear] >= 26.5 ){Ep <- -415.85 + 32.25*T[t,cyear] - .43* T[t,cyear]^2;}
      potEv[t,cyear] <- Ep;
      
      ##### Now calculate soil moisture according to the CPC Leaky Bucket model
      ##### (see J. Huang et al, 1996).
      
      if (t > 1){
        # evapotranspiration:
        Etrans <- Ep*M[t-1,cyear]*rootd/(Mmax*rootd);
        # groundwater loss via percolation:
        G <- mu.th*alph/(1+mu.th)*M[t-1,cyear]*rootd;
        # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
        R <- P[t,cyear]*(M[t-1,cyear]*rootd/(Mmax*rootd))^m.th +
        (alph/(1+mu.th))*M[t-1,cyear]*rootd;
        dWdt <- P[t,cyear] - Etrans - R - G;
        M[t,cyear] <- M[t-1,cyear] + dWdt/rootd;
      }
      
      if( t == 1 && cyear > 1){
        # evapotranspiration:
        Etrans <- Ep*M[12,cyear-1]*rootd/(Mmax*rootd);
        # groundwater loss via percolation:
        G <- mu.th*alph/(1+mu.th)*M[12,cyear-1]*rootd;
        # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
        R <- P[t,cyear]*(M[12,cyear-1]*rootd/(Mmax*rootd))^m.th +
          (alph/(1+mu.th))*M[12,cyear-1]*rootd;
        dWdt <- P[t,cyear] - Etrans - R - G;
        M[t,cyear] <- M[12,cyear-1] + dWdt/rootd;
      }
      
      if (t == 1 && cyear == 1){
        if (M0 < 0){ M0 <- .20;}
        # evapotranspiration (take initial soil moisture value to be 200 mm)
        Etrans <- Ep*M0*rootd/(Mmax*rootd);
        # groundwater loss via percolation:
        G <- mu.th*alph/(1+mu.th)*(M0*rootd);
        # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
        R <- P[t,cyear]*(M0*rootd/(Mmax*rootd))^m.th + (alph/(1+mu.th))*M0*rootd;
        dWdt <- P[t,cyear] - Etrans - R - G;
        M[t,cyear] <- M0 + dWdt/rootd;
      }
      
      # error-catching:
      if (M[t,cyear] <= Mmin) {M[t,cyear] <- Mmin;}
      if (M[t,cyear] >= Mmax) {M[t,cyear] <- Mmax;}
      if (is.na(M[t,cyear])==1){ M[t,cyear] <- Mmin;}
    } # end month (t) cycle
    #########################################################################################
  } # end year cycle
  
  return(M)
  
}
