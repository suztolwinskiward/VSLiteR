#' Run CPC Leaky Bucket model (submonthly version).
#' 
#' \code{leakybucket.submonthly} simulates soil moisture with submonthly time step.
#' 
#' Modifications of Suz Tolwinski-Ward's monthly time-step code by Nick Graham in 2011.
#' Ported to R by SETW in 2015. Implementation of CPC Leaky Bucket model as described in
#' Huang et al., 'Analysis of Model-Calculated Soil Moisture over the United States 
#' (1931-1993) and Applications to Long-Range Temperature Forecasts,' J. Clim. (1995)
#' #' Fixes to model consistent with Matlab VSLite version 2.5. 
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
#' @seealso \code{\link{leakybucket.monthly}}
#' 
#' @export

## LEAKY BUCKET WITH SUBSTEPPING ##
leakybucket.submonthly <- function(syear,eyear,phi,T,P,Mmax = 0.76,Mmin = 0.01,alph = 0.093,
                                   m.th = 4.886,mu.th = 5.8,rootd = 1000,M0 = .2){
  
  #  [M,potEv,ndl,cdays] <- leakybucket_submonthly(syear,eyear,phi,T,P,...
  #                                                      Mmax,Mmin,alph,m_th,mu_th,rootd,M0)
  
  # leackybucket_submonthly.m - Simulate soil moisture; substeps within monthly timesteps
  # to better capture nonlinearities and improve moisture estimates.
  ####################################################################################################
  # Usage: [M,potEv,ndl,cdays] <- leakybucket_submonthly(syear,eyear,phi,T,P,...
  #                   Mmax,Mmin,alph,m_th,mu_th,rootd,M0)
  #    outputs simulated soil moisture and potential evapotranspiration.
  #
  # Inputs:
  #   syear <- start year of simulation.
  #   eyear <- end year of simulation.
  #   phi <- latitude of site (in degrees N)
  #   T <- (12 x Nyrs) matrix of ordered mean monthly temperatures (in degEes C)
  #   P <- (12 x Nyrs) matrix of ordered accumulated monthly precipitation (in mm)
  #   Mmax <- scalar maximum soil moisture held by the soil (in v/v)
  #   Mmin <- scalar minimum soil moisture (for error-catching) (in v/v)
  #   alph <- scalar runoff parameter 1 (in inverse months)
  #   m_th <- scalar runoff parameter 3 (unitless)
  #   mu_th <- scalar runoff parameter 2 (unitless)
  #   rootd <- scalar root/"bucket" depth (in mm)
  #   M0 <- initial value for previous month's soil moisture at t <- 1 (in v/v)
  #
  # Outputs:
  #   M <- soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs)
  #   potEv <- potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm)
  #
  # SETW+ N. Graham and K. Georgakakos 2011
  
  # modified by Nick G. and K. Georgakakos - to sub-step the monthly steps. Also this version has added
  # soil moisture initial conditions for restarts, or spin-up.  Hands back monthly soil moisture
  # and summer soil moisture as well - see varargout.  Nick G. 2011/06
  
  ####################################################################################################
  iyear <- syear:eyear;
  nyrs <- length(iyear);
  
  # Storage for growth response output variables (size [12 x Nyears]):
  M <- potEv <- matrix(NA,12,nyrs);
  
  if(M0 < 0.){M0 <- 200/rootd;}
  
  # Compute normalized daylength (neglecting small difference in calculation for leap-years)
  L <- daylength.factor.from.lat(phi); 
  
  # Pre-calculation of istar and I, using input T to compute the climatology:
  Tm <- colMeans(T);
  istar <- (Tm/5)^1.514;
  istar[Tm < 0] <- 0;
  I <- sum(istar);

  # precalculation of the exponent alpha in the Thornwaite (1948) equation:
  a <- (6.75e-7)*I^3 - (7.71e-5)*I^2 + (1.79e-2)*I + 0.49;

  #########################################################################################
  #### -- year cycle -- ####
  # syear <- start (first) year of simulation
  # eyear <- end (last) year of simulation
  # cyear <- year the model is currently working on
  # iyear <- index of simulation year
  
  for (cyear in 1:length(iyear)) {      # begin cycling over years
    #########################################################################################
    for (t in 1:12){  # begin cycling over months in a year
      
      ##### Compute potential evapotranspiration for current month after Thornthwaite:
      if (T[t,cyear] < 0){Ep <- 0;}
      if (T(t,cyear) >= 0 && T(t,cyear) < 26.5){Ep <- 16*L[t]*(10*T[t,cyear]/I)^a;}
      if (T[t,cyear] >= 26.5){Ep <- -415.85 + 32.25*T[t,cyear] - .43* T[t,cyear]^2;}
      potEv[t,cyear] <- Ep;
      
      ##### Now calculate soil moisture according to the CPC Leaky Bucket model
      ##### (see J. Huang et al, 1996). Set n-steps according to 2 mm increments
      ##### have to update alpha and Ep as well - 2 mm increments came from
      ##### testing by K. Georgakakos, but one could use 5 or more, with less "accurate" results.
      ##### Stepping is necessary because the parametization is linearized around init condition.
      #################
      
      dp <- 2.0; # mm of precip per increment
      nstep <- floor(P[t,cyear]/dp)+1; # number of sub-monthly substeps
      Pinc <- P[t,cyear]/nstep; # precip per substep
      alphinc <- alph/nstep; # runoff rate per substep time interval
      Epinc <- Ep/nstep; # potential evapotrans per substep.
      #################
      
      # handling for sm_init
      if (t > 1){M0 <- M[t-1,cyear]}
      if (t == 1 && cyear > 1){M0 <- M[12,cyear-1]}
      sm0 <- M0;
      for(istep in 1:nstep){
        # evapotranspiration:
        Etrans <- Epinc*sm0*rootd/(Mmax*rootd);
        # groundwater loss via percolation:
        G <- mu.th*alphinc/(1+mu.th)*sm0*rootd;
        # runoff; contributions from surface flow (1st term) and subsurface (2nd term)
        R <- Pinc*(sm0*rootd/(Mmax*rootd))^m.th + (alphinc/(1+mu.th))*sm0*rootd;
        dWdt <- Pinc - Etrans - R - G;
        sm1 <- sm0 + dWdt/rootd;
        #
        sm0 <- max(sm1,Mmin);
        sm0 <- min(sm0,Mmax);
      }
      M[t,cyear] <- sm0;
      
      # error-catching:
      if (M[t,cyear] <= Mmin){M[t,cyear] <- Mmin;}
      if (M[t,cyear] >= Mmax){M(t,cyear) <- Mmax;}
      if (is.na(M[t,cyear])==1){M[t,cyear] <- Mmin;}
    
    } # end month (t) cycle
  } # end year cycle

  return(M)

}
