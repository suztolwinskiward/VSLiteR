#' VS-Lite model of tree ring width growth.
#' 
#' \code{VSLiteR} simulates tree ring width growth.
#' 
#' R port of VS-Lite Model of Tree Ring Width by Suz TOlwinski-Ward, 2015. For more references,
#' see xxxxyyyyyzzzz.
#' 
#' @param syear Start year of simulation.
#' @param eyear End year of simulation.
#' @param phi Latitude of site (in degrees N).
#' @param T (12 x Nyrs) Matrix of ordered mean monthly temperatures (in degEes C).
#' @param P (12 x Nyrs) Matrix of ordered accumulated monthly precipitation (in mm).
#' @param T1 Lower temperature threshold for growth to begin (scalar, deg. C).
#' @param T1 Upper temperature threshold for growth sensitivity to temp (scalar, deg. C).
#' @param M1 Lower moisture threshold for growth to begin (scalar, v.v).
#' @param M2 Upper moisture threshold for growth sensitivity to moisture (scalar, v/v).
#' @param Mmax Scalar maximum soil moisture held by the soil (in v/v).
#' @param Mmin Scalar minimum soil moisture (for error-catching) (in v/v).
#' @param alph Scalar runoff parameter 1 (in inverse months).
#' @param m.th Scalar runoff parameter 3 (unitless).
#' @param mu.th Scalar runoff parameter 2 (unitless).
#' @param rootd Scalar root/"bucket" depth (in mm).
#' @param M0 Initial value for previous month's soil moisture at t = 1 (in v/v).
#' @param subtep 
#' @param I_0 
#' @param I_f
#' @param hydroclim Switch; value is either "P" (default) or "M" depending on whether the 
#' second input climate variable is precipitation, in which case soil moisture is estimated
#' using the Leaky Bucket model of the CPC, or soil moisture, in which case the inputs are 
#' used directly to compute the growth response.
#' 
#' @return trw
#' @return gT
#' @return gM
#' @return gE
#' @return M
#' @return potEv
#' @return sample.mean.width 
#' @return sample.std.width 
#' 
#' @seealso \code{\link{compute.gE}},\code{\link{std.ramp}},\code{\link{leakybucket.monthly}},\code{\link{leakybucket.submonthly}}
#'
#' @export
####################################################################################################


VSLiteR <- function(syear,eyear,phi,T1,T2,M1,M2,T,P,varargin,
                        Mmax = 0.76,Mmin = 0.01,alph = 0.093,
                        m.th = 4.886,mu.th = 5.8,rootd = 1000,M0 = .2,
                        substep = 0,I_0 = 1,I_f = 12,hydroclim = "P"){
  #############################################################################
  nyrs <- length(syear:eyear)
  Gr <- gT <- gM <- M <- potEv <- matrix(NA,12,nyrs);
  #############################################################################
  
  ## Load in soil moisture, or estimate it with the Leaky Bucket model:
  if(hydroclim == "M"){
    ## Read in soil moisture:
    M = P;
  }else{# Compute soil moisture:
    if(substep == 1){
      M <- leakybucket_submonthly(syear,eyear,phi,T,P,
                                  Mmax,Mmin,alph,m.th,mu.th,rootd,M0);
    }else{
      M <- leakybucket.monthly(syear,eyear,phi,T,P,
                               Mmax,Mmin,alph,m.th,mu.th,rootd,M0);
    }
    if(substep !=1 && substep != 0){
      cat("'substep' param must either be set to 1 or 0.");
      return
    }
  }
  
  # Compute gE, the scaled monthly proxy for insolation:
  gE <- compute.gE(phi);
  
  #############################################################################
  ### Calculate Growth Response functions gT and gM
  
  # Temperature growth response:
  gT <- std.ramp(T,T1,T2)
  
  # Soil moisture growth response:
  gM <- std.ramp(M,M1,M2)
  
  # Compute overall growth rate:
  Gr <- kronecker(matrix(1,1,nyrs),gE)*pmin(gT,gM)
 
  ############## Compute proxy quantity from growth responses #################
  width <- matrix(NA,nyrs,1);
  if (phi>0){ # Site in Northern Hemisphere:
    if (I_0<0){ # if we include part of the previous year in each year's modeled growth:
      startmo <- 13+I_0;
      endmo <- I_f;
      # use average of growth data across modeled years to estimate first year's growth due
      # to previous year:
      width[1] <- sum(Gr[1:endmo,1]) + sum(rowMeans(Gr[startmo:12,]));
      for(cyear in 2:nyrs){
        width[cyear] <- colSums(Gr[startmo:12,cyear-1]) + colSums(Gr[1:endmo,cyear]);
      }
    }else{ # no inclusion of last year's growth conditions in estimates of this year's growth:
      startmo <- I_0+1;
      endmo <- I_f;
      width[cyear] <- colSums(Gr[startmo:endmo,])
    }
  }
  if(phi<0){ # if site is in the Southern Hemisphere:
    # (Note: in the Southern Hemisphere, ring widths are dated to the year in which growth began!)
    startmo <- 7+I_0; # (eg. I_0 = -4 in SH corresponds to starting integration in March of cyear)
    endmo <- I_f-6; # (eg. I_f = 12 in SH corresponds to ending integraion in June of next year)
    for (cyear in 1:(nyrs-1)){
      width(cyear) <- sum(Gr[startmo:12,cyear]) + sum(Gr[1:endmo,cyear+1]);
    }
    # use average of growth data across modeled years to estimate last year's growth due
    # to the next year:
    width[nyrs] <- sum(Gr[startmo:12,nyrs])+sum(rowMeans(Gr[1:endmo,]));
  }
  
  # Simulated proxy series standardized width:
  trw <- t((width-mean(width))/std(width)); 

  #############################################################################
  # Return output:
  out <- list(trw = trw, gT = gT, gM = gM, gE = gE, M = M, potEv = potEv,
              sample.mean.width = mean(width), sample.std.width = sd(width))
  return(out)

}


