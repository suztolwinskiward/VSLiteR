#' Test port of VSLite to R, and draw some figures.
#' 
#' Test VSLiteR code.  
#' 
#' @param site Index 1, 2, 3, or 4, specifies which of 4 test sites to run.
#
#' @seealso \code{\link{VSLiteR}}
#' 
#' @export

test_VSLiteR <- function(site = 1){
  # load the test data:
  data(sitecoords,Tmp,P,syear,eyear,trw.obs)
  
  # select climate and location data for the chosen test site:
  T <- Tmp[,,site];
  P <- P[,,site];
  phi <- sitecoords[site,1];
  
  #  # estimate the climate response parameters:
  #  cat('Performing Bayesian estimation of VS-Lite parameters for chosen site.')
  # [T1,T2,M1,M2] = estimate_vslite_params_v2_3(T,P,phi,trw_obs(:,site)','nsamp',2000);
  
  # Run VS-Lite.
  out <- VSLiteR(syear,eyear,phi,T,P); # all other parameters: use defaults.
  
  # Now make some plots of the output:
  par(mfrow=c(2,2))
  plot(rowMeans(out$gM),col='blue',type="o",ylim=c(0,1.1),
       xlab="months",ylab="Growth responses");
  lines(rowMeans(out$gT),col='red',type="o");
  title('Mean gT (red) and gM (blue)')
  #
  maps::map("state")
  points(sitecoords[site,2],sitecoords[site,1],col="red",pch = "*")
  title('Site location')
  #
  plot(syear:eyear,out$trw,col='red',type="o",
       ylim=c(1.1*min(c(out$trw,trw.obs[,site])),
              1.1*max(c(out$trw,trw.obs[,site]))));
  lines(syear:eyear,trw.obs[,site]);
  title('Simulated RW (red) and observed (black)');
  
}
