# ######### CONDITIONAL PARAMETER SAMPLING SUBROUTINES ###########
# 
# #' gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
# #' T is the matrix of temperature for all years and all months of the previous simulation.
# #' att is the lower bound on the support of the uniform prior distribution for Tt
# #' btt is the upper bound on the support of the uniform prior distribution for Tt
# #'
# #' SETW 6/10/2010
# #' @export
# 
# 
#  Tt_U_aux <- function(Ttcurr,T,To,gM,RW,errorpars,gE,Gterms,att,btt,intwindow,cyrs){
#    
#    # book-keeping matters:
#    Ny <- dim(Gterms)[2]
#    I_0 <- intwindow[1]
#    I_f <- intwindow[2]
#    
#    # Sample proposal from prior distribution:
#    Ttprop <- runif(1,att,btt);
#    
#   # Compute the growth responses that fall out:
#   gTprop <- matrix(NA,12,Ny);
#   gTprop[T<Ttprop] <- 0
#   gTprop[T>To] <- 1
#   gTprop[ T < To & T > Ttprop] <- (T[T < To & T >Ttprop]-Ttprop)/(To-Ttprop)
#   gprop <- diag(as.vector(gE))%*%pmin(gM,gTprop);
#   gcurr <- Gterms;
#   
#   ########## account for variable integration window:
#   if (I_0<0){ # if we include part of the previous year in each year's modeled growth:
#     startmo <- 13+I_0;
#     endmo <- I_f;
#     prevseas <- cbind(as.matrix(rowMeans(gprop[startmo:12,])),
#                       gprop[startmo:12,1:(Ny-1)])
#     gprop <- gprop[1:endmo,]
#     gprop <- rbind(prevseas,gprop)
#     prevseas <- cbind(as.matrix(rowMeans(gcurr[startmo:12,])),
#                       gcurr[startmo:12,1:(Ny-1)]);
#     gcurr <- gcurr[1:endmo,]
#     gcurr = rbind(prevseas,gcurr)
#   }else{ # no inclusion of last year's growth conditions in estimates of this year's growth:
#     startmo <- I_0+1;
#     endmo <- I_f;
#     gprop <- gprop[startmo:endmo,];
#     gcurr <- gcurr[startmo:endmo,];
#   }
#   mu.curr <- mean(colSums(gcurr))
#   sd.curr <- sd(colSums(gcurr))
#   
#   mu.prop <- mean(colSums(gprop))
#   sd.prop <- sd(colSums(gprop))
#   
#   if (length(errorpars) == 1){ # White noise error model:
#     sigma2rw <- errorpars;
#     Wcurr <- sqrt(1-sigma2rw) * (colSums(gcurr[,cyrs]) - mu.curr)/sd.curr
#     Wprop <- sqrt(1-sigma2rw) * (colSums(gprop[,cyrs]) - mu.prop)/sd.prop
#     expcurr <- sum((RW[cyrs] - Wcurr)^2)
#     expprop <- sum((RW[cyrs] - Wprop)^2)
#     HR = exp(-.5*(expprop-expcurr)/sigma2rw);
#   }
#   if (length(errorpars) == 2){ # AR(1) error model:
#     phi1 <- errorpars[1]; 
#     tau2 <- errorpars[2];
#     sigma2rw <- tau2/(1-phi1^2);
#     #
#     iSig <- makeAR1covmat(phi1,tau2,length(cyrs));
#     #
#     Wcurr = sqrt(1-sigma2rw)*t((sum(gcurr[,cyrs])-mu.curr)/sd.curr)
#     Wprop = sqrt(1-sigma2rw)*t((sum(gprop[,cyrs])-mu.prop)/sd.prop)
#     #
#     logLprop = -.5*t(RW[cyrs]-Wprop)%*%iSig%*%(RW[cyrs] - Wprop);
#     logLcurr = -.5*t(RW[cyrs]-Wcurr)%*%iSig%*%(RW[cyrs] - Wcurr);
#     HR = exp(logLprop-logLcurr);
#   }
# 
#   # accept or reject the proposal.
#   if (binornd(1,min(HR,1))==1){
#     Tt  <-  Ttprop;
#   }else{
#     Tt <- Ttcurr;
#   }
# }
# 
# return(Tt)
# }
# ################################################################
# # function [Tt] = Tt_lit_aux(Ttcurr,T,To,gM,RW,errorpars,gE,Gterms,att,btt,slp,int,intwindow,cyrs)
# # # gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
# # # T is the matrix of temperature for all years and all months of the previous simulation.
# # # att is the lower bound on the support of the uniform prior distribution for Tt
# # # btt is the upper bound on the support of the uniform prior distribution for Tt
# # #
# # # SETW 6/10/2010
# # Ny = size(Gterms,2);
# # I_0 = intwindow(1); I_f = intwindow(2);
# # #
# # if 1 # Sample from prior as proposal distribution!
# # Ttprop = slp*betarnd(att,btt)+int;
# # #     upperlim = normcdf(int+slp,Ttcurr,.5);
# # #     lowerlim = normcdf(int,Ttcurr,.5);
# # #     U = unifrnd(lowerlim,upperlim);
# # #     Ttprop = norminv(U,Ttcurr,.5);
# # #
# # gTprop = NaN*ones(12,Ny);
# # gTprop(T<Ttprop) = 0;
# # gTprop(T>To) = 1;
# # gTprop(T<To&T>Ttprop) = (T(T<To&T>Ttprop)-Ttprop)/(To-Ttprop);
# # gprop = diag(gE)*min(gM,gTprop);
# # gcurr = Gterms;
# # #
# # ########## account for variable integration window:
# # if I_0<0; # if we include part of the previous year in each year's modeled growth:
# # startmo = 13+I_0;
# # endmo = I_f;
# # prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
# # gprop = gprop(1:endmo,:);
# # gprop = [prevseas; gprop];
# # prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
# # gcurr = gcurr(1:endmo,:);
# # gcurr = [prevseas; gcurr];
# # else # no inclusion of last year's growth conditions in estimates of this year's growth:
# #   startmo = I_0+1;
# # endmo = I_f;
# # gprop = gprop(startmo:endmo,:);
# # gcurr = gcurr(startmo:endmo,:);
# # end
# # ############
# # if length(errorpars) == 1 # White noise error model:
# # sigma2rw = errorpars;
# # expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
# # expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
# # HR = exp(-.5*(expprop-expcurr)/sigma2rw);
# # #
# # #         # must also account for relative prior probs if using "geyer default" sampling:
# # #         liklicurr = betapdf((Ttcurr-int)/slp,att,btt);
# # #         likliprop = betapdf((Ttprop-int)/slp,att,btt);
# # #         HR = (likliprop/liklicurr)*exp(-.5*(expprop-expcurr)/sigma2rw);
# # elseif length(errorpars) == 2 # AR(1) error model:
# # phi1 = errorpars(1); tau2 = errorpars(2);
# # sigma2rw = tau2/(1-phi1^2);
# # #
# # [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
# # #
# # Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
# # Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
# # #
# # logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
# # logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
# # HR = exp(logLprop-logLcurr);
# # #
# # #         # must also account for relative prior probs if using "geyer default" sampling:
# # #         liklicurr = betapdf((Ttcurr-int)/slp,att,btt);
# # #         likliprop = betapdf((Ttprop-int)/slp,att,btt);
# # #         HR = (likliprop/liklicurr)*exp(logLprop-logLcurr);
# # end
# # end
# # # accept or reject the proposal.
# # if binornd(1,min(HR,1))==1
# # Tt = Ttprop;
# # else
# #   Tt = Ttcurr;
# # end
# # end
