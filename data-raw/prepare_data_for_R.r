


tmp <- R.matlab::readMat('~/vslite_testdata.mat')
sitecoords <- tmp$sitecoords
syear <- tmp$syear
eyear <- tmp$eyear
trw.obs <- tmp$trw.obs
Tmp <- P <- array(NA,dim = c(12,length(syear:eyear),4))
climate <- tmp$m08clim
for(i in 1:4){
  Tmp[,,i] <- climate[,,i]$T
  P[,,i] <- climate[,,i]$P
}

devtools::use_data(Tmp)
devtools::use_data(P)
devtools::use_data(syear)
devtools::use_data(eyear)
devtools::use_data(sitecoords)
devtools::use_data(trw.obs)

