
Behav_Class<-read.csv('Behav_Class.csv')


library(mgcv)
###model for foraging
###Only uses complete cases.
m1<-bam(cbind(forage, npics)~TRAIL+te(X, Y, Day, bs=c('cr', 'cr', 'cc'), k=c(5, 5, 20))+s(Camera, bs='re')+
          s(DailyEVI, k=20, bs='cr')+s(DailyResEVI, k=20, bs='cr')+s(IntEVI, k=20, bs='cr')+
          s(DailySnowDepth, k=20, bs='cr')+
          s(Nightlight, k=20, bs='cr')+
          s(EdgeDensity, k=20, bs='cr')+
          s(Grassland, k=20, bs='cr')+
          s(Forest, k=20, bs='cr')+
          s(RichnessL35km, k=20, bs='cr')+
          s(Cropland, k=20, bs='cr'), 
        data=Behav_Class, family=binomial(), discrete=TRUE)


###model for vigilance
###Only uses complete cases.

m2<-bam(cbind(vigilant, npics)~TRAIL+te(X, Y, Day, bs=c('cr', 'cr', 'cc'), k=c(5, 5, 20))+s(Camera, bs='re')+
          s(DailyEVI, k=20, bs='cr')+s(DailyResEVI, k=20, bs='cr')+s(IntEVI, k=20, bs='cr')+
          s(DailySnowDepth, k=20, bs='cr')+
          s(Nightlight, k=20, bs='cr')+
          s(EdgeDensity, k=20, bs='cr')+
          s(Grassland, k=20, bs='cr')+
          s(Forest, k=20, bs='cr')+
          s(RichnessL35km, k=20, bs='cr')+
          s(Cropland, k=20, bs='cr'), 
        data=Behav_Class, family=binomial(), discrete=TRUE)


###here, description of the creation of DailyEVI, DailyResEVI, IntEVI variables.
###Recall, double-logistic function used to estimate phenology parameters based upon calculated Veg Indices:
###Minimum VI, Maximum VI, SOS (start of rise), EOS (start of descent), RISE (slope of rise), FALL (slope of descent).
###*DailyEVI* is a predicted (smoothed) estimate of the EVI for a given pixel on a given day using the parameter estimates
###And the double logistic estimating equation (Beck et al. 2006).
###*IntEVI* is the sum of all DailyEVI values for a specific pixel b/w the estimated SOS and EOS.
###*DailyResEVI* is the difference between DailyEVI at a given pixel, and the concurrent mean Daily EVI of the the 8 surrounding pixels.

###Daily Snow Depth derived from Snodas layers--resampled to 500 m grid to align with existing pixels.


