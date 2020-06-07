
Behav_Class<-read.csv('Behav_Class.csv')



###model for foraging
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
