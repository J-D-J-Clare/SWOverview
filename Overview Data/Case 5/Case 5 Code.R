
library(dismo)
read.csv("DeerDetects2018.csv")

gbm_deer_MS_2018<-gbm.step(data=DeerDetects2018GBM_SP@data, 
                           gbm.x=c("MaxEVI","MinEVI",
                                   "SOSEVI","EOSEVI","LOSEVI","MaxNDVI","MinNDVI","SOSNDVI","EOSNDVI","LOSNDVI",
                                   "EdgeDensity", "Nightlight", "AgIndex", "CoreArea", "JanLST", "JanMaxSnow",
                                   "Aspen_Birch", "Broadleaf", "Coniferous", "ConiferousWetland", "Cropland", "Developed",
                                   "Forest", "ForestedWetland", "Grassland", "HardwoodWetland", "HerbWetland","Oak","Pine", 
                                   "propJPine", "propRPine", "propWater", "propSMaple", "propCHW", "propNHW", "Wetland", "propSpruceFir",
                                   "ShannonL3", "RichnessL3","IntEVI",
                                   "MaxEVI1km","MinEVI1km",
                                   "SOSEVI1km","EOSEVI1km","LOSEVI1km","MaxNDVI1km","MinNDVI1km","SOSNDVI1km","EOSNDVI1km","LOSNDVI1km",
                                   "EdgeDensity1km", "Nightlight1km", "AgIndex1km", "CoreArea1km", "JanLST1km", "JanMaxSnow1km",
                                   "Aspen_Birch1km", "Broadleaf1km", "Coniferous1km", "ConiferousWetland1km", "Cropland1km", "Developed1km",
                                   "Forest1km", "ForestedWetland1km", "Grassland1km", "HardwoodWetland1km", "HerbWetland1km","Oak1km","Pine1km", 
                                   "propJPine1km", "propRPine1km", "propWater1km", "propSMaple1km", "propCHW1km", "propNHW1km", "Wetland1km", "propSpruceFir1km",
                                   "ShannonL31km", "RichnessL31km","IntEVI1km",
                                   "MaxEVI5km","MinEVI5km",
                                   "SOSEVI5km","EOSEVI5km","LOSEVI5km","MaxNDVI5km","MinNDVI5km","SOSNDVI5km","EOSNDVI5km","LOSNDVI5km",
                                   "EdgeDensity5km", "Nightlight5km", "CoreArea5km", "JanLST5km", "JanMaxSnow5km",
                                   "Aspen_Birch5km", "Broadleaf5km", "Coniferous5km", "ConiferousWetland5km", "Cropland5km", "Developed5km",
                                   "Forest5km", "ForestedWetland5km", "Grassland5km", "HardwoodWetland5km", "HerbWetland5km","Oak5km","Pine5km", 
                                   "propJPine5km", "propRPine5km", "propWater5km", "propSMaple5km", "propCHW5km", "propNHW5km", "Wetland5km", "propSpruceFir5km",
                                   "ShannonL35km", "RichnessL35km","IntEVI5km",
                                   "MaxEVI10km","MinEVI10km",
                                   "SOSEVI10km","EOSEVI10km","LOSEVI10km","MaxNDVI10km","MinNDVI10km","SOSNDVI10km","EOSNDVI10km","LOSNDVI10km",
                                   "EdgeDensity10km", "Nightlight10km", "CoreArea10km", "JanLST10km", "JanMaxSnow10km",
                                   "Aspen_Birch10km", "Broadleaf10km", "Coniferous10km", "ConiferousWetland10km", "Cropland10km", "Developed10km",
                                   "Forest10km", "ForestedWetland10km", "Grassland10km", "HardwoodWetland10km", "HerbWetland10km","Oak10km","Pine10km", 
                                   "propJPine10km", "propRPine10km", "propWater10km", "propSMaple10km", "propCHW10km", "propNHW10km", "Wetland10km", "propSpruceFir10km",
                                   "ShannonL310km", "RichnessL310km","IntEVI10km", "X", "Y"), 
                           gbm.y="CPUE", 
                           family="gaussian",    tree.complexity=5, learning.rate=0.01, bag.fraction=.5)


Lyrs<-raster::stack('SW.lyr.grd')
DeerCPUE_MS2018<-predict(Lyrs, gbm_deer_MS_2018, n.trees=gbm_deer_MS_2018$gbm.call$best.trees, type='response', filename="DEER_CPUE_MS2018.tif", overwrite=TRUE)


library(maptools)
WI_DMU<-readShapePoly("DMU.shp")
WI_DMU<-extract(DeerCPUE_MS2018, WI_DMU, fun=mean, na.rm=T, sp=T)
cor(WI_DMU$DeerDens18[WI_DMU$DeerDens18!=-999], WI_DMU$DEER_CPUE_MS2018[WI_DMU$DeerDens18!=-999]) ###-999 coded for missing data
