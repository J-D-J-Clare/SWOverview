
Cam_Loc_Summary_Env<-read.csv("GBM_DATA.csv")

library(dismo)

gbm_SKUNKSTRIPED_MS<-gbm.step(data=Cam_Loc_Summary_Env, 
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
                                      "EdgeDensity5km", "Nightlight5km", "AgIndex5km", "CoreArea5km", "JanLST5km", "JanMaxSnow5km",
                                      "Aspen_Birch5km", "Broadleaf5km", "Coniferous5km", "ConiferousWetland5km", "Cropland5km", "Developed5km",
                                      "Forest5km", "ForestedWetland5km", "Grassland5km", "HardwoodWetland5km", "HerbWetland5km","Oak5km","Pine5km", 
                                      "propJPine5km", "propRPine5km", "propWater5km", "propSMaple5km", "propCHW5km", "propNHW5km", "Wetland5km", "propSpruceFir5km",
                                      "ShannonL35km", "RichnessL35km","IntEVI5km",
                                      "MaxEVI10km","MinEVI10km",
                                      "SOSEVI10km","EOSEVI10km","LOSEVI10km","MaxNDVI10km","MinNDVI10km","SOSNDVI10km","EOSNDVI10km","LOSNDVI10km",
                                      "EdgeDensity10km", "Nightlight10km", "AgIndex10km", "CoreArea10km", "JanLST10km", "JanMaxSnow10km",
                                      "Aspen_Birch10km", "Broadleaf10km", "Coniferous10km", "ConiferousWetland10km", "Cropland10km", "Developed10km",
                                      "Forest10km", "ForestedWetland10km", "Grassland10km", "HardwoodWetland10km", "HerbWetland10km","Oak10km","Pine10km", 
                                      "propJPine10km", "propRPine10km", "propWater10km", "propSMaple10km", "propCHW10km", "propNHW10km", "Wetland10km", "propSpruceFir10km",
                                      "ShannonL310km", "RichnessL310km","IntEVI10km", "trapnights", "X", "Y") 
                              , gbm.y="SKUNKSTRIPEDDetsPA", 
                              family="bernoulli",
                              tree.complexity=5, learning.rate=0.01, bag.fraction=.5)



gbm_OPOSSUM_MS<-gbm.step(data=Cam_Loc_Summary_Env, 
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
                                 "EdgeDensity5km", "Nightlight5km", "AgIndex5km", "CoreArea5km", "JanLST5km", "JanMaxSnow5km",
                                 "Aspen_Birch5km", "Broadleaf5km", "Coniferous5km", "ConiferousWetland5km", "Cropland5km", "Developed5km",
                                 "Forest5km", "ForestedWetland5km", "Grassland5km", "HardwoodWetland5km", "HerbWetland5km","Oak5km","Pine5km", 
                                 "propJPine5km", "propRPine5km", "propWater5km", "propSMaple5km", "propCHW5km", "propNHW5km", "Wetland5km", "propSpruceFir5km",
                                 "ShannonL35km", "RichnessL35km","IntEVI5km",
                                 "MaxEVI10km","MinEVI10km",
                                 "SOSEVI10km","EOSEVI10km","LOSEVI10km","MaxNDVI10km","MinNDVI10km","SOSNDVI10km","EOSNDVI10km","LOSNDVI10km",
                                 "EdgeDensity10km", "Nightlight10km", "AgIndex10km", "CoreArea10km", "JanLST10km", "JanMaxSnow10km",
                                 "Aspen_Birch10km", "Broadleaf10km", "Coniferous10km", "ConiferousWetland10km", "Cropland10km", "Developed10km",
                                 "Forest10km", "ForestedWetland10km", "Grassland10km", "HardwoodWetland10km", "HerbWetland10km","Oak10km","Pine10km", 
                                 "propJPine10km", "propRPine10km", "propWater10km", "propSMaple10km", "propCHW10km", "propNHW10km", "Wetland10km", "propSpruceFir10km",
                                 "ShannonL310km", "RichnessL310km","IntEVI10km", "trapnights", "X", "Y") 
                         , gbm.y="OPOSSUMDetsPA", 
                         family="bernoulli",
                         tree.complexity=5, learning.rate=0.01, bag.fraction=.5)

###example prediction
###Lyrs<-raster::stack('SW.lyr.grd')
###stacked grid of predictors way too large for github, 
###but will be accessible once moved to alternative archive.

pOpossum_MS<-predict(Lyrs, gbm_Opossum_MS, n.trees=gbm_Opossum_MS$gbm.call$best.trees, type='response', filename="Opossum_Distribution.tif", overwrite=TRUE)

plot(pOpossum_MS)



gbm_COYOTE_MS<-gbm.step(data=Cam_Loc_Summary_Env, 
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
                                "EdgeDensity5km", "Nightlight5km", "AgIndex5km", "CoreArea5km", "JanLST5km", "JanMaxSnow5km",
                                "Aspen_Birch5km", "Broadleaf5km", "Coniferous5km", "ConiferousWetland5km", "Cropland5km", "Developed5km",
                                "Forest5km", "ForestedWetland5km", "Grassland5km", "HardwoodWetland5km", "HerbWetland5km","Oak5km","Pine5km", 
                                "propJPine5km", "propRPine5km", "propWater5km", "propSMaple5km", "propCHW5km", "propNHW5km", "Wetland5km", "propSpruceFir5km",
                                "ShannonL35km", "RichnessL35km","IntEVI5km",
                                "MaxEVI10km","MinEVI10km",
                                "SOSEVI10km","EOSEVI10km","LOSEVI10km","MaxNDVI10km","MinNDVI10km","SOSNDVI10km","EOSNDVI10km","LOSNDVI10km",
                                "EdgeDensity10km", "Nightlight10km", "AgIndex10km", "CoreArea10km", "JanLST10km", "JanMaxSnow10km",
                                "Aspen_Birch10km", "Broadleaf10km", "Coniferous10km", "ConiferousWetland10km", "Cropland10km", "Developed10km",
                                "Forest10km", "ForestedWetland10km", "Grassland10km", "HardwoodWetland10km", "HerbWetland10km","Oak10km","Pine10km", 
                                "propJPine10km", "propRPine10km", "propWater10km", "propSMaple10km", "propCHW10km", "propNHW10km", "Wetland10km", "propSpruceFir10km",
                                "ShannonL310km", "RichnessL310km","IntEVI10km", "trapnights", "X", "Y") 
                        , gbm.y="COYOTEDetsPA", 
                        family="bernoulli",
                        tree.complexity=5, learning.rate=0.01, bag.fraction=.5)

###Briefly, creation of spatial layers.
###Original grid size off 500 x 500 m MODIS pixels used directly for EVI, NDVI, snow layers
###Nighttime lights cells originally 1 km resolution, resampled to 500m
###Land Cover proportions use a polygon (the 500m Modis grid) extraction from base layers (Wiscland2.0).
###Same extraction used to generate edge density, Aggregation,  Shannon, and richness metrics.
###Richness and Shannon's Diversity Index calculated based on level 3 land classes. 
###Edge Density, Core Area, Aggregation Index calculated based on level 1 land classes
###Larger proportional aggregations (e.g., 1km, 5km, 10km suffixes above) derived using raster::focal and raster::focalWeight....
###...e.g., fw1km<-focalweight(r, 1000, 'circle); Pine1km<-(Pine, fw1km, mean)
###Landscape metrics at different scales derived by reclipping the landscape (either original 500 m cell, or circular buffers
###surrounding centroid of original cell. Edge Density, Core Area, Aggregation Index are 'landscape stats', i.e., derived 
###(summed, averaged, etc.) across all patches and cover classes
###Richness and Shannon diversity acriss larger extents reflect averages (eg., mean richness within the 500 m pixels in 
###a 5 km buffer surrounding the centroid of the focal cell)

###Variable explanations:
###Response Variables
#COYOTEDetsPA-Coyote detected (1) or not (0)
#OPOSSUMDetsPA-Opossum detected (1) or not (0)
#SKUNKSTRIPEDDetsPA-Striped Skunk detected (1) or not (0)
###Phenology Metrics
###Max, Min, SOS (day growing season starts), LOS (day growing season ends) are parameters derived from double logistic model
###fit to vegetation indices (Enhanced Vegetation Index, Normalized Difference Vegetation Index) calculated from 16-d r. data
###Values here denote averaged estimates from years 2014 through 2017.
#MaxEVI
#MinEVI
#SOSEVI
#EOSEVI
#LOSEVI --Length of Season = EOS-SOS
#MaxNDVI
#MinNDVI
#SOSNDVI
#EOSNDVI
#LOSNDVI
#IntEVI--see Case 4 code for more detailed description.
###Climate/Snow
#JanLST--This is the average land surface--coverted to f--termperature across January from 2001-2017
#JanMaxSnow--Averaged 8day Maximum Snow Extent sampled across January from 2001-2017.
###Nightlight--Raw nighttime-lights reflectance
###Landscape Configuration Metrics (described above; Wiscland 2.0)
#EdgeDensity
#AgIndex
#CoreArea
###Landscape Composition (Wiscland 2.0)
## % Class --What percent each class (note some redudancy...forest include deciduous, which includes aspend_birch...)
#Aspen_Birch
#Broadleaf
#Coniferous
#ConiferousWetland
#Cropland
#Developed
#Forest
#ForestedWetland
#Grassland
#HardwoodWetland
#HerbWetland
#Oak
#Pine
#propJPine
#propRPine
#propWater
#propSMaple
#propCHW
#propNHW
#Wetland
#propSpruceFir
##Composition summaries (WiscLand 2.0)
#ShannonL3
#RichnessL3 
###Other
#trapnights-how long each point sampled
#X scaled X coordinate
#Y scaled Y coordinate

###Definition of a 'site' here is the specific camera location. However, note point-specific attributes calculated using a point-based extraction.
###In other words, 2 cameras 150 m apart that fall within the same Modis pixel will have the same attributes.