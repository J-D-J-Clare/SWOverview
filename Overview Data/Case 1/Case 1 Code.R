###Not run--how certain data created...
###Cell_Data<-read.csv('CELL_DATA.csv')
##library(fields)
#DummyX<-Cell_Data$X/100000
#DummyY<-Cell_Data$Y/100000
##knots<-cover.design(cbind(DummyX, DummyY), 30)
##omega.all <- (rdist(knots$design, knots$design)/5)^3
##svd.omega.all <- svd(omega.all)
##sqrt.omega.all <- t(svd.omega.all$v %*% (t(svd.omega.all$u)*sqrt(svd.omega.all$d)))
##Z.k <- (rdist(cbind(DummyX, DummyY), knots$design)/5)^3
##Z.matrix <- t(solve(sqrt.omega.all, t(Z.k))) ##Basis functions below have been re-transposed

## Generate neighborhood data and pairwise differences
## librray(spdep)
#CelLNB<-dnearneigh(as.matrix(Cell_Data[, c(65, 66)]), 0, 10000)
#CELL_NB_WB<-nb2WB(CelLNB)

#Neighinfo<-matrix(NA, 6050, 12)
#for (i in 1:6050){
#  Neighinfo[i, 1:CELL_NB_WB$num[i]]<-as.vector(CelLNB[[i]])
#}

#Cell_Data$Forest<-Cell_Data$NLCDConif+Cell_Data$NLCDDecid+Cell_Data$NLCDMixed

#Neighinfo<-matrix(NA, 6050, 12)
#for (i in 1:6050){
#  Neighinfo[i, 1:CELL_NB_WB$num[i]]<-as.vector(CelLNB[[i]])
#}


#diff.cov<-matrix(NA, 6050, 12)
#for (i in 1:6050){
#  diff.cov[i, 1:CELL_NB_WB$num[i]]<-as.vector(Cell_Data_For_RSE2$Forest[i]-Cell_Data_For_RSE2$Forest[CelLNB[[i]]])
#}

#diff.cov2<-matrix(NA, 6050, 12)
#for (i in 1:6050){
#  diff.cov2[i, 1:CELL_NB_WB$num[i]]<-as.vector(Cell_Data_For_RSE2$NLCDWwet[i]-Cell_Data_For_RSE2$NLCDWwet[CelLNB[[i]]])
#}

##Cell_Data$Forest<-scale(Cell_Data$Forest)
##Cell_Data$FW<-scale(Cell_Data$NLCDWwet)




###Read data...
diff.cov<-read.csv('diffusioncovForest.csv') ###the difference in  % Forest between a cell (column) and its neighbors
diff.cov2<-read.csv('diffusioncovFW.csv') ###the difference in  % Forested Wetland between a cell (column) and its neighbors
NN<-read.csv('NN.csv') ###the cell ID of a given cell's neighbors
numn<-read.csv( 'numn.csv') ###how many neighbors a cell has
B.space<-read.csv('BasisFunctions.csv') ###basis function for cells
Forest<-read.csv('Forest.csv') ###scaled % Forest in a cell
FW<-read.csv('FW.csv') ###scaled % Forested Wetland in a cell
Maintained<-read.csv('Maintained.csv') ###Camera placed along maintained trail or not?
Detects<-read.csv('DetectData.csv')
###year, number of 24-hr sampling occassions, number of sampling occassions w/ detections, the cell ID, and the camera ID




###JAGS code

cat("
    model{
    b0_p ~ dnorm (0,0.368) ##Base p(detection)
    b1_p ~ dnorm (0, 0.368) ##Effect of maintained trail
    var.p ~dunif(0, 3) ###variance for RE
    tau.p<-1/var.p^2
    
    
    
    a0_psi~dnorm(0, 0.368) ##intercept
    a1_psi~dnorm(0, 0.368)
    a2_psi~dnorm(0, 0.368)
    
    for (b in 1:30){
    a_space[b] ~ dnorm(0, taub)
    }
    sigb~dunif(0, 5)
    taub<-1/sigb^2
    
    a0.phi~dnorm(0, .368) ##Intercept for persistance
    a1.phi~dnorm(0, 0.368)  ##of forest on persistence
    a2.phi~dnorm(0, 0.368)  ##forested wetland
    a3.phi~dnorm(0, 0.368)  ##autocovariate effect on persistence
    a0.gamma~dnorm(0, 0.368)
    
    a1.gamma~dnorm(0, 0.368)
    a2.gamma~dnorm(0, 0.368)
    a0.dgam~dnorm(0, .368)
    a1.dgam~dnorm(0, .368) ###effect of surrounding diffusion gradient on pr colonization
    a2.dgam~dnorm(0, .368)
    
    for (i in 1:ncells){
    for (n in 1:numn[i]){
    logit(d.gam[i, n])<-a0.dgam+a1.dgam * diff.cov[n, i]+a2.dgam * diff.cov2[n, i]
    }
    }
    
    
    
    ###Likelihood
    for (i in 1:ncells){
    
    logit(psi[i, 1])<-a0_psi+a1_psi*Forest[i]+a2_psi*FW[i]+a_space[1]*B.space[i, 1]+a_space[2]*B.space[i, 2]+a_space[3]*B.space[i, 3]+a_space[4]*B.space[i, 4]+
    a_space[5]*B.space[i, 1]+a_space[6]*B.space[i, 6]+a_space[7]*B.space[i, 7]+a_space[8]*B.space[i, 8]+
    a_space[9]*B.space[i, 9]+a_space[10]*B.space[i, 10]+a_space[11]*B.space[i, 11]+a_space[12]*B.space[i, 1]+
    a_space[13]*B.space[i, 13]+a_space[14]*B.space[i, 14]+a_space[15]*B.space[i, 15]+a_space[16]*B.space[i, 16]+
    a_space[17]*B.space[i, 17]+a_space[18]*B.space[i, 18]+a_space[19]*B.space[i, 19]+a_space[20]*B.space[i, 20]+
    a_space[21]*B.space[i, 21]+a_space[22]*B.space[i, 22]+a_space[23]*B.space[i, 23]+a_space[24]*B.space[i, 24]+
    a_space[25]*B.space[i, 25]+a_space[26]*B.space[i, 26]+a_space[27]*B.space[i, 27]+a_space[28]*B.space[i, 28]+
    a_space[29]*B.space[i, 29]+a_space[30]*B.space[i, 30]
    
    z[i, 1]~dbern(psi[i, 1])
    for (t in 2:nyears){
    
    logit(phi[i, t-1])<-a0.phi+a1.phi*Forest[i]+a2.phi*FW[i]+a3.phi*(sum(z[NN[i,1:numn[i]], t-1])-sum(psi[NN[i,1:numn[i]], t-1]))     ##+eps.phi[i, t-1]
    logit(gamma[i, t-1])<-a0.gamma+a1.gamma*Forest[i]+a2.gamma*FW[i]#*eps.gamma[i, t-1]
    IndNeighOccMat[i, t-1]<-ifelse(sum(z[NN[i,1:numn[i]], t-1]) > 0, 1, 0)
    z.d.vec[i, t-1]<-sum(z[NN[i,1:numn[i]], t-1]*log(1-d.gam[i,1:numn[i]]))
    d.bar[i, t-1]<-1-exp(z.d.vec[i, t-1])
    z[i, t]~dbern(z[i, t-1]*phi[i, t-1]+
    (1-z[i, t-1])*IndNeighOccMat[i, t-1]*(d.bar[i, t-1])+
    (1-z[i, t-1])*(1-IndNeighOccMat[i, t-1])*gamma[i, t-1])
    psi[i, t]<-psi[i, t-1]*phi[i, t-1]+(1-psi[i, t-1])*(d.bar[i, t-1])*(1-gamma[i, t-1])+(1-psi[i, t-1])*gamma[i, t-1]*(1-d.bar[i, t-1])
    }
    }
    
    ###Obs likelihood
    for (t in 1:nyears){
    
    
    for (j in 1:nlocs){
    eps.p[j, t]~dnorm(0, tau.p)
    logit(p[j, t])<-b0_p+b1_p*Maintained[j]+eps.p[j, t]
    }
    }
    
    for (o in 1:nobs){
    mu.p[o]<-z[cell[o], year[o]]*p[cam[o], year[o]]
    y[o]~dbin(mu.p[o], K[o])
    }
    
    }
    ",fill=TRUE, file="beardiffus.jags")


zst<-matrix(1, 6050, 11)
inits =  function() {list(z=zst,a_space=rep(0, 30), 
                          a0.phi=3, a1.phi=0,a2.phi=0,  b0_p=-3,
                          b1_p=0, a0.gamma=-3, a1.gamma=0,
                          a0.dgam=-2, a1.dgam=0)}

params<-c('a0_psi', 'a1_psi', 'a2_psi', 'a_space', 'a0.dgam', 'a1.dgam', 'a2.dgam', 'b0_p', 'b1_p',  'a0.gamma', 'a1.gamma',
          'a2.gamma', 
          'a0.phi', 'a1.phi', 'a2.phi','a3.phi', 'var.p', 'sigb', 'psi', 'z')


###compile data
data<-list(ncells=6050, nyears=11, nobs=nrow(Detects), 
           nlocs=1770, 
           cell=Detects$ID , 
           year=Detects$Year-2013, ###add 2013 for real year
           cam=Detects$CamNum_Redux, 
           y=as.integer(Detects$Ndets), 
           K=as.integer(Detects$Ndays),
           Maintained=as.integer(unlist(Maintained)),
           diff.cov=as.matrix(diff.cov),
           diff.cov2=t(as.matrix(diff.cov2)),
           NN=as.matrix(NN),
           B.space=as.matrix(B.space),
           numn=as.integer(unlist(numn)),
           Forest=as.numeric(unlist(Forest)),
           FW=as.numeric(unlist(FW)))


library(jagsUI)

###fit model...
fitbearAll_2no_arriv2= jags(data=data, inits=inits, params, model.file="beardiffus.jags",
                            n.chains=4, n.iter=75000, n.adapt=5000, n.burnin=10000, n.thin=50, parallel=T)



