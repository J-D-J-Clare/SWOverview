###Load needed packages
library(scrbook)
library(raster)
library(boot)
library(gstat)
library(matrixcalc)
library(msm)
library(emdbook)
library(spdep)

###replicate simulation study....first, some likelihood functions to pass to nlm/optim
###great credit due to maintainers of the scrbook package--scr likelihood presented
###essentially follows that package's "intlik3Poisson"

###Beta-binomial likelihood, no weighting
SCROCC_BB<-function (start = NULL, y = y, K = NULL, X = traps, X2=NULL, y2=NULL, K2=NULL,
                          G=gr, z = NULL)
{
  
  require(scrbook)
  require(raster)
  require(boot)
  require(gstat)
  require(matrixcalc)
  require(msm)
  nG <- nrow(G)
  D <- e2dist(X, G)
  
  alpha0 <- start[1]
  sigma <- exp(start[2])
  Dens <- exp(start[3]+start[4]*z)
  theta<-exp(start[5])
  psi<-Dens/sum(Dens)
  loglam <- alpha0 - (1/(2 * sigma * sigma)) * D * D # + alpha2 * ztrap Ignore RSF
  probcap <- 1 - exp(-exp(loglam))
  Pm <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
  Pm2 <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
  ymat <- y
  ymat <- rbind(y, rep(0, ncol(y)))
  lik.marg <- rep(NA, nrow(ymat))
  
  if (!is.null(y2)){
    D_2<-e2dist(X2, G)
    loglamOcc<-alpha0 - (1/(2 * sigma * sigma)) * D_2 * D_2
    lamOcc<-exp(loglamOcc)
    lamOcc2<-matrix(NA, nrow(lamOcc), ncol(lamOcc))
    Pm3<-rep(NA, nrow(lamOcc))
    like.occ<-rep(NA,nrow(lamOcc))
    for (i in 1:nrow(lamOcc)){
      lamOcc2[i, ]<-lamOcc[i, ]*Dens
      Pm3[i]<-1-exp(-sum(lamOcc2[i, ]))
      like.occ[i]<-dbetabinom(y2[i],Pm3[i], K2,theta, log=TRUE)
    }
  }
  for (i in 1:nrow(ymat)) {
    Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K,nG), probcap[1:length(Pm)], log = TRUE))
    lik.cond <- exp(colSums(Pm))
    lik.marg[i] <- sum(lik.cond * psi)
  }
  
  nv<-c(rep(1, length(lik.marg)-1), 1)
  atheta<-1-lik.marg[nrow(ymat)]
  nind<-nrow(ymat)-1
  pixels<-rep(1, nG)
  part1 <- nind * log(sum(Dens*pixels)) - sum(Dens * pixels)*atheta
  part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
  
  if (!is.null(y2)){
    part3<-sum(like.occ)
    out <- -1 * (part1 + part2+part3)}
  else {
    out <- -1 * (part1 + part2)
  }
}




###Weighted beta-binomail pseudolikelihood for occupancy data

SCROCC_Weighted_BB_OSCR<-function (start = NULL, y = y, K = NULL, X = traps, X2=NULL, y2=NULL, K2=NULL,
                                   G=gr, z = NULL, trapweights=NULL)
{
  
  require(scrbook)
  require(raster)
  require(boot)
  require(gstat)
  require(matrixcalc)
  require(msm)
  nG <- nrow(G)
  D <- e2dist(X, G)
  
  alpha0 <- start[1]
  sigma <- exp(start[2])
  Dens <- exp(start[3]+start[4]*z)
  theta<-exp(start[5])
  psi<-Dens/sum(Dens)
  loglam <- alpha0 - (1/(2 * sigma * sigma)) * D * D # + alpha2 * ztrap Ignore RSF
  probcap <- 1 - exp(-exp(loglam))
  Pm <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
  Pm2 <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
  ymat <- y
  ymat <- rbind(y, rep(0, ncol(y)))
  lik.marg <- rep(NA, nrow(ymat))
  
  if (!is.null(y2)){
    D_2<-e2dist(X2, G)
    loglamOcc<-alpha0 - (1/(2 * sigma * sigma)) * D_2 * D_2
    lamOcc<-exp(loglamOcc)
    lamOcc2<-matrix(NA, nrow(lamOcc), ncol(lamOcc))
    Pm3<-rep(NA, nrow(lamOcc))
    like.occ<-rep(NA,nrow(lamOcc))
    if (is.null(trapweights)){
      trapweights<-rep(1, nrow(lamOcc))
    }
    for (i in 1:nrow(lamOcc)){
      lamOcc2[i, ]<-lamOcc[i, ]*Dens
      Pm3[i]<-1-exp(-sum(lamOcc2[i, ]))
      like.occ[i]<-trapweights[i]*dbetabinom(y2[i],Pm3[i], K2,theta, log=TRUE)
    }
  }
  for (i in 1:nrow(ymat)) {
    Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K,nG), probcap[1:length(Pm)], log = TRUE))
    lik.cond <- exp(colSums(Pm))
    lik.marg[i] <- sum(lik.cond * psi)
  }
  
  nv<-c(rep(1, length(lik.marg)-1), 1)
  atheta<-1-lik.marg[nrow(ymat)]
  nind<-nrow(ymat)-1
  pixels<-rep(1, nG)
  part1 <- nind * log(sum(Dens*pixels)) - sum(Dens * pixels)*atheta
  part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
  
  if (!is.null(y2)){
    part3<-sum(like.occ)
    out <- -1 * (part1 + part2+part3)}
  else {
    out <- -1 * (part1 + part2)
  }
}


###not considered in the text, but for comparison (if interested):
###likelihood below ignores any overdispersion in the occupancy data.

SCROCC_B<-function (start = NULL, y = y, K = NULL, X = traps, X2=NULL, y2=NULL, K2=NULL,
                         G=gr, z = NULL)
{
  
  require(scrbook)
  require(raster)
  require(boot)
  require(gstat)
  require(matrixcalc)
  require(msm)
  nG <- nrow(G)
  D <- e2dist(X, G)
  
  alpha0 <- start[1]
  sigma <- exp(start[2])
  Dens <- exp(start[3]+start[4]*z)
  #theta<-exp(start[5])
  psi<-Dens/sum(Dens)
  loglam <- alpha0 - (1/(2 * sigma * sigma)) * D * D # + alpha2 * ztrap Ignore RSF
  probcap <- 1 - exp(-exp(loglam))
  Pm <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
  Pm2 <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
  ymat <- y
  ymat <- rbind(y, rep(0, ncol(y)))
  lik.marg <- rep(NA, nrow(ymat))
  
  if (!is.null(y2)){
    D_2<-e2dist(X2, G)
    loglamOcc<-alpha0 - (1/(2 * sigma * sigma)) * D_2 * D_2
    lamOcc<-exp(loglamOcc)
    lamOcc2<-matrix(NA, nrow(lamOcc), ncol(lamOcc))
    Pm3<-rep(NA, nrow(lamOcc))
    like.occ<-rep(NA,nrow(lamOcc))
    for (i in 1:nrow(lamOcc)){
      lamOcc2[i, ]<-lamOcc[i, ]*Dens
      Pm3[i]<-1-exp(-sum(lamOcc2[i, ]))
      like.occ[i]<-dbinom(y2[i], K2,Pm3[i],log=TRUE)
    }
  }
  for (i in 1:nrow(ymat)) {
    Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K,nG), probcap[1:length(Pm)], log = TRUE))
    lik.cond <- exp(colSums(Pm))
    lik.marg[i] <- sum(lik.cond * psi)
  }
  
  nv<-c(rep(1, length(lik.marg)-1), 1)
  atheta<-1-lik.marg[nrow(ymat)]
  nind<-nrow(ymat)-1
  pixels<-rep(1, nG)
  part1 <- nind * log(sum(Dens*pixels)) - sum(Dens * pixels)*atheta
  part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
  
  if (!is.null(y2)){
    part3<-sum(like.occ)
    out <- -1 * (part1 + part2+part3)}
  else {
    out <- -1 * (part1 + part2)
  }
}


###simulation function...

Sim_fit_SCR_OCC_joint2mod<-function(nsims=10, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                                    X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                                    X2=cbind( sort(rep( seq(10,90, 10),9)), rep( seq(10,90,10),9)))
{
  gr<-expand.grid(1:100,1:100)
  Dmat<-as.matrix(dist(gr))
  names(gr)<-c("x", "y")
  
  
  ntraps1<-nrow(X1)
  ntraps2<-nrow(X2)
  
  
  SCR_nind<-rep(NA, nsims)
  SCR_ndet<-rep(NA, nsims)
  Occ_ndet<-rep(NA, nsims)
  
  JointA0<-matrix(NA, nsims, 4)
  JointSigma<-matrix(NA, nsims, 4)
  JointB0<-matrix(NA, nsims, 4)
  JointB1<-matrix(NA, nsims, 4)
  JointN<-matrix(NA, nsims, 3)
  JointChat<-rep(NA, nsims)
  varCovJoint<-array(NA, dim=c(2,2, nsims))
  set.seed(6000)
  g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=0, model=vgm(psill=0.25,model="Exp",range=10), nmax=20)
  z<- predict(g.dummy, newdata=gr, nsim=1)
  z<-z[,3]
  probs2<-exp(beta1*z)/sum(exp(beta1*z))
  
  for (s in 1:nsims){
    
    Sid<- sample(1:10000,N,replace=TRUE, prob=probs2)
    S<-gr[Sid,]
    D1<- e2dist(S,X1) ## N x ntraps
    D2<- e2dist(S,X2)
    loglam1<- alpha0 -(1/(2*sigma*sigma))*D1*D1 #+ alpha2*Zmat1 
    p1<- 1-exp(-exp(loglam1))
    loglam2<- alpha0 -(1/(2*sigma*sigma))*D2*D2 #+ alpha2*Zmat2 
    loglamOcc<-log(colSums(exp(loglam2))) 
    p2<- 1-exp(-exp(loglamOcc))
    y1<-matrix(NA,nrow=N,ncol=ntraps1)
    y2<-rep(NA,ntraps2)
    for (i in 1:ntraps2){ y2[i]<- rbinom(1,K,p2[i])}
    Occ_ndet[s]<-sum(y2)
    for(i in 1:N){
      y1[i,]<- rbinom(ntraps1,K,p1[i,])
    }
    ###which guys seen?
    cap1<-apply(y1,1,sum)>0
    y1<-y1[cap1,]
    SCR_nind[s]<-nrow(y1)
    SCR_ndet[s]<-sum(y1)
    ###fitting the BB likelihood without any trap specific weights
    tmp2<-nlm(SCROCC_BB,c(-1,0,-3, 0, 2),y=y1,K=K,X=X1,G=gr, z=z, y2=y2, K2=K,
              X2=X2, hessian = T)
    if (is.non.singular.matrix(tmp2$hessian)){
      fisher_info<-solve(tmp2$hessian)
      varCovJoint[,,s]<-fisher_info[3:4, 3:4]
      prop_sigma<-sqrt(diag(fisher_info))
      upper<-tmp2$estimate+1.96*prop_sigma
      lower<-tmp2$estimate-1.96*prop_sigma
      interval<-data.frame(value=tmp2$estimate, SE=prop_sigma, lower=lower, upper=upper)
      D_Pred<-tmp2$estimate[3]+tmp2$estimate[4]*z
      
    }
    
    else {
      interval<-data.frame(value=tmp2$estimate, SE=rep(NA, 4), lower=rep(NA, 4), upper=rep(NA, 4))  
      D_Pred<-tmp2$estimate[3]+tmp2$estimate[4]*z
      
    } 
    
    JointA0[s,]<-as.numeric(interval[1,])
    JointSigma[s,]<-as.numeric(interval[2,])
    JointB0[s,]<-as.numeric(interval[3,])
    JointB1[s,]<-as.numeric(interval[4,])
    JointN[s,1]<-sum(exp(D_Pred))
    JointChat[s]<-as.numeric(exp(interval[5,]))
    
  }
  
  out=list(
    SCR_nind=SCR_nind,
    SCR_ndet=SCR_ndet,
    Occ_ndet=Occ_ndet,
    JointA0=JointA0,
    JointSigma=JointSigma,
    JointB0=JointB0,
    JointB1=JointB1,
    JointN=JointN,
    varCovJoint=varCovJoint,
    JointChat=JointChat,
    z=z)
  return(out)
}





sim1<-Sim_fit_SCR_OCC_joint2mod(nsims=300, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                               X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                               X2=cbind( sort(rep( seq(10,90, 10),9)), rep( seq(10,90,10),9)))


sim2<-Sim_fit_SCR_OCC_joint2mod(nsims=300, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                                X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                                X2=cbind( sort(rep( seq(63,87, 3),9)), rep( seq(28,52,3),9)))



sim3<-Sim_fit_SCR_OCC_joint2mod(nsims=300, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                                X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                                X2=cbind( sort(rep(c(seq(10, 16, 3), seq(45, 51, 3), seq(80, 86, 3)), 9)), rep(c(seq(10, 16, 3), seq(45, 51, 3), seq(80, 86, 3)), 9)))



sim4<-Sim_fit_SCR_OCC_joint2mod(nsims=300, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                                X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                                X2=cbind( sort(rep( seq(10,90, 4),21)), rep( seq(10,90,4), 21)))







####simulations with trap-specific weighting
####Note, we ignore any 'sim1w' below b/c the likelihood
####is nearly exactly the same as the unweighted likelihood
###(The trap-specific weights are essentially all 0.995 rather than 1)


Sim_fit_SCR_OCC_joint2modweights<-function(nsims=10, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                                           X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                                           X2=cbind( sort(rep( seq(10,90, 10),9)), rep( seq(10,90,10),9)), trapweights=NULL)
{
  gr<-expand.grid(1:100,1:100)
  Dmat<-as.matrix(dist(gr))
  names(gr)<-c("x", "y")
  
  
  ntraps1<-nrow(X1)
  ntraps2<-nrow(X2)
  
  
  SCR_nind<-rep(NA, nsims)
  SCR_ndet<-rep(NA, nsims)
  Occ_ndet<-rep(NA, nsims)
  
  JointA0<-matrix(NA, nsims, 4)
  JointSigma<-matrix(NA, nsims, 4)
  JointB0<-matrix(NA, nsims, 4)
  JointB1<-matrix(NA, nsims, 4)
  JointN<-matrix(NA, nsims, 3)
  JointChat<-rep(NA, nsims)
  varCovJoint<-array(NA, dim=c(2,2, nsims))
  set.seed(6000)
  g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=0, model=vgm(psill=0.25,model="Exp",range=10), nmax=20)
  z<- predict(g.dummy, newdata=gr, nsim=1)
  z<-z[,3]
  probs2<-exp(beta1*z)/sum(exp(beta1*z))
  
  if (is.null(trapweights)){
    trapweights<-rep(1, nrow(X2))
  }
  
  for (s in 1:nsims){
    
    Sid<- sample(1:10000,N,replace=TRUE, prob=probs2)
    S<-gr[Sid,]
    D1<- e2dist(S,X1) ## N x ntraps
    D2<- e2dist(S,X2)
    
    
    #Zmat1<- matrix(z[raster.point1],nrow=N,ncol=ntraps1,byrow=TRUE) # note make dims the same
    loglam1<- alpha0 -(1/(2*sigma*sigma))*D1*D1 #+ alpha2*Zmat1 
    p1<- 1-exp(-exp(loglam1))
    #Zmat2<- matrix(z[raster.point2],nrow=N,ncol=ntraps2,byrow=TRUE) # note make dims the same
    loglam2<- alpha0 -(1/(2*sigma*sigma))*D2*D2 #+ alpha2*Zmat2 
    loglamOcc<-log(colSums(exp(loglam2))) 
    p2<- 1-exp(-exp(loglamOcc))
    y1<-matrix(NA,nrow=N,ncol=ntraps1)
    y2<-rep(NA,ntraps2)
    for (i in 1:ntraps2){ y2[i]<- rbinom(1,K,p2[i])}
    Occ_ndet[s]<-sum(y2)
    for(i in 1:N){
      y1[i,]<- rbinom(ntraps1,K,p1[i,])
    }
    ###which guys seen?
    cap1<-apply(y1,1,sum)>0
    y1<-y1[cap1,]
    SCR_nind[s]<-nrow(y1)
    SCR_ndet[s]<-sum(y1)
    tmp2<-nlm(SCROCC_Weighted_BB,c(-1,0,-3, 0, 2),y=y1,K=K,X=X1,G=gr, z=z, y2=y2, K2=K,
              X2=X2, hessian = T, trapweights=trapweights)
    if (is.non.singular.matrix(tmp2$hessian)){
      fisher_info<-solve(tmp2$hessian)
      varCovJoint[,,s]<-fisher_info[3:4, 3:4]
      prop_sigma<-sqrt(diag(fisher_info))
      upper<-tmp2$estimate+1.96*prop_sigma
      lower<-tmp2$estimate-1.96*prop_sigma
      interval<-data.frame(value=as.numeric(tmp2$estimate), SE=as.numeric(prop_sigma), lower=as.numeric(lower), upper=as.numeric(upper))
      D_Pred<-tmp2$estimate[3]+tmp2$estimate[4]*z
      
    } else {
      interval<-data.frame(value=as.numeric(tmp2$estimate), SE=rep(NA, 5), lower=rep(NA, 5), upper=rep(NA, 5))  
      D_Pred<-tmp2$estimate[3]+tmp2$estimate[4]*z
      
    } 
    
    JointA0[s,]<-as.numeric(interval[1,])
    JointSigma[s,]<-as.numeric(interval[2,])
    JointB0[s,]<-as.numeric(interval[3,])
    JointB1[s,]<-as.numeric(interval[4,])
    JointN[s,1]<-sum(exp(D_Pred))
    JointChat[s]<-as.numeric(exp(tmp2$estimate[5]))
    
  }
  
  out=list(
    SCR_nind=SCR_nind,
    SCR_ndet=SCR_ndet,
    Occ_ndet=Occ_ndet,
    JointA0=JointA0,
    JointSigma=JointSigma,
    JointB0=JointB0,
    JointB1=JointB1,
    JointN=JointN,
    varCovJoint=varCovJoint,
    JointChat=JointChat,
    z=z)
  return(out)
}


X2=cbind( sort(rep( seq(10,90, 10),9)), rep( seq(10,90,10),9))
blah2<-exp(-(e2dist(X2, X2)^2)/(2*2^2))
blah3<-blah2
diag(blah3)<-0
weightdummy<-1-(rowSums(blah3)/rowSums(blah2))


sim1W<-Sim_fit_SCR_OCC_joint2modweights(nsims=300, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                                        X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                                        X2=cbind( sort(rep( seq(10,90, 10),9)), rep( seq(10,90,10),9)), trapweights=weightdummy)


X2=cbind( sort(rep( seq(63,87, 3),9)), rep( seq(28,52,3),9))
blah2<-exp(-(e2dist(X2, X2)^2)/(2*2^2))
blah3<-blah2
diag(blah3)<-0
weightdummy<-1-(rowSums(blah3)/rowSums(blah2))


sim2w<-Sim_fit_SCR_OCC_joint2modweights(nsims=300, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                                        X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                                        X2=cbind( sort(rep( seq(63,87, 3),9)), rep( seq(28,52,3),9)), trapweights=weightdummy)


X2=cbind( sort(rep(c(seq(10, 16, 3), seq(45, 51, 3), seq(80, 86, 3)), 9)), 
          rep(c(seq(10, 16, 3), seq(45, 51, 3), seq(80, 86, 3)), 9))
sim3w<-Sim_fit_SCR_OCC_joint2modweights(nsims=300, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                                        X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                                        X2=cbind( sort(rep(c(seq(10, 16, 3), seq(45, 51, 3), seq(80, 86, 3)), 9)), 
                                                  rep(c(seq(10, 16, 3), seq(45, 51, 3), seq(80, 86, 3)), 9)), sigma.weights = as.numeric(exp(sim3$JointSigma[, 1])))


X2=cbind( sort(rep( seq(10,90, 4),21)), rep( seq(10,90,4), 21))
blah2<-exp(-(e2dist(X2, X2)^2)/(2*2^2))
blah3<-blah2
diag(blah3)<-0
weightdummy<-1-(rowSums(blah3)/rowSums(blah2))

sim4w<-Sim_fit_SCR_OCC_joint2modweights(nsims=300, K=20, K2=20, alpha0 = -2, sigma = 2, N= 400, beta1=1, 
                                        X1=cbind( sort(rep( seq(40,56,2),9)), rep( seq(20,36,2),9)),
                                        X2=cbind( sort(rep( seq(10,90, 4),21)), rep( seq(10,90,4), 21)), trapweights=weightdummy)




###Okay...consider bias, coverage.
### bias is relatively straightforward...

rb1<-mean((sim1$JointN[i, 1]-400)/400)
rb2<-mean((sim2$JointN[i, 1]-400)/400)
rb3<-mean((sim3$JointN[i, 1]-400)/400)
rb4<-mean((sim4$JointN[i, 1]-400)/400)


rb2w<-mean((sim2w$JointN[i, 1]-400)/400)
rb3w<-mean((sim3w$JointN[i, 1]-400)/400)
rb4w<-mean((sim4w$JointN[i, 1]-400)/400)



###Coverage requires a few steps.
###probably more convenient to embed these calculations within the simulation
###function, but so it goes.


VarNSim1<-rep(NA, 300)
VarNSim2<-rep(NA, 300)
VarNSim3<-rep(NA, 300)
VarNSim4<-rep(NA, 300)
VarNSim2w<-rep(NA, 300)
VarNSim3w<-rep(NA, 300)
VarNSim4w<-rep(NA, 300)


CprimeSim1<-rep(NA, 300)
CprimeSim2<-rep(NA, 300)
CprimeSim3<-rep(NA, 300)
CprimeSim4<-rep(NA, 300)
CprimeSim2w<-rep(NA, 300)
CprimeSim3w<-rep(NA, 300)
CprimeSim4w<-rep(NA, 300)


for (i in 1:300){
  VarNSim1[i]<-deltavar(sum(exp(b0+b1*z)), 
                        meanval=c(b0=sim1$JointB0[i, 1], b1=sim1$JointB1[i, 1]),
                        Sigma=sim1$varCovJoint[,,i])
  
  VarNSim2[i]<-deltavar(sum(exp(b0+b1*z)), 
                        meanval=c(b0=sim2$JointB0[i, 1], b1=sim2$JointB1[i, 1]),
                        Sigma=sim2$varCovJoint[,,i])
  
  VarNSim3[i]<-deltavar(sum(exp(b0+b1*z)), 
                        meanval=c(b0=sim3$JointB0[i, 1], b1=sim3$JointB1[i, 1]),
                        Sigma=sim3$varCovJoint[,,i])
  
  VarNSim4[i]<-deltavar(sum(exp(b0+b1*z)), 
                        meanval=c(b0=sim4$JointB0[i, 1], b1=sim4$JointB1[i, 1]),
                        Sigma=sim4$varCovJoint[,,i])
  
  VarNSim2w[i]<-deltavar(sum(exp(b0+b1*z)), 
                         meanval=c(b0=sim2w$JointB0[i, 1], b1=sim2w$JointB1[i, 1]),
                         Sigma=sim2w$varCovJoint[,,i])
  
  
  VarNSim3w[i]<-deltavar(sum(exp(b0+b1*z)), 
                         meanval=c(b0=sim3w$JointB0[i, 1], b1=sim3w$JointB1[i, 1]),
                         Sigma=sim3w$varCovJoint[,,i])
  
  VarNSim4w[i]<-deltavar(sum(exp(b0+b1*z)), 
                         meanval=c(b0=sim4w$JointB0[i, 1], b1=sim4w$JointB1[i, 1]),
                         Sigma=sim4w$varCovJoint[,,i])
  
  CprimeSim1[i]<-exp(1.96*sqrt(log(1+VarNSim1[i]/sim1$JointN[i, 1]^2)))  
  CprimeSim2[i]<-exp(1.96*sqrt(log(1+VarNSim2[i]/sim2$JointN[i, 1]^2)))  
  CprimeSim3[i]<-exp(1.96*sqrt(log(1+VarNSim3[i]/sim3$JointN[i, 1]^2)))  
  CprimeSim4[i]<-exp(1.96*sqrt(log(1+VarNSim4[i]/sim4$JointN[i, 1]^2)))  
  CprimeSim2w[i]<-exp(1.96*sqrt(log(1+VarNSim2w[i]/sim2w$JointN[i, 1]^2)))  
  CprimeSim3w[i]<-exp(1.96*sqrt(log(1+VarNSim3w[i]/sim3w$JointN[i, 1]^2)))  
  CprimeSim4w[i]<-exp(1.96*sqrt(log(1+VarNSim4w[i]/sim4w$JointN[i, 1]^2)))  
  
  
  sim1$JointN[i, 2]<-sim1$JointN[i, 1]/CprimeSim1[i] ###lower ci
  sim1$JointN[i, 3]<-sim1$JointN[i, 1]*CprimeSim1[i] ### upper ci
  
  sim2$JointN[i, 2]<-sim2$JointN[i, 1]/CprimeSim2[i]
  sim2$JointN[i, 3]<-sim2$JointN[i, 1]*CprimeSim2[i]
  
  sim3$JointN[i, 2]<-sim3$JointN[i, 1]/CprimeSim3[i]
  sim3$JointN[i, 3]<-sim3$JointN[i, 1]*CprimeSim3[i]
  
  sim4$JointN[i, 2]<-sim4$JointN[i, 1]/CprimeSim4[i]
  sim4$JointN[i, 3]<-sim4$JointN[i, 1]*CprimeSim4[i]
  
  sim2w$JointN[i, 2]<-sim2w$JointN[i, 1]/CprimeSim2w[i]
  sim2w$JointN[i, 3]<-sim2w$JointN[i, 1]*CprimeSim2w[i]
  
  sim3w$JointN[i, 2]<-sim3w$JointN[i, 1]/CprimeSim3w[i]
  sim3w$JointN[i, 3]<-sim3w$JointN[i, 1]*CprimeSim3w[i]
  
  sim4w$JointN[i, 2]<-sim4w$JointN[i, 1]/CprimeSim4w[i]
  sim4w$JointN[i, 3]<-sim4w$JointN[i, 1]*CprimeSim4w[i]
  
}





###Empirical case study here...
###not run--the creation of certain variables outlined below...
#blah2<-exp(-(e2dist(TraplocsOcc, TraplocsOcc)^2)/(2*1.53^2))
#blah3<-blah2
#diag(blah3)<-0
#weightbcat<-1-(rowSums(blah3)/rowSums(blah2)) ###creating trap specific weights




###Full_Mask is a much larger set of coordinates than neccessary
###define msk1 and msk2...


#Msk1<-Full_Mask[Idx_SCR$x, 1:2]/1000 ###set of integration points for scrdata
#Msk2<-Full_Mask[Idx_2018$x, 1:2]/1000 ### set of integration points for occupancy data


#D <- e2dist(TraplocsSCR, Msk1) creating 
#D_2<-e2dist(TraplocsOcc, Msk2_1)







#Woody covariate based on scaling the sum of wooded cover types--e.g....
#...All_Woody_Full_Mask<-as.numeric(scale(Full_Mask$PropDF+Full_Mask$PropEF+Full_Mask$PropMF+Full_Mask$PropWW+Full_Mask$PropSH))

#WoodyNieghs <-matrix(All_Woody_Full_Mask[as.numeric(N_3mat)], nrow(N_3mat), ncol(N_3mat))    

###create spline bases functions
#library(fields)
#knots2<-cover.design(cbind(scale(Full_Mask$X), scale(Full_Mask$Y)), 10)

# Define the omega and Z matrices for the random effects (see Crainiceanu et al. 2005)
#omega.all <- (rdist(knots2$design, knots2$design))^3
#svd.omega.all <- svd(omega.all)
#sqrt.omega.all <- t(svd.omega.all$v %*% (t(svd.omega.all$u)*sqrt(svd.omega.all$d)))
#Z.k <- (rdist(cbind(scale(Full_Mask$X), scale(Full_Mask$Y)), knots2$design))^3
#Z.matrix <- t(solve(sqrt.omega.all, t(Z.k)))
#using scaled x/y to create basis functions keeps basis function values on similar
#scale to covariates, etc.



###read data
LyRuCaps<-read.csv('LyRuCaps.csv') ###SCR data, 2-d form
LyRu_Traps<-read.csv("LyRu_Traps.csv") ###SCR trap data
LyruOcc<-read.csv('LyruOcc.csv')
D<-read.csv('DMat.csv')
#D_2 matrix too large to upload here...see line 616 for how created...will
#have to get this onto a different archiving site to make recreatable 
#D_2<-read.csv('D_2mat.csv')
All_Woody_Full_Mask<-read.csv('All_Woody_Full_Mask.csv')
Z.matrix<-as.matrix(read.csv('Z.matrix.csv'))
Idx_SCR<-read.csv('Idx_SCR.csv') ###which subset of full mask to use for calcs in scr likeliood
weightbcat<-read.csv('weightbcat.csv')
Idx_WI<-read.csv("Idx_WI.csv") ###for prediction/region.N type calculations,
### need to get rid of mask areas outside of WI, some of which are needed for integration
water<-read.csv('water.csv') ### some parts of the region of prediction are in open water...need to get rid of this, too




##D_3 matrix is wayyyy to large to store--need to recrete from full mask... 
##

Full_Mask<-read.csv('Full_Mask.csv')

D_3<-dnearneigh(as.matrix(Full_Mask[, 1:2]/1000), 0.5, 7)
N_3<-D_3
D_3<-list()
for(i in 1:length(N_3)){
  D_3[[i]]<-unlist(lapply(N_3[[i]], function(x) e2dist(Full_Mask[x, 1:2]/1000, Full_Mask[i, 1:2]/1000)))
}
###N_3<-refers to neighbors
###D_3<-refers to distances
rowMax <- max(sapply(N_3, length)) 
N_3mat<-do.call(rbind, lapply(N_3, function(x){ 
  length(x) <- rowMax 
  x
}))
D_3mat<-do.call(rbind, lapply(D_3, function(x){ 
  length(x) <- rowMax 
  x
}))

WoodyNieghs <-matrix(All_Woody_Full_Mask[as.numeric(N_3mat)], nrow(N_3mat), ncol(N_3mat))       


IntLik_Spline_NB_Woody<-function (start = rnorm(19, 0, 1), y = LyRuCaps, K = LyRu_Traps$Activedays,  
                                  y2=LyruOcc$Ndets, K2=LyruOcc$Ndays,
                                  D=D, D_2=D_2, z1 =All_Woody_Full_Mask , z2=WoodyNieghs,
                                  tcovSCR=LyRu_Traps$Bare_trail, tcovOcc=LyruOcc$TRAVEL_CORRIDOR,
                                  D_3=D_3mat, thinSCR=Idx_SCR$x, 
                                  Zbas1=Z.matrix[,1], Zbas2=Z.matrix[,2], Zbas3=Z.matrix[,3],
                                  Zbas4=Z.matrix[,4], Zbas5=Z.matrix[,5], Zbas6=Z.matrix[,6], Zbas7=Z.matrix[,7],
                                  Zbas8=Z.matrix[,8], Zbas9=Z.matrix[,9], Zbas10=Z.matrix[,10], weights=weightbcat)
{
  require(scrbook)
  nG1 <- nrow(G1)
  alpha0SCR <- start[1]
  alpha0occ<-start[2]
  alpha1<-start[3]
  sigma <- exp(start[4])
  sigmaN<-exp(start[5])
  twosigsq<-2*sigmaN^2
  z3<-rowSums(z2*exp(-D_3^2/twosigsq), na.rm=T)/max(rowSums(exp(-D_3^2/twosigsq), na.rm=T))
  Dens<-exp(start[6]+start[7]*z1+start[8]*z3+start[9]*Zbas1+start[10]*Zbas2+start[11]*Zbas3+start[12]*Zbas4+
              start[13]*Zbas5+start[14]*Zbas6+start[15]*Zbas7+start[16]*Zbas8+start[17]*Zbas9+start[18]*Zbas10)
  THETA<-exp(start[19])
  Dens1<-Dens[thinSCR]
  psi<-Dens1/sum(Dens1)
  loglam <- alpha0SCR+alpha1*tcovSCR - (1/(2 * sigma * sigma)) * D * D # + alpha2 * ztrap Ignore RSF
  probcap <- 1 - exp(-exp(loglam))
  Pm <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
  ymat <- y
  ymat <- rbind(y, rep(0, ncol(y)))
  lik.marg <- rep(NA, nrow(ymat))
  
  loglamOcc<-alpha0occ + alpha1*tcovOcc - (1/(2 * sigma * sigma)) * D_2 * D_2
  Dens2<-Dens
  lamOcc<-exp(loglamOcc)
  lamOcc2<-matrix(NA, nrow(lamOcc), ncol(lamOcc))
  Pm3<-rep(NA, nrow(lamOcc))
  like.occ<-rep(NA,nrow(lamOcc))
  for (i in 1:nrow(lamOcc)){
    z<-which(D_2[i, ] <= 7)
    lamOcc2[i, z]<-lamOcc[i, z]*Dens2[z]
    Pm3[i]<-1-exp(-sum(lamOcc2[i, z]))
    like.occ[i]<-dbetabinom(y2[i],Pm3[i], K2[i], theta=THETA, log=TRUE)*weights[i]
  }
  
  for (i in 1:nrow(ymat)) {
    Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG1), rep(K[i],nG1), probcap[1:length(Pm)], log = TRUE))
    lik.cond <- exp(colSums(Pm))
    lik.marg[i] <- sum(lik.cond * psi)
  }
  
  nv<-c(rep(1, length(lik.marg)-1), 1)
  atheta<-1-lik.marg[nrow(ymat)]
  nind<-nrow(ymat)-1
  pixels<-rep(1, nG1)
  part1 <- nind * log(sum(Dens1*pixels)) - sum(Dens1 * pixels)*atheta
  part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
  part3<-sum(like.occ)
  out <- -1 * (part1 + part2+part3) 
  
}   

start<-c(-4, -4.5, .8, .4, .7, -4, -.1, .4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3)
fit<-nlm(IntLik_Spline_NB_Woody, p=start, hessian = TRUE, print.level = 2, iterlim=500)


fisher_info<-solve(fit$hessian)
SE<-sqrt(diag(fisher_info))



dummypred<-exp(fit$estimate[6]+fit$estimate[7]*All_Woody_Full_Mask[Idx_WI$x]+
              fit$estimate[8]*rowSums(WoodyNieghs[Idx_WI$x,]*exp(-D_3mat[Idx_WI$x,]^2/
              (2*exp(fit$estimate[5])^2)), na.rm=T)/max(rowSums(exp(-D_3mat[Idx_WI$x,]^2/(2*exp(fit$estimate[5])^2)), na.rm=T))+
              Z.matrix[Idx_WI$x,]%*%fit_4_2_2020$estimate[9:18])

Nhat<-sum(dummypred[-water])

dummyvar<-deltavar(sum(exp(b0+b1*All_Woody_Full_Mask[Idx_WI$x[-water]]+
                      b2*rowSums(WoodyNieghs[Idx_WI$x[-water],]*exp(-D_3mat[Idx_WI$x[-water],]^2/(2*exp(a0)^2)), na.rm=T)/max(rowSums(exp(-D_3mat[Idx_WI$x[-water],]^2/(2*exp(a0)^2)), na.rm=T))+
                      b3*Z.matrix[Idx_WI$x[-water],1]+b4*Z.matrix[Idx_WI$x[-water],2]+b5*Z.matrix[Idx_WI$x[-water],3]+b6*Z.matrix[Idx_WI$x[-water],4]+b7*Z.matrix[Idx_WI$x[-water],5]+b8*Z.matrix[Idx_WI$x[-water],6]+
                      b9*Z.matrix[Idx_WI$x[-water],7]+b10*Z.matrix[Idx_WI$x[-water],8]+b11*Z.matrix[Idx_WI$x[-water],9]+b12*Z.matrix[Idx_WI$x[-water],10])), 
                  meanval=c(a0=fit$estimate[5], b0=fit$estimate[6], b1=fit$estimate[7], b2=fit$estimate[8],
                      b3=fit$estimate[9], b4=fit$estimate[10], b5=fit$estimate[11],
                      b6=fit$estimate[12], b7=fit$estimate[13], b8=fit$estimate[14],
                      b9=fit$estimate[15], b10=fit$estimate[16], b11=fit$estimate[17],
                      b12=fit$estimate[18]),
                      Sigma=fisher_info[5:18, 5:18])


Cprime<-exp(1.96*sqrt(log(1+dummyvar/Nhat)))
print(c(Nhat, Nhat/Cprime, Nhat*Cprime))


###report pixel-wise uncertainty
varpred<-deltavar(exp(b0+b1*All_Woody_Full_Mask[Idx_WI$x[-water]]+
                                 b2*rowSums(WoodyNieghs[Idx_WI$x[-water],]*exp(-D_3mat[Idx_WI$x[-water],]^2/(2*exp(a0)^2)), na.rm=T)/max(rowSums(exp(-D_3mat[Idx_WI$x[-water],]^2/(2*exp(a0)^2)), na.rm=T))+
                                 b3*Z.matrix[Idx_WI$x[-water],1]+b4*Z.matrix[Idx_WI$x[-water],2]+b5*Z.matrix[Idx_WI$x[-water],3]+b6*Z.matrix[Idx_WI$x[-water],4]+b7*Z.matrix[Idx_WI$x[-water],5]+b8*Z.matrix[Idx_WI$x[-water],6]+
                                 b9*Z.matrix[Idx_WI$x[-water],7]+b10*Z.matrix[Idx_WI$x[-water],8]+b11*Z.matrix[Idx_WI$x[-water],9]+b12*Z.matrix[Idx_WI$x[-water],10]), 
                  meanval=c(a0=fit$estimate[5], b0=fit$estimate[6], b1=fit$estimate[7], b2=fit$estimate[8],
                            b3=fit$estimate[9], b4=fit$estimate[10], b5=fit$estimate[11],
                            b6=fit$estimate[12], b7=fit$estimate[13], b8=fit$estimate[14],
                            b9=fit$estimate[15], b10=fit$estimate[16], b11=fit$estimate[17],
                            b12=fit$estimate[18]),
                  Sigma=fisher_info[5:18, 5:18])

SE_Pred<-sqrt(varpred)