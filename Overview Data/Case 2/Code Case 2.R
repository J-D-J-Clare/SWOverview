###import data into workspace





Trail<-read.csv('Trail.csv')
Forest5km<-read.csv('Forest5km.csv') ###scaled. mean = 0.521, sd = 0.248
Forest<-read.csv('Forest.csv')    ###scaled. mean = 0.609, sd = 0.316
Cropland<-read.csv('Cropland.csv') ###scaled. mean = 0.089, sd = 0.193
EVIDAY<-read.csv('EVIDAY.csv') ###scaled. mean = 0.340, sd = 0.160
BasXY<-read.csv('BasXY.csv')
knownocc<-read.csv('knownocc.csv')
knownavailable<-read.csv('knownavailable.csv')
camact<-read.csv('camact.csv')
basis_t<-read.csv('basis_t.csv')

DetectData<-read.csv('DetectData.csv')

library(dplyr)
tmp<-DetectData %>% group_by(PolyIDX_500m, Day) %>% mutate(ID=row_number())


y.arr<-array(NA, dim=c(1057, 365, 2))
T_idx<-array(NA, dim=c(1057, 365, 2))

for (i in 1:nrow(tmp)){
  y.arr[tmp$PolyIDX_500m[i], tmp$Day[i], tmp$ID[i]]<-tmp$BearPA[i]
  T_idx[tmp$PolyIDX_500m[i], tmp$Day[i], tmp$ID[i]]<-tmp$CamNUM[i]
}



y.arr[is.na(y.arr)]<-2 ###stan does not permit NA's--need to fill the arrays with something
T_idx[is.na(T_idx)]<-5000

###Compile data object for stan



data <- list(y = y.arr+1, # detection history
             T_idx=T_idx,
             ncells = 1057,
             ncams= 1103,
             Trail=as.numeric(unlist(Trail)),
             Forest5km=as.numeric(unlist(Forest5km)),
             Forest=as.numeric(unlist(Forest)),
             Cropland=as.numeric(unlist(Cropland)),
             EVIDAY=as.matrix(EVIDAY), 
             num_basisXY=30,
             num_bastime = 15,
             BasXY=as.matrix(BasXY),
             knownocc=as.numeric(unlist(knownocc)),
             knownavailable=as.matrix(knownavailable), 
             camact=as.matrix(camact),
             basis_t=as.matrix(basis_t))


###stan model

write(" // 
      
      data {
      int<lower=1> ncells;                // number of cells sampled--distinct from cell number across broader lattice
      int<lower=1> ncams;                 //number of distinct cameras
      int y[ncells, 365, 2];                 // observation vector
      int camact[ncells, 365];             //how many cameras active in given cell during particular day (ranges from 0 - 2)
      int T_idx[ncells, 365, 2];          //which camera correponds to a given y?
      real Forest5km[ncells];
      int num_basisXY;
      int num_bastime;
      matrix[num_basisXY, ncells] BasXY;
      matrix[num_bastime, 365] basis_t;
      real Forest[ncells];       //% forest locally (500 m pixel)
      real Cropland[ncells];    //% cropland locally (500 m pixel)
      real EVIDAY[ncells, 365]; //the EVI on the day of interest
      real Trail[ncams];         //camera on maintained trail or not
      int<lower=0, upper=1> knownocc[ncells]; //bear observed in cell ever
      int<lower=0, upper=2> knownavailable[ncells, 365]; //bear observed somewhere in cell on particular day?
      }
      
      
      parameters {
      real a0_psi; 
      real a1_psi; //the effect of forest
      row_vector[num_basisXY] a_XY;
      
      
      real a_EVI; //effect of evi
      real a_AL;  //effect of availability in previous time period
      row_vector[num_bastime] a_day_raw;
      row_vector[num_bastime] a_day_forest_raw;
      row_vector[num_bastime] a_day_cropland_raw;
      
      real<lower=0> var_a_day; //intercept varies by day
      real<lower=0> var_a_day_forest; //effect of forest varies by day
      real<lower=0> var_a_day_cropland; //effect of cropland varies by day
      
      
      
      real mu_b0_p; //this is the mean for the camera-specific random intercept
      real b1_p;
      real<lower=0> var_p;
      vector[ncams] b0_p_raw;
      
      
      
      
      }
      
      transformed parameters {
      vector<lower=0, upper=1>[ncams] p;   //prob of true detection per unit given covs
      vector<lower=0, upper=1>[ncells] psi;    //prob of occupancy 
      matrix[ncells, 365] logittheta;              //prob available
      
      
      simplex[3] ps[3, ncells, 364];  //transitions between states
      simplex[2] po[3, ncams]; //emissions--or p(y)|state
      
      
      //states: 1=unoccupied, 2= occupied, not available, 3= occupied and available
      //observations: 1=not seen, 2= seen;
      
      row_vector[num_bastime] a_day;
      row_vector[num_bastime] a_day_forest;
      row_vector[num_bastime] a_day_cropland;
      
      
      vector[ncams] b0_p;
      
      b0_p = mu_b0_p + var_p * b0_p_raw;
      
      a_day[1]=a_day_raw[1];
      a_day_forest[1]=a_day_forest_raw[1];
      a_day_cropland[1]=a_day_cropland_raw[1];
      
      for (t in 2:num_bastime){
        a_day[t]=a_day[t-1]+a_day_raw[t]*var_a_day;
        a_day_forest[t]=a_day_forest[t-1]+a_day_forest_raw[t]*var_a_day_forest;
        a_day_cropland[t]=a_day_cropland[t-1]+a_day_cropland_raw[t]*var_a_day_cropland;
      }
      
      for (c in 1:ncams){
        p[c] = inv_logit(b0_p[c]+b1_p*Trail[c]);       
        po[1, c, 1]=1;     
        po[1, c, 2]=0;
        po[2, c, 1]=1;
        po[2, c, 2]=0;
        po[3, c, 1]=1-p[c];
        po[3, c, 2]=p[c];
      }
      
      for (i in 1:ncells){
        psi[i]=inv_logit(a0_psi+a1_psi*Forest5km[i]+dot_product(a_XY, BasXY[,i]));
          for (t in 1:365){
            logittheta[i, t]=dot_product(a_day, basis_t[,t])+dot_product(a_day_forest, basis_t[,t])*Forest[i]+
            dot_product(a_day_cropland, basis_t[,t])*Cropland[i]+a_EVI*EVIDAY[i, t];       
          }
      }
      
      for (i in 1:ncells){
        for (t in 1:364){
          ps[1, i, t, 1]=1; //unoccupied remains unoccupied
          ps[1, i, t, 2]=0;
          ps[1, i, t, 3]=0;
          ps[2, i, t, 1]=0;
          ps[2, i, t, 2]=1-inv_logit(logittheta[i, t+1]); 
          ps[2, i, t, 3]=inv_logit(logittheta[i, t+1]);
          ps[3, i, t, 1]=0;
          ps[3, i, t, 2]=1-inv_logit(logittheta[i, t+1]+a_AL);
          ps[3, i, t, 3]=inv_logit(logittheta[i, t+1]+a_AL);
        }
      }
      
      }//end transformed
      
      model {
      real acc[3];
      vector[3] gam[365];
      // priors
      
      
      for (b in 1:num_basisXY){ //hiearchical approach here works poorly... 
        a_XY[b]~normal(0, 1.6);
      } 
      
      
      
      
      var_p ~ cauchy (0, 1); //random variance for obs. component
      mu_b0_p ~ normal (0, 1.6); 
      b0_p_raw ~ std_normal();
      b1_p ~ normal (0, 1.6);
      a0_psi ~ normal (0, 1.6); 
      a1_psi ~ normal (0, 1.6); 
      
      a_EVI ~ normal (0, 1.6); 
      a_AL ~normal (0, 1.6);
      
      
      a_day_raw~std_normal();
      a_day_forest_raw~std_normal();
      a_day_cropland_raw~std_normal();
      
      
      var_a_day ~ cauchy (0, 1);
      var_a_day_forest ~ cauchy (0, 1);
      var_a_day_cropland ~ cauchy (0, 1);
      
      
      //likelihood...
      
      
      for (i in 1:ncells) {
      // initial state
        if (camact[i, 1]==0){
          gam[1, 1] = 1-psi[i];
          gam[1, 2] = psi[i]*(1-inv_logit(logittheta[i, 1]));
          gam[1, 3] = psi[i]*inv_logit(logittheta[i, 1]);
      }     
      
      else{
        gam[1, 1] =1-psi[i]* po[1, T_idx[i, 1, 1], y[i, 1, 1]];
        gam[1, 2] =psi[i]*(1-inv_logit(logittheta[i, 1]))* po[2, T_idx[i, 1, 1], y[i, 1, 1]];
        gam[1, 3] =psi[i]*inv_logit(logittheta[i, 1])* po[3, T_idx[i, 1, 1], y[i, 1, 1]];
      }
      
        for (t in 2:365) {
          for (k in 1:3) {
            for (j in 1:3){
              if (camact[i, t]==0){
                acc[j] = gam[t - 1, j] * ps[j, i, t - 1, k];
              }
              else if (camact[i, t]==1){
                acc[j] = gam[t - 1, j] *ps[j, i, t - 1, k]*po[k, T_idx[i, t, 1], y[i, t, 1]];
              }
              else {
                acc[j] = gam[t - 1, j] *ps[j, i, t - 1, k]*po[k, T_idx[i, t, 1], y[i, t, 1]]*po[k, T_idx[i, t, 2], y[i, t, 2]];
              }
            gam[t, k] = sum(acc);
            }
          }
        }
        target += log(sum(gam[365]));
      
      } //end i loop
      
      } //end of model
      
      
      
      
      generated quantities {
      
      }
      
      
      ",
      
      "Bear_Intraannual.stan")





params <- c('a0_psi', 'a1_psi', "a_XY",  "var_p", "b1_p", "mu_b0_p",
            "a_day","a_day_forest", "a_day_cropland", "a_EVI", "a_AL",
            "var_a_day", "var_a_day_forest", "var_a_day_cropland")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

outstanHMM2 = stan("Bear_Intraannual.stan",
                   data = data,
                   pars = params,
                   chains = 2,
                   #init = inits,
                   iter = 2000,
                   warmup = 1000, 
                   thin = 1,
                   seed = 1,
                   open_progress = FALSE,
                   control=list(adapt_delta=0.8))
