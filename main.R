#Estimate PACE-PEG PK with PBPK modeling
#Owen Richfield
#1/9/2024
#Collaboration with Laura Bracaglia, Alex Piotrowski-Daspit

#PBPK Model structure: compartmental with compartments for:

#Brain
#Liver
#Spleen
#Heart
#Bone
#Kidneys
#Lungs

#Organism: Mouse (0.02 kg)

#Code is organized as follows:
  #First section: verify conservation of dose
  #Second section: rescale parameters using correlations
  #Third section: second round of rescaling parameters
  #Fourth section: fine tuning parameters

  #Fifth section: simulation runs with model
  
#############################################################
#First section: verify conservation of dose

#Clear all variables
rm(list=ls())

#Source scripts
source("libs.R")
source("parms.R")
source("funcs.R")
#source("FACS_IVIS_data_analysis_20230426.R")

#Set output parameters to 0. 

Kbile <- 0
Kurine <- 0
nt <- nt_Au

Kmax_t_Liver <- Kmax_t_Liver_Au
Kmax_t_Spleen <- Kmax_t_Spleen_Au
Kmax_t_Kidneys <- Kmax_t_Kidneys_Au
Kmax_t_Lungs <- Kmax_t_Lungs_Au

Krelease_t_Liver <- Krelease_t_Liver_Au
Krelease_t_Spleen <- Krelease_t_Spleen_Au
Krelease_t_Kidneys <- Krelease_t_Kidneys_Au
Krelease_t_Lungs <- Krelease_t_Lungs_Au

K50_t_Liver <- K50_t_Liver_Au
K50_t_Spleen <- K50_t_Spleen_Au
K50_t_Kidneys <- K50_t_Kidneys_Au
K50_t_Lungs <- K50_t_Lungs_Au

PACt_Liver <- PACt_Liver_Au
PACt_Lungs <- PACt_Lungs_Au
PACt_Kidneys <- PACt_Kidneys_Au
PACt_Spleen <- PACt_Spleen_Au
PACt_Bone <- PACt_Bone_Au
PACt_Heart <- PACt_Heart_Au
PACt_Rest <- PACt_Rest_Au
PACt_Brain <- PACt_Brain_Au

Init_condition <- rep(0,22)

OP <- ode(y=Init_condition,
          times=times,
          func = NP_PBPK2, 
          parms=list(),
          method='lsoda')
OP <- data.frame(OP)
OP_rename <- OP

print("Conservation of dose?")
#max(rowSums(OP[,2:dim(OP)[2]]))-min(rowSums(OP[,2:dim(OP)[2]])) < 1e-13

indx_samp <- seq(from = 1, to = Nt, by = 1)

plot(times[times<0.1],rowSums(OP_rename[,2:(dim(OP_rename)[2])])[times<0.1])

max(abs(rowSums(OP_rename[times > 1,2:(dim(OP_rename)[2])])-dose)/dose)*100

########################################################
#Second section: rescale parameters using correlations

#Clear all variables
rm(list=ls())

set.seed(7)

#Source scripts
source("libs.R")
source("parms.R")
alpha_temp_FACS <- 1e-2
alpha_temp_blood <- 1e-2
source("funcs.R")

Nit <- 1e4
PACt_mod <- 1#e-3
Kmax_Spleen_mod <- 1#e-2

for (j in 1:length(dose_rel)){
  
source("FACS_IVIS_data_analysis_20230426.R")
  
blood_conc <- blood_conc[blood_conc$dose==dose_rel[j],]
blood_summary <- blood_summary[blood_summary$dose==dose_rel[j],]
organ_FACS_summary <- organ_FACS_summary[organ_FACS_summary$Dose==dose_rel[j],]
organ_FACS_time48_summary <- organ_FACS_summary
organ_FACS_time3_summary <- organ_FACS_time_summary[organ_FACS_time_summary$time==3,]

dose_rate <- dose_rate_vec[j]
  
Init_condition <- rep(0,22)

OP_IS <- data.frame(
  dose = rep(dose_rel[j],Nit),
  
  nt = rep(0,Nit),
  Kbile = rep(0,Nit),
  Kurine = rep(0,Nit),
  
  PACt_Liver = rep(0,Nit),
  PACt_Kidneys = rep(0,Nit),
  PACt_Spleen = rep(0,Nit),
  PACt_Lungs = rep(0,Nit),
  PACt_Heart = rep(0,Nit),
  PACt_Bone = rep(0,Nit),
  PACt_Rest = rep(0,Nit),
  
  Kmax_t_Liver = rep(0,Nit),
  Kmax_t_Spleen  = rep(0,Nit),
  Kmax_t_Kidneys  = rep(0,Nit),
  Kmax_t_Lungs  = rep(0,Nit),
  
  Krelease_t_Liver  = rep(0,Nit),
  Krelease_t_Spleen  = rep(0,Nit),
  Krelease_t_Kidneys  = rep(0,Nit),
  Krelease_t_Lungs = rep(0,Nit),
  
  K50_t_Liver  = rep(0,Nit),
  K50_t_Spleen  = rep(0,Nit),
  K50_t_Kidneys  = rep(0,Nit),
  K50_t_Lungs = rep(0,Nit),
  
  blood_03 = rep(0,Nit),
  
  blood_2 = rep(0,Nit),
  
  blood_4 = rep(0,Nit),
  
  blood_8 = rep(0,Nit),
  
  blood_10 = rep(0,Nit),
  
  blood_24 = rep(0,Nit),
  
  blood_48 = rep(0,Nit),
  
  Heart_3 = rep(0,Nit),
  
  Kidneys_3 = rep(0,Nit),
  
  Spleen_3 = rep(0,Nit),
  
  Lungs_3 = rep(0,Nit),
  
  Bone_3 = rep(0,Nit),
  
  Heart_48 = rep(0,Nit),
  
  Kidneys_48 = rep(0,Nit),
  
  Spleen_48 = rep(0,Nit),
  
  Lungs_48 = rep(0,Nit),
  
  Bone_48 = rep(0,Nit)
)

 for (i in seq(Nit)){
   
   nt <- nt_Au*10^round(runif(1,0,1))
   Kbile <- Kbile_Au*10^round(runif(1,-2,0))
   Kurine <- 0
   PACt_Kidneys <- 0
   
   PACt_Lungs <- PACt_Lungs_Au*10^round(runif(1,-5,1))*PACt_mod
   PACt_Heart <- PACt_Heart_Au*10^round(runif(1,-5,1))*PACt_mod
   PACt_Bone <- PACt_Bone_Au*10^round(runif(1,-5,1))*PACt_mod
   PACt_Spleen <- PACt_Spleen_Au*10^round(runif(1,-5,1))*PACt_mod
   PACt_Rest <- PACt_Rest_Au*10^round(runif(1,-5,1))*PACt_mod
   PACt_Brain <- PACt_Rest
   
   Kmax_t_Liver <- 10^round(runif(1,-4,1))
   Kmax_t_Kidneys <- 10^round(runif(1,-4,1))
   Kmax_t_Spleen <- 10^round(runif(1,-4,1))
   Kmax_t_Lungs <- 10^round(runif(1,-4,1))
   
   K50_t_Liver <- K50_t_Liver_Au*10^round(runif(1,-2,1))
   K50_t_Spleen <- K50_t_Spleen_Au*10^round(runif(1,-2,1))
   K50_t_Kidneys <- K50_t_Kidneys_Au*10^round(runif(1,-2,1))
   K50_t_Lungs <- K50_t_Lungs_Au*10^round(runif(1,-2,1))
   
   Krelease_t_Liver <- Kmax_t_Liver*10^round(runif(1,-3,3))
   Krelease_t_Spleen <- Kmax_t_Spleen*10^round(runif(1,-3,3))
   Krelease_t_Kidneys <- Kmax_t_Kidneys*10^round(runif(1,-3,3))
   Krelease_t_Lungs <- Kmax_t_Lungs*10^round(runif(1,-3,3))
   
   PACt_Liver <- Kbile*10^round(runif(1,-3,3))
   
   
OP <- ode(y=Init_condition,
          times=times,
          func = NP_PBPK2, 
          parms=list(),
          method='lsoda')
OP <- data.frame(OP)
OP_rename <- OP

indx_samp <- seq(from = 1, to = Nt, by = 1)

OP_rename$M_Art <- OP[indx_samp,2]
OP_rename$M_Ven <- OP[indx_samp,3]

OP_rename$M_vasc_Liver <- OP[indx_samp,4]
OP_rename$M_extra_Liver <- OP[indx_samp,5]
OP_rename$M_phago_Liver <- OP[indx_samp,6]

OP_rename$M_vasc_Spleen <- OP[indx_samp,7]
OP_rename$M_extra_Spleen <- OP[indx_samp,8]
OP_rename$M_phago_Spleen <- OP[indx_samp,9]

OP_rename$M_vasc_Kidneys <- OP[indx_samp,10]
OP_rename$M_extra_Kidneys <- OP[indx_samp,11]
OP_rename$M_phago_Kidneys <- OP[indx_samp,12]

OP_rename$M_vasc_Lungs <- OP[indx_samp,13]
OP_rename$M_extra_Lungs <- OP[indx_samp,14]
OP_rename$M_phago_Lungs <- OP[indx_samp,15]

OP_rename$M_vasc_Brain <- OP[indx_samp,16]
OP_rename$M_extra_Brain <- OP[indx_samp,17]

OP_rename$M_vasc_Heart <- OP[indx_samp,18]
OP_rename$M_extra_Heart <- OP[indx_samp,19]

OP_rename$M_vasc_Bone <- OP[indx_samp,20]
OP_rename$M_extra_Bone <- OP[indx_samp,21]

OP_rename$M_vasc_Rest <- OP[indx_samp,22]
OP_rename$M_extra_Rest <- OP[indx_samp,23]

LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/LV
LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/LV
SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/LV
HV <- OP_rename$M_extra_Heart/LV
BoneV <- OP_rename$M_extra_Bone/LV
LV <- 1

blood <- (OP_rename$M_Ven)/(BW*VVenC)

plot(times,blood,type='l',ylim=c(0, 1500))
points(blood_summary$time, blood_summary$LCI,col='blue')
points(blood_summary$time, blood_summary$UCI,col='red')

OP_IS$PACt_Liver[i] <- PACt_Liver
OP_IS$Kbile[i] = Kbile

OP_IS$PACt_Kidneys[i] <- PACt_Kidneys
OP_IS$Kurine[i] = Kurine

OP_IS$PACt_Bone[i] <- PACt_Bone
OP_IS$PACt_Spleen[i] <- PACt_Spleen
OP_IS$PACt_Lungs[i] <- PACt_Lungs
OP_IS$PACt_Heart[i] <- PACt_Heart
OP_IS$PACt_Rest[i] <- PACt_Rest

OP_IS$nt[i] = nt

OP_IS$Kmax_t_Liver[i] = Kmax_t_Liver 
OP_IS$Kmax_t_Spleen[i]  = Kmax_t_Spleen
OP_IS$Kmax_t_Kidneys[i]  = Kmax_t_Kidneys 
OP_IS$Kmax_t_Lungs[i]  = Kmax_t_Lungs

OP_IS$Krelease_t_Liver[i]  = Krelease_t_Liver
OP_IS$Krelease_t_Spleen[i]  = Krelease_t_Spleen
OP_IS$Krelease_t_Kidneys[i]  = Krelease_t_Kidneys
OP_IS$Krelease_t_Lungs[i] = Krelease_t_Lungs

OP_IS$K50_t_Liver[i]  = K50_t_Liver
OP_IS$K50_t_Spleen[i]  = K50_t_Spleen
OP_IS$K50_t_Kidneys[i]  = K50_t_Kidneys
OP_IS$K50_t_Lungs[i] = K50_t_Lungs

OP_IS$blood_03[i] = blood[which.min(abs(times-0.03))]

OP_IS$blood_1[i] = blood[which.min(abs(times-1))]

OP_IS$blood_2[i] = blood[which.min(abs(times-2))]

OP_IS$blood_4[i] = blood[which.min(abs(times-4))]

OP_IS$blood_8[i] = blood[which.min(abs(times-8))]

OP_IS$blood_10[i] = blood[which.min(abs(times-10))]

OP_IS$blood_24[i] = blood[which.min(abs(times-24))]

OP_IS$blood_48[i] = blood[which.min(abs(times-48))]

OP_IS$Heart_3[i] = HV[which.min(abs(times-3))]

OP_IS$Kidneys_3[i] = KV[which.min(abs(times-3))]

OP_IS$Spleen_3[i] = SV[which.min(abs(times-3))]

OP_IS$Lungs_3[i] = LuV[which.min(abs(times-3))]

OP_IS$Bone_3[i] = BoneV[which.min(abs(times-3))]

OP_IS$Heart_48[i] = HV[which.min(abs(times-48))]

OP_IS$Kidneys_48[i] = KV[which.min(abs(times-48))]

OP_IS$Spleen_48[i] = SV[which.min(abs(times-48))]

OP_IS$Lungs_48[i] = LuV[which.min(abs(times-48))]

OP_IS$Bone_48[i] = BoneV[which.min(abs(times-48))]
}

OP_IS$cond.blood_03_yn <- 0*OP_IS$nt
OP_IS$cond.blood_03_yn[OP_IS$blood_03 < blood_summary$UCI[blood_summary$time==0.03] & 
                         OP_IS$blood_03 > blood_summary$LCI[blood_summary$time==0.03]] <- 1

OP_IS$cond.blood_1_yn <- 0*OP_IS$nt
OP_IS$cond.blood_1_yn[OP_IS$blood_1 < blood_summary$UCI[blood_summary$time==1] & 
                        OP_IS$blood_1 > blood_summary$LCI[blood_summary$time==1]] <- 1

OP_IS$cond.blood_2_yn <- 0*OP_IS$nt
OP_IS$cond.blood_2_yn[OP_IS$blood_2 < blood_summary$UCI[blood_summary$time==2] & 
                        OP_IS$blood_2 > blood_summary$LCI[blood_summary$time==2]] <- 1

OP_IS$cond.blood_4_yn <- 0*OP_IS$nt
OP_IS$cond.blood_4_yn[OP_IS$blood_4 < blood_summary$UCI[blood_summary$time==4] & 
                        OP_IS$blood_4 > blood_summary$LCI[blood_summary$time==4]] <- 1

OP_IS$cond.blood_8_yn <- 0*OP_IS$nt
OP_IS$cond.blood_8_yn[OP_IS$blood_8 < blood_summary$UCI[blood_summary$time==8] & 
                        OP_IS$blood_8 > blood_summary$LCI[blood_summary$time==8]] <- 1

OP_IS$cond.blood_10_yn <- 0*OP_IS$nt
OP_IS$cond.blood_10_yn[OP_IS$blood_10 < blood_summary$UCI[blood_summary$time==10] & 
                        OP_IS$blood_10 > blood_summary$LCI[blood_summary$time==10]] <- 1

OP_IS$cond.blood_24_yn <- 0*OP_IS$nt
OP_IS$cond.blood_24_yn[OP_IS$blood_24< blood_summary$UCI[blood_summary$time==24] & 
                         OP_IS$blood_24 > blood_summary$LCI[blood_summary$time==24]] <- 1

OP_IS$cond.blood_48_yn <- 0*OP_IS$nt
OP_IS$cond.blood_48_yn[OP_IS$blood_48 < blood_summary$UCI[blood_summary$time==48] & 
                         OP_IS$blood_48 > blood_summary$LCI[blood_summary$time==48]] <- 1

OP_IS$cond.blood.tot <- OP_IS$cond.blood_1_yn + OP_IS$cond.blood_8_yn + OP_IS$cond.blood_24_yn +
                        OP_IS$cond.blood_48_yn + OP_IS$cond.blood_03_yn  

OP_IS$cond.Heart_3_yn <- 0*OP_IS$nt
OP_IS$cond.Heart_3_yn[OP_IS$Heart_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Heart'] & 
                         OP_IS$Heart_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Heart']] <- 1

OP_IS$cond.Kidneys_3_yn <- 0*OP_IS$nt
OP_IS$cond.Kidneys_3_yn[OP_IS$Kidneys_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Kidneys'] & 
                           OP_IS$Kidneys_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Kidneys']] <- 1

OP_IS$cond.Lungs_3_yn <- 0*OP_IS$nt
OP_IS$cond.Lungs_3_yn[OP_IS$Lungs_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Lungs'] & 
                         OP_IS$Lungs_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Lungs']] <- 1

OP_IS$cond.Bone_3_yn <- 0*OP_IS$nt
OP_IS$cond.Bone_3_yn[OP_IS$Bone_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Bone'] & 
                        OP_IS$Bone_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Bone']] <- 1

OP_IS$cond.Spleen_3_yn <- 0*OP_IS$nt
OP_IS$cond.Spleen_3_yn[OP_IS$Spleen_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Spleen'] & 
                          OP_IS$Spleen_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Spleen']] <- 1

OP_IS$cond.organ.tot_3 <- OP_IS$cond.Spleen_3_yn + OP_IS$cond.Heart_3_yn + 
  OP_IS$cond.Kidneys_3_yn + OP_IS$cond.Lungs_3_yn + OP_IS$cond.Bone_3_yn

OP_IS$cond.Heart_48_yn <- 0*OP_IS$nt
OP_IS$cond.Heart_48_yn[OP_IS$Heart_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Heart'] & 
                      OP_IS$Heart_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Heart']] <- 1

OP_IS$cond.Kidneys_48_yn <- 0*OP_IS$nt
OP_IS$cond.Kidneys_48_yn[OP_IS$Kidneys_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Kidneys'] & 
                        OP_IS$Kidneys_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Kidneys']] <- 1

OP_IS$cond.Lungs_48_yn <- 0*OP_IS$nt
OP_IS$cond.Lungs_48_yn[OP_IS$Lungs_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Lungs'] & 
                      OP_IS$Lungs_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Lungs']] <- 1

OP_IS$cond.Bone_48_yn <- 0*OP_IS$nt
OP_IS$cond.Bone_48_yn[OP_IS$Bone_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Bone'] & 
                     OP_IS$Bone_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Bone']] <- 1

OP_IS$cond.Spleen_48_yn <- 0*OP_IS$nt
OP_IS$cond.Spleen_48_yn[OP_IS$Spleen_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Spleen'] & 
                       OP_IS$Spleen_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Spleen']] <- 1

OP_IS$cond.organ.tot_48 <- OP_IS$cond.Spleen_48_yn + OP_IS$cond.Heart_48_yn + 
  OP_IS$cond.Kidneys_48_yn + OP_IS$cond.Lungs_48_yn + OP_IS$cond.Bone_48_yn

OP_IS$PACt_Liver_scale <- log10(OP_IS$PACt_Liver/OP_IS$Kbile)
OP_IS$Kbile_scale <- log10(OP_IS$Kbile/Kbile_Au)

OP_IS$PACt_Bone_scale <- log10(OP_IS$PACt_Bone/PACt_Bone_Au)
OP_IS$PACt_Heart_scale <- log10(OP_IS$PACt_Heart/PACt_Heart_Au)
OP_IS$PACt_Spleen_scale <- log10(OP_IS$PACt_Spleen/PACt_Spleen_Au)
OP_IS$PACt_Lungs_scale <- log10(OP_IS$PACt_Lungs/PACt_Lungs_Au)
OP_IS$PACt_Rest_scale <- log10(OP_IS$PACt_Rest/PACt_Rest_Au)

OP_IS$Kmax_t_Liver_scale <- log10(OP_IS$Kmax_t_Liver)
OP_IS$Kmax_t_Spleen_scale <- log10(OP_IS$Kmax_t_Spleen)
OP_IS$Kmax_t_Kidneys_scale <- log10(OP_IS$Kmax_t_Kidneys)
OP_IS$Kmax_t_Lungs_scale <- log10(OP_IS$Kmax_t_Lungs)

OP_IS$Krelease_t_Liver_scale <- log10(OP_IS$Krelease_t_Liver/OP_IS$Kmax_t_Liver)
OP_IS$Krelease_t_Spleen_scale <- log10(OP_IS$Krelease_t_Spleen/OP_IS$Kmax_t_Spleen)
OP_IS$Krelease_t_Kidneys_scale <- log10(OP_IS$Krelease_t_Kidneys/OP_IS$Kmax_t_Kidneys)
OP_IS$Krelease_t_Lungs_scale <- log10(OP_IS$Krelease_t_Lungs/OP_IS$Kmax_t_Lungs)

OP_IS$K50_t_Liver_scale <- log10(OP_IS$K50_t_Liver/K50_t_Liver_Au)
OP_IS$K50_t_Spleen_scale <- log10(OP_IS$K50_t_Spleen/K50_t_Spleen_Au)
OP_IS$K50_t_Kidneys_scale <- log10(OP_IS$K50_t_Kidneys/K50_t_Kidneys_Au)
OP_IS$K50_t_Lungs_scale <- log10(OP_IS$K50_t_Lungs/K50_t_Lungs_Au)

OP_IS$nt_scale <- log10(OP_IS$nt/nt_Au)

OP_IS$dose <- rep(dose_rel[j],Nit)

OP_IS_nana <- OP_IS[complete.cases(OP_IS),]

if (j==1){
  
  OP_IS_tot <- OP_IS_nana
  
}else{
  
  OP_IS_tot <- rbind(OP_IS_tot,OP_IS_nana)
  
}

}

saveRDS(OP_IS_tot,"model_generated_data/scale_MC_20240118.rds")

OP_IS_max_cond <- OP_IS_tot[OP_IS_tot$cond.blood.tot%in%c(5),]
OP_IS_max_cond

OP_IS_max_cond_organs <- OP_IS_max_cond[OP_IS_max_cond$cond.organ.tot_3%in%c(3,4,5) &
                                          OP_IS_max_cond$cond.organ.tot_48%in%c(3,4,5),]

OPISmax_Heart <- OP_IS_max_cond[OP_IS_max_cond$cond.Heart_3_yn==1,]# & OP_IS_max_cond$cond.Heart_48_yn, ]
OPISmax_Lungs <- OP_IS_max_cond[OP_IS_max_cond$cond.Lungs_3_yn==1 & OP_IS_max_cond$cond.Lungs_48_yn==1, ]
OPISmax_Kidneys <- OP_IS_max_cond[OP_IS_max_cond$cond.Kidneys_3_yn==1 & OP_IS_max_cond$cond.Kidneys_48_yn==1, ]

scales <- data.frame(
  PACt_Liver=  mean(OP_IS_max_cond$PACt_Liver_scale),
  Kbile= mean(OP_IS_max_cond$Kbile_scale),

  PACt_Heart= mean(OPISmax_Heart$PACt_Heart_scale),
  PACt_Lungs= mean(OPISmax_Lungs$PACt_Lungs_scale),
  PACt_Bone= mean(OP_IS_max_cond$PACt_Bone_scale),
  PACt_Spleen= mean(OP_IS_max_cond$PACt_Spleen_scale),
  PACt_Rest= mean(OP_IS_max_cond$PACt_Rest_scale),

  Kmax_t_Liver= mean(OP_IS_max_cond$Kmax_t_Liver_scale),
  Kmax_t_Spleen= mean(OP_IS_max_cond$Kmax_t_Spleen_scale),
  Kmax_t_Lungs= mean(OPISmax_Lungs$Kmax_t_Lungs_scale),
  Kmax_t_Kidneys= mean(OPISmax_Kidneys$Kmax_t_Kidneys_scale),

  Krelease_t_Liver= mean(OP_IS_max_cond$Krelease_t_Liver_scale),
  Krelease_t_Spleen= mean(OP_IS_max_cond$Krelease_t_Spleen_scale),
  Krelease_t_Lungs= mean(OPISmax_Lungs$Krelease_t_Lungs_scale),
  Krelease_t_Kidneys= mean(OPISmax_Kidneys$Krelease_t_Kidneys_scale),

  K50_t_Liver= mean(OP_IS_max_cond$K50_t_Liver_scale),
  K50_t_Spleen= mean(OP_IS_max_cond$K50_t_Spleen_scale),
  K50_t_Lungs= mean(OPISmax_Lungs$K50_t_Lungs_scale),
  K50_t_Kidneys= mean(OPISmax_Kidneys$K50_t_Kidneys_scale),

  nt = mean(OP_IS_max_cond$nt_scale)
)

saveRDS(scales,"model_generated_data/scales_20240118.rds")


########################################################
#Third section: rescale, unrounded

#Clear all variables
rm(list=ls())

set.seed(7)

#Source scripts
source("libs.R")
source("parms.R")
alpha_temp_FACS <- 1e-2
alpha_temp_blood <- 1e-2
source("funcs.R")

OP_IS_tot <- readRDS("model_generated_data/scale_MC_20240118.rds")

OP_IS_max_cond <- OP_IS_tot[OP_IS_tot$cond.blood.tot%in%c(5),]
OP_IS_max_cond

OP_IS_max_cond_organs <- OP_IS_max_cond[OP_IS_max_cond$cond.organ.tot_3%in%c(3,4,5) &
                                          OP_IS_max_cond$cond.organ.tot_48%in%c(3,4,5),]
OP_IS_max_cond_organs
OPISmax_Heart <- OP_IS_max_cond[OP_IS_max_cond$cond.Heart_3_yn==1,]# & OP_IS_max_cond$cond.Heart_48_yn, ]
OPISmax_Lungs <- OP_IS_max_cond[OP_IS_max_cond$cond.Lungs_3_yn==1 & OP_IS_max_cond$cond.Lungs_48_yn==1, ]
OPISmax_Kidneys <- OP_IS_max_cond[OP_IS_max_cond$cond.Kidneys_3_yn==1 & OP_IS_max_cond$cond.Kidneys_48_yn==1, ]

scales <- data.frame(
  PACt_Liver=  mean(OP_IS_max_cond_organs$PACt_Liver_scale),
  Kbile= mean(OP_IS_max_cond_organs$Kbile_scale),
  
  PACt_Heart= mean(OPISmax_Heart$PACt_Heart_scale),
  PACt_Lungs= mean(OPISmax_Lungs$PACt_Lungs_scale),
  PACt_Bone= mean(OP_IS_max_cond_organs$PACt_Bone_scale),
  PACt_Spleen= mean(OP_IS_max_cond_organs$PACt_Spleen_scale),
  PACt_Rest= mean(OP_IS_max_cond_organs$PACt_Rest_scale),
  
  Kmax_t_Liver= mean(OP_IS_max_cond_organs$Kmax_t_Liver_scale),
  Kmax_t_Spleen= mean(OP_IS_max_cond_organs$Kmax_t_Spleen_scale),
  Kmax_t_Lungs= mean(OPISmax_Lungs$Kmax_t_Lungs_scale),
  Kmax_t_Kidneys= mean(OPISmax_Kidneys$Kmax_t_Kidneys_scale),
  
  Krelease_t_Liver= mean(OP_IS_max_cond_organs$Krelease_t_Liver_scale),
  Krelease_t_Spleen= mean(OP_IS_max_cond_organs$Krelease_t_Spleen_scale),
  Krelease_t_Lungs= mean(OPISmax_Lungs$Krelease_t_Lungs_scale),
  Krelease_t_Kidneys= mean(OPISmax_Kidneys$Krelease_t_Kidneys_scale),
  
  K50_t_Liver= mean(OP_IS_max_cond_organs$K50_t_Liver_scale),
  K50_t_Spleen= mean(OP_IS_max_cond_organs$K50_t_Spleen_scale),
  K50_t_Lungs= mean(OPISmax_Lungs$K50_t_Lungs_scale),
  K50_t_Kidneys= mean(OPISmax_Kidneys$K50_t_Kidneys_scale),
  
  nt = mean(OP_IS_max_cond_organs$nt_scale)
)


Nit <- 1e4
PACt_mod <- 1#e-3
Kmax_Spleen_mod <- 1#e-2

for (j in 1:length(dose_rel)){
  
  source("FACS_IVIS_data_analysis_20230426.R")
  
  blood_conc <- blood_conc[blood_conc$dose==dose_rel[j],]
  blood_summary <- blood_summary[blood_summary$dose==dose_rel[j],]
  organ_FACS_summary <- organ_FACS_summary[organ_FACS_summary$Dose==dose_rel[j],]
  organ_FACS_time48_summary <- organ_FACS_summary
  organ_FACS_time3_summary <- organ_FACS_time_summary[organ_FACS_time_summary$time==3,]
  
  dose_rate <- dose_rate_vec[j]
  
  Init_condition <- rep(0,22)
  
  OP_IS <- data.frame(
    dose = rep(dose_rel[j],Nit),
    
    nt = rep(0,Nit),
    Kbile = rep(0,Nit),
    Kurine = rep(0,Nit),
    
    PACt_Liver = rep(0,Nit),
    PACt_Kidneys = rep(0,Nit),
    PACt_Spleen = rep(0,Nit),
    PACt_Lungs = rep(0,Nit),
    PACt_Heart = rep(0,Nit),
    PACt_Bone = rep(0,Nit),
    PACt_Rest = rep(0,Nit),
    
    Kmax_t_Liver = rep(0,Nit),
    Kmax_t_Spleen  = rep(0,Nit),
    Kmax_t_Kidneys  = rep(0,Nit),
    Kmax_t_Lungs  = rep(0,Nit),
    
    Krelease_t_Liver  = rep(0,Nit),
    Krelease_t_Spleen  = rep(0,Nit),
    Krelease_t_Kidneys  = rep(0,Nit),
    Krelease_t_Lungs = rep(0,Nit),
    
    K50_t_Liver  = rep(0,Nit),
    K50_t_Spleen  = rep(0,Nit),
    K50_t_Kidneys  = rep(0,Nit),
    K50_t_Lungs = rep(0,Nit),
    
    blood_03 = rep(0,Nit),
    
    blood_2 = rep(0,Nit),
    
    blood_4 = rep(0,Nit),
    
    blood_8 = rep(0,Nit),
    
    blood_10 = rep(0,Nit),
    
    blood_24 = rep(0,Nit),
    
    blood_48 = rep(0,Nit),
    
    Heart_3 = rep(0,Nit),
    
    Kidneys_3 = rep(0,Nit),
    
    Spleen_3 = rep(0,Nit),
    
    Lungs_3 = rep(0,Nit),
    
    Bone_3 = rep(0,Nit),
    
    Heart_48 = rep(0,Nit),
    
    Kidneys_48 = rep(0,Nit),
    
    Spleen_48 = rep(0,Nit),
    
    Lungs_48 = rep(0,Nit),
    
    Bone_48 = rep(0,Nit)
  )
  
  for (i in seq(Nit)){
    
    nt <- nt_Au*10^scales$nt*10^runif(1,-1,1)
    Kbile <- Kbile_Au*10^scales$Kbile*10^rnorm(1,0,0.2)
    Kurine <- 0
    PACt_Kidneys <- 0
    
    PACt_Lungs <- PACt_Lungs_Au*10^scales$PACt_Lungs*10^runif(1,-1,1)
    PACt_Heart <- PACt_Heart_Au*10^scales$PACt_Heart*10^runif(1,-2,0)
    PACt_Bone <- PACt_Bone_Au*10^scales$PACt_Bone*10^runif(1,-1,1)
    PACt_Spleen <- PACt_Spleen_Au*10^scales$PACt_Spleen*10^runif(1,-1,1)
    PACt_Rest <- PACt_Rest_Au*10^scales$PACt_Rest*10^runif(1,-1,1)
    PACt_Brain <- PACt_Rest
    
    Kmax_t_Liver <- 10^scales$Kmax_t_Liver*10^rnorm(1,0,0.2)
    Kmax_t_Kidneys <- 10^scales$Kmax_t_Kidneys*10^runif(1,-1,1)
    Kmax_t_Spleen <- 10^scales$Kmax_t_Spleen*10^runif(1,-1,1)
    Kmax_t_Lungs <- 10^scales$Kmax_t_Lungs*10^runif(1,-1,1)
    
    K50_t_Liver <- K50_t_Liver_Au*10^scales$K50_t_Liver*10^rnorm(1,0,0.2)
    K50_t_Spleen <- K50_t_Spleen_Au*10^scales$K50_t_Spleen*10^runif(1,-1,1)
    K50_t_Kidneys <- K50_t_Kidneys_Au*10^scales$K50_t_Kidneys*10^runif(1,-1,1)
    K50_t_Lungs <- K50_t_Lungs_Au*10^scales$K50_t_Lungs*10^runif(1,-1,1)
    
    Krelease_t_Liver <- Kmax_t_Liver*10^scales$Krelease_t_Liver#*10^rnorm(1,0,0.0001)
    Krelease_t_Spleen <- Kmax_t_Spleen*10^scales$Krelease_t_Spleen#*10^rnorm(1,0,0.0001)
    Krelease_t_Kidneys <- Kmax_t_Kidneys*10^scales$Krelease_t_Kidneys#*10^rnorm(1,0,0.0001)
    Krelease_t_Lungs <- Kmax_t_Lungs*10^scales$Krelease_t_Lungs#*10^rnorm(1,0,0.0001)
    
    PACt_Liver <- Kbile*10^scales$PACt_Liver#*10^rnorm(1,0,0.0001)
    
    
    OP <- ode(y=Init_condition,
              times=times,
              func = NP_PBPK2, 
              parms=list(),
              method='lsoda')
    OP <- data.frame(OP)
    OP_rename <- OP
    
    indx_samp <- seq(from = 1, to = Nt, by = 1)
    
    OP_rename$M_Art <- OP[indx_samp,2]
    OP_rename$M_Ven <- OP[indx_samp,3]
    
    OP_rename$M_vasc_Liver <- OP[indx_samp,4]
    OP_rename$M_extra_Liver <- OP[indx_samp,5]
    OP_rename$M_phago_Liver <- OP[indx_samp,6]
    
    OP_rename$M_vasc_Spleen <- OP[indx_samp,7]
    OP_rename$M_extra_Spleen <- OP[indx_samp,8]
    OP_rename$M_phago_Spleen <- OP[indx_samp,9]
    
    OP_rename$M_vasc_Kidneys <- OP[indx_samp,10]
    OP_rename$M_extra_Kidneys <- OP[indx_samp,11]
    OP_rename$M_phago_Kidneys <- OP[indx_samp,12]
    
    OP_rename$M_vasc_Lungs <- OP[indx_samp,13]
    OP_rename$M_extra_Lungs <- OP[indx_samp,14]
    OP_rename$M_phago_Lungs <- OP[indx_samp,15]
    
    OP_rename$M_vasc_Brain <- OP[indx_samp,16]
    OP_rename$M_extra_Brain <- OP[indx_samp,17]
    
    OP_rename$M_vasc_Heart <- OP[indx_samp,18]
    OP_rename$M_extra_Heart <- OP[indx_samp,19]
    
    OP_rename$M_vasc_Bone <- OP[indx_samp,20]
    OP_rename$M_extra_Bone <- OP[indx_samp,21]
    
    OP_rename$M_vasc_Rest <- OP[indx_samp,22]
    OP_rename$M_extra_Rest <- OP[indx_samp,23]
    
    LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
    KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/LV
    LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/LV
    SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/LV
    HV <- OP_rename$M_extra_Heart/LV
    BoneV <- OP_rename$M_extra_Bone/LV
    LV <- 1
    
    blood <- (OP_rename$M_Ven)/(BW*VVenC)
    
    plot(times,blood,type='l',ylim=c(0, 1500))
    points(blood_summary$time, blood_summary$LCI,col='blue')
    points(blood_summary$time, blood_summary$UCI,col='red')
    
    OP_IS$PACt_Liver[i] <- PACt_Liver
    OP_IS$Kbile[i] = Kbile
    
    OP_IS$PACt_Kidneys[i] <- PACt_Kidneys
    OP_IS$Kurine[i] = Kurine
    
    OP_IS$PACt_Bone[i] <- PACt_Bone
    OP_IS$PACt_Spleen[i] <- PACt_Spleen
    OP_IS$PACt_Lungs[i] <- PACt_Lungs
    OP_IS$PACt_Heart[i] <- PACt_Heart
    OP_IS$PACt_Rest[i] <- PACt_Rest
    
    OP_IS$nt[i] = nt
    
    OP_IS$Kmax_t_Liver[i] = Kmax_t_Liver 
    OP_IS$Kmax_t_Spleen[i]  = Kmax_t_Spleen
    OP_IS$Kmax_t_Kidneys[i]  = Kmax_t_Kidneys 
    OP_IS$Kmax_t_Lungs[i]  = Kmax_t_Lungs
    
    OP_IS$Krelease_t_Liver[i]  = Krelease_t_Liver
    OP_IS$Krelease_t_Spleen[i]  = Krelease_t_Spleen
    OP_IS$Krelease_t_Kidneys[i]  = Krelease_t_Kidneys
    OP_IS$Krelease_t_Lungs[i] = Krelease_t_Lungs
    
    OP_IS$K50_t_Liver[i]  = K50_t_Liver
    OP_IS$K50_t_Spleen[i]  = K50_t_Spleen
    OP_IS$K50_t_Kidneys[i]  = K50_t_Kidneys
    OP_IS$K50_t_Lungs[i] = K50_t_Lungs
    
    OP_IS$blood_03[i] = blood[which.min(abs(times-0.03))]
    
    OP_IS$blood_1[i] = blood[which.min(abs(times-1))]
    
    OP_IS$blood_2[i] = blood[which.min(abs(times-2))]
    
    OP_IS$blood_4[i] = blood[which.min(abs(times-4))]
    
    OP_IS$blood_8[i] = blood[which.min(abs(times-8))]
    
    OP_IS$blood_10[i] = blood[which.min(abs(times-10))]
    
    OP_IS$blood_24[i] = blood[which.min(abs(times-24))]
    
    OP_IS$blood_48[i] = blood[which.min(abs(times-48))]
    
    OP_IS$Heart_3[i] = HV[which.min(abs(times-3))]
    
    OP_IS$Kidneys_3[i] = KV[which.min(abs(times-3))]
    
    OP_IS$Spleen_3[i] = SV[which.min(abs(times-3))]
    
    OP_IS$Lungs_3[i] = LuV[which.min(abs(times-3))]
    
    OP_IS$Bone_3[i] = BoneV[which.min(abs(times-3))]
    
    OP_IS$Heart_48[i] = HV[which.min(abs(times-48))]
    
    OP_IS$Kidneys_48[i] = KV[which.min(abs(times-48))]
    
    OP_IS$Spleen_48[i] = SV[which.min(abs(times-48))]
    
    OP_IS$Lungs_48[i] = LuV[which.min(abs(times-48))]
    
    OP_IS$Bone_48[i] = BoneV[which.min(abs(times-48))]
  }
  
  OP_IS$cond.blood_03_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_03_yn[OP_IS$blood_03 < blood_summary$UCI[blood_summary$time==0.03] & 
                           OP_IS$blood_03 > blood_summary$LCI[blood_summary$time==0.03]] <- 1
  
  OP_IS$cond.blood_1_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_1_yn[OP_IS$blood_1 < blood_summary$UCI[blood_summary$time==1] & 
                          OP_IS$blood_1 > blood_summary$LCI[blood_summary$time==1]] <- 1
  
  OP_IS$cond.blood_2_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_2_yn[OP_IS$blood_2 < blood_summary$UCI[blood_summary$time==2] & 
                          OP_IS$blood_2 > blood_summary$LCI[blood_summary$time==2]] <- 1
  
  OP_IS$cond.blood_4_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_4_yn[OP_IS$blood_4 < blood_summary$UCI[blood_summary$time==4] & 
                          OP_IS$blood_4 > blood_summary$LCI[blood_summary$time==4]] <- 1
  
  OP_IS$cond.blood_8_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_8_yn[OP_IS$blood_8 < blood_summary$UCI[blood_summary$time==8] & 
                          OP_IS$blood_8 > blood_summary$LCI[blood_summary$time==8]] <- 1
  
  OP_IS$cond.blood_10_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_10_yn[OP_IS$blood_10 < blood_summary$UCI[blood_summary$time==10] & 
                           OP_IS$blood_10 > blood_summary$LCI[blood_summary$time==10]] <- 1
  
  OP_IS$cond.blood_24_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_24_yn[OP_IS$blood_24< blood_summary$UCI[blood_summary$time==24] & 
                           OP_IS$blood_24 > blood_summary$LCI[blood_summary$time==24]] <- 1
  
  OP_IS$cond.blood_48_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_48_yn[OP_IS$blood_48 < blood_summary$UCI[blood_summary$time==48] & 
                           OP_IS$blood_48 > blood_summary$LCI[blood_summary$time==48]] <- 1
  
  OP_IS$cond.blood.tot <- OP_IS$cond.blood_1_yn + OP_IS$cond.blood_8_yn + OP_IS$cond.blood_24_yn +
    OP_IS$cond.blood_48_yn + OP_IS$cond.blood_03_yn  
  
  OP_IS$cond.Heart_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Heart_3_yn[OP_IS$Heart_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Heart'] & 
                          OP_IS$Heart_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Heart']] <- 1
  
  OP_IS$cond.Kidneys_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Kidneys_3_yn[OP_IS$Kidneys_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Kidneys'] & 
                            OP_IS$Kidneys_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Kidneys']] <- 1
  
  OP_IS$cond.Lungs_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Lungs_3_yn[OP_IS$Lungs_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Lungs'] & 
                          OP_IS$Lungs_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Lungs']] <- 1
  
  OP_IS$cond.Bone_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Bone_3_yn[OP_IS$Bone_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Bone'] & 
                         OP_IS$Bone_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Bone']] <- 1
  
  OP_IS$cond.Spleen_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Spleen_3_yn[OP_IS$Spleen_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Spleen'] & 
                           OP_IS$Spleen_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Spleen']] <- 1
  
  OP_IS$cond.organ.tot_3 <- OP_IS$cond.Spleen_3_yn + OP_IS$cond.Heart_3_yn + 
    OP_IS$cond.Kidneys_3_yn + OP_IS$cond.Lungs_3_yn + OP_IS$cond.Bone_3_yn
  
  OP_IS$cond.Heart_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Heart_48_yn[OP_IS$Heart_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Heart'] & 
                           OP_IS$Heart_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Heart']] <- 1
  
  OP_IS$cond.Kidneys_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Kidneys_48_yn[OP_IS$Kidneys_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Kidneys'] & 
                             OP_IS$Kidneys_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Kidneys']] <- 1
  
  OP_IS$cond.Lungs_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Lungs_48_yn[OP_IS$Lungs_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Lungs'] & 
                           OP_IS$Lungs_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Lungs']] <- 1
  
  OP_IS$cond.Bone_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Bone_48_yn[OP_IS$Bone_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Bone'] & 
                          OP_IS$Bone_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Bone']] <- 1
  
  OP_IS$cond.Spleen_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Spleen_48_yn[OP_IS$Spleen_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Spleen'] & 
                            OP_IS$Spleen_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Spleen']] <- 1
  
  OP_IS$cond.organ.tot_48 <- OP_IS$cond.Spleen_48_yn + OP_IS$cond.Heart_48_yn + 
    OP_IS$cond.Kidneys_48_yn + OP_IS$cond.Lungs_48_yn + OP_IS$cond.Bone_48_yn
  
  OP_IS$PACt_Liver_scale <- log10(OP_IS$PACt_Liver/OP_IS$Kbile)
  OP_IS$Kbile_scale <- log10(OP_IS$Kbile/Kbile_Au)
  
  OP_IS$PACt_Bone_scale <- log10(OP_IS$PACt_Bone/PACt_Bone_Au)
  OP_IS$PACt_Heart_scale <- log10(OP_IS$PACt_Heart/PACt_Heart_Au)
  OP_IS$PACt_Spleen_scale <- log10(OP_IS$PACt_Spleen/PACt_Spleen_Au)
  OP_IS$PACt_Lungs_scale <- log10(OP_IS$PACt_Lungs/PACt_Lungs_Au)
  OP_IS$PACt_Rest_scale <- log10(OP_IS$PACt_Rest/PACt_Rest_Au)
  
  OP_IS$Kmax_t_Liver_scale <- log10(OP_IS$Kmax_t_Liver)
  OP_IS$Kmax_t_Spleen_scale <- log10(OP_IS$Kmax_t_Spleen)
  OP_IS$Kmax_t_Kidneys_scale <- log10(OP_IS$Kmax_t_Kidneys)
  OP_IS$Kmax_t_Lungs_scale <- log10(OP_IS$Kmax_t_Lungs)
  
  OP_IS$Krelease_t_Liver_scale <- log10(OP_IS$Krelease_t_Liver/OP_IS$Kmax_t_Liver)
  OP_IS$Krelease_t_Spleen_scale <- log10(OP_IS$Krelease_t_Spleen/OP_IS$Kmax_t_Spleen)
  OP_IS$Krelease_t_Kidneys_scale <- log10(OP_IS$Krelease_t_Kidneys/OP_IS$Kmax_t_Kidneys)
  OP_IS$Krelease_t_Lungs_scale <- log10(OP_IS$Krelease_t_Lungs/OP_IS$Kmax_t_Lungs)
  
  OP_IS$K50_t_Liver_scale <- log10(OP_IS$K50_t_Liver/K50_t_Liver_Au)
  OP_IS$K50_t_Spleen_scale <- log10(OP_IS$K50_t_Spleen/K50_t_Spleen_Au)
  OP_IS$K50_t_Kidneys_scale <- log10(OP_IS$K50_t_Kidneys/K50_t_Kidneys_Au)
  OP_IS$K50_t_Lungs_scale <- log10(OP_IS$K50_t_Lungs/K50_t_Lungs_Au)
  
  OP_IS$nt_scale <- log10(OP_IS$nt/nt_Au)
  
  OP_IS$dose <- rep(dose_rel[j],Nit)
  
  OP_IS_nana <- OP_IS[complete.cases(OP_IS),]
  
  if (j==1){
    
    OP_IS_tot <- OP_IS_nana
    
  }else{
    
    OP_IS_tot <- rbind(OP_IS_tot,OP_IS_nana)
    
  }
  
}

saveRDS(OP_IS_tot,"model_generated_data/parm_MC_20240119.rds")

OP_IS_max_cond <- OP_IS_tot[OP_IS_tot$cond.blood.tot%in%c(5),]
OP_IS_max_cond

OP_IS_max_cond_organs <- OP_IS_max_cond[OP_IS_max_cond$cond.organ.tot_3%in%c(3,4,5) &
                                          OP_IS_max_cond$cond.organ.tot_48%in%c(4,5),]
OPISmax_Heart <- OP_IS_max_cond[OP_IS_max_cond$cond.Heart_3_yn==1,]# & OP_IS_max_cond$cond.Heart_48_yn, ]
OPISmax_Lungs <- OP_IS_max_cond[OP_IS_max_cond$cond.Lungs_3_yn==1 & OP_IS_max_cond$cond.Lungs_48_yn==1, ]
OPISmax_Kidneys <- OP_IS_max_cond[OP_IS_max_cond$cond.Kidneys_3_yn==1 & OP_IS_max_cond$cond.Kidneys_48_yn==1, ]

scales <- data.frame(
  PACt_Liver=  mean(OP_IS_max_cond_organs$PACt_Liver_scale),
  Kbile= mean(OP_IS_max_cond_organs$Kbile_scale),
  
  PACt_Heart= mean(OPISmax_Heart$PACt_Heart_scale),
  PACt_Lungs= mean(OPISmax_Lungs$PACt_Lungs_scale),
  PACt_Bone= mean(OP_IS_max_cond_organs$PACt_Bone_scale),
  PACt_Spleen= mean(OP_IS_max_cond_organs$PACt_Spleen_scale),
  PACt_Rest= mean(OP_IS_max_cond_organs$PACt_Rest_scale),
  
  Kmax_t_Liver= mean(OP_IS_max_cond_organs$Kmax_t_Liver_scale),
  Kmax_t_Spleen= mean(OP_IS_max_cond_organs$Kmax_t_Spleen_scale),
  Kmax_t_Lungs= mean(OPISmax_Lungs$Kmax_t_Lungs_scale),
  Kmax_t_Kidneys= mean(OPISmax_Kidneys$Kmax_t_Kidneys_scale),
  
  Krelease_t_Liver= mean(OP_IS_max_cond_organs$Krelease_t_Liver_scale),
  Krelease_t_Spleen= mean(OP_IS_max_cond_organs$Krelease_t_Spleen_scale),
  Krelease_t_Lungs= mean(OPISmax_Lungs$Krelease_t_Lungs_scale),
  Krelease_t_Kidneys= mean(OPISmax_Kidneys$Krelease_t_Kidneys_scale),
  
  K50_t_Liver= mean(OP_IS_max_cond_organs$K50_t_Liver_scale),
  K50_t_Spleen= mean(OP_IS_max_cond_organs$K50_t_Spleen_scale),
  K50_t_Lungs= mean(OPISmax_Lungs$K50_t_Lungs_scale),
  K50_t_Kidneys= mean(OPISmax_Kidneys$K50_t_Kidneys_scale),
  
  nt = mean(OP_IS_max_cond_organs$nt_scale)
)

saveRDS(scales,"model_generated_data/scales_20240119.rds")

########################################################
#Fourth section: rescale, with normal distributions

#Clear all variables
rm(list=ls())

set.seed(21)

#Source scripts
source("libs.R")
source("parms.R")
alpha_temp_FACS <- 1e-2
alpha_temp_blood <- 1e-2
source("funcs.R")

scales <- readRDS("model_generated_data/scales_20240119.rds")


Nit <- 1e4
PACt_mod <- 1#e-3
Kmax_Spleen_mod <- 1#e-2

for (j in 1:length(dose_rel)){
  
  source("FACS_IVIS_data_analysis_20230426.R")
  
  blood_conc <- blood_conc[blood_conc$dose==dose_rel[j],]
  blood_summary <- blood_summary[blood_summary$dose==dose_rel[j],]
  organ_FACS_summary <- organ_FACS_summary[organ_FACS_summary$Dose==dose_rel[j],]
  organ_FACS_time48_summary <- organ_FACS_summary
  organ_FACS_time3_summary <- organ_FACS_time_summary[organ_FACS_time_summary$time==3,]
  
  dose_rate <- dose_rate_vec[j]
  
  Init_condition <- rep(0,22)
  
  OP_IS <- data.frame(
    dose = rep(dose_rel[j],Nit),
    
    nt = rep(0,Nit),
    Kbile = rep(0,Nit),
    Kurine = rep(0,Nit),
    
    PACt_Liver = rep(0,Nit),
    PACt_Kidneys = rep(0,Nit),
    PACt_Spleen = rep(0,Nit),
    PACt_Lungs = rep(0,Nit),
    PACt_Heart = rep(0,Nit),
    PACt_Bone = rep(0,Nit),
    PACt_Rest = rep(0,Nit),
    
    Kmax_t_Liver = rep(0,Nit),
    Kmax_t_Spleen  = rep(0,Nit),
    Kmax_t_Kidneys  = rep(0,Nit),
    Kmax_t_Lungs  = rep(0,Nit),
    
    Krelease_t_Liver  = rep(0,Nit),
    Krelease_t_Spleen  = rep(0,Nit),
    Krelease_t_Kidneys  = rep(0,Nit),
    Krelease_t_Lungs = rep(0,Nit),
    
    K50_t_Liver  = rep(0,Nit),
    K50_t_Spleen  = rep(0,Nit),
    K50_t_Kidneys  = rep(0,Nit),
    K50_t_Lungs = rep(0,Nit),
    
    blood_03 = rep(0,Nit),
    
    blood_2 = rep(0,Nit),
    
    blood_4 = rep(0,Nit),
    
    blood_8 = rep(0,Nit),
    
    blood_10 = rep(0,Nit),
    
    blood_24 = rep(0,Nit),
    
    blood_48 = rep(0,Nit),
    
    Heart_3 = rep(0,Nit),
    
    Kidneys_3 = rep(0,Nit),
    
    Spleen_3 = rep(0,Nit),
    
    Lungs_3 = rep(0,Nit),
    
    Bone_3 = rep(0,Nit),
    
    Heart_48 = rep(0,Nit),
    
    Kidneys_48 = rep(0,Nit),
    
    Spleen_48 = rep(0,Nit),
    
    Lungs_48 = rep(0,Nit),
    
    Bone_48 = rep(0,Nit)
  )
  
  for (i in seq(Nit)){
    
    nt <- nt_Au*10^scales$nt*10^rnorm(1,0,0.2)
    Kbile <- Kbile_Au*10^scales$Kbile*10^rnorm(1,0,0.2)
    Kurine <- 0
    PACt_Kidneys <- 0
    
    PACt_Lungs <- PACt_Lungs_Au*10^scales$PACt_Lungs*10^rnorm(1,0,0.2)
    PACt_Heart <- PACt_Heart_Au*10^scales$PACt_Heart*10^rnorm(1,0,0.2)
    PACt_Bone <- PACt_Bone_Au*10^scales$PACt_Bone*10^rnorm(1,0,0.2)
    PACt_Spleen <- PACt_Spleen_Au*10^scales$PACt_Spleen*10^rnorm(1,0,0.2)
    PACt_Rest <- PACt_Rest_Au*10^scales$PACt_Rest*10^rnorm(1,0,0.2)
    PACt_Brain <- PACt_Rest
    
    Kmax_t_Liver <- 10^scales$Kmax_t_Liver*10^rnorm(1,0,0.2)
    Kmax_t_Kidneys <- 10^scales$Kmax_t_Kidneys*10^rnorm(1,0,0.2)
    Kmax_t_Spleen <- 10^scales$Kmax_t_Spleen*10^rnorm(1,0,0.2)
    Kmax_t_Lungs <- 10^scales$Kmax_t_Lungs*10^rnorm(1,0,0.2)
    
    K50_t_Liver <- K50_t_Liver_Au*10^scales$K50_t_Liver*10^rnorm(1,0,0.2)
    K50_t_Spleen <- K50_t_Spleen_Au*10^scales$K50_t_Spleen*10^rnorm(1,0,0.2)
    K50_t_Kidneys <- K50_t_Kidneys_Au*10^scales$K50_t_Kidneys*10^rnorm(1,0,0.2)
    K50_t_Lungs <- K50_t_Lungs_Au*10^scales$K50_t_Lungs*10^rnorm(1,0,0.2)
    
    Krelease_t_Liver <- Kmax_t_Liver*10^scales$Krelease_t_Liver*10^rnorm(1,0,0.01)
    Krelease_t_Spleen <- Kmax_t_Spleen*10^scales$Krelease_t_Spleen*10^rnorm(1,0,0.01)
    Krelease_t_Kidneys <- Kmax_t_Kidneys*10^scales$Krelease_t_Kidneys*10^rnorm(1,0,0.01)
    Krelease_t_Lungs <- Kmax_t_Lungs*10^scales$Krelease_t_Lungs*10^rnorm(1,0,0.01)
    
    PACt_Liver <- Kbile*10^scales$PACt_Liver*10^rnorm(1,0,0.01)
    
    
    OP <- ode(y=Init_condition,
              times=times,
              func = NP_PBPK2, 
              parms=list(),
              method='lsoda')
    OP <- data.frame(OP)
    OP_rename <- OP
    
    indx_samp <- seq(from = 1, to = Nt, by = 1)
    
    OP_rename$M_Art <- OP[indx_samp,2]
    OP_rename$M_Ven <- OP[indx_samp,3]
    
    OP_rename$M_vasc_Liver <- OP[indx_samp,4]
    OP_rename$M_extra_Liver <- OP[indx_samp,5]
    OP_rename$M_phago_Liver <- OP[indx_samp,6]
    
    OP_rename$M_vasc_Spleen <- OP[indx_samp,7]
    OP_rename$M_extra_Spleen <- OP[indx_samp,8]
    OP_rename$M_phago_Spleen <- OP[indx_samp,9]
    
    OP_rename$M_vasc_Kidneys <- OP[indx_samp,10]
    OP_rename$M_extra_Kidneys <- OP[indx_samp,11]
    OP_rename$M_phago_Kidneys <- OP[indx_samp,12]
    
    OP_rename$M_vasc_Lungs <- OP[indx_samp,13]
    OP_rename$M_extra_Lungs <- OP[indx_samp,14]
    OP_rename$M_phago_Lungs <- OP[indx_samp,15]
    
    OP_rename$M_vasc_Brain <- OP[indx_samp,16]
    OP_rename$M_extra_Brain <- OP[indx_samp,17]
    
    OP_rename$M_vasc_Heart <- OP[indx_samp,18]
    OP_rename$M_extra_Heart <- OP[indx_samp,19]
    
    OP_rename$M_vasc_Bone <- OP[indx_samp,20]
    OP_rename$M_extra_Bone <- OP[indx_samp,21]
    
    OP_rename$M_vasc_Rest <- OP[indx_samp,22]
    OP_rename$M_extra_Rest <- OP[indx_samp,23]
    
    LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
    KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/LV
    LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/LV
    SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/LV
    HV <- OP_rename$M_extra_Heart/LV
    BoneV <- OP_rename$M_extra_Bone/LV
    LV <- 1
    
    blood <- (OP_rename$M_Ven)/(BW*VVenC)
    
    plot(times,blood,type='l',ylim=c(0, 1500))
    points(blood_summary$time, blood_summary$LCI,col='blue')
    points(blood_summary$time, blood_summary$UCI,col='red')
    
    OP_IS$PACt_Liver[i] <- PACt_Liver
    OP_IS$Kbile[i] = Kbile
    
    OP_IS$PACt_Kidneys[i] <- PACt_Kidneys
    OP_IS$Kurine[i] = Kurine
    
    OP_IS$PACt_Bone[i] <- PACt_Bone
    OP_IS$PACt_Spleen[i] <- PACt_Spleen
    OP_IS$PACt_Lungs[i] <- PACt_Lungs
    OP_IS$PACt_Heart[i] <- PACt_Heart
    OP_IS$PACt_Rest[i] <- PACt_Rest
    
    OP_IS$nt[i] = nt
    
    OP_IS$Kmax_t_Liver[i] = Kmax_t_Liver 
    OP_IS$Kmax_t_Spleen[i]  = Kmax_t_Spleen
    OP_IS$Kmax_t_Kidneys[i]  = Kmax_t_Kidneys 
    OP_IS$Kmax_t_Lungs[i]  = Kmax_t_Lungs
    
    OP_IS$Krelease_t_Liver[i]  = Krelease_t_Liver
    OP_IS$Krelease_t_Spleen[i]  = Krelease_t_Spleen
    OP_IS$Krelease_t_Kidneys[i]  = Krelease_t_Kidneys
    OP_IS$Krelease_t_Lungs[i] = Krelease_t_Lungs
    
    OP_IS$K50_t_Liver[i]  = K50_t_Liver
    OP_IS$K50_t_Spleen[i]  = K50_t_Spleen
    OP_IS$K50_t_Kidneys[i]  = K50_t_Kidneys
    OP_IS$K50_t_Lungs[i] = K50_t_Lungs
    
    OP_IS$blood_03[i] = blood[which.min(abs(times-0.03))]
    
    OP_IS$blood_1[i] = blood[which.min(abs(times-1))]
    
    OP_IS$blood_2[i] = blood[which.min(abs(times-2))]
    
    OP_IS$blood_4[i] = blood[which.min(abs(times-4))]
    
    OP_IS$blood_8[i] = blood[which.min(abs(times-8))]
    
    OP_IS$blood_10[i] = blood[which.min(abs(times-10))]
    
    OP_IS$blood_24[i] = blood[which.min(abs(times-24))]
    
    OP_IS$blood_48[i] = blood[which.min(abs(times-48))]
    
    OP_IS$Heart_3[i] = HV[which.min(abs(times-3))]
    
    OP_IS$Kidneys_3[i] = KV[which.min(abs(times-3))]
    
    OP_IS$Spleen_3[i] = SV[which.min(abs(times-3))]
    
    OP_IS$Lungs_3[i] = LuV[which.min(abs(times-3))]
    
    OP_IS$Bone_3[i] = BoneV[which.min(abs(times-3))]
    
    OP_IS$Heart_48[i] = HV[which.min(abs(times-48))]
    
    OP_IS$Kidneys_48[i] = KV[which.min(abs(times-48))]
    
    OP_IS$Spleen_48[i] = SV[which.min(abs(times-48))]
    
    OP_IS$Lungs_48[i] = LuV[which.min(abs(times-48))]
    
    OP_IS$Bone_48[i] = BoneV[which.min(abs(times-48))]
  }
  
  OP_IS$cond.blood_03_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_03_yn[OP_IS$blood_03 < blood_summary$UCI[blood_summary$time==0.03] & 
                           OP_IS$blood_03 > blood_summary$LCI[blood_summary$time==0.03]] <- 1
  
  OP_IS$cond.blood_1_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_1_yn[OP_IS$blood_1 < blood_summary$UCI[blood_summary$time==1] & 
                          OP_IS$blood_1 > blood_summary$LCI[blood_summary$time==1]] <- 1
  
  OP_IS$cond.blood_2_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_2_yn[OP_IS$blood_2 < blood_summary$UCI[blood_summary$time==2] & 
                          OP_IS$blood_2 > blood_summary$LCI[blood_summary$time==2]] <- 1
  
  OP_IS$cond.blood_4_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_4_yn[OP_IS$blood_4 < blood_summary$UCI[blood_summary$time==4] & 
                          OP_IS$blood_4 > blood_summary$LCI[blood_summary$time==4]] <- 1
  
  OP_IS$cond.blood_8_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_8_yn[OP_IS$blood_8 < blood_summary$UCI[blood_summary$time==8] & 
                          OP_IS$blood_8 > blood_summary$LCI[blood_summary$time==8]] <- 1
  
  OP_IS$cond.blood_10_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_10_yn[OP_IS$blood_10 < blood_summary$UCI[blood_summary$time==10] & 
                           OP_IS$blood_10 > blood_summary$LCI[blood_summary$time==10]] <- 1
  
  OP_IS$cond.blood_24_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_24_yn[OP_IS$blood_24< blood_summary$UCI[blood_summary$time==24] & 
                           OP_IS$blood_24 > blood_summary$LCI[blood_summary$time==24]] <- 1
  
  OP_IS$cond.blood_48_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_48_yn[OP_IS$blood_48 < blood_summary$UCI[blood_summary$time==48] & 
                           OP_IS$blood_48 > blood_summary$LCI[blood_summary$time==48]] <- 1
  
  OP_IS$cond.blood.tot <- OP_IS$cond.blood_1_yn + OP_IS$cond.blood_8_yn + OP_IS$cond.blood_24_yn +
    OP_IS$cond.blood_48_yn + OP_IS$cond.blood_03_yn  
  
  OP_IS$cond.Heart_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Heart_3_yn[OP_IS$Heart_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Heart'] & 
                          OP_IS$Heart_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Heart']] <- 1
  
  OP_IS$cond.Kidneys_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Kidneys_3_yn[OP_IS$Kidneys_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Kidneys'] & 
                            OP_IS$Kidneys_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Kidneys']] <- 1
  
  OP_IS$cond.Lungs_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Lungs_3_yn[OP_IS$Lungs_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Lungs'] & 
                          OP_IS$Lungs_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Lungs']] <- 1
  
  OP_IS$cond.Bone_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Bone_3_yn[OP_IS$Bone_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Bone'] & 
                         OP_IS$Bone_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Bone']] <- 1
  
  OP_IS$cond.Spleen_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Spleen_3_yn[OP_IS$Spleen_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Spleen'] & 
                           OP_IS$Spleen_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Spleen']] <- 1
  
  OP_IS$cond.organ.tot_3 <- OP_IS$cond.Spleen_3_yn + OP_IS$cond.Heart_3_yn + 
    OP_IS$cond.Kidneys_3_yn + OP_IS$cond.Lungs_3_yn + OP_IS$cond.Bone_3_yn
  
  OP_IS$cond.Heart_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Heart_48_yn[OP_IS$Heart_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Heart'] & 
                           OP_IS$Heart_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Heart']] <- 1
  
  OP_IS$cond.Kidneys_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Kidneys_48_yn[OP_IS$Kidneys_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Kidneys'] & 
                             OP_IS$Kidneys_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Kidneys']] <- 1
  
  OP_IS$cond.Lungs_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Lungs_48_yn[OP_IS$Lungs_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Lungs'] & 
                           OP_IS$Lungs_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Lungs']] <- 1
  
  OP_IS$cond.Bone_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Bone_48_yn[OP_IS$Bone_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Bone'] & 
                          OP_IS$Bone_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Bone']] <- 1
  
  OP_IS$cond.Spleen_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Spleen_48_yn[OP_IS$Spleen_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Spleen'] & 
                            OP_IS$Spleen_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Spleen']] <- 1
  
  OP_IS$cond.organ.tot_48 <- OP_IS$cond.Spleen_48_yn + OP_IS$cond.Heart_48_yn + 
    OP_IS$cond.Kidneys_48_yn + OP_IS$cond.Lungs_48_yn + OP_IS$cond.Bone_48_yn
  
  OP_IS$PACt_Liver_scale <- log10(OP_IS$PACt_Liver/OP_IS$Kbile)
  OP_IS$Kbile_scale <- log10(OP_IS$Kbile/Kbile_Au)
  
  OP_IS$PACt_Bone_scale <- log10(OP_IS$PACt_Bone/PACt_Bone_Au)
  OP_IS$PACt_Heart_scale <- log10(OP_IS$PACt_Heart/PACt_Heart_Au)
  OP_IS$PACt_Spleen_scale <- log10(OP_IS$PACt_Spleen/PACt_Spleen_Au)
  OP_IS$PACt_Lungs_scale <- log10(OP_IS$PACt_Lungs/PACt_Lungs_Au)
  OP_IS$PACt_Rest_scale <- log10(OP_IS$PACt_Rest/PACt_Rest_Au)
  
  OP_IS$Kmax_t_Liver_scale <- log10(OP_IS$Kmax_t_Liver)
  OP_IS$Kmax_t_Spleen_scale <- log10(OP_IS$Kmax_t_Spleen)
  OP_IS$Kmax_t_Kidneys_scale <- log10(OP_IS$Kmax_t_Kidneys)
  OP_IS$Kmax_t_Lungs_scale <- log10(OP_IS$Kmax_t_Lungs)
  
  OP_IS$Krelease_t_Liver_scale <- log10(OP_IS$Krelease_t_Liver/OP_IS$Kmax_t_Liver)
  OP_IS$Krelease_t_Spleen_scale <- log10(OP_IS$Krelease_t_Spleen/OP_IS$Kmax_t_Spleen)
  OP_IS$Krelease_t_Kidneys_scale <- log10(OP_IS$Krelease_t_Kidneys/OP_IS$Kmax_t_Kidneys)
  OP_IS$Krelease_t_Lungs_scale <- log10(OP_IS$Krelease_t_Lungs/OP_IS$Kmax_t_Lungs)
  
  OP_IS$K50_t_Liver_scale <- log10(OP_IS$K50_t_Liver/K50_t_Liver_Au)
  OP_IS$K50_t_Spleen_scale <- log10(OP_IS$K50_t_Spleen/K50_t_Spleen_Au)
  OP_IS$K50_t_Kidneys_scale <- log10(OP_IS$K50_t_Kidneys/K50_t_Kidneys_Au)
  OP_IS$K50_t_Lungs_scale <- log10(OP_IS$K50_t_Lungs/K50_t_Lungs_Au)
  
  OP_IS$nt_scale <- log10(OP_IS$nt/nt_Au)
  
  OP_IS$dose <- rep(dose_rel[j],Nit)
  
  OP_IS_nana <- OP_IS[complete.cases(OP_IS),]
  
  if (j==1){
    
    OP_IS_tot <- OP_IS_nana
    
  }else{
    
    OP_IS_tot <- rbind(OP_IS_tot,OP_IS_nana)
    
  }
  
}

saveRDS(OP_IS_tot,"model_generated_data/parm2_MC_20240122.rds")

OP_IS_max_cond <- OP_IS_tot[OP_IS_tot$cond.blood.tot%in%c(5),]
OP_IS_max_cond

OP_IS_max_cond_organs <- OP_IS_max_cond[OP_IS_max_cond$cond.organ.tot_3%in%c(3,4,5) &
                                          OP_IS_max_cond$cond.organ.tot_48%in%c(4,5),]
OP_IS_max_cond_organs

OPISmax_Heart <- OP_IS_max_cond[OP_IS_max_cond$cond.Heart_3_yn==1,]# & OP_IS_max_cond$cond.Heart_48_yn, ]
OPISmax_Lungs <- OP_IS_max_cond[OP_IS_max_cond$cond.Lungs_3_yn==1 & OP_IS_max_cond$cond.Lungs_48_yn==1, ]
OPISmax_Kidneys <- OP_IS_max_cond[OP_IS_max_cond$cond.Kidneys_3_yn==1 & OP_IS_max_cond$cond.Kidneys_48_yn==1, ]

scales <- data.frame(
  PACt_Liver=  mean(OP_IS_max_cond_organs$PACt_Liver_scale),
  Kbile= mean(OP_IS_max_cond_organs$Kbile_scale),
  
  PACt_Heart= mean(OPISmax_Heart$PACt_Heart_scale),
  PACt_Lungs= mean(OPISmax_Lungs$PACt_Lungs_scale),
  PACt_Bone= mean(OP_IS_max_cond_organs$PACt_Bone_scale),
  PACt_Spleen= mean(OP_IS_max_cond_organs$PACt_Spleen_scale),
  PACt_Rest= mean(OP_IS_max_cond_organs$PACt_Rest_scale),
  
  Kmax_t_Liver= mean(OP_IS_max_cond_organs$Kmax_t_Liver_scale),
  Kmax_t_Spleen= mean(OP_IS_max_cond_organs$Kmax_t_Spleen_scale),
  Kmax_t_Lungs= mean(OPISmax_Lungs$Kmax_t_Lungs_scale),
  Kmax_t_Kidneys= mean(OPISmax_Kidneys$Kmax_t_Kidneys_scale),
  
  Krelease_t_Liver= mean(OP_IS_max_cond_organs$Krelease_t_Liver_scale),
  Krelease_t_Spleen= mean(OP_IS_max_cond_organs$Krelease_t_Spleen_scale),
  Krelease_t_Lungs= mean(OPISmax_Lungs$Krelease_t_Lungs_scale),
  Krelease_t_Kidneys= mean(OPISmax_Kidneys$Krelease_t_Kidneys_scale),
  
  K50_t_Liver= mean(OP_IS_max_cond_organs$K50_t_Liver_scale),
  K50_t_Spleen= mean(OP_IS_max_cond_organs$K50_t_Spleen_scale),
  K50_t_Lungs= mean(OPISmax_Lungs$K50_t_Lungs_scale),
  K50_t_Kidneys= mean(OPISmax_Kidneys$K50_t_Kidneys_scale),
  
  nt = mean(OP_IS_max_cond_organs$nt_scale)
)

saveRDS(scales,"model_generated_data/scales_20240122.rds")


##########################################
#Parameter analysis - distribution determination

#PACt_Heart - uniform
ff_PACt_Heart <- fitdistr(OP_IS_max_cond_organs$PACt_Heart*1e5, 'gamma')
ff_PACt_Heart_sim <- qgamma(ppoints(OP_IS_max_cond_organs$PACt_Heart), 
                           shape = ff_PACt_Heart$estimate[1],
                           rate = ff_PACt_Heart$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Heart, ff_PACt_Heart_sim)

ff_PACt_Heart <- fitdistr(OP_IS_max_cond_organs$PACt_Heart*1e5, 'exponential')
ff_PACt_Heart_sim <- qexp(ppoints(OP_IS_max_cond_organs$PACt_Heart), 
                         rate = ff_PACt_Heart$estimate)
qqplot(OP_IS_max_cond_organs$PACt_Heart, ff_PACt_Heart_sim)

ff_PACt_Heart <- fitdistr(OP_IS_max_cond_organs$PACt_Heart*1e5, 'normal')
ff_PACt_Heart_sim <- qnorm(ppoints(OP_IS_max_cond_organs$PACt_Heart), 
                           mean = ff_PACt_Heart$estimate[1],
                           sd = ff_PACt_Heart$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Heart, ff_PACt_Heart_sim)

ff_PACt_Heart_sim <- qunif(ppoints(OP_IS_max_cond_organs$PACt_Heart),
                           min(OP_IS_max_cond_organs$PACt_Heart),
                           max(OP_IS_max_cond_organs$PACt_Heart))
qqplot(OP_IS_max_cond_organs$PACt_Heart, ff_PACt_Heart_sim)

#PACt_Bone - gamma
ff_PACt_Bone <- fitdistr(OP_IS_max_cond_organs$PACt_Bone*1e5, 'gamma')
ff_PACt_Bone_sim <- qgamma(ppoints(OP_IS_max_cond_organs$PACt_Bone), 
                            shape = ff_PACt_Bone$estimate[1],
                            rate = ff_PACt_Bone$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Bone, ff_PACt_Bone_sim)

ff_PACt_Bone <- fitdistr(OP_IS_max_cond_organs$PACt_Bone*1e5, 'exponential')
ff_PACt_Bone_sim <- qexp(ppoints(OP_IS_max_cond_organs$PACt_Bone), 
                          rate = ff_PACt_Bone$estimate)
qqplot(OP_IS_max_cond_organs$PACt_Bone, ff_PACt_Bone_sim)

ff_PACt_Bone <- fitdistr(OP_IS_max_cond_organs$PACt_Bone*1e5, 'normal')
ff_PACt_Bone_sim <- qnorm(ppoints(OP_IS_max_cond_organs$PACt_Bone), 
                           mean = ff_PACt_Bone$estimate[1],
                           sd = ff_PACt_Bone$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Bone, ff_PACt_Bone_sim)

ff_PACt_Bone_sim <- qunif(ppoints(OP_IS_max_cond_organs$PACt_Bone),
                           min(OP_IS_max_cond_organs$PACt_Bone),
                           max(OP_IS_max_cond_organs$PACt_Bone))
qqplot(OP_IS_max_cond_organs$PACt_Bone, ff_PACt_Bone_sim)

#PACt_Spleen - gamma
ff_PACt_Spleen <- fitdistr(OP_IS_max_cond_organs$PACt_Spleen*1e5, 'gamma')
ff_PACt_Spleen_sim <- qgamma(ppoints(OP_IS_max_cond_organs$PACt_Spleen), 
                            shape = ff_PACt_Spleen$estimate[1],
                            rate = ff_PACt_Spleen$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Spleen, ff_PACt_Spleen_sim)

ff_PACt_Spleen <- fitdistr(OP_IS_max_cond_organs$PACt_Spleen*1e5, 'exponential')
ff_PACt_Spleen_sim <- qexp(ppoints(OP_IS_max_cond_organs$PACt_Spleen), 
                          rate = ff_PACt_Spleen$estimate)
qqplot(OP_IS_max_cond_organs$PACt_Spleen, ff_PACt_Spleen_sim)

ff_PACt_Spleen <- fitdistr(OP_IS_max_cond_organs$PACt_Spleen*1e5, 'normal')
ff_PACt_Spleen_sim <- qnorm(ppoints(OP_IS_max_cond_organs$PACt_Spleen), 
                           mean = ff_PACt_Spleen$estimate[1],
                           sd = ff_PACt_Spleen$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Spleen, ff_PACt_Spleen_sim)

ff_PACt_Spleen_sim <- qunif(ppoints(OP_IS_max_cond_organs$PACt_Spleen),
                           min(OP_IS_max_cond_organs$PACt_Spleen),
                           max(OP_IS_max_cond_organs$PACt_Spleen))
qqplot(OP_IS_max_cond_organs$PACt_Spleen, ff_PACt_Spleen_sim)



#PACt_Lungs - gamma
ff_PACt_Lungs <- fitdistr(OP_IS_max_cond_organs$PACt_Lungs*1e5, 'gamma')
ff_PACt_Lungs_sim <- qgamma(ppoints(OP_IS_max_cond_organs$PACt_Lungs), 
                           shape = ff_PACt_Lungs$estimate[1],
                           rate = ff_PACt_Lungs$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Lungs, ff_PACt_Lungs_sim)

ff_PACt_Lungs <- fitdistr(OP_IS_max_cond_organs$PACt_Lungs*1e5, 'exponential')
ff_PACt_Lungs_sim <- qexp(ppoints(OP_IS_max_cond_organs$PACt_Lungs), 
                         rate = ff_PACt_Lungs$estimate)
qqplot(OP_IS_max_cond_organs$PACt_Lungs, ff_PACt_Lungs_sim)

ff_PACt_Lungs <- fitdistr(OP_IS_max_cond_organs$PACt_Lungs*1e5, 'normal')
ff_PACt_Lungs_sim <- qnorm(ppoints(OP_IS_max_cond_organs$PACt_Lungs), 
                          mean = ff_PACt_Lungs$estimate[1],
                          sd = ff_PACt_Lungs$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Lungs, ff_PACt_Lungs_sim)

ff_PACt_Lungs_sim <- qunif(ppoints(OP_IS_max_cond_organs$PACt_Lungs),
                          min(OP_IS_max_cond_organs$PACt_Lungs),
                          max(OP_IS_max_cond_organs$PACt_Lungs))
qqplot(OP_IS_max_cond_organs$PACt_Lungs, ff_PACt_Lungs_sim)

#PACt_Rest - gamma
ff_PACt_Rest <- fitdistr(OP_IS_max_cond_organs$PACt_Rest*1e5, 'gamma')
ff_PACt_Rest_sim <- qgamma(ppoints(OP_IS_max_cond_organs$PACt_Rest), 
                            shape = ff_PACt_Rest$estimate[1],
                            rate = ff_PACt_Rest$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Rest, ff_PACt_Rest_sim)

ff_PACt_Rest <- fitdistr(OP_IS_max_cond_organs$PACt_Rest*1e5, 'exponential')
ff_PACt_Rest_sim <- qexp(ppoints(OP_IS_max_cond_organs$PACt_Rest), 
                          rate = ff_PACt_Rest$estimate)
qqplot(OP_IS_max_cond_organs$PACt_Rest, ff_PACt_Rest_sim)

ff_PACt_Rest <- fitdistr(OP_IS_max_cond_organs$PACt_Rest*1e5, 'normal')
ff_PACt_Rest_sim <- qnorm(ppoints(OP_IS_max_cond_organs$PACt_Rest), 
                           mean = ff_PACt_Rest$estimate[1],
                           sd = ff_PACt_Rest$estimate[2])
qqplot(OP_IS_max_cond_organs$PACt_Rest, ff_PACt_Rest_sim)

ff_PACt_Rest_sim <- qunif(ppoints(OP_IS_max_cond_organs$PACt_Rest),
                           min(OP_IS_max_cond_organs$PACt_Rest),
                           max(OP_IS_max_cond_organs$PACt_Rest))
qqplot(OP_IS_max_cond_organs$PACt_Rest, ff_PACt_Rest_sim)

#Kbile - uniform
ff_Kbile <- fitdistr(OP_IS_max_cond_organs$Kbile*1e5, 'gamma')
ff_Kbile_sim <- qgamma(ppoints(OP_IS_max_cond_organs$Kbile), 
                           shape = ff_Kbile$estimate[1],
                           rate = ff_Kbile$estimate[2])
qqplot(OP_IS_max_cond_organs$Kbile, ff_Kbile_sim)

ff_Kbile <- fitdistr(OP_IS_max_cond_organs$Kbile*1e5, 'exponential')
ff_Kbile_sim <- qexp(ppoints(OP_IS_max_cond_organs$Kbile), 
                         rate = ff_Kbile$estimate)
qqplot(OP_IS_max_cond_organs$Kbile, ff_Kbile_sim)

ff_Kbile <- fitdistr(OP_IS_max_cond_organs$Kbile*1e5, 'normal')
ff_Kbile_sim <- qnorm(ppoints(OP_IS_max_cond_organs$Kbile), 
                          mean = ff_Kbile$estimate[1],
                          sd = ff_Kbile$estimate[2])
qqplot(OP_IS_max_cond_organs$Kbile, ff_Kbile_sim)

ff_Kbile_sim <- qunif(ppoints(OP_IS_max_cond_organs$Kbile),
                          min(OP_IS_max_cond_organs$Kbile),
                          max(OP_IS_max_cond_organs$Kbile))
qqplot(OP_IS_max_cond_organs$Kbile, ff_Kbile_sim)

#nt - gamma
ff_nt <- fitdistr(OP_IS_max_cond_organs$nt, 'gamma')
ff_nt_sim <- qgamma(ppoints(OP_IS_max_cond_organs$nt), 
                       shape = ff_nt$estimate[1],
                       rate = ff_nt$estimate[2])
qqplot(OP_IS_max_cond_organs$nt, ff_nt_sim)

ff_nt <- fitdistr(OP_IS_max_cond_organs$nt, 'exponential')
ff_nt_sim <- qexp(ppoints(OP_IS_max_cond_organs$nt), 
                     rate = ff_nt$estimate)
qqplot(OP_IS_max_cond_organs$nt, ff_nt_sim)

ff_nt <- fitdistr(OP_IS_max_cond_organs$nt, 'normal')
ff_nt_sim <- qnorm(ppoints(OP_IS_max_cond_organs$nt), 
                      mean = ff_nt$estimate[1],
                      sd = ff_nt$estimate[2])
qqplot(OP_IS_max_cond_organs$nt, ff_nt_sim)

ff_nt_sim <- qunif(ppoints(OP_IS_max_cond_organs$nt),
                      min(OP_IS_max_cond_organs$nt),
                      max(OP_IS_max_cond_organs$nt))
qqplot(OP_IS_max_cond_organs$nt, ff_nt_sim)

#Kmax_t_Liver - gamma
ff_Kmax_t_Liver <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Liver, 'gamma')
ff_Kmax_t_Liver_sim <- qgamma(ppoints(OP_IS_max_cond_organs$Kmax_t_Liver), 
                    shape = ff_Kmax_t_Liver$estimate[1],
                    rate = ff_Kmax_t_Liver$estimate[2])
qqplot(OP_IS_max_cond_organs$Kmax_t_Liver, ff_Kmax_t_Liver_sim)

ff_Kmax_t_Liver <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Liver, 'exponential')
ff_Kmax_t_Liver_sim <- qexp(ppoints(OP_IS_max_cond_organs$Kmax_t_Liver), 
                  rate = ff_Kmax_t_Liver$estimate)
qqplot(OP_IS_max_cond_organs$Kmax_t_Liver, ff_Kmax_t_Liver_sim)

ff_Kmax_t_Liver <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Liver, 'normal')
ff_Kmax_t_Liver_sim <- qnorm(ppoints(OP_IS_max_cond_organs$Kmax_t_Liver), 
                   mean = ff_Kmax_t_Liver$estimate[1],
                   sd = ff_Kmax_t_Liver$estimate[2])
qqplot(OP_IS_max_cond_organs$Kmax_t_Liver, ff_Kmax_t_Liver_sim)

ff_Kmax_t_Liver_sim <- qunif(ppoints(OP_IS_max_cond_organs$Kmax_t_Liver),
                   min(OP_IS_max_cond_organs$Kmax_t_Liver),
                   max(OP_IS_max_cond_organs$Kmax_t_Liver))
qqplot(OP_IS_max_cond_organs$Kmax_t_Liver, ff_Kmax_t_Liver_sim)

#Kmax_t_Lungs - gamma
ff_Kmax_t_Lungs <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Lungs, 'gamma')
ff_Kmax_t_Lungs_sim <- qgamma(ppoints(OP_IS_max_cond_organs$Kmax_t_Lungs), 
                              shape = ff_Kmax_t_Lungs$estimate[1],
                              rate = ff_Kmax_t_Lungs$estimate[2])
qqplot(OP_IS_max_cond_organs$Kmax_t_Lungs, ff_Kmax_t_Lungs_sim)

ff_Kmax_t_Lungs <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Lungs, 'exponential')
ff_Kmax_t_Lungs_sim <- qexp(ppoints(OP_IS_max_cond_organs$Kmax_t_Lungs), 
                            rate = ff_Kmax_t_Lungs$estimate)
qqplot(OP_IS_max_cond_organs$Kmax_t_Lungs, ff_Kmax_t_Lungs_sim)

ff_Kmax_t_Lungs <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Lungs, 'normal')
ff_Kmax_t_Lungs_sim <- qnorm(ppoints(OP_IS_max_cond_organs$Kmax_t_Lungs), 
                             mean = ff_Kmax_t_Lungs$estimate[1],
                             sd = ff_Kmax_t_Lungs$estimate[2])
qqplot(OP_IS_max_cond_organs$Kmax_t_Lungs, ff_Kmax_t_Lungs_sim)

ff_Kmax_t_Lungs_sim <- qunif(ppoints(OP_IS_max_cond_organs$Kmax_t_Lungs),
                             min(OP_IS_max_cond_organs$Kmax_t_Lungs),
                             max(OP_IS_max_cond_organs$Kmax_t_Lungs))
qqplot(OP_IS_max_cond_organs$Kmax_t_Lungs, ff_Kmax_t_Lungs_sim)

#Kmax_t_Kidneys - gamma
ff_Kmax_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Kidneys, 'gamma')
ff_Kmax_t_Kidneys_sim <- qgamma(ppoints(OP_IS_max_cond_organs$Kmax_t_Kidneys), 
                              shape = ff_Kmax_t_Kidneys$estimate[1],
                              rate = ff_Kmax_t_Kidneys$estimate[2])
qqplot(OP_IS_max_cond_organs$Kmax_t_Kidneys, ff_Kmax_t_Kidneys_sim)

ff_Kmax_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Kidneys, 'exponential')
ff_Kmax_t_Kidneys_sim <- qexp(ppoints(OP_IS_max_cond_organs$Kmax_t_Kidneys), 
                            rate = ff_Kmax_t_Kidneys$estimate)
qqplot(OP_IS_max_cond_organs$Kmax_t_Kidneys, ff_Kmax_t_Kidneys_sim)

ff_Kmax_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Kidneys, 'normal')
ff_Kmax_t_Kidneys_sim <- qnorm(ppoints(OP_IS_max_cond_organs$Kmax_t_Kidneys), 
                             mean = ff_Kmax_t_Kidneys$estimate[1],
                             sd = ff_Kmax_t_Kidneys$estimate[2])
qqplot(OP_IS_max_cond_organs$Kmax_t_Kidneys, ff_Kmax_t_Kidneys_sim)

ff_Kmax_t_Kidneys_sim <- qunif(ppoints(OP_IS_max_cond_organs$Kmax_t_Kidneys),
                             min(OP_IS_max_cond_organs$Kmax_t_Kidneys),
                             max(OP_IS_max_cond_organs$Kmax_t_Kidneys))
qqplot(OP_IS_max_cond_organs$Kmax_t_Kidneys, ff_Kmax_t_Kidneys_sim)

#Kmax_t_Spleen - gamma
ff_Kmax_t_Spleen <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Spleen, 'gamma')
ff_Kmax_t_Spleen_sim <- qgamma(ppoints(OP_IS_max_cond_organs$Kmax_t_Spleen), 
                                shape = ff_Kmax_t_Spleen$estimate[1],
                                rate = ff_Kmax_t_Spleen$estimate[2])
qqplot(OP_IS_max_cond_organs$Kmax_t_Spleen, ff_Kmax_t_Spleen_sim)

ff_Kmax_t_Spleen <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Spleen, 'exponential')
ff_Kmax_t_Spleen_sim <- qexp(ppoints(OP_IS_max_cond_organs$Kmax_t_Spleen), 
                              rate = ff_Kmax_t_Spleen$estimate)
qqplot(OP_IS_max_cond_organs$Kmax_t_Spleen, ff_Kmax_t_Spleen_sim)

ff_Kmax_t_Spleen <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Spleen, 'normal')
ff_Kmax_t_Spleen_sim <- qnorm(ppoints(OP_IS_max_cond_organs$Kmax_t_Spleen), 
                               mean = ff_Kmax_t_Spleen$estimate[1],
                               sd = ff_Kmax_t_Spleen$estimate[2])
qqplot(OP_IS_max_cond_organs$Kmax_t_Spleen, ff_Kmax_t_Spleen_sim)

ff_Kmax_t_Spleen_sim <- qunif(ppoints(OP_IS_max_cond_organs$Kmax_t_Spleen),
                               min(OP_IS_max_cond_organs$Kmax_t_Spleen),
                               max(OP_IS_max_cond_organs$Kmax_t_Spleen))
qqplot(OP_IS_max_cond_organs$Kmax_t_Spleen, ff_Kmax_t_Spleen_sim)

#K50_t_Liver - gamma
ff_K50_t_Liver <- fitdistr(OP_IS_max_cond_organs$K50_t_Liver, 'gamma')
ff_K50_t_Liver_sim <- qgamma(ppoints(OP_IS_max_cond_organs$K50_t_Liver), 
                              shape = ff_K50_t_Liver$estimate[1],
                              rate = ff_K50_t_Liver$estimate[2])
qqplot(OP_IS_max_cond_organs$K50_t_Liver, ff_K50_t_Liver_sim)

ff_K50_t_Liver <- fitdistr(OP_IS_max_cond_organs$K50_t_Liver, 'exponential')
ff_K50_t_Liver_sim <- qexp(ppoints(OP_IS_max_cond_organs$K50_t_Liver), 
                            rate = ff_K50_t_Liver$estimate)
qqplot(OP_IS_max_cond_organs$K50_t_Liver, ff_K50_t_Liver_sim)

ff_K50_t_Liver <- fitdistr(OP_IS_max_cond_organs$K50_t_Liver, 'normal')
ff_K50_t_Liver_sim <- qnorm(ppoints(OP_IS_max_cond_organs$K50_t_Liver), 
                             mean = ff_K50_t_Liver$estimate[1],
                             sd = ff_K50_t_Liver$estimate[2])
qqplot(OP_IS_max_cond_organs$K50_t_Liver, ff_K50_t_Liver_sim)

ff_K50_t_Liver_sim <- qunif(ppoints(OP_IS_max_cond_organs$K50_t_Liver),
                             min(OP_IS_max_cond_organs$K50_t_Liver),
                             max(OP_IS_max_cond_organs$K50_t_Liver))
qqplot(OP_IS_max_cond_organs$K50_t_Liver, ff_K50_t_Liver_sim)

#K50_t_Lungs - gamma
ff_K50_t_Lungs <- fitdistr(OP_IS_max_cond_organs$K50_t_Lungs, 'gamma')
ff_K50_t_Lungs_sim <- qgamma(ppoints(OP_IS_max_cond_organs$K50_t_Lungs), 
                              shape = ff_K50_t_Lungs$estimate[1],
                              rate = ff_K50_t_Lungs$estimate[2])
qqplot(OP_IS_max_cond_organs$K50_t_Lungs, ff_K50_t_Lungs_sim)

ff_K50_t_Lungs <- fitdistr(OP_IS_max_cond_organs$K50_t_Lungs, 'exponential')
ff_K50_t_Lungs_sim <- qexp(ppoints(OP_IS_max_cond_organs$K50_t_Lungs), 
                            rate = ff_K50_t_Lungs$estimate)
qqplot(OP_IS_max_cond_organs$K50_t_Lungs, ff_K50_t_Lungs_sim)

ff_K50_t_Lungs <- fitdistr(OP_IS_max_cond_organs$K50_t_Lungs, 'normal')
ff_K50_t_Lungs_sim <- qnorm(ppoints(OP_IS_max_cond_organs$K50_t_Lungs), 
                             mean = ff_K50_t_Lungs$estimate[1],
                             sd = ff_K50_t_Lungs$estimate[2])
qqplot(OP_IS_max_cond_organs$K50_t_Lungs, ff_K50_t_Lungs_sim)

ff_K50_t_Lungs_sim <- qunif(ppoints(OP_IS_max_cond_organs$K50_t_Lungs),
                             min(OP_IS_max_cond_organs$K50_t_Lungs),
                             max(OP_IS_max_cond_organs$K50_t_Lungs))
qqplot(OP_IS_max_cond_organs$K50_t_Lungs, ff_K50_t_Lungs_sim)

#K50_t_Kidneys - gamma
ff_K50_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$K50_t_Kidneys, 'gamma')
ff_K50_t_Kidneys_sim <- qgamma(ppoints(OP_IS_max_cond_organs$K50_t_Kidneys), 
                                shape = ff_K50_t_Kidneys$estimate[1],
                                rate = ff_K50_t_Kidneys$estimate[2])
qqplot(OP_IS_max_cond_organs$K50_t_Kidneys, ff_K50_t_Kidneys_sim)

ff_K50_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$K50_t_Kidneys, 'exponential')
ff_K50_t_Kidneys_sim <- qexp(ppoints(OP_IS_max_cond_organs$K50_t_Kidneys), 
                              rate = ff_K50_t_Kidneys$estimate)
qqplot(OP_IS_max_cond_organs$K50_t_Kidneys, ff_K50_t_Kidneys_sim)

ff_K50_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$K50_t_Kidneys, 'normal')
ff_K50_t_Kidneys_sim <- qnorm(ppoints(OP_IS_max_cond_organs$K50_t_Kidneys), 
                               mean = ff_K50_t_Kidneys$estimate[1],
                               sd = ff_K50_t_Kidneys$estimate[2])
qqplot(OP_IS_max_cond_organs$K50_t_Kidneys, ff_K50_t_Kidneys_sim)

ff_K50_t_Kidneys_sim <- qunif(ppoints(OP_IS_max_cond_organs$K50_t_Kidneys),
                               min(OP_IS_max_cond_organs$K50_t_Kidneys),
                               max(OP_IS_max_cond_organs$K50_t_Kidneys))
qqplot(OP_IS_max_cond_organs$K50_t_Kidneys, ff_K50_t_Kidneys_sim)

#K50_t_Spleen - gamma
ff_K50_t_Spleen <- fitdistr(OP_IS_max_cond_organs$K50_t_Spleen, 'gamma')
ff_K50_t_Spleen_sim <- qgamma(ppoints(OP_IS_max_cond_organs$K50_t_Spleen), 
                               shape = ff_K50_t_Spleen$estimate[1],
                               rate = ff_K50_t_Spleen$estimate[2])
qqplot(OP_IS_max_cond_organs$K50_t_Spleen, ff_K50_t_Spleen_sim)

ff_K50_t_Spleen <- fitdistr(OP_IS_max_cond_organs$K50_t_Spleen, 'exponential')
ff_K50_t_Spleen_sim <- qexp(ppoints(OP_IS_max_cond_organs$K50_t_Spleen), 
                             rate = ff_K50_t_Spleen$estimate)
qqplot(OP_IS_max_cond_organs$K50_t_Spleen, ff_K50_t_Spleen_sim)

ff_K50_t_Spleen <- fitdistr(OP_IS_max_cond_organs$K50_t_Spleen, 'normal')
ff_K50_t_Spleen_sim <- qnorm(ppoints(OP_IS_max_cond_organs$K50_t_Spleen), 
                              mean = ff_K50_t_Spleen$estimate[1],
                              sd = ff_K50_t_Spleen$estimate[2])
qqplot(OP_IS_max_cond_organs$K50_t_Spleen, ff_K50_t_Spleen_sim)

ff_K50_t_Spleen_sim <- qunif(ppoints(OP_IS_max_cond_organs$K50_t_Spleen),
                              min(OP_IS_max_cond_organs$K50_t_Spleen),
                              max(OP_IS_max_cond_organs$K50_t_Spleen))
qqplot(OP_IS_max_cond_organs$K50_t_Spleen, ff_K50_t_Spleen_sim)

###################################################
#Fifth section: generate simulation runs at 0.5mg

#Clear all variables
rm(list=ls())

set.seed(126)

#Source scripts
source("libs.R")
source("parms.R")
alpha_temp_FACS <- 1e-2
alpha_temp_blood <- 3e-2
source("funcs.R")

OP_IS_tot <- readRDS("model_generated_data/parm2_MC_20240122.rds")
scales <- readRDS("model_generated_data/scales_20240122.rds")

OP_IS_max_cond <- OP_IS_tot[OP_IS_tot$cond.blood.tot%in%c(5),]
OP_IS_max_cond

OP_IS_max_cond_organs <- OP_IS_max_cond[OP_IS_max_cond$cond.organ.tot_3%in%c(3,4,5) &
                                   OP_IS_max_cond$cond.organ.tot_48%in%c(4,5),]
OPISmax_Heart <- OP_IS_max_cond[OP_IS_max_cond$cond.Heart_3_yn==1,]# & OP_IS_max_cond$cond.Heart_48_yn, ]
OPISmax_Lungs <- OP_IS_max_cond[OP_IS_max_cond$cond.Lungs_3_yn==1 & OP_IS_max_cond$cond.Lungs_48_yn==1, ]
OPISmax_Kidneys <- OP_IS_max_cond[OP_IS_max_cond$cond.Kidneys_3_yn==1 & OP_IS_max_cond$cond.Kidneys_48_yn==1, ]

scales <- data.frame(
  PACt_Liver=  mean(OP_IS_max_cond_organs$PACt_Liver_scale),
  Kbile= mean(OP_IS_max_cond_organs$Kbile_scale),
  
  PACt_Heart= mean(OPISmax_Heart$PACt_Heart_scale),
  PACt_Lungs= mean(OP_IS_max_cond_organs$PACt_Lungs_scale),
  PACt_Bone= mean(OP_IS_max_cond_organs$PACt_Bone_scale),
  PACt_Spleen= mean(OP_IS_max_cond_organs$PACt_Spleen_scale),
  PACt_Rest= mean(OP_IS_max_cond_organs$PACt_Rest_scale),
  
  Kmax_t_Liver= mean(OP_IS_max_cond_organs$Kmax_t_Liver_scale),
  Kmax_t_Spleen= mean(OP_IS_max_cond_organs$Kmax_t_Spleen_scale),
  Kmax_t_Lungs= mean(OP_IS_max_cond_organs$Kmax_t_Lungs_scale),
  Kmax_t_Kidneys= mean(OP_IS_max_cond_organs$Kmax_t_Kidneys_scale),
  
  Krelease_t_Liver= mean(OP_IS_max_cond_organs$Krelease_t_Liver_scale),
  Krelease_t_Spleen= mean(OP_IS_max_cond_organs$Krelease_t_Spleen_scale),
  Krelease_t_Lungs= mean(OP_IS_max_cond_organs$Krelease_t_Lungs_scale),
  Krelease_t_Kidneys= mean(OP_IS_max_cond_organs$Krelease_t_Kidneys_scale),
  
  K50_t_Liver= mean(OP_IS_max_cond_organs$K50_t_Liver_scale),
  K50_t_Spleen= mean(OP_IS_max_cond_organs$K50_t_Spleen_scale),
  K50_t_Lungs= mean(OP_IS_max_cond_organs$K50_t_Lungs_scale),
  K50_t_Kidneys= mean(OP_IS_max_cond_organs$K50_t_Kidneys_scale),
  
  nt = mean(OP_IS_max_cond_organs$nt_scale)
)

dist_setmod <- 1e5

ff_PACt_Bone <- fitdistr(OP_IS_max_cond_organs$PACt_Bone*dist_setmod, 'gamma')
ff_PACt_Spleen <- fitdistr(OP_IS_max_cond_organs$PACt_Spleen*dist_setmod, 'gamma')
ff_PACt_Lungs <- fitdistr(OPISmax_Lungs$PACt_Lungs*dist_setmod, 'gamma')
ff_PACt_Rest <- fitdistr(OP_IS_max_cond_organs$PACt_Rest*dist_setmod, 'gamma')

ff_Kbile <- fitdistr(OP_IS_max_cond_organs$Kbile*dist_setmod, 'gamma')
ff_nt <- fitdistr(OP_IS_max_cond_organs$nt, 'gamma')

ff_Kmax_t_Liver <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Liver, 'gamma')
ff_Kmax_t_Lungs <- fitdistr(OPISmax_Lungs$Kmax_t_Lungs, 'gamma')
ff_Kmax_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Kidneys, 'gamma')
ff_Kmax_t_Spleen <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Spleen, 'gamma')

ff_K50_t_Liver <- fitdistr(OP_IS_max_cond_organs$K50_t_Liver, 'gamma')
ff_K50_t_Lungs <- fitdistr(OP_IS_max_cond_organs$K50_t_Lungs, 'gamma')
ff_K50_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$K50_t_Kidneys, 'gamma')
ff_K50_t_Spleen <- fitdistr(OP_IS_max_cond_organs$K50_t_Spleen, 'gamma')

OP_list <- list()

Nit <- 1e2

for (j in 1:length(dose_rel)){
  
  source("FACS_IVIS_data_analysis_20230426.R")
  
  blood_conc <- blood_conc[blood_conc$dose==dose_rel[j],]
  blood_summary <- blood_summary[blood_summary$dose==dose_rel[j],]
  organ_FACS_time48_summary <- organ_FACS_summary[organ_FACS_summary$Dose==dose_rel[j],]
  organ_FACS_time24_summary <- organ_FACS_time_summary[organ_FACS_time_summary$time==24,]
  organ_FACS_time3_summary <- organ_FACS_time_summary[organ_FACS_time_summary$time==3,]
  
  dose_rate <- dose_rate_vec[j]
  
  Init_condition <- rep(0,22)
  
  OP_IS <- data.frame(
    dose = rep(dose_rel[j],Nit),
    
    nt = rep(0,Nit),
    Kbile = rep(0,Nit),
    Kurine = rep(0,Nit),
    
    PACt_Liver = rep(0,Nit),
    PACt_Kidneys = rep(0,Nit),
    PACt_Spleen = rep(0,Nit),
    PACt_Lungs = rep(0,Nit),
    PACt_Heart = rep(0,Nit),
    PACt_Bone = rep(0,Nit),
    PACt_Rest = rep(0,Nit),
    
    Kmax_t_Liver = rep(0,Nit),
    Kmax_t_Spleen  = rep(0,Nit),
    Kmax_t_Kidneys  = rep(0,Nit),
    Kmax_t_Lungs  = rep(0,Nit),
    
    Krelease_t_Liver  = rep(0,Nit),
    Krelease_t_Spleen  = rep(0,Nit),
    Krelease_t_Kidneys  = rep(0,Nit),
    Krelease_t_Lungs = rep(0,Nit),
    
    K50_t_Liver  = rep(0,Nit),
    K50_t_Spleen  = rep(0,Nit),
    K50_t_Kidneys  = rep(0,Nit),
    K50_t_Lungs = rep(0,Nit),
    
    blood_03 = rep(0,Nit),
    
    blood_1 = rep(0,Nit),
    
    blood_2 = rep(0,Nit),
    
    blood_4 = rep(0,Nit),
    
    blood_8 = rep(0,Nit),
    
    blood_10 = rep(0,Nit),
    
    blood_24 = rep(0,Nit),
    
    blood_48 = rep(0,Nit),
    
    Heart_3 = rep(0,Nit),
    
    Kidneys_3 = rep(0,Nit),
    
    Spleen_3 = rep(0,Nit),
    
    Lungs_3 = rep(0,Nit),
    
    Bone_3 = rep(0,Nit),
    
    Heart_48 = rep(0,Nit),
    
    Kidneys_48 = rep(0,Nit),
    
    Spleen_48 = rep(0,Nit),
    
    Lungs_48 = rep(0,Nit),
    
    Bone_48 = rep(0,Nit)
  )
  
  for (i in seq(Nit)){
    
    nt <- rgamma(1,shape=ff_nt$estimate[1],
                 rate=ff_nt$estimate[2])
    Kbile <- runif(1,min(OP_IS_max_cond_organs$Kbile),
                        max(OP_IS_max_cond_organs$Kbile))
    
    PACt_Lungs <- rgamma(1,shape=ff_PACt_Lungs$estimate[1],
                          rate=ff_PACt_Lungs$estimate[2])/dist_setmod*PACt_Lungs_Mod
    PACt_Heart <- runif(1,min(OP_IS_max_cond_organs$PACt_Heart),
                        max(OP_IS_max_cond_organs$PACt_Heart))*PACt_Heart_Mod
    PACt_Bone <- rgamma(1,shape=ff_PACt_Bone$estimate[1],
                        rate=ff_PACt_Bone$estimate[2])/dist_setmod
    PACt_Spleen <- rgamma(1,shape=ff_PACt_Spleen$estimate[1],
                          rate=ff_PACt_Spleen$estimate[2])/dist_setmod*PACt_Spleen_Mod
    PACt_Rest <- rgamma(1,shape=ff_PACt_Rest$estimate[1],
                        rate=ff_PACt_Rest$estimate[2])/dist_setmod
    
    Kmax_t_Liver <- rgamma(1,shape=ff_Kmax_t_Liver$estimate[1],
                           rate=ff_Kmax_t_Liver$estimate[2])
    Kmax_t_Kidneys <- rgamma(1,shape=ff_Kmax_t_Kidneys$estimate[1],
                             rate=ff_Kmax_t_Kidneys$estimate[2])*Kmax_t_Kidneys_Mod
    Kmax_t_Lungs <- rgamma(1,shape=ff_Kmax_t_Lungs$estimate[1],
                           rate=ff_Kmax_t_Lungs$estimate[2])*Kmax_t_Lungs_Mod
    Kmax_t_Spleen <- rgamma(1,shape=ff_Kmax_t_Spleen$estimate[1],
                            rate=ff_Kmax_t_Spleen$estimate[2])
    
    K50_t_Liver <- rgamma(1,shape=ff_K50_t_Liver$estimate[1],
                          rate=ff_K50_t_Liver$estimate[2])
    K50_t_Spleen <- rgamma(1,shape=ff_K50_t_Spleen$estimate[1],
                           rate=ff_K50_t_Spleen$estimate[2])
    K50_t_Kidneys <- rgamma(1,shape=ff_K50_t_Kidneys$estimate[1],
                            rate=ff_K50_t_Kidneys$estimate[2])
    K50_t_Lungs <- rgamma(1,shape=ff_K50_t_Lungs$estimate[1],
                          rate=ff_K50_t_Lungs$estimate[2])
    
    PACt_Liver <- Kbile*10^scales$PACt_Liver*10^rnorm(1,0,0.01)
    
    Krelease_t_Liver <- Kmax_t_Liver*10^scales$Krelease_t_Liver*10^rnorm(1,0,0.01)
    Krelease_t_Spleen <- Kmax_t_Spleen*10^scales$Krelease_t_Spleen*10^rnorm(1,0,0.01)
    Krelease_t_Kidneys <- Kmax_t_Kidneys*10^scales$Krelease_t_Kidneys*10^rnorm(1,0,0.01)
    Krelease_t_Lungs <- Kmax_t_Lungs*10^scales$Krelease_t_Lungs*10^rnorm(1,0,0.01)
    
    Kurine <- 0
    PACt_Kidneys <- PACt_Kidneys_Mod
    PACt_Brain <- PACt_Rest
    
    
    OP <- ode(y=Init_condition,
              times=times,
              func = NP_PBPK2, 
              parms=list(),
              method='lsoda')
    OP <- data.frame(OP)
    OP_rename <- OP
    
    indx_samp <- seq(from = 1, to = Nt, by = 1)
    
    OP_rename$M_Art <- OP[indx_samp,2]
    OP_rename$M_Ven <- OP[indx_samp,3]
    
    OP_rename$M_vasc_Liver <- OP[indx_samp,4]
    OP_rename$M_extra_Liver <- OP[indx_samp,5]
    OP_rename$M_phago_Liver <- OP[indx_samp,6]
    
    OP_rename$M_vasc_Spleen <- OP[indx_samp,7]
    OP_rename$M_extra_Spleen <- OP[indx_samp,8]
    OP_rename$M_phago_Spleen <- OP[indx_samp,9]
    
    OP_rename$M_vasc_Kidneys <- OP[indx_samp,10]
    OP_rename$M_extra_Kidneys <- OP[indx_samp,11]
    OP_rename$M_phago_Kidneys <- OP[indx_samp,12]
    
    OP_rename$M_vasc_Lungs <- OP[indx_samp,13]
    OP_rename$M_extra_Lungs <- OP[indx_samp,14]
    OP_rename$M_phago_Lungs <- OP[indx_samp,15]
    
    OP_rename$M_vasc_Brain <- OP[indx_samp,16]
    OP_rename$M_extra_Brain <- OP[indx_samp,17]
    
    OP_rename$M_vasc_Heart <- OP[indx_samp,18]
    OP_rename$M_extra_Heart <- OP[indx_samp,19]
    
    OP_rename$M_vasc_Bone <- OP[indx_samp,20]
    OP_rename$M_extra_Bone <- OP[indx_samp,21]
    
    OP_rename$M_vasc_Rest <- OP[indx_samp,22]
    OP_rename$M_extra_Rest <- OP[indx_samp,23]
    
    LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
    KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/LV
    LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/LV
    SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/LV
    HV <- OP_rename$M_extra_Heart/LV
    BoneV <- OP_rename$M_extra_Bone/LV
    LV <- 1
    
    blood <- (OP_rename$M_Ven)/(BW*VVenC)
    
    OP_rename$LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
    OP_rename$KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/OP_rename$LV
    OP_rename$LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/OP_rename$LV
    OP_rename$SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/OP_rename$LV
    OP_rename$HV <- OP_rename$M_extra_Heart/OP_rename$LV
    OP_rename$BoneV <- OP_rename$M_extra_Bone/OP_rename$LV
    
    OP_rename$Blood <- (OP_rename$M_Ven+OP_rename$M_Art)/(BW*VBloodC)
    OP_rename$Liver <- (OP_rename$M_extra_Liver+OP_rename$M_phago_Liver)/(BW*VLC)
    OP_rename$Spleen <- (OP_rename$M_extra_Spleen+OP_rename$M_phago_Spleen)/(BW*VSC)
    OP_rename$Kidneys <- (OP_rename$M_extra_Kidneys+OP_rename$M_phago_Kidneys)/(BW*VKC)
    OP_rename$Lungs <- (OP_rename$M_extra_Lungs+OP_rename$M_phago_Lungs)/(BW*VLuC)
    OP_rename$Heart <- (OP_rename$M_extra_Heart)/(BW*VHC)
    OP_rename$Bone <- (OP_rename$M_extra_Bone)/(BW*VBoneC)  
    OP_rename$Rest <- (OP_rename$M_extra_Rest)/(BW*VrestC) 
    OP_rename$Brain <- (OP_rename$M_extra_Brain)/(BW*VBRC) 
    
    
    OP_list[[i]] <- OP_rename
    
    plot(times,blood,type='l',ylim=c(0, 1500))
    points(blood_summary$time, blood_summary$LCI,col='blue')
    points(blood_summary$time, blood_summary$UCI,col='red')
    
    OP_IS$PACt_Liver[i] <- PACt_Liver
    OP_IS$Kbile[i] = Kbile
    
    OP_IS$PACt_Kidneys[i] <- PACt_Kidneys
    OP_IS$Kurine[i] = Kurine
    
    OP_IS$PACt_Bone[i] <- PACt_Bone
    OP_IS$PACt_Spleen[i] <- PACt_Spleen
    OP_IS$PACt_Lungs[i] <- PACt_Lungs
    OP_IS$PACt_Heart[i] <- PACt_Heart
    OP_IS$PACt_Rest[i] <- PACt_Rest
    
    OP_IS$nt[i] = nt
    
    OP_IS$Kmax_t_Liver[i] = Kmax_t_Liver 
    OP_IS$Kmax_t_Spleen[i]  = Kmax_t_Spleen
    OP_IS$Kmax_t_Kidneys[i]  = Kmax_t_Kidneys 
    OP_IS$Kmax_t_Lungs[i]  = Kmax_t_Lungs
    
    OP_IS$Krelease_t_Liver[i]  = Krelease_t_Liver
    OP_IS$Krelease_t_Spleen[i]  = Krelease_t_Spleen
    OP_IS$Krelease_t_Kidneys[i]  = Krelease_t_Kidneys
    OP_IS$Krelease_t_Lungs[i] = Krelease_t_Lungs
    
    OP_IS$K50_t_Liver[i]  = K50_t_Liver
    OP_IS$K50_t_Spleen[i]  = K50_t_Spleen
    OP_IS$K50_t_Kidneys[i]  = K50_t_Kidneys
    OP_IS$K50_t_Lungs[i] = K50_t_Lungs
    
    OP_IS$blood_03[i] = blood[which.min(abs(times-0.03))]
    
    OP_IS$blood_1[i] = blood[which.min(abs(times-1))]
    
    OP_IS$blood_2[i] = blood[which.min(abs(times-2))]
    
    OP_IS$blood_4[i] = blood[which.min(abs(times-4))]
    
    OP_IS$blood_8[i] = blood[which.min(abs(times-8))]
    
    OP_IS$blood_10[i] = blood[which.min(abs(times-10))]
    
    OP_IS$blood_24[i] = blood[which.min(abs(times-24))]
    
    OP_IS$blood_48[i] = blood[which.min(abs(times-48))]
    
    OP_IS$Heart_3[i] = HV[which.min(abs(times-3))]
    
    OP_IS$Kidneys_3[i] = KV[which.min(abs(times-3))]
    
    OP_IS$Spleen_3[i] = SV[which.min(abs(times-3))]
    
    OP_IS$Lungs_3[i] = LuV[which.min(abs(times-3))]
    
    OP_IS$Bone_3[i] = BoneV[which.min(abs(times-3))]
    
    OP_IS$Heart_48[i] = HV[which.min(abs(times-48))]
    
    OP_IS$Kidneys_48[i] = KV[which.min(abs(times-48))]
    
    OP_IS$Spleen_48[i] = SV[which.min(abs(times-48))]
    
    OP_IS$Lungs_48[i] = LuV[which.min(abs(times-48))]
    
    OP_IS$Bone_48[i] = BoneV[which.min(abs(times-48))]
  }
  
  OP_IS$cond.blood_03_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_03_yn[OP_IS$blood_03 < blood_summary$UCI[blood_summary$time==0.03] & 
                           OP_IS$blood_03 > blood_summary$LCI[blood_summary$time==0.03]] <- 1
  
  OP_IS$cond.blood_1_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_1_yn[OP_IS$blood_1 < blood_summary$UCI[blood_summary$time==1] & 
                          OP_IS$blood_1 > blood_summary$LCI[blood_summary$time==1]] <- 1
  
  OP_IS$cond.blood_2_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_2_yn[OP_IS$blood_2 < blood_summary$UCI[blood_summary$time==2] & 
                          OP_IS$blood_2 > blood_summary$LCI[blood_summary$time==2]] <- 1
  
  OP_IS$cond.blood_4_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_4_yn[OP_IS$blood_4 < blood_summary$UCI[blood_summary$time==4] & 
                          OP_IS$blood_4 > blood_summary$LCI[blood_summary$time==4]] <- 1
  
  OP_IS$cond.blood_8_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_8_yn[OP_IS$blood_8 < blood_summary$UCI[blood_summary$time==8] & 
                          OP_IS$blood_8 > blood_summary$LCI[blood_summary$time==8]] <- 1
  
  OP_IS$cond.blood_10_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_10_yn[OP_IS$blood_10 < blood_summary$UCI[blood_summary$time==10] & 
                           OP_IS$blood_10 > blood_summary$LCI[blood_summary$time==10]] <- 1
  
  OP_IS$cond.blood_48_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_48_yn[OP_IS$blood_48< blood_summary$UCI[blood_summary$time==48] & 
                           OP_IS$blood_48 > blood_summary$LCI[blood_summary$time==48]] <- 1
  
  OP_IS$cond.blood_48_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_48_yn[OP_IS$blood_48 < blood_summary$UCI[blood_summary$time==48] & 
                           OP_IS$blood_48 > blood_summary$LCI[blood_summary$time==48]] <- 1
  
  OP_IS$cond.blood.tot <- OP_IS$cond.blood_1_yn + OP_IS$cond.blood_8_yn + 
    OP_IS$cond.blood_48_yn + OP_IS$cond.blood_03_yn  
  
  OP_IS$cond.Heart_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Heart_3_yn[OP_IS$Heart_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Heart'] & 
                          OP_IS$Heart_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Heart']] <- 1
  
  OP_IS$cond.Kidneys_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Kidneys_3_yn[OP_IS$Kidneys_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Kidneys'] & 
                            OP_IS$Kidneys_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Kidneys']] <- 1
  
  OP_IS$cond.Lungs_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Lungs_3_yn[OP_IS$Lungs_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Lungs'] & 
                          OP_IS$Lungs_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Lungs']] <- 1
  
  OP_IS$cond.Bone_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Bone_3_yn[OP_IS$Bone_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Bone'] & 
                         OP_IS$Bone_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Bone']] <- 1
  
  OP_IS$cond.Spleen_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Spleen_3_yn[OP_IS$Spleen_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Spleen'] & 
                           OP_IS$Spleen_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Spleen']] <- 1
  
  OP_IS$cond.organ.tot_3 <- OP_IS$cond.Spleen_3_yn + OP_IS$cond.Heart_3_yn + 
    OP_IS$cond.Kidneys_3_yn + OP_IS$cond.Lungs_3_yn + OP_IS$cond.Bone_3_yn
  
  OP_IS$cond.Heart_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Heart_48_yn[OP_IS$Heart_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Heart'] & 
                           OP_IS$Heart_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Heart']] <- 1
  
  OP_IS$cond.Kidneys_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Kidneys_48_yn[OP_IS$Kidneys_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Kidneys'] & 
                             OP_IS$Kidneys_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Kidneys']] <- 1
  
  OP_IS$cond.Lungs_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Lungs_48_yn[OP_IS$Lungs_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Lungs'] & 
                           OP_IS$Lungs_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Lungs']] <- 1
  
  OP_IS$cond.Bone_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Bone_48_yn[OP_IS$Bone_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Bone'] & 
                          OP_IS$Bone_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Bone']] <- 1
  
  OP_IS$cond.Spleen_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Spleen_48_yn[OP_IS$Spleen_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Spleen'] & 
                            OP_IS$Spleen_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Spleen']] <- 1
  
  OP_IS$cond.organ.tot_48 <- OP_IS$cond.Spleen_48_yn + OP_IS$cond.Heart_48_yn + 
    OP_IS$cond.Kidneys_48_yn + OP_IS$cond.Lungs_48_yn + OP_IS$cond.Bone_48_yn
  
  OP_IS$PACt_Liver_scale <- log10(OP_IS$PACt_Liver/OP_IS$Kbile)
  OP_IS$Kbile_scale <- log10(OP_IS$Kbile/Kbile_Au)
  
  OP_IS$PACt_Bone_scale <- log10(OP_IS$PACt_Bone/PACt_Bone_Au)
  OP_IS$PACt_Heart_scale <- log10(OP_IS$PACt_Heart/PACt_Heart_Au)
  OP_IS$PACt_Spleen_scale <- log10(OP_IS$PACt_Spleen/PACt_Spleen_Au)
  OP_IS$PACt_Lungs_scale <- log10(OP_IS$PACt_Lungs/PACt_Lungs_Au)
  OP_IS$PACt_Rest_scale <- log10(OP_IS$PACt_Rest/PACt_Rest_Au)
  
  OP_IS$Kmax_t_Liver_scale <- log10(OP_IS$Kmax_t_Liver/Kmax_t_Liver_Au)
  OP_IS$Kmax_t_Spleen_scale <- log10(OP_IS$Kmax_t_Spleen/Kmax_t_Spleen_Au)
  OP_IS$Kmax_t_Kidneys_scale <- log10(OP_IS$Kmax_t_Kidneys/Kmax_t_Kidneys_Au)
  OP_IS$Kmax_t_Lungs_scale <- log10(OP_IS$Kmax_t_Lungs/Kmax_t_Lungs_Au)
  
  OP_IS$Krelease_t_Liver_scale <- log10(OP_IS$Krelease_t_Liver/OP_IS$Kmax_t_Liver)
  OP_IS$Krelease_t_Spleen_scale <- log10(OP_IS$Krelease_t_Spleen/OP_IS$Kmax_t_Spleen)
  OP_IS$Krelease_t_Kidneys_scale <- log10(OP_IS$Krelease_t_Kidneys/OP_IS$Kmax_t_Kidneys)
  OP_IS$Krelease_t_Lungs_scale <- log10(OP_IS$Krelease_t_Lungs/OP_IS$Kmax_t_Lungs)
  
  OP_IS$K50_t_Liver_scale <- log10(OP_IS$K50_t_Liver/K50_t_Liver_Au)
  OP_IS$K50_t_Spleen_scale <- log10(OP_IS$K50_t_Spleen/K50_t_Spleen_Au)
  OP_IS$K50_t_Kidneys_scale <- log10(OP_IS$K50_t_Kidneys/K50_t_Kidneys_Au)
  OP_IS$K50_t_Lungs_scale <- log10(OP_IS$K50_t_Lungs/K50_t_Lungs_Au)
  
  OP_IS$nt_scale <- log10(OP_IS$nt/nt_Au)
  
  OP_IS$dose <- rep(dose_rel[j],Nit)
  
  OP_IS_nana <- OP_IS[complete.cases(OP_IS),]
  
  if (j==1){
    
    OP_IS_tot <- OP_IS_nana
    
  }else{
    
    OP_IS_tot <- rbind(OP_IS_tot,OP_IS_nana)
    
  }
  
}

Blood <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Blood)})),nrow=Nt)

Bone <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Bone)})),nrow=Nt)

Kidneys <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Kidneys)})),nrow=Nt)

Liver <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Liver)})),nrow=Nt)

Spleen <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Spleen)})),nrow=Nt)

Lungs <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Lungs)})),nrow=Nt)

Heart <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Heart)})),nrow=Nt)

Brain <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Brain)})),nrow=Nt)

Rest <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Rest)})),nrow=Nt)

OP_mn_sd <- data.frame(Blood_mn = rep(0,Nt),
                       Blood_sd = rep(0,Nt),
                       Heart_mn = rep(0,Nt),
                       Heart_sd = rep(0,Nt),
                       Liver_mn = rep(0,Nt),
                       Liver_sd = rep(0,Nt),
                       Bone_mn = rep(0,Nt),
                       Bone_sd = rep(0,Nt),
                       Kidneys_mn = rep(0,Nt),
                       Kidneys_sd = rep(0,Nt),
                       Spleen_mn = rep(0,Nt),
                       Spleen_sd = rep(0,Nt),
                       Lungs_mn = rep(0,Nt),
                       Lungs_sd = rep(0,Nt),
                       Rest_mn = rep(0,Nt),
                       Rest_sd = rep(0,Nt),
                       Brain_mn = rep(0,Nt),
                       Brain_sd = rep(0,Nt))

OP_mn_sd$Blood_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Blood[ind,]))
})
OP_mn_sd$Blood_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Blood[ind,]))
})

OP_mn_sd$Heart_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Heart[ind,]))
})
OP_mn_sd$Heart_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Heart[ind,]))
})

OP_mn_sd$Lungs_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Lungs[ind,]))
})
OP_mn_sd$Lungs_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Lungs[ind,]))
})

OP_mn_sd$Bone_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Bone[ind,]))
})
OP_mn_sd$Bone_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Bone[ind,]))
})

OP_mn_sd$Spleen_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Spleen[ind,]))
})
OP_mn_sd$Spleen_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Spleen[ind,]))
})

OP_mn_sd$Kidneys_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Kidneys[ind,]))
})
OP_mn_sd$Kidneys_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Kidneys[ind,]))
})

OP_mn_sd$Liver_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Liver[ind,]))
})
OP_mn_sd$Liver_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Liver[ind,]))
})

OP_mn_sd$Rest_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Rest[ind,]))
})
OP_mn_sd$Rest_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Rest[ind,]))
})
OP_mn_sd$Brain_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Brain[ind,]))
})
OP_mn_sd$Brain_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Brain[ind,]))
})

indx_temp <- times==3

OP_organ_sum3 <- data.frame(Heart_mn=mean(Heart[indx_temp,]*VHC/Liver[indx_temp,]/VLC),
                            Heart_sd=sd(Heart[indx_temp,]*VHC/Liver[indx_temp,]/VLC),
                            Lungs_mn=mean(Lungs[indx_temp,]*VLuC/Liver[indx_temp,]/VLC),
                            Lungs_sd=sd(Lungs[indx_temp,]*VLuC/Liver[indx_temp,]/VLC),
                            Kidneys_mn=mean(Kidneys[indx_temp,]*VKC/Liver[indx_temp,]/VLC),
                            Kidneys_sd=sd(Kidneys[indx_temp,]*VKC/Liver[indx_temp,]/VLC),
                            Spleen_mn=mean(Spleen[indx_temp,]*VSC/Liver[indx_temp,]/VLC),
                            Spleen_sd=sd(Spleen[indx_temp,]*VSC/Liver[indx_temp,]/VLC),
                            Bone_mn=mean(Bone[indx_temp,]*VBoneC/Liver[indx_temp,]/VLC),
                            Bone_sd=sd(Bone[indx_temp,]*VBoneC/Liver[indx_temp,]/VLC))

indx_temp <- times==6

OP_organ_sum6 <- data.frame(Heart_mn=mean(Heart[indx_temp,]*VHC/Liver[indx_temp,]/VLC),
                            Heart_sd=sd(Heart[indx_temp,]*VHC/Liver[indx_temp,]/VLC),
                            Lungs_mn=mean(Lungs[indx_temp,]*VLuC/Liver[indx_temp,]/VLC),
                            Lungs_sd=sd(Lungs[indx_temp,]*VLuC/Liver[indx_temp,]/VLC),
                            Kidneys_mn=mean(Kidneys[indx_temp,]*VKC/Liver[indx_temp,]/VLC),
                            Kidneys_sd=sd(Kidneys[indx_temp,]*VKC/Liver[indx_temp,]/VLC),
                            Spleen_mn=mean(Spleen[indx_temp,]*VSC/Liver[indx_temp,]/VLC),
                            Spleen_sd=sd(Spleen[indx_temp,]*VSC/Liver[indx_temp,]/VLC),
                            Bone_mn=mean(Bone[indx_temp,]*VBoneC/Liver[indx_temp,]/VLC),
                            Bone_sd=sd(Bone[indx_temp,]*VBoneC/Liver[indx_temp,]/VLC))

indx_temp <- times==24

OP_organ_sum24 <- data.frame(Heart_mn=mean(Heart[indx_temp,]*VHC/Liver[indx_temp,]/VLC),
                             Heart_sd=sd(Heart[indx_temp,]*VHC/Liver[indx_temp,]/VLC),
                             Lungs_mn=mean(Lungs[indx_temp,]*VLuC/Liver[indx_temp,]/VLC),
                             Lungs_sd=sd(Lungs[indx_temp,]*VLuC/Liver[indx_temp,]/VLC),
                             Kidneys_mn=mean(Kidneys[indx_temp,]*VKC/Liver[indx_temp,]/VLC),
                             Kidneys_sd=sd(Kidneys[indx_temp,]*VKC/Liver[indx_temp,]/VLC),
                             Spleen_mn=mean(Spleen[indx_temp,]*VSC/Liver[indx_temp,]/VLC),
                             Spleen_sd=sd(Spleen[indx_temp,]*VSC/Liver[indx_temp,]/VLC),
                             Bone_mn=mean(Bone[indx_temp,]*VBoneC/Liver[indx_temp,]/VLC),
                             Bone_sd=sd(Bone[indx_temp,]*VBoneC/Liver[indx_temp,]/VLC))
indx_temp <- times==48

OP_organ_sum48 <- data.frame(Heart_mn=mean(Heart[indx_temp,]*VHC/Liver[indx_temp,]/VLC),
                             Heart_sd=sd(Heart[indx_temp,]*VHC/Liver[indx_temp,]/VLC),
                             Lungs_mn=mean(Lungs[indx_temp,]*VLuC/Liver[indx_temp,]/VLC),
                             Lungs_sd=sd(Lungs[indx_temp,]*VLuC/Liver[indx_temp,]/VLC),
                             Kidneys_mn=mean(Kidneys[indx_temp,]*VKC/Liver[indx_temp,]/VLC),
                             Kidneys_sd=sd(Kidneys[indx_temp,]*VKC/Liver[indx_temp,]/VLC),
                             Spleen_mn=mean(Spleen[indx_temp,]*VSC/Liver[indx_temp,]/VLC),
                             Spleen_sd=sd(Spleen[indx_temp,]*VSC/Liver[indx_temp,]/VLC),
                             Bone_mn=mean(Bone[indx_temp,]*VBoneC/Liver[indx_temp,]/VLC),
                             Bone_sd=sd(Bone[indx_temp,]*VBoneC/Liver[indx_temp,]/VLC))



writeMat("model_generated_data/sim_runs_20240128.mat",OP=OP_mn_sd)
writeMat("model_generated_data/sim_runs_20240128_organ_sum.mat",
         OP_organ_sum3=OP_organ_sum3, OP_organ_sum6=OP_organ_sum6, 
         OP_organ_sum24=OP_organ_sum24,OP_organ_sum48=OP_organ_sum48)

###################################################
#Sixth section: generate simulation runs for all doses

#Clear all variables
rm(list=ls())

set.seed(87)

#Source scripts
source("libs.R")
source("parms.R")
alpha_temp_FACS <- 1e-2
alpha_temp_blood <- 3e-2
source("funcs.R")

OP_IS_tot <- readRDS("model_generated_data/parm2_MC_20240122.rds")
scales <- readRDS("model_generated_data/scales_20240122.rds")

OP_IS_max_cond <- OP_IS_tot[OP_IS_tot$cond.blood.tot%in%c(5),]
OP_IS_max_cond

OP_IS_max_cond_organs <- OP_IS_max_cond[OP_IS_max_cond$cond.organ.tot_3%in%c(3,4,5) &
                                          OP_IS_max_cond$cond.organ.tot_48%in%c(4,5),]
OPISmax_Heart <- OP_IS_max_cond[OP_IS_max_cond$cond.Heart_3_yn==1,]# & OP_IS_max_cond$cond.Heart_48_yn, ]
OPISmax_Lungs <- OP_IS_max_cond[OP_IS_max_cond$cond.Lungs_3_yn==1 & OP_IS_max_cond$cond.Lungs_48_yn==1, ]
OPISmax_Kidneys <- OP_IS_max_cond[OP_IS_max_cond$cond.Kidneys_3_yn==1 & OP_IS_max_cond$cond.Kidneys_48_yn==1, ]

scales <- data.frame(
  PACt_Liver=  mean(OP_IS_max_cond_organs$PACt_Liver_scale),
  Kbile= mean(OP_IS_max_cond_organs$Kbile_scale),
  
  PACt_Heart= mean(OPISmax_Heart$PACt_Heart_scale),
  PACt_Lungs= mean(OP_IS_max_cond_organs$PACt_Lungs_scale),
  PACt_Bone= mean(OP_IS_max_cond_organs$PACt_Bone_scale),
  PACt_Spleen= mean(OP_IS_max_cond_organs$PACt_Spleen_scale),
  PACt_Rest= mean(OP_IS_max_cond_organs$PACt_Rest_scale),
  
  Kmax_t_Liver= mean(OP_IS_max_cond_organs$Kmax_t_Liver_scale),
  Kmax_t_Spleen= mean(OP_IS_max_cond_organs$Kmax_t_Spleen_scale),
  Kmax_t_Lungs= mean(OP_IS_max_cond_organs$Kmax_t_Lungs_scale),
  Kmax_t_Kidneys= mean(OP_IS_max_cond_organs$Kmax_t_Kidneys_scale),
  
  Krelease_t_Liver= mean(OP_IS_max_cond_organs$Krelease_t_Liver_scale),
  Krelease_t_Spleen= mean(OP_IS_max_cond_organs$Krelease_t_Spleen_scale),
  Krelease_t_Lungs= mean(OP_IS_max_cond_organs$Krelease_t_Lungs_scale),
  Krelease_t_Kidneys= mean(OP_IS_max_cond_organs$Krelease_t_Kidneys_scale),
  
  K50_t_Liver= mean(OP_IS_max_cond_organs$K50_t_Liver_scale),
  K50_t_Spleen= mean(OP_IS_max_cond_organs$K50_t_Spleen_scale),
  K50_t_Lungs= mean(OP_IS_max_cond_organs$K50_t_Lungs_scale),
  K50_t_Kidneys= mean(OP_IS_max_cond_organs$K50_t_Kidneys_scale),
  
  nt = mean(OP_IS_max_cond_organs$nt_scale)
)

dist_setmod <- 1e5

ff_PACt_Bone <- fitdistr(OP_IS_max_cond_organs$PACt_Bone*dist_setmod, 'gamma')
ff_PACt_Spleen <- fitdistr(OP_IS_max_cond_organs$PACt_Spleen*dist_setmod, 'gamma')
ff_PACt_Lungs <- fitdistr(OPISmax_Lungs$PACt_Lungs*dist_setmod, 'gamma')
ff_PACt_Rest <- fitdistr(OP_IS_max_cond_organs$PACt_Rest*dist_setmod, 'gamma')

ff_Kbile <- fitdistr(OP_IS_max_cond_organs$Kbile*dist_setmod, 'gamma')
ff_nt <- fitdistr(OP_IS_max_cond_organs$nt, 'gamma')

ff_Kmax_t_Liver <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Liver, 'gamma')
ff_Kmax_t_Lungs <- fitdistr(OPISmax_Lungs$Kmax_t_Lungs, 'gamma')
ff_Kmax_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Kidneys, 'gamma')
ff_Kmax_t_Spleen <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Spleen, 'gamma')

ff_K50_t_Liver <- fitdistr(OP_IS_max_cond_organs$K50_t_Liver, 'gamma')
ff_K50_t_Lungs <- fitdistr(OP_IS_max_cond_organs$K50_t_Lungs, 'gamma')
ff_K50_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$K50_t_Kidneys, 'gamma')
ff_K50_t_Spleen <- fitdistr(OP_IS_max_cond_organs$K50_t_Spleen, 'gamma')


OP_list <- list()

Nit <- 1e2

for (j in 1:length(dose_rel)){
  
  source("FACS_IVIS_data_analysis_20230426.R")
  
  blood_conc <- blood_conc[blood_conc$dose==dose_rel[j],]
  blood_summary <- blood_summary[blood_summary$dose==dose_rel[j],]
  organ_FACS_time48_summary <- organ_FACS_summary[organ_FACS_summary$Dose==dose_rel[j],]
  organ_FACS_time24_summary <- organ_FACS_time_summary[organ_FACS_time_summary$time==24,]
  organ_FACS_time3_summary <- organ_FACS_time_summary[organ_FACS_time_summary$time==3,]
  
  dose_rate <- dose_rate_vec[j]
  
  Init_condition <- rep(0,22)
  
  OP_IS <- data.frame(
    dose = rep(dose_rel[j],Nit),
    
    nt = rep(0,Nit),
    Kbile = rep(0,Nit),
    Kurine = rep(0,Nit),
    
    PACt_Liver = rep(0,Nit),
    PACt_Kidneys = rep(0,Nit),
    PACt_Spleen = rep(0,Nit),
    PACt_Lungs = rep(0,Nit),
    PACt_Heart = rep(0,Nit),
    PACt_Bone = rep(0,Nit),
    PACt_Rest = rep(0,Nit),
    
    Kmax_t_Liver = rep(0,Nit),
    Kmax_t_Spleen  = rep(0,Nit),
    Kmax_t_Kidneys  = rep(0,Nit),
    Kmax_t_Lungs  = rep(0,Nit),
    
    Krelease_t_Liver  = rep(0,Nit),
    Krelease_t_Spleen  = rep(0,Nit),
    Krelease_t_Kidneys  = rep(0,Nit),
    Krelease_t_Lungs = rep(0,Nit),
    
    K50_t_Liver  = rep(0,Nit),
    K50_t_Spleen  = rep(0,Nit),
    K50_t_Kidneys  = rep(0,Nit),
    K50_t_Lungs = rep(0,Nit),
    
    blood_03 = rep(0,Nit),
    
    blood_1 = rep(0,Nit),
    
    blood_2 = rep(0,Nit),
    
    blood_4 = rep(0,Nit),
    
    blood_8 = rep(0,Nit),
    
    blood_10 = rep(0,Nit),
    
    blood_24 = rep(0,Nit),
    
    blood_48 = rep(0,Nit),
    
    Heart_3 = rep(0,Nit),
    
    Kidneys_3 = rep(0,Nit),
    
    Spleen_3 = rep(0,Nit),
    
    Lungs_3 = rep(0,Nit),
    
    Bone_3 = rep(0,Nit),
    
    Heart_48 = rep(0,Nit),
    
    Kidneys_48 = rep(0,Nit),
    
    Spleen_48 = rep(0,Nit),
    
    Lungs_48 = rep(0,Nit),
    
    Bone_48 = rep(0,Nit)
  )
  
  for (i in seq(Nit)){
    
    nt <- rgamma(1,shape=ff_nt$estimate[1],
                 rate=ff_nt$estimate[2])
    Kbile <- runif(1,min(OP_IS_max_cond_organs$Kbile),
                   max(OP_IS_max_cond_organs$Kbile))
    
    PACt_Lungs <- rgamma(1,shape=ff_PACt_Lungs$estimate[1],
                         rate=ff_PACt_Lungs$estimate[2])/dist_setmod*PACt_Lungs_Mod
    PACt_Heart <- runif(1,min(OP_IS_max_cond_organs$PACt_Heart),
                        max(OP_IS_max_cond_organs$PACt_Heart))*PACt_Heart_Mod
    PACt_Bone <- rgamma(1,shape=ff_PACt_Bone$estimate[1],
                        rate=ff_PACt_Bone$estimate[2])/dist_setmod
    PACt_Spleen <- rgamma(1,shape=ff_PACt_Spleen$estimate[1],
                          rate=ff_PACt_Spleen$estimate[2])/dist_setmod
    PACt_Rest <- rgamma(1,shape=ff_PACt_Rest$estimate[1],
                        rate=ff_PACt_Rest$estimate[2])/dist_setmod
    
    Kmax_t_Liver <- rgamma(1,shape=ff_Kmax_t_Liver$estimate[1],
                           rate=ff_Kmax_t_Liver$estimate[2])
    Kmax_t_Kidneys <- rgamma(1,shape=ff_Kmax_t_Kidneys$estimate[1],
                             rate=ff_Kmax_t_Kidneys$estimate[2])*Kmax_t_Kidneys_Mod
    Kmax_t_Lungs <- rgamma(1,shape=ff_Kmax_t_Lungs$estimate[1],
                           rate=ff_Kmax_t_Lungs$estimate[2])*Kmax_t_Lungs_Mod
    Kmax_t_Spleen <- rgamma(1,shape=ff_Kmax_t_Spleen$estimate[1],
                            rate=ff_Kmax_t_Spleen$estimate[2])
    
    K50_t_Liver <- rgamma(1,shape=ff_K50_t_Liver$estimate[1],
                          rate=ff_K50_t_Liver$estimate[2])
    K50_t_Spleen <- rgamma(1,shape=ff_K50_t_Spleen$estimate[1],
                           rate=ff_K50_t_Spleen$estimate[2])
    K50_t_Kidneys <- rgamma(1,shape=ff_K50_t_Kidneys$estimate[1],
                            rate=ff_K50_t_Kidneys$estimate[2])
    K50_t_Lungs <- rgamma(1,shape=ff_K50_t_Lungs$estimate[1],
                          rate=ff_K50_t_Lungs$estimate[2])
    
    PACt_Liver <- Kbile*10^scales$PACt_Liver*10^rnorm(1,0,0.01)
    
    Krelease_t_Liver <- Kmax_t_Liver*10^scales$Krelease_t_Liver*10^rnorm(1,0,0.01)
    Krelease_t_Spleen <- Kmax_t_Spleen*10^scales$Krelease_t_Spleen*10^rnorm(1,0,0.01)
    Krelease_t_Kidneys <- Kmax_t_Kidneys*10^scales$Krelease_t_Kidneys*10^rnorm(1,0,0.01)
    Krelease_t_Lungs <- Kmax_t_Lungs*10^scales$Krelease_t_Lungs*10^rnorm(1,0,0.01)
    
    Kurine <- 0
    PACt_Kidneys <- PACt_Kidneys_Mod
    PACt_Brain <- PACt_Rest
    
    
    OP <- ode(y=Init_condition,
              times=times,
              func = NP_PBPK2, 
              parms=list(),
              method='lsoda')
    OP <- data.frame(OP)
    OP_rename <- OP
    
    indx_samp <- seq(from = 1, to = Nt, by = 1)
    
    OP_rename$M_Art <- OP[indx_samp,2]
    OP_rename$M_Ven <- OP[indx_samp,3]
    
    OP_rename$M_vasc_Liver <- OP[indx_samp,4]
    OP_rename$M_extra_Liver <- OP[indx_samp,5]
    OP_rename$M_phago_Liver <- OP[indx_samp,6]
    
    OP_rename$M_vasc_Spleen <- OP[indx_samp,7]
    OP_rename$M_extra_Spleen <- OP[indx_samp,8]
    OP_rename$M_phago_Spleen <- OP[indx_samp,9]
    
    OP_rename$M_vasc_Kidneys <- OP[indx_samp,10]
    OP_rename$M_extra_Kidneys <- OP[indx_samp,11]
    OP_rename$M_phago_Kidneys <- OP[indx_samp,12]
    
    OP_rename$M_vasc_Lungs <- OP[indx_samp,13]
    OP_rename$M_extra_Lungs <- OP[indx_samp,14]
    OP_rename$M_phago_Lungs <- OP[indx_samp,15]
    
    OP_rename$M_vasc_Brain <- OP[indx_samp,16]
    OP_rename$M_extra_Brain <- OP[indx_samp,17]
    
    OP_rename$M_vasc_Heart <- OP[indx_samp,18]
    OP_rename$M_extra_Heart <- OP[indx_samp,19]
    
    OP_rename$M_vasc_Bone <- OP[indx_samp,20]
    OP_rename$M_extra_Bone <- OP[indx_samp,21]
    
    OP_rename$M_vasc_Rest <- OP[indx_samp,22]
    OP_rename$M_extra_Rest <- OP[indx_samp,23]
    
    LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
    KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/LV
    LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/LV
    SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/LV
    HV <- OP_rename$M_extra_Heart/LV
    BoneV <- OP_rename$M_extra_Bone/LV
    LV <- 1
    
    blood <- (OP_rename$M_Ven)/(BW*VVenC)
    
    OP_rename$LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
    OP_rename$KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/OP_rename$LV
    OP_rename$LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/OP_rename$LV
    OP_rename$SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/OP_rename$LV
    OP_rename$HV <- OP_rename$M_extra_Heart/OP_rename$LV
    OP_rename$BoneV <- OP_rename$M_extra_Bone/OP_rename$LV
    
    OP_rename$Blood <- (OP_rename$M_Ven+OP_rename$M_Art)/(BW*VBloodC)
    OP_rename$Liver <- (OP_rename$M_extra_Liver+OP_rename$M_phago_Liver)/(BW*VLC)
    OP_rename$Spleen <- (OP_rename$M_extra_Spleen+OP_rename$M_phago_Spleen)/(BW*VSC)
    OP_rename$Kidneys <- (OP_rename$M_extra_Kidneys+OP_rename$M_phago_Kidneys)/(BW*VKC)
    OP_rename$Lungs <- (OP_rename$M_extra_Lungs+OP_rename$M_phago_Lungs)/(BW*VLuC)
    OP_rename$Heart <- (OP_rename$M_extra_Heart)/(BW*VHC)
    OP_rename$Bone <- (OP_rename$M_extra_Bone)/(BW*VBoneC)  
    OP_rename$Rest <- (OP_rename$M_extra_Rest)/(BW*VrestC) 
    OP_rename$Brain <- (OP_rename$M_extra_Brain)/(BW*VBRC) 
    
    
    OP_list[[i]] <- OP_rename
    
    plot(times,blood,type='l',ylim=c(0, 1500))
    points(blood_summary$time, blood_summary$LCI,col='blue')
    points(blood_summary$time, blood_summary$UCI,col='red')
    
    OP_IS$PACt_Liver[i] <- PACt_Liver
    OP_IS$Kbile[i] = Kbile
    
    OP_IS$PACt_Kidneys[i] <- PACt_Kidneys
    OP_IS$Kurine[i] = Kurine
    
    OP_IS$PACt_Bone[i] <- PACt_Bone
    OP_IS$PACt_Spleen[i] <- PACt_Spleen
    OP_IS$PACt_Lungs[i] <- PACt_Lungs
    OP_IS$PACt_Heart[i] <- PACt_Heart
    OP_IS$PACt_Rest[i] <- PACt_Rest
    
    OP_IS$nt[i] = nt
    
    OP_IS$Kmax_t_Liver[i] = Kmax_t_Liver 
    OP_IS$Kmax_t_Spleen[i]  = Kmax_t_Spleen
    OP_IS$Kmax_t_Kidneys[i]  = Kmax_t_Kidneys 
    OP_IS$Kmax_t_Lungs[i]  = Kmax_t_Lungs
    
    OP_IS$Krelease_t_Liver[i]  = Krelease_t_Liver
    OP_IS$Krelease_t_Spleen[i]  = Krelease_t_Spleen
    OP_IS$Krelease_t_Kidneys[i]  = Krelease_t_Kidneys
    OP_IS$Krelease_t_Lungs[i] = Krelease_t_Lungs
    
    OP_IS$K50_t_Liver[i]  = K50_t_Liver
    OP_IS$K50_t_Spleen[i]  = K50_t_Spleen
    OP_IS$K50_t_Kidneys[i]  = K50_t_Kidneys
    OP_IS$K50_t_Lungs[i] = K50_t_Lungs
    
    OP_IS$blood_03[i] = blood[which.min(abs(times-0.03))]
    
    OP_IS$blood_1[i] = blood[which.min(abs(times-1))]
    
    OP_IS$blood_2[i] = blood[which.min(abs(times-2))]
    
    OP_IS$blood_4[i] = blood[which.min(abs(times-4))]
    
    OP_IS$blood_8[i] = blood[which.min(abs(times-8))]
    
    OP_IS$blood_10[i] = blood[which.min(abs(times-10))]
    
    OP_IS$blood_24[i] = blood[which.min(abs(times-24))]
    
    OP_IS$blood_48[i] = blood[which.min(abs(times-48))]
    
    OP_IS$Heart_3[i] = HV[which.min(abs(times-3))]
    
    OP_IS$Kidneys_3[i] = KV[which.min(abs(times-3))]
    
    OP_IS$Spleen_3[i] = SV[which.min(abs(times-3))]
    
    OP_IS$Lungs_3[i] = LuV[which.min(abs(times-3))]
    
    OP_IS$Bone_3[i] = BoneV[which.min(abs(times-3))]
    
    OP_IS$Heart_48[i] = HV[which.min(abs(times-48))]
    
    OP_IS$Kidneys_48[i] = KV[which.min(abs(times-48))]
    
    OP_IS$Spleen_48[i] = SV[which.min(abs(times-48))]
    
    OP_IS$Lungs_48[i] = LuV[which.min(abs(times-48))]
    
    OP_IS$Bone_48[i] = BoneV[which.min(abs(times-48))]
  }
  
  OP_IS$cond.blood_03_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_03_yn[OP_IS$blood_03 < blood_summary$UCI[blood_summary$time==0.03] & 
                           OP_IS$blood_03 > blood_summary$LCI[blood_summary$time==0.03]] <- 1
  
  OP_IS$cond.blood_1_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_1_yn[OP_IS$blood_1 < blood_summary$UCI[blood_summary$time==1] & 
                          OP_IS$blood_1 > blood_summary$LCI[blood_summary$time==1]] <- 1
  
  OP_IS$cond.blood_2_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_2_yn[OP_IS$blood_2 < blood_summary$UCI[blood_summary$time==2] & 
                          OP_IS$blood_2 > blood_summary$LCI[blood_summary$time==2]] <- 1
  
  OP_IS$cond.blood_4_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_4_yn[OP_IS$blood_4 < blood_summary$UCI[blood_summary$time==4] & 
                          OP_IS$blood_4 > blood_summary$LCI[blood_summary$time==4]] <- 1
  
  OP_IS$cond.blood_8_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_8_yn[OP_IS$blood_8 < blood_summary$UCI[blood_summary$time==8] & 
                          OP_IS$blood_8 > blood_summary$LCI[blood_summary$time==8]] <- 1
  
  OP_IS$cond.blood_10_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_10_yn[OP_IS$blood_10 < blood_summary$UCI[blood_summary$time==10] & 
                           OP_IS$blood_10 > blood_summary$LCI[blood_summary$time==10]] <- 1
  
  OP_IS$cond.blood_48_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_48_yn[OP_IS$blood_48< blood_summary$UCI[blood_summary$time==48] & 
                           OP_IS$blood_48 > blood_summary$LCI[blood_summary$time==48]] <- 1
  
  OP_IS$cond.blood_48_yn <- 0*OP_IS$nt
  OP_IS$cond.blood_48_yn[OP_IS$blood_48 < blood_summary$UCI[blood_summary$time==48] & 
                           OP_IS$blood_48 > blood_summary$LCI[blood_summary$time==48]] <- 1
  
  OP_IS$cond.blood.tot <- OP_IS$cond.blood_1_yn + OP_IS$cond.blood_8_yn + 
    OP_IS$cond.blood_48_yn + OP_IS$cond.blood_03_yn  
  
  OP_IS$cond.Heart_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Heart_3_yn[OP_IS$Heart_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Heart'] & 
                          OP_IS$Heart_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Heart']] <- 1
  
  OP_IS$cond.Kidneys_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Kidneys_3_yn[OP_IS$Kidneys_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Kidneys'] & 
                            OP_IS$Kidneys_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Kidneys']] <- 1
  
  OP_IS$cond.Lungs_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Lungs_3_yn[OP_IS$Lungs_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Lungs'] & 
                          OP_IS$Lungs_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Lungs']] <- 1
  
  OP_IS$cond.Bone_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Bone_3_yn[OP_IS$Bone_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Bone'] & 
                         OP_IS$Bone_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Bone']] <- 1
  
  OP_IS$cond.Spleen_3_yn <- 0*OP_IS$nt
  OP_IS$cond.Spleen_3_yn[OP_IS$Spleen_3 < organ_FACS_time3_summary$UCI[organ_FACS_time3_summary$Organ=='Spleen'] & 
                           OP_IS$Spleen_3 > organ_FACS_time3_summary$LCI[organ_FACS_time3_summary$Organ=='Spleen']] <- 1
  
  OP_IS$cond.organ.tot_3 <- OP_IS$cond.Spleen_3_yn + OP_IS$cond.Heart_3_yn + 
    OP_IS$cond.Kidneys_3_yn + OP_IS$cond.Lungs_3_yn + OP_IS$cond.Bone_3_yn
  
  OP_IS$cond.Heart_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Heart_48_yn[OP_IS$Heart_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Heart'] & 
                           OP_IS$Heart_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Heart']] <- 1
  
  OP_IS$cond.Kidneys_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Kidneys_48_yn[OP_IS$Kidneys_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Kidneys'] & 
                             OP_IS$Kidneys_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Kidneys']] <- 1
  
  OP_IS$cond.Lungs_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Lungs_48_yn[OP_IS$Lungs_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Lungs'] & 
                           OP_IS$Lungs_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Lungs']] <- 1
  
  OP_IS$cond.Bone_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Bone_48_yn[OP_IS$Bone_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Bone'] & 
                          OP_IS$Bone_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Bone']] <- 1
  
  OP_IS$cond.Spleen_48_yn <- 0*OP_IS$nt
  OP_IS$cond.Spleen_48_yn[OP_IS$Spleen_48 < organ_FACS_time48_summary$UCI[organ_FACS_time48_summary$Organ=='Spleen'] & 
                            OP_IS$Spleen_48 > organ_FACS_time48_summary$LCI[organ_FACS_time48_summary$Organ=='Spleen']] <- 1
  
  OP_IS$cond.organ.tot_48 <- OP_IS$cond.Spleen_48_yn + OP_IS$cond.Heart_48_yn + 
    OP_IS$cond.Kidneys_48_yn + OP_IS$cond.Lungs_48_yn + OP_IS$cond.Bone_48_yn
  
  OP_IS$PACt_Liver_scale <- log10(OP_IS$PACt_Liver/OP_IS$Kbile)
  OP_IS$Kbile_scale <- log10(OP_IS$Kbile/Kbile_Au)
  
  OP_IS$PACt_Bone_scale <- log10(OP_IS$PACt_Bone/PACt_Bone_Au)
  OP_IS$PACt_Heart_scale <- log10(OP_IS$PACt_Heart/PACt_Heart_Au)
  OP_IS$PACt_Spleen_scale <- log10(OP_IS$PACt_Spleen/PACt_Spleen_Au)
  OP_IS$PACt_Lungs_scale <- log10(OP_IS$PACt_Lungs/PACt_Lungs_Au)
  OP_IS$PACt_Rest_scale <- log10(OP_IS$PACt_Rest/PACt_Rest_Au)
  
  OP_IS$Kmax_t_Liver_scale <- log10(OP_IS$Kmax_t_Liver/Kmax_t_Liver_Au)
  OP_IS$Kmax_t_Spleen_scale <- log10(OP_IS$Kmax_t_Spleen/Kmax_t_Spleen_Au)
  OP_IS$Kmax_t_Kidneys_scale <- log10(OP_IS$Kmax_t_Kidneys/Kmax_t_Kidneys_Au)
  OP_IS$Kmax_t_Lungs_scale <- log10(OP_IS$Kmax_t_Lungs/Kmax_t_Lungs_Au)
  
  OP_IS$Krelease_t_Liver_scale <- log10(OP_IS$Krelease_t_Liver/OP_IS$Kmax_t_Liver)
  OP_IS$Krelease_t_Spleen_scale <- log10(OP_IS$Krelease_t_Spleen/OP_IS$Kmax_t_Spleen)
  OP_IS$Krelease_t_Kidneys_scale <- log10(OP_IS$Krelease_t_Kidneys/OP_IS$Kmax_t_Kidneys)
  OP_IS$Krelease_t_Lungs_scale <- log10(OP_IS$Krelease_t_Lungs/OP_IS$Kmax_t_Lungs)
  
  OP_IS$K50_t_Liver_scale <- log10(OP_IS$K50_t_Liver/K50_t_Liver_Au)
  OP_IS$K50_t_Spleen_scale <- log10(OP_IS$K50_t_Spleen/K50_t_Spleen_Au)
  OP_IS$K50_t_Kidneys_scale <- log10(OP_IS$K50_t_Kidneys/K50_t_Kidneys_Au)
  OP_IS$K50_t_Lungs_scale <- log10(OP_IS$K50_t_Lungs/K50_t_Lungs_Au)
  
  OP_IS$nt_scale <- log10(OP_IS$nt/nt_Au)
  
  OP_IS$dose <- rep(dose_rel[j],Nit)
  
  OP_IS_nana <- OP_IS[complete.cases(OP_IS),]
  
  if (j==1){
    
    OP_IS_tot <- OP_IS_nana
    
  }else{
    
    OP_IS_tot <- rbind(OP_IS_tot,OP_IS_nana)
    
  }
  
Blood <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Blood)})),nrow=Nt)

Bone <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Bone)})),nrow=Nt)

Kidneys <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Kidneys)})),nrow=Nt)

Liver <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Liver)})),nrow=Nt)

Spleen <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Spleen)})),nrow=Nt)

Lungs <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Lungs)})),nrow=Nt)

Heart <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Heart)})),nrow=Nt)

Brain <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Brain)})),nrow=Nt)

Rest <- matrix(unlist(lapply(OP_list,FUN=function(x){
  return(x$Rest)})),nrow=Nt)

OP_mn_sd <- data.frame(Blood_mn = rep(0,Nt),
                       Blood_sd = rep(0,Nt),
                       Heart_mn = rep(0,Nt),
                       Heart_sd = rep(0,Nt),
                       Liver_mn = rep(0,Nt),
                       Liver_sd = rep(0,Nt),
                       Bone_mn = rep(0,Nt),
                       Bone_sd = rep(0,Nt),
                       Kidneys_mn = rep(0,Nt),
                       Kidneys_sd = rep(0,Nt),
                       Spleen_mn = rep(0,Nt),
                       Spleen_sd = rep(0,Nt),
                       Lungs_mn = rep(0,Nt),
                       Lungs_sd = rep(0,Nt),
                       Rest_mn = rep(0,Nt),
                       Rest_sd = rep(0,Nt),
                       Brain_mn = rep(0,Nt),
                       Brain_sd = rep(0,Nt))

OP_mn_sd$Blood_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Blood[ind,]))
})
OP_mn_sd$Blood_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Blood[ind,]))
})

OP_mn_sd$Heart_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Heart[ind,]))
})
OP_mn_sd$Heart_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Heart[ind,]))
})

OP_mn_sd$Lungs_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Lungs[ind,]))
})
OP_mn_sd$Lungs_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Lungs[ind,]))
})

OP_mn_sd$Bone_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Bone[ind,]))
})
OP_mn_sd$Bone_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Bone[ind,]))
})

OP_mn_sd$Spleen_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Spleen[ind,]))
})
OP_mn_sd$Spleen_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Spleen[ind,]))
})

OP_mn_sd$Kidneys_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Kidneys[ind,]))
})
OP_mn_sd$Kidneys_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Kidneys[ind,]))
})

OP_mn_sd$Liver_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Liver[ind,]))
})
OP_mn_sd$Liver_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Liver[ind,]))
})

OP_mn_sd$Rest_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Rest[ind,]))
})
OP_mn_sd$Rest_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Rest[ind,]))
})
OP_mn_sd$Brain_mn <- sapply(seq(Nt),FUN=function(ind){
  return(mean(Brain[ind,]))
})
OP_mn_sd$Brain_sd <- sapply(seq(Nt),FUN=function(ind){
  return(sd(Brain[ind,]))
})

writeMat(paste('model_generated_data/sim_runs_',
               as.character(dose_rel[j]),'_20240118.mat'),OP=OP_mn_sd)
}


##########################################################
#Seventh section: Sensitivity analysis

#Clear all variables
rm(list=ls())

#Source scripts
source("libs.R")
source("parms.R")
alpha_temp_FACS <- 1e-2
alpha_temp_blood <- 3e-2
source("funcs.R")

OP_IS_tot <- readRDS("model_generated_data/parm2_MC_20240122.rds")
scales <- readRDS("model_generated_data/scales_20240122.rds")

OP_IS_max_cond <- OP_IS_tot[OP_IS_tot$cond.blood.tot%in%c(5),]
OP_IS_max_cond

OP_IS_max_cond_organs <- OP_IS_max_cond[OP_IS_max_cond$cond.organ.tot_3%in%c(3,4,5) &
                                          OP_IS_max_cond$cond.organ.tot_48%in%c(4,5),]
OPISmax_Heart <- OP_IS_max_cond[OP_IS_max_cond$cond.Heart_3_yn==1,]# & OP_IS_max_cond$cond.Heart_48_yn, ]
OPISmax_Lungs <- OP_IS_max_cond[OP_IS_max_cond$cond.Lungs_3_yn==1 & OP_IS_max_cond$cond.Lungs_48_yn==1, ]
OPISmax_Kidneys <- OP_IS_max_cond[OP_IS_max_cond$cond.Kidneys_3_yn==1 & OP_IS_max_cond$cond.Kidneys_48_yn==1, ]

scales <- data.frame(
  PACt_Liver=  mean(OP_IS_max_cond_organs$PACt_Liver_scale),
  Kbile= mean(OP_IS_max_cond_organs$Kbile_scale),
  
  PACt_Heart= mean(OPISmax_Heart$PACt_Heart_scale),
  PACt_Lungs= mean(OP_IS_max_cond_organs$PACt_Lungs_scale),
  PACt_Bone= mean(OP_IS_max_cond_organs$PACt_Bone_scale),
  PACt_Spleen= mean(OP_IS_max_cond_organs$PACt_Spleen_scale),
  PACt_Rest= mean(OP_IS_max_cond_organs$PACt_Rest_scale),
  
  Kmax_t_Liver= mean(OP_IS_max_cond_organs$Kmax_t_Liver_scale),
  Kmax_t_Spleen= mean(OP_IS_max_cond_organs$Kmax_t_Spleen_scale),
  Kmax_t_Lungs= mean(OP_IS_max_cond_organs$Kmax_t_Lungs_scale),
  Kmax_t_Kidneys= mean(OP_IS_max_cond_organs$Kmax_t_Kidneys_scale),
  
  Krelease_t_Liver= mean(OP_IS_max_cond_organs$Krelease_t_Liver_scale),
  Krelease_t_Spleen= mean(OP_IS_max_cond_organs$Krelease_t_Spleen_scale),
  Krelease_t_Lungs= mean(OP_IS_max_cond_organs$Krelease_t_Lungs_scale),
  Krelease_t_Kidneys= mean(OP_IS_max_cond_organs$Krelease_t_Kidneys_scale),
  
  K50_t_Liver= mean(OP_IS_max_cond_organs$K50_t_Liver_scale),
  K50_t_Spleen= mean(OP_IS_max_cond_organs$K50_t_Spleen_scale),
  K50_t_Lungs= mean(OP_IS_max_cond_organs$K50_t_Lungs_scale),
  K50_t_Kidneys= mean(OP_IS_max_cond_organs$K50_t_Kidneys_scale),
  
  nt = mean(OP_IS_max_cond_organs$nt_scale)
)

dist_setmod <- 1e5

ff_PACt_Bone <- fitdistr(OP_IS_max_cond_organs$PACt_Bone*dist_setmod, 'gamma')
ff_PACt_Spleen <- fitdistr(OP_IS_max_cond_organs$PACt_Spleen*dist_setmod, 'gamma')
ff_PACt_Lungs <- fitdistr(OPISmax_Lungs$PACt_Lungs*dist_setmod, 'gamma')
ff_PACt_Rest <- fitdistr(OP_IS_max_cond_organs$PACt_Rest*dist_setmod, 'gamma')

ff_Kbile <- fitdistr(OP_IS_max_cond_organs$Kbile*dist_setmod, 'gamma')
ff_nt <- fitdistr(OP_IS_max_cond_organs$nt, 'gamma')

ff_Kmax_t_Liver <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Liver, 'gamma')
ff_Kmax_t_Lungs <- fitdistr(OPISmax_Lungs$Kmax_t_Lungs, 'gamma')
ff_Kmax_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Kidneys, 'gamma')
ff_Kmax_t_Spleen <- fitdistr(OP_IS_max_cond_organs$Kmax_t_Spleen, 'gamma')

ff_K50_t_Liver <- fitdistr(OP_IS_max_cond_organs$K50_t_Liver, 'gamma')
ff_K50_t_Lungs <- fitdistr(OP_IS_max_cond_organs$K50_t_Lungs, 'gamma')
ff_K50_t_Kidneys <- fitdistr(OP_IS_max_cond_organs$K50_t_Kidneys, 'gamma')
ff_K50_t_Spleen <- fitdistr(OP_IS_max_cond_organs$K50_t_Spleen, 'gamma')


#Calculate parameter means
nt_mean <- ff_nt$estimate[1]/ff_nt$estimate[2]
Kbile_mean <- (min(OP_IS_max_cond_organs$Kbile)+max(OP_IS_max_cond_organs$Kbile))/2

PACt_Lungs_mean <- ff_PACt_Lungs$estimate[1]/ff_PACt_Lungs$estimate[2]/dist_setmod*PACt_Lungs_Mod
PACt_Heart_mean <- (min(OP_IS_max_cond_organs$PACt_Heart)+max(OP_IS_max_cond_organs$PACt_Heart))/2*PACt_Heart_Mod
PACt_Bone_mean <- ff_PACt_Bone$estimate[1]/ff_PACt_Bone$estimate[2]/dist_setmod
PACt_Spleen_mean <- ff_PACt_Spleen$estimate[1]/ff_PACt_Spleen$estimate[2]/dist_setmod*PACt_Spleen_Mod
PACt_Rest_mean <- ff_PACt_Rest$estimate[1]/ff_PACt_Rest$estimate[2]/dist_setmod

Kmax_t_Liver_mean <- ff_Kmax_t_Liver$estimate[1]/ff_Kmax_t_Liver$estimate[2]
Kmax_t_Kidneys_mean <- ff_Kmax_t_Kidneys$estimate[1]/ff_Kmax_t_Kidneys$estimate[2]*Kmax_t_Kidneys_Mod
Kmax_t_Lungs_mean <- ff_Kmax_t_Lungs$estimate[1]/ff_Kmax_t_Lungs$estimate[2]*Kmax_t_Lungs_Mod
Kmax_t_Spleen_mean <- ff_Kmax_t_Spleen$estimate[1]/ff_Kmax_t_Spleen$estimate[2]

K50_t_Liver_mean <- ff_K50_t_Liver$estimate[1]/ff_K50_t_Liver$estimate[2]
K50_t_Spleen_mean <- ff_K50_t_Spleen$estimate[1]/ff_K50_t_Spleen$estimate[2]
K50_t_Kidneys_mean <- ff_K50_t_Kidneys$estimate[1]/ff_K50_t_Kidneys$estimate[2]
K50_t_Lungs_mean <- ff_K50_t_Lungs$estimate[1]/ff_K50_t_Lungs$estimate[2]

PACt_Liver_mean <- Kbile_mean*10^scales$PACt_Liver

Krelease_t_Liver_mean <- Kmax_t_Liver_mean*10^scales$Krelease_t_Liver
Krelease_t_Spleen_mean <- Kmax_t_Spleen_mean*10^scales$Krelease_t_Spleen
Krelease_t_Kidneys_mean <- Kmax_t_Kidneys_mean*10^scales$Krelease_t_Kidneys
Krelease_t_Lungs_mean <- Kmax_t_Lungs_mean*10^scales$Krelease_t_Lungs

###Calculate standard deviations
nt_sd <- sqrt(ff_nt$estimate[1])/ff_nt$estimate[2]
Kbile_sd <- (-min(OP_IS_max_cond_organs$Kbile)+max(OP_IS_max_cond_organs$Kbile))/sqrt(12)

PACt_Lungs_sd <- sqrt(ff_PACt_Lungs$estimate[1])/ff_PACt_Lungs$estimate[2]/dist_setmod*PACt_Lungs_Mod
PACt_Heart_sd <- (-min(OP_IS_max_cond_organs$PACt_Heart)+max(OP_IS_max_cond_organs$PACt_Heart))/sqrt(12)*PACt_Heart_Mod
PACt_Bone_sd <- sqrt(ff_PACt_Bone$estimate[1])/ff_PACt_Bone$estimate[2]/dist_setmod
PACt_Spleen_sd <- sqrt(ff_PACt_Spleen$estimate[1])/ff_PACt_Spleen$estimate[2]/dist_setmod*PACt_Spleen_Mod
PACt_Rest_sd <- sqrt(ff_PACt_Rest$estimate[1])/ff_PACt_Rest$estimate[2]/dist_setmod

Kmax_t_Liver_sd <- sqrt(ff_Kmax_t_Liver$estimate[1])/ff_Kmax_t_Liver$estimate[2]
Kmax_t_Kidneys_sd <- sqrt(ff_Kmax_t_Kidneys$estimate[1])/ff_Kmax_t_Kidneys$estimate[2]*Kmax_t_Kidneys_Mod
Kmax_t_Lungs_sd <- sqrt(ff_Kmax_t_Lungs$estimate[1])/ff_Kmax_t_Lungs$estimate[2]*Kmax_t_Lungs_Mod
Kmax_t_Spleen_sd <- sqrt(ff_Kmax_t_Spleen$estimate[1])/ff_Kmax_t_Spleen$estimate[2]

K50_t_Liver_sd <- sqrt(ff_K50_t_Liver$estimate[1])/ff_K50_t_Liver$estimate[2]
K50_t_Spleen_sd <- sqrt(ff_K50_t_Spleen$estimate[1])/ff_K50_t_Spleen$estimate[2]
K50_t_Kidneys_sd <- sqrt(ff_K50_t_Kidneys$estimate[1])/ff_K50_t_Kidneys$estimate[2]
K50_t_Lungs_sd <- sqrt(ff_K50_t_Lungs$estimate[1])/ff_K50_t_Lungs$estimate[2]

PACt_Liver_sd <- Kbile_sd*10^scales$PACt_Liver

Krelease_t_Liver_sd <- Kmax_t_Liver_sd*10^scales$Krelease_t_Liver
Krelease_t_Spleen_sd <- Kmax_t_Spleen_sd*10^scales$Krelease_t_Spleen
Krelease_t_Kidneys_sd <- Kmax_t_Kidneys_sd*10^scales$Krelease_t_Kidneys
Krelease_t_Lungs_sd <- Kmax_t_Lungs_sd*10^scales$Krelease_t_Lungs

####Report Mean values
nt_mean 
Kbile_mean 

PACt_Lungs_mean 
PACt_Heart_mean 
PACt_Bone_mean 
PACt_Spleen_mean 
PACt_Rest_mean 

Kmax_t_Liver_mean 
Kmax_t_Kidneys_mean 
Kmax_t_Lungs_mean 
Kmax_t_Spleen_mean 

K50_t_Liver_mean 
K50_t_Spleen_mean 
K50_t_Kidneys_mean 
K50_t_Lungs_mean 

PACt_Liver_mean 

Krelease_t_Liver_mean 
Krelease_t_Spleen_mean 
Krelease_t_Kidneys_mean 
Krelease_t_Lungs_mean 

###Report Standard deviations
nt_sd 
Kbile_sd 

PACt_Lungs_sd 
PACt_Heart_sd 
PACt_Bone_sd 
PACt_Spleen_sd 
PACt_Rest_sd 

Kmax_t_Liver_sd 
Kmax_t_Kidneys_sd 
Kmax_t_Lungs_sd 
Kmax_t_Spleen_sd 

K50_t_Liver_sd 
K50_t_Spleen_sd 
K50_t_Kidneys_sd 
K50_t_Lungs_sd 

PACt_Liver_sd 

Krelease_t_Liver_sd 
Krelease_t_Spleen_sd 
Krelease_t_Kidneys_sd 
Krelease_t_Lungs_sd

###Calculate Coefficients of Variation
nt_sd/nt_mean 
Kbile_sd/Kbile_mean 

PACt_Lungs_sd/PACt_Lungs_mean 
PACt_Heart_sd/PACt_Heart_mean 
PACt_Bone_sd/PACt_Bone_mean 
PACt_Spleen_sd/PACt_Spleen_mean 
PACt_Rest_sd/PACt_Rest_mean 

Kmax_t_Liver_sd/Kmax_t_Liver_mean 
Kmax_t_Kidneys_sd/Kmax_t_Kidneys_mean 
Kmax_t_Lungs_sd/Kmax_t_Lungs_mean 
Kmax_t_Spleen_sd/Kmax_t_Spleen_mean 

K50_t_Liver_sd/K50_t_Liver_mean 
K50_t_Spleen_sd /K50_t_Spleen_mean 
K50_t_Kidneys_sd/K50_t_Kidneys_mean 
K50_t_Lungs_sd/K50_t_Lungs_mean 

Init_condition <- rep(0,22)

OP_list <- list()

nt <- nt_mean
Kbile <- Kbile_mean

PACt_Lungs <- PACt_Lungs_mean
PACt_Heart <- PACt_Heart_mean
PACt_Bone <- PACt_Bone_mean
PACt_Spleen <- PACt_Spleen_mean
PACt_Rest <- PACt_Rest_mean

Kmax_t_Liver <- Kmax_t_Liver_mean
Kmax_t_Kidneys <- Kmax_t_Kidneys_mean
Kmax_t_Spleen <- Kmax_t_Spleen_mean
Kmax_t_Lungs <- Kmax_t_Lungs_mean

K50_t_Liver <- K50_t_Liver_mean
K50_t_Spleen <- K50_t_Spleen_mean
K50_t_Kidneys <- K50_t_Kidneys_mean
K50_t_Lungs <- K50_t_Lungs_mean

PACt_Liver <- PACt_Liver_mean

Krelease_t_Liver <- Krelease_t_Liver_mean
Krelease_t_Spleen <- Krelease_t_Spleen_mean
Krelease_t_Kidneys <- Krelease_t_Kidneys_mean
Krelease_t_Lungs <- Krelease_t_Lungs_mean

Kurine <- 0
PACt_Kidneys <- PACt_Kidneys_Mod
PACt_Brain <- PACt_Rest_mean

OP <- ode(y=Init_condition,
          times=times,
          func = NP_PBPK2, 
          parms=list(),
          method='lsoda')
OP <- data.frame(OP)
OP_rename <- OP

indx_samp <- seq(from = 1, to = Nt, by = 1)

OP_rename$M_Art <- OP[indx_samp,2]
OP_rename$M_Ven <- OP[indx_samp,3]

OP_rename$M_vasc_Liver <- OP[indx_samp,4]
OP_rename$M_extra_Liver <- OP[indx_samp,5]
OP_rename$M_phago_Liver <- OP[indx_samp,6]

OP_rename$M_vasc_Spleen <- OP[indx_samp,7]
OP_rename$M_extra_Spleen <- OP[indx_samp,8]
OP_rename$M_phago_Spleen <- OP[indx_samp,9]

OP_rename$M_vasc_Kidneys <- OP[indx_samp,10]
OP_rename$M_extra_Kidneys <- OP[indx_samp,11]
OP_rename$M_phago_Kidneys <- OP[indx_samp,12]

OP_rename$M_vasc_Lungs <- OP[indx_samp,13]
OP_rename$M_extra_Lungs <- OP[indx_samp,14]
OP_rename$M_phago_Lungs <- OP[indx_samp,15]

OP_rename$M_vasc_Brain <- OP[indx_samp,16]
OP_rename$M_extra_Brain <- OP[indx_samp,17]

OP_rename$M_vasc_Heart <- OP[indx_samp,18]
OP_rename$M_extra_Heart <- OP[indx_samp,19]

OP_rename$M_vasc_Bone <- OP[indx_samp,20]
OP_rename$M_extra_Bone <- OP[indx_samp,21]

OP_rename$M_vasc_Rest <- OP[indx_samp,22]
OP_rename$M_extra_Rest <- OP[indx_samp,23]

LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/LV
LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/LV
SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/LV
HV <- OP_rename$M_extra_Heart/LV
BoneV <- OP_rename$M_extra_Bone/LV
LV <- 1

OP_rename$LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
OP_rename$KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/OP_rename$LV
OP_rename$LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/OP_rename$LV
OP_rename$SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/OP_rename$LV
OP_rename$HV <- OP_rename$M_extra_Heart/OP_rename$LV
OP_rename$BoneV <- OP_rename$M_extra_Bone/OP_rename$LV

blood<- (OP_rename$M_Ven)/(BW*VVenC)
OP_rename$Blood <- (OP_rename$M_Ven)/(BW*VVenC)
OP_rename$Liver <- (OP_rename$M_extra_Liver+OP_rename$M_phago_Liver)/(BW*VLC)
OP_rename$Spleen <- (OP_rename$M_extra_Spleen+OP_rename$M_phago_Spleen)/(BW*VSC)
OP_rename$Kidneys <- (OP_rename$M_extra_Kidneys+OP_rename$M_phago_Kidneys)/(BW*VKC)
OP_rename$Lungs <- (OP_rename$M_extra_Lungs+OP_rename$M_phago_Lungs)/(BW*VLuC)
OP_rename$Heart <- (OP_rename$M_extra_Heart)/(BW*VHC)
OP_rename$Bone <- (OP_rename$M_extra_Bone)/(BW*VBoneC)  

AUC <- dt/2*(OP_rename$Blood[1]+OP_rename$Blood[Nt]) + dt*sum(OP_rename$Blood[2:(Nt-1)])

parm_mod <- rep(1,20)
AUC_keep <- 0*parm_mod

for (i in seq(parm_mod)){
  parm_mod <- rep(1,20)
  parm_mod[i] <- 1.01
  
  nt <- nt_mean*parm_mod[1]
  Kbile <- Kbile_mean*parm_mod[2]
  
  PACt_Lungs <- PACt_Lungs_mean*parm_mod[3]
  PACt_Heart <- PACt_Heart_mean*parm_mod[4]
  PACt_Bone <- PACt_Bone_mean*parm_mod[5]
  PACt_Spleen <- PACt_Spleen_mean*parm_mod[6]
  PACt_Rest <- PACt_Rest_mean*parm_mod[7]
  
  Kmax_t_Liver <- Kmax_t_Liver_mean*parm_mod[8]
  Kmax_t_Kidneys <- Kmax_t_Kidneys_mean*parm_mod[9]
  Kmax_t_Spleen <- Kmax_t_Spleen_mean*parm_mod[10]
  Kmax_t_Lungs <- Kmax_t_Lungs_mean*parm_mod[11]
  
  K50_t_Liver <- K50_t_Liver_mean*parm_mod[12]
  K50_t_Spleen <- K50_t_Spleen_mean*parm_mod[13]
  K50_t_Kidneys <- K50_t_Kidneys_mean*parm_mod[14]
  K50_t_Lungs <- K50_t_Lungs_mean*parm_mod[15]
  
  PACt_Liver <- PACt_Liver_mean*parm_mod[16]

  Krelease_t_Liver <- Krelease_t_Liver_mean*parm_mod[17]
  Krelease_t_Spleen <- Krelease_t_Spleen_mean*parm_mod[18]
  Krelease_t_Kidneys <- Krelease_t_Kidneys_mean*parm_mod[19]
  Krelease_t_Lungs <- Krelease_t_Lungs_mean*parm_mod[20]
  
  Kurine <- 0
  PACt_Kidneys <- PACt_Kidneys_Mod
  PACt_Brain <- PACt_Rest_mean
  
  OP <- ode(y=Init_condition,
            times=times,
            func = NP_PBPK2, 
            parms=list(),
            method='lsoda')
  OP <- data.frame(OP)
  OP_rename <- OP
  
  indx_samp <- seq(from = 1, to = Nt, by = 1)
  
  OP_rename$M_Art <- OP[indx_samp,2]
  OP_rename$M_Ven <- OP[indx_samp,3]
  
  OP_rename$M_vasc_Liver <- OP[indx_samp,4]
  OP_rename$M_extra_Liver <- OP[indx_samp,5]
  OP_rename$M_phago_Liver <- OP[indx_samp,6]
  
  OP_rename$M_vasc_Spleen <- OP[indx_samp,7]
  OP_rename$M_extra_Spleen <- OP[indx_samp,8]
  OP_rename$M_phago_Spleen <- OP[indx_samp,9]
  
  OP_rename$M_vasc_Kidneys <- OP[indx_samp,10]
  OP_rename$M_extra_Kidneys <- OP[indx_samp,11]
  OP_rename$M_phago_Kidneys <- OP[indx_samp,12]
  
  OP_rename$M_vasc_Lungs <- OP[indx_samp,13]
  OP_rename$M_extra_Lungs <- OP[indx_samp,14]
  OP_rename$M_phago_Lungs <- OP[indx_samp,15]
  
  OP_rename$M_vasc_Brain <- OP[indx_samp,16]
  OP_rename$M_extra_Brain <- OP[indx_samp,17]
  
  OP_rename$M_vasc_Heart <- OP[indx_samp,18]
  OP_rename$M_extra_Heart <- OP[indx_samp,19]
  
  OP_rename$M_vasc_Bone <- OP[indx_samp,20]
  OP_rename$M_extra_Bone <- OP[indx_samp,21]
  
  OP_rename$M_vasc_Rest <- OP[indx_samp,22]
  OP_rename$M_extra_Rest <- OP[indx_samp,23]
  
  LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
  KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/LV
  LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/LV
  SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/LV
  HV <- OP_rename$M_extra_Heart/LV
  BoneV <- OP_rename$M_extra_Bone/LV
  LV <- 1
  
  OP_rename$LV <- (OP_rename$M_phago_Liver+OP_rename$M_extra_Liver)
  OP_rename$KV <- (OP_rename$M_phago_Kidneys+OP_rename$M_extra_Kidneys)/OP_rename$LV
  OP_rename$LuV <- (OP_rename$M_phago_Lungs+OP_rename$M_extra_Lungs)/OP_rename$LV
  OP_rename$SV <- (OP_rename$M_phago_Spleen+OP_rename$M_extra_Spleen)/OP_rename$LV
  OP_rename$HV <- OP_rename$M_extra_Heart/OP_rename$LV
  OP_rename$BoneV <- OP_rename$M_extra_Bone/OP_rename$LV
  
  blood<- (OP_rename$M_Ven)/(BW*VVenC)
  OP_rename$Blood <- (OP_rename$M_Ven+OP_rename$M_Art)/(BW*VBloodC)
  OP_rename$Liver <- (OP_rename$M_extra_Liver+OP_rename$M_phago_Liver)/(BW*VLC)
  OP_rename$Spleen <- (OP_rename$M_extra_Spleen+OP_rename$M_phago_Spleen)/(BW*VSC)
  OP_rename$Kidneys <- (OP_rename$M_extra_Kidneys+OP_rename$M_phago_Kidneys)/(BW*VKC)
  OP_rename$Lungs <- (OP_rename$M_extra_Lungs+OP_rename$M_phago_Lungs)/(BW*VLuC)
  OP_rename$Heart <- (OP_rename$M_extra_Heart)/(BW*VHC)
  OP_rename$Bone <- (OP_rename$M_extra_Bone)/(BW*VBoneC)  
  
  AUC_i <- dt/2*(OP_rename$Blood[1]+OP_rename$Blood[Nt]) + dt*sum(OP_rename$Blood[2:(Nt-1)])
  AUC_keep[i] <- AUC_i
}

NSC <- (AUC_keep-AUC)/AUC/0.01
which(abs(NSC)>0.2)

NSC

