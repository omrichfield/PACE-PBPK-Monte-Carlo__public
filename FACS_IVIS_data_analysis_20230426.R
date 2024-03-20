#FACS and IVIS data analysis
#Owen Richfield
#12/5/2022

#Clear all variables
#rm(list=ls())

#source("parms.R")
#source("funcs.R")

n <- 3
tol_fac <- 2

##############################################################################################
#Organ NP FACS
alpha <- alpha_temp_FACS

organ_FACS_temp <- data.frame(read_excel("data_for_model/organ_FACS_MFI_proc_all_forms.xlsx"))

organ_FACS <- organ_FACS_temp[organ_FACS_temp$Form=='PACE_PEG',]

organ_vol_vec <- c(rep(CNH,n),
                   rep(CNLu,n),
                   rep(CNL,n),
                   rep(CNS,n),
                   rep(CNK,n),
                   rep(CNBone,n))
# c(rep(VHC,n),
#   rep(VLuC,n),
#   rep(VLC,n),
#   rep(VSC,n),
#   rep(VKC,n),
#   rep(VBoneC,n))

organ_FACS$Fluor_vol <- organ_FACS$MFI*organ_vol_vec
organ_FACS$Liver_norm <- organ_FACS$Fluor_vol*0
organ_FACS_Liver_Fluor <- organ_FACS$Fluor_vol[organ_FACS$Organ=='Liver']
organ_FACS_Doses <- unique(organ_FACS$Dose)
num_Doses <- length(unique(organ_FACS$Dose))


for (i in seq(num_Doses)){
  indx_temp <- organ_FACS$Dose==organ_FACS_Doses[i]
  organ_FACS_di <- organ_FACS[indx_temp,]
  liver_indx <- organ_FACS_di$Organ=='Liver'
  organ_FACS$Liver_norm[indx_temp] <- organ_FACS_di$Fluor_vol/organ_FACS_di$Fluor_vol[liver_indx]
}

FACS_Organ <- aggregate(organ_FACS$Liver_norm, list(Organ = organ_FACS$Organ, Dose=organ_FACS$Dose),mean)$Organ
FACS_dose <- aggregate(organ_FACS$Liver_norm, list(Organ = organ_FACS$Organ, Dose=organ_FACS$Dose),mean)$Dose

FACS_mean <- aggregate(organ_FACS$Liver_norm, list(Organ = organ_FACS$Organ, Dose=organ_FACS$Dose),mean)$x
FACS_sd <- aggregate(organ_FACS$Liver_norm, list(Organ = organ_FACS$Organ, Dose=organ_FACS$Dose),sd)$x

organ_FACS_summary <- data.frame(Dose=FACS_dose, Organ=FACS_Organ, mean=FACS_mean, sd=FACS_sd)

organ_FACS_summary$std_err <- organ_FACS_summary$sd/sqrt(n)
organ_FACS_summary$t_score <- qt(p=alpha/2, df=n-1,lower.tail=F)
organ_FACS_summary$margin_err <- organ_FACS_summary$std_err*organ_FACS_summary$t_score
organ_FACS_summary$LCI <- organ_FACS_summary$mean -  organ_FACS_summary$margin_err
organ_FACS_summary$UCI <- organ_FACS_summary$mean +  organ_FACS_summary$margin_err

writeMat('model_generated_data/PACE_FACS_summary_20240109.mat',FACS_summary=organ_FACS_summary)

# ##############################################################################################
#Organ NP FACS in time

organ_FACS <- data.frame(read_excel("data_for_model/organ_FACS_MFI_time.xlsx"))

organ_vol_vec <- c(rep(CNH,n),
                   rep(CNLu,n),
                   rep(CNL,n),
                   rep(CNS,n),
                   rep(CNK,n),
                   rep(CNBone,n))
# c(rep(VHC,n),
#   rep(VLuC,n),
#   rep(VLC,n),
#   rep(VSC,n),
#   rep(VKC,n),
#   rep(VBoneC,n))

organ_FACS$Fluor_vol <- organ_FACS$MFI*organ_vol_vec
organ_FACS$Liver_norm <- organ_FACS$Fluor_vol*0
organ_FACS_Liver_Fluor <- organ_FACS$Fluor_vol[organ_FACS$Organ=='Liver']
organ_FACS_times <- unique(organ_FACS$time)
num_times <- length(unique(organ_FACS$time))


for (i in seq(num_times)){
  indx_temp <- organ_FACS$time==organ_FACS_times[i]
  organ_FACS_di <- organ_FACS[indx_temp,]
  liver_indx <- organ_FACS_di$Organ=='Liver'
  organ_FACS$Liver_norm[indx_temp] <- organ_FACS_di$Fluor_vol/organ_FACS_di$Fluor_vol[liver_indx]
}

FACS_Organ <- aggregate(organ_FACS$Liver_norm, list(Organ = organ_FACS$Organ, time=organ_FACS$time),mean)$Organ
FACS_time <- aggregate(organ_FACS$Liver_norm, list(Organ = organ_FACS$Organ, time=organ_FACS$time),mean)$time

FACS_mean <- aggregate(organ_FACS$Liver_norm, list(Organ = organ_FACS$Organ, time=organ_FACS$time),mean)$x
FACS_sd <- aggregate(organ_FACS$Liver_norm, list(Organ = organ_FACS$Organ, time=organ_FACS$time),sd)$x

organ_FACS_time_summary <- data.frame(time=FACS_time, Organ=FACS_Organ, mean=FACS_mean, sd=FACS_sd)

organ_FACS_time_summary$std_err <- organ_FACS_time_summary$sd/sqrt(n)
organ_FACS_time_summary$t_score <- qt(p=alpha/2, df=n-1,lower.tail=F)
organ_FACS_time_summary$margin_err <- organ_FACS_time_summary$std_err*organ_FACS_time_summary$t_score
organ_FACS_time_summary$LCI <- organ_FACS_time_summary$mean -  organ_FACS_time_summary$margin_err
organ_FACS_time_summary$UCI <- organ_FACS_time_summary$mean +  organ_FACS_time_summary$margin_err

writeMat('model_generated_data/PACE_FACS_summary_time_20240109.mat',FACS_summary_time=organ_FACS_time_summary)

##############################################################################################
#Blood Concentration of NPs
alpha <- alpha_temp_blood

blood_conc_temp <- data.frame(read_excel("data_for_model/blood_PACEPEG_conc.xlsx"))

blood_conc <- data.frame(time=rep(blood_conc_temp$time,5),
                         dose=c(rep(0.1,length(blood_conc_temp$time)),
                                rep(0.5,length(blood_conc_temp$time)),
                                rep(2,length(blood_conc_temp$time)),
                                rep(2.5,length(blood_conc_temp$time)),
                                rep(3.5,length(blood_conc_temp$time))),
                         conc=c(blood_conc_temp$PACEPEG_01mg,
                                blood_conc_temp$PACEPEG_05mg,
                                blood_conc_temp$PACEPEG_2mg,
                                blood_conc_temp$PACEPEG_25mg,
                                blood_conc_temp$PACEPEG_35mg))

blood_time <- aggregate(blood_conc, list(time = blood_conc$time, dose = blood_conc$dose),mean)$time
blood_dose <- aggregate(blood_conc, list(time = blood_conc$time, dose = blood_conc$dose),mean)$dose

blood_mean <- aggregate(blood_conc$conc, list(time = blood_conc$time, dose = blood_conc$dose),mean,na.rm=T)$x
blood_sd <- aggregate(blood_conc$conc, list(time = blood_conc$time, dose = blood_conc$dose),sd,na.rm=T)$x

blood_summary <- data.frame(dose = blood_dose, time=blood_time, mean=blood_mean, sd=blood_sd)

blood_summary$std_err <- blood_summary$sd/sqrt(n)
blood_summary$t_score <- qt(p=alpha/2, df=n-1,lower.tail=F)
blood_summary$margin_err <- blood_summary$std_err*blood_summary$t_score
blood_summary$LCI <- blood_summary$mean -  blood_summary$margin_err
blood_summary$UCI <- blood_summary$mean +  blood_summary$margin_err
blood_summary <- na.omit(blood_summary)

writeMat('model_generated_data/PACE_blood_summary_20240109.mat',blood_summary=blood_summary)

##############################################################################################
#Plot the 
dd=data.frame(time=c(organ_FACS_time_summary$time,
                     rep(48,length(organ_FACS_summary$mean[organ_FACS_summary$Dose==0.5]))), 
              mean=c(organ_FACS_time_summary$mean, 
                     organ_FACS_summary$mean[organ_FACS_summary$Dose==0.5]),
              Organ=c(organ_FACS_time_summary$Organ, 
                      organ_FACS_summary$Organ[organ_FACS_summary$Dose==0.5]))

ggplot(data=dd,aes(x=time,y=mean,color=Organ))+
              geom_line()+ylim(c(0,0.3))+
              ylab('NP mass (rel to Liver, by cell #)')



#Liver NP FACS

# liver_FACS <- data.frame(read_excel("liver_lung_FACS_proc.xlsx"))
# 
# liver_cellfrac_vec <- c(rep(0.18,n),rep(0.22,n))*V_Liver
# 
# liver_FACS$Fluor_vol <- liver_FACS$MFI*liver_cellfrac_vec
# 
# liver_FACS_endo <- liver_FACS[liver_FACS$Organ.Cell=="Liver.Endo",]
# liver_FACS_phago <- liver_FACS[liver_FACS$Organ.Cell=="Liver.Phago",]
# 
# FACS_liver <- data.frame(Dose=liver_FACS_endo$Dose,
#                          endo.phago.norm = liver_FACS_endo$Fluor_vol/liver_FACS_phago$Fluor_vol)
# 
# FACS_dose <- aggregate(FACS_liver$endo.phago.norm, list(Dose=FACS_liver$Dose),mean)$Dose
# FACS_mean <- aggregate(FACS_liver$endo.phago.norm, list(Dose = FACS_liver$Dose),mean)$x
# FACS_sd <- aggregate(FACS_liver$endo.phago.norm, list(Dose = FACS_liver$Dose),sd)$x
# 
# liver_FACS_summary <- data.frame(Dose=FACS_dose, mean=FACS_mean, sd=FACS_sd)
# 
# liver_FACS_summary$std_err <- liver_FACS_summary$sd/sqrt(n)
# liver_FACS_summary$t_score <- qt(p=alpha/2, df=n-1,lower.tail=F)
# liver_FACS_summary$margin_err <- liver_FACS_summary$std_err*liver_FACS_summary$t_score
# liver_FACS_summary$LCI <- liver_FACS_summary$mean -  liver_FACS_summary$margin_err
# liver_FACS_summary$UCI <- liver_FACS_summary$mean +  liver_FACS_summary$margin_err

##################################################################################
  
# 
# 
# OPT <- readRDS("OPT_20221209_2.rds")
# #OPT <- readRDS("OPT_20221209.rds")
# 
# compare_sim_data <- data.frame()
# 
# 
# cond_Heart_UCI <- OPT$Heart_norm < organ_FACS_summary$UCI[organ_FACS_summary$Organ=='Heart']  
# cond_Heart_LCI <- OPT$Heart_norm > organ_FACS_summary$LCI[organ_FACS_summary$Organ=='Heart']  
# 
# which(cond_Heart_LCI & cond_Heart_UCI)
# OPT$cond_Heart <- cond_Heart_LCI & cond_Heart_UCI
# 
# 
# cond_Kidney_UCI <- OPT$Kidney_norm < organ_FACS_summary$UCI[organ_FACS_summary$Organ=='Kidneys']  
# cond_Kidney_LCI <- OPT$Kidney_norm > organ_FACS_summary$LCI[organ_FACS_summary$Organ=='Kidneys']  
# 
# which(cond_Kidney_LCI & cond_Kidney_UCI)
# OPT$cond_Kidney <- cond_Kidney_LCI & cond_Kidney_UCI
# 
# 
# cond_Lung_UCI <- OPT$Lung_norm < organ_FACS_summary$UCI[organ_FACS_summary$Organ=='Lungs']  
# cond_Lung_LCI <- OPT$Lung_norm > organ_FACS_summary$LCI[organ_FACS_summary$Organ=='Lungs']  
# 
# which(cond_Lung_LCI & cond_Lung_UCI)
# OPT$cond_Lung <- cond_Lung_LCI & cond_Lung_UCI
# 
# 
# cond_Spleen_UCI <- OPT$Spleen_norm < organ_FACS_summary$UCI[organ_FACS_summary$Organ=='Spleen']  
# cond_Spleen_LCI <- OPT$Spleen_norm > organ_FACS_summary$LCI[organ_FACS_summary$Organ=='Spleen']  
# 
# which(cond_Spleen_LCI & cond_Spleen_UCI)
# OPT$cond_Spleen <- cond_Spleen_LCI & cond_Spleen_UCI
# 
# cond_Other_UCI <- OPT$Other_norm < organ_FACS_summary$UCI[organ_FACS_summary$Organ=='Bone']  
# cond_Other_LCI <- OPT$Other_norm > organ_FACS_summary$LCI[organ_FACS_summary$Organ=='Bone']  
# 
# which(cond_Other_LCI & cond_Other_UCI)
# OPT$cond_Other <- cond_Other_LCI & cond_Other_UCI
# 
# cond_blood_UCI <- OPT$blood_4 < blood_summary$UCI[blood_summary$time==4] 
# cond_blood_LCI <- OPT$blood_4 > blood_summary$LCI[blood_summary$time==4] 
# 
# which(cond_blood_LCI & cond_blood_UCI)
# OPT$cond_blood_4 <- cond_blood_LCI & cond_blood_UCI
# 
# cond_blood_UCI <- OPT$blood_48 < blood_summary$UCI[blood_summary$time==48] 
# cond_blood_LCI <- OPT$blood_48 > blood_summary$LCI[blood_summary$time==48] 
# 
# which(cond_blood_LCI & cond_blood_UCI)
# OPT$cond_blood_48 <- cond_blood_LCI & cond_blood_UCI
# 
# OPT$cond <- OPT$cond_blood_4 & OPT$cond_blood_48 & OPT$cond_Heart & OPT$cond_Kidney & OPT$cond_Lung & OPT$cond_Other & OPT$cond_Spleen
# 
# OPT_cond <- OPT[OPT$cond,]
# 
# ################################################################################
# set.seed(6)
# code_time_start <- proc.time()
# dose_vec <- c(0.1,0.5,2,2.5)*1000
# 
# kve_Heart_keep <- rep(0,N)
# kve_Kidney_keep <- rep(0,N)
# kve_Liver_keep <- rep(0,N)
# kve_Lung_keep <- rep(0,N)
# kve_Spleen_keep <- rep(0,N)
# kve_Other_keep <- rep(0,N)
# kexc_Kidney_keep <- rep(0,N)
# kexc_Liver_keep <- rep(0,N)
# kvp_Liver_keep <- rep(0,N)
# kvp_Spleen_keep <- rep(0,N)
# kvp_Lung_keep <- rep(0,N)
# MPmax_Liver_keep <- rep(0,N)
# MPmax_Spleen_keep <- rep(0,N)
# MPmax_Lung_keep <- rep(0,N)
# kdeg_Liver_keep <- rep(0,N)
# kdeg_Spleen_keep <- rep(0,N)
# kdeg_Lung_keep <- rep(0,N)
# 
# OP_keep <- list()
# 
# OP_blood_4_keep <- rep(0,N)
# OP_blood_48_keep <- rep(0,N)
# OP_Heart_Liver_norm <- rep(0,N)
# OP_Kidney_Liver_norm <- rep(0,N)
# OP_Lung_Liver_norm <- rep(0,N)
# OP_Spleen_Liver_norm <- rep(0,N)
# OP_Other_Liver_norm <- rep(0,N)
# 
# for (i in seq(N)){
#   
#   kve_Heart <- abs(rnorm(1,mean(OPT_cond$kve_Heart),var(OPT_cond$kve_Heart)))     #kve
#   kve_Kidney <- abs(rnorm(1,mean(OPT_cond$kve_Kidney),var(OPT_cond$kve_Kidney))) 
#   kve_Liver <- abs(rnorm(1,mean(OPT_cond$kve_Liver),var(OPT_cond$kve_Liver))) 
#   kve_Lung <- abs(rnorm(1,mean(OPT_cond$kve_Lung),var(OPT_cond$kve_Lung))) 
#   kve_Spleen <- abs(rnorm(1,mean(OPT_cond$kve_Spleen),var(OPT_cond$kve_Spleen))) 
#   kve_Other <- abs(rnorm(1,mean(OPT_cond$kve_Other),var(OPT_cond$kve_Other))) 
#   
#   kexc_Kidney <- abs(rnorm(1,mean(OPT_cond$kexc_Kidney),var(OPT_cond$kexc_Kidney))) 
#   kexc_Liver <- abs(rnorm(1,mean(OPT_cond$kexc_Liver),var(OPT_cond$kexc_Liver))) 
#   
#   kvp_Liver <- abs(rnorm(1,mean(OPT_cond$kvp_Liver),var(OPT_cond$kvp_Liver)))        #kvp
#   kvp_Spleen <- abs(rnorm(1,mean(OPT_cond$kvp_Spleen),var(OPT_cond$kvp_Spleen)))
#   kvp_Lung <- abs(rnorm(1,mean(OPT_cond$kvp_Lung),var(OPT_cond$kvp_Lung)))
#   
#   MPmax_Liver <- abs(rnorm(1,mean(OPT_cond$MPmax_Liver),var(OPT_cond$MPmax_Liver))) #MPMax
#   MPmax_Spleen <- abs(rnorm(1,mean(OPT_cond$MPmax_Spleen),var(OPT_cond$MPmax_Spleen)))
#   MPmax_Lung <- abs(rnorm(1,mean(OPT_cond$MPmax_Lung),var(OPT_cond$MPmax_Lung)))
#   
#   kdeg_Liver <- abs(rnorm(1,mean(OPT_cond$kdeg_Liver),var(OPT_cond$kdeg_Liver)))      #kdeg
#   kdeg_Spleen <- abs(rnorm(1,mean(OPT_cond$kdeg_Spleen),var(OPT_cond$kdeg_Spleen)))
#   kdeg_Lung <- abs(rnorm(1,mean(OPT_cond$kdeg_Lung),var(OPT_cond$kdeg_Lung)))
#   
#   OP <- ode(y=Init_condition,
#             times=times,
#             func = NP_PBPK, 
#             parms=list(),
#             method='lsoda')
#   OP <- data.frame(OP)
#   OP_keep[[i]] <- OP
#   
#   kve_Heart_keep[i] <- kve_Heart
#   kve_Kidney_keep[i] <- kve_Kidney
#   kve_Liver_keep[i] <- kve_Liver
#   kve_Lung_keep[i] <- kve_Lung
#   kve_Spleen_keep[i] <- kve_Spleen
#   kve_Other_keep[i] <- kve_Other
#   kexc_Kidney_keep[i] <- kexc_Kidney
#   kexc_Liver_keep[i] <- kexc_Liver
#   kvp_Liver_keep[i] <- kvp_Liver
#   kvp_Spleen_keep[i] <- kvp_Spleen
#   kvp_Lung_keep[i] <- kvp_Lung
#   MPmax_Liver_keep[i] <- MPmax_Liver
#   MPmax_Spleen_keep[i] <- MPmax_Spleen
#   MPmax_Lung_keep[i] <- MPmax_Lung
#   kdeg_Liver_keep[i] <- kdeg_Liver
#   kdeg_Spleen_keep[i] <- kdeg_Spleen
#   kdeg_Lung_keep[i] <- kdeg_Lung
#   
#   OP_blood_4_keep[i] <- OP$X2[OP$time==4]/V_Ven
#   OP_blood_48_keep[i] <- OP$X2[OP$time==48]/V_Ven
#   
#   OP_48 <- OP[OP$time==48,]
#   
#   OP_Heart_Liver_norm[i] <- OP_48$X4/(OP_48$X8+OP_48$X9)
#   
#   OP_Kidney_Liver_norm[i] <- OP_48$X6/(OP_48$X8+OP_48$X9)
#   
#   OP_Spleen_Liver_norm[i] <- (OP_48$X11+OP_48$X12)/(OP_48$X8+OP_48$X9)
#   
#   OP_Lung_Liver_norm[i] <- (OP_48$X14+OP_48$X15)/(OP_48$X8+OP_48$X9)
#   
#   OP_Other_Liver_norm[i] <- OP_48$X17/(OP_48$X8+OP_48$X9)
#   
# }
# 
# code_time <- proc.time() - code_time_start
# 
# OPT2 <- data.frame(run=seq(N),
#                   kve_Heart=kve_Heart_keep,
#                   kve_Kidney=kve_Kidney_keep,
#                   kve_Liver=kve_Liver_keep,
#                   kve_Lung=kve_Lung_keep,
#                   kve_Spleen=kve_Spleen_keep,
#                   kve_Other=kve_Other_keep,
#                   kexc_Kidney=kexc_Kidney_keep,
#                   kexc_Liver=kexc_Liver_keep,
#                   kvp_Liver=kvp_Liver_keep,
#                   kvp_Spleen=kvp_Spleen_keep,
#                   kvp_Lung=kvp_Lung_keep,
#                   MPmax_Liver=MPmax_Liver_keep,
#                   MPmax_Spleen=MPmax_Spleen_keep,
#                   MPmax_Lung=MPmax_Lung_keep,
#                   kdeg_Liver=kdeg_Liver_keep,
#                   kdeg_Spleen=kdeg_Spleen_keep,
#                   kdeg_Lung=kdeg_Lung_keep,
#                   
#                   blood_4 = OP_blood_4_keep,
#                   blood_48 = OP_blood_48_keep,
#                   
#                   Heart_norm=OP_Heart_Liver_norm,
#                   Kidney_norm=OP_Kidney_Liver_norm,
#                   Spleen_norm=OP_Spleen_Liver_norm,
#                   Lung_norm=OP_Lung_Liver_norm,
#                   Other_norm=OP_Other_Liver_norm)
# 
# 
# Z_blood_4 <- (mean(OPT2$blood_4)-blood_summary$mean[blood_summary$time==4])/
#                     sqrt(var(OPT2$blood_4)/n+(blood_summary$sd[blood_summary$time==4])^2/n)
# 
# Z_blood_48 <- (mean(OPT2$blood_48)-blood_summary$mean[blood_summary$time==48])/
#   sqrt(var(OPT2$blood_48)/n+(blood_summary$sd[blood_summary$time==48])^2/n)
# 
# Z_Heart <- (mean(OPT2$Heart_norm)-organ_FACS_summary$mean[organ_FACS_summary$Organ=='Heart'])/
#   sqrt(var(OPT2$Heart_norm)/n+(organ_FACS_summary$sd[organ_FACS_summary$Organ=='Heart'])^2/n)
# 
# Z_Kidney <- (mean(OPT2$Kidney_norm)-organ_FACS_summary$mean[organ_FACS_summary$Organ=='Kidneys'])/
#   sqrt(var(OPT2$Kidney_norm)/n+(organ_FACS_summary$sd[organ_FACS_summary$Organ=='Kidneys'])^2/n)
# 
# Z_Lung <- (mean(OPT2$Lung_norm)-organ_FACS_summary$mean[organ_FACS_summary$Organ=='Lungs'])/
#   sqrt(var(OPT2$Lung_norm)/n+(organ_FACS_summary$sd[organ_FACS_summary$Organ=='Lungs'])^2/n)
# 
# Z_Spleen <- (mean(OPT2$Spleen_norm)-organ_FACS_summary$mean[organ_FACS_summary$Organ=='Spleen'])/
#   sqrt(var(OPT2$Spleen_norm)/n+(organ_FACS_summary$sd[organ_FACS_summary$Organ=='Spleen'])^2/n)
# 
# Z_Other <- (mean(OPT2$Other_norm)-organ_FACS_summary$mean[organ_FACS_summary$Organ=='Bone'])/
#   sqrt(var(OPT2$Other_norm)/n+(organ_FACS_summary$sd[organ_FACS_summary$Organ=='Bone'])^2/n)
# 
# Z_scores <- c(Z_blood_4, Z_blood_48, Z_Heart,Z_Kidney,Z_Lung,Z_Spleen,Z_Other)
# 
# parms_keep <- list(kve_Heart_mean = mean(OPT_cond$kve_Heart),
#                     kve_Heart_var = var(OPT_cond$kve_Heart),
#                     kve_Kidney_mean = mean(OPT_cond$kve_Kidney),
#                     kve_Kidney_var = var(OPT_cond$kve_Kidney),
#                     kve_Lung_mean = mean(OPT_cond$kve_Lung),
#                     kve_Lung_var = var(OPT_cond$kve_Lung),
#                     kve_Liver_mean = mean(OPT_cond$kve_Liver),
#                     kve_Liver_var = var(OPT_cond$kve_Liver),
#                     kve_Other_mean = mean(OPT_cond$kve_Other),
#                     kve_Other_var = var(OPT_cond$kve_Other),
#                     kve_Spleen_mean = mean(OPT_cond$kve_Spleen),
#                     kve_Spleen_var = var(OPT_cond$kve_Spleen),
#                    
#                    kvp_Lung_mean = mean(OPT_cond$kvp_Lung),
#                    kvp_Lung_var = var(OPT_cond$kvp_Lung),
#                    kvp_Liver_mean = mean(OPT_cond$kvp_Liver),
#                    kvp_Liver_var = var(OPT_cond$kvp_Liver),
#                    kvp_Spleen_mean = mean(OPT_cond$kvp_Spleen),
#                    kvp_Spleen_var = var(OPT_cond$kvp_Spleen),
#                    
#                    MPmax_Lung_mean = mean(OPT_cond$MPmax_Lung),
#                    MPmax_Lung_var = var(OPT_cond$MPmax_Lung),
#                    MPmax_Liver_mean = mean(OPT_cond$MPmax_Liver),
#                    MPmax_Liver_var = var(OPT_cond$MPmax_Liver),
#                    MPmax_Spleen_mean = mean(OPT_cond$MPmax_Spleen),
#                    MPmax_Spleen_var = var(OPT_cond$MPmax_Spleen),
#                    
#                    kdeg_Lung_mean = mean(OPT_cond$kdeg_Lung),
#                    kdeg_Lung_var = var(OPT_cond$kdeg_Lung),
#                    kdeg_Liver_mean = mean(OPT_cond$kdeg_Liver),
#                    kdeg_Liver_var = var(OPT_cond$kdeg_Liver),
#                    kdeg_Spleen_mean = mean(OPT_cond$kdeg_Spleen),
#                    kdeg_Spleen_var = var(OPT_cond$kdeg_Spleen),
# 
#                    kexc_Liver_mean = mean(OPT_cond$kexc_Liver),
#                    kexc_Liver_var = var(OPT_cond$kexc_Liver),
#                    kexc_Kidney_mean = mean(OPT_cond$kexc_Kidney),
#                    kexc_Kidney_var = var(OPT_cond$kexc_Kidney))
# 
# saveRDS(parms_keep,"parms_keep_20221209.rds")
# saveRDS(OP_keep,"OP_keep_20221209.rds")


