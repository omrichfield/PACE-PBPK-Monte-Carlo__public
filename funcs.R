#Functions for PACE NP PBPK Model
#Owen Richfield, Oct. 26 2023

#Functional form from (Lin, 2016)
#Citation: Zhoumeng Lin, Nancy A. Monteiro-Riviere & Jim E. Riviere (2016) A
#physiologically based pharmacokinetic model for polyethylene glycol-coated gold nanoparticles of
#different sizes in adult mice, Nanotoxicology, 10:2, 162-172, DOI: 10.3109/17435390.2015.1027314

NP_PBPK3 <- function(t,y,parms){
  
  #Assign variables
  M_Art <- y[1]
  M_Ven <- y[2]
  
  M_vasc_Liver <- y[3]
  M_extra_Liver <- y[4]
  M_phago_Liver <- y[5]
  
  M_vasc_Spleen <- y[6]
  M_extra_Spleen <- y[7]
  M_phago_Spleen <- y[8]
  
  M_vasc_Kidneys <- y[9]
  M_extra_Kidneys <- y[10]
  M_phago_Kidneys <- y[11]
  
  M_vasc_Lungs <- y[12]
  M_extra_Lungs <- y[13]
  M_phago_Lungs <- y[14]
  
  M_vasc_Brain <- y[15]
  M_extra_Brain <- y[16]
  
  M_vasc_Heart <- y[17]
  M_extra_Heart <- y[18]
  
  M_vasc_Bone <- y[19]
  M_extra_Bone <- y[20]
  
  M_vasc_Rest <- y[21]
  M_extra_Rest <- y[22]
  
  #Calculate concentrations
  
  C_Art <- M_Art/(VArtC*BW)
  C_Ven <- M_Ven/(VVenC*BW)
  
  C_vasc_Liver <- M_vasc_Liver/(VLC*BVL*BW)
  C_extra_Liver <- M_extra_Liver/(VLC*(1-BVL)*BW)
  
  C_vasc_Spleen <- M_vasc_Spleen/(VSC*BVS*BW)  
  C_extra_Spleen <- M_extra_Spleen/(VSC*(1-BVS)*BW) 
  
  C_vasc_Kidneys <- M_vasc_Kidneys/(VKC*BVK*BW)
  C_extra_Kidneys <- M_extra_Kidneys/(VKC*(1-BVK)*BW) 
  
  C_vasc_Lungs <- M_vasc_Lungs/(VLuC*BVLu*BW) 
  C_extra_Lungs <- M_extra_Lungs/(VLuC*(1-BVLu)*BW) 
  
  C_vasc_Heart <- M_vasc_Heart/(VHC*BVH*BW)
  C_extra_Heart <- M_extra_Heart/(VHC*(1-BVH)*BW) 
  
  C_vasc_Bone <- M_vasc_Bone/(VBoneC*BVBone*BW)
  C_extra_Bone <- M_extra_Bone/(VBoneC*(1-BVBone)*BW) 
  
  C_vasc_Brain <- M_vasc_Brain/(VBRC*BVBR*BW)
  C_extra_Brain <- M_extra_Brain/(VBRC*(1-BVBR)*BW) 
  
  C_vasc_Rest <- M_vasc_Rest/(VrestC*BVrest*BW) 
  C_extra_Rest <- M_extra_Rest/(VrestC*(1-BVrest)*BW) 
  
  #Arterial
  
  dM_Art <- QC*C_vasc_Heart - QC*C_Art
  
  #Venous
  
  if (t <= dose_time){
    R_dose <- dose_rate
  }else{
    R_dose <- 0
  }
  
  dM_Ven <- QC*QLC*C_vasc_Liver + QC*QBRC*C_vasc_Brain + QC*QBoneC*C_vasc_Bone +
    QC*QKC*C_vasc_Kidneys + QC*QrestC*C_vasc_Rest - QC*C_Ven + R_dose
  
  #Liver
  
  K_up_Liver <- (Kmax_t_Liver*t^nt)/(K50_t_Liver^nt + t^nt)
  
  R_up_Liver <- K_up_Liver*M_vasc_Liver
  R_release_Liver <- Krelease_t_Liver*M_phago_Liver
  R_bile <- Kbile*C_extra_Liver
  R_deg_Liver <- Kdeg_t_Liver*M_phago_Liver
  
  dM_phago_Liver <- R_up_Liver - R_release_Liver - R_deg_Liver
  
  dM_vasc_Liver <- QC*QLC*(C_Art - C_vasc_Liver) + QC*QSC*C_vasc_Spleen - PACt_Liver*C_vasc_Liver + 
    (PACt_Liver*C_extra_Liver)/Pt_Liver + R_release_Liver - R_up_Liver 
  
  dM_extra_Liver <- PACt_Liver*C_vasc_Liver - (PACt_Liver*C_extra_Liver)/Pt_Liver - R_bile
  
  #Spleen
  
  K_up_Spleen <- (Kmax_t_Spleen*t^nt)/(K50_t_Spleen^nt + t^nt)
  
  R_up_Spleen <- K_up_Spleen*M_vasc_Spleen
  R_release_Spleen <- Krelease_t_Spleen*M_phago_Spleen
  R_deg_Spleen <- Kdeg_t_Spleen*M_phago_Spleen
  
  dM_phago_Spleen <- R_up_Spleen - R_release_Spleen - R_deg_Spleen 
  
  dM_vasc_Spleen <- QC*QSC*(C_Art - C_vasc_Spleen) - PACt_Spleen*C_vasc_Spleen + 
    (PACt_Spleen*C_extra_Spleen)/Pt_Spleen + R_release_Spleen - R_up_Spleen
  
  dM_extra_Spleen <- PACt_Spleen*C_vasc_Spleen - (PACt_Spleen*C_extra_Spleen)/Pt_Spleen
  
  #Kidneys
  
  K_up_Kidneys <- (Kmax_t_Kidneys*t^nt)/(K50_t_Kidneys^nt + t^nt)
  
  R_up_Kidneys <- K_up_Kidneys*M_vasc_Kidneys
  R_release_Kidneys <- Krelease_t_Kidneys*M_phago_Kidneys
  R_urine <- Kurine*C_extra_Kidneys
  R_deg_Kidneys <- Kdeg_t_Kidneys*M_phago_Kidneys
  
  dM_phago_Kidneys <- R_up_Kidneys - R_release_Kidneys - R_deg_Kidneys
  
  dM_vasc_Kidneys <- QC*QKC*(C_Art - C_vasc_Kidneys) - PACt_Kidneys*C_vasc_Kidneys + 
    (PACt_Kidneys*C_extra_Kidneys)/Pt_Kidneys + R_release_Kidneys - R_up_Kidneys 
  
  dM_extra_Kidneys <- PACt_Kidneys*C_vasc_Kidneys - 
    (PACt_Kidneys*C_extra_Kidneys)/Pt_Kidneys - R_urine
  
  #Lungs
  
  K_up_Lungs <- (Kmax_t_Lungs*t^nt)/(K50_t_Lungs^nt + t^nt)
  
  R_up_Lungs <- K_up_Lungs*M_vasc_Lungs
  R_release_Lungs <- Krelease_t_Lungs*M_phago_Lungs
  R_deg_Lungs <- Kdeg_t_Lungs*M_phago_Lungs
  
  dM_phago_Lungs <- R_up_Lungs - R_release_Lungs 
  
  dM_vasc_Lungs <- QC*QLuC*(C_vasc_Heart - C_vasc_Lungs) - PACt_Lungs*C_vasc_Lungs + 
    (PACt_Lungs*C_extra_Lungs)/Pt_Lungs + R_release_Lungs - R_up_Lungs
  
  dM_extra_Lungs <- PACt_Lungs*C_vasc_Lungs - (PACt_Lungs*C_extra_Lungs)/Pt_Lungs
  
  #Heart
  
  dM_vasc_Heart <- QC*QHC*(C_Ven + C_vasc_Lungs - 2*C_vasc_Heart) - PACt_Heart*C_vasc_Heart + 
    (PACt_Heart*C_extra_Heart)/Pt_Heart 
  
  dM_extra_Heart <- PACt_Heart*C_vasc_Heart - (PACt_Heart*C_extra_Heart)/Pt_Heart
  
  #Brain
  
  dM_vasc_Brain <- QC*QBRC*(C_Art - C_vasc_Brain) - PACt_Brain*C_vasc_Brain + 
    (PACt_Brain*C_extra_Brain)/Pt_Brain 
  
  dM_extra_Brain <- PACt_Brain*C_vasc_Brain - (PACt_Brain*C_extra_Brain)/Pt_Brain
  
  #Bone
  
  dM_vasc_Bone <- QC*QBoneC*(C_Art - C_vasc_Bone) - PACt_Bone*C_vasc_Bone + 
    (PACt_Bone*C_extra_Bone)/Pt_Bone 
  
  dM_extra_Bone <- PACt_Bone*C_vasc_Bone - (PACt_Bone*C_extra_Bone)/Pt_Bone
  
  #Rest of body
  
  dM_vasc_Rest <- QC*QrestC*(C_Art - C_vasc_Rest) - PACt_Rest*C_vasc_Rest + 
    (PACt_Rest*C_extra_Rest)/Pt_Rest 
  
  dM_extra_Rest <- PACt_Rest*C_vasc_Rest - (PACt_Rest*C_extra_Rest)/Pt_Rest
  
  
  
  
  output <- list(c(dM_Art,
                   dM_Ven,
                   
                   dM_vasc_Liver,
                   dM_extra_Liver,
                   dM_phago_Liver,
                   
                   dM_vasc_Spleen,
                   dM_extra_Spleen,
                   dM_phago_Spleen,
                   
                   dM_vasc_Kidneys,
                   dM_extra_Kidneys,
                   dM_phago_Kidneys,
                   
                   dM_vasc_Lungs,
                   dM_extra_Lungs,
                   dM_phago_Lungs,
                   
                   dM_vasc_Brain,
                   dM_extra_Brain,
                   
                   dM_vasc_Heart,
                   dM_extra_Heart,
                   
                   dM_vasc_Bone,
                   dM_extra_Bone,
                   
                   dM_vasc_Rest,
                   dM_extra_Rest))
  
  return(output)
  
}

NP_PBPK2 <- function(t,y,parms){
  
  #Assign variables
  M_Art <- y[1]
  M_Ven <- y[2]
  
  M_vasc_Liver <- y[3]
  M_extra_Liver <- y[4]
  M_phago_Liver <- y[5]
  
  M_vasc_Spleen <- y[6]
  M_extra_Spleen <- y[7]
  M_phago_Spleen <- y[8]
  
  M_vasc_Kidneys <- y[9]
  M_extra_Kidneys <- y[10]
  M_phago_Kidneys <- y[11]
  
  M_vasc_Lungs <- y[12]
  M_extra_Lungs <- y[13]
  M_phago_Lungs <- y[14]
  
  M_vasc_Brain <- y[15]
  M_extra_Brain <- y[16]
  
  M_vasc_Heart <- y[17]
  M_extra_Heart <- y[18]
  
  M_vasc_Bone <- y[19]
  M_extra_Bone <- y[20]
  
  M_vasc_Rest <- y[21]
  M_extra_Rest <- y[22]
  
  #Calculate concentrations
  
  C_Art <- M_Art/(VArtC*BW)
  C_Ven <- M_Ven/(VVenC*BW)
  
  C_vasc_Liver <- M_vasc_Liver/(VLC*BVL*BW)
  C_extra_Liver <- M_extra_Liver/(VLC*(1-BVL)*BW)
  
  C_vasc_Spleen <- M_vasc_Spleen/(VSC*BVS*BW)  
  C_extra_Spleen <- M_extra_Spleen/(VSC*(1-BVS)*BW) 
  
  C_vasc_Kidneys <- M_vasc_Kidneys/(VKC*BVK*BW)
  C_extra_Kidneys <- M_extra_Kidneys/(VKC*(1-BVK)*BW) 
  
  C_vasc_Lungs <- M_vasc_Lungs/(VLuC*BVLu*BW) 
  C_extra_Lungs <- M_extra_Lungs/(VLuC*(1-BVLu)*BW) 
  
  C_vasc_Heart <- M_vasc_Heart/(VHC*BVH*BW)
  C_extra_Heart <- M_extra_Heart/(VHC*(1-BVH)*BW) 
  
  C_vasc_Bone <- M_vasc_Bone/(VBoneC*BVBone*BW)
  C_extra_Bone <- M_extra_Bone/(VBoneC*(1-BVBone)*BW) 
  
  C_vasc_Brain <- M_vasc_Brain/(VBRC*BVBR*BW)
  C_extra_Brain <- M_extra_Brain/(VBRC*(1-BVBR)*BW) 
  
  C_vasc_Rest <- M_vasc_Rest/(VrestC*BVrest*BW) 
  C_extra_Rest <- M_extra_Rest/(VrestC*(1-BVrest)*BW) 
  
  #Arterial
  
  dM_Art <- QC*C_vasc_Heart - QC*C_Art
  
  #Venous
  
  if (t <= dose_time){
    R_dose <- dose_rate
  }else{
    R_dose <- 0
  }
  
  dM_Ven <- QC*QLC*C_vasc_Liver + QC*QBRC*C_vasc_Brain + QC*QBoneC*C_vasc_Bone +
    QC*QKC*C_vasc_Kidneys + QC*QrestC*C_vasc_Rest - QC*C_Ven + R_dose
  
  #Liver
  
  K_up_Liver <- (Kmax_t_Liver*t^nt)/(K50_t_Liver^nt + t^nt)
  
  R_up_Liver <- K_up_Liver*M_vasc_Liver
  R_release_Liver <- Krelease_t_Liver*M_phago_Liver
  R_bile <- Kbile*C_extra_Liver
  
  dM_phago_Liver <- R_up_Liver - R_release_Liver 
  
  dM_vasc_Liver <- QC*QLC*(C_Art - C_vasc_Liver) + QC*QSC*C_vasc_Spleen - PACt_Liver*C_vasc_Liver + 
    (PACt_Liver*C_extra_Liver)/Pt_Liver + R_release_Liver - R_up_Liver 
  
  dM_extra_Liver <- PACt_Liver*C_vasc_Liver - (PACt_Liver*C_extra_Liver)/Pt_Liver - R_bile
  
  #Spleen
  
  K_up_Spleen <- (Kmax_t_Spleen*t^nt)/(K50_t_Spleen^nt + t^nt)
  
  R_up_Spleen <- K_up_Spleen*M_vasc_Spleen
  R_release_Spleen <- Krelease_t_Spleen*M_phago_Spleen
  
  dM_phago_Spleen <- R_up_Spleen - R_release_Spleen 
  
  dM_vasc_Spleen <- QC*QSC*(C_Art - C_vasc_Spleen) - PACt_Spleen*C_vasc_Spleen + 
    (PACt_Spleen*C_extra_Spleen)/Pt_Spleen + R_release_Spleen - R_up_Spleen
  
  dM_extra_Spleen <- PACt_Spleen*C_vasc_Spleen - (PACt_Spleen*C_extra_Spleen)/Pt_Spleen
  
  #Kidneys
  
  K_up_Kidneys <- (Kmax_t_Kidneys*t^nt)/(K50_t_Kidneys^nt + t^nt)
  
  R_up_Kidneys <- K_up_Kidneys*M_vasc_Kidneys
  R_release_Kidneys <- Krelease_t_Kidneys*M_phago_Kidneys
  R_urine <- Kurine*C_extra_Kidneys
  
  dM_phago_Kidneys <- R_up_Kidneys - R_release_Kidneys 
  
  dM_vasc_Kidneys <- QC*QKC*(C_Art - C_vasc_Kidneys) - PACt_Kidneys*C_vasc_Kidneys + 
                        (PACt_Kidneys*C_extra_Kidneys)/Pt_Kidneys + R_release_Kidneys - R_up_Kidneys 
  
  dM_extra_Kidneys <- PACt_Kidneys*C_vasc_Kidneys - 
                            (PACt_Kidneys*C_extra_Kidneys)/Pt_Kidneys - R_urine
  
  #Lungs
  
  K_up_Lungs <- (Kmax_t_Lungs*t^nt)/(K50_t_Lungs^nt + t^nt)
  
  R_up_Lungs <- K_up_Lungs*M_vasc_Lungs
  R_release_Lungs <- Krelease_t_Lungs*M_phago_Lungs
  
  dM_phago_Lungs <- R_up_Lungs - R_release_Lungs 
  
  dM_vasc_Lungs <- QC*QLuC*(C_vasc_Heart - C_vasc_Lungs) - PACt_Lungs*C_vasc_Lungs + 
    (PACt_Lungs*C_extra_Lungs)/Pt_Lungs + R_release_Lungs - R_up_Lungs
  
  dM_extra_Lungs <- PACt_Lungs*C_vasc_Lungs - (PACt_Lungs*C_extra_Lungs)/Pt_Lungs
  
  #Heart
  
  dM_vasc_Heart <- QC*QHC*(C_Ven + C_vasc_Lungs - 2*C_vasc_Heart) - PACt_Heart*C_vasc_Heart + 
    (PACt_Heart*C_extra_Heart)/Pt_Heart 
  
  dM_extra_Heart <- PACt_Heart*C_vasc_Heart - (PACt_Heart*C_extra_Heart)/Pt_Heart
  
  #Brain
  
  dM_vasc_Brain <- QC*QBRC*(C_Art - C_vasc_Brain) - PACt_Brain*C_vasc_Brain + 
    (PACt_Brain*C_extra_Brain)/Pt_Brain 
  
  dM_extra_Brain <- PACt_Brain*C_vasc_Brain - (PACt_Brain*C_extra_Brain)/Pt_Brain
  
  #Bone
  
  dM_vasc_Bone <- QC*QBoneC*(C_Art - C_vasc_Bone) - PACt_Bone*C_vasc_Bone + 
    (PACt_Bone*C_extra_Bone)/Pt_Bone 
  
  dM_extra_Bone <- PACt_Bone*C_vasc_Bone - (PACt_Bone*C_extra_Bone)/Pt_Bone
  
  #Rest of body
  
  dM_vasc_Rest <- QC*QrestC*(C_Art - C_vasc_Rest) - PACt_Rest*C_vasc_Rest + 
    (PACt_Rest*C_extra_Rest)/Pt_Rest 
  
  dM_extra_Rest <- PACt_Rest*C_vasc_Rest - (PACt_Rest*C_extra_Rest)/Pt_Rest
  
  
  
  
  output <- list(c(dM_Art,
                   dM_Ven,
                   
                   dM_vasc_Liver,
                   dM_extra_Liver,
                   dM_phago_Liver,
                   
                   dM_vasc_Spleen,
                   dM_extra_Spleen,
                   dM_phago_Spleen,
                   
                   dM_vasc_Kidneys,
                   dM_extra_Kidneys,
                   dM_phago_Kidneys,
                   
                   dM_vasc_Lungs,
                   dM_extra_Lungs,
                   dM_phago_Lungs,
                   
                   dM_vasc_Brain,
                   dM_extra_Brain,
                   
                   dM_vasc_Heart,
                   dM_extra_Heart,
                   
                   dM_vasc_Bone,
                   dM_extra_Bone,
                   
                   dM_vasc_Rest,
                   dM_extra_Rest))
  
  return(output)
  
}

NP_PBPK <- function(t,y,parms){
  
  #Assign variables
  M_Art <- y[1]
  M_Ven <- y[2]
  
  M_vasc_Liver <- y[3]
  M_extra_Liver <- y[4]
  M_phago_Liver <- y[5]

  M_vasc_Spleen <- y[6]
  M_extra_Spleen <- y[7]
  M_phago_Spleen <- y[8]
  
  M_vasc_Kidneys <- y[9]
  M_extra_Kidneys <- y[10]
  M_phago_Kidneys <- y[11]
  
  M_vasc_Lungs <- y[12]
  M_extra_Lungs <- y[13]
  M_phago_Lungs <- y[14]
  
  M_vasc_Brain <- y[15]
  M_extra_Brain <- y[16]
  
  M_vasc_Rest <- y[17]
  M_extra_Rest <- y[18]
  
  #Calculate concentrations
  
  C_Art <- M_Art/(VArtC*BW)
  C_Ven <- M_Ven/(VVenC*BW)
  
  C_vasc_Liver <- M_vasc_Liver/(VLC*BVL*BW)
  C_extra_Liver <- M_extra_Liver/(VLC*(1-BVL)*BW)

  C_vasc_Spleen <- M_vasc_Spleen/(VSC*BVS*BW)  
  C_extra_Spleen <- M_extra_Spleen/(VSC*(1-BVS)*BW) 

  C_vasc_Kidneys <- M_vasc_Kidneys/(VKC*BVK*BW)
  C_extra_Kidneys <- M_extra_Kidneys/(VKC*(1-BVK)*BW) 

  C_vasc_Lungs <- M_vasc_Lungs/(VLuC*BVLu*BW) 
  C_extra_Lungs <- M_extra_Lungs/(VLuC*(1-BVLu)*BW) 

  C_vasc_Brain <- M_vasc_Brain/(VBRC*BVBR*BW)
  C_extra_Brain <- M_extra_Brain/(VBRC*(1-BVBR)*BW) 
  
  C_vasc_Rest <- M_vasc_Rest/(VrestC*BVrest*BW) 
  C_extra_Rest <- M_extra_Rest/(VrestC*(1-BVrest)*BW) 

  #Arterial
  
  dM_Art <- QC*C_vasc_Lungs - QC*C_Art
  
  #Venous
  
  if (t > 0.03 & t <= dose_time + 0.03){
    R_dose <- dose_rate*dtphr
  }else{
    R_dose <- 0
  }
  
  dM_Ven <- QC*QLC*C_vasc_Liver + QC*QBRC*C_vasc_Brain + 
                  QC*QKC*C_vasc_Kidneys + QC*QrestC*C_vasc_Rest - QC*C_Ven + R_dose
  
  #Liver

  K_up_Liver <- (Kmax_t_Liver*t^nt)/(K50_t_Liver^nt + t^nt)

  R_up_Liver <- K_up_Liver*M_vasc_Liver
  R_release_Liver <- Krelease_t_Liver*M_phago_Liver
  R_bile <- Kbile*C_extra_Liver
  
  dM_phago_Liver <- R_up_Liver - R_release_Liver 
    
  dM_vasc_Liver <- QC*QLC*(C_Art - C_vasc_Liver)+ QC*QSC*C_vasc_Spleen - PACt_Liver*C_vasc_Liver + 
                    (PACt_Liver*C_extra_Liver)/Pt_Liver + R_release_Liver - R_up_Liver - R_bile
  
  dM_extra_Liver <- PACt_Liver*C_vasc_Liver - (PACt_Liver*C_extra_Liver)/Pt_Liver
  
  #Spleen
  
  K_up_Spleen <- (Kmax_t_Spleen*t^nt)/(K50_t_Spleen^nt + t^nt)
  
  R_up_Spleen <- K_up_Spleen*M_vasc_Spleen
  R_release_Spleen <- Krelease_t_Spleen*M_phago_Spleen
  
  dM_phago_Spleen <- R_up_Spleen - R_release_Spleen 
  
  dM_vasc_Spleen <- QC*QSC*(C_Art - C_vasc_Spleen) - PACt_Spleen*C_vasc_Spleen + 
    (PACt_Spleen*C_extra_Spleen)/Pt_Spleen + R_release_Spleen - R_up_Spleen
  
  dM_extra_Spleen <- PACt_Spleen*C_vasc_Spleen - (PACt_Spleen*C_extra_Spleen)/Pt_Spleen

  #Kidneys
  
  K_up_Kidneys <- (Kmax_t_Kidneys*t^nt)/(K50_t_Kidneys^nt + t^nt)
  
  R_up_Kidneys <- K_up_Kidneys*M_vasc_Kidneys
  R_release_Kidneys <- Krelease_t_Kidneys*M_phago_Kidneys
  R_urine <- Kurine*C_extra_Kidneys
  
  dM_phago_Kidneys <- R_up_Kidneys - R_release_Kidneys 
  
  dM_vasc_Kidneys <- QC*QKC*(C_Art - C_vasc_Kidneys) - PACt_Kidneys*C_vasc_Kidneys + 
    (PACt_Kidneys*C_extra_Kidneys)/Pt_Kidneys + R_release_Kidneys - R_up_Kidneys - R_urine
  
  dM_extra_Kidneys <- PACt_Kidneys*C_vasc_Kidneys - (PACt_Kidneys*C_extra_Kidneys)/Pt_Kidneys
  
  #Lungs
  
  K_up_Lungs <- (Kmax_t_Lungs*t^nt)/(K50_t_Lungs^nt + t^nt)
  
  R_up_Lungs <- K_up_Lungs*M_vasc_Lungs
  R_release_Lungs <- Krelease_t_Lungs*M_phago_Lungs
  
  dM_phago_Lungs <- R_up_Lungs - R_release_Lungs 
  
  dM_vasc_Lungs <- QC*QLuC*(C_Ven - C_vasc_Lungs) - PACt_Lungs*C_vasc_Lungs + 
    (PACt_Lungs*C_extra_Lungs)/Pt_Lungs + R_release_Lungs - R_up_Lungs
  
  dM_extra_Lungs <- PACt_Lungs*C_vasc_Lungs - (PACt_Lungs*C_extra_Lungs)/Pt_Lungs
  
  #Brain
  
  dM_vasc_Brain <- QC*QBRC*(C_Art - C_vasc_Brain) - PACt_Brain*C_vasc_Brain + 
    (PACt_Brain*C_extra_Brain)/Pt_Brain 
  
  dM_extra_Brain <- PACt_Brain*C_vasc_Brain - (PACt_Brain*C_extra_Brain)/Pt_Brain
  
  #Rest of body
  
  dM_vasc_Rest <- QC*QrestC*(C_Art - C_vasc_Rest) - PACt_Rest*C_vasc_Rest + 
    (PACt_Rest*C_extra_Rest)/Pt_Rest 
  
  dM_extra_Rest <- PACt_Rest*C_vasc_Rest - (PACt_Rest*C_extra_Rest)/Pt_Rest
  

  
  
  output <- list(c(dM_Art,
                   dM_Ven,
                   
                   dM_vasc_Liver,
                   dM_extra_Liver,
                   dM_phago_Liver,
                   
                   dM_vasc_Spleen,
                   dM_extra_Spleen,
                   dM_phago_Spleen,
                   
                   dM_vasc_Kidneys,
                   dM_extra_Kidneys,
                   dM_phago_Kidneys,
                   
                   dM_vasc_Lungs,
                   dM_extra_Lungs,
                   dM_phago_Lungs,
                   
                   dM_vasc_Brain,
                   dM_extra_Brain,
                   
                   dM_vasc_Rest,
                   dM_extra_Rest))
  
  return(output)
  
}