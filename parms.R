#Parameters for PACE NP PBPK Model
#Owen Richfield, Oct. 26 2023

#Most parameters initially derived from (Lin, 2016) unless otherwise noted
#Citation: Zhoumeng Lin, Nancy A. Monteiro-Riviere & Jim E. Riviere (2016) A
#physiologically based pharmacokinetic model for polyethylene glycol-coated gold nanoparticles of
#different sizes in adult mice, Nanotoxicology, 10:2, 162-172, DOI: 10.3109/17435390.2015.1027314

#All NP PK parameters correspond to the 100nm diameter NP case in (Lin, 2016)

#All parameters pertain to organs with the order: 
#Liver (L)
#Spleen (S)
#Kidneys (K)
#Lungs (Lu)
#Brain (BR)
#Rest of Body (Rest)

#Time parameters
time_dum <- 2
Ttot <- 48                              #hrs
dtphr <- 2.5e2
Nt <- Ttot*dtphr       
dt <- Ttot/Nt                           #hrs
times <- seq(from=1, to=Nt, by=1)*dt    #hrs

#Data analysis parameters
alpha_scale <- 1e-2
alpha_parm <- 3e-2

#Physiological Parameters - units: L, kg -> g, h
BW <- 0.02 #body weight (kg)
QC <- 16.5*BW #cardiac output (L/h)

#Dose parameters
dose_rel <- 0.5#c(0.1,0.5,2,2.5,3.5) #mg
dose <- dose_rel             #kg
dose_time <- 20/3600         #hr
dose_rate <- dose/dose_time/0.14*0.5  #kg/hr
dose_rate_vec <- dose/dose_time/0.14*0.5  #kg/hr

#Blood flow to organ (fraction of CO)
QLC <- 0.161
QSC <- 0.011
QKC <- 0.091
QLuC <- 1
QHC <- 1
QBRC <- 0.033
QBoneC <- 0.072 #relative to QLC, based on Davies
QrestC <- 0.704 - QBoneC

#Organ volumes (fraction of BW)
VLC <- 0.055
VSC <- 0.005
VKC <- 0.017
VLuC <- 0.007
VBRC <- 0.017
VHC <- 0.005
VBoneC <- 0.0175
VrestC <- 0.85-VHC-VBoneC
VBloodC <- 0.049
VPlasmaC <- 0.029
VVenC <- 0.8*VBloodC
VArtC <- 0.2*VBloodC

#Volume fraction of blood in organs
BVL <- 0.31
BVS <- 0.17
BVK <- 0.24
BVLu <- 0.50
BVBone <- 0.275
BVH <- 0.24 #assumed
BVBR <- 0.03
BVrest <- 0.04

#################################################

#Pharmacokinetic parameters

#Tissue:plasma distribution coefficient Pt (unitless)
Pt_Liver <- 0.08
Pt_Spleen <- 0.15
Pt_Kidneys <- 0.15
Pt_Heart <- 0.15
Pt_Bone <- 0.15
Pt_Lungs <- 0.15
Pt_Brain <- 0.15
Pt_Rest <- 0.15

#Permeability coefficient PACt (unitless)
PACt_Liver_Au <- 0.001

PACt_Spleen_Au <- 0.001
PACt_Kidneys_Au <- 0.001 
PACt_Lungs_Au <- 0.001
PACt_Heart_Au <- 0.001  
PACt_Bone_Au <- 0.001
PACt_Brain_Au <- 0.000001
PACt_Rest_Au <- 0.000001

#Maximum uptake rate constant of phagocytic cells (1/h)
Kmax_t_Liver_Au <- 4
Kmax_t_Spleen_Au <- 10
Kmax_t_Kidneys_Au <- 0.1 
Kmax_t_Lungs_Au <- 0.1

#Time for reaching half maximum uptake rate (h)
K50_t_Liver_Au <- 24
K50_t_Spleen_Au <- 24
K50_t_Kidneys_Au <- 24
K50_t_Lungs_Au <- 24

#Hill coefficient (unitless)
nt_Au <- 0.1

#Release rate constant of phagocytic cells
Krelease_t_Liver_Au <- 0.0075
Krelease_t_Spleen_Au <- 0.003
Krelease_t_Kidneys_Au <- 0.01
Krelease_t_Lungs_Au <- 0.005

#Elimination rate constants
Kbile_Au <- 0.0012
Kurine_Au <- 0.00012

PACt_Heart_Mod <- 1/10
PACt_Bone_Mod <- 1/10#1.5e-2
PACt_Spleen_Mod <- 1#3e-2
PACt_Lungs_Mod <- 1/100#1.25e-4
PACt_Rest_Mod <- 1
PACt_Kidneys_Mod <- 1e-6
Kmax_t_Kidneys_Mod <- 1#3.5e0
Kmax_t_Lungs_Mod <- 1/100#3.5e0

CNH <- 3.3e4 #https://journals.physiology.org/doi/full/10.1152/ajpheart.00137.2019#:~:text=4%2C%20the%20absolute%20number%20of,evaluation%20of%20flow%20cytometry%20data.
CNL <- 264.8e6 #https://www.cambridge.org/core/services/aop-cambridge-core/content/view/FDAE3F61AD2CCC363AB88B677BE3708C/S0016672300018061a.pdf/cell_numbers_and_cell_sizes_in_organs_of_mice_selected_for_large_and_small_body_size.pdf
CNLu <- 86.1e6 #^
CNS <- 193.9e6 #^
CNK <- 142.0e6 #^
CNBone <- 300e6 #Nombela




