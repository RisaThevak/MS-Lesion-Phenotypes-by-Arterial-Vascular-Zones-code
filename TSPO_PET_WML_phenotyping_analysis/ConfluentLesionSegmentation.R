library(neurobase)
library(extrantsr)
library(ANTsRCore)
library(ANTsR)
library(parallel)
library(pdist)
library(plyr)
library(mclust)

#Loading functions for J. Dworkin's lesiontools toolbox 
source('/export02/data/risa/lesiontools-master/R/gradient.R')
source('/export02/data/risa/lesiontools-master/R/hessian.R')
source('/export02/data/risa/lesiontools-master/R/lesioncenters.R')
source('/export02/data/risa/lesiontools-master/R/lesioncenters_c3d.R')
source('/export02/data/risa/lesiontools-master/R/lesionclusters.R')
source('/export02/data/risa/lesiontools-master/R/labelreg.R')
source('/export02/data/risa/lesiontools-master/R/frangi.R')
source('/export02/data/risa/lesiontools-master/R/dtboundary.R')
source('/export02/data/risa/lesiontools-master/R/centralveins.R')

#subj7TDirs <- c("MS_AN_093/MS_AN_093_25-07-24_7T", "MS_BG_025/MS_BG_025_12-07-22_7T", "MS_CB_020/MS_CB_020_16-11-21_7T", "MS_CB_060/MS_CB_060_06-05-23_7T", "MS_CF_024/MS_CF_024_23-03-23_7T", "MS_CK_082/MS_CK_082_19-04-24_7T", "MS_CR_067/MS_CR_067_03-08-23_7T", "MS_DG_074/MS_DG_074_01-09-23_7T", "MS_DK_013/MS_DK_013_30-08-22_7T", "MS_EB_006/MS_EB_006_22-03-23_7T", "MS_ED_087/MS_ED_087_17-04-24_7T", "MS_FG_079/MS_FG_079_22-02-24_7T", "MS_GC_076/MS_GC_076_27-10-23_7T", "MS_JB_004/MS_JB_004_01-09-22_7T", "MS_JB_008/MS_JB_008_09-02-22_7T", "MS_JNP_069/MS_JNP_069_01-09-23_7T", "MS_KSF_084/MS_KSF_084_24-04-24_7T", "MS_LG_051/MS_LG_051_12-05-22_7T", "MS_LM_023/MS_LM_023_30-11-21_7T", "MS_LV_094/MS_LV_094_15-07-24_7T", "MS_MC_018/MS_MC_018_09-11-22_7T", "MS_MF_068/MS_MF_068_02-08-23_7T", "MS_MP_014/MS_MP_014_28-06-22_7T", "MS_MPG_039/MS_MPG_039_22-11-22_7T", "MS_NE_095/MS_NE_095_02-08-24_7T", "MS_PL_096/MS_PL_086_24-07-24_7T", "MS_PM_026/MS_PM_026_21-10-22_7T", "MS_PM_056/MS_PM_056_21-12-21_7T", "MS_PT_089/MS_PT_089_18-04-24_7T", "MS_RC_066/MS_RC_066_21-08-23_7T", "MS_SC_002/MS_SC_002_07-12-22_7T", "MS_SG_034/MS_SG_034_31-08-22_7T", "MS_SH_019/MS_SH_019_24-08-22_7T", "MS_SK_022/MS_SK_022_29-08-22_7T", "MS_SP_091/MS_SP_091_13-06-24_7T", "MS_VMA_090/MS_VMA_090_24_04_24_7T", "MS_WEH_023/MS_WEH_023_19-10-22_7T")
subj7TDirs <- c("MS_CB_020/MS_CB_020_16-11-21_7T", "MS_CF_024/MS_CF_024_23-03-23_7T", "MS_CK_082/MS_CK_082_19-04-24_7T", "MS_CR_067/MS_CR_067_03-08-23_7T", "MS_DG_074/MS_DG_074_01-09-23_7T", "MS_DK_013/MS_DK_013_30-08-22_7T", "MS_EA_081/MS_EA_081_19-04-24_7T", "MS_FG_079/MS_FG_079_22-02-24_7T", "MS_GC_076/MS_GC_076_27-10-23_7T", "MS_JB_004/MS_JB_004_01-09-22_7T", "MS_JB_008/MS_JB_008_09-02-22_7T", "MS_JNP_069/MS_JNP_069_01-09-23_7T", "MS_KSF_084/MS_KSF_084_24-04-24_7T", "MS_LG_051/MS_LG_051_12-05-22_7T", "MS_LM_023/MS_LM_023_30-11-21_7T", "MS_LV_094/MS_LV_094_15-07-24_7T", "MS_MC_018/MS_MC_018_09-11-22_7T", "MS_MF_068/MS_MF_068_02-08-23_7T", "MS_MP_014/MS_MP_014_28-06-22_7T", "MS_MPG_039/MS_MPG_039_22-11-22_7T", "MS_NE_095/MS_NE_095_02-08-24_7T", "MS_PL_096/MS_PL_086_24-07-24_7T", "MS_PM_026/MS_PM_026_21-10-22_7T", "MS_PM_056/MS_PM_056_21-12-21_7T", "MS_RC_066/MS_RC_066_21-08-23_7T", "MS_SC_002/MS_SC_002_07-12-22_7T", "MS_SG_034/MS_SG_034_31-08-22_7T", "MS_SH_019/MS_SH_019_24-08-22_7T", "MS_VMA_090/MS_VMA_090_24_04_24_7T")
homeDir <- "/export02/data/risa/01_DOD_Study/"
subj7TDirs <- paste(homeDir,subj7TDirs,sep="")
for (subj7TDir in subj7TDirs)
{
#subj7TDir <- "MS_EA_081/MS_EA_081_19-04-24_7T"
probmapPath <- sprintf("%s/lesion_prob_semi.nii",subj7TDir)
print(probmapPath)
binmapPath <- sprintf("%s/WM_lesions.nii",subj7TDir)
lesion.probs <- readnii(probmapPath)
lesion.mask <- readnii(binmapPath)

clusters <- lesionclusters(probmap=lesion.probs, binmap=lesion.mask, 
                           parallel=TRUE, cores=4, minCenterSize=30, 
                           c3d_path="/export02/data/risa/c3d-1.1.0-Linux-x86_64/bin/c3d")
writenii(nim=clusters[["lesionclusters.nn"]], filename=sprintf("%s/lesion_clusters_semi.nii",subj7TDir))

}

subj7TDir <- "/export02/data/risa/01_DOD_Study/MS_PT_089/MS_PT_089_18-04-24_7T"
probmapPath <- sprintf("%s/lesion_prob_semi.nii",subj7TDir)
binmapPath <- sprintf("%s/WM_lesions.nii",subj7TDir)
lesion.probs <- readnii(probmapPath)
lesion.mask <- readnii(binmapPath)

clusters <- lesionclusters(probmap=lesion.probs, binmap=lesion.mask, 
                           parallel=TRUE, cores=8, minCenterSize=90, 
                           c3d_path="/export02/data/risa/c3d-1.1.0-Linux-x86_64/bin/c3d")
writenii(nim=clusters[["lesionclusters.nn"]], filename=sprintf("%s/lesion_clusters_semi.nii",subj7TDir))
