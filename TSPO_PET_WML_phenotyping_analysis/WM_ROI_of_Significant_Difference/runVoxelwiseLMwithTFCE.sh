#!/bin/bash

'filelist=""
for dir in HC_AM_071 HC_AW_073 HC_BJP_037 HC_CA_012 HC_CR_072 HC_DR_032 HC_EDL_039 HC_EM_092 HC_FR_065 HC_GK_070 HC_IC_062 HC_JH_085 HC_KMD_040 HC_KT_033 HC_LB_083 HC_MB_036 HC_MB_075 HC_NA_045 HC_SG_078 HC_TK_086 MS_AN_093 MS_BG_025 MS_CB_020 MS_CF_024 MS_DG_074 MS_DK_013 MS_EA_081 MS_EB_006 MS_ED_087 MS_FG_079 MS_GC_076 MS_JB_004 MS_JB_008 MS_KSF_084 MS_LG_051 MS_LM_023 MS_MC_018 MS_MF_068 MS_MPG_039 MS_MP_014 MS_NE_095 MS_PL_096 MS_PM_056 MS_PT_089 MS_RC_066 MS_SC_002 MS_SG_034 MS_SH_019 MS_SK_022 MS_SP_091 MS_VMA_090 MS_WEH_023; do 
	SUVR_MNI152=/export02/data/risa/01_DOD_Study/${dir}/*/ResampledInMNI152Space/PET_SUVR_MNI152.nii
	
	filelist="${filelist}${SUVR_MNI152} "
done 
	mrcat ${filelist} /export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/CerebralWMAnalysis/Input_4D_Image.nii  -force'
	
dir=/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/CerebralWMAnalysis
randomise -i ${dir}/Input_4D_Image.nii -o ${dir}/Output -d ${dir}/clinicalStatus_genotype_age.mat -t ${dir}/clinicalStatus_contrasts.con -m ${dir}/MNI152_T1_1mm_cerebralWM_mask.nii -n 5000 -D -T
#Trial 1 - de-meaning data, no smoothing PET SUVR MNI152 images beforehand, 2 HCs are NOT omitted, LABs are excluded
