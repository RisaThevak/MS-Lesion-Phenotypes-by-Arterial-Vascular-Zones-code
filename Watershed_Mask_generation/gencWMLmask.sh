#!/bin/bash 

NIFTYREG_INSTALL=/export02/data/risa/niftyreg_install/bin
fslDir=/export02/data/risa/fsl-6.0/bin
MNI152_T1_1mm_orig=/export02/data/risa/01_DOD_Study/MNI152_T1_1mm_brain

sumWMLs_ws=/export02/data/risa/01_DOD_Study/sumwsWMLs
binaryWMLs_ws=/export02/data/risa/01_DOD_Study/binarywsWMLs
sumWMLs_mixed=/export02/data/risa/01_DOD_Study/summixedWMLs
binaryWMLs_mixed=/export02/data/risa/01_DOD_Study/binarymixedWMLs
sumWMLs_nws=/export02/data/risa/01_DOD_Study/sumnwsWMLs
binaryWMLs_nws=/export02/data/risa/01_DOD_Study/binarynwsWMLs
${fslDir}/fslmaths ${MNI152_T1_1mm_orig}.nii.gz -mul 0 ${sumWMLs_ws}.nii.gz
${fslDir}/fslmaths ${MNI152_T1_1mm_orig}.nii.gz -mul 0 ${sumWMLs_mixed}.nii.gz
${fslDir}/fslmaths ${MNI152_T1_1mm_orig}.nii.gz -mul 0 ${sumWMLs_nws}.nii.gz

# Define an array of MS patient folder names 
#MSvolumes=(MS_AN_093/MS_AN_093_25-07-24_7T MS_BG_025/MS_BG_025_12-07-22_7T MS_CB_020/MS_CB_020_16-11-21_7T MS_CB_060/MS_CB_060_06-05-23_7T MS_CF_024/MS_CF_024_23-03-23_7T MS_CK_082/MS_CK_082_19-04-24_7T MS_CR_067/MS_CR_067_03-08-23_7T MS_DG_074/MS_DG_074_01-09-23_7T MS_DK_013/MS_DK_013_30-08-22_7T MS_EA_081/MS_EA_081_19-04-24_7T MS_EB_006/MS_EB_006_22-03-23_7T MS_ED_087/MS_ED_087_17-04-24_7T MS_FG_079/MS_FG_079_22-02-24_7T MS_GC_076/MS_GC_076_27-10-23_7T MS_JB_004/MS_JB_004_01-09-22_7T MS_JB_008/MS_JB_008_09-02-22_7T MS_JNP_069/MS_JNP_069_01-09-23_7T MS_KSF_084/MS_KSF_084_24-04-24_7T MS_LG_051/MS_LG_051_12-05-22_7T MS_LM_023/MS_LM_023_30-11-21_7T MS_LV_094/MS_LV_094_15-07-24_7T MS_MC_018/MS_MC_018_09-11-22_7T MS_MF_068/MS_MF_068_02-08-23_7T MS_MP_014/MS_MP_014_28-06-22_7T MS_MPG_039/MS_MPG_039_22-11-22_7T MS_NE_095/MS_NE_095_02-08-24_7T MS_PL_096/MS_PL_086_24-07-24_7T MS_PM_026/MS_PM_026_21-10-22_7T MS_PM_056/MS_PM_056_21-12-21_7T MS_PT_089/MS_PT_089_18-04-24_7T MS_RC_066/MS_RC_066_21-08-23_7T MS_SC_002/MS_SC_002_07-12-22_7T MS_SG_034/MS_SG_034_31-08-22_7T MS_SH_019/MS_SH_019_24-08-22_7T MS_SK_022/MS_SK_022_29-08-22_7T MS_SP_091/MS_SP_091_13-06-24_7T MS_VMA_090/MS_VMA_090_24_04_24_7T MS_WEH_023/MS_WEH_023_19-10-22_7T)
MSvolumes=(MS_AN_093/MS_AN_093_25-07-24_7T MS_BG_025/MS_BG_025_12-07-22_7T MS_CB_020/MS_CB_020_16-11-21_7T MS_CB_060/MS_CB_060_06-05-23_7T MS_CF_024/MS_CF_024_23-03-23_7T MS_CK_082/MS_CK_082_19-04-24_7T MS_CR_067/MS_CR_067_03-08-23_7T MS_DG_074/MS_DG_074_01-09-23_7T MS_DK_013/MS_DK_013_30-08-22_7T MS_EA_081/MS_EA_081_19-04-24_7T MS_EB_006/MS_EB_006_22-03-23_7T MS_FG_079/MS_FG_079_22-02-24_7T MS_GC_076/MS_GC_076_27-10-23_7T MS_JB_004/MS_JB_004_01-09-22_7T MS_JB_008/MS_JB_008_09-02-22_7T MS_JNP_069/MS_JNP_069_01-09-23_7T MS_KSF_084/MS_KSF_084_24-04-24_7T MS_LG_051/MS_LG_051_12-05-22_7T MS_LM_023/MS_LM_023_30-11-21_7T MS_LV_094/MS_LV_094_15-07-24_7T MS_MC_018/MS_MC_018_09-11-22_7T MS_MF_068/MS_MF_068_02-08-23_7T MS_MP_014/MS_MP_014_28-06-22_7T MS_MPG_039/MS_MPG_039_22-11-22_7T MS_NE_095/MS_NE_095_02-08-24_7T MS_PL_096/MS_PL_086_24-07-24_7T MS_PM_026/MS_PM_026_21-10-22_7T MS_PM_056/MS_PM_056_21-12-21_7T MS_RC_066/MS_RC_066_21-08-23_7T MS_SC_002/MS_SC_002_07-12-22_7T MS_SG_034/MS_SG_034_31-08-22_7T MS_SH_019/MS_SH_019_24-08-22_7T MS_SK_022/MS_SK_022_29-08-22_7T MS_SP_091/MS_SP_091_13-06-24_7T MS_VMA_090/MS_VMA_090_24_04_24_7T MS_WEH_023/MS_WEH_023_19-10-22_7T)

HCvolumes=(HC_AM_071/HC_AM_071_15-11-23_7T HC_AW_073/HC_AW_073_15-09-23_7T HC_BJP_037/HC_BJP_037_19-07-22_7T HC_CA_012/HC_CA_012_12-10-21_7T HC_CR_072/HC_CR_072_07-08-23_7T HC_DR_032/HC_DR_032_09-03-22_7T HC_EDL_039/HC_EDL_039_27-07-22_7T HC_EM_092/HC_EM_092_28-06-24_7T HC_FR_065/HC_FR_065_18-08-23_7T HC_GK_070/HC_GK_070_15-08-23_7T HC_IC_062/HC_IC_062_01-08-23_7T HC_JH_085/HC_JH_085_18-04-24_7T HC_KMD_040/HC_KMD_040_28-06-22_7T HC_KT_033/HC_KT_033_21-12-21_7T HC_LB_083/HC_LB_083_19-04-24_7T HC_MB_036/HC_MB_036_17-05-22_7T HC_MB_075/HC_MB_075_27-10-23_7T HC_MJ_061/HC_MJ_061_11-08-23_7T HC_NA_045/HC_NA_045_07-03-22_7T HC_SG_078/HC_SG_078_20-10-23_7T HC_TK_086/HC_TK_086_03-05-24_7T)

#STEP 1: Generate binary WML masks, specific to ws-, mixed- & nws- arterial vascular zones, in MNI152 space 
if ! [ -s ${binaryWMLs_ws}.nii ]; then 

for dir in ${MSvolumes[@]}; do #Going through all MS patient folder names 

	SUBJECT_7TDIR=/export02/data/risa/01_DOD_Study/$dir
	UNIDEN_lin=${SUBJECT_7TDIR}/ResampledInMNI152Space/MP2RAGE_UNIDEN_brain_lin
	WMLs_vascularlbl=${SUBJECT_7TDIR}/ResampledInMNI152Space/WM_lesions_vascularlbl
	WMLs_ws=${SUBJECT_7TDIR}/ResampledInMNI152Space/WM_lesions_ws
	WMLs_mixed=${SUBJECT_7TDIR}/ResampledInMNI152Space/WM_lesions_mixed
	WMLs_nws=${SUBJECT_7TDIR}/ResampledInMNI152Space/WM_lesions_nws
	MNI152_T1_1mm=${SUBJECT_7TDIR}/ResampledInMNI152Space/MNI152_T1_1mm_brain
	
	${fslDir}/fslmaths ${WMLs_vascularlbl}.nii -thr 1 -uthr 1 -bin ${WMLs_ws}.nii
	${fslDir}/fslmaths ${WMLs_vascularlbl}.nii -thr 2 -uthr 2 -bin ${WMLs_mixed}.nii
	${fslDir}/fslmaths ${WMLs_vascularlbl}.nii -thr 3 -uthr 3 -bin ${WMLs_nws}.nii

	nii2mnc ${MNI152_T1_1mm}.nii ${MNI152_T1_1mm}.mnc

	for WMLs in ${WMLs_ws} ${WMLs_mixed} ${WMLs_nws}; do 
		WMLsMNI152=${WMLs}_MNI152
		nii2mnc ${WMLs}.nii ${WMLs}.mnc

		#Initial linear registration step + non-linear registration 
		mincresample -like ${MNI152_T1_1mm}.mnc ${WMLs}.mnc ${WMLsMNI152}_lin.mnc -transformation ${UNIDEN_lin}0_GenericAffine.xfm -clobber 
		mnc2nii ${WMLsMNI152}_lin.mnc ${WMLsMNI152}_lin.nii 
		${NIFTYREG_INSTALL}/reg_resample -ref ${MNI152_T1_1mm}.nii -flo ${WMLsMNI152}_lin.nii -cpp ${UNIDEN_lin}_transform_cpp.nii -res ${WMLsMNI152}.nii -inter 1 #Check interpolation method 

		# Threshold to isolate lesions (assuming lesions are non-zero)
		${fslDir}/fslmaths ${WMLsMNI152}.nii -thr 0 -bin ${WMLsMNI152}.nii; gunzip ${WMLsMNI152}.nii.gz --force

		rm -f ${WMLs}.mnc; rm -f ${WMLsMNI152}_lin.mnc; rm -f ${WMLsMNI152}_lin.nii; rm -f ${WMLs}.nii 

	done 
	rm -f ${MNI152_T1_1mm}.mnc

	${fslDir}/fslmaths ${sumWMLs_ws}.nii.gz -add ${WMLs_ws}_MNI152.nii ${sumWMLs_ws}.nii.gz
        ${fslDir}/fslmaths ${sumWMLs_mixed}.nii.gz -add ${WMLs_mixed}_MNI152.nii ${sumWMLs_mixed}.nii.gz
        ${fslDir}/fslmaths ${sumWMLs_nws}.nii.gz -add ${WMLs_nws}_MNI152.nii ${sumWMLs_nws}.nii.gz

	rm -f ${WMLs_ws}_MNI152.nii; rm -f ${WMLs_mixed}_MNI152.nii; rm -f ${WMLs_nws}_MNI152.nii
done 

# Binarized summed WML masks (MNI152 space)
${fslDir}/fslmaths ${sumWMLs_ws}.nii.gz -thr 0 -bin ${binaryWMLs_ws}.nii
${fslDir}/fslmaths ${sumWMLs_mixed}.nii.gz -thr 0 -bin ${binaryWMLs_mixed}.nii
${fslDir}/fslmaths ${sumWMLs_nws}.nii.gz -thr 0 -bin ${binaryWMLs_nws}.nii
gunzip ${binaryWMLs_ws}.nii.gz --force; gunzip ${binaryWMLs_mixed}.nii.gz --force; gunzip ${binaryWMLs_nws}.nii.gz --force

fi 

#STEP 2: Transform binary WML mask in MNI152 space, specific to ws-, mixed- & nws- arterial vascular zones, to each healthy control native sapce
for dir in ${HCvolumes[@]}; do #Going through all healthy control folder names 

	SUBJECT_7TDIR=/export02/data/risa/01_DOD_Study/$dir
	UNIDEN_lin=${SUBJECT_7TDIR}/ResampledInMNI152Space/MP2RAGE_UNIDEN_brain_lin
	UNIDEN=${SUBJECT_7TDIR}/MP2RAGE_UNIDEN_PREGD
	MNI152_T1_1mm=${SUBJECT_7TDIR}/ResampledInMNI152Space/MNI152_T1_1mm_brain
	watershedMask=${SUBJECT_7TDIR}/ResampledInMNI152Space/watershed_mask
	nonwatershedMask=${SUBJECT_7TDIR}/ResampledInMNI152Space/nonwatershed_mask
	
	for vascular_zone in ws mixed nws; do 
		WMLs=/export02/data/risa/01_DOD_Study/binary${vascular_zone}WMLs
		CounterWMLs=${SUBJECT_7TDIR}/ResampledInMNI152Space/CounterWMLs_${vascular_zone}

		${NIFTYREG_INSTALL}/reg_resample -ref ${UNIDEN_lin}.nii -flo ${WMLs}.nii -cpp ${UNIDEN_lin}_transform_cpp_inv.nii -res ${CounterWMLs}_nonlin_inv.nii -inter 1 #Check interpolation method 
		nii2mnc ${CounterWMLs}_nonlin_inv.nii ${CounterWMLs}_nonlin_inv.mnc

		nii2mnc ${UNIDEN}.nii ${UNIDEN}.mnc 
		mincresample -like ${UNIDEN}.mnc ${CounterWMLs}_nonlin_inv.mnc ${CounterWMLs}.mnc -transformation ${UNIDEN_lin}0_GenericAffine_inv.xfm -clobber -trilinear
		mnc2nii ${CounterWMLs}.mnc ${CounterWMLs}.nii 

		#Binarize transformed counter WML image 
		${fslDir}/fslmaths ${CounterWMLs}.nii -thr 0 -bin ${CounterWMLs}.nii.gz

		if [ "${vascular_zone}" == "ws" ]; then #Zero out portions outside of subject-specific ws region 
			${fslDir}/fslmaths ${CounterWMLs}.nii.gz -mas ${watershedMask}.nii ${CounterWMLs}.nii.gz
			
		elif [ "${vascular_zone}" == "nws" ]; then #Zero out portion outside of subject-specific nws region 
			${fslDir}/fslmaths ${watershedMask}.nii -binv ${nonwatershedMask}.nii.gz
			${fslDir}/fslmaths ${CounterWMLs}.nii.gz -mas ${nonwatershedMask}.nii.gz ${CounterWMLs}.nii.gz
		fi 

		#Zero out non-WM values in counter WML image 
		${fslDir}/fslmaths ${CounterWMLs}.nii.gz -mas ${SUBJECT_7TDIR}/cerebralNAWM.nii ${CounterWMLs}.nii.gz

		gunzip ${CounterWMLs}.nii.gz --force; gunzip ${nonwatershedMask}.nii.gz --force
		rm -f ${CounterWMLs}.nii.gz; rm -f ${CounterWMLs}.mnc; rm -f ${CounterWMLs}_nonlin_inv.mnc; rm -f ${CounterWMLs}_nonlin_inv.nii; rm -f ${UNIDEN}.mnc; rm -f ${nonwatershedMask}.nii.gz
	done 

done 

