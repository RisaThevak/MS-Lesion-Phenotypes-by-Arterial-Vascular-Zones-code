#!/bin/bash

SUBJECT_MAIN_FOL=$1
SUBJECT_7TFOL=${SUBJECT_MAIN_FOL}/$2
#SUBJECT_PETFOL=${SUBJECT_MAIN_FOL}/$3
SUBJECT_PETFOL=${SUBJECT_MAIN_FOL} #Temporary 
SUBJECT_7TDIR=/export02/data/risa/01_DOD_Study/$SUBJECT_7TFOL
PET_MAIN_DIR=/export02/data/risa/01_DOD_Study/$SUBJECT_PETFOL
PET_TRANS_DIR=${SUBJECT_7TDIR}/PETOutput
SPM12_DIR=/usr/local/spm12/bin/

PET_emi_alltpts=${PET_TRANS_DIR}/PET_avg_native_alltpts
T1w=${SUBJECT_7TDIR}/MP2RAGE_UNIDEN_PREGD_corr
MRI_Template=${SUBJECT_7TDIR}/MP2RAGE_UNIDEN_brain_forPETreg 

tmp=(${PET_MAIN_DIR}/*_4D_MC01.mnc)
tmp=${tmp[0]}

if ! [ -s ${PET_TRANS_DIR}/ ]; then 
	mkdir ${PET_TRANS_DIR}/
fi 

if ! [ -s ${SUBJECT_7TDIR}/FIRST_subcort/ ]; then 
	mkdir ${SUBJECT_7TDIR}/FIRST_subcort/
fi 

UNIDEN_corr_brain=${SUBJECT_7TDIR}/FIRST_subcort/MP2RAGE_UNIDEN_brain_forPETreg
FSLDIR_tmp=/export02/data/risa/fsl-6.0

if ! [ -s ${SUBJECT_7TDIR}/FIRST_subcort ]; then 
	mkdir ${SUBJECT_7TDIR}/FIRST_subcort
fi 

#1) Compute 11C-PBR28 PET SUVR map in MATLAB
if ! [ -s ${PET_TRANS_DIR}/PET_SUVR_T1w_reg_noPVC.nii ]; then 
	read -p "Enter body weight of study participant in kg: " bw
	read -p "Enter injected radiotracer dose of study participant in mCi: " injd
	string3=$(printf "addpath(genpath(\"/export02/data/risa/01_DOD_Study/0_Scripts_PET\")); calSUVR_noPVC(%s, %s, \"%s\"); exit;" $bw $injd ${SUBJECT_7TFOL}) 
	matlab20b -nodesktop -nosplash -batch "$string3"
fi 

#2) Blur 11C-PBR28 PET SUVR map using SPM12 "smooth" function 
if ! [ -s ${PET_TRANS_DIR}/blur_PET_SUVR_T1w_reg_noPVC.nii ]; then 
	echo "Filter PET_SUVR_T1w_reg image with Gaussian Kernel using SPM12 'smooth':"
	echo "Rationale: This allows for localization of PET SUVR differences between MS patients and Healthy controls."

	read -p "Enter prefix (usually 'blur_'): " prefix 

	sentinel="n"
	while [ "${sentinel}" = "n" ]; do 
		read -p "Enter FHWM (mm) of Gaussian Kernel for smoothing (optimal = 4mm FWHM): " FWHM 

		string3=$(printf "addpath(genpath(\"/export02/data/risa/01_DOD_Study/0_Scripts_PET/\")); SPM12smooth(\"%s\", %s, \"%s\"); exit;" ${PET_TRANS_DIR}/PET_SUVR_T1w_reg_noPVC.nii $FWHM $prefix) 
		matlab20b -nodesktop -nosplash -batch "$string3"

		remote fsleyes ${SUBJECT_7TDIR}/MP2RAGE_UNIDEN_PREGD_corr.nii -dr 0 1000000 -a 100 ${PET_TRANS_DIR}/${prefix}PET_SUVR_T1w_reg_noPVC.nii -dr 0 2 -a 40 -cm brain_colours_actc ${PET_TRANS_DIR}/PET_SUVR_T1w_reg_noPVC.nii -dr 0 2 -a 40 -cm brain_colours_actc ${PET_TRANS_DIR}/RefMaskSVCA.nii -dr 0 1 -a 1 -cm blue
		read -p "Are you satisfied gaussian filtering of PET SUVR image? (y/n): " sentinel
	done  
fi 

rm -f ${MRI_Template}.mnc

