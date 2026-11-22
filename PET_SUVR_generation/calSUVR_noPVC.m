function [] = calSUVR_noPVC(bw, id, subjMRIDir)
%INPUTS: 
%1) bw: patient body weight in kg 
%2) id: injected radiotracer dose in mCi 
%3) subjMRIDir: name of 7 T MRI directory of subject 

addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 

bw=bw*1.e3; %kg to g
id=id*1.e6; %mCi to nCi

subjMRIDir = char(subjMRIDir); 
home = '/export02/data/risa/01_DOD_Study/'; 
subjMRIDir = [home subjMRIDir]; 
outputFol = [subjMRIDir '/PETOutput/']; 

%Create binary mask for pseudoreference region
refMask = double(load_untouch_nii([outputFol 'RefMaskSVCA.nii']).img); 
m_refROI=single(zeros(size(refMask))); 
m_refROI(refMask~=0)=1; 

%Determining number of voxels in ref ROI mask
s_refROI=sum(m_refROI(:)); 

%Read in averaged+registered PET emission image that was registered to MP2RAGE UNI space earlier: 
volem = load_untouch_nii([outputFol 'PETPVC/PET_avg_T1w_reg_onlybrain.nii']); %read in file via nifti
tmp = volem.hdr; volem = double(volem.img); 

%Compute mean SUV for reference ROI
suv_refROI=sum(volem(:).*m_refROI(:))/s_refROI/id*bw;

%Save PET SUV image as nifti file 
suv_tmp = volem/(id/bw); 
suv_nii = make_nii(suv_tmp); 
suv_nii.hdr.dime.pixdim = tmp.dime.pixdim; suv_nii.hdr.hist = tmp.hist; 
save_nii(suv_nii,[outputFol 'PET_SUV_T1w_reg_noPVC.nii']); %save through nifti 

%IMPORTANT: Compute PET SUVR map 
suvr = volem/(id*suv_refROI/bw);

%Save PET SUVR map as nifti file and minc file 
suvr_nii = make_nii(suvr); 
suvr_nii.hdr.dime.pixdim = tmp.dime.pixdim; suvr_nii.hdr.hist = tmp.hist; 
save_nii(suvr_nii,[outputFol 'PET_SUVR_T1w_reg_noPVC.nii']); %save through nifti 
end 
