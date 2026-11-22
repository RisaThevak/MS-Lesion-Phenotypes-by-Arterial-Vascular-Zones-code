%% Resample watershed regions mask to each patient's native MP2RAGE space
clear; clc; close all; 
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/infoPET.mat"); 
NIFTYREG_INSTALL="/export02/data/risa/niftyreg_install/bin";

DirsMS = DirsPET(string(groupsPET)=="MS"); 
NamesMS = NamesPET(string(groupsPET)=="MS"); 
noMS = length(DirsMS); 
for id = 1:noMS
    subjDir = [DirsMS{id} '/' NamesMS{id}];
    
    fname = "/export02/data/risa/Atlas_MNI152/Atlas_MNI152/watershed_regions.nii"; 
    
    UNIDEN_brain_lin=sprintf("%s/ResampledInMNI152Space/MP2RAGE_UNIDEN_brain_lin",subjDir); 
    UNIDEN_corr_brain=sprintf("%s/MP2RAGE_UNIDEN_PREGD",subjDir); 
    prefix = sprintf("%s/ResampledInMNI152Space/watershed_mask",subjDir); 
    com = sprintf("%s/reg_resample -ref %s.nii -flo %s -cpp %s_transform_cpp_inv.nii -res %s_nonlin_inv.nii -inter 1",...
    NIFTYREG_INSTALL, UNIDEN_brain_lin, fname, UNIDEN_brain_lin, prefix); 
    system(char(com)); 
    
    com = sprintf("nii2mnc %s_nonlin_inv.nii %s_nonlin_inv.mnc", prefix, prefix); 
    system(char(com)); 
    
    com = sprintf("nii2mnc %s.nii %s.mnc", UNIDEN_corr_brain, UNIDEN_corr_brain); 
    system(char(com)); 
    
    com = sprintf("mincresample -like %s.mnc %s_nonlin_inv.mnc %s.mnc -transformation %s0_GenericAffine_inv.xfm -clobber -trilinear",...
        UNIDEN_corr_brain, prefix, prefix, UNIDEN_brain_lin); 
    system(char(com)); 
    
    com = sprintf("mnc2nii %s.mnc %s.nii", prefix, prefix); 
    system(char(com)); 
    
    com = sprintf("rm -f %s_nonlin_inv.nii; rm -f %s_nonlin_inv.mnc; rm -f %s.mnc; rm -f %s.mnc",...
        prefix, prefix, prefix, UNIDEN_corr_brain); 
    system(char(com)); 
end 

%% Resample watershed regions mask to each control's native MP2RAGE space
clear; clc; close all; 
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/infoPET.mat"); 
NIFTYREG_INSTALL="/export02/data/risa/niftyreg_install/bin";

DirsHC = DirsPET(string(groupsPET)=="HC"); 
NamesHC = NamesPET(string(groupsPET)=="HC"); 
noHC = length(DirsHC); 
for id = 1:noHC
    subjDir = [DirsHC{id} '/' NamesHC{id}];
    
    fname = "/export02/data/risa/Atlas_MNI152/Atlas_MNI152/watershed_regions.nii"; 
    
    UNIDEN_brain_lin=sprintf("%s/ResampledInMNI152Space/MP2RAGE_UNIDEN_brain_lin",subjDir); 
    UNIDEN_corr_brain=sprintf("%s/MP2RAGE_UNIDEN_PREGD",subjDir); 
    prefix = sprintf("%s/ResampledInMNI152Space/watershed_mask",subjDir); 
    com = sprintf("%s/reg_resample -ref %s.nii -flo %s -cpp %s_transform_cpp_inv.nii -res %s_nonlin_inv.nii -inter 1",...
    NIFTYREG_INSTALL, UNIDEN_brain_lin, fname, UNIDEN_brain_lin, prefix); 
    system(char(com)); 
    
    com = sprintf("nii2mnc %s_nonlin_inv.nii %s_nonlin_inv.mnc", prefix, prefix); 
    system(char(com)); 
    
    com = sprintf("nii2mnc %s.nii %s.mnc", UNIDEN_corr_brain, UNIDEN_corr_brain); 
    system(char(com)); 
    
    com = sprintf("mincresample -like %s.mnc %s_nonlin_inv.mnc %s.mnc -transformation %s0_GenericAffine_inv.xfm -clobber -trilinear",...
        UNIDEN_corr_brain, prefix, prefix, UNIDEN_brain_lin); 
    system(char(com)); 
    
    com = sprintf("mnc2nii %s.mnc %s.nii", prefix, prefix); 
    system(char(com)); 
    
    com = sprintf("rm -f %s_nonlin_inv.nii; rm -f %s_nonlin_inv.mnc; rm -f %s.mnc; rm -f %s.mnc",...
        prefix, prefix, prefix, UNIDEN_corr_brain); 
    system(char(com)); 
end 
