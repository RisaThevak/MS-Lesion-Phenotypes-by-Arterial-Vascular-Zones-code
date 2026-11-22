%% 1) Determine 'lesion volume' threshold to define confluent lesion units

WMLsVolcollapsed = [vertcat(WMLsVol.watershed{:});...
                    vertcat(WMLsVol.mixed{:});...
                    vertcat(WMLsVol.nonWatershed{:})]; 
                
% WMLsVolcollapsed = [vertcat(WMLsVol.watershed{:});...
%                     vertcat(WMLsVol.mixed{:})];                 

histogram(WMLsVolcollapsed(:),'binWidth',0.5); 
title('Frequency of WML volumes'); 
xlabel('Individual WML volume (cm3)'); 
ylabel('Lesion count'); 
xlim([0 60]);
ylim([0 20]); 
Q = quantile(WMLsVolcollapsed(:),[0.025 0.25 0.5 0.75 0.950 0.975]); 

%% 2) Determine 'minimum lesion volume' threshold to filter out random T1 hypointensities or T2* hyperintensities 

%FWHM of PET system PSF along transaxial dimension is ~2.8mm -> minimum
%lesion volume = (3)^3 = 27mm3 to avoid PET PVEs (along if object is sphere
%of size > 2.5-3x the FWHM, then the maximum PET emission value approaches
%the true mean PET emission intensity of the object, assuming homogeneity
%and noise free) 

% WMLsVolcollapsed = [vertcat(WMLsVol.watershed{:});...
%                     vertcat(WMLsVol.mixed{:});...
%                     vertcat(WMLsVol.nonWatershed{:})]; 
                
WMLsVolcollapsed = [vertcat(WMLsVol.mixed{:})];                 

histogram(WMLsVolcollapsed(:),'binWidth',0.01); 
title('Frequency of WML volumes'); 
xlabel('Individual WML volume (cm3)'); 
ylabel('Lesion count');  
xlim([0 1.2]); 
Q = quantile(WMLsVolcollapsed(:),[0.025 0.25 0.5 0.75 0.950 0.975]); 
median(WMLsVolcollapsed(:))
range(WMLsVolcollapsed(:))
%Given this being a left-skewed distribution, the median seems to occur at
%around 50mm3 

%%
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/HC_data_WMLphenotyping.mat");
histogram(tmpblur,'binWidth',0.05); 
title('Frequency of WML volumes'); 
xlabel('Individual WML volume (cm3)'); 
ylabel('Lesion count'); 
ylim([0 7]); 
Q = quantile(WMLsVolcollapsed(:),[0.025 0.25 0.5 0.75 0.950 0.975]);

%% 3) Exporting Resampled + Smoothed PET SUVR Data 
clear; close all; 
addpath(genpath("/export02/data/risa/NIfTI_20140122/"));
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/infoPET.mat"); 
FWHM=2.5; %6.0 blurring kernel 

noSubjs = length(DirsPET); 
for id = 1:noSubjs
    fprintf(['Processing subject: ' NamesPET{id} '\n']);
    subj7TDir = [DirsPET{id} '/' NamesPET{id}];  
    
    %blurring resampled PET suvr image using 3D geodesic kernel of given FWHM
    imgPath = [subj7TDir '/PETOutput/PET_SUVR_T1w_reg.nii'];  
%     imgPath = [subj7TDir '/PETOutput/PET_SUVR_T1w_reg_noPVC.nii']; %TEMPORARY
    SPM12smooth(imgPath,FWHM,sprintf("blur%.1fmm_",FWHM)); 
end

%% 4a) Set thresholds for TSPO phenotyping using healthy control group - why Nylund specifies in her 2026 JN paper 
clc; clear; close all; 
addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/infoPET.mat"); 

% DirsHC = DirsPET(string(groupsPET)=="HC"); 
% NamesHC = NamesPET(string(groupsPET)=="HC"); 
DirsHC = DirsPET(string(groupsPET)=="HC" & string(bindAff')=="HAB"); %CHANGE
NamesHC = NamesPET(string(groupsPET)=="HC" & string(bindAff')=="HAB"); %CHANGE
noHC = length(DirsHC);

WMsuvr = cell(1,noHC); WMsuvrblur = WMsuvr;
meanWMsuvrHC = zeros(1,noHC); meanWMsuvrHCblur = meanWMsuvrHC; stdWMsuvrHC = meanWMsuvrHC; stdWMsuvrHCblur = meanWMsuvrHCblur; 
thresh1 = meanWMsuvrHC; thresh2 = meanWMsuvrHCblur; 
for id = 1:noHC 
    CAT12mrisDir = [DirsHC{id} '/' NamesHC{id} '/CAT12_CSE/mri']; 
    Dir = [DirsHC{id} '/' NamesHC{id}]; 
    
    WM = double(load_untouch_nii([Dir '/cerebralNAWM.nii']).img); %okay since healthy controls have no WMLs
 
    SUVR = double(load_untouch_nii([Dir '/PETOutput/PET_SUVR_T1w_reg.nii']).img);
    WMsuvr{id} = SUVR(WM==1)'; 
    SUVR = double(load_untouch_nii([Dir '/PETOutput/blur2.5mm_PET_SUVR_T1w_reg.nii']).img);
    WMsuvrblur{id} = SUVR(WM==1)'; 
end 

for id = 1:noHC
    meanWMsuvrHC(id) = mean(WMsuvr{id},'all'); stdWMsuvrHC(id) = std(WMsuvr{id},0,'all');
    meanWMsuvrHCblur(id) = mean(WMsuvrblur{id},'all'); stdWMsuvrHCblur(id) = std(WMsuvrblur{id},0,'all');
    thresh1(id) = meanWMsuvrHC(id) + 1.96*stdWMsuvrHC(id);
    thresh2(id) = meanWMsuvrHCblur(id) + 1.96*stdWMsuvrHCblur(id);

    fprintf('id = %d\n',id)
end 

threshSUVR = mean(thresh1,'all'); threshSUVRblur = mean(thresh2,'all'); 

%CHANGE NAME WITH 'HAB', 'MAB' or no addition
save('/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/HC_HABdata_WMLphenotyping.mat','threshSUVR', 'threshSUVRblur');

%% 4b) Set thresholds for TSPO phenotyping using healthy control group - What I initially did
clc; clear; close all; 
addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/infoPET.mat"); 

% DirsHC = DirsPET(string(groupsPET)=="HC"); 
% NamesHC = NamesPET(string(groupsPET)=="HC"); 
DirsHC = DirsPET(string(groupsPET)=="HC" & string(bindAff')=="MAB"); %CHANGE
NamesHC = NamesPET(string(groupsPET)=="HC" & string(bindAff')=="MAB"); %CHANGE
noHC = length(DirsHC);

WMsuvr = cell(1,noHC); WMsuvrblur = WMsuvr; 
for id = 1:noHC 
    CAT12mrisDir = [DirsHC{id} '/' NamesHC{id} '/CAT12_CSE/mri']; 
    Dir = [DirsHC{id} '/' NamesHC{id}]; 
    
    WM = double(load_untouch_nii([Dir '/cerebralNAWM.nii']).img); %okay since healthy controls have no WMLs
 
    SUVR = double(load_untouch_nii([Dir '/PETOutput/PET_SUVR_T1w_reg.nii']).img);
    WMsuvr{id} = SUVR(WM==1)'; 
    SUVR = double(load_untouch_nii([Dir '/PETOutput/blur2.5mm_PET_SUVR_T1w_reg.nii']).img);
    WMsuvrblur{id} = SUVR(WM==1)'; 

end 

meanWMsuvrHC = zeros(1,noHC); meanWMsuvrHCblur = meanWMsuvrHC; 
for i = 1:noHC
    meanWMsuvrHC(i) = mean(WMsuvr{i},'all'); 
    meanWMsuvrHCblur(i) = mean(WMsuvrblur{i},'all'); 
end 
tmpblur = meanWMsuvrHCblur; tmp = meanWMsuvrHC; 
stdWMsuvrHC = std(meanWMsuvrHC,0,'all'); meanWMsuvrHC = mean(meanWMsuvrHC,'all'); 
stdWMsuvrHCblur = std(meanWMsuvrHCblur,0,'all'); meanWMsuvrHCblur = mean(meanWMsuvrHCblur,'all');  

threshSUVR = meanWMsuvrHC + 1.96*stdWMsuvrHC;
threshSUVRblur = meanWMsuvrHCblur + 1.96*stdWMsuvrHCblur; %Based on Nylund (2020) paper

%CHANGE NAME WITH 'HAB', 'MAB' or no addition
save('/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/HC_MABdata_WMLphenotyping.mat','threshSUVR', 'threshSUVRblur');

%% 4c) Set thresholds for TSPO phenotyping using healthy control group - Treaba (2021) method
clc; clear; close all; 
addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/infoPET.mat"); 

DirsHC = DirsPET(string(groupsPET)=="HC" & string(bindAff')=="LAB"); %CHANGE
NamesHC = NamesPET(string(groupsPET)=="HC" & string(bindAff')=="LAB"); %CHANGE
noHC = length(DirsHC);
DirsMS = DirsPET(string(groupsPET)=="MS" & string(bindAff')=="LAB"); %CHANGE
NamesMS = NamesPET(string(groupsPET)=="MS" & string(bindAff')=="LAB"); %CHANGE
noMS = length(DirsMS);

meanROISUVRMS = zeros(1,noMS);
for id = 1:noMS 
    CAT12mrisDir = [DirsMS{id} '/' NamesMS{id} '/CAT12_CSE/mri']; 
    Dir = [DirsMS{id} '/' NamesMS{id}]; 

    % SUVRMNI152 = double(load_untouch_nii([Dir '/ResampledInMNI152Space/PET_SUVR_MNI152.nii']).img);
    SUVRMNI152 = double(load_untouch_nii([Dir '/ResampledInMNI152Space/blur2.500000e+00mm_PET_SUVR_MNI152.nii']).img);
    ROIDir='/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/CerebralWMAnalysis'; 
    ROI = double(load_untouch_nii([ROIDir '/Output_tfce_corrp_tstat1.nii']).img);
    meanROISUVRMS(id) = mean(SUVRMNI152(ROI>0.95),'all'); 

end 

WMsuvr = cell(1,noHC); meanWMsuvrHC = zeros(1,noHC); meanROISUVRHC = zeros(1,noHC); 
for id = 1:noHC 
    CAT12mrisDir = [DirsHC{id} '/' NamesHC{id} '/CAT12_CSE/mri']; 
    Dir = [DirsHC{id} '/' NamesHC{id}]; 
    
    WM = double(load_untouch_nii([Dir '/cerebralNAWM.nii']).img); %okay since healthy controls have no WMLs
 
    % SUVR = double(load_untouch_nii([Dir '/PETOutput/PET_SUVR_T1w_reg.nii']).img);
    SUVR = double(load_untouch_nii([Dir '/PETOutput/blur2.5mm_PET_SUVR_T1w_reg.nii']).img);
    meanWMsuvrHC(id) = mean(SUVR(WM==1),'all'); 

    % SUVRMNI152 = double(load_untouch_nii([Dir '/ResampledInMNI152Space/PET_SUVR_MNI152.nii']).img);
    SUVRMNI152 = double(load_untouch_nii([Dir '/ResampledInMNI152Space/blur2.500000e+00mm_PET_SUVR_MNI152.nii']).img);
    ROIDir='/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/CerebralWMAnalysis'; 
    ROI = double(load_untouch_nii([ROIDir '/Output_tfce_corrp_tstat1.nii']).img);
    meanROISUVRHC(id) = mean(SUVRMNI152(ROI>0.95),'all'); 
    
end 

meanWMsuvrHC = mean(meanWMsuvrHC,'all'); %group-level mean 
meanROISUVRHC = mean(meanROISUVRHC,'all'); %group-level mean
meanROISUVRMS = mean(meanROISUVRMS,'all'); %group-level mean
 
threshSUVRblur = (meanROISUVRMS - meanROISUVRHC) + meanWMsuvrHC; 
%CHANGE NAME WITH 'HAB', 'MAB' or no addition
save('/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/HC_LABdata_WMLphenotyping.mat', 'threshSUVRblur');

%% 5) Export measures at individual WML level specific to ws/mixed/nws classification
clear; close all; clc; 
tic
addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/infoPET.mat"); 
conn = 6; %CHECK

%CHANGE this definition if omitting PVWM or not 
confluentThresh = 14; %cm3 %TESTING
minLesionVol = 0.027; %cm3

SPlesions = struct('MS_AN_093',nan, 'MS_BG_025', nan, 'MS_CB_020', nan, 'MS_CB_060', nan, 'MS_CF_024', [28, 58, 95, 70],...
'MS_CK_082', [158], 'MS_CR_067', [219, 104, 206, 185, 146, 164, 172, 34, 107, 119],'MS_DG_074', [10],...
'MS_DK_013', [1, 2], 'MS_EA_081', [141, 175, 215, 158, 199, 50, 34], 'MS_EB_006', [192, 169, 221, 226],... 
'MS_ED_087', nan, 'MS_FG_079', [103, 50, 36, 40, 41],'MS_GC_076', [48, 61, 54], 'MS_JB_004', [15], 'MS_JB_008', [23],...
'MS_JNP_069', [1], 'MS_KSF_084', [66], 'MS_LG_051', [29, 37, 30, 1], 'MS_LM_023', nan,...
'MS_LV_094', [191, 208, 254, 161, 160], 'MS_MC_018', [8, 15, 12], 'MS_MF_068', [101, 67, 57], 'MS_MPG_039', [25, 30],... 
'MS_MP_014', [14], 'MS_NE_095', [158, 170], 'MS_PL_096', [102, 149, 157], 'MS_PM_026', [44, 127, 101],...
'MS_PM_056', [81, 69], 'MS_PT_089', nan, 'MS_RC_066', [128, 155, 108], 'MS_SC_002', [1, 117], 'MS_SG_034', [74, 80, 84],... 
'MS_SH_019', [66, 50, 54, 46, 41], 'MS_SK_022', [2], 'MS_SP_051', nan, 'MS_VMA_090', [130, 170, 184], 'MS_WEH_023', [1]); 

DirsMS = DirsPET(string(groupsPET)=="MS"); 
NamesMS = NamesPET(string(groupsPET)=="MS"); 
TIVMS = TIV(string(groupsPET)=="MS");
noMS = length(DirsMS);
bindAffMS = bindAff(string(groupsPET)=="MS"); 

% Initialize a structure to store the MTsat values for each patient and lesion classification
WMLsMTsat = cell(1, noMS); 

% Similarly,...
WMLsSUVR = WMLsMTsat; WMLsVol = WMLsMTsat; WMLscount = WMLsMTsat; 
RimActive = WMLsMTsat; Active = WMLsMTsat; Inactive = WMLsMTsat;
RimActiveVol = WMLsMTsat; ActiveVol = WMLsMTsat; InactiveVol = WMLsMTsat;
RimActiveMTsat = WMLsMTsat; ActiveMTsat = WMLsMTsat; InactiveMTsat = WMLsMTsat;
RimActiveSUVR = WMLsMTsat; ActiveSUVR = WMLsMTsat; InactiveSUVR = WMLsMTsat;
RimActivePerCore = WMLsMTsat; ActivePerCore = WMLsMTsat; 

%OMIT ED_087!!!!!
% Iterate through each patient's directory 
parfor id = 1:noMS
% for id = 4:4
    mainDir = DirsMS{id}; 
    subj7TDir = [DirsMS{id} '/' NamesMS{id}]; 

    name = [subj7TDir '/WM_lesions.nii'];  
    WMLs = load_untouch_nii(name); 
    DimDat = WMLs.hdr.dime.pixdim; Hist = WMLs.hdr.hist; 
    WMLs = double(WMLs.img); WMLs(WMLs~=0)=1; 
    WMLvascularlbl = zeros(size(WMLs)); %Initialize WML mask with arterial vascular zone classification
    TSPOActClass = zeros(size(WMLs)); %Initialize WML mask with TSPO PET phenotype classification 

    % Temporary variables for the current patient
    MTsatwsl = []; MTsatmixedl = []; MTsatnwsl = [];
    SUVRwsl = []; SUVRmixedl = []; SUVRnwsl= []; 
    Volwsl = []; Volmixedl = []; Volnwsl = []; 
    
    % Initialize counts for TSPO PET-based WML phenotypes 
    noInactive = zeros(1,3); noActive = zeros(1,3); noRimActive = zeros(1,3); noWMLs = zeros(1,3); 
    IAVol = cell(1,3); AVol = cell(1,3); RAVol = cell(1,3); 
    IAMTsat = cell(1,3); AMTsat = cell(1,3); RAMTsat = cell(1,3);
    IASUVR = cell(1,3); ASUVR = cell(1,3); RASUVR = cell(1,3);
    APerCore = cell(1,3); RAPerCore = cell(1,3); 

    CC = bwconncomp(WMLs,conn); dim = DimDat(2:4); 
    for i = 1:length(CC.PixelIdxList) % Iterate through lesion units
        %Logical array specifying current WM lesion
        lesionMask = zeros(size(WMLs)); 
        lesionMask(CC.PixelIdxList{i}) = 1;   
        
        %Calculate volume of current lesion
        %In the case of a confluent lesion unit, this would be the volume of the entire confluent unit 
        volLesion = dim(1)*dim(2)*dim(3)*length(CC.PixelIdxList{i})/1000; %mm3 to cm3
        
        name = [subj7TDir '/lesion_clusters_semi.nii.gz'];
        if (isfile(name)) %If current lesion is one to be excluded in septum (but doesn't meet volume for CONFLUENT LESION UNIT)
           %Load lesion instance mask from J.J. Dworkin method of confluent
           %lesion segmentation
           DworkinLIM =  double(load_untouch_nii(name).img); 

           sentinel = 0; 
           %Remove confounding lesions in septum pellucidum 
           SP = strsplit(mainDir, '/'); 
           SPidxs = getfield(SPlesions,SP{end}); 
           if ~isnan(SPidxs)
               for i = 1:length(SPidxs)
                   tmp = DworkinLIM==SPidxs(i); 
                   if all(lesionMask==tmp)
                       sentinel = 1; 
                   end 
               end 
           end 

           if sentinel, continue; end 
        end 
           
        if (isfile(name) && volLesion > confluentThresh) %CONFLUENT LESION UNIT
           %Load lesion instance mask from J.J. Dworkin method of confluent
           %lesion segmentation
           DworkinLIM =  double(load_untouch_nii(name).img); 
           
           %Remove confounding lesions in septum pellucidum 
           SP = strsplit(mainDir, '/'); 
           SPidxs = getfield(SPlesions,SP{end}); 
           if ~isnan(SPidxs)
               for i = 1:length(SPidxs)
                DworkinLIM(DworkinLIM==SPidxs(i)) = 0; 
               end 
           end 
           
           %Zero out any lesion instances except for those in current confluent unit (lesionMask==1)
           DworkinLIM(lesionMask==0) = 0; 
           %Indices identifying lesion instances 
           C = unique(DworkinLIM,'sorted'); 
           C(C==0) = []; %Remove zero (it is background intensity, not lesion instance)
           
%            disp("Confluent unit detected and MTsatwsl size before: "); disp(size(MTsatwsl));
           for ii = 1:length(C) %Iterate through each lesion in confluent unit 
                lesionMaskCUI = (DworkinLIM == C(ii)); 
                
                [TSPOActClass, val, Perc_ActiveCore] = TSPOPhenotypeNylund(subj7TDir, bindAffMS{id}, lesionMaskCUI, lesionMask,...
                    minLesionVol, TSPOActClass, dim); %Pass here binding affinity phenotype
                
                [MTsatwsl, MTsatmixedl, MTsatnwsl,...
                SUVRwsl, SUVRmixedl, SUVRnwsl,...
                Volwsl, Volmixedl, Volnwsl,...
                noInactive, noActive, noRimActive, noWMLs,...
                IAVol, AVol, RAVol,...
                IAMTsat, AMTsat, RAMTsat,... 
                IASUVR, ASUVR, RASUVR,...
                APerCore, RAPerCore,...
                WMLvascularlbl] = ...
                addWMLdata(subj7TDir, lesionMaskCUI, minLesionVol,...
                MTsatwsl, MTsatmixedl, MTsatnwsl,...
                SUVRwsl, SUVRmixedl, SUVRnwsl,...
                Volwsl, Volmixedl, Volnwsl,...
                noInactive, noActive, noRimActive, noWMLs,...
                IAVol, AVol, RAVol,...
                IAMTsat, AMTsat, RAMTsat,...
                IASUVR, ASUVR, RASUVR,...
                APerCore, RAPerCore,...
                WMLvascularlbl, val, Perc_ActiveCore, dim);
           end 
           
%            disp("Confluent unit detected and MTsatwsl size after: "); disp(size(MTsatwsl));
        else %ISOLATED, NONCONFLUENT LESION
%             disp("Confluent unit NOT detected and MTsatwsl size before: "); disp(size(MTsatwsl));
           
            [TSPOActClass, val, Perc_ActiveCore] = TSPOPhenotypeNylund(subj7TDir, bindAffMS{id}, lesionMask, lesionMask,...
                minLesionVol, TSPOActClass, dim); %pass binding affinity phenotype here
                
            [MTsatwsl, MTsatmixedl, MTsatnwsl,...
            SUVRwsl, SUVRmixedl, SUVRnwsl,...
            Volwsl, Volmixedl, Volnwsl,...
            noInactive, noActive, noRimActive, noWMLs,...
            IAVol, AVol, RAVol,...
            IAMTsat, AMTsat, RAMTsat,...  
            IASUVR, ASUVR, RASUVR,...
            APerCore, RAPerCore,...
            WMLvascularlbl] = ...
            addWMLdata(subj7TDir, lesionMask, minLesionVol,...
            MTsatwsl, MTsatmixedl, MTsatnwsl,...
            SUVRwsl, SUVRmixedl, SUVRnwsl,...
            Volwsl, Volmixedl, Volnwsl,...
            noInactive, noActive, noRimActive, noWMLs,...
            IAVol, AVol, RAVol,...
            IAMTsat, AMTsat, RAMTsat,... 
            IASUVR, ASUVR, RASUVR,...
            APerCore, RAPerCore,...
            WMLvascularlbl, val, Perc_ActiveCore, dim); 
        
%             disp("Confluent unit NOT detected and MTsatwsl size after: "); disp(size(MTsatwsl));
        end   
      
    end   

    WMLsMTsat{id}.watershed = MTsatwsl;
    WMLsMTsat{id}.mixed = MTsatmixedl;
    WMLsMTsat{id}.nonWatershed = MTsatnwsl;
    
    WMLsSUVR{id}.watershed = SUVRwsl;
    WMLsSUVR{id}.mixed = SUVRmixedl;
    WMLsSUVR{id}.nonWatershed = SUVRnwsl;
    
    WMLsVol{id}.watershed = Volwsl;
    WMLsVol{id}.mixed = Volmixedl;
    WMLsVol{id}.nonWatershed = Volnwsl;
    
    WMLscount{id}.watershed = noWMLs(1); 
    WMLscount{id}.mixed = noWMLs(2); 
    WMLscount{id}.nonWatershed = noWMLs(3); 
    
    RimActive{id}.watershed = noRimActive(1);
    RimActive{id}.mixed = noRimActive(2);
    RimActive{id}.nonWatershed = noRimActive(3);
    Active{id}.watershed = noActive(1);
    Active{id}.mixed = noActive(2);
    Active{id}.nonWatershed = noActive(3);
    Inactive{id}.watershed = noInactive(1);
    Inactive{id}.mixed = noInactive(2);
    Inactive{id}.nonWatershed = noInactive(3); 
    
    RimActiveVol{id}.watershed = RAVol{1};
    RimActiveVol{id}.mixed = RAVol{2};
    RimActiveVol{id}.nonWatershed = RAVol{3};
    ActiveVol{id}.watershed = AVol{1};
    ActiveVol{id}.mixed = AVol{2};
    ActiveVol{id}.nonWatershed = AVol{3};
    InactiveVol{id}.watershed = IAVol{1};
    InactiveVol{id}.mixed = IAVol{2};
    InactiveVol{id}.nonWatershed = IAVol{3};  
    
    RimActiveMTsat{id}.watershed = RAMTsat{1};
    RimActiveMTsat{id}.mixed = RAMTsat{2};
    RimActiveMTsat{id}.nonWatershed = RAMTsat{3};
    ActiveMTsat{id}.watershed = AMTsat{1};
    ActiveMTsat{id}.mixed = AMTsat{2};
    ActiveMTsat{id}.nonWatershed = AMTsat{3};
    InactiveMTsat{id}.watershed = IAMTsat{1};
    InactiveMTsat{id}.mixed = IAMTsat{2};
    InactiveMTsat{id}.nonWatershed = IAMTsat{3};   

    RimActiveSUVR{id}.watershed = RASUVR{1};
    RimActiveSUVR{id}.mixed = RASUVR{2};
    RimActiveSUVR{id}.nonWatershed = RASUVR{3};
    ActiveSUVR{id}.watershed = ASUVR{1};
    ActiveSUVR{id}.mixed = ASUVR{2};
    ActiveSUVR{id}.nonWatershed = ASUVR{3};
    InactiveSUVR{id}.watershed = IASUVR{1};
    InactiveSUVR{id}.mixed = IASUVR{2};
    InactiveSUVR{id}.nonWatershed = IASUVR{3};

    RimActivePerCore{id}.watershed = RAPerCore{1};
    RimActivePerCore{id}.mixed = RAPerCore{2};
    RimActivePerCore{id}.nonWatershed = RAPerCore{3};
    ActivePerCore{id}.watershed = APerCore{1};
    ActivePerCore{id}.mixed = APerCore{2};
    ActivePerCore{id}.nonWatershed = APerCore{3};
end
delete(gcp('nocreate')); 

% Save the structure in a single MAT file
outDir = '/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis';

WMLsMTsat = reformStruct(WMLsMTsat); 
WMLsSUVR = reformStruct(WMLsSUVR); 
WMLsVol = reformStruct(WMLsVol); 
RimActive = reformStruct(RimActive); 
Active = reformStruct(Active); 
Inactive = reformStruct(Inactive); 
WMLscount = reformStruct(WMLscount);

RimActiveVol = reformStruct(RimActiveVol); 
ActiveVol = reformStruct(ActiveVol); 
InactiveVol = reformStruct(InactiveVol); 
RimActiveMTsat = reformStruct(RimActiveMTsat); 
ActiveMTsat = reformStruct(ActiveMTsat); 
InactiveMTsat = reformStruct(InactiveMTsat);
RimActiveSUVR = reformStruct(RimActiveSUVR); 
ActiveSUVR = reformStruct(ActiveSUVR); 
InactiveSUVR = reformStruct(InactiveSUVR);
RimActivePerCore = reformStruct(RimActivePerCore); 
ActivePerCore = reformStruct(ActivePerCore); 

save([outDir '/WMLsMeasures.mat'],'WMLsMTsat', 'WMLsSUVR', 'WMLsVol',...
    'RimActive', 'Active', 'Inactive', 'WMLscount', 'RimActiveVol', 'ActiveVol', 'InactiveVol',...
    'RimActiveMTsat', 'ActiveMTsat', 'InactiveMTsat', 'RimActiveSUVR', 'ActiveSUVR', 'InactiveSUVR',...
    'RimActivePerCore', 'ActivePerCore');

% Display the total number of lesions in each category for all patients
fprintf('Total number of lesions in watershed regions only: %d\n', sum(cellfun(@numel, WMLsMTsat.watershed)));
fprintf('Total number of lesions in mixed regions (watershed + non-watershed): %d\n', sum(cellfun(@numel, WMLsMTsat.mixed)));
fprintf('Total number of lesions in non-watershed regions only: %d\n', sum(cellfun(@numel, WMLsMTsat.nonWatershed)));
toc

%% 6) Output TSPO PET SUVR maps masked only for lesion + lesional rims 
clear; close all; clc; 
tic
addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/infoPET.mat"); 

DirsMS = DirsPET(string(groupsPET)=="MS"); 
NamesMS = NamesPET(string(groupsPET)=="MS"); 
noMS = length(DirsMS); 

% Iterate through each patient's directory 
for id = 1:noMS
    mainDir = DirsMS{id}; 
    subj7TDir = [DirsMS{id} '/' NamesMS{id}]; 

    name = [subj7TDir '/WM_lesions.nii'];  
    WMLs = load_untouch_nii(name); 
    DimDat = WMLs.hdr.dime.pixdim; Hist = WMLs.hdr.hist; WMLs = double(WMLs.img); 
%     name = [subj7TDir '/cerebralNAWM.nii'];  
%     NAWM = double(load_untouch_nii(name).img); 
%     NAWM(WMLs==1)=1;  

    name = sprintf("%s/ResampledInMNI152Space/WM_lesions_TSPOPhenotype.nii",subj7TDir);
    TSPOWMLphenotype = double(load_untouch_nii(char(name)).img); 
    TSPOWMLphenotype(TSPOWMLphenotype~=0) = 1; 
    
    LV = double(load_untouch_nii([subj7TDir '/LV_geo_seg_output/LateralVentricles_mask.nii']).img);
    
    se = strel('sphere',1);
    %Dilate lesion ('lesion core') by 3 voxels to get 2.1mm perilesional rim
    dilatedMask = imdilate(TSPOWMLphenotype,se); %dilate lesion edge by approx. 1 voxel 
    dilatedMask = imdilate(dilatedMask,se); %dilate lesion edge by approx. 1 voxel 
    dilatedMask = imdilate(dilatedMask,se); %dilate lesion edge by approx. 1 voxel

    name = sprintf("%s/PETOutput/PET_SUVR_T1w_reg.nii",subj7TDir);
    SUVR = double(load_untouch_nii(char(name)).img); 
    SUVR(dilatedMask==0) = 0; SUVR(LV~=0) = 0;  
    name = sprintf("%s/ResampledInMNI152Space/SUVR_WMLcorerims_masked.nii",subj7TDir);
    imag_nii = make_nii(double(SUVR)); 
    if(~isempty(DimDat)), imag_nii.hdr.dime.pixdim = DimDat; end 
    if(~isempty(Hist)),imag_nii.hdr.hist = Hist; end 
    save_nii(imag_nii,char(name));
    
    %blurring resampled PET suvr image using 3D geodesic kernel of given FWHM
    imgPath = [subj7TDir '/PETOutput/PET_SUVR_T1w_reg.nii']; FWHM = 3; 
%     SPM12smooth(imgPath,FWHM,sprintf("blur%dmm_",FWHM)); 
    name = sprintf("%s/PETOutput/blur%.1fmm_PET_SUVR_T1w_reg.nii",subj7TDir,FWHM);
    SUVR = double(load_untouch_nii(char(name)).img); 
    SUVR(dilatedMask==0) = 0; SUVR(LV~=0) = 0; 
    name = sprintf("%s/ResampledInMNI152Space/blurSUVR_WMLcorerims_masked.nii",subj7TDir);
    imag_nii = make_nii(double(SUVR)); 
    if(~isempty(DimDat)), imag_nii.hdr.dime.pixdim = DimDat; end 
    if(~isempty(Hist)),imag_nii.hdr.hist = Hist; end 
    save_nii(imag_nii,char(name));    
end

%% Generate "baseline" HC data for MS patient WM lesion measures
clear; close all; clc; 
addpath(genpath("/export02/data/risa/NIfTI_20140122/")); 
load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/infoPET.mat"); 
conn = 6; %CHECK

DirsHC = DirsPET(string(groupsPET)=="HC"); 
NamesHC = NamesPET(string(groupsPET)=="HC"); 
TIVHC = TIV(string(groupsPET)=="HC");
noHC = length(DirsHC);

% Initialize a structure to store the MTsat values for each patient and lesion classification
cWMLsMTsat = struct();
cWMLsMTsat.watershed = cell(1, noHC);     % Watershed lesions per patient
cWMLsMTsat.mixed = cell(1, noHC);         % Mixed lesions (watershed + non-watershed) per patient
cWMLsMTsat.nonWatershed = cell(1, noHC);  % Non-watershed lesions per patient
WMMTsat = zeros(1, noHC); 

% Similarly,...
cWMLsSUVR = cWMLsMTsat; cWMLsVol = cWMLsMTsat; WMSUVR = WMMTsat; WMVol = WMMTsat; 

% Iterate through each control's directory 
for id = 1:noHC
    mainDir = DirsHC{id}; 
    subj7TDir = [DirsHC{id} '/' NamesHC{id}]; 

    %Load data volumes 
    MTsatVolume = double(load_untouch_nii([subj7TDir '/MTSat/MTSat.nii']).img);
    name = sprintf("%s/PETOutput/PET_SUVR_T1w_reg.nii",subj7TDir); 
%     name = sprintf("%s/PETOutput/PET_SUVR_T1w_reg_noPVC.nii",subj7TDir); %TEMPORARY
    SUVRVolume = double(load_untouch_nii(char(name)).img); 
    
    % Whole cerebral WM specific measures 
    WM = double(load_untouch_nii([subj7TDir '/cerebralNAWM.nii']).img);
    MTsat = MTsatVolume(WM==1); MTsat(MTsat<0 | MTsat>5) = []; 
    SUVR = SUVRVolume(WM==1);  
    WMMTsat(id) = mean(MTsat,'omitnan'); WMSUVR(id) = mean(SUVR,'omitnan'); 
    
    % 1) Counter ws-specific WM lesions 
    % Temporary variables for the current patient
    MTsatarr = []; SUVRarr = []; Volarr = [];   
    name = [subj7TDir '/ResampledInMNI152Space/CounterWMLs_ws.nii'];  
    cWMLs = load_untouch_nii(name); DimDat = cWMLs.hdr.dime.pixdim; cWMLs = double(cWMLs.img);  
    %Exclusion criterion 1: Omit Counter WM lesion if having volume <= 27mm3
    CC = bwconncomp(cWMLs,conn); dim = DimDat(2:4);  
    for i = 1:length(CC.PixelIdxList)
        volLesion = dim(1)*dim(2)*dim(3)*length(CC.PixelIdxList{i}); 
        if (volLesion <= 27) 
            cWMLs(CC.PixelIdxList{i}) = 0; 
        end 
    end
    CC = bwconncomp(cWMLs,conn);     
    for i = 1:length(CC.PixelIdxList) % Iterate through each counter WM lesion
        %Logical array specifying current counter WM lesion
        lesionMask = zeros(size(cWMLs)); 
        lesionMask(CC.PixelIdxList{i}) = 1;  
        
        lesionMask = lesionMask==1; 
        
        % Extract measures 
        MTsat = MTsatVolume(lesionMask); MTsat(MTsat<0 | MTsat>5) = []; 
        SUVR = SUVRVolume(lesionMask);  
        volLesion = dim(1)*dim(2)*dim(3)*length(CC.PixelIdxList{i})/1000; %mm3 to cm3
%         volLesion = volLesion*100/TIVHC(id); %Normalize WML volume by TIV
      
        MTsatarr = [MTsatarr; mean(MTsat,'omitnan')]; 
        SUVRarr = [SUVRarr; mean(SUVR,'omitnan')]; 
        Volarr = [Volarr; volLesion]; 
    end
    % Store the results for this patient in the structure
    cWMLsMTsat.watershed{id} = MTsatarr;
    cWMLsSUVR.watershed{id} = SUVRarr;
    cWMLsVol.watershed{id} = Volarr; 
        
    % 2) Counter mixed-specific WM lesions 
    % Temporary variables for the current patient
    MTsatarr = []; SUVRarr = []; Volarr = []; 
    name = [subj7TDir '/ResampledInMNI152Space/CounterWMLs_mixed.nii'];  
    cWMLs = load_untouch_nii(name); DimDat = cWMLs.hdr.dime.pixdim; cWMLs = double(cWMLs.img);  
    %Exclusion criterion 1: Omit Counter WM lesion if having volume <= 27mm3
    CC = bwconncomp(cWMLs,conn); dim = DimDat(2:4);  
    for i = 1:length(CC.PixelIdxList)
        volLesion = dim(1)*dim(2)*dim(3)*length(CC.PixelIdxList{i}); 
        if (volLesion <= 27) 
            cWMLs(CC.PixelIdxList{i}) = 0; 
        end 
    end
    CC = bwconncomp(cWMLs,conn);     
    for i = 1:length(CC.PixelIdxList) % Iterate through each counter WM lesion
        %Logical array specifying current counter WM lesion
        lesionMask = zeros(size(cWMLs)); 
        lesionMask(CC.PixelIdxList{i}) = 1;  
        
        lesionMask = lesionMask==1; 
        
        % Extract measures 
        MTsat = MTsatVolume(lesionMask); MTsat(MTsat<0 | MTsat>5) = []; 
        SUVR = SUVRVolume(lesionMask);  
        volLesion = dim(1)*dim(2)*dim(3)*length(CC.PixelIdxList{i})/1000; %mm3 to cm3
%         volLesion = volLesion*100/TIVHC(id); %Normalize WML volume by TIV
      
        MTsatarr = [MTsatarr; mean(MTsat,'omitnan')]; 
        SUVRarr = [SUVRarr; mean(SUVR,'omitnan')]; 
        Volarr = [Volarr; volLesion];
    end
    % Store the results for this patient in the structure
    cWMLsMTsat.mixed{id} = MTsatarr;
    cWMLsSUVR.mixed{id} = SUVRarr;
    cWMLsVol.mixed{id} = Volarr; 

    % 3) Counter nws-specific WM lesions 
    % Temporary variables for the current patient
    MTsatarr = []; SUVRarr = []; Volarr = []; 
    name = [subj7TDir '/ResampledInMNI152Space/CounterWMLs_nws.nii'];  
    cWMLs = load_untouch_nii(name); DimDat = cWMLs.hdr.dime.pixdim; cWMLs = double(cWMLs.img);  
    %Exclusion criterion 1: Omit Counter WM lesion if having volume <= 27mm3
    CC = bwconncomp(cWMLs,conn); dim = DimDat(2:4);  
    for i = 1:length(CC.PixelIdxList)
        volLesion = dim(1)*dim(2)*dim(3)*length(CC.PixelIdxList{i}); 
        if (volLesion <= 27) 
            cWMLs(CC.PixelIdxList{i}) = 0; 
        end 
    end
    CC = bwconncomp(cWMLs,conn);     
    for i = 1:length(CC.PixelIdxList) % Iterate through each counter WM lesion
        %Logical array specifying current counter WM lesion
        lesionMask = zeros(size(cWMLs)); 
        lesionMask(CC.PixelIdxList{i}) = 1;  
        
        lesionMask = lesionMask==1; 
        
        % Extract measures 
        MTsat = MTsatVolume(lesionMask); MTsat(MTsat<0 | MTsat>5) = []; 
        SUVR = SUVRVolume(lesionMask);  
        volLesion = dim(1)*dim(2)*dim(3)*length(CC.PixelIdxList{i})/1000; %mm3 to cm3
%         volLesion = volLesion*100/TIVHC(id); %Normalize WML volume by TIV
      
        MTsatarr = [MTsatarr; mean(MTsat,'omitnan')]; 
        SUVRarr = [SUVRarr; mean(SUVR,'omitnan')]; 
        Volarr = [Volarr; volLesion];
    end
    % Store the results for this patient in the structure
    cWMLsMTsat.nonWatershed{id} = MTsatarr;
    cWMLsSUVR.nonWatershed{id} = SUVRarr;
    cWMLsVol.nonWatershed{id} = Volarr; 
    
    WMVol(id)= dim(1)*dim(2)*dim(3)*sum(WM==1,'all')/1000; %mm3 to cm3
end

% Save the structure in a single MAT file
outDir = '/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis'; 
save([outDir '/counterWMLsMeasures.mat'],'cWMLsMTsat', 'WMMTsat', 'cWMLsSUVR', 'WMSUVR', 'cWMLsVol', 'WMVol'); 
 
%% addWMLdata() function 
function [MTsatwsl, MTsatmixedl, MTsatnwsl,...
    SUVRwsl, SUVRmixedl, SUVRnwsl,...
    Volwsl, Volmixedl, Volnwsl,...
    noInactive, noActive, noRimActive, noWMLs,...
    IAVol, AVol, RAVol,...
    IAMTsat, AMTsat, RAMTsat,...
    IASUVR, ASUVR, RASUVR,...
    APerCore, RAPerCore,...
    WMLvascularlbl] = ...
    addWMLdata(subj7TDir, lesionMask, minLesionVol,...
    MTsatwsl, MTsatmixedl, MTsatnwsl,...
    SUVRwsl, SUVRmixedl, SUVRnwsl,...
    Volwsl, Volmixedl, Volnwsl,...
    noInactive, noActive, noRimActive, noWMLs,...
    IAVol, AVol, RAVol,...
    IAMTsat, AMTsat, RAMTsat,...  
    IASUVR, ASUVR, RASUVR,...
    APerCore, RAPerCore,...
    WMLvascularlbl, val, Perc_ActiveCore, dim)

% 1) Load the shared region mask (watershed and non-watershed regions)
name = sprintf("%s/ResampledInMNI152Space/watershed_mask.nii",subj7TDir);
regionMask = double(load_untouch_nii(char(name)).img);  

% 2) Load data volumes
MTsatVolume = double(load_untouch_nii([subj7TDir '/MTSat/MTSat.nii']).img);
name = sprintf("%s/PETOutput/PET_SUVR_T1w_reg.nii",subj7TDir); 
% name = sprintf("%s/PETOutput/PET_SUVR_T1w_reg_noPVC.nii",subj7TDir); %TEMPORARY
SUVRVolume = double(load_untouch_nii(char(name)).img); 
WMLs = double(load_untouch_nii([subj7TDir '/WM_lesions.nii']).img); WMLs(WMLs~=0)=1; 
WM = double(load_untouch_nii([subj7TDir '/cerebralNAWM.nii']).img); WM(WMLs==1) = 1; 

% 4) Calculate lesion volume and see if it meets minimum lesion volume (> 27mm3?)
volLesion = dim(1)*dim(2)*dim(3)*sum(lesionMask(:),'omitnan')/1000; %mm3 to cm3
%If volume of current lesion is <=27mm3, omit it 
if (volLesion <= minLesionVol), return; end %CRITERION
                
% 5) Nylund (2022) definition for TSPO PET-based phenotyping 
se = strel('sphere',1);
%Dilate lesion ('lesion core') by 3 voxels to get 2.1mm perilesional rim
lesionDilated = imdilate(lesionMask,se); %dilate lesion edge by approx. 1 voxel 
lesionDilated = imdilate(lesionDilated,se); %dilate lesion edge by approx. 1 voxel 
lesionDilated = imdilate(lesionDilated,se); %dilate lesion edge by approx. 1 voxel

lesionMask = lesionMask==1;

% 6) Extract measures for this lesion
MTsat = MTsatVolume(lesionMask); MTsat(MTsat<0 | MTsat>5) = []; 
SUVR = SUVRVolume(lesionMask); 

% 7) Check if the lesion overlaps with watershed or non-watershed regions
lesionInWatershed = any(regionMask(lesionMask) == 1);
lesionInNonWatershed = any(regionMask(lesionMask) == 0);
        
% 9) Classify the lesion and store the MTsat mean values
if lesionInWatershed && ~lesionInNonWatershed
    % Lesion falls in watershed regions only
    MTsatwsl = [MTsatwsl; mean(MTsat,'omitnan')];
    SUVRwsl = [SUVRwsl; mean(SUVR,'omitnan')]; 
    Volwsl = [Volwsl; volLesion]; 

    WMLvascularlbl(lesionMask) = 1; %ws classification
    noWMLs(1) = noWMLs(1) + 1; 
    if (val==1) 
        noRimActive(1) = noRimActive(1) + 1; 
        RAVol{1} = [RAVol{1}; volLesion]; 
        RAMTsat{1} = [RAMTsat{1}; mean(MTsat,'omitnan')];
        RASUVR{1} = [RASUVR{1}; mean(SUVR,'omitnan')]; 
        RAPerCore{1} = [RAPerCore{1}; Perc_ActiveCore]; 
    elseif (val==2)
        noInactive(1) = noInactive(1) + 1; 
        IAVol{1} = [IAVol{1}; volLesion]; 
        IAMTsat{1} = [IAMTsat{1}; mean(MTsat,'omitnan')];
        IASUVR{1} = [IASUVR{1}; mean(SUVR,'omitnan')]; 
    elseif (val==3)
        noActive(1) = noActive(1) + 1; 
        AVol{1} = [AVol{1}; volLesion]; 
        AMTsat{1} = [AMTsat{1}; mean(MTsat,'omitnan')];
        ASUVR{1} = [ASUVR{1}; mean(SUVR,'omitnan')];
        APerCore{1} = [APerCore{1}; Perc_ActiveCore]; 
    end 
elseif lesionInWatershed && lesionInNonWatershed
    % Lesion falls in both watershed and non-watershed regions
    
    %TEMPORARY: Omit Mixed Lesions with greater than 50% volume in watershed zones%%%%%%%%%%%%%%%% 
    tmpVal = regionMask(lesionMask) == 1; 
    volLesionws = dim(1)*dim(2)*dim(3)*sum(tmpVal(:),'omitnan')/1000; %mm3 to cm3
    if (volLesionws*100/volLesion > 50), return; end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MTsatmixedl = [MTsatmixedl; mean(MTsat,'omitnan')];
    SUVRmixedl = [SUVRmixedl; mean(SUVR,'omitnan')];
    Volmixedl = [Volmixedl; volLesion]; 

    WMLvascularlbl(lesionMask) = 2; %mixed classification
    noWMLs(2) = noWMLs(2) + 1;
    if (val==1) 
        noRimActive(2) = noRimActive(2) + 1; 
        RAVol{2} = [RAVol{2}; volLesion];
        RAMTsat{2} = [RAMTsat{2}; mean(MTsat,'omitnan')];
        RASUVR{2} = [RASUVR{2}; mean(SUVR,'omitnan')]; 
        RAPerCore{2} = [RAPerCore{2}; Perc_ActiveCore]; 
    elseif (val==2)
        noInactive(2) = noInactive(2) + 1;
        IAVol{2} = [IAVol{2}; volLesion]; 
        IAMTsat{2} = [IAMTsat{2}; mean(MTsat,'omitnan')];
        IASUVR{2} = [IASUVR{2}; mean(SUVR,'omitnan')]; 
    elseif (val==3)
        noActive(2) = noActive(2) + 1;
        AVol{2} = [AVol{2}; volLesion];
        AMTsat{2} = [AMTsat{2}; mean(MTsat,'omitnan')];
        ASUVR{2} = [ASUVR{2}; mean(SUVR,'omitnan')]; 
        APerCore{2} = [APerCore{2}; Perc_ActiveCore];
    end 
elseif ~lesionInWatershed && lesionInNonWatershed
    % Lesion falls in non-watershed regions only
    MTsatnwsl = [MTsatnwsl; mean(MTsat,'omitnan')];
    SUVRnwsl = [SUVRnwsl; mean(SUVR,'omitnan')]; 
    Volnwsl = [Volnwsl; volLesion]; 

    WMLvascularlbl(lesionMask) = 3; %nws classification
    noWMLs(3) = noWMLs(3) + 1;
    if (val==1) 
        noRimActive(3) = noRimActive(3) + 1; 
        RAVol{3} = [RAVol{3}; volLesion];
        RAMTsat{3} = [RAMTsat{3}; mean(MTsat,'omitnan')];
        RASUVR{3} = [RASUVR{3}; mean(SUVR,'omitnan')]; 
        RAPerCore{3} = [RAPerCore{3}; Perc_ActiveCore];
    elseif (val==2)
        noInactive(3) = noInactive(3) + 1; 
        IAVol{3} = [IAVol{3}; volLesion]; 
        IAMTsat{3} = [IAMTsat{3}; mean(MTsat,'omitnan')];
        IASUVR{3} = [IASUVR{3}; mean(SUVR,'omitnan')]; 
    elseif (val==3)
        noActive(3) = noActive(3) + 1; 
        AVol{3} = [AVol{3}; volLesion];
        AMTsat{3} = [AMTsat{3}; mean(MTsat,'omitnan')];
        ASUVR{3} = [ASUVR{3}; mean(SUVR,'omitnan')]; 
        APerCore{3} = [APerCore{3}; Perc_ActiveCore];
    end 
end

end 

%% Perform TSPO phenotyping for all MS patients (Purely Nylund 2020 method)
%But not separating according to TSPO binding affinity unlike Nylund 2020,
%Hamzaoui 2023 and Poirion 2021

function [TSPOActClass, val, percActive_core] = TSPOPhenotypeNylund(subj7TDir, bindAff, lesionMask, lesionMaskCU,...
    minLesionVol, TSPOActClass, dim)

addpath(genpath("/export02/data/risa/NIfTI_20140122/"));
if (bindAff=="HAB")
    load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/HC_HABdata_WMLphenotyping.mat"); 
elseif (bindAff=="MAB") 
    load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/HC_MABdata_WMLphenotyping.mat");
else 
    load("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/HC_LABdata_WMLphenotyping.mat"); 
end 
%Use threshold from 2.5mm Gaussian kernel blurred PET SUVR map
%NOTE: we are using 2.5mm blurred images and threshold calculated from
%those images, to enhance detectability of "active voxel CLUSTERS" and
%reduce noise, like the Airas group does 
thresh = threshSUVRblur; 

% 1) Calculate lesion volume and see if it meets minimum lesion volume (> 27mm3?)
volLesion = dim(1)*dim(2)*dim(3)*sum(lesionMask(:),'omitnan')/1000; %mm3 to cm3
%If volume of current lesion is <=27mm3, omit it 
if (volLesion <= minLesionVol), val = 0; percActive_core = 0; return; end %CRITERION

WMLs = double(load_untouch_nii([subj7TDir '/WM_lesions.nii']).img); WMLs(WMLs~=0)=1; 
WM = double(load_untouch_nii([subj7TDir '/cerebralNAWM.nii']).img); WM(WMLs==1) = 1; 
                
% 2) Nylund (2022) definition for TSPO PET-based phenotyping 
se = strel('sphere',1);
%Dilate lesion ('lesion core') by 3 voxels to get 2.1mm perilesional rim
lesionDilated = imdilate(lesionMask,se); %dilate lesion edge by approx. 1 voxel 
lesionDilated = imdilate(lesionDilated,se); %dilate lesion edge by approx. 1 voxel 
lesionDilated = imdilate(lesionDilated,se); %dilate lesion edge by approx. 1 voxel

% 3) Define perilesional rim 
rim = lesionDilated; rim(lesionMask==1) = 0; %get rim
rim(WM==0) = 0; %remove non-WM voxels from rim 
rim(lesionMaskCU==1) = 0; %remove shared part of rim (if part of confluent unit)
LV = double(load_untouch_nii([subj7TDir '/LV_geo_seg_output/LateralVentricles_mask.nii']).img);
rim(LV==1) = 0; %remove LV voxels from rim 

lesionMask = lesionMask==1; lesionRimMask = rim==1; %logical arrays
% 4) Calculate rim volume and see if it meets minimum lesion volume (> 27mm3?)
volRim = dim(1)*dim(2)*dim(3)*sum(lesionRimMask(:),'omitnan')/1000; %mm3 to cm3
%If volume of current rim is <=27mm3, omit it 
if (volRim <= 0.027), val = 0; percActive_core = 0; return; end %CRITERION

% name = char(sprintf("%s/PETOutput/PET_SUVR_T1w_reg.nii",subj7TDir));
%NOTE: we are using 2.5mm blurred images and threshold calculated from
%those images, to enhance detectability of "active voxel CLUSTERS" and
%reduce noise, like the Airas group does
name = char(sprintf("%s/PETOutput/blur2.5mm_PET_SUVR_T1w_reg.nii",subj7TDir)); 
% name = char(sprintf("%s/PETOutput/blur2.5mm_PET_SUVR_T1w_reg_noPVC.nii",subj7TDir)); %TEMPORARY
SUVR = double(load_untouch_nii(name).img); %need to associate WML rim and core 
SUVR_active_mask = SUVR>=thresh; %general mask, including all tissue compartments 
CC = bwconncomp(SUVR_active_mask, 6); 
for i = 1:length(CC.PixelIdxList)
    if length(CC.PixelIdxList{i}) < 8 %omit active clusters of less than 8 voxels to remove random peak values
        SUVR_active_mask(CC.PixelIdxList{i}) = 0; 
    end     
end 

%Percentage proportion of voxels in lesional rim that are TSPO-active
active_rim_idxs = lesionRimMask & SUVR_active_mask; 
noPx_active_rim = sum(lesionRimMask(active_rim_idxs),'all','omitnan'); 
noPx_rim = sum(lesionRimMask,'all','omitnan'); 
percActive_rim = 100*noPx_active_rim/noPx_rim; 

%Percentage proportion of voxels in lesional core that are TSPO-active
active_core_idxs = lesionMask & SUVR_active_mask;    
noPx_active_core= sum(lesionMask(active_core_idxs),'all','omitnan'); 
noPx_core = sum(lesionMask,'all','omitnan'); 
percActive_core = 100*noPx_active_core/noPx_core;

if (percActive_core < 5.0) && (percActive_rim >= (percActive_core+5.0))
    val = 1; %rim-active
elseif (percActive_core <= 20) && (percActive_core >= 5.0) && (percActive_rim >= 2*percActive_core)
    val = 1; %rim-active
elseif (percActive_core == 0) %&& (percActive_rim == 0) #TEMPORARY
    val = 2; %inactive
else %There are active voxels in the lesion core, along with or NOT in lesion rim
    val = 3; %overall-active 
end 
TSPOActClass(lesionMask) = val; 

end 

%% cell of structs to struct of cells 
function [output] = reformStruct(dat)

tmp = cell2mat(dat); 
output.watershed = {tmp.watershed}; 
output.mixed = {tmp.mixed}; 
output.nonWatershed = {tmp.nonWatershed}; 

end 
