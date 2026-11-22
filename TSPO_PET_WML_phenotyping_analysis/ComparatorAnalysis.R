library(lme4)
library(lmerTest)
library(ggsignif)
library(emmeans)
library(ggprism)
library(ggpubr)

mainDir <- "/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/"
outDir <- "/export02/data/risa/DODManuscript2_Revisions/"

#FOR FORMATTING DATA FRAMES
idxsHC <- explaVars$subj.type=="HC" 
subject.name <- unlist(explaVars$subj.name[idxsHC]) 
clinical.status <- as.character(unlist(explaVars$clinicalStatus[idxsHC])); names(clinical.status) <- subject.name
sexList <- as.character(unlist(explaVars$sex[idxsHC])); names(sexList) <- subject.name
ageList <- unlist(explaVars$age[idxsHC]); names(ageList) <- subject.name
bindAffList <- as.character(unlist(explaVars$bindAff[idxsHC])); names(bindAffList) <- subject.name
DMTeffList <- as.character(unlist(explaVars$DMT.effect[idxsHC])); names(DMTeffList) <- subject.name
edssList <- unlist(explaVars$edss[idxsHC]); names(edssList) <- subject.name
gctList <- unlist(explaVars$gct[idxsHC]); names(gctList) <- subject.name
# meanMTsatList <- unlist(explaVars$mean.cxMTsat[idxsHC]); names(meanMTsatList) <- subject.name #CHANGE
disDurList <- unlist(explaVars$diseaseDuration[idxsHC]); names(disDurList) <- subject.name
ageOnsetList <- unlist(explaVars$age.at.onset[idxsHC]); names(ageOnsetList) <- subject.name
WBVList <- unlist(explaVars$WBV[idxsHC]); names(WBVList) <- subject.name
TIVList <- unlist(explaVars$TIV[idxsHC]); names(TIVList) <- subject.name
WMLTotalVolList <- unlist(explaVars$WMLTotalVol[idxsHC]); names(WMLTotalVolList) <- subject.name
pyramidalList <- unlist(explaVars$pyramidal[idxsHC]); names(pyramidalList) <- subject.name
cerebellarList <- unlist(explaVars$cerebellar[idxsHC]); names(cerebellarList) <- subject.name
brainstemList <- unlist(explaVars$brainstem[idxsHC]); names(brainstemList) <- subject.name
sensoryList <- unlist(explaVars$sensory[idxsHC]); names(sensoryList) <- subject.name
bbList <- unlist(explaVars$bowel.bladder[idxsHC]); names(bbList) <- subject.name
visualList <- unlist(explaVars$visual[idxsHC]); names(visualList) <- subject.name
mentalList <- unlist(explaVars$mental[idxsHC]); names(mentalList) <- subject.name
format.WMLs.df <- function(WMLs.temp)
{ 
  WMLs.temp <- within(WMLs.temp, region <- relevel(factor(region), ref = "ws"))
  WMLs.temp <- WMLs.temp %>% mutate(age=dplyr::recode(subject.name,!!!ageList))
  WMLs.temp <- WMLs.temp %>% mutate(sex=dplyr::recode(subject.name,!!!sexList))
  WMLs.temp <- WMLs.temp %>% mutate(clinical.status=dplyr::recode(subject.name,!!!clinical.status))
  WMLs.temp <- WMLs.temp %>% mutate(bindAff=dplyr::recode(subject.name,!!!bindAffList))
  WMLs.temp <- WMLs.temp %>% mutate(mental=dplyr::recode(subject.name,!!!mentalList))
  WMLs.temp <- WMLs.temp %>% mutate(edss=dplyr::recode(subject.name,!!!edssList))
  WMLs.temp <- WMLs.temp %>% mutate(gct=dplyr::recode(subject.name,!!!gctList)) 
  # WMLs.temp <- WMLs.temp %>% mutate(meanø.cxMTsat=dplyr::recode(subject.name,!!!meanMTsatList)) #CHANGE
  WMLs.temp <- WMLs.temp %>% mutate(disDur=dplyr::recode(subject.name,!!!disDurList))
  WMLs.temp <- WMLs.temp %>% mutate(ageOnset=dplyr::recode(subject.name,!!!ageOnsetList))
  WMLs.temp <- WMLs.temp %>% mutate(WBV=dplyr::recode(subject.name,!!!WBVList))
  WMLs.temp <- WMLs.temp %>% mutate(TIV=dplyr::recode(subject.name,!!!TIVList))
  WMLs.temp <- WMLs.temp %>% mutate(WMLTotalVol=dplyr::recode(subject.name,!!!WMLTotalVolList))
  WMLs.temp <- WMLs.temp %>% mutate(pyramidal=dplyr::recode(subject.name,!!!pyramidalList))
  WMLs.temp <- WMLs.temp %>% mutate(cerebellar=dplyr::recode(subject.name,!!!cerebellarList))
  WMLs.temp <- WMLs.temp %>% mutate(brainstem=dplyr::recode(subject.name,!!!brainstemList))
  WMLs.temp <- WMLs.temp %>% mutate(sensory=dplyr::recode(subject.name,!!!sensoryList))
  WMLs.temp <- WMLs.temp %>% mutate(bb=dplyr::recode(subject.name,!!!bbList))
  WMLs.temp <- WMLs.temp %>% mutate(visual=dplyr::recode(subject.name,!!!visualList))
  WMLs.temp <- WMLs.temp %>% mutate(mental=dplyr::recode(subject.name,!!!mentalList))
  WMLs.temp <- WMLs.temp %>% mutate(DMTeff=dplyr::recode(subject.name,!!!DMTeffList))
  
  return(WMLs.temp)
}

# WM LESIONS -------------------------------------------------------------- -------------------------------------------------------------
# volume at individual WML level ----------------------
dat <- readMat(sprintf("%scounterWMLsMeasures.mat",mainDir))

ws <- map(dat$cWMLsVol[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.Vol <- stack(tmp)
ws.Vol <- cbind(ws.Vol,rep("ws",nrow(ws.Vol)))
colnames(ws.Vol) <- c("Vol", "subject.name", "region")

mixed <- map(dat$cWMLsVol[[2]],unlist);
tmp = list()
tmp[subject.name] <- mixed
mixed.Vol <- stack(tmp)
mixed.Vol <- cbind(mixed.Vol,rep("mixed",nrow(mixed.Vol)))
colnames(mixed.Vol) <- c("Vol", "subject.name", "region")

nws <- map(dat$cWMLsVol[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.Vol <- stack(tmp)
nws.Vol <- cbind(nws.Vol,rep("nws",nrow(nws.Vol)))
colnames(nws.Vol) <- c("Vol", "subject.name", "region")

cWMLs.Vol <- rbind(ws.Vol,mixed.Vol,nws.Vol)
cWMLs.Vol <- format.WMLs.df(cWMLs.Vol)
cWMLs.Vol$logVol <- log(cWMLs.Vol$Vol)

#Volume of aggregated binary lesion mask matched to each HC's MRI space
cWMLs.TotalVol <- cWMLs.Vol %>%
  dplyr::group_by(subject.name) %>%
  dplyr::summarise(sumVol = sum(Vol), WBV = max(WBV), TIV = max(TIV))

p3 <- ggplot(cWMLs.Vol, aes(region, logVol)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimWMLVol)

mmod3 <- lmer(logVol~(1|subject.name)+region+age, cWMLs.Vol,REML=TRUE)
summary(mmod3)
step(mmod3, direction="backward")
anova(mmod3)
#Fixed effect term 'age' and random effect term '(1 | subject.name:region)'
#found to be insignificant. Omitting them. 

#Final model 
mmod3 <- lmer(logVol ~ region + (1 | subject.name), cWMLs.Vol,REML=TRUE)
summary(mmod3)
tmp <- emmeans_test(cWMLs.Vol,logVol~region,p.adjust.method="fdr",model=mmod3)
ggsave(sprintf("%sModel3_plot.jpg",outDir),p3,width=16,height=16,units="cm",dpi=1200)

# mean MTsat at individual WML level -------------------------
dat <- readMat(sprintf("%scounterWMLsMeasures.mat",mainDir))

ws <- map(dat$cWMLsMTsat[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.mtsat <- stack(tmp)
ws.mtsat <- cbind(ws.mtsat,rep("ws",nrow(ws.mtsat)))
colnames(ws.mtsat) <- c("mtsat", "subject.name", "region")

mixed <- map(dat$cWMLsMTsat[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.mtsat <- stack(tmp)
mixed.mtsat <- cbind(mixed.mtsat,rep("mixed",nrow(mixed.mtsat)))
colnames(mixed.mtsat) <- c("mtsat", "subject.name", "region")

nws <- map(dat$cWMLsMTsat[[3]],unlist);
tmp = list()
tmp[subject.name] <- nws
nws.mtsat <- stack(tmp)
nws.mtsat <- cbind(nws.mtsat,rep("nws",nrow(nws.mtsat)))
colnames(nws.mtsat) <- c("mtsat", "subject.name", "region")

cWMLs.mtsat <- rbind(ws.mtsat,mixed.mtsat,nws.mtsat)
cWMLs.mtsat <- format.WMLs.df(cWMLs.mtsat)
cWMLs.mtsat <- cWMLs.mtsat %>% mutate(WML.volume=cWMLs.Vol$Vol)
cWMLs.mtsat$mtsatpervol <- cWMLs.mtsat$mtsat/(cWMLs.mtsat$WML.volume*1000)

p1 <- ggplot(cWMLs.mtsat, aes(region, mtsat)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimMTsat)

mmod1 <- lmer(mtsat~1+(1|subject.name)+region+age+WML.volume,cWMLs.mtsat,REML=TRUE)
summary(mmod1)
step(mmod1, direction="backward")
#Fixed effect terms 'age' and 'region' found to be insignificant, but keeping 
#'region' as regressor

#Final model 
mmod1 <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region ,cWMLs.mtsat,REML=TRUE)
summary(mmod1)
tmp <- emmeans_test(cWMLs.mtsat,mtsat~region,p.adjust.method="fdr",model=mmod1)
ggsave(sprintf("%sModel1_plot.jpg",outDir),p1,width=16,height=16,units="cm",dpi=1200)

# mean SUVR at individual WML level ----------------------
dat <- readMat(sprintf("%scounterWMLsMeasures.mat",mainDir))

ws <- map(dat$cWMLsSUVR[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws 
ws.SUVR <- stack(tmp)
ws.SUVR <- cbind(ws.SUVR,rep("ws",nrow(ws.SUVR)))
colnames(ws.SUVR) <- c("SUVR", "subject.name", "region")

mixed <- map(dat$cWMLsSUVR[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.SUVR <- stack(tmp)
mixed.SUVR <- cbind(mixed.SUVR,rep("mixed",nrow(mixed.SUVR)))
colnames(mixed.SUVR) <- c("SUVR", "subject.name", "region")

nws <- map(dat$cWMLsSUVR[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.SUVR <- stack(tmp)
nws.SUVR <- cbind(nws.SUVR,rep("nws",nrow(nws.SUVR)))
colnames(nws.SUVR) <- c("SUVR", "subject.name", "region")

cWMLs.SUVR <- rbind(ws.SUVR,mixed.SUVR,nws.SUVR)
cWMLs.SUVR <- format.WMLs.df(cWMLs.SUVR)
cWMLs.SUVR <- cWMLs.SUVR %>% mutate(WML.volume=cWMLs.Vol$Vol)

p2 <- ggplot(cWMLs.SUVR, aes(region, SUVR)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimSUVR)

mmod2 <- lmer(SUVR~1+(1|subject.name)+region+age+WML.volume+bindAff,cWMLs.SUVR,REML=TRUE)
summary(mmod2)
step(mmod2, direction="backward")
#Fixed effect term 'age' found to be insignificant 

#Final model 
mmod2 <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume+bindAff,cWMLs.SUVR,REML=TRUE)
summary(mmod2)
tmp <- emmeans_test(cWMLs.SUVR,SUVR~region,p.adjust.method="fdr",model=mmod2)
ggsave(sprintf("%sModel2_plot.jpg",outDir),p2,width=16,height=16,units="cm",dpi=1200)

# EMPTY SPACE -------------------------------------------------------------
# Comparison for WML MTsat  --------------------------------------------------------------
cWMLs.mtsat <- cWMLs.mtsat %>% mutate(subj.type="HC")
WMLs.mtsat <- WMLs.mtsat %>% mutate(subj.type="MS")
MTsat <- rbind(cWMLs.mtsat,WMLs.mtsat)

#Test for normality: paired samples t-test or paired samples wilcoxon test, with FDR mcp?
tmp2 <- MTsat$mtsat[MTsat$subj.type=="HC" & MTsat$region=="ws"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)

tmp2 <- MTsat$mtsat[MTsat$subj.type=="MS" & MTsat$region=="ws"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
#Conclusion: HC distribution not normal for watershed regions by visual inspection

tmp2 <- MTsat$mtsat[MTsat$subj.type=="HC" & MTsat$region=="mixed"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)

tmp2 <- MTsat$mtsat[MTsat$subj.type=="MS" & MTsat$region=="mixed"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
#Conclusion: MS distribution not normal for mixed regions by visual inspection 

tmp2 <- MTsat$mtsat[MTsat$subj.type=="HC" & MTsat$region=="nws"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)

tmp2 <- MTsat$mtsat[MTsat$subj.type=="MS" & MTsat$region=="nws"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
#Conclusion: Not both normal distributions for non-watershed regions by visual inspection

MTsat %>%
  group_by(subj.type, region) %>%
  shapiro_test(mtsat)

#Conclusion from normality test:
#Visual inspection, combined with Shapiro's test for normality, indicate that the 
#sample distributions for mtsat in ws, nws and mixed regions are 
#significantly different from a normal distribution -> perform Wilcoxon rank test with 
#unpaired samples 

pairwise.wilcox.test(MTsat$mtsat[MTsat$region=="ws"], MTsat$subj.type[MTsat$region=="ws"],
                     p.adjust.method = "fdr",paired=FALSE,exact=FALSE)

MTsat %>% 
  group_by(subj.type, region) %>%
  get_summary_stats(mtsat, type="mean_ci")

# Comparison for WML SUVR  --------------------------------------------------------------
cWMLs.SUVR <- cWMLs.SUVR %>% mutate(subj.type="HC")
WMLs.SUVR <- WMLs.SUVR %>% mutate(subj.type="MS")
SUVR <- rbind(cWMLs.SUVR,WMLs.SUVR)

#Test for normality: paired samples t-test or paired samples wilcoxon test, with FDR mcp?
tmp2 <- SUVR$SUVR[SUVR$subj.type=="HC" & SUVR$region=="ws"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)

tmp2 <- SUVR$SUVR[SUVR$subj.type=="MS" & SUVR$region=="ws"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
#Conclusion: HC distribution not normal for watershed regions by visual inspection

tmp2 <- SUVR$SUVR[SUVR$subj.type=="HC" & SUVR$region=="mixed"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)

tmp2 <- SUVR$SUVR[SUVR$subj.type=="MS" & SUVR$region=="mixed"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
#Conclusion: HC distribution not normal for mixed regions by visual inspection 

tmp2 <- SUVR$SUVR[SUVR$subj.type=="HC" & SUVR$region=="nws"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)

tmp2 <- SUVR$SUVR[SUVR$subj.type=="MS" & SUVR$region=="nws"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
#Conclusion: Not both normal distributions for non-watershed regions by visual inspection

SUVR %>%
  group_by(subj.type, region) %>%
  shapiro_test(SUVR)

#Conclusion from normality test:
#Visual inspection, combined with Shapiro's test for normality, indicate that the 
#sample distributions for SUVR in ws, nws and mixed regions are 
#significantly different from a normal distribution -> perform Wilcoxon rank test with 
#unpaired samples 

pairwise.wilcox.test(SUVR$SUVR[SUVR$region=="nws"], SUVR$subj.type[SUVR$region=="nws"],
                     p.adjust.method = "fdr",paired=FALSE,exact=FALSE)

SUVR %>% 
  group_by(subj.type, region) %>%
  get_summary_stats(SUVR, type="mean_ci")

# EMPTY SPACE -------------------------------------------------------------
# Saving objects ----------------------------------------------------------
#Remove MS patients 'PT_089' and 'ED_087' -> confluent lesion segmentation using J.J Dworkin
#lesiontools algorithm did not work; also omitted from gencWMLmask.sh script 
rmvSubjs <- c("MS_PT_089","MS_ED_087")
#WM lesions
WMLs.mtsat <- WMLs.mtsat[!(WMLs.mtsat$subject.name %in% rmvSubjs),]
WMLs.SUVR <- WMLs.SUVR[!(WMLs.SUVR$subject.name %in% rmvSubjs),]
WMLs.Vol <- WMLs.Vol[!(WMLs.Vol$subject.name %in% rmvSubjs),]
WMLs.Count <- WMLs.Count[!(WMLs.Count$subject.name %in% rmvSubjs),]

#TSPO PET WML phenotypes
RimActive.Vol <- RimActive.Vol[!(RimActive.Vol$subject.name %in% rmvSubjs),]
Active.Vol <- Active.Vol[!(Active.Vol$subject.name %in% rmvSubjs),]
Inactive.Vol <- Inactive.Vol[!(Inactive.Vol$subject.name %in% rmvSubjs),]
RimActive.TotalVol <- RimActive.TotalVol[!(RimActive.TotalVol$subject.name %in% rmvSubjs),]
Active.TotalVol <- Active.TotalVol[!(Active.TotalVol$subject.name %in% rmvSubjs),]
Inactive.TotalVol <- Inactive.TotalVol[!(Inactive.TotalVol$subject.name %in% rmvSubjs),]
RimActive.Count <- RimActive.Count[!(RimActive.Count$subject.name %in% rmvSubjs),]
Active.Count <- Active.Count[!(Active.Count$subject.name %in% rmvSubjs),]
Inactive.Count <- Inactive.Count[!(Inactive.Count$subject.name %in% rmvSubjs),]
RimActive.mtsat <- RimActive.mtsat[!(RimActive.mtsat$subject.name %in% rmvSubjs),]
Active.mtsat <- Active.mtsat[!(Active.mtsat$subject.name %in% rmvSubjs),]
Inactive.mtsat <- Inactive.mtsat[!(Inactive.mtsat$subject.name %in% rmvSubjs),]
RimActive.SUVR <- RimActive.SUVR[!(RimActive.SUVR$subject.name %in% rmvSubjs),]
Active.SUVR <- Active.SUVR[!(Active.SUVR$subject.name %in% rmvSubjs),]
Inactive.SUVR <- Inactive.SUVR[!(Inactive.SUVR$subject.name %in% rmvSubjs),]
RimActive.PerCore <- RimActive.PerCore[!(RimActive.PerCore$subject.name %in% rmvSubjs),]
Active.PerCore <- Active.PerCore[!(Active.PerCore$subject.name %in% rmvSubjs),]

#Setting up Data Frames
cWMLs.mtsat <- cWMLs.mtsat %>% mutate(subj.type="HC")
WMLs.mtsat <- WMLs.mtsat %>% mutate(subj.type="MS")
Fin.WMLs.mtsat <- rbind(cWMLs.mtsat,WMLs.mtsat)

cWMLs.SUVR <- cWMLs.SUVR %>% mutate(subj.type="HC")
WMLs.SUVR <- WMLs.SUVR %>% mutate(subj.type="MS")
Fin.WMLs.SUVR <- rbind(cWMLs.SUVR,WMLs.SUVR)

WMLs.Vol <- WMLs.Vol %>% mutate(subj.type="MS")
WMLs.Count <- WMLs.Count %>% mutate(subj.type="MS")

#TSPO WM lesion phenotypes 
#Individual lesion volume per TSPO PET phenotype
RimActive.Vol <- RimActive.Vol %>% mutate(subj.type="MS")
Active.Vol <- Active.Vol %>% mutate(subj.type="MS")
Inactive.Vol <- Inactive.Vol %>% mutate(subj.type="MS")
#Total lesion volume per TSPO PET phenotype 
RimActive.TotalVol <- RimActive.TotalVol %>% mutate(subj.type="MS")
Active.TotalVol <- Active.TotalVol %>% mutate(subj.type="MS")
Inactive.TotalVol <- Inactive.TotalVol %>% mutate(subj.type="MS")
#Lesion count per TSPO PET phenotype 
RimActive.Count <- RimActive.Count %>% mutate(subj.type="MS")
Active.Count <- Active.Count %>% mutate(subj.type="MS")
Inactive.Count <- Inactive.Count %>% mutate(subj.type="MS")
#Individual mean lesion MTsat per TSPO PET phenotype 
RimActive.mtsat <- RimActive.mtsat %>% mutate(subj.type="MS")
Active.mtsat <- Active.mtsat %>% mutate(subj.type="MS")
Inactive.mtsat <- Inactive.mtsat %>% mutate(subj.type="MS")
#Individual mean lesion SUVR per TSPO PET phenotype 
RimActive.SUVR <- RimActive.SUVR %>% mutate(subj.type="MS")
Active.SUVR <- Active.SUVR %>% mutate(subj.type="MS")
Inactive.SUVR <- Inactive.SUVR %>% mutate(subj.type="MS")
#Individual mean lesion SUVR per TSPO PET phenotype 
RimActive.PerCore <- RimActive.PerCore %>% mutate(subj.type="MS")
Active.PerCore <- Active.PerCore %>% mutate(subj.type="MS")

explaVars <- explaVars[!(unlist(explaVars$subj.name) %in% rmvSubjs),]

save(WMLs.Vol,cWMLs.TotalVol,Fin.WMLs.mtsat,Fin.WMLs.SUVR, WMLs.Count,
     RimActive.Vol, Active.Vol, Inactive.Vol,
     RimActive.TotalVol, Active.TotalVol, Inactive.TotalVol,
     RimActive.Count, Active.Count, Inactive.Count,
     RimActive.mtsat, Active.mtsat, Inactive.mtsat,
     RimActive.SUVR, Active.SUVR, Inactive.SUVR,
     RimActive.PerCore, Active.PerCore,
     explaVars, file=sprintf("%sfinalizedData_withLABs.RData",mainDir))

# DOD Brain Comms Paper #2 Revisions
RimActive.Vol <- RimActive.Vol[RimActive.Vol$bindAff!="LAB",]
Active.Vol <- Active.Vol[Active.Vol$bindAff!="LAB",]
Inactive.Vol <- Inactive.Vol[Inactive.Vol$bindAff!="LAB",]
RimActive.TotalVol <- RimActive.TotalVol[RimActive.TotalVol$bindAff!="LAB",]
Active.TotalVol <- Active.TotalVol[Active.TotalVol$bindAff!="LAB",]
Inactive.TotalVol <- Inactive.TotalVol[Inactive.TotalVol$bindAff!="LAB",]
RimActive.Count <- RimActive.Count[RimActive.Count$bindAff!="LAB",]
Active.Count <- Active.Count[Active.Count$bindAff!="LAB",]
Inactive.Count <- Inactive.Count[Inactive.Count$bindAff!="LAB",]
RimActive.mtsat <- RimActive.mtsat[RimActive.mtsat$bindAff!="LAB",]
Active.mtsat <- Active.mtsat[Active.mtsat$bindAff!="LAB",]
Inactive.mtsat <- Inactive.mtsat[Inactive.mtsat$bindAff!="LAB",]
RimActive.SUVR <- RimActive.SUVR[RimActive.SUVR$bindAff!="LAB",]
Active.SUVR <- Active.SUVR[Active.SUVR$bindAff!="LAB",]
Inactive.SUVR <- Inactive.SUVR[Inactive.SUVR$bindAff!="LAB",]
RimActive.PerCore <- RimActive.PerCore[RimActive.PerCore$bindAff!="LAB",]
Active.PerCore <- Active.PerCore[Active.PerCore$bindAff!="LAB",]
Fin.WMLs.mtsat <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$bindAff!="LAB",]
Fin.WMLs.SUVR <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$bindAff!="LAB",]
explaVars <- explaVars[explaVars$bindAff!="LAB",]

save(WMLs.Vol,Fin.WMLs.mtsat,Fin.WMLs.SUVR, WMLs.Count,
     RimActive.Vol, Active.Vol, Inactive.Vol,
     RimActive.TotalVol, Active.TotalVol, Inactive.TotalVol,
     RimActive.Count, Active.Count, Inactive.Count,
     RimActive.mtsat, Active.mtsat, Inactive.mtsat,
     RimActive.SUVR, Active.SUVR, Inactive.SUVR,
     RimActive.PerCore, Active.PerCore,
     explaVars, file=sprintf("%sfinalizedData_noLABs.RData",mainDir))

# Comparison of region-specific MS data with matched region-specific HC data --------
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggsignif)
library(emmeans)
library(ggprism)
library(ggpubr)
library(rstatix)
library(scales)

mainDir <- "/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/"
outDir <- "/export02/data/risa/DODManuscript2_Preliminary/"
load(sprintf("%sfinalizedData.RData",mainDir))
threshWMLVol <- 6 #threshold [cm3] for WML to avoid boundary singular fit in LMM

#WML mtsat - model testing 
tmp <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$region=="nws",] #switch between "ws", "mixed" & "nws"
#test model 
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + clinical.status + age,
             tmp, REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
#final model 
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + clinical.status,
             tmp, REML=TRUE)
summary(mmod)

#WML SUVR - model testing 
tmp <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$region=="ws",] #switch between "ws", "mixed" & "nws"
#test model 
mmod <- lmer(SUVR ~ (1 | subject.name) + WML.volume + clinical.status + age + bindAff,
             tmp, REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
#final model 
mmod <- lmer(SUVR ~ (1 | subject.name) + WML.volume + clinical.status + bindAff,
             tmp, REML=TRUE)
summary(mmod)

