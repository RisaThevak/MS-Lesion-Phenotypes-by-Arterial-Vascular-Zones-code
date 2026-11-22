# Loading packages --------------------------------------------------------
library(R.matlab)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
library(ggstatsplot)
library(stats)
library(purrr)
library(patchwork)
library(MASS)
library(ggcorrplot)
library(ppcor)
library(scales)
library(broom)
library(rstatix)
library(Hmisc)
library(olsrr)
library(car)
library(reshape)
library(emmeans)
library(FDRestimation)
library(lme4)
library(lmerTest)
library(ggsignif)
library(emmeans)
library(ggprism)
library(ggpubr)

# Loading data ------------------------------------------------------------
statDir <- "/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/"
dat <- readMat(sprintf("%sinfoPET.mat",statDir))
explaVars <-  data.frame(age=t(dat$age),
                         sex=factor(unlist(t(dat$sex)),levels=c('F','M')),
                         subj.type=factor(unlist(t(dat$type)),levels=c('HC','MS')),
                         clinicalStatus=factor(unlist(t(dat$clinicalStatus)),levels=c('HC','RRMS','PMS')),
                         edss=t(dat$edss),
                         pyramidal=t(dat$pyramidal),
                         cerebellar=t(dat$cerebellar),
                         brainstem=t(dat$brainstem),
                         sensory=t(dat$sensory),
                         bowel.bladder=t(dat$bowel.bladder),
                         visual=t(dat$visual),
                         mental=t(dat$mental),
                         diseaseDuration=t(dat$diseaseDuration),
                         age.at.onset=t(dat$age)-t(dat$diseaseDuration),
                         DMT.effect=factor(unlist(t(dat$DMT.effect)),levels=c('HC','None','M-DMT','H-DMT')),
                         bindAff=factor(unlist(t(dat$bindAff)),levels=c('LAB','MAB','HAB')),
                         subj.id=c(1:length(dat$age)),
                         subj.name=t(dat$Names),
                         TIV=t(dat$TIV),
                         DMT=unlist(t(dat$DMT)),
                         WMLTotalVol=unlist(t(dat$WMLvolume)), 
                         WBV=t(dat$WBV))

#C) Global Mean Cortical Thickness
csdat <- readMat(sprintf("%s/Vertexdata_Cohort2/AllVertex_fwhmQMRI%d_cxTh%d_SUVR%d.mat",statDir,5,5,5))
cxTh.lh <- csdat$cxTh.lh
cxTh.rh <- csdat$cxTh.rh
cxTh <- cbind(cxTh.lh,cxTh.rh)
gct <- rowMeans(cxTh,na.rm=TRUE)
# #Reading in whole cortical MTsat data
# mtsat <- readMat("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/mtsat.mat")
# mean.cxMTsat <- rowMeans(t(mtsat$mean.corticalGM),na.rm=TRUE)
# explaVars <- cbind(explaVars,data.frame(gct),data.frame(mean.cxMTsat))
explaVars <- cbind(explaVars,data.frame(gct))

mainDir <- "/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/"
# mainDir <- "/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/NoPVWMLs/"
outDir <- "/export02/data/risa/DODManuscript2_Preliminary/"

#FOR FORMATTING DATA FRAMES
idxsMS <- explaVars$subj.type=="MS" 
subject.name <- unlist(explaVars$subj.name[idxsMS]) 
clinical.status <- as.character(unlist(explaVars$clinicalStatus[idxsMS])); names(clinical.status) <- subject.name
sexList <- as.character(unlist(explaVars$sex[idxsMS])); names(sexList) <- subject.name
ageList <- unlist(explaVars$age[idxsMS]); names(ageList) <- subject.name
bindAffList <- as.character(unlist(explaVars$bindAff[idxsMS])); names(bindAffList) <- subject.name
DMTeffList <- as.character(unlist(explaVars$DMT.effect[idxsMS])); names(DMTeffList) <- subject.name
edssList <- unlist(explaVars$edss[idxsMS]); names(edssList) <- subject.name
gctList <- unlist(explaVars$gct[idxsMS]); names(gctList) <- subject.name
# meanMTsatList <- unlist(explaVars$mean.cxMTsat[idxsMS]); names(meanMTsatList) <- subject.name #CHANGE
disDurList <- unlist(explaVars$diseaseDuration[idxsMS]); names(disDurList) <- subject.name
ageOnsetList <- unlist(explaVars$age.at.onset[idxsMS]); names(ageOnsetList) <- subject.name
WBVList <- unlist(explaVars$WBV[idxsMS]); names(WBVList) <- subject.name
TIVList <- unlist(explaVars$TIV[idxsMS]); names(TIVList) <- subject.name
WMLTotalVolList <- unlist(explaVars$WMLTotalVol[idxsMS]); names(WMLTotalVolList) <- subject.name
pyramidalList <- unlist(explaVars$pyramidal[idxsMS]); names(pyramidalList) <- subject.name
cerebellarList <- unlist(explaVars$cerebellar[idxsMS]); names(cerebellarList) <- subject.name
brainstemList <- unlist(explaVars$brainstem[idxsMS]); names(brainstemList) <- subject.name
sensoryList <- unlist(explaVars$sensory[idxsMS]); names(sensoryList) <- subject.name
bbList <- unlist(explaVars$bowel.bladder[idxsMS]); names(bbList) <- subject.name
visualList <- unlist(explaVars$visual[idxsMS]); names(visualList) <- subject.name
mentalList <- unlist(explaVars$mental[idxsMS]); names(mentalList) <- subject.name
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
  # WMLs.temp <- WMLs.temp %>% mutate(mean.cxMTsat=dplyr::recode(subject.name,!!!meanMTsatList)) #CHANGE
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

ylimMTsat <- c(0,4.2)
ylimPara <- c(0,0.118)*1000
ylimSUVR <- c(0.25,2.3)
ylimWMLVol <- c(-3.5,5)

# 1) Used emmeans() function and FDR for mcp of p-values from t-test of regression coefficient
#(Sattertwaite's approximation for degrees of freedom given random effect terms) 

#2) Find sample means (NOT marginal means used for mcp)

#Test Kenward-Roger's approximation and determine if similar to Sattertwaite approximation 

#4) Is t-statistic square root of F-statistic? 
#Answer: Yes in case of fixed effects terms that are numerical variables 

# WM LESIONS -------------------------------------------------------------- -------------------------------------------------------------
# volume at individual WML level* ----------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$WMLsVol[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.Vol <- stack(tmp)
ws.Vol <- cbind(ws.Vol,rep("ws",nrow(ws.Vol)))
colnames(ws.Vol) <- c("Vol", "subject.name", "region")

mixed <- map(dat$WMLsVol[[2]],unlist);
tmp = list()
tmp[subject.name] <- mixed
mixed.Vol <- stack(tmp)
mixed.Vol <- cbind(mixed.Vol,rep("mixed",nrow(mixed.Vol)))
colnames(mixed.Vol) <- c("Vol", "subject.name", "region")

nws <- map(dat$WMLsVol[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.Vol <- stack(tmp)
nws.Vol <- cbind(nws.Vol,rep("nws",nrow(nws.Vol)))
colnames(nws.Vol) <- c("Vol", "subject.name", "region")

WMLs.Vol <- rbind(ws.Vol,mixed.Vol,nws.Vol)
WMLs.Vol <- format.WMLs.df(WMLs.Vol)
WMLs.Vol$VolCubeRoot <- WMLs.Vol$Vol^(1/3) #Cubic root transformation for normality (reduce skewness and kurtosis)
WMLs.Vol$logVol <- log(WMLs.Vol$Vol)

p1 <- ggplot(WMLs.Vol, aes(region, logVol)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimWMLVol)

mmod <- lmer(logVol ~ (1|subject.name)+region+age, WMLs.Vol,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
anova(mmod)
#Random effect term '(1 | subject.name:region)'
#found to be insignificant. Omitting them. 

# tent.mdl <- lm(logVol ~ region + age, WMLs.Vol)
# ggplot(augment(tent.mdl), aes(x=.fitted,y=.resid)) + geom_point() #residual plot
# plot(tent.mdl) #observations numbered according to full dataset 
# crPlots(tent.mdl) #partial residual plots

#Final model 
mmod <- lmer(logVol ~ region + (1 | subject.name) + age,
              WMLs.Vol,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(WMLs.Vol,logVol~region,p.adjust.method="fdr",model=mmod)

# volume at individual WML level (RRMS)* ----------------------
WMLs.Vol.RRMS <- WMLs.Vol[WMLs.Vol$clinical.status=="RRMS",]

p1b <- ggplot(WMLs.Vol.RRMS, aes(region, logVol)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimWMLVol)

mmodb <- lmer(logVol~(1|subject.name)+region+age,
              WMLs.Vol.RRMS,REML=TRUE)
summary(mmodb)
step(mmodb, direction="backward")
#Fixed effect term 'age' and random effect term '(1 | subject.name:region)'
#found to be insignificant. Omitting them. 

#Final model 
mmodb <- lmer(logVol ~ region + (1 | subject.name),
               WMLs.Vol.RRMS,REML=TRUE)
summary(mmodb)
tmp <- emmeans_test(WMLs.Vol.RRMS,logVol~region,p.adjust.method="fdr",model=mmodb)

# volume at individual WML level (PMS)* ----------------------
WMLs.Vol.PMS <- WMLs.Vol[WMLs.Vol$clinical.status=="PMS",]

p1c <- ggplot(WMLs.Vol.PMS, aes(region, logVol)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimWMLVol)

mmodc <- lmer(logVol~(1|subject.name)+region+age,
               WMLs.Vol.PMS,REML=TRUE)
summary(mmodc)
step(mmodc, direction="backward")
#Fixed effect term 'age' and random effect term '(1 | subject.name:region)'
#found to be insignificant. Omitting them. 

#Final model 
mmodc <- lmer(logVol ~ region + (1 | subject.name),
               WMLs.Vol.PMS,REML=TRUE)
summary(mmodc)
tmp <- emmeans_test(WMLs.Vol.PMS,logVol~region,p.adjust.method="fdr",model=mmodc)

# EMPTY SPACE -------------------------------------------------------------
# WM lesion count -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$WMLscount[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.count <- stack(tmp)
ws.count <- cbind(ws.count,rep("ws",nrow(ws.count)))
colnames(ws.count) <- c("count", "subject.name", "region")

mixed <- map(dat$WMLscount[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.count <- stack(tmp)
mixed.count <- cbind(mixed.count,rep("mixed",nrow(mixed.count)))
colnames(mixed.count) <- c("count", "subject.name", "region")

nws <- map(dat$WMLscount[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.count <- stack(tmp)
nws.count <- cbind(nws.count,rep("nws",nrow(nws.count)))
colnames(nws.count) <- c("count", "subject.name", "region")

WMLs.Count <- rbind(ws.count,mixed.count,nws.count)
WMLs.Count <- format.WMLs.df(WMLs.Count)

pViolinBack <- geom_violin(draw_quantiles = c(0.25, 0.75), trim=TRUE, 
                           linetype="solid",linewidth=0.6)
pViolinFront <- geom_violin(aes(fill=region),draw_quantiles = c(0.5), trim=TRUE, 
                            linewidth=1.2, fill="transparent")
p2 <- ggplot(WMLs.Count, aes(region, count)) + 
  pViolinBack + pViolinFront + geom_jitter(aes(color=subject.name), height=0) 

#Assuming different intercept according to region (arterial vascular zone), 
#same slope w.r.t. WMLcount and age, and subject random-effect 
# mmod <- lmer(count ~ (1 | subject.name) + region + age, WMLs.Count,REML=TRUE)
# summary(mmod)
# step(mmod, direction="backward")
#Fixed effect term 'age' found to be insignificant

#Final model 
mmod <- lmer(count ~ (1 | subject.name) + region,
              WMLs.Count,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(WMLs.Count,count~region,p.adjust.method="fdr",model=mmod)

ggplot(WMLs.Count, aes(x=reorder(subject.name,count),y=count,fill=region)) + 
  geom_bar(stat="identity",position="stack",linewidth=1.2) + 
  facet_wrap(~clinical.status,nrow=1,scales="free_x")

# EMPTY SPACE -------------------------------------------------------------
# mean MTsat at individual WML level -------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$WMLsMTsat[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.mtsat <- stack(tmp)
ws.mtsat <- cbind(ws.mtsat,rep("ws",nrow(ws.mtsat)))
colnames(ws.mtsat) <- c("mtsat", "subject.name", "region")

mixed <- map(dat$WMLsMTsat[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.mtsat <- stack(tmp)
mixed.mtsat <- cbind(mixed.mtsat,rep("mixed",nrow(mixed.mtsat)))
colnames(mixed.mtsat) <- c("mtsat", "subject.name", "region")

nws <- map(dat$WMLsMTsat[[3]],unlist);
tmp = list()
tmp[subject.name] <- nws
nws.mtsat <- stack(tmp)
nws.mtsat <- cbind(nws.mtsat,rep("nws",nrow(nws.mtsat)))
colnames(nws.mtsat) <- c("mtsat", "subject.name", "region")

WMLs.mtsat <- rbind(ws.mtsat,mixed.mtsat,nws.mtsat)
WMLs.mtsat <- format.WMLs.df(WMLs.mtsat)
WMLs.mtsat <- WMLs.mtsat %>% mutate(WML.volume=WMLs.Vol$Vol)
WMLs.mtsat$mtsatpervol <- WMLs.mtsat$mtsat/(WMLs.mtsat$WML.volume*1000)

p3 <- ggplot(WMLs.mtsat, aes(region, mtsat)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimMTsat)

mmod <- lmer(mtsat~(1|subject.name)+region+age+WML.volume,WMLs.mtsat,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
#Fixed effect terms 'age' and 'region' found to be insignificant, but keeping 
#'region' as regressor

#Final model 
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, WMLs.mtsat,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(WMLs.mtsat,mtsat~region,p.adjust.method="fdr",model=mmod)
# ggsave(sprintf("%sModel1_plot.jpg",outDir),p3,width=16,height=16,units="cm",dpi=1200)

# mean MTsat at individual WML level (RRMS) -------------------------
WMLs.mtsat.RRMS <- WMLs.mtsat[WMLs.mtsat$clinical.status=="RRMS",]

p3b <- ggplot(WMLs.mtsat.RRMS, aes(region, mtsat)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimMTsat)

mmodb <- lmer(mtsat~(1|subject.name)+region+age+WML.volume,
               WMLs.mtsat.RRMS,REML=TRUE)
summary(mmodb)
step(mmodb, direction="backward")
#Fixed effect terms 'age', 'region' and 'WML.volume' found to be insignificant, 
#but for consistency keeping 'region' and 'WML.volume'

#Final model 
mmodb <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region,
               WMLs.mtsat.RRMS,REML=TRUE)
summary(mmodb)
tmp <- emmeans_test(WMLs.mtsat.RRMS,mtsat~region,p.adjust.method="fdr",model=mmodb)

# ggsave(sprintf("%sModel1b_plot.jpg",outDir),p3b,width=16,height=16,units="cm",dpi=1200)

# mean MTsat at individual WML level (PMS) -------------------------
WMLs.mtsat.PMS <- WMLs.mtsat[WMLs.mtsat$clinical.status=="PMS",]

p3c <- ggplot(WMLs.mtsat.PMS, aes(region, mtsat)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimMTsat)

mmodc <- lmer(mtsat~(1|subject.name)+region+age+WML.volume,
               WMLs.mtsat.PMS,REML=TRUE)
summary(mmodc)
step(mmodc, direction="backward")
#Fixed effect terms "age" and "region" found to be insignificant, but 
#keeping "region" as a regressor 

#Final model
mmodc <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region,
               WMLs.mtsat.PMS,REML=TRUE)
summary(mmodc)
tmp <- emmeans_test(WMLs.mtsat.PMS,mtsat~region,p.adjust.method="fdr",model=mmodc)
# ggsave(sprintf("%sModel1c_plot.jpg",outDir),p3c,width=16,height=16,units="cm",dpi=1200)

# EMPTY SPACE -------------------------------------------------------------
# mean SUVR at individual WML level* ----------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$WMLsSUVR[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws 
ws.SUVR <- stack(tmp)
ws.SUVR <- cbind(ws.SUVR,rep("ws",nrow(ws.SUVR)))
colnames(ws.SUVR) <- c("SUVR", "subject.name", "region")

mixed <- map(dat$WMLsSUVR[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.SUVR <- stack(tmp)
mixed.SUVR <- cbind(mixed.SUVR,rep("mixed",nrow(mixed.SUVR)))
colnames(mixed.SUVR) <- c("SUVR", "subject.name", "region")

nws <- map(dat$WMLsSUVR[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.SUVR <- stack(tmp)
nws.SUVR <- cbind(nws.SUVR,rep("nws",nrow(nws.SUVR)))
colnames(nws.SUVR) <- c("SUVR", "subject.name", "region")

WMLs.SUVR <- rbind(ws.SUVR,mixed.SUVR,nws.SUVR)
WMLs.SUVR <- format.WMLs.df(WMLs.SUVR)
WMLs.SUVR <- WMLs.SUVR %>% mutate(WML.volume=WMLs.Vol$Vol)

p4 <- ggplot(WMLs.SUVR, aes(region, SUVR)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimSUVR)

mmod <- lmer(SUVR~(1|subject.name)+region+age+WML.volume+bindAff,WMLs.SUVR,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
#Fixed effect term 'age' found to be insignificant 

#Final model 
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff,WMLs.SUVR,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(WMLs.SUVR,SUVR~region,p.adjust.method="fdr",model=mmod)
# ggsave(sprintf("%sModel2_plot.jpg",outDir),p4,width=16,height=16,units="cm",dpi=1200)

# mean SUVR at individual WML level (RRMS)* ----------------------
WMLs.SUVR.RRMS <- WMLs.SUVR[WMLs.SUVR$clinical.status=="RRMS",]

p4b <- ggplot(WMLs.SUVR.RRMS, aes(region, SUVR)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimSUVR)

mmodb <- lmer(SUVR~(1|subject.name)+region+age+WML.volume+bindAff,
               WMLs.SUVR.RRMS,REML=TRUE)
summary(mmodb)
step(mmodb, direction="backward")

#Fixed effect terms 'age' and 'WML.volume' and random effect 
#term (1 | subject.name:region) found to be insignificant. But omitting 
#only 'age' as regressor for consistency 

#Final model 
mmodb <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff,
               WMLs.SUVR.RRMS,REML=TRUE)
summary(mmodb)
tmp <- emmeans_test(WMLs.SUVR.RRMS,SUVR~region,p.adjust.method="fdr",model=mmodb)
# ggsave(sprintf("%sModel2b_plot.jpg",outDir),p4b,width=16,height=16,units="cm",dpi=1200)

# mean SUVR at individual WML level (PMS)* ----------------------
WMLs.SUVR.PMS <- WMLs.SUVR[WMLs.SUVR$clinical.status=="PMS",]

p4c <- ggplot(WMLs.SUVR.PMS, aes(region, SUVR)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimSUVR)

mmodc <- lmer(SUVR~(1|subject.name)+region+age+WML.volume+bindAff,
               WMLs.SUVR.PMS,REML=TRUE)
summary(mmodc)
step(mmodc, direction="backward")
#Fixed effect term 'age' and random effect term '(1 | subject.name:region)
#found to be insignificant. For consistency, only omitting 'age' as regressor 

#Final model 
mmodc <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume+bindAff,
               WMLs.SUVR.PMS,REML=TRUE)
summary(mmodc)
tmp <- emmeans_test(WMLs.SUVR.PMS,SUVR~region,p.adjust.method="fdr",model=mmodc)
# ggsave(sprintf("%sModel2c_plot.jpg",outDir),p4c,width=16,height=16,units="cm",dpi=1200)

# TSPO WM lesion phenotyping ----------------------------------------------- -------------------------------------------------------------
# Rim active lesion count -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$RimActive[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.count <- stack(tmp)
ws.count <- cbind(ws.count,rep("ws",nrow(ws.count)))
colnames(ws.count) <- c("count", "subject.name", "region")

mixed <- map(dat$RimActive[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.count <- stack(tmp)
mixed.count <- cbind(mixed.count,rep("mixed",nrow(mixed.count)))
colnames(mixed.count) <- c("count", "subject.name", "region")

nws <- map(dat$RimActive[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.count <- stack(tmp)
nws.count <- cbind(nws.count,rep("nws",nrow(nws.count)))
colnames(nws.count) <- c("count", "subject.name", "region")

RimActive.Count <- rbind(ws.count,mixed.count,nws.count)
RimActive.Count <- format.WMLs.df(RimActive.Count)

# Active lesion count -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$Active[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.count <- stack(tmp)
ws.count <- cbind(ws.count,rep("ws",nrow(ws.count)))
colnames(ws.count) <- c("count", "subject.name", "region")

mixed <- map(dat$Active[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.count <- stack(tmp)
mixed.count <- cbind(mixed.count,rep("mixed",nrow(mixed.count)))
colnames(mixed.count) <- c("count", "subject.name", "region")

nws <- map(dat$Active[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.count <- stack(tmp)
nws.count <- cbind(nws.count,rep("nws",nrow(nws.count)))
colnames(nws.count) <- c("count", "subject.name", "region")

Active.Count <- rbind(ws.count,mixed.count,nws.count)
Active.Count <- format.WMLs.df(Active.Count)

# Inactive lesion count -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$Inactive[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.count <- stack(tmp)
ws.count <- cbind(ws.count,rep("ws",nrow(ws.count)))
colnames(ws.count) <- c("count", "subject.name", "region")

mixed <- map(dat$Inactive[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.count <- stack(tmp)
mixed.count <- cbind(mixed.count,rep("mixed",nrow(mixed.count)))
colnames(mixed.count) <- c("count", "subject.name", "region")

nws <- map(dat$Inactive[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.count <- stack(tmp)
nws.count <- cbind(nws.count,rep("nws",nrow(nws.count)))
colnames(nws.count) <- c("count", "subject.name", "region")

Inactive.Count <- rbind(ws.count,mixed.count,nws.count)
Inactive.Count <- format.WMLs.df(Inactive.Count)

# Patient-wise Proportions of TSPO WML Phenotypes (Count) -----------------------------------------------------
TSPOWMLs.Count <- RimActive.Count$count + Active.Count$count + Inactive.Count$count
RimActive.Count$Prop <- RimActive.Count$count*100/TSPOWMLs.Count
Active.Count$Prop <- Active.Count$count*100/TSPOWMLs.Count
Inactive.Count$Prop <- Inactive.Count$count*100/TSPOWMLs.Count

# Remove patients where TSPO PET WML phenotyping criteria returned no WMLs at all 
RimActive.Count <- RimActive.Count[TSPOWMLs.Count!=0,]
Active.Count <- Active.Count[TSPOWMLs.Count!=0,]
Inactive.Count <- Inactive.Count[TSPOWMLs.Count!=0,]

# Plotting Rim-Active Lesion Count 
pViolinBack <- geom_violin(draw_quantiles = c(0.25, 0.75), trim=TRUE, 
                           linetype="solid",linewidth=0.6)
pViolinFront <- geom_violin(aes(fill=region),draw_quantiles = c(0.5), trim=TRUE, 
                            linewidth=1.2, fill="transparent")
ggplot(RimActive.Count, aes(region, count)) + 
  pViolinBack + pViolinFront + geom_jitter(aes(color=subject.name), height=0) 

#Assuming different intercept according to region (arterial vascular zone), 
#same slope w.r.t. WMLcount and age, and subject random-effect 
mmod <- lmer(count ~ (1 | subject.name) + region + age, RimActive.Count,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")

#Final model 
mmod <- lmer(count ~ (1 | subject.name) + region,
             RimActive.Count,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(RimActive.Count,count~region,p.adjust.method="fdr",model=mmod)

# Plotting Active Lesion Count 
ggplot(Active.Count, aes(region, count)) + 
  pViolinBack + pViolinFront + geom_jitter(aes(color=subject.name), height=0) + ylim(0,40) 

#Assuming different intercept according to region (arterial vascular zone), 
#same slope w.r.t. WMLcount and age, and subject random-effect 
mmod <- lmer(count ~ (1 | subject.name) + region + age, Active.Count,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")

#Final model 
mmod <- lmer(count ~ (1 | subject.name) + region,
             Active.Count,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(Active.Count,count~region,p.adjust.method="fdr",model=mmod)

# Plotting Inactive Lesion Count 
ggplot(Inactive.Count, aes(region, count)) + 
  pViolinBack + pViolinFront + geom_jitter(aes(color=subject.name), height=0) 

#Assuming different intercept according to region (arterial vascular zone), 
#same slope w.r.t. WMLcount and age, and subject random-effect 
mmod <- lmer(count ~ (1 | subject.name) + region + age, Inactive.Count,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")

#Final model 
mmod <- lmer(count ~ (1 | subject.name) + region,
             Inactive.Count,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(Inactive.Count,count~region,p.adjust.method="fdr",model=mmod)

# EMPTY SPACE -------------------------------------------------------------
# Rim active lesion Vol -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$RimActiveVol[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.Vol <- stack(tmp)
ws.Vol <- cbind(ws.Vol,rep("ws",nrow(ws.Vol)))
colnames(ws.Vol) <- c("Vol", "subject.name", "region")

mixed <- map(dat$RimActiveVol[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.Vol <- stack(tmp)
mixed.Vol <- cbind(mixed.Vol,rep("mixed",nrow(mixed.Vol)))
colnames(mixed.Vol) <- c("Vol", "subject.name", "region")

nws <- map(dat$RimActiveVol[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.Vol <- stack(tmp)
nws.Vol <- cbind(nws.Vol,rep("nws",nrow(nws.Vol)))
colnames(nws.Vol) <- c("Vol", "subject.name", "region")

RimActive.Vol <- rbind(ws.Vol,mixed.Vol,nws.Vol)
RimActive.Vol <- format.WMLs.df(RimActive.Vol)
RimActive.Vol$VolCubeRoot <- RimActive.Vol$Vol^(1/3) #Cubic root transformation for normality (reduce skewness and kurtosis)
RimActive.Vol$logVol <- log(RimActive.Vol$Vol)

p1 <- ggplot(RimActive.Vol, aes(region, logVol)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimWMLVol)

mmod <- lmer(logVol~(1|subject.name)+region+age,
             RimActive.Vol,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
anova(mmod)
#Random effect term '(1 | subject.name:region)'
#found to be insignificant. Omitting them. 

#Final model 
mmod <- lmer(logVol ~ region + (1 | subject.name) + age,
             RimActive.Vol,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(RimActive.Vol,logVol~region,p.adjust.method="fdr",model=mmod)

# Active lesion Vol -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$ActiveVol[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.Vol <- stack(tmp)
ws.Vol <- cbind(ws.Vol,rep("ws",nrow(ws.Vol)))
colnames(ws.Vol) <- c("Vol", "subject.name", "region")

mixed <- map(dat$ActiveVol[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.Vol <- stack(tmp)
mixed.Vol <- cbind(mixed.Vol,rep("mixed",nrow(mixed.Vol)))
colnames(mixed.Vol) <- c("Vol", "subject.name", "region")

nws <- map(dat$ActiveVol[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.Vol <- stack(tmp)
nws.Vol <- cbind(nws.Vol,rep("nws",nrow(nws.Vol)))
colnames(nws.Vol) <- c("Vol", "subject.name", "region")

Active.Vol <- rbind(ws.Vol,mixed.Vol,nws.Vol)
Active.Vol <- format.WMLs.df(Active.Vol)
Active.Vol$VolCubeRoot <- Active.Vol$Vol^(1/3) #Cubic root transformation for normality (reduce skewness and kurtosis)
Active.Vol$logVol <- log(Active.Vol$Vol)

p2 <- ggplot(Active.Vol, aes(region, logVol)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimWMLVol)

mmod <- lmer(logVol~(1|subject.name)+region+age,
             Active.Vol,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
anova(mmod)
#Random effect term '(1 | subject.name:region)'
#found to be insignificant. Omitting them. 

#Final model 
mmod <- lmer(logVol ~ region + (1 | subject.name) + age,
             Active.Vol,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(Active.Vol,logVol~region,p.adjust.method="fdr",model=mmod)

# Inactive lesion Vol -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$InactiveVol[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.Vol <- stack(tmp)
ws.Vol <- cbind(ws.Vol,rep("ws",nrow(ws.Vol)))
colnames(ws.Vol) <- c("Vol", "subject.name", "region")

mixed <- map(dat$InactiveVol[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.Vol <- stack(tmp)
mixed.Vol <- cbind(mixed.Vol,rep("mixed",nrow(mixed.Vol)))
colnames(mixed.Vol) <- c("Vol", "subject.name", "region")

nws <- map(dat$InactiveVol[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.Vol <- stack(tmp)
nws.Vol <- cbind(nws.Vol,rep("nws",nrow(nws.Vol)))
colnames(nws.Vol) <- c("Vol", "subject.name", "region")

Inactive.Vol <- rbind(ws.Vol,mixed.Vol,nws.Vol)
Inactive.Vol <- format.WMLs.df(Inactive.Vol)
Inactive.Vol$VolCubeRoot <- Inactive.Vol$Vol^(1/3) #Cubic root transformation for normality (reduce skewness and kurtosis)
Inactive.Vol$logVol <- log(Inactive.Vol$Vol)

p3 <- ggplot(Inactive.Vol, aes(region, logVol)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimWMLVol)

mmod <- lmer(logVol~(1|subject.name)+region+age,
             Inactive.Vol,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
anova(mmod)
#Random effect term '(1 | subject.name:region)'
#found to be insignificant. Omitting them. 

#Final model 
mmod <- lmer(logVol ~ region + (1 | subject.name) + age,
             Inactive.Vol,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(Inactive.Vol,logVol~region,p.adjust.method="fdr",model=mmod)

# Patient-wise Proportions of TSPO WML Phenotypes (Volume) -----------------------------------------------------
RimActive.TotalVol <- RimActive.Vol %>%
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(sumVol = sum(Vol), age = max(age), edss=max(edss), bindAff=first(bindAff), disDur=max(disDur), gct=max(gct))

Active.TotalVol <- Active.Vol %>%
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(sumVol = sum(Vol), age = max(age), edss=max(edss), bindAff=first(bindAff), disDur=max(disDur), gct=max(gct))

Inactive.TotalVol <- Inactive.Vol %>%
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(sumVol = sum(Vol), age = max(age), edss=max(edss), bindAff=first(bindAff), disDur=max(disDur), gct=max(gct))

tmp <- rbind(RimActive.TotalVol, Active.TotalVol, Inactive.TotalVol)
TSPOWMLs.TotalVol <- tmp %>%
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(sumVol = sum(sumVol))

for (i in levels(TSPOWMLs.TotalVol$subject.name))
{
  for (j in levels(TSPOWMLs.TotalVol$region)) {
    idxs2 <- TSPOWMLs.TotalVol$subject.name==i & TSPOWMLs.TotalVol$region==j 
    
    idxs1 <- RimActive.TotalVol$subject.name==i & RimActive.TotalVol$region==j
    if (any(idxs1) && any(idxs2)) {
      RimActive.TotalVol$Prop[idxs1] <- RimActive.TotalVol$sumVol[idxs1]*100/TSPOWMLs.TotalVol$sumVol[idxs2]
    }
    idxs1 <- Active.TotalVol$subject.name==i & Active.TotalVol$region==j
    if (any(idxs1) && any(idxs2)) {
      Active.TotalVol$Prop[idxs1] <- Active.TotalVol$sumVol[idxs1]*100/TSPOWMLs.TotalVol$sumVol[idxs2]
    }
    idxs1 <- Inactive.TotalVol$subject.name==i & Inactive.TotalVol$region==j
    if (any(idxs1) && any(idxs2)) {
      Inactive.TotalVol$Prop[idxs1] <- Inactive.TotalVol$sumVol[idxs1]*100/TSPOWMLs.TotalVol$sumVol[idxs2]
    }
  }
} 

mmod <- lmer(Prop ~ (1 | subject.name) + region + age, RimActive.TotalVol,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")

# EMPTY SPACE -------------------------------------------------------------
# Rim active lesion MTsat -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$RimActiveMTsat[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.mtsat <- stack(tmp)
ws.mtsat <- cbind(ws.mtsat,rep("ws",nrow(ws.mtsat)))
colnames(ws.mtsat) <- c("mtsat", "subject.name", "region")

mixed <- map(dat$RimActiveMTsat[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.mtsat <- stack(tmp)
mixed.mtsat <- cbind(mixed.mtsat,rep("mixed",nrow(mixed.mtsat)))
colnames(mixed.mtsat) <- c("mtsat", "subject.name", "region")

nws <- map(dat$RimActiveMTsat[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.mtsat <- stack(tmp)
nws.mtsat <- cbind(nws.mtsat,rep("nws",nrow(nws.mtsat)))
colnames(nws.mtsat) <- c("mtsat", "subject.name", "region")

RimActive.mtsat <- rbind(ws.mtsat,mixed.mtsat,nws.mtsat)
RimActive.mtsat <- format.WMLs.df(RimActive.mtsat)
RimActive.mtsat <- RimActive.mtsat %>% mutate(WML.volume=RimActive.Vol$Vol)
RimActive.mtsat$mtsatpervol <- RimActive.mtsat$mtsat/(RimActive.mtsat$WML.volume*1000)

p1 <- ggplot(RimActive.mtsat, aes(region, mtsat)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimMTsat)

mmod <- lmer(mtsat~(1|subject.name)+region+age+WML.volume,RimActive.mtsat,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
#Fixed effect terms 'age' and 'region' found to be insignificant, but keeping 
#'region' as regressor

#Final model 
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, RimActive.mtsat,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(RimActive.mtsat,mtsat~region,p.adjust.method="fdr",model=mmod)
# ggsave(sprintf("%sModel1_plot.jpg",outDir),p3,width=16,height=16,units="cm",dpi=1200)

# Active lesion MTsat -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$ActiveMTsat[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.mtsat <- stack(tmp)
ws.mtsat <- cbind(ws.mtsat,rep("ws",nrow(ws.mtsat)))
colnames(ws.mtsat) <- c("mtsat", "subject.name", "region")

mixed <- map(dat$ActiveMTsat[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.mtsat <- stack(tmp)
mixed.mtsat <- cbind(mixed.mtsat,rep("mixed",nrow(mixed.mtsat)))
colnames(mixed.mtsat) <- c("mtsat", "subject.name", "region")

nws <- map(dat$ActiveMTsat[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.mtsat <- stack(tmp)
nws.mtsat <- cbind(nws.mtsat,rep("nws",nrow(nws.mtsat)))
colnames(nws.mtsat) <- c("mtsat", "subject.name", "region")

Active.mtsat <- rbind(ws.mtsat,mixed.mtsat,nws.mtsat)
Active.mtsat <- format.WMLs.df(Active.mtsat)
Active.mtsat <- Active.mtsat %>% mutate(WML.volume=Active.Vol$Vol)
Active.mtsat$mtsatpervol <- Active.mtsat$mtsat/(Active.mtsat$WML.volume*1000)

p2 <- ggplot(Active.mtsat, aes(region, mtsat)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimMTsat)

mmod <- lmer(mtsat~(1|subject.name)+region+age+WML.volume,Active.mtsat,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
#Fixed effect terms 'age' and 'region' found to be insignificant, but keeping 
#'region' as regressor

#Final model 
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, Active.mtsat,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(Active.mtsat,mtsat~region,p.adjust.method="fdr",model=mmod)
# ggsave(sprintf("%sModel1_plot.jpg",outDir),p3,width=16,height=16,units="cm",dpi=1200)

# Inactive lesion MTsat -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$InactiveMTsat[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.mtsat <- stack(tmp)
ws.mtsat <- cbind(ws.mtsat,rep("ws",nrow(ws.mtsat)))
colnames(ws.mtsat) <- c("mtsat", "subject.name", "region")

mixed <- map(dat$InactiveMTsat[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.mtsat <- stack(tmp)
mixed.mtsat <- cbind(mixed.mtsat,rep("mixed",nrow(mixed.mtsat)))
colnames(mixed.mtsat) <- c("mtsat", "subject.name", "region")

nws <- map(dat$InactiveMTsat[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.mtsat <- stack(tmp)
nws.mtsat <- cbind(nws.mtsat,rep("nws",nrow(nws.mtsat)))
colnames(nws.mtsat) <- c("mtsat", "subject.name", "region")

Inactive.mtsat <- rbind(ws.mtsat,mixed.mtsat,nws.mtsat)
Inactive.mtsat <- format.WMLs.df(Inactive.mtsat)
Inactive.mtsat <- Inactive.mtsat %>% mutate(WML.volume=Inactive.Vol$Vol)
Inactive.mtsat$mtsatpervol <- Inactive.mtsat$mtsat/(Inactive.mtsat$WML.volume*1000)

p3 <- ggplot(Inactive.mtsat, aes(region, mtsat)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimMTsat)

mmod <- lmer(mtsat~(1|subject.name)+region+age+WML.volume,Inactive.mtsat,REML=TRUE)
summary(mmod)
step(mmod, direction="backward")
#Fixed effect terms 'age' and 'region' found to be insignificant, but keeping 
#'region' as regressor

#Final model 
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, Inactive.mtsat,REML=TRUE)
summary(mmod)
tmp <- emmeans_test(Inactive.mtsat,mtsat~region,p.adjust.method="fdr",model=mmod)
# ggsave(sprintf("%sModel1_plot.jpg",outDir),p3,width=16,height=16,units="cm",dpi=1200)

# Rim active lesion SUVR -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$RimActiveSUVR[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.SUVR <- stack(tmp)
ws.SUVR <- cbind(ws.SUVR,rep("ws",nrow(ws.SUVR)))
colnames(ws.SUVR) <- c("SUVR", "subject.name", "region")

mixed <- map(dat$RimActiveSUVR[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.SUVR <- stack(tmp)
mixed.SUVR <- cbind(mixed.SUVR,rep("mixed",nrow(mixed.SUVR)))
colnames(mixed.SUVR) <- c("SUVR", "subject.name", "region")

nws <- map(dat$RimActiveSUVR[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.SUVR <- stack(tmp)
nws.SUVR <- cbind(nws.SUVR,rep("nws",nrow(nws.SUVR)))
colnames(nws.SUVR) <- c("SUVR", "subject.name", "region")

RimActive.SUVR <- rbind(ws.SUVR,mixed.SUVR,nws.SUVR)
RimActive.SUVR <- format.WMLs.df(RimActive.SUVR)
RimActive.SUVR <- RimActive.SUVR %>% mutate(WML.volume=RimActive.Vol$Vol)

p1 <- ggplot(RimActive.SUVR, aes(region, SUVR)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimSUVR)

# Active lesion SUVR -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$ActiveSUVR[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.SUVR <- stack(tmp)
ws.SUVR <- cbind(ws.SUVR,rep("ws",nrow(ws.SUVR)))
colnames(ws.SUVR) <- c("SUVR", "subject.name", "region")

mixed <- map(dat$ActiveSUVR[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.SUVR <- stack(tmp)
mixed.SUVR <- cbind(mixed.SUVR,rep("mixed",nrow(mixed.SUVR)))
colnames(mixed.SUVR) <- c("SUVR", "subject.name", "region")

nws <- map(dat$ActiveSUVR[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.SUVR <- stack(tmp)
nws.SUVR <- cbind(nws.SUVR,rep("nws",nrow(nws.SUVR)))
colnames(nws.SUVR) <- c("SUVR", "subject.name", "region")

Active.SUVR <- rbind(ws.SUVR,mixed.SUVR,nws.SUVR)
Active.SUVR <- format.WMLs.df(Active.SUVR)
Active.SUVR <- Active.SUVR %>% mutate(WML.volume=Active.Vol$Vol)

p2 <- ggplot(Active.SUVR, aes(region, SUVR)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimSUVR)

# Inactive lesion SUVR -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$InactiveSUVR[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.SUVR <- stack(tmp)
ws.SUVR <- cbind(ws.SUVR,rep("ws",nrow(ws.SUVR)))
colnames(ws.SUVR) <- c("SUVR", "subject.name", "region")

mixed <- map(dat$InactiveSUVR[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.SUVR <- stack(tmp)
mixed.SUVR <- cbind(mixed.SUVR,rep("mixed",nrow(mixed.SUVR)))
colnames(mixed.SUVR) <- c("SUVR", "subject.name", "region")

nws <- map(dat$InactiveSUVR[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.SUVR <- stack(tmp)
nws.SUVR <- cbind(nws.SUVR,rep("nws",nrow(nws.SUVR)))
colnames(nws.SUVR) <- c("SUVR", "subject.name", "region")

Inactive.SUVR <- rbind(ws.SUVR,mixed.SUVR,nws.SUVR)
Inactive.SUVR <- format.WMLs.df(Inactive.SUVR)
Inactive.SUVR <- Inactive.SUVR %>% mutate(WML.volume=Inactive.Vol$Vol)

p3 <- ggplot(Inactive.SUVR, aes(region, SUVR)) + 
  geom_violin() + geom_boxplot() + geom_jitter(aes(color=subject.name)) + ylim(ylimSUVR)

# Rim active lesion PerCore -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$RimActivePerCore[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.PerCore <- stack(tmp)
ws.PerCore <- cbind(ws.PerCore,rep("ws",nrow(ws.PerCore)))
colnames(ws.PerCore) <- c("PerCore", "subject.name", "region")

mixed <- map(dat$RimActivePerCore[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.PerCore <- stack(tmp)
mixed.PerCore <- cbind(mixed.PerCore,rep("mixed",nrow(mixed.PerCore)))
colnames(mixed.PerCore) <- c("PerCore", "subject.name", "region")

nws <- map(dat$RimActivePerCore[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.PerCore <- stack(tmp)
nws.PerCore <- cbind(nws.PerCore,rep("nws",nrow(nws.PerCore)))
colnames(nws.PerCore) <- c("PerCore", "subject.name", "region")

RimActive.PerCore <- rbind(ws.PerCore,mixed.PerCore,nws.PerCore)
RimActive.PerCore <- format.WMLs.df(RimActive.PerCore)
RimActive.PerCore <- RimActive.PerCore %>% mutate(WML.volume=RimActive.Vol$Vol)

# Active lesion PerCore -----------------------------------------------------
dat <- readMat(sprintf("%sWMLsMeasures.mat",mainDir))

ws <- map(dat$ActivePerCore[[1]],unlist); 
tmp = list()
tmp[subject.name] <- ws
ws.PerCore <- stack(tmp)
ws.PerCore <- cbind(ws.PerCore,rep("ws",nrow(ws.PerCore)))
colnames(ws.PerCore) <- c("PerCore", "subject.name", "region")

mixed <- map(dat$ActivePerCore[[2]],unlist); 
tmp = list()
tmp[subject.name] <- mixed
mixed.PerCore <- stack(tmp)
mixed.PerCore <- cbind(mixed.PerCore,rep("mixed",nrow(mixed.PerCore)))
colnames(mixed.PerCore) <- c("PerCore", "subject.name", "region")

nws <- map(dat$ActivePerCore[[3]],unlist); 
tmp = list()
tmp[subject.name] <- nws
nws.PerCore <- stack(tmp)
nws.PerCore <- cbind(nws.PerCore,rep("nws",nrow(nws.PerCore)))
colnames(nws.PerCore) <- c("PerCore", "subject.name", "region")

Active.PerCore <- rbind(ws.PerCore,mixed.PerCore,nws.PerCore)
Active.PerCore <- format.WMLs.df(Active.PerCore)
Active.PerCore <- Active.PerCore %>% mutate(WML.volume=Active.Vol$Vol)







