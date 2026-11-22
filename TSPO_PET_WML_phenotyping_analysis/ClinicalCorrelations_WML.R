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
library(patchwork)

mainDir <- "/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/"
outDir <- "/export02/data/risa/DODManuscript2_Revisions/"
load(sprintf("%sfinalizedData_noLABs.RData",mainDir))

# WM lesion count  ----------------------
tmp <- WMLs.Count 
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "count", y = "gct", #CHANGE
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "count", ylab = "gct") #CHANGE
result <- cor.test(tmp$count, tmp$gct, use = "complete.obs",method = "spearman",exact=FALSE)
#Count, mixed 
#Count, nws

# WM lesion average volume ----------------------
tmp <- WMLs.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct), omsss=max(omsss))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "mean.logVol", y = "gct", #CHANGE
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "mean log volume over all lesions", ylab = "gct") #CHANGE
#Volume cube root, mixed

# WM lesion total volume ----------------------
tmp <- WMLs.Vol %>% 
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(total.Vol = sum(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct), omsss=max(omsss))
selIdxs <- tmp$region=="mixed" #CHANGE
tmp <- tmp[selIdxs,]
ggscatter(tmp, x = "total.Vol", y = "gct", #CHANGE
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Total lesion volume", ylab = "gct") #CHANGE

#Total lesion volume, mixed

# WM lesion average MTsat -------------------------
WMLs.mtsat <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$subj.type=="MS",] #CUSTOM
tmp <- WMLs.mtsat %>% 
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.mtsat = mean(mtsat), age = max(age), edss=max(edss), 
                   gct=max(gct), omsss=max(omsss), disDur=max(disDur))
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "mean.mtsat", y = "gct", #CHANGE
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "mean mtsat over all lesions", ylab = "gct") #CHANGE

#MTsat, nws

# WM lesion average SUVR ----------------------
WMLs.SUVR <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$subj.type=="MS",] #CUSTOM
tmp <- WMLs.SUVR %>% 
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.SUVR = mean(SUVR), age = max(age), edss=max(edss), 
                   gct=max(gct), disDur=max(disDur))
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "mean.SUVR", y = "gct", #CHANGE
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "mean SUVR over all lesions", ylab = "gct") #CHANGE
#SUVR, ws
#SUVR, nws

# WM lesion average MTsat per unit lesion volume -------------------------
WMLs.mtsat <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$subj.type=="MS",] #CUSTOM
tmp <- WMLs.mtsat %>% 
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.mtsatpervol = mean(mtsatpervol), age = max(age), edss=max(edss), 
                   gct=max(gct), omsss=max(omsss), disDur=max(disDur))
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "mean.mtsatpervol", y = "gct", #CHANGE
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "mean mtsat over all lesions", ylab = "gct") #CHANGE

#MTsat, nws

# TSPO PET WML phenotypes --------------------------------------- -------------------------------------------------

# TSPO PET WML phenotype count ----------------------
tmp <- Active.Count #CHANGE
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "count", y = "gct", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "count", ylab = "gct")
result <- cor.test(tmp$count, tmp$gct, use = "complete.obs",method = "spearman",exact=FALSE)

#Inactive.Count, nws (R = 0.42)
#Inactive.Count, mixed (R = 0.55)
#InInactive.Count, ws (R = -0.43)

# TSPO PET WML phenotype proportion by count ----------------------
tmp <- RimActive.Count #CHANGE
selIdxs <- tmp$region=="ws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "Prop", y = "gct", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Proportion of Total WM Count", ylab = "gct")
result <- cor.test(tmp$Prop, tmp$gct, use = "complete.obs",method = "spearman",exact=FALSE)

#InInactive, Proportion by count, ws (R = -0.33)
#InInactive, Proportion by count, mixed (R = -0.33)

# TSPO PET WML phenotype average volume ----------------------
tmp <- RimActive.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
# tmp <- RimInactive.Vol #CHANGE
selIdxs <- tmp$region=="ws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "mean.logVol", y = "edss", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "mean log volume over all lesions", ylab = "edss")
result <- cor.test(tmp$mean.logVol, tmp$edss, use = "complete.obs",method = "spearman",exact=FALSE)

#RimInactive.Vol, ws
#Inactive.Vol, mixed

# TSPO PET WML phenotype proportion by total volume ----------------------
tmp <- RimActive.TotalVol #CHANGE
selIdxs <- tmp$region=="ws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "Prop", y = "edss", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Proportion of Total WM Volume", ylab = "edss")
result <- cor.test(tmp$Prop, tmp$edss, use = "complete.obs",method = "spearman",exact=FALSE)

#RimInactive, Proportion by Total Volume, ws (R = 0.54)
#RimInactive, Proportion by Total Volume, nws (R = 0.52)

# TSPO PET WML phenotype average MTsat ----------------------
tmp <- Active.mtsat %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.mtsat = mean(mtsat), age = max(age), edss=max(edss), 
                   gct=max(gct), disDur=max(disDur))
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
0

# TSPO PET WML phenotype average MTsat per unit lesion volume ----------------------
tmp <- Active.mtsat %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.mtsatpervol = mean(mtsatpervol), age = max(age), edss=max(edss), 
                   gct=max(gct), disDur=max(disDur))
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
ggscatter(tmp, x = "mean.mtsatpervol", y = "gct", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "mean mtsat over all lesions", ylab = "gct")
result <- cor.test(tmp$mean.mtsatpervol, tmp$gct, use = "complete.obs",method = "spearman",exact=FALSE)

# # TSPO PET WML phenotype total volumes ----------------------
# tmp <- Active.TotalVol #CHANGE
# selIdxs <- tmp$region=="mixed" #CHANGE
# tmp <- tmp[selIdxs,]
# rownames(tmp) <- c(1:nrow(tmp))
# ggscatter(tmp, x = "sumVol", y = "edss", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "spearman",
#           xlab = "Total WM Volume", ylab = "EDSS")
# result <- cor.test(tmp$sumVol, tmp$edss, use = "complete.obs",method = "spearman",exact=FALSE)
# 
# #RimActive.TotalVol, ws (R = 0.59)
# #Active.TotalVol, mixed (R = 0.62)

# # mean para at individual WML level ----------------------
# WMLs.para <- Fin.WMLs.para[Fin.WMLs.para$subj.type=="MS",] #CUSTOM
# # WMLs.para <- Fin.WMLs.paraBHs[Fin.WMLs.paraBHs$subj.type=="MS",] #CUSTOM
# # WMLs.para <- Fin.WMLs.paranonBHs[Fin.WMLs.paranonBHs$subj.type=="MS",] #CUSTOM
# selIdxs <- WMLs.para$region=="mixed" #CUSTOM
# tmp <- WMLs.para[selIdxs,]
# rownames(tmp) <- c(1:nrow(tmp))
# mdl <- lm(para ~ age + WML.volume, tmp)
# summary(mdl)
# #are original row numbers from data frame "tmp"
# tmpIdxs <- as.integer(names(resid(mdl))) 
# tmp$residX[tmpIdxs] <- resid(mdl)
# tmp <- tmp %>% #CHANGE
#   dplyr::group_by(subject.name, region) %>%
#   dplyr::summarise(residX = mean(residX), age = max(age), edss=max(edss))
# ggscatter(tmp, x = "residX", y = "edss", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "spearman",
#           xlab = "para", ylab = "EDSS")
# result <- cor.test(tmp$para, tmp$edss, use = "complete.obs",method = "spearman",exact=FALSE)

# Multiple Comparisons Correction -----------------------------------------
pvals <- c(0.6, 0.07, 0.2, 0.01)

res <- p.adjust(pvals, method="fdr")

# # Figure 10a, b ----------------------------------------------------------
# #ggplot2 common components for customized plots
# # colors.TSPOPhenotype <- setNames(c('purple4','violet','pink'),c("RA", "A", "IA"))
# # colors.zone <- setNames(c('indianred2','orange','royalblue3'),levels(WMLs.Vol$region))
# # pScaleClrs <- scale_fill_manual(values=colors.zone)
# pPoint <- geom_point(fill='black',size=1.8,shape=21)
# totalWidth <- 26.67
# ph <- 8
# ptheme1 <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
#                  axis.text.x = element_text(size=9),
#                  axis.text.y = element_text(size=12),
#                  axis.title.y = element_text(size=12),
#                  axis.title.x = element_text(size=12),
#                  text = element_text(family="Times New Roman"),
#                  legend.position="none")
# ptheme2 <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
#                  axis.text.x = element_text(size=9),
#                  axis.text.y = element_blank(),
#                  axis.title.y = element_blank(),
#                  axis.title.x = element_text(size=12),
#                  text = element_text(family="Times New Roman"),
#                  legend.position="none")
# 
# #Average WML SUVR in ws zone panel #1
# WMLs.SUVR <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$subj.type=="MS",] #CUSTOM
# tmp <- WMLs.SUVR %>% 
#   dplyr::group_by(subject.name, region) %>%
#   dplyr::summarise(mean.SUVR = mean(SUVR), age = max(age), edss=max(edss), 
#                    gct=max(gct), omsss=max(omsss), disDur=max(disDur))
# selIdxs <- tmp$region=="ws" #CHANGE
# tmp <- tmp[selIdxs,]
# rownames(tmp) <- c(1:nrow(tmp))
# yScale <- scale_y_continuous(limits=c(0.5,1.25),breaks=c(0.5,0.75,1,1.25))
# p1a <- ggplot(tmp, aes(age,mean.SUVR))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme1 +
#   labs(title="", y = ~paste(bold("Average lesion SUVR")), x = "Age (years)") + 
#   yScale
# p1b <- ggplot(tmp, aes(edss,mean.SUVR))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "EDSS") + yScale
# p1c <- ggplot(tmp, aes(disDur,mean.SUVR))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "Disease Duration (years)") + yScale
# p1d <- ggplot(tmp, aes(omsss,mean.SUVR))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "MSSS") + yScale
# p1e <- ggplot(tmp, aes(gct,mean.SUVR))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# # plot1 <- p1a + p1b + p1c + p1d + p1e + plot_layout(nrow=1,ncol=5)
# 
# #WM lesion count in mixed zone panel #2
# tmp <- WMLs.Count 
# selIdxs <- tmp$region=="mixed" #CHANGE
# tmp <- tmp[selIdxs,]
# rownames(tmp) <- c(1:nrow(tmp))
# yScale <- scale_y_continuous(limits=c(0,65),breaks=c(0,10,20,30,40,50,60))
# p2a <- ggplot(tmp, aes(age,count))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme1 +
#   labs(title="", y = ~paste(bold("Lesion count")), x = "Age (years)") + 
#   yScale
# p2b <- ggplot(tmp, aes(edss,count))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "EDSS") + yScale
# p2c <- ggplot(tmp, aes(disDur,count))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Disease Duration (years)") + yScale
# p2d <- ggplot(tmp, aes(omsss,count))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "MSSS") + yScale
# p2e <- ggplot(tmp, aes(gct,count))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# # plot2 <- p2a + p2b + p2c + p2d + p2e + plot_layout(nrow=1,ncol=5)
# 
# #Average lesion volume in mixed zone panel #1 
# tmp <- WMLs.Vol %>% #CHANGE
#   dplyr::group_by(subject.name, region) %>%
#   dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
#                    disDur=max(disDur), gct=max(gct), omsss=max(omsss))
# tmp$mean.logVol <- log(tmp$mean.Vol)
# selIdxs <- tmp$region=="mixed" #CHANGE
# tmp <- tmp[selIdxs,]
# rownames(tmp) <- c(1:nrow(tmp))
# yScale <- scale_y_continuous(limits=c(-4,2),breaks=c(-4,-3,-2,-1,0,1,2))
# p3a <- ggplot(tmp, aes(age,mean.logVol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme1 +
#   labs(title="", y = ~paste(bold("Average lesion volume ("),bold("cm"^"3"),
#                             bold(") "),bold("(log-trans)")), x = "Age (years)") + 
#   yScale
# p3b <- ggplot(tmp, aes(edss,mean.logVol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "EDSS") + yScale
# p3c <- ggplot(tmp, aes(disDur,mean.logVol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Disease Duration (years)") + yScale
# p3d <- ggplot(tmp, aes(omsss,mean.logVol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "MSSS") + yScale
# p3e <- ggplot(tmp, aes(gct,mean.logVol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# # plot3 <- p3a + p3b + p3c + p3d + p3e + plot_layout(nrow=1,ncol=5)
# 
# #Total WM lesion volume in mixed zone panel #2 
# tmp <- WMLs.Vol %>% 
#   dplyr::group_by(subject.name, region) %>%
#   dplyr::summarise(total.Vol = sum(Vol), age = max(age), edss=max(edss), 
#                    disDur=max(disDur), gct=max(gct), omsss=max(omsss))
# tmp$total.Vol <- log(tmp$total.Vol)
# selIdxs <- tmp$region=="mixed" #CHANGE
# tmp <- tmp[selIdxs,]
# yScale <- scale_y_continuous(limits=c(-4,4.5),breaks=c(-4,-3,-2,-1,0,1,2,3,4))
# p4a <- ggplot(tmp, aes(age,total.Vol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme1 +
#   labs(title="", y = ~paste(bold("Total lesion volume ("),bold("cm"^"3"),
#                             bold(") "),bold("(log-trans)")), x = "Age (years)") + 
#   yScale
# p4b <- ggplot(tmp, aes(edss,total.Vol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "EDSS") + yScale
# p4c <- ggplot(tmp, aes(disDur,total.Vol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Disease Duration (years)") + yScale
# p4d <- ggplot(tmp, aes(omsss,total.Vol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "MSSS") + yScale
# p4e <- ggplot(tmp, aes(gct,total.Vol))  +
#   pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# # plot4 <- p4a + p4b + p4c + p4d + p4e + plot_layout(nrow=1,ncol=5)
# 
# plot1 <- p1a + p1b + p1c + p1d + p1e + 
#   p2a + p2b + p2c + p2d + p2e +  plot_layout(nrow=2,ncol=5)
# 
# ggsave(sprintf("%sFigureG.jpg",outDir),plot1,width=totalWidth,height=ph*2+1.6,units="cm",dpi=1200)
# 
# plot2 <- p3a + p3b + p3c + p3d + p3e +
#   p4a + p4b + p4c + p4d + p4e + plot_layout(nrow=2,ncol=5)
# 
# ggsave(sprintf("%sFigureH.jpg",outDir),plot2,width=totalWidth,height=ph*2+1.6,units="cm",dpi=1200)
# 
# # Figure 11a, b ----------------------------------------------------------
# #ggplot2 common components for customized plots
# # colors.TSPOPhenotype <- setNames(c('purple4','violet','pink'),c("RA", "A", "IA"))
# # colors.zone <- setNames(c('indianred2','orange','royalblue3'),levels(WMLs.Vol$region))
# # pScaleClrs <- scale_fill_manual(values=colors.zone)
# pPoint12 <- geom_point(fill='purple4',size=1.8,shape=21)
# pPoint34 <- geom_point(fill='violet',size=1.8,shape=21)
# totalWidth <- 26.67
# ph <- 8
# ptheme1 <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
#                  axis.text.x = element_text(size=9),
#                  axis.text.y = element_text(size=12),
#                  axis.title.y = element_text(size=12),
#                  axis.title.x = element_text(size=12),
#                  text = element_text(family="Times New Roman"),
#                  legend.position="none")
# ptheme2 <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
#                  axis.text.x = element_text(size=9),
#                  axis.text.y = element_blank(),
#                  axis.title.y = element_blank(),
#                  axis.title.x = element_text(size=12),
#                  text = element_text(family="Times New Roman"),
#                  legend.position="none")
# 
# #Rim-Active lesions in ws zone panel #1
# tmp <- RimActive.Vol %>% #CHANGE
#   dplyr::group_by(subject.name, region) %>%
#   dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
#                    disDur=max(disDur), gct=max(gct), omsss=max(omsss))
# tmp$mean.logVol <- log(tmp$mean.Vol)
# selIdxs <- tmp$region=="ws" #CHANGE
# tmp <- tmp[selIdxs,]
# rownames(tmp) <- c(1:nrow(tmp))
# yScale <- scale_y_continuous(limits=c(-4,2),breaks=c(-4,-3,-2,-1,0,1,2))
# p1a <- ggplot(tmp, aes(age,mean.logVol))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme1 +
#   labs(title="", y = ~paste(bold("Average lesion volume ("),bold("cm"^"3"),
#                             bold(") "),bold("(log-trans)")), x = "Age (years)") + 
#   yScale
# p1b <- ggplot(tmp, aes(edss,mean.logVol))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "EDSS") + yScale
# p1c <- ggplot(tmp, aes(disDur,mean.logVol))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "Disease Duration (years)") + yScale
# p1d <- ggplot(tmp, aes(omsss,mean.logVol))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "MSSS") + yScale
# p1e <- ggplot(tmp, aes(gct,mean.logVol))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# # plot1 <- p1a + p1b + p1c + p1d + p1e + plot_layout(nrow=1,ncol=5)
# 
# #Rim-Active lesions in ws zone panel #2
# tmp <- RimActive.TotalVol #CHANGE
# selIdxs <- tmp$region=="ws" #CHANGE
# tmp <- tmp[selIdxs,]
# rownames(tmp) <- c(1:nrow(tmp))
# yScale <- scale_y_continuous(limits=c(0,90),breaks=c(0,10,20,30,40,50,60,70,80,90))
# p2a <- ggplot(tmp, aes(age,Prop))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme1 +
#   labs(title="", y = ~paste(bold("Proportion of total lesion volume (%)")), x = "Age (years)") + 
#   yScale
# p2b <- ggplot(tmp, aes(edss,Prop))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "EDSS") + yScale
# p2c <- ggplot(tmp, aes(disDur,Prop))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "Disease Duration (years)") + yScale
# p2d <- ggplot(tmp, aes(omsss,Prop))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "MSSS") + yScale
# p2e <- ggplot(tmp, aes(gct,Prop))  +
#   pPoint12 + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
#   labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# # plot2 <- p2a + p2b + p2c + p2d + p2e + plot_layout(nrow=1,ncol=5)
# 
# #Active lesions in mixed zone panel #1 
# tmp <- Active.Count #CHANGE
# selIdxs <- tmp$region=="mixed" #CHANGE
# tmp <- tmp[selIdxs,]
# rownames(tmp) <- c(1:nrow(tmp))
# yScale <- scale_y_continuous(limits=c(0,65),breaks=c(0,10,20,30,40,50,60))
# p3a <- ggplot(tmp, aes(age,count))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme1 +
#   labs(title="", y = ~paste(bold("Lesion count")), x = "Age (years)") + 
#   yScale
# p3b <- ggplot(tmp, aes(edss,count))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "EDSS") + yScale
# p3c <- ggplot(tmp, aes(disDur,count))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Disease Duration (years)") + yScale
# p3d <- ggplot(tmp, aes(omsss,count))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "MSSS") + yScale
# p3e <- ggplot(tmp, aes(gct,count))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# # plot3 <- p3a + p3b + p3c + p3d + p3e + plot_layout(nrow=1,ncol=5)
# 
# #Active lesions in mixed zone panel #2 
# tmp <- Active.Vol %>% #CHANGE
#   dplyr::group_by(subject.name, region) %>%
#   dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
#                    disDur=max(disDur), gct=max(gct), omsss=max(omsss))
# tmp$mean.logVol <- log(tmp$mean.Vol)
# selIdxs <- tmp$region=="mixed" #CHANGE
# tmp <- tmp[selIdxs,]
# rownames(tmp) <- c(1:nrow(tmp))
# yScale <- scale_y_continuous(limits=c(-4,2),breaks=c(-4,-3,-2,-1,0,1,2))
# p4a <- ggplot(tmp, aes(age,mean.logVol))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme1 +
#   labs(title="", y = ~paste(bold("Average lesion volume ("),bold("cm"^"3"),
#                             bold(") "),bold("(log-trans)")), x = "Age (years)") + 
#   yScale
# p4b <- ggplot(tmp, aes(edss,mean.logVol))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "EDSS") + yScale
# p4c <- ggplot(tmp, aes(disDur,mean.logVol))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Disease Duration (years)") + yScale
# p4d <- ggplot(tmp, aes(omsss,mean.logVol))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "MSSS") + yScale
# p4e <- ggplot(tmp, aes(gct,mean.logVol))  +
#   pPoint34 + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
#   labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# # plot4 <- p4a + p4b + p4c + p4d + p4e + plot_layout(nrow=1,ncol=5)
# 
# plot1 <- p1a + p1b + p1c + p1d + p1e + 
#   p2a + p2b + p2c + p2d + p2e +  plot_layout(nrow=2,ncol=5)
# ggsave(sprintf("%sFigureI.jpg",outDir),plot1,width=totalWidth,height=ph*2+1.6,units="cm",dpi=1200)
# 
# plot2 <- p3a + p3b + p3c + p3d + p3e +
#   p4a + p4b + p4c + p4d + p4e + plot_layout(nrow=2,ncol=5)
# ggsave(sprintf("%sFigureJ.jpg",outDir),plot2,width=totalWidth,height=ph*2+1.6,units="cm",dpi=1200)

# Figure 10a, b ----------------------------------------------------------
#ggplot2 common components for customized plots
# colors.TSPOPhenotype <- setNames(c('purple4','violet','pink'),c("RA", "A", "IA"))
# colors.zone <- setNames(c('indianred2','orange','royalblue3'),levels(WMLs.Vol$region))
# pScaleClrs <- scale_fill_manual(values=colors.zone)
pPoint <- geom_point(fill='black',size=1.8,shape=21)
totalWidth <- 26.67
ph <- 8
ptheme1 <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
                 axis.text.x = element_text(size=9),
                 axis.text.y = element_text(size=12),
                 axis.title.y = element_text(size=12),
                 axis.title.x = element_text(size=12),
                 text = element_text(family="Times New Roman"),
                 legend.position="none")
ptheme2 <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
                 axis.text.x = element_text(size=9),
                 axis.text.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.title.x = element_text(size=12),
                 text = element_text(family="Times New Roman"),
                 legend.position="none")

#Average lesion volume in ws zone panel #1 
tmp <- WMLs.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="ws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(-4,1.5),breaks=c(-4,-3,-2,-1,0,1))
p1a <- ggplot(tmp, aes(age,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme1 +
  labs(title="", y = "", x = "Age (years)") + 
  yScale
p1b <- ggplot(tmp, aes(edss,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  labs(title="", x = "EDSS") + yScale
p1c <- ggplot(tmp, aes(disDur,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  labs(title="", x = "Disease Duration (years)") + yScale
p1d <- ggplot(tmp, aes(gct,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot1 <- p1a + p1b + p1c + p1d + plot_layout(nrow=1,ncol=4)

#Average lesion volume in mixed zone panel #2
tmp <- WMLs.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="mixed" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(-4,1.5),breaks=c(-4,-3,-2,-1,0,1))
p2a <- ggplot(tmp, aes(age,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme1 +
  labs(title="", y = "", x = "Age (years)") + 
  yScale
p2b <- ggplot(tmp, aes(edss,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  labs(title="", x = "EDSS") + yScale
p2c <- ggplot(tmp, aes(disDur,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  labs(title="", x = "Disease Duration (years)") + yScale
p2d <- ggplot(tmp, aes(gct,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot2 <- p2a + p2b + p2c + p2d + plot_layout(nrow=1,ncol=4)

#Average lesion volume in nws zone panel #3
tmp <- WMLs.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(-4,1.5),breaks=c(-4,-3,-2,-1,0,1))
p3a <- ggplot(tmp, aes(age,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme1 +
  labs(title="", y = "", x = "Age (years)") + 
  yScale
p3b <- ggplot(tmp, aes(edss,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  labs(title="", x = "EDSS") + yScale
p3c <- ggplot(tmp, aes(disDur,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  labs(title="", x = "Disease Duration (years)") + yScale
p3d <- ggplot(tmp, aes(gct,mean.logVol))  +
  pPoint + geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot3 <- p3a + p3b + p3c + p3d + plot_layout(nrow=1,ncol=4)

#Average WML SUVR in ws zone panel #1
WMLs.SUVR <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$subj.type=="MS",] #CUSTOM
tmp <- WMLs.SUVR %>% 
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.SUVR = mean(SUVR), age = max(age), edss=max(edss), 
                   gct=max(gct), disDur=max(disDur))
selIdxs <- tmp$region=="ws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(0.5,1.25),breaks=c(0.5,0.75,1,1.25))
p4a <- ggplot(tmp, aes(age,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme1 +
  labs(title="", y = "", x = "Age (years)") + 
  yScale
p4b <- ggplot(tmp, aes(edss,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
  labs(title="", x = "EDSS") + yScale
p4c <- ggplot(tmp, aes(disDur,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
  labs(title="", x = "Disease Duration (years)") + yScale
p4d <- ggplot(tmp, aes(gct,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='indianred2') +  ptheme2 +
  labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot4 <- p4a + p4b + p4c + p4d + plot_layout(nrow=1,ncol=4)

#Average WML SUVR in mixed zone panel #2
WMLs.SUVR <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$subj.type=="MS",] #CUSTOM
tmp <- WMLs.SUVR %>% 
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.SUVR = mean(SUVR), age = max(age), edss=max(edss), 
                   gct=max(gct), disDur=max(disDur))
selIdxs <- tmp$region=="mixed" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(0.5,1.25),breaks=c(0.5,0.75,1,1.25))
p5a <- ggplot(tmp, aes(age,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme1 +
  labs(title="", y="", x = "Age (years)") + 
  yScale
p5b <- ggplot(tmp, aes(edss,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
  labs(title="", x = "EDSS") + yScale
p5c <- ggplot(tmp, aes(disDur,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
  labs(title="", x = "Disease Duration (years)") + yScale
p5d <- ggplot(tmp, aes(gct,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='orange') +  ptheme2 +
  labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot5 <- p5a + p5b + p5c + p5d + plot_layout(nrow=1,ncol=4)

#Average WML SUVR in nws zone panel #3
WMLs.SUVR <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$subj.type=="MS",] #CUSTOM
tmp <- WMLs.SUVR %>% 
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.SUVR = mean(SUVR), age = max(age), edss=max(edss), 
                   gct=max(gct), disDur=max(disDur))
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(0.5,1.25),breaks=c(0.5,0.75,1,1.25))
p6a <- ggplot(tmp, aes(age,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='royalblue3') +  ptheme1 +
  labs(title="", y="", x = "Age (years)") + 
  yScale
p6b <- ggplot(tmp, aes(edss,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='royalblue3') +  ptheme2 +
  labs(title="", x = "EDSS") + yScale
p6c <- ggplot(tmp, aes(disDur,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='royalblue3') +  ptheme2 +
  labs(title="", x = "Disease Duration (years)") + yScale
p6d <- ggplot(tmp, aes(gct,mean.SUVR))  +
  pPoint + geom_smooth(method="lm", color="black",fill='royalblue3') +  ptheme2 +
  labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot6 <- p6a + p6b + p6c + p6d + plot_layout(nrow=1,ncol=4)

plot1 <- p1a + p1b + p1c + p1d + 
  p2a + p2b + p2c + p2d + 
  p3a + p3b + p3c + p3d + plot_layout(nrow=3,ncol=4) 
ggsave(sprintf("%sSupp_Figure_1_Revised.jpg",outDir),plot1,width=totalWidth-1,height=ph*3+1.6,units="cm",dpi=1200)

# plot2 <- p4a + p4b + p4c + p4d + 
#   p5a + p5b + p5c + p5d + 
#   p6a + p6b + p6c + p6d + plot_layout(nrow=3,ncol=4)
# ggsave(sprintf("%sFigureH.jpg",outDir),plot2,width=totalWidth-1,height=ph*3+1.6,units="cm",dpi=1200)

# Figure 11a, b ----------------------------------------------------------
#ggplot2 common components for customized plots
# colors.TSPOPhenotype <- setNames(c('purple4','violet','pink'),c("RA", "A", "IA"))
# colors.zone <- setNames(c('indianred2','orange','royalblue3'),levels(WMLs.Vol$region))
# pScaleClrs <- scale_fill_manual(values=colors.zone)
pPoint1 <- geom_point(fill='purple4',size=1.8,shape=21)
pPoint2 <- geom_point(fill='violet',size=1.8,shape=21)
totalWidth <- 26.67
ph <- 8
ptheme1 <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
                 axis.text.x = element_text(size=9),
                 axis.text.y = element_text(size=12),
                 axis.title.y = element_text(size=12),
                 axis.title.x = element_text(size=12),
                 text = element_text(family="Times New Roman"),
                 legend.position="none")
ptheme2 <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
                 axis.text.x = element_text(size=9),
                 axis.text.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.title.x = element_text(size=12),
                 text = element_text(family="Times New Roman"),
                 legend.position="none")

#Rim-Active lesions in ws zone panel #1
tmp <- RimActive.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="ws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(-4,0),breaks=c(-4,-3,-2,-1,0))
p1a <- ggplot(tmp, aes(age,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme1 +
  pPoint1 + labs(title="", y = "", x = "Age (years)") + 
  yScale
p1b <- ggplot(tmp, aes(edss,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint1 + labs(title="", x = "EDSS") + yScale
p1c <- ggplot(tmp, aes(disDur,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint1 + labs(title="", x = "Disease Duration (years)") + yScale
p1d <- ggplot(tmp, aes(gct,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint1 + labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot1 <- p1a + p1b + p1c + p1d + plot_layout(nrow=1,ncol=4)

#Rim-Active lesions in mixed zone panel #2
tmp <- RimActive.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="mixed" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(-4,0),breaks=c(-4,-3,-2,-1,0))
p2a <- ggplot(tmp, aes(age,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme1 +
  pPoint1 + labs(title="", y = "", x = "Age (years)") + 
  yScale
p2b <- ggplot(tmp, aes(edss,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint1 + labs(title="", x = "EDSS") + yScale
p2c <- ggplot(tmp, aes(disDur,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint1 + labs(title="", x = "Disease Duration (years)") + yScale
p2d <- ggplot(tmp, aes(gct,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint1 + labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot2 <- p2a + p2b + p2c + p2d + plot_layout(nrow=1,ncol=4)

#Rim-Active lesions in nws zone panel #3
tmp <- RimActive.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(-4,0),breaks=c(-4,-3,-2,-1,0))
p3a <- ggplot(tmp, aes(age,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme1 +
  pPoint1 + labs(title="", y = "", x = "Age (years)") + 
  yScale
p3b <- ggplot(tmp, aes(edss,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint1 + labs(title="", x = "EDSS") + yScale
p3c <- ggplot(tmp, aes(disDur,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint1 + labs(title="", x = "Disease Duration (years)") + yScale
p3d <- ggplot(tmp, aes(gct,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint1 + labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot3 <- p3a + p3b + p3c + p3d + plot_layout(nrow=1,ncol=4)

#Active lesions in ws zone panel #1
tmp <- Active.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="ws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(-4,1.5),breaks=c(-4,-3,-2,-1,0,1))
p4a <- ggplot(tmp, aes(age,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme1 +
  pPoint2 + labs(title="", y = "", x = "Age (years)") + 
  yScale
p4b <- ggplot(tmp, aes(edss,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint2 + labs(title="", x = "EDSS") + yScale
p4c <- ggplot(tmp, aes(disDur,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint2 + labs(title="", x = "Disease Duration (years)") + yScale
p4d <- ggplot(tmp, aes(gct,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint2 + labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot4 <- p4a + p4b + p4c + p4d + plot_layout(nrow=1,ncol=4)

#Active lesions in mixed zone panel #2
tmp <- Active.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="mixed" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(-4,1.5),breaks=c(-4,-3,-2,-1,0,1))
p5a <- ggplot(tmp, aes(age,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme1 +
  pPoint2 + labs(title="", y = "", x = "Age (years)") + 
  yScale
p5b <- ggplot(tmp, aes(edss,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint2 + labs(title="", x = "EDSS") + yScale
p5c <- ggplot(tmp, aes(disDur,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint2 + labs(title="", x = "Disease Duration (years)") + yScale
p5d <- ggplot(tmp, aes(gct,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint2 + labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot5 <- p5a + p5b + p5c + p5d + plot_layout(nrow=1,ncol=4)

#Active lesions in nws zone panel #3
tmp <- Active.Vol %>% #CHANGE
  dplyr::group_by(subject.name, region) %>%
  dplyr::summarise(mean.Vol = mean(Vol), age = max(age), edss=max(edss), 
                   disDur=max(disDur), gct=max(gct))
tmp$mean.logVol <- log(tmp$mean.Vol)
selIdxs <- tmp$region=="nws" #CHANGE
tmp <- tmp[selIdxs,]
rownames(tmp) <- c(1:nrow(tmp))
yScale <- scale_y_continuous(limits=c(-4,1.5),breaks=c(-4,-3,-2,-1,0,1))
p6a <- ggplot(tmp, aes(age,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme1 +
  pPoint2 + labs(title="", y = "", x = "Age (years)") + 
  yScale
p6b <- ggplot(tmp, aes(edss,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint2 + labs(title="", x = "EDSS") + yScale
p6c <- ggplot(tmp, aes(disDur,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint2 + labs(title="", x = "Disease Duration (years)") + yScale
p6d <- ggplot(tmp, aes(gct,mean.logVol))  +
  geom_smooth(method="lm", color="black",se=FALSE) +  ptheme2 +
  pPoint2 + labs(title="", x = "Global Cortical Thickness (mm)") + yScale
# plot6 <- p6a + p6b + p6c + p6d + plot_layout(nrow=1,ncol=4)

# ylab = ~paste(bold("Average lesion volume ("),bold("cm"^"3"),
#            bold(") "),bold("(log-trans)"))
plot1 <- p1a + p1b + p1c + p1d + 
  p2a + p2b + p2c + p2d + 
  p3a + p3b + p3c + p3d + plot_layout(nrow=3,ncol=4) 
ggsave(sprintf("%sFigure7_Revised.jpg",outDir),plot1,width=totalWidth-1,height=ph*3+1.6,units="cm",dpi=1200)

plot2 <- p4a + p4b + p4c + p4d + 
  p5a + p5b + p5c + p5d + 
  p6a + p6b + p6c + p6d + plot_layout(nrow=3,ncol=4)
ggsave(sprintf("%sSupp_Figure_2_Revised.jpg",outDir),plot2,width=totalWidth-1,height=ph*3+1.6,units="cm",dpi=1200)

# Official Plot #2: Clinical severity bar/violin plots ---------------------
clrs <- setNames(c('white','black','gray'),levels(explaVars$clinicalStatus))
pScaleFill <- scale_fill_manual(values=clrs)
pScaleClr <- scale_color_manual(values=clrs)
pPointB <- geom_point(aes(fill=clinical.status),size=1.2,shape=21,stroke=0.8,
                      position=position_jitterdodge(jitter.height=0,jitter.width=1.3))
pViolinBack <- geom_boxplot(aes(color=clinical.status,alpha=0.6),
                            linetype="solid",linewidth=0.6, outlier.shape = NA, coef = 0)
pViolinFront <- geom_violin(aes(color=clinical.status),trim=TRUE,
                            linewidth=1.0, fill="transparent")
pthemeB <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold", vjust=0),
                 axis.text.x = element_text(size=8),
                 axis.text.y = element_text(size=12),
                 axis.title.y = element_text(size=12),
                 axis.title.x = element_blank(),
                 text = element_text(family="Times New Roman"),
                 legend.position="none")
add_stats <- function(df,...) {
  add_pvalue(data=df, 
             label = "{p.signif}",
             fontfamily="Times New Roman",
             label.size=5,
             tip.length = 0.05, ...)
}

#1) Rim-active ws lesion volume (RRMS/PMS)
tmp <- RimActive.Vol[RimActive.Vol$region=="ws",] %>% #CHANGE
  dplyr::group_by(subject.name) %>%
  dplyr::summarise(mean.Vol = mean(Vol), clinical.status=first(clinical.status), age=first(age))
tmp$mean.logVol <- log(tmp$mean.Vol)
df.form <- tmp %>%
  rstatix::t_test(mean.logVol ~ clinical.status) %>%
  rstatix::add_xy_position(x = "clinical.status") # important for positioning!
df.form$y.position <- df.form$y.position + 0.10

#MUST DO: Perform normality test
tmp2 <- tmp$mean.logVol[tmp$clinical.status=="RRMS"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
tmp2 <- tmp$mean.logVol[tmp$clinical.status=="PMS"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
tmp %>%
  group_by(clinical.status) %>%
  shapiro_test(mean.logVol)
#Conclusion from normality test:
#Visual inspection, combined with Shapiro's test for normality, indicate that the 
#sample distributions for log average rim-active ws lesions in RRMS and PMS patients 
#DO NOT significantly diverge from normal distributions 

#Linear regression adjusting for age, then get P-value
mmod <- lm(mean.logVol ~ age + clinical.status, tmp) #CUSTOM
summary(mmod)
df.form$p.signif <- "ns" 

p1 <- ggplot(tmp, aes(clinical.status,mean.logVol)) + pViolinBack + pViolinFront + 
  pPointB + pScaleClr + pScaleFill + pthemeB + 
  labs(y="Average rim-active ws lesion volume",title="A") + add_stats(df.form)

#2) Rim-active ws lesion count (RRMS/PMS)
tmp <- RimActive.Count[RimActive.Count$region=="ws",]
df.form <- tmp %>%
  rstatix::t_test(count ~ clinical.status) %>%
  rstatix::add_xy_position(x = "clinical.status") # important for positioning!
df.form$y.position <- df.form$y.position + 0.10

#MUST DO: Perform normality test
tmp2 <- tmp$count[tmp$clinical.status=="RRMS"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
tmp2 <- tmp$count[tmp$clinical.status=="PMS"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
tmp %>%
  group_by(clinical.status) %>%
  shapiro_test(count)
#Conclusion from normality test:
#Visual inspection, combined with Shapiro's test for normality, indicate that the 
#sample distributions for rim-active ws lesion count in RRMS and PMS patients 
#significantly diverge from normal distributions 

#Mann-Whitney U Test
pairwise.wilcox.test(tmp$count,tmp$clinical.status,
                     p.adjust.method = "fdr",paired=FALSE,exact=FALSE)
df.form$p.signif <- "*" #P=0.05 

p2 <- ggplot(tmp, aes(clinical.status,count)) + pViolinBack + pViolinFront + 
  pPointB + pScaleClr + pScaleFill + pthemeB + 
  labs(y="Rim-active ws lesion count",title="B") + add_stats(df.form)

#5) Rim-active ws lesion volume (cGML High/cGML Low)
med.cxMTsat.HC <- median(explaVars$mean.cxMTsat[explaVars$subj.type=="HC"],na.rm=TRUE)
IQR.cxMTsat.HC <- IQR(explaVars$mean.cxMTsat[explaVars$subj.type=="HC"],na.rm=TRUE)
thresh <- med.cxMTsat.HC - 2*IQR.cxMTsat.HC

tmp <- RimActive.Vol
tmp$cxMTsat.status[tmp$mean.cxMTsat<thresh] <- "cGML High"
tmp$cxMTsat.status[tmp$mean.cxMTsat>=thresh] <- "cGML Low"
tmp$cxMTsat.status <- as.factor(tmp$cxMTsat.status)
tmp <- tmp[tmp$region=="ws",] %>% #CHANGE
  dplyr::group_by(subject.name) %>%
  dplyr::summarise(mean.Vol = mean(Vol), cxMTsat.status=first(cxMTsat.status), age=first(age))
tmp$mean.logVol <- log(tmp$mean.Vol)
df.form <- tmp %>%
  rstatix::t_test(mean.logVol ~ cxMTsat.status) %>%
  rstatix::add_xy_position(x = "cxMTsat.status") # important for positioning!
df.form$y.position <- df.form$y.position + 0.10

#MUST DO: Perform normality test
tmp2 <- tmp$mean.logVol[tmp$cxMTsat.status=="cGML High"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
tmp2 <- tmp$mean.logVol[tmp$cxMTsat.status=="cGML Low"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
tmp %>%
  group_by(cxMTsat.status) %>%
  shapiro_test(mean.logVol)
#Conclusion from normality test:
#Visual inspection, combined with Shapiro's test for normality, indicate that the 
#sample distributions for rim-active ws lesion count in RRMS and PMS patients 
#DO NOT significantly diverge from normal distributions

#Linear regression adjusting for age and get P-value
mmod <- lm(mean.logVol ~ age + cxMTsat.status, tmp) #CUSTOM
summary(mmod)
test.result <- emmeans_test(tmp,logVol~cxMTsat.status,p.adjust.method="fdr",model=mmod) #CUSTOM
df.form$p.signif <- "**" #P = 0.007

clrs <- setNames(c('firebrick','orangered'),levels(tmp$cxMTsat.status))
pScaleFill <- scale_fill_manual(values=clrs)
pScaleClr <- scale_color_manual(values=clrs)
pPointB <- geom_point(aes(fill=cxMTsat.status),size=1.2,shape=21,stroke=0.8,
                      position=position_jitterdodge(jitter.height=0,jitter.width=1.3))
pViolinBack <- geom_boxplot(aes(color=cxMTsat.status,alpha=0.6),
                            linetype="solid",linewidth=0.6, outlier.shape = NA, coef = 0)
pViolinFront <- geom_violin(aes(color=cxMTsat.status),trim=TRUE,
                            linewidth=1.0, fill="transparent")
p3 <- ggplot(tmp, aes(cxMTsat.status,mean.logVol)) + pViolinBack + pViolinFront +
  pPointB + pScaleFill + pScaleClr + pthemeB + labs(y="Average rim-active ws lesion volume",title="C") + add_stats(df.form)

#6) Rim-active ws lesion count (cGML High/cGML Low)
tmp <- RimActive.Count
tmp$cxMTsat.status[tmp$mean.cxMTsat<thresh] <- "cGML High"
tmp$cxMTsat.status[tmp$mean.cxMTsat>=thresh] <- "cGML Low"
tmp$cxMTsat.status <- as.factor(tmp$cxMTsat.status)
tmp <- tmp[tmp$region=="ws",] 
df.form <- tmp %>%
  rstatix::t_test(count ~ cxMTsat.status) %>%
  rstatix::add_xy_position(x = "cxMTsat.status") # important for positioning!
df.form$y.position <- df.form$y.position + 0.10

#MUST DO: Perform normality test
tmp2 <- tmp$count[tmp$cxMTsat.status=="cGML High"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
tmp2 <- tmp$count[tmp$cxMTsat.status=="cGML Low"] #CHANGE
ggdensity(tmp2, fill="lightpink")
ggqqplot(tmp2)
tmp %>%
  group_by(cxMTsat.status) %>%
  shapiro_test(count)
#Conclusion from normality test:
#Visual inspection, combined with Shapiro's test for normality, indicate that the 
#sample distributions for rim-active ws lesion count in RRMS and PMS patients 
#significantly diverge from normal distributions 

#Mann-Whitney U Test
pairwise.wilcox.test(tmp$count,tmp$cxMTsat.status,
                     p.adjust.method = "fdr",paired=FALSE,exact=FALSE)
df.form$p.signif <- "***" #P=0.0009

clrs <- setNames(c('firebrick','orangered'),levels(tmp$cxMTsat.status))
pScaleFill <- scale_fill_manual(values=clrs)
pScaleClr <- scale_color_manual(values=clrs)
pPointB <- geom_point(aes(fill=cxMTsat.status),size=1.2,shape=21,stroke=0.8,
                      position=position_jitterdodge(jitter.height=0,jitter.width=1.3))
pViolinBack <- geom_boxplot(aes(color=cxMTsat.status,alpha=0.6),
                            linetype="solid",linewidth=0.6, outlier.shape = NA, coef = 0)
pViolinFront <- geom_violin(aes(color=cxMTsat.status),trim=TRUE,
                            linewidth=1.0, fill="transparent")
p4 <- ggplot(tmp, aes(cxMTsat.status,count)) + pViolinBack + pViolinFront +
  pPointB + pScaleFill + pScaleClr + pthemeB + labs(y="Rim-active ws lesion count",title="D") + add_stats(df.form)


plot2 <- p1+p2+p3+p4+plot_layout(ncol=4,nrow=1)
ggsave(sprintf("%sSupplementaryFigure6.jpg",outDir),plot2,width=21,height=8,units="cm",dpi=1200)

# #3) Rim-active ws lesion volume (Male/Female) - Male sample size is too small
# tmp <- RimActive.Vol[RimActive.Vol$region=="ws",] %>% #CHANGE
#   dplyr::group_by(subject.name) %>%
#   dplyr::summarise(mean.Vol = mean(Vol), sex=first(sex))
# tmp$mean.logVol <- log(tmp$mean.Vol)
# df.form <- tmp %>%
#   rstatix::t_test(mean.logVol ~ sex) %>%
#   rstatix::add_xy_position(x = "sex") # important for positioning!
# df.form$y.position <- df.form$y.position + 0.10
# #MUST DO: Perform normality test, then linear regression adjusting for age and sex, then get P-value and input here
# df.form$p.signif <- "**" #P=0.009 (Linear regression adjusting for age)
# ggplot(tmp, aes(sex,mean.logVol)) + geom_boxplot() + geom_jitter() + 
#   labs(y="Rim-active ws lesion volume",title="A") + add_stats(df.form)
# 
# #4) Rim-active ws lesion count (Male/Female) - insignificant result  
# tmp <- RimActive.Count[RimActive.Count$region=="ws",]
# df.form <- tmp %>%
#   rstatix::t_test(count ~ sex) %>%
#   rstatix::add_xy_position(x = "sex") # important for positioning!
# df.form$y.position <- df.form$y.position + 0.10
# #MUST DO: Perform normality test, then linear regression adjusting for age and sex, then get P-value and input here
# df.form$p.signif <- "**" #P=0.009 (Linear regression adjusting for age)
# ggplot(tmp, aes(sex,count)) + geom_boxplot() + geom_jitter()+ 
#   labs(y="Rim-active ws lesion count",title="B") + add_stats(df.form)
# Stats for Table #1 (Demographics) ---------------------------------------
library(ggpubr)
library(ggsignif)
library(emmeans)
library(ggprism)
library(R.matlab)

#Load the .RData with LABs for Table 1 summary measures and statistics, since LABs were 
#ONLY omitted in the case of TSPO PET WML phenotyping
mainDir <- "/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/"
outDir <- "/export02/data/risa/DODManuscript2_Revisions/"
load(sprintf("%sfinalizedData_withLABs.RData",mainDir))

#1) Age, mean (std)
tmp <- data.frame(age=explaVars$age, #CHANGE
                  subj.type=explaVars$subj.type,
                  clinicalStatus=explaVars$clinicalStatus)
useMetric = "mean_sd"
tmp %>% 
  group_by(clinicalStatus) %>% #CHANGE
  get_summary_stats(type=useMetric)
tmp %>% 
  group_by(subj.type) %>% #CHANGE
  get_summary_stats(type=useMetric)

#2) EDSS and Disease Duration, median (IQR)
tmp <- data.frame(edss=explaVars$edss,
                  disDur=explaVars$diseaseDuration,
                  subj.type=explaVars$subj.type,
                  clinicalStatus=explaVars$clinicalStatus)
useMetric = "median_iqr"
tmp %>% 
  group_by(clinicalStatus) %>% #CHANGE
  get_summary_stats(type=useMetric)
tmp %>% 
  group_by(subj.type) %>% #CHANGE
  get_summary_stats(type=useMetric)

#3) Total WML volume, WBV PF %, global cortical thickness, thalamic volume
tmp <- data.frame(WMLTotalVol=explaVars$WMLTotalVol,
                  WMLTotalVol_PF=(explaVars$WMLTotalVol/explaVars$TIV),
                  WBV_PF=(explaVars$WBV/explaVars$TIV),
                  gct=explaVars$gct,
                  subj.type=explaVars$subj.type,
                  clinicalStatus=explaVars$clinicalStatus)
useMetric = "mean_sd"
tmp %>% 
  group_by(clinicalStatus) %>% #CHANGE
  get_summary_stats(type=useMetric)
tmp %>% 
  group_by(subj.type) %>% #CHANGE
  get_summary_stats(type=useMetric)

#4) Aggregated binary lesion mask in healthy controls - is this in cm3?
tmp <- (cWMLs.TotalVol$sumVol/cWMLs.TotalVol$WBV)*100
mean(tmp)
sd(tmp,na.rm=TRUE)

# volume <- readMat("/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/volume.mat")
# gc()
# dfvol <- data.frame(cortex=unlist(t(volume$corticalGM))/explaVars$TIV,
#                     cerebralWM=unlist(t(volume$cerebralWM))/explaVars$TIV,
#                     WMLs=explaVars$WMLvolume)
# dfvol <- dfvol[!rmvIdxs,]

#3) Stats: WML volume 
#FDR mcp, for pairwise comparisons 
datStat <- cbind(dfvol,explaVars)
mdl <- lm(WMLs/TIV~age+clinicalStatus,datStat)
means <- emmeans(mdl,~clinicalStatus)
contrast(means, method="pairwise", adjust="fdr")

#4) Stats: Global cortical thickness 
#FDR mcp, for pairwise comparisons 
datStat <- explaVars
mdl <- lm(gct~age+clinicalStatus,datStat)
means <- emmeans(mdl,~clinicalStatus)
contrast(means, method="pairwise", adjust="fdr")

datStat <- explaVars
mdl <- lm(gct~age+subj.type,datStat)
means <- emmeans(mdl,~subj.type)
contrast(means, method="pairwise", adjust="fdr")

#4) Stats: Whole Brain Volume PF% 
#FDR mcp, for pairwise comparisons 
datStat <- explaVars
mdl <- lm(WBV/TIV~age+clinicalStatus,datStat)
means <- emmeans(mdl,~clinicalStatus)
contrast(means, method="pairwise", adjust="fdr")

datStat <- explaVars
mdl <- lm(WBV/TIV~age+subj.type,datStat)
means <- emmeans(mdl,~subj.type)
contrast(means, method="pairwise", adjust="fdr")

#5) Stats: Fisher's exact tests for binding affinity and 
fisher.test(table(explaVars$bindAff,explaVars$clinicalStatus))

#6) Stats: Mann-Whitney U Test with FDR mcp correction 
idxHC <- explaVars$subj.type=="HC"
wilcox.test(diseaseDuration~clinicalStatus,explaVars,subset=!idxHC,exact=FALSE)

wilcox.test(edss~clinicalStatus,explaVars,subset=!idxHC,exact=FALSE)

pairwise.wilcox.test(explaVars$age, explaVars$clinicalStatus,
                     p.adjust.method = "fdr",exact=FALSE)

#7) Numbers of Active/Inactive/Rim-active WM lesions 
#Load the .RData with LABs for Table 1 summary measures and statistics, since LABs WERE 
#OMITTED in the case of TSPO PET WML phenotyping
mainDir <- "/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/"
outDir <- "/export02/data/risa/DODManuscript2_Revisions/"
load(sprintf("%sfinalizedData_noLABs.RData",mainDir))

RimActive.RRMS.count <- sum(RimActive.Count$count[RimActive.Count$clinical.status=="RRMS"])
Active.RRMS.count <- sum(Active.Count$count[Active.Count$clinical.status=="RRMS"])
Inactive.RRMS.count <- sum(Inactive.Count$count[Inactive.Count$clinical.status=="RRMS"])

RimActive.PMS.count <- sum(RimActive.Count$count[RimActive.Count$clinical.status=="PMS"])
Active.PMS.count <- sum(Active.Count$count[Active.Count$clinical.status=="PMS"])
Inactive.PMS.count <- sum(Inactive.Count$count[Inactive.Count$clinical.status=="PMS"])

PMS.count <- RimActive.PMS.count + Active.PMS.count + Inactive.PMS.count 
RRMS.count <- RimActive.RRMS.count + Active.RRMS.count + Inactive.RRMS.count 

RimActive.RRMS.count <- (RimActive.RRMS.count/RRMS.count)*100
Active.RRMS.count <- (Active.RRMS.count/RRMS.count)*100
Inactive.RRMS.count <- (Inactive.RRMS.count/RRMS.count)*100

RimActive.PMS.count <- (RimActive.PMS.count/PMS.count)*100
Active.PMS.count <- (Active.PMS.count/PMS.count)*100
Inactive.PMS.count <- (Inactive.PMS.count/PMS.count)*100
