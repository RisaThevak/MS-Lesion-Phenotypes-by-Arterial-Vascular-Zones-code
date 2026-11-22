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
library(epitools)
library(LMMstar)

mainDir <- "/export02/data/risa/01_DOD_Study/0_Scripts_StatisticalAnalysis/WMLphenotypeAnalysis/"
outDir <- "/export02/data/risa/DODManuscript2_Revisions/" #CHANGE
load(sprintf("%sfinalizedData_noLABs.RData",mainDir)) #CHANGE

#ggplot2 common components for customized plots
colors.TSPOPhenotype <- setNames(c('purple4','violet','pink'),c("RA", "A", "IA"))
colors.zone <- setNames(c('indianred2','orange','royalblue3'),levels(WMLs.Vol$region))
pScaleColor <- scale_color_manual(values=colors.zone)
pScaleFill <- scale_fill_manual(values=colors.zone)
pPoint <- geom_jitter(aes(color="gray10"),width=0.25, height=0, size=0.5)
#draw_quantiles() argument in geom_violin() should not be used as geom_violin() calculates quantiles 
#using the density distribution, where boxplot() and geom_boxplot() calculates quantiles based on 
#original data
pViolinBack <- geom_boxplot(aes(color=region,alpha=0.6),linetype="solid",linewidth=0.6, outlier.shape = NA, coef = 0)
pViolinFront <- geom_violin(aes(color=region),trim=TRUE,
                            linewidth=1.0, fill="transparent")
add_stats <- function(df,...) {
  add_pvalue(data=df, 
             xmin = "xmin", 
             xmax = "xmax",
             label = "{p.adj.signif}",
             label.size = 4.0,
             tip.length = 0.015,
             bracket.size = 0.85,...)
}
pwCol1 <- 11.38
pwCol2 <- 15.29
ph <- 17.28
pStepIncrease <- 0.11
ptheme <- theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
                axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12),
                axis.title.y = element_text(face="bold",size=12),
                # axis.title.x = element_text(face="bold",size=11),
                axis.title.x = element_blank(),
                text = element_text(family="Liberation Serif"),
                legend.position="none")
pXlabels <- scale_x_discrete(labels=c("ws", "mixed", "nws"))
pSUVRYlabels <- scale_y_continuous(limits=c(0.25,2.0), 
                                   breaks=c(0.5,1.0,1.5,2.0),
                                   oob=rescale_none, 
                                   labels=c("0.5","1.0","1.5","2.0"))
pMTsatYlabels <- scale_y_continuous(limits=c(0,3.0), 
                                   breaks=c(0,1,2,3),
                                   oob=rescale_none, 
                                   labels=c("0","1.0","2.0","3.0"))
pVolumeYlabels <- scale_y_continuous(limits=c(-3.75,4.5), 
                                     breaks=c(-3,-2,-1,0,1,2,3,4),
                                     oob=rescale_none, 
                                     labels=c("-3.0","-2.0","-1.0","0","1.0","2.0","3.0","4.0"))
pPerCoreYlabels <- scale_y_continuous(limits=c(0,100), 
                                      breaks=c(0,25,50,75,100),
                                      oob=rescale_none, 
                                      labels=c("0","25","50","75","100"))

pWMPropy <- ylim(0,100)
pthemepatientwise <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
                           axis.text.x = element_blank(),
                           axis.text.y = element_text(size=12),
                           axis.title.y = element_text(face="bold",size=12),
                           axis.title.x = element_blank(),
                           text = element_text(family="Liberation Serif"),
                           legend.position="bottom", 
                           legend.title=element_blank(),
                           legend.box.margin = margin(0,0,0,0),
                           legend.margin = margin(0,0,0,0),
                           legend.text = element_text(size=12),
                           strip.placement = "outside",   # format to look like title
                           strip.background = element_blank(),
                           strip.text = element_text(face="bold", size=13))
pCountlabels <- scale_y_continuous(limits=c(0,170), 
                                   breaks=c(0, 25, 50, 75, 100, 125, 150), 
                                   expand = expansion(mult = c(0, 0)))

pthemeProp <- theme(plot.title = element_text(hjust = 0.5, size = 11.5, face="bold"),
                     axis.text.x = element_text(size=12),
                     axis.text.y = element_text(size=12),
                     axis.title.y = element_text(face="bold",size=12),
                     axis.title.x = element_blank(),
                     text = element_text(family="Liberation Serif"),
                     legend.position="none", 
                     strip.placement = "outside",   # format to look like title
                     strip.background = element_blank(),
                     strip.text = element_text(face="bold", size=13))
pProplabels <- scale_y_continuous(limits=c(0,115),
                                   breaks=c(0,10,20,30,40,50,60,70,80,90,100), 
                                   expand = expansion(mult = c(0, 0)))

cal.Par.Resid <- function(datPlot, mmod)
{ 
  addEffects <- c()
  coefmmod <- coef(mmod)$subject.name #Make into matrix 
  for (i in 1:nrow(datPlot))
  {
    subjStr <- toString(datPlot$subject.name[i])
    if (datPlot$region[i]=="ws"){
      addEffects[i] <- coefmmod[subjStr,"(Intercept)"]
    } else if (datPlot$region[i]=="mixed"){
      addEffects[i] <- coefmmod[subjStr,"(Intercept)"] + coefmmod[subjStr,"regionmixed"] 
    }else if (datPlot$region[i]=="nws") {
      addEffects[i] <- coefmmod[subjStr,"(Intercept)"] + coefmmod[subjStr,"regionnws"]  
    }
  }
  parResid <- resid(mmod) + addEffects 
  
  return(parResid)
}

fisher.tests.mcp <- function(A)
{ 
  fisher.pvals <- numeric(length=9)
  fisher.odds <- numeric(length=9)
  fisher.odds.confInt <- list()
  
  #1
  mat <- matrix(nrow=2, ncol=2, byrow=TRUE,
                      c(A[1,1], A[1,2]+A[1,3], A[2,1], A[2,2]+A[2,3]))
  res <- fisher.test(mat, y = NULL, alternative = "two.sided",
              conf.int = TRUE, conf.level = 0.95)
  fisher.pvals[1] <- res$p.value
  fisher.odds[1] <- res$estimate
  fisher.odds.confInt[[1]] <- res$conf.int
  
  
  #2
  mat <- matrix(nrow=2, ncol=2, byrow=TRUE,
                   c(A[1,1], A[1,2]+A[1,3], A[3,1], A[3,2]+A[3,3]))
  res <- fisher.test(mat, y = NULL, alternative = "two.sided",
                     conf.int = TRUE, conf.level = 0.95)
  fisher.pvals[2] <- res$p.value
  fisher.odds[2] <- res$estimate
  fisher.odds.confInt[[2]] <- res$conf.int
  
  #3
  mat <- matrix(nrow=2, ncol=2, byrow=TRUE,
                   c(A[2,1], A[2,2]+A[2,3], A[3,1], A[3,2]+A[3,3]))
  res <- fisher.test(mat, y = NULL, alternative = "two.sided",
                     conf.int = TRUE, conf.level = 0.95)
  fisher.pvals[3] <- res$p.value
  fisher.odds[3] <- res$estimate
  fisher.odds.confInt[[3]] <- res$conf.int
  
  #4
  mat <- matrix(nrow=2, ncol=2, byrow=TRUE,
                   c(A[1,2], A[1,1]+A[1,3], A[2,2], A[2,1]+A[2,3]))
  res <- fisher.test(mat, y = NULL, alternative = "two.sided",
                     conf.int = TRUE, conf.level = 0.95)
  fisher.pvals[4] <- res$p.value
  fisher.odds[4] <- res$estimate
  fisher.odds.confInt[[4]] <- res$conf.int
  
  #5
  mat <- matrix(nrow=2, ncol=2, byrow=TRUE,
                   c(A[1,2], A[1,1]+A[1,3], A[3,2], A[3,1]+A[3,3]))
  res <- fisher.test(mat, y = NULL, alternative = "two.sided",
                     conf.int = TRUE, conf.level = 0.95)
  fisher.pvals[5] <- res$p.value
  fisher.odds[5] <- res$estimate
  fisher.odds.confInt[[5]] <- res$conf.int
  
  #6
  mat <- matrix(nrow=2, ncol=2, byrow=TRUE,
                   c(A[2,2], A[2,1]+A[2,3], A[3,2], A[3,1]+A[3,3]))
  res <- fisher.test(mat, y = NULL, alternative = "two.sided",
                     conf.int = TRUE, conf.level = 0.95)
  fisher.pvals[6] <- res$p.value
  fisher.odds[6] <- res$estimate
  fisher.odds.confInt[[6]] <- res$conf.int
  
  #7
  mat <- matrix(nrow=2, ncol=2, byrow=TRUE,
                   c(A[1,3], A[1,1]+A[1,2], A[2,3], A[2,1]+A[2,2]))
  res <- fisher.test(mat, y = NULL, alternative = "two.sided",
                     conf.int = TRUE, conf.level = 0.95)
  fisher.pvals[7] <- res$p.value
  fisher.odds[7] <- res$estimate
  fisher.odds.confInt[[7]] <- res$conf.int
  
  #8
  mat <- matrix(nrow=2, ncol=2, byrow=TRUE,
                   c(A[1,3], A[1,1]+A[1,2], A[3,3], A[3,1]+A[3,2]))
  res <- fisher.test(mat, y = NULL, alternative = "two.sided",
                     conf.int = TRUE, conf.level = 0.95)
  fisher.pvals[8] <- res$p.value
  fisher.odds[8] <- res$estimate
  fisher.odds.confInt[[8]] <- res$conf.int
  
  #9
  mat <- matrix(nrow=2, ncol=2, byrow=TRUE,
                   c(A[2,3], A[2,1]+A[2,2], A[3,3], A[3,1]+A[3,2]))
  res <- fisher.test(mat, y = NULL, alternative = "two.sided",
                     conf.int = TRUE, conf.level = 0.95)
  fisher.pvals[9] <- res$p.value
  fisher.odds[9] <- res$estimate
  fisher.odds.confInt[[9]] <- res$conf.int
  
  pvals <- p.adjust(fisher.pvals, method="fdr")
  res <- list(fisher.pvals,pvals,fisher.odds,fisher.odds.confInt)
  
  return(res)
}

#Use ML fitting method to get AIC for model 
#RATIONALE: Faraway's (2006) Extending the linear model with R (p. 156): 
#"The reason is that REML estimates the random effects by considering linear 
#combinations of the data that remove the fixed effects. 
#If these fixed effects are changed, the likelihoods of the two models will not
#be directly comparable."
#RATIONALE: Zuur et al. (2009; p. 122) suggest that "To compare models with nested
#fixed effects (but with the same random structure), ML estimation must be used and
#not REML." This indicates to me that I ought to use ML since my random effects are
#the same in both models, but my fixed effects differ. [Zuur et al. 2009. 
#Mixed Effect Models and Extensions in Ecology with R. Springer.]

#Sage artwork guidelines: asterized based files (i.e. with .tiff or .jpeg extension) require a resolution 
#of at least 300 dpi (dots per inch). Line art should be supplied with a minimum resolution of 800 dpi. 

# Figure 1: WM lesion count ----------------------
#Ordinary means and standard deviations 
tmp <- data.frame(count=WMLs.Count$count,
                  region=WMLs.Count$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  get_summary_stats(type="mean_sd")

tmp <- data.frame(count=WMLs.Count$count,
                  clinical.status=WMLs.Count$clinical.status,
                  region=WMLs.Count$region) #CUSTOM
tmp %>% 
  group_by(clinical.status,region) %>%
  get_summary_stats(type="mean_sd")

#Plot 1b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot1b <- WMLs.Count #CUSTOM
mmod <- lmer(count ~ (1 | subject.name) + region, datPlot1b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot1b,count~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot1b, c('region','count')) #CUSTOM
temp$region <- factor(datPlot1b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(count ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val1b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 1c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot1c <- WMLs.Count[WMLs.Count$clinical.status=="RRMS",] #CUSTOM
mmod <- lmer(count ~ (1 | subject.name) + region, datPlot1c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot1c,count~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot1c, c('region','count')) #CUSTOM
temp$region <- factor(datPlot1c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(count ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val1c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 1d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot1d <- WMLs.Count[WMLs.Count$clinical.status=="PMS",] #CUSTOM
mmod <- lmer(count ~ (1 | subject.name) + region, datPlot1d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot1d,count~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot1d, c('region','count')) #CUSTOM
temp$region <- factor(datPlot1d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(count ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val1d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#MANUAL
df_p_val1b$y.position <- df_p_val1b$y.position-10
df_p_val1d$y.position <- df_p_val1d$y.position-10

#Plotting
nonFirstyScale <- scale_y_continuous(limits=c(0,85), 
                                     breaks=c(0,10,20,30,40,50,60,70,80),
                                     oob=rescale_none, 
                                     labels=c("","","","","","","","",""),
                                     name=element_blank())
plot1b <- ggplot(datPlot1b, aes(region,count))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels +
  labs(y = 'Lesion count', title = "MS") + #CUSTOM
  scale_y_continuous(limits=c(0,85), 
                     breaks=c(0,10,20,30,40,50,60,70,80))
plot1c <- ggplot(datPlot1c, aes(region,count))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "RRMS") + nonFirstyScale
plot1d <- ggplot(datPlot1d, aes(region,count))  +
  pPoint + pViolinBack + pViolinFront +
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Progressive MS") + nonFirstyScale

if (nrow(df_p_val1b)!=0) { plot1b <- plot1b + add_stats(df_p_val1b) }
if (nrow(df_p_val1c)!=0) { plot1c <- plot1c + add_stats(df_p_val1c) }
if (nrow(df_p_val1d)!=0) { plot1d <- plot1d + add_stats(df_p_val1d) }

p1 <- plot1b + plot1c + plot1d + plot_layout(ncol=3)
# ggsave(sprintf("%sFigure1.jpg",outDir),p1,width=0.75*(pwCol1+pwCol2),height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot1c, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(count ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot1c, region <- relevel(region, ref = "nws"))
# mmod <- lmer(count ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(count ~ (1 | subject.name) + region, datPlot1c,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# Figure 2: WML Volume ----------------------
#Ordinary means and standard deviations 
tmp <- data.frame(Vol=WMLs.Vol$logVol,
                  region=WMLs.Vol$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  get_summary_stats(type="mean_sd")

tmp <- data.frame(Vol=WMLs.Vol$logVol,
                  clinical.status=WMLs.Vol$clinical.status,
                  region=WMLs.Vol$region) #CUSTOM
tmp %>% 
  group_by(clinical.status,region) %>%
  get_summary_stats(type="mean_sd")

#Plot 2b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot2b <- WMLs.Vol #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot2b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot2b,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot2b, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot2b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val2b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
  
#Plot 2c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot2c <- WMLs.Vol[WMLs.Vol$clinical.status=="RRMS",] #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot2c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot2c,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot2c, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot2c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val2c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 2d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot2d <- WMLs.Vol[WMLs.Vol$clinical.status=="PMS",] #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot2d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot2d,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot2d, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot2d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val2d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#MANUAL
df_p_val2b$y.position[3] <- df_p_val2b$y.position[3] + 0.2
df_p_val2d$y.position[2] <- df_p_val2d$y.position[2] + 0.2

# Plotting 
nonFirstyScale <- scale_y_continuous(limits=c(-3.75,4.5), 
                                     breaks=c(-3,-2,-1,0,1,2,3,4),
                                     oob=rescale_none, 
                                     labels=c("","","","","","","",""),
                                     name=element_blank())
plot2b <- ggplot(datPlot2b, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels +  
  labs(y = ~paste(bold("Log-transformed lesion volume ("),bold("cm"^"3"),bold(") "))) + pVolumeYlabels 
plot2c <- ggplot(datPlot2c, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale
plot2d <- ggplot(datPlot2d, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale

if (nrow(df_p_val2b)!=0) { plot2b <- plot2b + add_stats(df_p_val2b) } 
if (nrow(df_p_val2c)!=0) { plot2c <- plot2c + add_stats(df_p_val2c) } 
if (nrow(df_p_val2d)!=0) { plot2d <- plot2d + add_stats(df_p_val2d) } 

p2 <- plot2b + plot2c + plot2d + plot_layout(ncol=3)
# ggsave(sprintf("%sFigure2.jpg",outDir),p2,width=0.75*(pwCol1+pwCol2),height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot2d, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(logVol ~ region + (1 | subject.name) + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot2d, region <- relevel(region, ref = "nws"))
# mmod <- lmer(logVol ~ region + (1 | subject.name) + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot2d,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# Figure 3: WM lesion count per patient -----------------------------------
p2b <- ggplot(WMLs.Count, aes(x=reorder(subject.name,count),y=count,fill=region)) + 
  geom_bar(stat="identity",position="stack",linewidth=0.4, width=0.65, color="black") + 
  pScaleFill + pthemepatientwise + pCountlabels + labs(y = 'WM lesion count') + 
  facet_wrap(~clinical.status,nrow=1,scales="free_x", strip.position = "bottom") 

# FIGURE A: Composite Figure ----------------------------------------------
figureA <- wrap_elements(full=p1)+wrap_elements(full=p2)+
  wrap_elements(p2b)+plot_layout(nrow=3,ncol=1)
ggsave(sprintf("%sFigureA.jpg",outDir),figureA,
       width=0.75*(pwCol1+pwCol2),height=26.5,units="cm",dpi=1200)

# EMPTY SPACE -------------------------------------------------------------
# Figure 4: WML MTsat ----------------------
#Ordinary means and standard deviations 
tmp <- data.frame(mtsat=Fin.WMLs.mtsat$mtsat,
                  subj.type=Fin.WMLs.mtsat$subj.type,
                  region=Fin.WMLs.mtsat$region) #CUSTOM
tmp %>% 
  group_by(subj.type,region) %>%
  get_summary_stats(type="mean_sd")

tmp <- data.frame(mtsat=Fin.WMLs.mtsat$mtsat,
                  clinical.status=Fin.WMLs.mtsat$clinical.status,
                  region=Fin.WMLs.mtsat$region) #CUSTOM
tmp %>% 
  group_by(clinical.status,region) %>%
  get_summary_stats(type="mean_sd")

#Plot 4a GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot4a <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$subj.type=="HC",] #CUSTOM
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, datPlot4a,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot4a,mtsat~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot4a, c('region','mtsat')) #CUSTOM
temp$region <- factor(datPlot4a$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsat ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val4a <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot4a$parResid <- cal.Par.Resid(datPlot4a,mmod)

#Plot 4b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot4b <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$subj.type=="MS",] #CUSTOM
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, datPlot4b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot4b,mtsat~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot4b, c('region','mtsat')) #CUSTOM
temp$region <- factor(datPlot4b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsat ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val4b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot4b$parResid <- cal.Par.Resid(datPlot4b,mmod)

#Plot 4c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot4c <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$clinical.status=="RRMS",] #CUSTOM
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, datPlot4c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot4c,mtsat~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot4c, c('region','mtsat')) #CUSTOM
temp$region <- factor(datPlot4c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsat ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val4c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot4c$parResid <- cal.Par.Resid(datPlot4c,mmod)

#Plot 4d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot4d <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$clinical.status=="PMS",] #CUSTOM
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, datPlot4d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot4d,mtsat~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot4d, c('region','mtsat')) #CUSTOM
temp$region <- factor(datPlot4d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsat ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val4d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot4d$parResid <- cal.Par.Resid(datPlot4d,mmod)

#MANUAL CHANGES
df_p_val4a <- data.frame()
df_p_val4b <- data.frame()
df_p_val4c <- data.frame()
df_p_val4d <- data.frame()

# Plotting 
nonFirstyScale <- scale_y_continuous(limits=c(0,3.0), 
                     breaks=c(0,1,2,3),
                     oob=rescale_none, 
                     labels=c("","","",""),
                     name=element_blank())
plot4a <- ggplot(datPlot4a, aes(region,parResid)) + pPoint +
  pViolinBack + pViolinFront +
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + pMTsatYlabels + 
  labs(y = "MTsat (p.u.)") #CUSTOM
plot4b <- ggplot(datPlot4b, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale
plot4c <- ggplot(datPlot4c, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale
plot4d <- ggplot(datPlot4d, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale

if (nrow(df_p_val4a)!=0) { plot4a <- plot4a + add_stats(df_p_val4a) }
if (nrow(df_p_val4b)!=0) { plot4b <- plot4b + add_stats(df_p_val4b) }
if (nrow(df_p_val4c)!=0) { plot4c <- plot4c + add_stats(df_p_val4c) }
if (nrow(df_p_val4d)!=0) { plot4d <- plot4d + add_stats(df_p_val4d) }

p4 <- plot4a + plot4b + plot4c + plot4d + plot_layout(ncol=4,nrow=1)
# ggsave(sprintf("%sFigure4.jpg",outDir),p4,width=pwCol1+pwCol2,height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot4c, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot4c, region <- relevel(region, ref = "nws"))
# mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(mtsat ~ (1 | subject.name) + WML.volume + region, datPlot4c,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# Figure 6: WML MTsat (per mm3) ----------------------
#Ordinary means and standard deviations 
tmp <- data.frame(mtsatpervol=Fin.WMLs.mtsat$mtsatpervol,
                  subj.type=Fin.WMLs.mtsat$subj.type,
                  region=Fin.WMLs.mtsat$region) #CUSTOM
tmp %>% 
  group_by(subj.type,region) %>%
  get_summary_stats(type="mean_sd")

tmp <- data.frame(mtsatpervol=Fin.WMLs.mtsat$mtsatpervol,
                  clinical.status=Fin.WMLs.mtsat$clinical.status,
                  region=Fin.WMLs.mtsat$region) #CUSTOM
tmp %>% 
  group_by(clinical.status,region) %>%
  get_summary_stats(type="mean_sd")

#Plot 6b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot6b <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$subj.type=="MS",] #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot6b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot6b,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot6b, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot6b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val6b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot6b$parResid <- cal.Par.Resid(datPlot6b,mmod)

#Plot 6c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot6c <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$clinical.status=="RRMS",] #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot6c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot6c,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot6c, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot6c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val6c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot6c$parResid <- cal.Par.Resid(datPlot6c,mmod)

#Plot 6d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot6d <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$clinical.status=="PMS",] #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot6d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot6d,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot6d, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot6d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val6d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot6d$parResid <- cal.Par.Resid(datPlot6d,mmod)

#To fit stat brackets within y-scale 
df_p_val6b$y.position[2] <- df_p_val6b$y.position[1] + 0.003
df_p_val6c$y.position[2] <- df_p_val6c$y.position[1] + 0.003
df_p_val6d$y.position[2] <- df_p_val6d$y.position[1] + 0.003

# Plotting 
pMTsatYlabelsv2 <- scale_y_continuous(limits=c(0,0.060), 
                                    breaks=c(0,0.015,0.030,0.045,0.060),
                                    oob=rescale_none, 
                                    labels=c("0","0.015","0.030","0.045","0.060"))
nonFirstyScale <- scale_y_continuous(limits=c(0,0.060), 
                                     breaks=c(0,0.015,0.030,0.045,0.060),
                                     oob=rescale_none, 
                                     labels=c("","","","",""),
                                     name=element_blank())
plot6b <- ggplot(datPlot6b, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + labs(y=expression(bold(paste('MTsat per unit volume (p.u./',mm^'3',')')))) + pMTsatYlabelsv2
plot6c <- ggplot(datPlot6c, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale
plot6d <- ggplot(datPlot6d, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale

if (nrow(df_p_val6b)!=0) { plot6b <- plot6b + add_stats(df_p_val6b) }
if (nrow(df_p_val6c)!=0) { plot6c <- plot6c + add_stats(df_p_val6c) }
if (nrow(df_p_val6d)!=0) { plot6d <- plot6d + add_stats(df_p_val6d) }

p6 <- plot_spacer() + plot6b + plot6c + plot6d + plot_layout(ncol=4,nrow=1)
ggsave(sprintf("%sSuppFigure_WML_MTsatpermm3.jpg",outDir),p6,width=pwCol1+pwCol2,height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot6d, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot6d, region <- relevel(region, ref = "nws"))
# mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot6d,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# Figure 5: WML SUVR ----------------------
#Ordinary means and standard deviations 
tmp <- data.frame(SUVR=Fin.WMLs.SUVR$SUVR,
                  subj.type=Fin.WMLs.SUVR$subj.type,
                  region=Fin.WMLs.SUVR$region) #CUSTOM
tmp %>% 
  group_by(subj.type,region) %>%
  get_summary_stats(type="mean_sd")

tmp <- data.frame(SUVR=Fin.WMLs.SUVR$SUVR,
                  clinical.status=Fin.WMLs.SUVR$clinical.status,
                  region=Fin.WMLs.SUVR$region) #CUSTOM
tmp %>% 
  group_by(clinical.status,region) %>%
  get_summary_stats(type="mean_sd")

#Plot 5a GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot5a <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$subj.type=="HC",] #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot5a,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot5a,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot5a, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot5a$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val5a <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot5a$parResid <- cal.Par.Resid(datPlot5a,mmod)

#Plot 5b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot5b <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$subj.type=="MS",] #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot5b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot5b,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot5b, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot5b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val5b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot5b$parResid <- cal.Par.Resid(datPlot5b,mmod)

#Plot 5c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot5c <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$clinical.status=="RRMS",] #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot5c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot5c,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot5c, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot5c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val5c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot5c$parResid <- cal.Par.Resid(datPlot5c,mmod)

#Plot 5d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot5d <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$clinical.status=="PMS",] #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot5d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot5d,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot5d, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot5d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val5d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot5d$parResid <- cal.Par.Resid(datPlot5d,mmod)

#MANUAL CHANGES
df_p_val5a <- data.frame()
df_p_val5b <- df_p_val5b[c(2,3),]
df_p_val5c <- df_p_val5c[c(2,3),]
df_p_val5d <- df_p_val5d[2,]

#To make stat brackets fall within y-axis range
df_p_val5b$y.position[1] <- df_p_val5b$y.position[1]-0.30
df_p_val5b$y.position[2] <- df_p_val5b$y.position[2]-0.24
df_p_val5c$y.position[1] <- df_p_val5c$y.position[1]-0.30
df_p_val5c$y.position[2] <- df_p_val5c$y.position[2]-0.24
df_p_val5d$y.position <- df_p_val5d$y.position-0.15

# Plotting 
nonFirstyScale <- scale_y_continuous(limits=c(0.25,2.0), 
                                     breaks=c(0.5,1.0,1.5,2.0),
                                     oob=rescale_none, 
                                     labels=c("","","",""),
                                     name=element_blank())
plot5a <- ggplot(datPlot5a, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + pSUVRYlabels + 
  labs(y = "[11C]PBR28 PET SUVR", title = "Healthy Control") #CUSTOM
plot5b <- ggplot(datPlot5b, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + labs(y = "", title = "MS") + nonFirstyScale 
plot5c <- ggplot(datPlot5c, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "RRMS") + nonFirstyScale 
plot5d <- ggplot(datPlot5d, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Progressive MS") + nonFirstyScale 

if (nrow(df_p_val5a)!=0) { plot5a <- plot5a + add_stats(df_p_val5a) }
if (nrow(df_p_val5b)!=0) { plot5b <- plot5b + add_stats(df_p_val5b) }
if (nrow(df_p_val5c)!=0) { plot5c <- plot5c + add_stats(df_p_val5c) }
if (nrow(df_p_val5d)!=0) { plot5d <- plot5d + add_stats(df_p_val5d) }

p5 <- plot5a + plot5b + plot5c + plot5d + plot_layout(ncol=4,nrow=1)
# ggsave(sprintf("%sFigure5.jpg",outDir),p5,width=pwCol1+pwCol2,height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot5d, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot5d, region <- relevel(region, ref = "nws"))
# mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]


# FIGURE B: Composite Figure ----------------------------------------------
# figureB <- wrap_elements(full=p4)+wrap_elements(full=p5)+plot_layout(nrow=2,ncol=1)
# ggsave(sprintf("%sFigureB.jpg",outDir),figureB,
#        width=pwCol1+pwCol2,height=17.28,units="cm",dpi=1200)

figureB <- wrap_elements(full=p5)+wrap_elements(full=p4)+
  wrap_elements(full=p6)+plot_layout(nrow=3,ncol=1,heights=c(2,1.5,2))
ggsave(sprintf("%sFigure5.jpg",outDir),figureB,
       width=0.75*(pwCol1+pwCol2),height=23,units="cm",dpi=1200)

#Try only area between 3rd and 4th quaritle of each violin plot a transparent white! 
# EMPTY SPACE -------------------------------------------------------------


# Figure 19: TSPO WML Phenotype Proportions (By Count), Violin Plots ------------------------------------
#Switch between "RimActive", "Active" and "Inactive"
#Ordinary means and standard deviations 
tmp <- data.frame(Prop=Inactive.Count$Prop,
                  region=Inactive.Count$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  # get_summary_stats(type="median_iqr")
  get_summary_stats(type="mean_sd")

#Plot 19b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot19b <- RimActive.Count #CUSTOM
mmod <- lmer(Prop ~ (1 | subject.name) + region + age, datPlot19b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot19b,Prop~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot19b, c('region','Prop')) #CUSTOM
temp$region <- factor(datPlot19b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(Prop ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val19b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 19c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot19c <- Active.Count #CUSTOM
mmod <- lmer(Prop ~ (1 | subject.name) + region + age, datPlot19c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot19c,Prop~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot19c, c('region','Prop')) #CUSTOM
temp$region <- factor(datPlot19c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(Prop ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val19c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 19d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot19d <- Inactive.Count #CUSTOM
mmod <- lmer(Prop ~ (1 | subject.name) + region + age, datPlot19d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot19d,Prop~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot19d, c('region','Prop')) #CUSTOM
temp$region <- factor(datPlot19d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(Prop ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val19d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#MANUAL #CHANGE
# df_p_val19c$y.position <- df_p_val19d$y.position 
df_p_val19c$y.position[2] <- df_p_val19c$y.position[2] + 6
# df_p_val19d$y.position[2] <- df_p_val19d$y.position[2] + 2

#Plotting
nonFirstyScale <- scale_y_continuous(limits=c(0,110), 
                                     breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                                     oob=rescale_none, 
                                     labels=c("","","","","","","","","","",""),
                                     name=element_blank())
plot19b <- ggplot(datPlot19b, aes(region,Prop))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels +
  labs(y = 'Proportion of WM lesions (%)', title = "Rim-Active") + #CUSTOM
  scale_y_continuous(limits=c(0,110), 
                     breaks=c(0,10,20,30,40,50,60,70,80,90,100))
plot19c <- ggplot(datPlot19c, aes(region,Prop))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Active") + nonFirstyScale
plot19d <- ggplot(datPlot19d, aes(region,Prop))  +
  pPoint + pViolinBack + pViolinFront +
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Inactive") + nonFirstyScale

if (nrow(df_p_val19b)!=0) { plot19b <- plot19b + add_stats(df_p_val19b) }
if (nrow(df_p_val19c)!=0) { plot19c <- plot19c + add_stats(df_p_val19c) }
if (nrow(df_p_val19d)!=0) { plot19d <- plot19d + add_stats(df_p_val19d) }

p19 <- plot19b + plot19c + plot19d + plot_layout(ncol=3)
# ggsave(sprintf("%sFigure19.jpg",outDir),p19,width=0.75*(pwCol1+pwCol2),height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot19d, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(Prop ~ (1 | subject.name) + region + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot19d, region <- relevel(region, ref = "nws"))
# mmod <- lmer(Prop ~ (1 | subject.name) + region + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]


# Figure 20: TSPO WML Phenotype Volumes -----------------------------------
#Switch between "RimActive", "Active" and "Inactive"
#Ordinary means and standard deviations 
tmp <- data.frame(Vol=Inactive.Vol$logVol,
                  region=Inactive.Vol$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  get_summary_stats(type="mean_sd")

#Plot 20b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot20b <- RimActive.Vol #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot20b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot20b,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot20b, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot20b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val20b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 20c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot20c <- Active.Vol #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot20c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot20c,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot20c, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot20c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val20c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 20d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot20d <- Inactive.Vol #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot20d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot20d,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot20d, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot20d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
# df.form <- temp %>%
#   rstatix::t_test(logVol ~ region) %>% #CUSTOM
#   rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val20d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#MANUAL #CHANGE
df_p_val20b$y.position[2] <- df_p_val20b$y.position[2] + 0.5
df_p_val20c$y.position[2] <- df_p_val20c$y.position[2] + 0.3
df_p_val20c$y.position[3] <- df_p_val20c$y.position[3] + 0.7

# Plotting 
nonFirstyScale <- scale_y_continuous(limits=c(-3.75,4.5), 
                                     breaks=c(-3,-2,-1,0,1,2,3,4),
                                     oob=rescale_none, 
                                     labels=c("","","","","","","",""),
                                     name=element_blank())
plot20b <- ggplot(datPlot20b, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels +
  labs(y = ~paste(bold("Log-transformed lesion volume ("),bold("cm"^"3"),bold(") "))) + pVolumeYlabels
plot20c <- ggplot(datPlot20c, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "") + nonFirstyScale
plot20d <- ggplot(datPlot20d, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "") + nonFirstyScale

if (nrow(df_p_val20b)!=0) { plot20b <- plot20b + add_stats(df_p_val20b) } 
if (nrow(df_p_val20c)!=0) { plot20c <- plot20c + add_stats(df_p_val20c) } 
if (nrow(df_p_val20d)!=0) { plot20d <- plot20d + add_stats(df_p_val20d) }

p20 <- plot20b + plot20c + plot20d + plot_layout(ncol=3)
# ggsave(sprintf("%sFigure20.jpg",outDir),p20,width=0.75*(pwCol1+pwCol2),height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot20c, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(logVol ~ (1 | subject.name) + region + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot20c, region <- relevel(region, ref = "nws"))
# mmod <- lmer(logVol ~ (1 | subject.name) + region + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]


# Figure 21: TSPO WML Phenotype Proportions (By Count), Stacked Bars -------------------------------
RimActive.ws.count <- sum(RimActive.Count$count[RimActive.Count$region=="ws"])
Active.ws.count <- sum(Active.Count$count[Active.Count$region=="ws"])
Inactive.ws.count <- sum(Inactive.Count$count[Inactive.Count$region=="ws"])

RimActive.mixed.count <- sum(RimActive.Count$count[RimActive.Count$region=="mixed"])
Active.mixed.count <- sum(Active.Count$count[Active.Count$region=="mixed"])
Inactive.mixed.count <- sum(Inactive.Count$count[Inactive.Count$region=="mixed"])

RimActive.nws.count <- sum(RimActive.Count$count[RimActive.Count$region=="nws"])
Active.nws.count <- sum(Active.Count$count[Active.Count$region=="nws"])
Inactive.nws.count <- sum(Inactive.Count$count[Inactive.Count$region=="nws"])

#Not using original WML counts, since for TSPO phenotyping, lesions with rim volume < 27mm3 were discounted 
ws.count <- RimActive.ws.count + Active.ws.count + Inactive.ws.count 
mixed.count <- RimActive.mixed.count + Active.mixed.count + Inactive.mixed.count 
nws.count <- RimActive.nws.count + Active.nws.count + Inactive.nws.count 

prop1 <- c(RimActive.ws.count/ws.count, Active.ws.count/ws.count, Inactive.ws.count/ws.count, 
          RimActive.mixed.count/mixed.count, Active.mixed.count/mixed.count, Inactive.mixed.count/mixed.count, 
          RimActive.nws.count/nws.count, Active.nws.count/nws.count, Inactive.nws.count/nws.count)*100 
df1 <- data.frame(region=factor(c("ws", "ws", "ws", "mixed", "mixed", "mixed", "nws", "nws", "nws"),levels=c("ws", "mixed", "nws")),
                 TSPOphenotype=factor(c("RA", "A", "IA", "RA", "A", "IA", "RA", "A", "IA"),levels=c("IA", "A", "RA")),
                 prop=prop1)
p21 <- ggplot(df1, aes(x=region,y=prop1,fill=TSPOphenotype)) + 
  geom_bar(stat="identity",position="stack",color="black",linewidth=0.2,width=0.7) + 
  scale_fill_manual(values=colors.TSPOPhenotype) + 
  labs(y="Proportion of total lesion count (%)",title="") + 
  pProplabels + pthemeProp

# Figure 22: TSPO WML Phenotype Proportions (By Total Volume), Stacked Bars -------------------------------
RimActive.ws.Vol <- sum(RimActive.Vol$Vol[RimActive.Vol$region=="ws"])
Active.ws.Vol <- sum(Active.Vol$Vol[Active.Vol$region=="ws"])
Inactive.ws.Vol <- sum(Inactive.Vol$Vol[Inactive.Vol$region=="ws"])

RimActive.mixed.Vol <- sum(RimActive.Vol$Vol[RimActive.Vol$region=="mixed"])
Active.mixed.Vol <- sum(Active.Vol$Vol[Active.Vol$region=="mixed"])
Inactive.mixed.Vol <- sum(Inactive.Vol$Vol[Inactive.Vol$region=="mixed"])

RimActive.nws.Vol <- sum(RimActive.Vol$Vol[RimActive.Vol$region=="nws"])
Active.nws.Vol <- sum(Active.Vol$Vol[Active.Vol$region=="nws"])
Inactive.nws.Vol <- sum(Inactive.Vol$Vol[Inactive.Vol$region=="nws"])

#Not using original WML volumes, since for TSPO phenotyping, lesions with rim volume < 27mm3 were discounted 
ws.Vol <- RimActive.ws.Vol + Active.ws.Vol + Inactive.ws.Vol 
mixed.Vol <- RimActive.mixed.Vol + Active.mixed.Vol + Inactive.mixed.Vol 
nws.Vol <- RimActive.nws.Vol + Active.nws.Vol + Inactive.nws.Vol 

prop2 <- c(RimActive.ws.Vol/ws.Vol, Active.ws.Vol/ws.Vol, Inactive.ws.Vol/ws.Vol, 
          RimActive.mixed.Vol/mixed.Vol, Active.mixed.Vol/mixed.Vol, Inactive.mixed.Vol/mixed.Vol, 
          RimActive.nws.Vol/nws.Vol, Active.nws.Vol/nws.Vol, Inactive.nws.Vol/nws.Vol)*100 
df2 <- data.frame(region=factor(c("ws", "ws", "ws", "mixed", "mixed", "mixed", "nws", "nws", "nws"),levels=c("ws", "mixed", "nws")),
                 TSPOphenotype=factor(c("RA", "A", "IA", "RA", "A", "IA", "RA", "A", "IA"),levels=c("IA", "A", "RA")),
                 prop=prop2)
p22 <- ggplot(df2, aes(x=region,y=prop2,fill=TSPOphenotype)) + 
  geom_bar(stat="identity",position="stack",color="black",linewidth=0.2,width=0.7) + 
  scale_fill_manual(values=colors.TSPOPhenotype) + 
  labs(y="Proportion of total lesion volume (%)",title="") + 
  pProplabels + pthemeProp 
 
# Figure 23: TSPO WML Phenotype SUVR -----------------------------------
#Switch between "RimActive", "Active" and "Inactive"
#Ordinary means and standard deviations 
tmp <- data.frame(SUVR=Active.SUVR$SUVR,
                  region=Active.SUVR$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  get_summary_stats(type="mean_sd")

#Plot 23b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot23b <- RimActive.SUVR #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot23b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot23b,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot23b, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot23b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val23b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot23b$parResid <- cal.Par.Resid(datPlot23b,mmod)

#Plot 23c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot23c <- Active.SUVR #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot23c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot23c,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot23c, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot23c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val23c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot23c$parResid <- cal.Par.Resid(datPlot23c,mmod)

#Plot 23d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot23d <- Inactive.SUVR #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot23d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot23d,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot23d, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot23d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
# df.form <- temp %>%
#   rstatix::t_test(SUVR ~ region) %>% #CUSTOM
#   rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val23d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot23d$parResid <- cal.Par.Resid(datPlot23d,mmod)

#MANUAL #CHANGE
#To make stat brackets fall within y-axis range
df_p_val23b$y.position <- NA
df_p_val23c$y.position[1] <- NA
df_p_val23c$y.position[2] <- df_p_val23c$y.position[2]-0.27
df_p_val23c$y.position[3] <- df_p_val23c$y.position[3]-0.17

# Plotting 
nonFirstyScale <- scale_y_continuous(limits=c(0.25,2.0), 
                                     breaks=c(0.5,1.0,1.5,2.0),
                                     oob=rescale_none, 
                                     labels=c("","","",""),
                                     name=element_blank())
plot23b <- ggplot(datPlot23b, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + pSUVRYlabels + 
  labs(y = "[11C]PBR28 PET SUVR", title = "Rim-Active") #CUSTOM
plot23c <- ggplot(datPlot23c, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Active") + nonFirstyScale
plot23d <- ggplot(datPlot23d, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Inactive") + nonFirstyScale

if (nrow(df_p_val23b)!=0) { plot23b <- plot23b + add_stats(df_p_val23b) }
if (nrow(df_p_val23c)!=0) { plot23c <- plot23c + add_stats(df_p_val23c) }
if (nrow(df_p_val23d)!=0) { plot23d <- plot23d + add_stats(df_p_val23d) }

p23 <- plot23b + plot23c + plot23d + plot_layout(ncol=3,nrow=1)
ggsave(sprintf("%sResponse_Figure_3.jpg",outDir),p23,width=pwCol1+pwCol2,height=ph,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot23c, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(SUVR ~ (1 | subject.name) + WML.volume + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot23c, region <- relevel(region, ref = "nws"))
# mmod <- lmer(SUVR ~ (1 | subject.name) + WML.volume + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)

# Figure 24: TSPO WML Phenotype MTsat (per mm3) -----------------------------------
#Switch between "RimActive", "Active" and "Inactive"
#Ordinary means and standard deviations 
tmp <- data.frame(mtsatpervol=Inactive.mtsat$mtsatpervol,
                  region=Inactive.mtsat$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  get_summary_stats(type="mean_sd")

#Plot 24b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot24b <- RimActive.mtsat #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot24b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot24b,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot24b, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot24b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val24b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot24b$parResid <- cal.Par.Resid(datPlot24b,mmod)

#Plot 24c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot24c <- Active.mtsat #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot24c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot24c,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot24c, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot24c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val24c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot24c$parResid <- cal.Par.Resid(datPlot24c,mmod)

#Plot 24d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot24d <- Inactive.mtsat #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot24d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot24d,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot24d, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot24d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
# df.form <- temp %>%
#   rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
#   rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val24d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot24d$parResid <- cal.Par.Resid(datPlot24d,mmod)

#To fit stat brackets within y-scale 
df_p_val24b$y.position[1] <- df_p_val24b$y.position[1] - 0.003
df_p_val24b$y.position[2] <- df_p_val24b$y.position[1] + 0.004
df_p_val24c$y.position[2] <- df_p_val24c$y.position[1] + 0.004

# Plotting 
pMTsatYlabelsv2 <- scale_y_continuous(limits=c(0,0.060), 
                                      breaks=c(0,0.015,0.030,0.045,0.060),
                                      oob=rescale_none, 
                                      labels=c("0","0.015","0.030","0.045","0.060"))
nonFirstyScale <- scale_y_continuous(limits=c(0,0.060), 
                                     breaks=c(0,0.015,0.030,0.045,0.060),
                                     oob=rescale_none, 
                                     labels=c("","","","",""),
                                     name=element_blank())
plot24b <- ggplot(datPlot24b, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + pMTsatYlabelsv2 + 
  labs(y=expression(bold(paste('MTsat per unit volume (p.u./',mm^'3',')')))) #CUSTOM
plot24c <- ggplot(datPlot24c, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale
plot24d <- ggplot(datPlot24d, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale

if (nrow(df_p_val24b)!=0) { plot24b <- plot24b + add_stats(df_p_val24b) }
if (nrow(df_p_val24c)!=0) { plot24c <- plot24c + add_stats(df_p_val24c) }
if (nrow(df_p_val24d)!=0) { plot24d <- plot24d + add_stats(df_p_val24d) }

p24 <- plot24b + plot24c + plot24d + plot_layout(ncol=3,nrow=1)
ggsave(sprintf("%sSuppFigure_TSPOWML_MTsatpermm3.jpg",outDir),p24,width=pwCol1+pwCol2,height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot24c, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot24c, region <- relevel(region, ref = "nws"))
# mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)

# Figure 25: TSPO WML Phenotype PerCore -----------------------------------
#Switch between "RimActive", "Active" and "Inactive"
#Ordinary means and standard deviations 

#Plot 25b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot25b <- RimActive.PerCore #CUSTOM
mmod <- lmer(PerCore ~ (1 | subject.name) + region + WML.volume, datPlot25b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot25b,PerCore~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot25b, c('region','PerCore')) #CUSTOM
temp$region <- factor(datPlot25b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(PerCore ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val25b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
# datPlot25b$parResid <- cal.Par.Resid(datPlot25b,mmod)

#Plot 25c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot25c <- Active.PerCore #CUSTOM
mmod <- lmer(PerCore ~ (1 | subject.name) + region + WML.volume, datPlot25c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot25c,PerCore~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot25c, c('region','PerCore')) #CUSTOM
temp$region <- factor(datPlot25c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(PerCore ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val25c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
# datPlot25c$parResid <- cal.Par.Resid(datPlot25c,mmod)

#MANUAL #CHANGE
#To make stat brackets fall within y-axis range
df_p_val25b$y.position[2] <- df_p_val25b$y.position[2] + 5
df_p_val25c$y.position[1] <- df_p_val25c$y.position[1] - 10
df_p_val25c$y.position[2] <- df_p_val25c$y.position[2] - 5

# Plotting 
nonFirstyScale <- scale_y_continuous(limits=c(0,100), 
                                     breaks=c(0,25,50,75,100),
                                     oob=rescale_none, 
                                     labels=c("0","25","50","75","100"),
                                     name=element_blank())
plot25b <- ggplot(datPlot25b, aes(region,PerCore))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + pPerCoreYlabels + 
  labs(y = "% Active voxels in lesion core", title = "Rim-Active") #CUSTOM
plot25c <- ggplot(datPlot25c, aes(region,PerCore))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Active") + nonFirstyScale

if (nrow(df_p_val25b)!=0) { plot25b <- plot25b + add_stats(df_p_val25b) }
if (nrow(df_p_val25c)!=0) { plot25c <- plot25c + add_stats(df_p_val25c) }

p25 <- plot25b + plot25c + plot_layout(ncol=2,nrow=1)
ggsave(sprintf("%sFigure25.jpg",outDir),p25,width=pwCol1+pwCol2,height=ph,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot25b, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(PerCore ~ (1 | subject.name) + WML.volume + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot25b, region <- relevel(region, ref = "nws"))
# mmod <- lmer(PerCore ~ (1 | subject.name) + WML.volume + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)

# FIGURE F: Composite Figure ----------------------------------------------
# design <- "AAAAAB
#            CCCCCD
#            FFFFFG"
# figureF <- wrap_elements(full=p19) + p21 + 
#   wrap_elements(full=p20) + p22 + wrap_elements(full=p23) + plot_layout(design=design)
# ggsave(sprintf("%sFigureF.jpg",outDir),figureF,
#        width=0.80*(pwCol1+pwCol2),height=25.92,units="cm",dpi=1200)
design <- "AAAAAB
           CCCCCD
           EEEEEF"
figureF <- wrap_elements(full=p19) + p21 + 
  wrap_elements(full=p20) + p22  + wrap_elements(full=p24) + plot_layout(design=design)
ggsave(sprintf("%sFigureF.jpg",outDir),figureF,
       width=0.80*(pwCol1+pwCol2),height=25.92,units="cm",dpi=1200)

# EMPTY SPACE -------------------------------------------------------------
# REPORTING STATS (BELOW) -------------------------------------------------------------


# Within-region and between-group comparisons: WML MTsat -----------------------------
#WML mtsat (ws)
# tmpDat <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$region=="ws",] #switch between "ws", "mixed" & "nws"
# tmpDat <- Fin.WMLs.mtsatBHs[Fin.WMLs.mtsatBHs$region=="nws",] #switch between "ws", "mixed" & "nws"
tmpDat <- Fin.WMLs.mtsatnonBHs[Fin.WMLs.mtsatnonBHs$region=="nws",] #switch between "ws", "mixed" & "nws"

#A) clinical status 
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + clinical.status,
             tmpDat, REML=TRUE)
summary(mmod)
tmp <- emmeans_test(tmpDat,mtsat~clinical.status,p.adjust.method="fdr",model=mmod)

#LINEAR MODEL STATS FOR REPORTING
#RRMS vs HC, PMS vs HC
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + clinical.status, tmpDat,REML=TRUE) #CUSTOM
coef(mmod)$subject.name["clinical.statusRRMS"][1,]
confint(mmod)["clinical.statusRRMS",]
coef(mmod)$subject.name["clinical.statusPMS"][1,]
confint(mmod)["clinical.statusPMS",]
#PMS vs RRMS
infotmpDat <- within(tmpDat, clinical.status <- relevel(clinical.status, ref = "RRMS"))
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + clinical.status, infotmpDat,REML=TRUE) #CUSTOM
coef(mmod)$subject.name["clinical.statusPMS"][1,]
confint(mmod)["clinical.statusPMS",]
#Criterions to report
REMLcrit(mmod)
#Use ML fitting method to get AIC for model
MLtmpDat <- lmer(mtsat ~ (1 | subject.name) + WML.volume + clinical.status, tmpDat,REML=FALSE) #CUSTOM
AIC(MLtmpDat)

#B) subject type 
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + subj.type,
             tmpDat, REML=TRUE)
summary(mmod)
tmp <- emmeans_test(tmpDat,mtsat~subj.type,p.adjust.method="fdr",model=mmod)

#LINEAR MODEL STATS FOR REPORTING
#MS vs HC
mmod <- lmer(mtsat ~ (1 | subject.name) + WML.volume + subj.type, tmpDat,REML=TRUE) #CUSTOM
coef(mmod)$subject.name["subj.typeMS"][1,]
confint(mmod)["subj.typeMS",]
#Criterions to report
REMLcrit(mmod)
#Use ML fitting method to get AIC for model
MLtmpDat <- lmer(mtsat ~ (1 | subject.name) + WML.volume + subj.type, tmpDat,REML=FALSE) #CUSTOM
AIC(MLtmpDat)

# Within-region and between-group comparisons: WML SUVR -----------------------------
#WML SUVR (ws)
# tmpDat <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$region=="nws",] #switch between "ws", "mixed" & "nws"
# tmpDat <- Fin.WMLs.SUVRBHs[Fin.WMLs.SUVRBHs$region=="nws",] #switch between "ws", "mixed" & "nws"
tmpDat <- Fin.WMLs.SUVRnonBHs[Fin.WMLs.SUVRnonBHs$region=="nws",] #switch between "ws", "mixed" & "nws"

#A) clinical status 
mmod <- lmer(SUVR ~ (1 | subject.name) + WML.volume + clinical.status + bindAff,
             tmpDat, REML=TRUE)
summary(mmod)
tmp <- emmeans_test(tmpDat,SUVR~clinical.status,p.adjust.method="fdr",model=mmod)

#LINEAR MODEL STATS FOR REPORTING
#RRMS vs HC, PMS vs HC
mmod <- lmer(SUVR ~ (1 | subject.name) + WML.volume + clinical.status + bindAff, tmpDat,REML=TRUE) #CUSTOM
coef(mmod)$subject.name["clinical.statusRRMS"][1,]
confint(mmod)["clinical.statusRRMS",]
coef(mmod)$subject.name["clinical.statusPMS"][1,]
confint(mmod)["clinical.statusPMS",]
#PMS vs RRMS
infotmpDat <- within(tmpDat, clinical.status <- relevel(clinical.status, ref = "RRMS"))
mmod <- lmer(SUVR ~ (1 | subject.name) + WML.volume + clinical.status + bindAff, infotmpDat,REML=TRUE) #CUSTOM
coef(mmod)$subject.name["clinical.statusPMS"][1,]
confint(mmod)["clinical.statusPMS",]
#Criterions to report
REMLcrit(mmod)
#Use ML fitting method to get AIC for model
MLtmpDat <- lmer(SUVR ~ (1 | subject.name) + WML.volume + clinical.status + bindAff, tmpDat,REML=FALSE) #CUSTOM
AIC(MLtmpDat)

#A) subject type 
mmod <- lmer(SUVR ~ (1 | subject.name) + WML.volume + subj.type + bindAff,
             tmpDat, REML=TRUE)
summary(mmod)
tmp <- emmeans_test(tmpDat,SUVR~subj.type,p.adjust.method="fdr",model=mmod)

#LINEAR MODEL STATS FOR REPORTING
#MS vs HC
mmod <- lmer(SUVR ~ (1 | subject.name) + WML.volume + subj.type + bindAff, tmpDat,REML=TRUE) #CUSTOM
coef(mmod)$subject.name["subj.typeMS"][1,]
confint(mmod)["subj.typeMS",]
#Criterions to report
REMLcrit(mmod)
#Use ML fitting method to get AIC for model
MLtmpDat <- lmer(SUVR ~ (1 | subject.name) + WML.volume + subj.type + bindAff, tmpDat,REML=FALSE) #CUSTOM
AIC(MLtmpDat)

# Within-region and between-group comparisons: WML Volume -----------------------------
#Smol mtsat (ws)
tmpDat <- WMLs.Vol[WMLs.Vol$Vol<=threshWMLVol & 
                           WMLs.Vol$region=="mixed",] #switch between "ws", "mixed" & "nws"

#A) clinical status 
mmod <- lmer(logVol ~ (1 | subject.name)  + age + clinical.status,
             tmpDat, REML=TRUE)
summary(mmod)
tmp <- emmeans_test(tmpDat,mtsat~clinical.status,p.adjust.method="fdr",model=mmod)

#LINEAR MODEL STATS FOR REPORTING
#PMS vs RRMS
infotmpDat <- within(tmpDat, clinical.status <- relevel(clinical.status, ref = "RRMS"))
mmod <- lmer(logVol ~ (1 | subject.name) + age + clinical.status, infotmpDat,REML=TRUE) #CUSTOM
coef(mmod)$subject.name["clinical.statusPMS"][1,]
confint(mmod)["clinical.statusPMS",]
#Criterions to report
REMLcrit(mmod)
#Use ML fitting method to get AIC for model
MLtmpDat <- lmer(logVol ~ (1 | subject.name) + clinical.status, tmpDat,REML=FALSE) #CUSTOM
AIC(MLtmpDat)


# TSPO WML Phenotype Proportions Fisher's test (By Count) ------------------------------------------
RimActive.ws.count <- sum(RimActive.Count$count[RimActive.Count$region=="ws"])
Active.ws.count <- sum(Active.Count$count[Active.Count$region=="ws"])
Inactive.ws.count <- sum(Inactive.Count$count[Inactive.Count$region=="ws"])

RimActive.mixed.count <- sum(RimActive.Count$count[RimActive.Count$region=="mixed"])
Active.mixed.count <- sum(Active.Count$count[Active.Count$region=="mixed"])
Inactive.mixed.count <- sum(Inactive.Count$count[Inactive.Count$region=="mixed"])

RimActive.nws.count <- sum(RimActive.Count$count[RimActive.Count$region=="nws"])
Active.nws.count <- sum(Active.Count$count[Active.Count$region=="nws"])
Inactive.nws.count <- sum(Inactive.Count$count[Inactive.Count$region=="nws"])

A = matrix(
  c(RimActive.ws.count, Active.ws.count, Inactive.ws.count, 
    RimActive.mixed.count, Active.mixed.count, Inactive.mixed.count, 
    RimActive.nws.count, Active.nws.count, Inactive.nws.count), 
  nrow = 3,   
  ncol = 3,         
  byrow = TRUE          
)

colnames(A) = c("RA","Act","InA")
rownames(A) = c("ws","mixed","nws")

chisq <- chisq.test(A)
res <- fisher.tests.mcp(A)

# Supplementary Figure: WML SUVR (mixed lesions with 50% ws volume omitted) ----------------------
#Ordinary means and standard deviations 
tmp <- data.frame(SUVR=Fin.WMLs.SUVR$SUVR,
                  subj.type=Fin.WMLs.SUVR$subj.type,
                  region=Fin.WMLs.SUVR$region) #CUSTOM
tmp %>% 
  group_by(subj.type,region) %>%
  get_summary_stats(type="mean_sd")

tmp <- data.frame(SUVR=Fin.WMLs.SUVR$SUVR,
                  clinical.status=Fin.WMLs.SUVR$clinical.status,
                  region=Fin.WMLs.SUVR$region) #CUSTOM
tmp %>% 
  group_by(clinical.status,region) %>%
  get_summary_stats(type="mean_sd")

#Plot 5a GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot5a <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$subj.type=="HC",] #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot5a,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot5a,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot5a, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot5a$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val5a <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot5a$parResid <- cal.Par.Resid(datPlot5a,mmod)

#Plot 5b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot5b <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$subj.type=="MS",] #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot5b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot5b,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot5b, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot5b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val5b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot5b$parResid <- cal.Par.Resid(datPlot5b,mmod)

#Plot 5c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot5c <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$clinical.status=="RRMS",] #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot5c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot5c,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot5c, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot5c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val5c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot5c$parResid <- cal.Par.Resid(datPlot5c,mmod)

#Plot 5d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot5d <- Fin.WMLs.SUVR[Fin.WMLs.SUVR$clinical.status=="PMS",] #CUSTOM
mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot5d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot5d,SUVR~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot5d, c('region','SUVR')) #CUSTOM
temp$region <- factor(datPlot5d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(SUVR ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val5d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot5d$parResid <- cal.Par.Resid(datPlot5d,mmod)

# Plotting 
nonFirstyScale <- scale_y_continuous(limits=c(0.25,2.3), 
                                     breaks=c(0.5,1.0,1.5,2.0),
                                     oob=rescale_none, 
                                     labels=c("","","",""),
                                     name=element_blank())
plot5a <- ggplot(datPlot5a, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + pSUVRYlabels + 
  labs(y = "[11C]PBR28 PET SUVR") #CUSTOM
plot5b <- ggplot(datPlot5b, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale 
plot5c <- ggplot(datPlot5c, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale 
plot5d <- ggplot(datPlot5d, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale 

if (nrow(df_p_val5a)!=0) { plot5a <- plot5a + add_stats(df_p_val5a) }
if (nrow(df_p_val5b)!=0) { plot5b <- plot5b + add_stats(df_p_val5b) }
if (nrow(df_p_val5c)!=0) { plot5c <- plot5c + add_stats(df_p_val5c) }
if (nrow(df_p_val5d)!=0) { plot5d <- plot5d + add_stats(df_p_val5d) }

p5 <- plot5a + plot5b + plot5c + plot5d + plot_layout(ncol=4,nrow=1)
ggsave(sprintf("%sSuppFigure_Mixed50prctOmit_TSPOWML_SUVR.jpg",outDir),p5,width=pwCol1+pwCol2,height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot5d, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot5d, region <- relevel(region, ref = "nws"))
# mmod <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(SUVR ~ (1 | subject.name) + region + WML.volume + bindAff, datPlot5d,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# Supplementary Figure: Composite Figure ----------------------------------------------
# figureB <- wrap_elements(full=p4)+wrap_elements(full=p5)+
#   wrap_elements(full=p6)+wrap_elements(full=p7)+plot_layout(nrow=4,ncol=1)
# ggsave(sprintf("%sFigureB.jpg",outDir),figureB,
#        width=pwCol1+pwCol2,height=34.56,units="cm",dpi=1200)
figureB <- wrap_elements(full=p4)+wrap_elements(full=p5)+plot_layout(nrow=2,ncol=1)
ggsave(sprintf("%sFigureB.jpg",outDir),figureB,
       width=pwCol1+pwCol2,height=17.28,units="cm",dpi=1200)

#Try only area between 3rd and 4th quaritle of each violin plot a transparent white! 

# SUPPLEMENTARY FIGURES -------------------------------------------------------------


# Figure 19: TSPO WML Phenotype Proportions (By Count), Violin Plots ------------------------------------
#Switch between "RimActive", "Active" and "Inactive"
#Ordinary means and standard deviations 
tmp <- data.frame(Prop=Active.Count$Prop,
                  region=Active.Count$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  # get_summary_stats(type="median_iqr")
  get_summary_stats(type="mean_sd")

#Plot 19b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot19b <- RimActive.Count #CUSTOM
mmod <- lmer(Prop ~ (1 | subject.name) + region + age, datPlot19b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot19b,Prop~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot19b, c('region','Prop')) #CUSTOM
temp$region <- factor(datPlot19b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(Prop ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val19b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 19c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot19c <- Active.Count #CUSTOM
mmod <- lmer(Prop ~ (1 | subject.name) + region + age, datPlot19c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot19c,Prop~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot19c, c('region','Prop')) #CUSTOM
temp$region <- factor(datPlot19c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(Prop ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val19c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 19d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot19d <- Inactive.Count #CUSTOM
mmod <- lmer(Prop ~ (1 | subject.name) + region + age, datPlot19d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot19d,Prop~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot19d, c('region','Prop')) #CUSTOM
temp$region <- factor(datPlot19d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(Prop ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val19d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#MANUAL
df_p_val19c$y.position <- df_p_val19d$y.position 
df_p_val19c$y.position[1] <- df_p_val19c$y.position[1] + 0.5
df_p_val19d$y.position[1] <- df_p_val19d$y.position[1] + 0.5
df_p_val19c$y.position[2] <- df_p_val19c$y.position[1] + 5
df_p_val19d$y.position[2] <- df_p_val19d$y.position[1] + 5

#Plotting
nonFirstyScale <- scale_y_continuous(limits=c(0,110), 
                                     breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                                     oob=rescale_none, 
                                     labels=c("","","","","","","","","","",""),
                                     name=element_blank())
plot19b <- ggplot(datPlot19b, aes(region,Prop))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels +
  labs(y = 'Proportion of WM lesions (%)', title = "Rim-Active") + #CUSTOM
  scale_y_continuous(limits=c(0,110), 
                     breaks=c(0,10,20,30,40,50,60,70,80,90,100))
plot19c <- ggplot(datPlot19c, aes(region,Prop))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Active") + nonFirstyScale
plot19d <- ggplot(datPlot19d, aes(region,Prop))  +
  pPoint + pViolinBack + pViolinFront +
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Inactive") + nonFirstyScale

if (nrow(df_p_val19b)!=0) { plot19b <- plot19b + add_stats(df_p_val19b) }
if (nrow(df_p_val19c)!=0) { plot19c <- plot19c + add_stats(df_p_val19c) }
if (nrow(df_p_val19d)!=0) { plot19d <- plot19d + add_stats(df_p_val19d) }

p19 <- plot19b + plot19c + plot19d + plot_layout(ncol=3)
# ggsave(sprintf("%sFigure19.jpg",outDir),p19,width=0.75*(pwCol1+pwCol2),height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot19d, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(Prop ~ (1 | subject.name) + region + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot19d, region <- relevel(region, ref = "nws"))
# mmod <- lmer(Prop ~ (1 | subject.name) + region + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(Prop ~ (1 | subject.name) + region + age, datPlot19d,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# Figure 20: TSPO WML Phenotype Volumes -----------------------------------
#Switch between "RimActive", "Active" and "Inactive"
#Ordinary means and standard deviations 
tmp <- data.frame(Vol=Inactive.Vol$logVol,
                  region=Inactive.Vol$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  get_summary_stats(type="mean_sd")

#Plot 20b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot20b <- RimActive.Vol #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot20b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot20b,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot20b, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot20b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val20b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 20c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot20c <- Active.Vol #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot20c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot20c,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot20c, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot20c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val20c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 20d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot20d <- Inactive.Vol #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot20d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot20d,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot20d, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot20d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val20d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#MANUAL
df_p_val20b$y.position[1] <- df_p_val20b$y.position[1] + 0.1
df_p_val20b$y.position[2] <- df_p_val20b$y.position[1] + 0.5
df_p_val20b$y.position[3] <- df_p_val20b$y.position[2] + 0.5
df_p_val20c$y.position[1] <- df_p_val20c$y.position[1] + 0.1
df_p_val20c$y.position[2] <- df_p_val20c$y.position[1] + 0.5
df_p_val20c$y.position[3] <- df_p_val20c$y.position[2] + 0.5

# Plotting 
nonFirstyScale <- scale_y_continuous(limits=c(-3.75,4.5), 
                                     breaks=c(-3,-2,-1,0,1,2,3,4),
                                     oob=rescale_none, 
                                     labels=c("","","","","","","",""),
                                     name=element_blank())
plot20b <- ggplot(datPlot20b, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels +
  labs(y = ~paste(bold("Log-transformed lesion volume ("),bold("cm"^"3"),bold(") "))) + pVolumeYlabels
plot20c <- ggplot(datPlot20c, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "") + nonFirstyScale
plot20d <- ggplot(datPlot20d, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "") + nonFirstyScale

if (nrow(df_p_val20b)!=0) { plot20b <- plot20b + add_stats(df_p_val20b) } 
if (nrow(df_p_val20c)!=0) { plot20c <- plot20c + add_stats(df_p_val20c) } 
if (nrow(df_p_val20d)!=0) { plot20d <- plot20d + add_stats(df_p_val20d) }

p20 <- plot20b + plot20c + plot20d + plot_layout(ncol=3)
# ggsave(sprintf("%sFigure20.jpg",outDir),p20,width=0.75*(pwCol1+pwCol2),height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot20c, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(logVol ~ (1 | subject.name) + region + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot20c, region <- relevel(region, ref = "nws"))
# mmod <- lmer(logVol ~ (1 | subject.name) + region + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(logVol ~ (1 | subject.name) + region + age, datPlot20c,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# Figure 21: TSPO WML Phenotype Proportions (By Count), Stacked Bars -------------------------------
RimActive.ws.count <- sum(RimActive.Count$count[RimActive.Count$region=="ws"])
Active.ws.count <- sum(Active.Count$count[Active.Count$region=="ws"])
Inactive.ws.count <- sum(Inactive.Count$count[Inactive.Count$region=="ws"])

RimActive.mixed.count <- sum(RimActive.Count$count[RimActive.Count$region=="mixed"])
Active.mixed.count <- sum(Active.Count$count[Active.Count$region=="mixed"])
Inactive.mixed.count <- sum(Inactive.Count$count[Inactive.Count$region=="mixed"])

RimActive.nws.count <- sum(RimActive.Count$count[RimActive.Count$region=="nws"])
Active.nws.count <- sum(Active.Count$count[Active.Count$region=="nws"])
Inactive.nws.count <- sum(Inactive.Count$count[Inactive.Count$region=="nws"])

#Not using original WML counts, since for TSPO phenotyping, lesions with rim volume < 27mm3 were discounted 
ws.count <- RimActive.ws.count + Active.ws.count + Inactive.ws.count 
mixed.count <- RimActive.mixed.count + Active.mixed.count + Inactive.mixed.count 
nws.count <- RimActive.nws.count + Active.nws.count + Inactive.nws.count 

prop1 <- c(RimActive.ws.count/ws.count, Active.ws.count/ws.count, Inactive.ws.count/ws.count, 
           RimActive.mixed.count/mixed.count, Active.mixed.count/mixed.count, Inactive.mixed.count/mixed.count, 
           RimActive.nws.count/nws.count, Active.nws.count/nws.count, Inactive.nws.count/nws.count)*100 
df1 <- data.frame(region=factor(c("ws", "ws", "ws", "mixed", "mixed", "mixed", "nws", "nws", "nws"),levels=c("ws", "mixed", "nws")),
                  TSPOphenotype=factor(c("RA", "A", "IA", "RA", "A", "IA", "RA", "A", "IA"),levels=c("IA", "A", "RA")),
                  prop=prop1)
p21 <- ggplot(df1, aes(x=region,y=prop1,fill=TSPOphenotype)) + 
  geom_bar(stat="identity",position="stack",color="black",linewidth=0.2,width=0.7) + 
  scale_fill_manual(values=colors.TSPOPhenotype) + 
  labs(y="Proportion of total lesion count (%)",title="") + 
  pProplabels + pthemeProp

# Figure 22: TSPO WML Phenotype Proportions (By Total Volume), Stacked Bars -------------------------------
RimActive.ws.Vol <- sum(RimActive.Vol$Vol[RimActive.Vol$region=="ws"])
Active.ws.Vol <- sum(Active.Vol$Vol[Active.Vol$region=="ws"])
Inactive.ws.Vol <- sum(Inactive.Vol$Vol[Inactive.Vol$region=="ws"])

RimActive.mixed.Vol <- sum(RimActive.Vol$Vol[RimActive.Vol$region=="mixed"])
Active.mixed.Vol <- sum(Active.Vol$Vol[Active.Vol$region=="mixed"])
Inactive.mixed.Vol <- sum(Inactive.Vol$Vol[Inactive.Vol$region=="mixed"])

RimActive.nws.Vol <- sum(RimActive.Vol$Vol[RimActive.Vol$region=="nws"])
Active.nws.Vol <- sum(Active.Vol$Vol[Active.Vol$region=="nws"])
Inactive.nws.Vol <- sum(Inactive.Vol$Vol[Inactive.Vol$region=="nws"])

#Not using original WML volumes, since for TSPO phenotyping, lesions with rim volume < 27mm3 were discounted 
ws.Vol <- RimActive.ws.Vol + Active.ws.Vol + Inactive.ws.Vol 
mixed.Vol <- RimActive.mixed.Vol + Active.mixed.Vol + Inactive.mixed.Vol 
nws.Vol <- RimActive.nws.Vol + Active.nws.Vol + Inactive.nws.Vol 

prop2 <- c(RimActive.ws.Vol/ws.Vol, Active.ws.Vol/ws.Vol, Inactive.ws.Vol/ws.Vol, 
           RimActive.mixed.Vol/mixed.Vol, Active.mixed.Vol/mixed.Vol, Inactive.mixed.Vol/mixed.Vol, 
           RimActive.nws.Vol/nws.Vol, Active.nws.Vol/nws.Vol, Inactive.nws.Vol/nws.Vol)*100 
df2 <- data.frame(region=factor(c("ws", "ws", "ws", "mixed", "mixed", "mixed", "nws", "nws", "nws"),levels=c("ws", "mixed", "nws")),
                  TSPOphenotype=factor(c("RA", "A", "IA", "RA", "A", "IA", "RA", "A", "IA"),levels=c("IA", "A", "RA")),
                  prop=prop2)
p22 <- ggplot(df2, aes(x=region,y=prop2,fill=TSPOphenotype)) + 
  geom_bar(stat="identity",position="stack",color="black",linewidth=0.2,width=0.7) + 
  scale_fill_manual(values=colors.TSPOPhenotype) + 
  labs(y="Proportion of total lesion volume (%)",title="") + 
  pProplabels + pthemeProp 

# Figure 24: TSPO WML Phenotype MTsat (per mm3) -----------------------------------
#Switch between "RimActive", "Active" and "Inactive"
#Ordinary means and standard deviations 
tmp <- data.frame(mtsatpervol=Inactive.mtsat$mtsatpervol,
                  region=Inactive.mtsat$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  get_summary_stats(type="mean_sd")

#Plot 24b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot24b <- RimActive.mtsat #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot24b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot24b,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot24b, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot24b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val24b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot24b$parResid <- cal.Par.Resid(datPlot24b,mmod)

#Plot 24c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot24c <- Active.mtsat #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot24c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot24c,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot24c, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot24c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val24c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot24c$parResid <- cal.Par.Resid(datPlot24c,mmod)

#Plot 24d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot24d <- Inactive.mtsat #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot24d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot24d,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot24d, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot24d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val24d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot24d$parResid <- cal.Par.Resid(datPlot24d,mmod)

#To fit stat brackets within y-scale 
df_p_val24b$y.position[1] <- df_p_val24b$y.position[1] - 0.003
df_p_val24b$y.position[2] <- df_p_val24b$y.position[1] + 0.004
df_p_val24b$y.position[3] <- df_p_val24b$y.position[2] + 0.004
df_p_val24c$y.position[2] <- df_p_val24c$y.position[1] + 0.004
df_p_val24c$y.position[3] <- df_p_val24c$y.position[2] + 0.004

# Plotting 
pMTsatYlabelsv2 <- scale_y_continuous(limits=c(0,0.060), 
                                      breaks=c(0,0.015,0.030,0.045,0.060),
                                      oob=rescale_none, 
                                      labels=c("0","0.015","0.030","0.045","0.060"))
nonFirstyScale <- scale_y_continuous(limits=c(0,0.060), 
                                     breaks=c(0,0.015,0.030,0.045,0.060),
                                     oob=rescale_none, 
                                     labels=c("","","","",""),
                                     name=element_blank())
plot24b <- ggplot(datPlot24b, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + pMTsatYlabelsv2 + 
  labs(y=expression(bold(paste('MTsat per unit volume (p.u./',mm^'3',')')))) #CUSTOM
plot24c <- ggplot(datPlot24c, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale
plot24d <- ggplot(datPlot24d, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale

if (nrow(df_p_val24b)!=0) { plot24b <- plot24b + add_stats(df_p_val24b) }
if (nrow(df_p_val24c)!=0) { plot24c <- plot24c + add_stats(df_p_val24c) }
if (nrow(df_p_val24d)!=0) { plot24d <- plot24d + add_stats(df_p_val24d) }

p24 <- plot24b + plot24c + plot24d + plot_layout(ncol=3,nrow=1)
ggsave(sprintf("%sSuppFigure_TSPOWML_MTsatpermm3.jpg",outDir),p24,width=pwCol1+pwCol2,height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot24d, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot24d, region <- relevel(region, ref = "nws"))
# mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot24d,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# FIGURE F: Composite Figure ----------------------------------------------
# design <- "AAAAAB
#            CCCCCD
#            FFFFFG"
# figureF <- wrap_elements(full=p19) + p21 + 
#   wrap_elements(full=p20) + p22 + wrap_elements(full=p23) + plot_layout(design=design)
# ggsave(sprintf("%sFigureF.jpg",outDir),figureF,
#        width=0.80*(pwCol1+pwCol2),height=25.92,units="cm",dpi=1200)
design <- "AAAAAB
           CCCCCD
           EEEEEF"
figureF <- wrap_elements(full=p19) + p21 + 
  wrap_elements(full=p20) + p22  + wrap_elements(full=p24) + plot_layout(design=design)
ggsave(sprintf("%sSupplementaryFigure3.jpg",outDir),figureF,
       width=0.80*(pwCol1+pwCol2),height=25.92,units="cm",dpi=1200)

# Figure 1: WM lesion count ----------------------
#Ordinary means and standard deviations 
tmp <- data.frame(count=WMLs.Count$count,
                  region=WMLs.Count$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  get_summary_stats(type="mean_sd")

tmp <- data.frame(count=WMLs.Count$count,
                  clinical.status=WMLs.Count$clinical.status,
                  region=WMLs.Count$region) #CUSTOM
tmp %>% 
  group_by(clinical.status,region) %>%
  get_summary_stats(type="mean_sd")

#Plot 1b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot1b <- WMLs.Count #CUSTOM
mmod <- lmer(count ~ (1 | subject.name) + region, datPlot1b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot1b,count~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot1b, c('region','count')) #CUSTOM
temp$region <- factor(datPlot1b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(count ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val1b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 1c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot1c <- WMLs.Count[WMLs.Count$clinical.status=="RRMS",] #CUSTOM
mmod <- lmer(count ~ (1 | subject.name) + region, datPlot1c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot1c,count~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot1c, c('region','count')) #CUSTOM
temp$region <- factor(datPlot1c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(count ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val1c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 1d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot1d <- WMLs.Count[WMLs.Count$clinical.status=="PMS",] #CUSTOM
mmod <- lmer(count ~ (1 | subject.name) + region, datPlot1d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot1d,count~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot1d, c('region','count')) #CUSTOM
temp$region <- factor(datPlot1d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(count ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val1d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#MANUAL
# df_p_val1b$y.position <- df_p_val1b$y.position-10
# df_p_val1d$y.position <- df_p_val1d$y.position-10

#Plotting
nonFirstyScale <- scale_y_continuous(limits=c(0,85), 
                                     breaks=c(0,10,20,30,40,50,60,70,80),
                                     oob=rescale_none, 
                                     labels=c("","","","","","","","",""),
                                     name=element_blank())
plot1b <- ggplot(datPlot1b, aes(region,count))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels +
  labs(y = 'Lesion count', title = "MS") + #CUSTOM
  scale_y_continuous(limits=c(0,85), 
                     breaks=c(0,10,20,30,40,50,60,70,80))
plot1c <- ggplot(datPlot1c, aes(region,count))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "RRMS") + nonFirstyScale
plot1d <- ggplot(datPlot1d, aes(region,count))  +
  pPoint + pViolinBack + pViolinFront +
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title = "Progressive MS") + nonFirstyScale

if (nrow(df_p_val1b)!=0) { plot1b <- plot1b + add_stats(df_p_val1b) }
if (nrow(df_p_val1c)!=0) { plot1c <- plot1c + add_stats(df_p_val1c) }
if (nrow(df_p_val1d)!=0) { plot1d <- plot1d + add_stats(df_p_val1d) }

p1 <- plot1b + plot1c + plot1d + plot_layout(ncol=3)
# ggsave(sprintf("%sFigure1.jpg",outDir),p1,width=0.75*(pwCol1+pwCol2),height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot1d, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(count ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot1d, region <- relevel(region, ref = "nws"))
# mmod <- lmer(count ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(count ~ (1 | subject.name) + region, datPlot1d,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# Figure 2: WML Volume ----------------------
#Ordinary means and standard deviations 
tmp <- data.frame(Vol=WMLs.Vol$logVol,
                  region=WMLs.Vol$region) #CUSTOM
tmp %>%
  group_by(region) %>%
  get_summary_stats(type="mean_sd")

tmp <- data.frame(Vol=WMLs.Vol$logVol,
                  clinical.status=WMLs.Vol$clinical.status,
                  region=WMLs.Vol$region) #CUSTOM
tmp %>% 
  group_by(clinical.status,region) %>%
  get_summary_stats(type="mean_sd")

#Plot 2b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot2b <- WMLs.Vol #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot2b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot2b,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot2b, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot2b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val2b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 2c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot2c <- WMLs.Vol[WMLs.Vol$clinical.status=="RRMS",] #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot2c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot2c,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot2c, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot2c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val2c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#Plot 2d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot2d <- WMLs.Vol[WMLs.Vol$clinical.status=="PMS",] #CUSTOM
mmod <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot2d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot2d,logVol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot2d, c('region','logVol')) #CUSTOM
temp$region <- factor(datPlot2d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(logVol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val2d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)

#MANUAL
df_p_val2b$y.position[1] <- df_p_val2b$y.position[1] + 0.2
df_p_val2b$y.position[2] <- df_p_val2b$y.position[1] + 0.4
df_p_val2b$y.position[3] <- df_p_val2b$y.position[2] + 0.4
df_p_val2c$y.position[3] <- df_p_val2c$y.position[2] + 0.4
df_p_val2d$y.position[1] <- df_p_val2d$y.position[1] + 0.2
df_p_val2d$y.position[2] <- df_p_val2d$y.position[1] + 0.4
df_p_val2d$y.position[3] <- df_p_val2d$y.position[2] + 0.4

# Plotting 
nonFirstyScale <- scale_y_continuous(limits=c(-3.75,4.5), 
                                     breaks=c(-3,-2,-1,0,1,2,3,4),
                                     oob=rescale_none, 
                                     labels=c("","","","","","","",""),
                                     name=element_blank())
plot2b <- ggplot(datPlot2b, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels +  
  labs(y = ~paste(bold("Log-transformed lesion volume ("),bold("cm"^"3"),bold(") "))) + pVolumeYlabels 
plot2c <- ggplot(datPlot2c, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale
plot2d <- ggplot(datPlot2d, aes(region,logVol))  +
  pPoint + pViolinBack + pViolinFront +  
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "") + nonFirstyScale

if (nrow(df_p_val2b)!=0) { plot2b <- plot2b + add_stats(df_p_val2b) } 
if (nrow(df_p_val2c)!=0) { plot2c <- plot2c + add_stats(df_p_val2c) } 
if (nrow(df_p_val2d)!=0) { plot2d <- plot2d + add_stats(df_p_val2d) } 

p2 <- plot2b + plot2c + plot2d + plot_layout(ncol=3)
# ggsave(sprintf("%sFigure2.jpg",outDir),p2,width=0.75*(pwCol1+pwCol2),height=ph/2,units="cm",dpi=1200)

# #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot2d, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(logVol ~ region + (1 | subject.name) + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot2d, region <- relevel(region, ref = "nws"))
# mmod <- lmer(logVol ~ region + (1 | subject.name) + age, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)["regionws",]
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)["regionmixed",]
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(logVol ~ region + (1 | subject.name) + age, datPlot2d,REML=FALSE) #CUSTOM
# AIC(MLtmp)

# Figure 3: WM lesion count per patient (stratified by vascular zone) -----------------------------------
p2b <- ggplot(WMLs.Count, aes(x=reorder(subject.name,count),y=count,fill=region)) + 
  geom_bar(stat="identity",position="stack",linewidth=0.4, width=0.65, color="black") + 
  pScaleFill + pthemepatientwise + pCountlabels + labs(y = 'WM lesion count') + 
  facet_wrap(~clinical.status,nrow=1,scales="free_x", strip.position = "bottom") 

# FIGURE A: Composite Figure ----------------------------------------------
figureA <- wrap_elements(full=p1)+wrap_elements(full=p2)+
  wrap_elements(p2b)+plot_layout(nrow=3,ncol=1)
ggsave(sprintf("%sSupplementaryFigure3.jpg",outDir),figureA,
       width=0.75*(pwCol1+pwCol2),height=26.5,units="cm",dpi=1200)

# Figure 6: WML MTsat (per mm3) ----------------------
#Ordinary means and standard deviations 
tmp <- data.frame(mtsatpervol=Fin.WMLs.mtsat$mtsatpervol,
                  subj.type=Fin.WMLs.mtsat$subj.type,
                  region=Fin.WMLs.mtsat$region) #CUSTOM
tmp %>% 
  group_by(subj.type,region) %>%
  get_summary_stats(type="mean_sd")

tmp <- data.frame(mtsatpervol=Fin.WMLs.mtsat$mtsatpervol,
                  clinical.status=Fin.WMLs.mtsat$clinical.status,
                  region=Fin.WMLs.mtsat$region) #CUSTOM
tmp %>% 
  group_by(clinical.status,region) %>%
  get_summary_stats(type="mean_sd")

#Plot 6b GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot6b <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$subj.type=="MS",] #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot6b,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot6b,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot6b, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot6b$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val6b <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot6b$parResid <- cal.Par.Resid(datPlot6b,mmod)

#Plot 6c GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot6c <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$clinical.status=="RRMS",] #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot6c,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot6c,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot6c, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot6c$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val6c <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot6c$parResid <- cal.Par.Resid(datPlot6c,mmod)

#Plot 6d GLM Stats
#FDR mcp, for pairwise comparisons
df <- data.frame()
datPlot6d <- Fin.WMLs.mtsat[Fin.WMLs.mtsat$clinical.status=="PMS",] #CUSTOM
mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot6d,REML=TRUE) #CUSTOM
summary(mmod)
tmp <- emmeans_test(datPlot6d,mtsatpervol~region,p.adjust.method="fdr",model=mmod) #CUSTOM
df <- rbind(df,tmp)
df$p.adj.signif[df$p.adj < 0.10 & df$p.adj > 0.05] <- "#"
#Plotting Stat Brackets
temp <- dplyr::select(datPlot6d, c('region','mtsatpervol')) #CUSTOM
temp$region <- factor(datPlot6d$region, levels=c("ws", "mixed", "nws"))
#Need to set t-test condition for detecting "not enough y observations"
df.form <- temp %>%
  rstatix::t_test(mtsatpervol ~ region) %>% #CUSTOM
  rstatix::add_xy_position(x = "region", dodge = 0.8, step.increase = pStepIncrease) # important for positioning!
repCols <- c("statistic", "df", "p", "p.adj", "p.adj.signif")
df.form[,repCols] <- df[,repCols]
df_p_val6d <- df.form[df.form$p.adj.signif != "ns",] #CUSTOM
rm(df)
datPlot6d$parResid <- cal.Par.Resid(datPlot6d,mmod)

#To fit stat brackets within y-scale 
df_p_val6b$y.position[2] <- df_p_val6b$y.position[1] + 0.003
df_p_val6b$y.position[3] <- df_p_val6b$y.position[2] + 0.004
df_p_val6c$y.position[2] <- df_p_val6c$y.position[1] + 0.003
df_p_val6c$y.position[3] <- df_p_val6c$y.position[2] + 0.004
df_p_val6d$y.position[2] <- df_p_val6d$y.position[1] + 0.003
df_p_val6d$y.position[3] <- df_p_val6d$y.position[2] + 0.004

# Plotting 
pMTsatYlabelsv2 <- scale_y_continuous(limits=c(0,0.060), 
                                      breaks=c(0,0.015,0.030,0.045,0.060),
                                      oob=rescale_none, 
                                      labels=c("0","0.015","0.030","0.045","0.060"))
nonFirstyScale <- scale_y_continuous(limits=c(0,0.060), 
                                     breaks=c(0,0.015,0.030,0.045,0.060),
                                     oob=rescale_none, 
                                     labels=c("","","","",""),
                                     name=element_blank())
plot6b <- ggplot(datPlot6b, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill  + ptheme +
  pXlabels + 
  labs(y=expression(bold(paste('MTsat per unit volume (p.u./',mm^'3',')'))), 
       title="MS") + pMTsatYlabelsv2
plot6c <- ggplot(datPlot6c, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title="RRMS") + nonFirstyScale
plot6d <- ggplot(datPlot6d, aes(region,parResid))  +
  pPoint + pViolinBack + pViolinFront + 
  pScaleColor + pScaleFill + ptheme +
  pXlabels + labs(y = "", title="Progressive MS") + nonFirstyScale

if (nrow(df_p_val6b)!=0) { plot6b <- plot6b + add_stats(df_p_val6b) }
if (nrow(df_p_val6c)!=0) { plot6c <- plot6c + add_stats(df_p_val6c) }
if (nrow(df_p_val6d)!=0) { plot6d <- plot6d + add_stats(df_p_val6d) }

p6 <- plot6b + plot6c + plot6d + plot_layout(ncol=3,nrow=1)
ggsave(sprintf("%sSupplementaryFigure2.jpg",outDir),p6,width=20.0025,height=ph/2,units="cm",dpi=1200)

# # #LINEAR MODEL STATS FOR REPORTING
# #Watershed versus mixed
# infotmp <- within(datPlot6b, region <- relevel(region, ref = "mixed"))
# mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# #Watershed versus non-watershed, Mixed versus non-watershed
# infotmp <- within(datPlot6b, region <- relevel(region, ref = "nws"))
# mmod <- lmer(mtsatpervol ~ (1 | subject.name) + region, infotmp,REML=TRUE) #CUSTOM
# coef(mmod)$subject.name["regionws"][1,]
# confint(mmod)
# coef(mmod)$subject.name["regionmixed"][1,]
# confint(mmod)
# #Criterions to report
# REMLcrit(mmod)
# #Use ML fitting method to get AIC for model
# MLtmp <- lmer(mtsatpervol ~ (1 | subject.name) + region, datPlot6b,REML=FALSE) #CUSTOM
# AIC(MLtmp)




# Figure 7: WM lesion count per patient (stratified by TSPO PET WML phenotype) -----------------------------------
tmp1 <- RimActive.Count %>% #For RA lesions
  dplyr::group_by(subject.name, clinical.status) %>% #Regarding of vascular zone of lesion (ws,mixed,nws)
  dplyr::summarize(total_count = sum(count, na.rm = TRUE))
tmp2 <- Active.Count %>% #For Active lesions
  dplyr::group_by(subject.name, clinical.status) %>% #Regarding of vascular zone of lesion (ws,mixed,nws)
  dplyr::summarize(total_count = sum(count, na.rm = TRUE))
tmp3 <- Inactive.Count %>% #For Inactive lesions
  dplyr::group_by(subject.name, clinical.status) %>% #Regarding of vascular zone of lesion (ws,mixed,nws)
  dplyr::summarize(total_count = sum(count, na.rm = TRUE))
WMLs.Count.TSPOPET = stack(list(RA=tmp1$total_count,A=tmp2$total_count,IA=tmp3$total_count))
WMLs.Count.TSPOPET = cbind(WMLs.Count.TSPOPET,tmp1$subject.name,tmp1$clinical.status)
names(WMLs.Count.TSPOPET) <- c("count","TSPOPhenotype","subject.name","clinical.status")

p7 <- ggplot(WMLs.Count.TSPOPET, aes(x=reorder(subject.name,count),y=count,fill=TSPOPhenotype)) + 
  geom_bar(stat="identity",position="stack",linewidth=0.4, width=0.65, color="black") + 
  scale_fill_manual(values=colors.TSPOPhenotype,labels=c("RA"="Rim-Active","A"="Active","IA"="Inactive")) + 
  pthemepatientwise + pCountlabels + labs(y = 'WM lesion count') + 
  facet_wrap(~clinical.status,nrow=1,scales="free_x", strip.position = "bottom") 
ggsave(sprintf("%sResponseFigure1.jpg",outDir),p7,
       width=0.80*(pwCol1+pwCol2),height=25.92/2,units="cm",dpi=1200)

# Figure 8: WM lesion volume per patient (stratified by TSPO PET WML phenotype) -----------------------------------
tmp1 <- RimActive.Vol %>% #For RA lesions
  dplyr::group_by(subject.name, clinical.status) %>% #Regarding of vascular zone of lesion (ws,mixed,nws)
  dplyr::summarize(avg.vol = mean(Vol, na.rm = TRUE))
tmp2 <- Active.Vol %>% #For Active lesions
  dplyr::group_by(subject.name, clinical.status) %>% #Regarding of vascular zone of lesion (ws,mixed,nws)
  dplyr::summarize(avg.vol = mean(Vol, na.rm = TRUE))
tmp3 <- Inactive.Vol %>% #For Inactive lesions
  dplyr::group_by(subject.name, clinical.status) %>% #Regarding of vascular zone of lesion (ws,mixed,nws)
  dplyr::summarize(avg.vol = mean(Vol, na.rm = TRUE))
WMLs.Vol.TSPOPET = stack(list(RA=tmp1$avg.vol,A=tmp2$avg.vol,IA=tmp3$avg.vol))
WMLs.Vol.TSPOPET = cbind(WMLs.Vol.TSPOPET,
                         c(tmp1$subject.name,tmp2$subject.name,tmp3$subject.name),
                         c(tmp1$clinical.status,tmp2$clinical.status,tmp3$clinical.status))
names(WMLs.Vol.TSPOPET) <- c("avg.vol","TSPOPhenotype","subject.name","clinical.status")

# Example data where some labels might be missing certain sub-categories
WMLs.Vol.TSPOPET <- WMLs.Vol.TSPOPET %>%
  complete(nesting(subject.name, clinical.status), TSPOPhenotype, fill = list(avg.vol = 0))

p8 <- ggplot(WMLs.Vol.TSPOPET, aes(x=reorder(subject.name,avg.vol),y=avg.vol,fill=TSPOPhenotype)) + 
  geom_bar(stat="identity",position="dodge",linewidth=0.4, width=0.65, color="black") + 
  scale_fill_manual(values=colors.TSPOPhenotype,labels=c("RA"="Rim-Active","A"="Active","IA"="Inactive")) + 
  pthemepatientwise + 
  # scale_y_continuous(limits=c(0,1.5), 
  #                    breaks=c(0,0.5,1,1.5),
  #                    oob=rescale_none, 
  #                    labels=c("0","0.5","1.0","1.5")) + 
  labs(y = 'Average lesion volume (cm3)') + 
  facet_wrap(~clinical.status,nrow=1,scales="free_x", strip.position = "bottom") 
#How to set a "fixed number of" stacked bars per individual
ggsave(sprintf("%sResponseFigure2.jpg",outDir),p8,
       width=0.80*(pwCol1+pwCol2),height=25.92/2,units="cm",dpi=1200)



