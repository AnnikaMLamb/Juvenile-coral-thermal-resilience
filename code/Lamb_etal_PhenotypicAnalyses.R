##R script for analyses of hybrid and purebred heat wave performance
#Data obtained using methods described in manuscript -  Purebred and hybrid coral juveniles demonstrate thermal resilience - methods
#Written by Annika Lamb 

#Load packages and functions
library(lme4)
library(nlme)
library(ggplot2)
library(multcomp)
library(glmmTMB)
library(sjmisc)
library(ggrepel)
library(lsmeans)


#Summary function

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  require(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  
  length2 <- function (x, na.rm=FALSE) {
    
    if (na.rm) sum(!is.na(x))
    
    else       length(x)
    
  }
  

  # This does the summary. For each group's data frame, return a vector with
  
  # N, mean, and sd
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 
                 .fun = function(xx, col) {
                   
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                     
                   )
                   
                 },
                 
                 measurevar
                 
  )
  # Rename the "mean" column    
  
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  
  # Calculate t-statistic for confidence interval: 
  
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
#General graphics
cols<- c("#fc8d62","#66c2a5","#8da0cb")
TP_names <- as_labeller(c('Ambient' = "Ambient", 'Elevated'="Elevated",`30` = "T30", `65` = "T65"))
axes_marks_settings <- element_text(color = "black", size =10)

###SURVIVORSHIP###
#Load data 
Survival <- read.table('HeatWaveExperiment_Results_Singles_Survival.csv', header=T, sep=',')
str(Survival)

#exclude NAs
Survival<-na.omit(Survival)

##Summary, averaged by tile 
summarySurvival_Tile<-summarySE(Survival, measurevar="Survival", groupvars=c("Day","Cross","Tank","TileID","Treatment","DHW"),conf.interval=0.95)
summarySurvival_Tile

#Averaged by tile and then plotted
SummarySurvival_Tile_Crunched<-summarySE(summarySurvival_Tile, measurevar="Survival", groupvars=c("Day","Cross","Treatment","DHW"),conf.interval=0.95)
SummarySurvival_Tile_Crunched

#Plot by treatment and timepoint
survival_graph=ggplot(summarySurvival_Tile, aes(y=Survival, x=Cross, fill=Cross))+
  facet_grid(Day~Treatment, labeller=TP_names)+
  geom_boxplot() + 
  scale_fill_manual(values=cols, labels = c(expression(italic("A. loripes")), "KL hybrid",expression(italic("A. kenti"))))+
  theme_classic()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(size = 0), axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),)+
  theme(legend.key = element_rect(fill="white",colour = NA))+
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  labs(x="", y="Proportion surviving", fill = "Offspring Group")+
  geom_text(data=SummarySurvival_Tile_Crunched, size = 6,(aes(x=factor(Cross), y= min(Survival)*0.8, label=paste("n =",N))))
survival_graph





##GLMM
Survival$TimePoint<-as.factor(Survival$TimePoint)
Survival_mod1<-glmer(Survival~Cross*Treatment*TimePoint+ (1|Tank/TileID), data=Survival, family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(Survival_mod1)

#Model with cross * treatment interaction 
Survival_mod2<-glmer(Survival~Cross*Treatment+TimePoint+ (1|Tank/TileID), data=Survival, family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(Survival_mod2)
anova(Survival_mod1,Survival_mod2)

#Model with cross * TimePoint interaction 
Survival_mod3<-glmer(Survival~Cross*TimePoint+Treatment+ (1|Tank/TileID), data=Survival, family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(Survival_mod3)
anova(Survival_mod1,Survival_mod3)
anova(Survival_mod2,Survival_mod3)

#Model with Treatment * TimePoint interaction 
Survival_mod4<-glmer(Survival~Cross+TimePoint*Treatment+ (1|Tank/TileID), data=Survival, family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(Survival_mod4)
anova(Survival_mod1,Survival_mod4)
anova(Survival_mod2,Survival_mod4)
anova(Survival_mod3,Survival_mod4)

#Model with no interactions
Survival_mod5<-glmer(Survival~Cross+Treatment+TimePoint+ (1|Tank/TileID), data=Survival, family=binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(Survival_mod5)
anova(Survival_mod1,Survival_mod5)
anova(Survival_mod2,Survival_mod5)
anova(Survival_mod3,Survival_mod5)
anova(Survival_mod4,Survival_mod5)

summary(glht(Survival_mod5, mcp(Cross="Tukey")))


###GROWTH###
Growth <- read.table('HeatWaveExperiment_Results_Singles_Growth.csv', header=T, sep=',')
str(Growth)

#exclude NAs
Growth<-na.omit(Growth)

##Summary, averaged by tile 
summaryGrowth_Tile<-summarySE(Growth, measurevar="Size", groupvars=c("Day","Cross","Tank","TileID","Treatment", "DHW"),conf.interval=0.95)
summaryGrowth_Tile
#Averaged by tile, then tank and then plotted
SummaryGrowth_Tile_Crunched<-summarySE(summaryGrowth_Tile, measurevar="Size", groupvars=c("Day","Cross","Treatment","DHW"),conf.interval=0.95)
SummaryGrowth_Tile_CrunchedAcrossTreatments<-summarySE(summaryGrowth_Tile, measurevar="Size", groupvars=c("Day","Cross","DHW"),conf.interval=0.95)
SummaryGrowth_Tile_CrunchedAcrossTreatments
summaryGrowth<-summarySE(Growth, measurevar="Size", groupvars=c("Day","Cross","Treatment","DHW"),conf.interval=0.95)
summaryGrowth

##Graph growth
yl <- expression(Area ~ mm^2)
#Graph by day
axes_marks_settings <- element_text(color = "black", size =10)
GrowthGraph=ggplot(summaryGrowth, aes(y=Size, x=Day, color=Cross, fill = Cross))+
  facet_grid(Treatment~.)+
  geom_ribbon(aes(ymin=Size-se, ymax=Size+se), position=position_dodge(0),alpha=0.8)+
  scale_fill_manual(values = cols, ,labels = c(expression(italic("A. loripes")), "KL hybrid",expression(italic("A. kenti")))) +
  scale_color_manual(values = c("#d95f02","#1b9e77","#7570b3")) +
  geom_point(size=2, fill = "black")+
  theme_classic()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),)+
  theme(legend.key = element_rect(fill="white",colour = NA))+
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  labs(x="Time point (days)", y=yl, color = "Offspring Group") +
  theme(axis.text = axes_marks_settings)+
  guides(fill = guide_legend("Offspring Group"), color = 'none') 
GrowthGraph
GrowthGraph + 
  geom_label_repel(label=summaryGrowth$N, fill = NA, size = 6, direction=c("both"), box.padding = 0.1, label.padding = 0.25, point.padding = 3, label.size = NA, arrow = NULL,show.legend = FALSE, fontface="bold" )

#Treat time point as a factor
Growth$TimePoint<- as.factor(Growth$TimePoint)

###lme using the nlme function and allowing for differences in variances between time points and crosses
Growth_mod1<-lme(Size~Cross*Treatment*TimePoint, random= ~1|Tank/TileID, data=Growth,weights=varIdent(form=~1|TimePoint*Cross))
summary(Growth_mod1)

#Model with Cross * Time Point interaction
Growth_mod2<-lme(Size~Cross*TimePoint+Treatment, random= ~1|Tank/TileID, data=Growth,weights=varIdent(form=~1|TimePoint*Cross))
summary(Growth_mod2)

#Model with Cross * Treatment interaction
Growth_mod3<-lme(Size~Cross*Treatment+TimePoint, random= ~1|Tank/TileID, data=Growth,weights=varIdent(form=~1|TimePoint*Cross))
summary(Growth_mod3)

#Model with Treatment * Time Point interaction
Growth_mod4<-lme(Size~Cross+TimePoint*Treatment, random= ~1|Tank/TileID, data=Growth,weights=varIdent(form=~1|TimePoint*Cross))
summary(Growth_mod4)

#Model with no interaction
Growth_mod5<-lme(Size~Cross+TimePoint+Treatment, random= ~1|Tank/TileID, data=Growth,weights=varIdent(form=~1|TimePoint*Cross))
summary(Growth_mod5)

##Modify models using Maximum Likelihood method to compare (instead of REML)
Growth_mod1_ml <- update(Growth_mod1, method="ML")
summary(Growth_mod1_ml)

Growth_mod2_ml <- update(Growth_mod2, method="ML")
summary(Growth_mod2_ml)
anova(Growth_mod1_ml,Growth_mod2_ml)

Growth_mod3_ml <- update(Growth_mod3, method="ML")
summary(Growth_mod3_ml)
anova(Growth_mod1_ml,Growth_mod3_ml)
anova(Growth_mod2_ml,Growth_mod3_ml)

Growth_mod4_ml <- update(Growth_mod4, method="ML")
summary(Growth_mod4_ml)
anova(Growth_mod1_ml,Growth_mod4_ml)
anova(Growth_mod2_ml,Growth_mod4_ml)
anova(Growth_mod3_ml,Growth_mod4_ml)

Growth_mod5_ml <- update(Growth_mod5, method="ML")
summary(Growth_mod5_ml)
anova(Growth_mod1_ml,Growth_mod5_ml)
anova(Growth_mod2_ml,Growth_mod5_ml)
anova(Growth_mod3_ml,Growth_mod5_ml)
anova(Growth_mod4_ml,Growth_mod5_ml)

##1 and 2 are equal best so simplified model --- Growth_mod2<-lme(Size~Cross*TimePoint+Treatment, random= ~1|Tank/TileID, data=Growth,weights=varIdent(form=~1|TimePoint*Cross)) --- is used
summary(Growth_mod2_ml)
summary(Growth_mod2)
  
leastsquare = lsmeans(Growth_mod2_ml,
                      pairwise ~ Cross:TimePoint,
                      adjust = "tukey")
leastsquare$contrasts

###Colour###
Colour <- read.table('HeatWaveExperiment_Results_Singles_Colour.csv', header=T, sep=',')
head(Colour)

#exclude NAs
Colour<-na.omit(Colour)

##Summary, averaged by tile 
summaryColour_Tile<-summarySE(Colour, measurevar="Colour", groupvars=c("Day","Cross","Tank","TileID","Treatment", "DHW"),conf.interval=0.95)
summaryColour_Tile
#Averaged by tile, then tank and then plotted
SummaryColour_Tile_Crunched<-summarySE(summaryColour_Tile, measurevar="Colour", groupvars=c("Day","Cross","Treatment","DHW"),conf.interval=0.95)
SummaryColour_Tile_Crunched

##Graph by day
summaryColour<-summarySE(Colour, measurevar="Colour", groupvars=c("Day","Cross", "Treatment","DHW"),conf.interval=0.95)
summaryColour

ColourGraph<-ggplot(summaryColour, aes(y=Colour, x=Day, color=Cross, fill = Cross))+
  facet_grid(Treatment~.)+
  geom_ribbon(aes(ymin=Colour-se, ymax=Colour+se), position=position_dodge(0),alpha=0.8)+
  scale_fill_manual(values = cols, labels = c(expression(italic("A. loripes")), "KL hybrid",expression(italic("A. kenti")))) +
  scale_color_manual(values = c("#d95f02","#1b9e77","#7570b3")) +
  geom_point(size=2, fill = "black")+
  theme_classic()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=16))+
  theme(axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),)+
  theme(legend.key = element_rect(fill="white",colour = NA))+
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  labs(x="Time point (days)", y="Colour score", color = "Offspring Group") +
  theme(axis.text = axes_marks_settings)+
  guides(fill = guide_legend("Offspring Group"), color = 'none') 
ColourGraph
ColourGraph + geom_label_repel(label=summaryColour$N, fill = NA, size = 6, direction=c("both"), box.padding = 0.1, label.padding = 0.5, point.padding = 3, label.size = NA, arrow = NULL,show.legend = FALSE, fontface="bold" )


##lme
Colour$TimePoint <- as.factor(Colour$TimePoint)
Colour_mod1<-lme(Colour~Cross*Treatment*TimePoint, random= ~1|Tank/TileID, data=Colour)
summary(Colour_mod1)

Colour_mod2<-lme(Colour~Cross*Treatment+TimePoint, random= ~1|Tank/TileID, data=Colour)
summary(Colour_mod2)

Colour_mod3<-lme(Colour~Cross+Treatment*TimePoint, random= ~1|Tank/TileID, data=Colour)
summary(Colour_mod3)

Colour_mod4<-lme(Colour~Cross*TimePoint+Treatment, random= ~1|Tank/TileID, data=Colour)
summary(Colour_mod4)

Colour_mod5<-lme(Colour~Cross+Treatment+TimePoint, random= ~1|Tank/TileID, data=Colour)
summary(Colour_mod5)

##Modify models using Maximum Likelihood method to compare (instead of REML)
Colour_mod1_ml <- update(Colour_mod1, method="ML")
summary(Colour_mod1_ml)

Colour_mod2_ml <- update(Colour_mod2, method="ML")
summary(Colour_mod2_ml)
anova(Colour_mod1_ml,Colour_mod2_ml)

Colour_mod3_ml <- update(Colour_mod3, method="ML")
summary(Colour_mod3_ml)
anova(Colour_mod1_ml,Colour_mod3_ml)
anova(Colour_mod2_ml,Colour_mod3_ml)

Colour_mod4_ml <- update(Colour_mod4, method="ML")
summary(Colour_mod4_ml)
anova(Colour_mod1_ml,Colour_mod4_ml)
anova(Colour_mod2_ml,Colour_mod4_ml)
anova(Colour_mod3_ml,Colour_mod4_ml)

Colour_mod5_ml <- update(Colour_mod5, method="ML")
summary(Colour_mod5_ml)
anova(Colour_mod1_ml,Colour_mod5_ml)
anova(Colour_mod2_ml,Colour_mod5_ml)
anova(Colour_mod3_ml,Colour_mod5_ml)
anova(Colour_mod4_ml,Colour_mod5_ml)

#Results for optimal model
summary(Colour_mod3_ml)
leastsquare = lsmeans(Colour_mod3_ml,
                      pairwise ~ Treatment:TimePoint,
                      adjust = "tukey")
leastsquare

###IPAM###
IPAM <- read.table('HybridPAM_PrelimResults_RecruitCentres_T4_Long.csv', header=T, sep=',')
head(IPAM)
#exclude NAs
IPAM<-na.omit(IPAM)
IPAM$Cross <- as.factor(IPAM$Cross)

#Summarise
summaryIPAM<-summarySE(IPAM, measurevar="Y", groupvars=c("TimePoint","Cross","Treatment","DHW"),conf.interval=0.95)

##Graph
IPAMGraph<-ggplot(summaryIPAM, aes(y=Y, x=TimePoint, color=Cross, fill = Cross))+
  facet_grid(Treatment~.)+
  geom_ribbon(aes(ymin=Y-se, ymax=Y+se), position=position_dodge(0),alpha=0.8)+
  scale_fill_manual(values = cols, labels = c(expression(italic("A. loripes")), "KL hybrid",expression(italic("A. kenti")))) +
  scale_color_manual(values = c("#d95f02","#1b9e77","#7570b3")) +
  geom_point(size=2, fill = "black")+
  theme_classic()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=16))+
  theme(axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),)+
  theme(legend.key = element_rect(fill="white",colour = NA))+
  theme(axis.text = element_text(colour = "black")) + theme(panel.background = element_rect(colour = "black", size=1)) + 
  labs(x="Time point (days)", y="Fv/Fm", color = "Offspring Group") +
  theme(axis.text = axes_marks_settings)+
  guides(fill = guide_legend("Offspring Group"), color = 'none') 
IPAMGraph
IPAMGraph + geom_label_repel(label=summaryIPAM$N, fill = NA, size = 6, direction=c("both"), box.padding = 0.1, label.padding = 0.25, point.padding = 3, label.size = NA, arrow = NULL,show.legend = FALSE, fontface="bold" )


##lme
IPAM_mod1<-lme(Y~Cross*Treatment*TimePoint, random= ~1|Tank/TileID, data=IPAM)
summary(IPAM_mod1)

IPAM_mod2<-lme(Y~Cross*Treatment+TimePoint, random= ~1|Tank/TileID, data=IPAM)
summary(IPAM_mod2)

IPAM_mod3<-lme(Y~Cross+Treatment*TimePoint, random= ~1|Tank/TileID, data=IPAM)
summary(IPAM_mod3)

IPAM_mod4<-lme(Y~Cross*TimePoint+Treatment, random= ~1|Tank/TileID, data=IPAM)
summary(IPAM_mod4)

IPAM_mod5<-lme(Y~Cross+Treatment+TimePoint, random= ~1|Tank/TileID, data=IPAM)
summary(IPAM_mod5)

#Compare models
IPAM_mod1_ml <- update(IPAM_mod1, method="ML")
summary(IPAM_mod1_ml)

IPAM_mod2_ml <- update(IPAM_mod2, method="ML")
anova(IPAM_mod1_ml,IPAM_mod2_ml)
summary(IPAM_mod2_ml)

IPAM_mod3_ml <- update(IPAM_mod3, method="ML")
anova(IPAM_mod1_ml,IPAM_mod3_ml)
anova(IPAM_mod2_ml,IPAM_mod3_ml)
summary(IPAM_mod3_ml)

IPAM_mod4_ml <- update(IPAM_mod4, method="ML")
summary(IPAM_mod4_ml)
anova(IPAM_mod1_ml,IPAM_mod4_ml)
anova(IPAM_mod4_ml,IPAM_mod2_ml)
anova(IPAM_mod3_ml,IPAM_mod4_ml)


IPAM_mod5_ml <- update(IPAM_mod5, method="ML")
summary(IPAM_mod5_ml)
anova(IPAM_mod1_ml,IPAM_mod5_ml)
anova(IPAM_mod2_ml,IPAM_mod5_ml)
anova(IPAM_mod3_ml,IPAM_mod5_ml)
anova(IPAM_mod4_ml,IPAM_mod5_ml)

#Summarise best model (mod4)
summary(IPAM_mod4_ml)
summ(IPAM_mod4_ml)
summary(glht(IPAM_mod4_ml, mcp(Cross="Tukey")))
