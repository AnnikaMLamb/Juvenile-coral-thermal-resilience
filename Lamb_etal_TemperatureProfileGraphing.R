##R script for graphing temperature profile of Lamb et al. 
#Written by Annika Lamb 

#Load packages
library(RcppRoll)
library(ggplot2)
library(dplyr)
library(vctrs)

#Import data
Dir = "C:/Users/alamb/OneDrive - Australian Institute of Marine Science/PhD/Heatwave/Stats/TemperatureData"
setwd(Dir)
Temp <- read.table('FinalTemperatureProfile_DHWsIncluded.csv', header=T, sep=',', stringsAsFactors = FALSE)
str(Temp)

#Replace 'non-stressful' temperatures with hotspot values of 0
Temp$hotspot<-replace(Temp$DegreesAboveDaviesFebruaryAverage, Temp$DegreesAboveDaviesFebruaryAverage<1,0) 

#visualise data
plot(Temp$hotspot~Temp$Day)

#calculate DHW estimates
data.frame(Temp, c2=RcppRoll::roll_sum(Temp$hotspot,14, fill=NA))
Temp$c3<-cumsum(Temp$hotspot)
Temp$DHW<-Temp$c3/7

#Plot average daily tank temperatures
ggplot(aes(x=Day, y=ElevatedTemperature), data=Temp)+ ylab("Mean daily temperature") + scale_x_continuous (n.breaks = 12 , expand = c(0,0),limits = c ( 1 , 66 ) ) + geom_hline(yintercept=27.8, linetype="dashed", color = "black", size = 1.2)+ geom_line(size=1.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
                                                                                                                                                                                                                                                                           axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),                                                                                                                                                                                                                                                               legend.text = element_text(size = 14), legend.title = element_text(size = 14))
#Plot cumulative DHWs
ggplot(aes(x=Day, y=DHW), data=Temp)+ ylab("Cummulative degree heating weeks") + scale_x_continuous (n.breaks = 12 , expand = c(0,0),limits = c ( 1 , 66 ) ) + scale_y_continuous (n.breaks = 13 , expand = c(0,0),limits = c ( 0 , 13 ) ) + geom_line(size=1.2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
                                                                                                                                                                                                                                                                                                                                                                    axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),                                                                                                                                                                                                                                                                                                                                                         legend.text = element_text(size = 14), legend.title = element_text(size = 14))

#Combined plot of average daily temperatures and DHWs
# Plot both lines on one set of axes
combined_plot <- ggplot(Temp, aes(x = Day)) +
  geom_line(aes(y = ElevatedTemperature), color = "black", size = 1.2) +
  geom_line(aes(y = RampTemp), color = "blue", size = 1.2) +
  geom_line(aes(y = (DHW*0.308) + 27), color = "red", size = 1.2) +  # Multiply DHW by 5 to demonstrate the secondary axis
  scale_y_continuous(name = "Mean daily temperature", n.breaks = 6 , expand = c(0,0),limits = c(27, 31), sec.axis = sec_axis(~(.*3.25)-87.75, name = "Cumulative degree heating weeks")) +
  geom_hline(yintercept=27.8, linetype="dashed", color = "black", size = 1.2)+
   xlab("Day") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.line = element_line(color = "black"),
        axis.line.y.right = element_line(color = "red"), 
        axis.ticks.y.right = element_line(color = "red"),
        axis.text.y.right = element_text(color = "red"))+
  scale_x_continuous (n.breaks = 12 , expand = c(0,0),limits = c ( 1 , 66 ))

print(combined_plot)

