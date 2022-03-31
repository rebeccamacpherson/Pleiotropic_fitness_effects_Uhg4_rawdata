#Graphs for MacPherson et al., Pleiotropic Fitness Effects of the lncRNA Uhg4 in Drosophila melanogaster

#Last edited 2022_03-04 by Rebecca MacPherson
#Created by Rebecca MacPherson


#load packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

#### Stress Response Phenotypes ############

#Define objects used in ggplot command

Ulabels <- c("208-XA", "208-XF", "208-XG", "Control") #xaxis lables
my_palette = c(brewer.pal(9, "Set1")[c(8,2)]) #blue and pink colors
bold.italic<-element_text(face="bold.italic", size=10, colour = "black") #graph text
axestext<-element_text(size=12, color = "black") #graph text - axes
axestitle<-element_text(size=12) #graph text - axes title

########    Chill Coma ####

#import data
setwd("G:path\\to\\working\\directory")


gene <- read.csv("cc_ForSAS.csv", header=TRUE) 
gene1<-gene
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))

#create graph object
cc_u<-ggplot(data=gene1, aes(x=Line, y=(Time), fill=Sex)) + 
  geom_boxplot(outlier.size=0.5)+
  labs(x="Fly Line", y="Time (s)")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_discrete(limit = c(0, 500, 1000, 1500, 2000))+
  expand_limits(y=c(0,2040))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none", axis.text = axestext,
                          axis.title.x = element_blank(), axis.title.y = element_text(margin=margin(r=10)), 
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle,
                          axis.text.x = bold.italic)


########    Heat Shock ####

#import data
setwd("G:path\\to\\working\\directory")


gene <- read.csv("Full_Raw_Data23.csv", header=TRUE)
gene$Percent_Alive<-(gene$Percent_Alive*100)
gene1<-gene
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))

#create graph object
hs_u<-ggplot(data=gene1, aes(x=Line, y=(Percent_Alive), fill=Sex)) + 
  geom_bar(position = position_dodge(width=0.8), stat = "summary", fun = "mean", width=0.8, color='black', size=0.5)+
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width=0.8),width=0.2)+
  labs(x="Fly Line", y="Percent Surviving")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_discrete(expand = c(0,0),limit = c(0, 25, 50, 75, 100))+
  expand_limits(y=c(0, 102))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(margin=margin(r=10)), plot.title = element_text(hjust=0.5,size=30),
                          axis.title = axestitle ,axis.text.x = bold.italic)


########    Ethanol Sensitivity  ####

#import object
setwd("G:path\\to\\working\\directory")

gene <- read.csv("Raw_Data_Scored_EtOHSens.csv", header=TRUE)
gene1<-gene
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))

#create graph object
etoh_u<-ggplot(data=gene1, aes(x=Line, y=(Time), fill=Sex)) + 
  geom_boxplot(outlier.size=0.5)+
  # geom_violin()+
  labs(x="Fly Line", y="Time (s)")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_discrete(expand=c(0,0), limit = c(0, 250, 500, 750, 1000))+
  expand_limits(y=c(0,1020))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position= 'none', axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle,
                          axis.text.x = bold.italic)



####arranging into a panel #########

stress_u<-plot_grid(cc_u, hs_u, etoh_u,
                    nrow=3, align = "hv", vjust=1.5,
                    labels = c("A", "B", "C"))


############ Sleep and Activity Phenotypes #############

########    Night Sleep       ####################
#import data
setwd("G:path\\to\\working\\directory")


gene <- read.csv(file= "Ind_night_sleep_nodead.csv", header=TRUE)

gene1<-gene
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))

#create graph object
ns_u<-ggplot(data=gene1, aes(x=Line, y=(mean_sleep_per_ind), fill=Sex)) + 
  geom_boxplot(outlier.size=0.5)+
  labs(x="Fly Line", y="Proportion of Time\nAsleep, Night")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_discrete(expand = c(0,0), limit = c(0, 0.25, 0.5, 0.75, 1))+
  expand_limits(y=c(0, 1.15))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle ,axis.text.x = bold.italic)
ns_u


########    Day Sleep       ####################
#import data
setwd("G:path\\to\\working\\directory")


gene <- read.csv(file= "Ind_day_sleep_nodead.csv", header=TRUE)

gene1<-gene
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))

#create graph object
ds_u<-ggplot(data=gene1, aes(x=Line, y=(mean_sleep_per_ind), fill=Sex)) + 
  geom_boxplot(outlier.size=0.5)+
  labs(x="Fly Line", y="Proportion of Time\n Asleep, Day")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_discrete(expand = c(0,0), limit = c(0, 0.25, 0.5, 0.75, 1))+
  expand_limits(y=c(0, 1.15))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle ,axis.text.x = bold.italic)
ds_u


########    Sleep Bout Count, Night   ###########
#import data
setwd("G:path\\to\\working\\directory")

gene <- read.csv("Ind_sleep_bout_nodead_bout_counts_night_compiled.csv", header=TRUE)
gene1<-gene
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))

#create graph object
sbc_n_u<-ggplot(data=gene1, aes(x=Line, y=(bout_count), fill=Sex)) + 
  geom_boxplot(outlier.size=0.5)+
  labs(x="Fly Line", y="Sleep Bout Count,\nNight")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_discrete(expand = c(0,0), limit = c(0, 5, 10, 15, 20))+
  expand_limits(y=c(0,23))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle ,axis.text.x = bold.italic)
sbc_n_u



########    Sleep Bout Count, Day   ###########
#import data
setwd("G:path\\to\\working\\directory")

gene <- read.csv("Ind_sleep_bout_nodead_bout_counts_day_compiled.csv", header=TRUE)
gene1<-gene
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))
#create graph object
sbc_d_u<-ggplot(data=gene1, aes(x=Line, y=(bout_count), fill=Sex)) + 
  geom_boxplot(outlier.size=0.5)+
  labs(x="Fly Line", y="Sleep Bout Count,\nDay")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_discrete(expand = c(0,0), limit = c(0, 5, 10, 15, 20))+
  expand_limits(y=c(0,23))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle ,axis.text.x = bold.italic)
sbc_d_u


########    Activity Bout Length, Day    #######
#import data
setwd("G:path\\to\\working\\directory")

gene <- read.csv("Ind_activity_bout_nodead_day_compiled.csv", header=TRUE)

gene1<-gene
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))

#create graph object
abl_d_u<-ggplot(data=gene1, aes(x=Line, y=(bout_length), fill=Sex)) + 
  geom_boxplot(outlier.size=0.5)+
  labs(x="Fly Line", y="Activity Bout Length,\nDay (min)")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_discrete(expand=c(0,0),limit = c(0, 50, 100, 150, 200))+
  expand_limits(y=c(0,230))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle ,axis.text.x = bold.italic)
abl_d_u


########    Locomotor Activity #######################
#import data
setwd("G:path\\to\\working\\directory")

gene <- read.csv("Ind_daily_locomotor_activity_data_nodead_compiled.csv", header=TRUE)
gene1<-gene
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))

#create graph object
act_u<-ggplot(data=gene1, aes(x=Line, y=(Activity), fill=Sex)) + 
  geom_boxplot(outlier.size=0.5)+
  labs(x="Fly Line", y="Total Activity\n(Counts)")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_discrete(expand=c(0,0),limit = c(0, 500, 1000, 1500, 2000, 2500))+
  expand_limits(y=c(0,2850))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle ,axis.text.x = bold.italic)
act_u



#### plotting on a panel     ##########################

s_act_u<-plot_grid(ds_u, ns_u, sbc_d_u,  sbc_n_u, abl_d_u, act_u,
                   ncol=2, nrow=3, align = "hv",
                   vjust = 1, labels = c("A", "B", "C", "D", "E", "F"))


########## Development Viability Phenotypes ############

#load packageds
library(ggplot2)
library(dplyr)
library(readxl)
library(RColorBrewer)
library(cowplot)

#######    Development####
#import data
setwd("G:path\\to\\working\\directory")


data_line <- read_excel("linechart.xlsx")

#defining objects for use in the ggplot command
Ulabels <- c("208-XA", "208-XF", "208-XG", "Control") #xaxis lables
axestext<-element_text(size=12, color = "black") #graph text - axes
axestitle<-element_text(size=12) #graph text - axes title
subtitle<-element_text(hjust=0.5,size=12, face="bold.italic")
my_palette_viability = c((brewer.pal(9, "Set1")[c(8)]),(brewer.pal(9, "Pastel1")[c(8)]),
                         (brewer.pal(9, "Set1")[c(2)]),(brewer.pal(9, "Pastel1")[c(2)]))


#subsetting data to only include homozygous flies with straight wings for use in graphs
datawtF<-data_line[which(data_line$Line == "208wt" | data_line$Line == "208F"),1:7]
datawtA<-data_line[which(data_line$Line == "208wt" | data_line$Line == "208A"),1:7]
datawtG<-data_line[which(data_line$Line == "208wt" | data_line$Line == "208G"),1:7]

datawtF_S<-datawtF[which(datawtF$Wing == "S"),]
datawtA_S<-datawtA[which(datawtA$Wing == "S"),]
datawtG_S<-datawtG[which(datawtG$Wing == "S"),]

dataF_S<-datawtF_S[which(datawtF_S$Line == "208F"),]
dataA_S<-datawtA_S[which(datawtA_S$Line == "208A"),]
dataG_S<-datawtG_S[which(datawtG_S$Line == "208G"),]
datawt_S<-datawtF_S[which(datawtF_S$Line == "208wt"),]

#making graphs
#208wt
d_u_208wt<-ggplot(data=datawtG_S[which(datawtG_S$Line=="208wt"),], aes(x=Day, y=Number, fill=Sex)) +
  geom_col(color="black") +
  scale_fill_manual(values=my_palette_viability[c(1,3)])+
  scale_y_continuous(expand=c(0,0),breaks=c(0,2,4,6,8,10))+
  expand_limits(y=c(0,10))+
  labs(x="Day", y="Number of Flies", subtitle = "Control")+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = axestitle, axis.title.y = element_text(margin=margin(r=10)),
                          plot.subtitle = subtitle, axis.title = axestitle ,axis.text.x = axestext)

#208A
d_u_208A<-ggplot(data=datawtA_S[which(datawtA_S$Line=="208A"),], aes(x=Day, y=Number, fill=Sex)) +
  geom_col(color="black") +
  scale_fill_manual(values=my_palette_viability[c(1,3)])+
  scale_y_continuous(expand=c(0,0),breaks=c(0,2,4,6,8,10))+
  expand_limits(y=c(0,10))+
  labs(x="Day", y="Number of Flies", subtitle = "208-XA")+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = axestitle, axis.title.y = element_text(margin=margin(r=10)),
                          plot.subtitle = subtitle, axis.title = axestitle ,axis.text.x = axestext)

#208F
d_u_208F<-ggplot(data=datawtF_S[which(datawtF_S$Line=="208F"),], aes(x=Day, y=Number, fill=Sex)) +
  geom_col(color="black") +
  scale_fill_manual(values=my_palette_viability[c(1,3)])+
  scale_y_continuous(expand=c(0,0),breaks=c(0,2,4,6,8,10))+
  expand_limits(y=c(0,10))+
  labs(x="Day", y="Number of Flies", subtitle="208-XF")+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = axestitle, axis.title.y = element_text(margin=margin(r=10)),
                          plot.subtitle = subtitle, axis.title = axestitle ,axis.text.x = axestext)

#208G
  d_u_208G<-ggplot(data=datawtG_S[which(datawtG_S$Line=="208G"),], aes(x=Day, y=Number, fill=Sex)) +
  geom_col(color="black") +
  scale_fill_manual(values=my_palette_viability[c(1,3)])+
  scale_y_continuous(expand=c(0,0),breaks=c(0,2,4,6,8,10))+
  expand_limits(y=c(0,10))+
  labs(x="Day", y="Number of Flies", subtitle="208-XG")+
    theme_classic() + theme(legend.position = "none",axis.text = axestext,
                            axis.title.x = axestitle, axis.title.y = element_text(margin=margin(r=10)),
                            plot.subtitle = subtitle, axis.title = axestitle ,axis.text.x = axestext)

#development final plot
d_u_byline<-plot_grid(d_u_208A, d_u_208F, d_u_208G,d_u_208wt, align="hv", ncol=4)





####        Viability #####
#import data
setwd("G:path\\to\\working\\directory")


gene <- read.csv("Viability_values_forgraphing.csv", header=TRUE)

#defining variables for use in graphing plots
gene1<-gene
Ulabels <- c("208-XA", "208-XF", "208-XG", "Control") #xaxis lables
my_palette = c(brewer.pal(9, "Set1")[c(8,2)]) #blue and pink colors
bold.italic<-element_text(face="bold.italic", size=10, colour = "black") #graph text

my_p_males = c(brewer.pal(9, "Set1")[c(9)])
my_p_females = c(brewer.pal(9, "Set1")[c(4)])
gene1$Line = factor(gene1$Line, 
                    levels = c("208A", "208F", "208G", "208wt"))


# viability both sexes (total)
v_u<-ggplot(data=gene1, aes(x=Line, y=(v_total))) + 
  geom_boxplot(outlier.size=0.5)+
  labs(x="Fly Line", y="Viability")+
  scale_x_discrete(labels = Ulabels)+
  scale_y_continuous(expand=c(0,0))+
    expand_limits(y=c(0,2.1))+
  scale_fill_manual(values=my_palette)+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle ,axis.text.x = bold.italic)


# viability - females
v_u_f<-ggplot(data=gene1, aes(x=Line, y=(v_female))) + 
  geom_boxplot(fill=my_palette[1],outlier.size=0.5)+
  labs(x="Fly Line", y="Viability")+
  scale_x_discrete(labels = Ulabels)+
  expand_limits(y=c(0,2.1))+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle ,axis.text.x = bold.italic)


# viability - males
v_u_m<-ggplot(data=gene1, aes(x=Line, y=(v_male))) + 
  geom_boxplot(fill=my_palette[2], outlier.size=0.5)+
  labs(x="Fly Line", y="Viability")+
  scale_x_discrete(labels = Ulabels)+
  expand_limits(y=c(0,2.1))+
  theme_classic() + theme(legend.position = "none",axis.text = axestext,
                          axis.title.x = element_blank(),axis.title.y = element_text(margin=margin(r=10)),
                          plot.title = element_text(hjust=0.5,size=30), axis.title = axestitle ,axis.text.x = bold.italic)


v_u_bothsexes<-plot_grid(v_u_f, v_u_m, align="hv", nrow=1, ncol=2, labels = c("B", "C"))


######### final plotting ####################


plot_grid(d_u_byline, v_u_bothsexes, nrow=2, labels = "A")    #viability,development time, deletion


