#RNA seq analysis
#Last Edited by Rebecca MacPherson 03-31-22


#####Filtering and normalization of counts data prior to statistical analysis#####

#load libraries
library(dplyr)
library(matrixStats)
library(data.table)
library(tidyr)

#import data
combined_counts <- read.table(file= "combined_ntr_counts.txt", header=TRUE)


# Import combined counts dataset
combined_counts <- as.data.frame(combined_counts)


# Adding rownames to combined_counts and renaming as Data

Data<-combined_counts
Data<-Data[,-1]
row.names(Data)<-combined_counts[,1]

#Remove the length column from Data and call it Data_nolength
Data_nolength<-Data[,-1]


#Data filtered by Median count>2
Datamatrix<-data.matrix(Data_nolength)
Medians<-rowMedians(Datamatrix)
Medians<-data.frame(Medians)
row.names(Medians)<-row.names(Data_nolength)
DatacbindMedian<-cbind(Data_nolength,Medians)
filterMedian<-filter(DatacbindMedian,Medians>=2)
filterMedian<-subset(filterMedian,select=-c(Medians))


#Filter out genes with proportion of non zero samples greater than 0.25
#18581 is the no of rows in filterMedian; "16" is the number of columns. 

colnotzero<-rep(0,18581)
colprop<-rep(0,18581)
for(i in 1:18581) {
  ctr<-0
  for(j in 1:16) {
    if(filterMedian[i,j] > 0) ctr=ctr+1
  }
  colnotzero[i]<-ctr
  colprop[i]<-colnotzero[i]/16.0
}
colprop<-data.frame(colprop)
filterMediancolprop<-cbind(filterMedian,colprop)
filtercolprop<-filter(filterMediancolprop,colprop>0.25)
filteredData<-subset(filtercolprop,select=-c(colprop))


#Merging based on Geneid
Geneid<-rownames(filteredData)
filteredData1<-cbind(filteredData,Geneid)
merged<-merge(filteredData1,combined_counts[c("Geneid","length")],by="Geneid")
rownames(merged)<-merged$Geneid
merged<-subset(merged,select=-c(Geneid))



# Prep data for normalization
Datafornormalization<-relocate(merged,"length", .before = "RAM.uhg4.208G.M.2_S16_L001_sorted.bam")
head(Datafornormalization)


write.csv(Datafornormalization, file='data_for_norm_ntr.csv', quote=FALSE)
#Datafornormalization includes Geneid (rownames) and Length (1st column of dataframe)



#####GeTMM normalization#####
##"The raw readcount matrix (tab-delimited text file) was used, in which the first column 
#holds the geneID from Ensembl that are used as row names in data matrix (x) in R, the second column 
#of the text file (thus the first column in x) 
#holds the gene length in kb and the remaining columns contain read counts of each sample." - (Smid et al., 2018)

# calculate RPK

x<-Datafornormalization
head(x)

rpk <- (x[,2:ncol(x)]/x[,1])

# remove length col in x
x <- x[,-1]
# for normalization purposes, no grouping of samples
group <- c(rep("A",ncol(x)))


#GeTMM
rpk.norm <- DGEList(counts=rpk,group=group)
rpk.norm <- calcNormFactors(rpk.norm)
norm.counts.rpk <- cpm(rpk.norm)

write.csv(norm.counts.rpk,file="norm_counts_GeTMM_uhg4_ntr.csv")


#####Transpose for input into SAS#####
a<-as.data.frame(norm.counts.rpk)
head(a)

#have to retain sample information within the df, not as column or row names
cnA<-colnames(a)
y<-rbind(cnA,a)
rny<-rownames(y)
z<-cbind(rny,y)
b<-transpose(z)

colnames(b)<-rownames(z)
b<-b[-1,] #remove the first row; that info is now in column names for pivot_longer
norm_counts_transposed<-b



norm_counts_forSAS<-pivot_longer(norm_counts_transposed,cols=starts_with("XLOC"),names_to = "geneid",values_to = "norm_counts")
write.csv(norm_counts_forSAS,file="norm_counts_forSAS_RNAseq_Uhg4_ntr.csv")
#separate metadata column into individual model terms (e.g. sex, genotypes) prior to running in SAS








##### SAS analysis ###########

#Modified by Rebecca MacPherson March 2022
#data below are in SAS notation
##########################################################################

/*Importing the data*/
  PROC IMPORT DATAFILE='/file/path.csv'
DBMS=CSV
OUT=norm_counts_forSAS replace;
GETNAMES=YES;
GUESSINGROWS=max;
RUN;

/* sort by gene id*/
PROC SORT DATA = norm_counts_forSAS OUT = WORK.norm_counts_sorted;
BY geneid;

/* format output tables*/
PROC template;
edit Stat.GLM.ProbF; format=E12.; end;
run;
PROC template;
edit Stat.GLM.LSMSlice; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Contrasts; format=E12.; end;
run;
PROC template;
edit Stat.GLM.Tests; format=E12.; end;
run;
proc template;
edit Stat.Mixed.Tests3; parent = Stat.Mixed.FTests; end;
run;


/*no html output or note writing to decrease run time*/
  ODS RESULTS OFF;
  options nonotes;

/*Naming output datasets*/
  ods noproctitle;
ods graphics / imagemap=off;
ods output FitStatistics=FS_ntr_01;
ods output LSMeans=LSMeans_ntr_01;
ods output ClassLevels= ClassLevels_ntr_01;
ods output NObs= NObs_ntr_01;
ods output OverallANOVA= OverallANOVA_ntr_01;
ods output ModelANOVA= ModelANOVA_ntr_01;
ods output SlicedANOVA= SlicedANOVA_ntr_01;
ods output Contrasts=Contrasts_ntr_01;


/*Defining and running the model for each gene (geneid)*/
proc glm data=WORK.norm_counts_sorted plots=none;
class line sex;
model norm_counts = line|sex;
by geneid;
lsmeans line*sex / slice= sex; 
contrast '208A - 208wt' line 1 0 0 -1;
contrast '208F - 208wt' line 0 1 0 -1;
contrast '208G - 208wt' line 0 0 1 -1;
contrast '208A - 208F' line 1 -1 0 0;
contrast '208A - 208G' line 1 0 -1 0;
contrast '208F - 208G' line 0 1 -1 0;
contrast '208A - 208wt Female' line 1 0 0 -1 line*sex 1 0 0 0 0 0 -1 0;
contrast '208F - 208wt Female' line 0 1 0 -1 line*sex 0 0 1 0 0 0 -1 0;
contrast '208G - 208wt Female' line 0 0 1 -1 line*sex 0 0 0 0 1 0 -1 0;
contrast '208A - 208F Female' line 1 -1 0 0 line*sex 1 0 -1 0 0 0 0 0;
contrast '208A - 208G Female' line 1 0 -1 0 line*sex 1 0 0 0 -1 0 0 0;
contrast '208F - 208G Female' line 0 1 -1 0 line*sex 0 0 1 0 -1 0 0 0;
contrast '208A - 208wt Male' line 1 0 0 -1 line*sex 0 1 0 0 0 0 0 -1;
contrast '208F - 208wt Male' line 0 1 0 -1 line*sex 0 0 0 1 0 0 0 -1;
contrast '208G - 208wt Male' line 0 0 1 -1 line*sex 0 0 0 0 0 1 0 -1;
contrast '208A - 208F Male' line 1 -1 0 0 line*sex 0 1 0 -1 0 0 0 0;
contrast '208A - 208G Male' line 1 0 -1 0 line*sex 0 1 0 0 0 -1 0 0;
contrast '208F - 208G Male' line 0 1 -1 0 line*sex 0 0 0 1 0 -1 0 0;
run;


ods trace output;
proc glm.....
run;
ods trace off;



##### Post SAS FDR correction and differential expression analyses ####
# Created by Rebecca MacPherson 2022-03-17
# Last edited by Rebecca MacPherson on 2022-03-17

#grab the initial ANOVA output file from the glm SAS analysis
setwd("file/path")

#output file from SAS is stacked improperly and includes Type I and Type III results
# 1-  Retain just the Type III results
# 2   Create lists of p values for each gene that correspond to each model term (Line, Sex, LinexSex)

model.anova.df<-read.csv("MODELANOVA_NTR_01.csv")
model.line<-model.anova.df[which(model.anova.df$Source=='line' & model.anova.df$HypothesisType =="3"),]
model.sex<-model.anova.df[which(model.anova.df$Source=='sex' & model.anova.df$HypothesisType =="3"),]
model.linexsex<-model.anova.df[which(model.anova.df$Source=='line*sex' & model.anova.df$HypothesisType =="3"),]


# Apply FDR correction to pvalues for each model term, separately

FDR_line<-p.adjust(model.line$ProbF,method="fdr")
FDR_sex<-p.adjust(model.sex$ProbF,method="fdr")
FDR_linexsex<-p.adjust(model.linexsex$ProbF,method="fdr")

# Combine new fdr p value columns with the rest of the metadata information (taken from model.line, since the first column is the same for "model.***")
gid<-as.data.frame(model.line$ï..geneid)

#make sure the order is the same
table(model.line$ï..geneid == model.linexsex$ï..geneid)
table(model.line$ï..geneid == model.sex$ï..geneid)

#combine into one df

id_pvalues<-cbind(gid, model.line$ProbF, FDR_line, model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex)
id_pvalues_export<-cbind(gid, model.line$ProbF, FDR_line, model.sex$ProbF, FDR_sex, model.linexsex$ProbF, FDR_linexsex)

colnames(id_pvalues)<-c("geneid","ProbF_line", "FDR_line", "ProbF_sex", "FDR_sex", "ProbF_linexsex", "FDR_linexsex")     
colnames(id_pvalues_export)<-c("geneid","ProbF_line_ntr", "FDR_line_ntr", "ProbF_sex", "FDR_sex", "ProbF_linexsex_ntr", "FDR_linexsex_ntr")     

#write table
setwd("file\path")
write.csv(id_pvalues_export, file="01_pvalues_FDR_allgenes_ntrs.csv")


#############################################################################

