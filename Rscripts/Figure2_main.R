# Figure 2 Genetic Load Downsampled to 10 x -
# only sites present in historical (removed exclusively modern sites)
# also other test not included in the paper. 
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(ggpubr)
library(bootstrap)

setwd("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/")
metadata <- read_excel("Metadata_crane.xlsx", sheet = 1)

samples<- read.table("Files/snpEff/Samples_down")
samples_cov <- c("EB31_S103_15x","EB31_S103_5x", "EB31_S103",
                 "EB33_S105_15x","EB33_S105_5x","EB33_S105")
xxx <- "RmExclusivModern"
#xxx <- "OnlySharedVariants"

#Change in Frequency of deleterious variation through time-----
listfreq3 <- read.csv2("Files/Change_frequency.csv") # removed sites exclusively modern
#listfreq3<- listfreq3[which(listfreq3$FreqHist>0&listfreq3$FreqModern>0),]  # only sites exclusivly shared

listfreq3$Type <- factor(listfreq3$Type, levels=c("Low","Moderate","High"), ordered = T)

# Plot distribution
freq_dist <- ggplot(listfreq3, aes(Difference, color=Type))
freq_dist <- freq_dist + geom_density() + theme_classic() + 
  scale_color_manual(values=c("#4f7cac","#7dde92","#696047"))+ xlab("Delta Frequency (Modern - Historical)")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")+ theme(legend.position = "top", legend.title = element_blank())

freq_dist

#percentage with higher frequency in modern
x<-length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="High",]$Variant)/length(listfreq3[listfreq3$Type=="High",]$Variant)
y<-length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="Moderate",]$Variant)/length(listfreq3[listfreq3$Type=="Moderate",]$Variant)
z<-length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="Low",]$Variant)/length(listfreq3[listfreq3$Type=="Low",]$Variant)

x<- length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="High",]$Variant)
xxx<- length(listfreq3[listfreq3$Difference<0 & listfreq3$Type=="High",]$Variant)
xx <-length(listfreq3[listfreq3$Type=="High",]$Variant)

y<-length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="Moderate",]$Variant)
yyy<-length(listfreq3[listfreq3$Difference<0 & listfreq3$Type=="Moderate",]$Variant)
yy <-length(listfreq3[listfreq3$Type=="Moderate",]$Variant)

z<-length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="Low",]$Variant)
zzz<-length(listfreq3[listfreq3$Difference<0 & listfreq3$Type=="Low",]$Variant)
zz<-length(listfreq3[listfreq3$Type=="Low",]$Variant)



df <- rbind.data.frame(c("High",x,xxx,xx),c("Moderate",y,yyy,yy),c("Low",z,zzz,zz))
colnames(df) <- c("Type","VariantsUP","VariantsDown","TotalVariants")

chisq <- chisq.test(as.numeric(df$VariantsUP), as.numeric(df$VariantsDown))
chisq



mean(listfreq3[listfreq3$Type=="High",]$Difference)
sd(listfreq3[listfreq3$Type=="High",]$Difference)
mean(listfreq3[listfreq3$Type=="Moderate",]$Difference)
mean(listfreq3[listfreq3$Type=="Low",]$Difference)

# Rxy ------
#read files
high_freq <- read.csv2("Files/Frequency_derived_high.csv", sep=",") # obtained in Load_frequencychange.R
moderate_freq <- read.csv2("Files/Frequency_derived_moderate.csv", sep=",")
low_freq <- read.csv2("Files/Frequency_derived_low.csv", sep=",")


#exclusivly present in both groups, not just not in modern - NOT INCLUDED.
#high_freq <- read.csv2("Files/Frequency_derived_high_nosingletons.csv", sep=",") # obtained in Load_frequencychange.R
#moderate_freq <- read.csv2("Files/Frequency_derived_moderate_nosingletons.csv", sep=",")
#low_freq <- read.csv2("Files/Frequency_derived_low_nosingletons.csv", sep=",")

high_freq$Type <- "High"
moderate_freq$Type <- "Moderate"
low_freq$Type <- "Low"

Frequencies_all <- rbind(high_freq,moderate_freq,low_freq)
Frequencies_all$FreqHist <- as.numeric(Frequencies_all$FreqHist)
Frequencies_all$FreqModern <- as.numeric(Frequencies_all$FreqModern)

# This only to use for the 1st case
Frequencies_all <- Frequencies_all[-which(Frequencies_all$FreqHist==0 & Frequencies_all$FreqModern>0),]

##Rxy
types <- c("High","Moderate", "Low")

Rxy_frame<-data.frame(matrix(ncol = 4, nrow = 0))
Rxy_jn<-data.frame(matrix(ncol = 3, nrow = 0))

y<-Frequencies_all$FreqHist*(1-Frequencies_all$FreqModern)
x<-Frequencies_all$FreqModern*(1-Frequencies_all$FreqHist)
xdata_total<-cbind.data.frame(x,y, Type=Frequencies_all$Type)

for (j in types){
  print(j)
  xdata <- xdata_total[xdata_total$Type==j,]
 
  size <- floor(length(xdata$x)*0.01) # take 1% of the sites for the jacknife
  #size <- ceiling(length(xdata$x)*0.01) # take 1% of the sites for the jacknife ONLY FOR THE SHARED FILES
  rxy<-function(a,xdata,b){
    sum(xdata[seq(a,a+(size-1),1),b])
  } 
  rx_ry<-data.frame()
  for(i in c(1:100)) {
    rx_ry[i,1]<-rxy(size*(i-1)+1,xdata,1)
    rx_ry[i,2]<-rxy(size*(i-1)+1,xdata,2)
  } 
  
  rat<-function(c,rx_ry){sum(rx_ry[c,1])/sum(rx_ry[c,2])}
  jn_h<-jackknife(1:100,rat,rx_ry) #do jackknife for dataframe rx_ry with the function of rat for 100 times
  Rxy_jn<-rbind(Rxy_jn, cbind(jn_h$jack.values, rep(j, 100)))
  Rxy_frame<-rbind(Rxy_frame, cbind(j, sum(xdata$x)/sum(xdata$y), jn_h$jack.se))
} 

colnames(Rxy_frame)<-c("Type","Rxy", "SE")
Rxy_frame$Rxy <- as.numeric(Rxy_frame$Rxy)
Rxy_frame$SE <- as.numeric(Rxy_frame$SE)

# Normalize by low 
Rxy_frame_low <- Rxy_frame[Rxy_frame$Type=="Low",]
Rxy_frame_notlow <- Rxy_frame[Rxy_frame$Type!="Low",]

Rxy_frame_normalized<- data.frame(Type=c("High","Moderate"),
                                  Rxy_norm=as.numeric(Rxy_frame_notlow$Rxy)/as.numeric(Rxy_frame_low$Rxy),
                                  SE_norm=as.numeric(Rxy_frame_notlow$S)/as.numeric(Rxy_frame_low$Rxy))

cut.off = 1
Rxy_frame_normalized$Rxy_norm_2 <- Rxy_frame_normalized$Rxy_norm -1

rxy_plot <- ggplot(Rxy_frame_normalized, aes(Type, Rxy_norm_2, fill=Type))
rxy_plot <- rxy_plot + geom_col(position = "dodge") + theme_classic()+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) +scale_fill_manual(values=c("#696047","#7dde92")) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + ylab("Rxy (Modern/Historical)") + theme(legend.position = "none")
rxy_plot

#HIGH --------------
high <- read.table("Files/snpEff/high_down10x.txt")
colnames(high) <- c("Chrom","Pos","REF","ALT", samples$V1)

high$chrompos<- paste0(high$Chrom,":",high$Pos)


high_fixed0 <- high %>% filter(grepl("0/0:[0-9]+", EB31_S103_15x)) %>% 
  filter(grepl("0/0:[0-9]+", EB31_S103_5x)) %>% filter(grepl("0/0:[0-9]+", EB33_S105_15x)) %>% 
  filter(grepl("0/0:[0-9]+", EB33_S105_5x))  %>% 
  filter(grepl("0/0:[0-9]+", EB11_S84)) %>%   filter(grepl("0/0:[0-9]+", EB12_S85)) %>% 
  filter(grepl("0/0:[0-9]+", EB13_S86)) %>%  filter(grepl("0/0:[0-9]+", EB14_S87)) %>% 
  filter(grepl("0/0:[0-9]+", EB15_S88))  %>%  filter(grepl("0/0:[0-9]+", EB16_S89)) %>% 
  filter(grepl("0/0:[0-9]+", EB17_S90))  %>%  filter(grepl("0/0:[0-9]+",EB18_S91 )) %>% 
  filter(grepl("0/0:[0-9]+", EB19_S92))  %>%  filter(grepl("0/0:[0-9]+",EB1_S114 )) %>% 
  filter(grepl("0/0:[0-9]+", EB20_S93))  %>%  filter(grepl("0/0:[0-9]+", EB21_S94)) %>% 
  filter(grepl("0/0:[0-9]+", EB22_S95))  %>%  filter(grepl("0/0:[0-9]+", EB24_S96)) %>% 
  filter(grepl("0/0:[0-9]+", EB25_S97))  %>%  filter(grepl("0/0:[0-9]+", EB26_S98)) %>% 
  filter(grepl("0/0:[0-9]+", EB27_S99))  %>%  filter(grepl("0/0:[0-9]+",EB28_S100)) %>% 
  filter(grepl("0/0:[0-9]+", EB29_S101))  %>%  filter(grepl("0/0:[0-9]+", EB2_S115)) %>% 
  filter(grepl("0/0:[0-9]+", EB31_S103)) %>% 
  filter(grepl("0/0:[0-9]+", EB32_S104))  %>%  filter(grepl("0/0:[0-9]+", EB33_S105)) %>% 
  filter(grepl("0/0:[0-9]+", EB34_S106))  %>%  filter(grepl("0/0:[0-9]+", EB35_S107)) %>% 
  filter(grepl("0/0:[0-9]+", EB37_S108))  %>%  filter(grepl("0/0:[0-9]+", EB38_S109)) %>% 
  filter(grepl("0/0:[0-9]+", EB39_S110))  %>%  filter(grepl("0/0:[0-9]+", EB3_S116)) %>% 
  filter(grepl("0/0:[0-9]+", EB40_S111))  %>%  filter(grepl("0/0:[0-9]+", EB41_S112)) %>% 
  filter(grepl("0/0:[0-9]+", EB42_S113))  %>%  filter(grepl("0/0:[0-9]+",EB6_S117 )) %>% 
  filter(grepl("0/0:[0-9]+", EB7_S118))  %>%  filter(grepl("0/0:[0-9]+", EB8_S119)) %>% 
  filter(grepl("0/0:[0-9]+", EB9_S83))  

fixed <- high_fixed0$chrompos

high <- high[which(!high$chrompos %in% fixed),]


highcounts<-list()
highcounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- high[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,high[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

highcounts_df <- do.call(rbind, highcounts)

#MODERATE --------------
moderate <- read.table("Files/snpEff/moderate_down10x.txt")
colnames(moderate) <- c("Chrom","Pos","REF","ALT", samples$V1)
moderate$chrompos<- paste0(moderate$Chrom,":",moderate$Pos)


moderate_fixed0 <- moderate %>% filter(grepl("0/0:[0-9]+", EB31_S103_15x)) %>% 
  filter(grepl("0/0:[0-9]+", EB31_S103_5x)) %>% filter(grepl("0/0:[0-9]+", EB33_S105_15x)) %>% 
  filter(grepl("0/0:[0-9]+", EB33_S105_5x))  %>% 
  filter(grepl("0/0:[0-9]+", EB11_S84)) %>%   filter(grepl("0/0:[0-9]+", EB12_S85)) %>% 
  filter(grepl("0/0:[0-9]+", EB13_S86)) %>%  filter(grepl("0/0:[0-9]+", EB14_S87)) %>% 
  filter(grepl("0/0:[0-9]+", EB15_S88))  %>%  filter(grepl("0/0:[0-9]+", EB16_S89)) %>% 
  filter(grepl("0/0:[0-9]+", EB17_S90))  %>%  filter(grepl("0/0:[0-9]+",EB18_S91 )) %>% 
  filter(grepl("0/0:[0-9]+", EB19_S92))  %>%  filter(grepl("0/0:[0-9]+",EB1_S114 )) %>% 
  filter(grepl("0/0:[0-9]+", EB20_S93))  %>%  filter(grepl("0/0:[0-9]+", EB21_S94)) %>% 
  filter(grepl("0/0:[0-9]+", EB22_S95))  %>%  filter(grepl("0/0:[0-9]+", EB24_S96)) %>% 
  filter(grepl("0/0:[0-9]+", EB25_S97))  %>%  filter(grepl("0/0:[0-9]+", EB26_S98)) %>% 
  filter(grepl("0/0:[0-9]+", EB27_S99))  %>%  filter(grepl("0/0:[0-9]+",EB28_S100)) %>% 
  filter(grepl("0/0:[0-9]+", EB29_S101))  %>%  filter(grepl("0/0:[0-9]+", EB2_S115)) %>% 
  filter(grepl("0/0:[0-9]+", EB31_S103)) %>% 
  filter(grepl("0/0:[0-9]+", EB32_S104))  %>%  filter(grepl("0/0:[0-9]+", EB33_S105)) %>% 
  filter(grepl("0/0:[0-9]+", EB34_S106))  %>%  filter(grepl("0/0:[0-9]+", EB35_S107)) %>% 
  filter(grepl("0/0:[0-9]+", EB37_S108))  %>%  filter(grepl("0/0:[0-9]+", EB38_S109)) %>% 
  filter(grepl("0/0:[0-9]+", EB39_S110))  %>%  filter(grepl("0/0:[0-9]+", EB3_S116)) %>% 
  filter(grepl("0/0:[0-9]+", EB40_S111))  %>%  filter(grepl("0/0:[0-9]+", EB41_S112)) %>% 
  filter(grepl("0/0:[0-9]+", EB42_S113))  %>%  filter(grepl("0/0:[0-9]+",EB6_S117 )) %>% 
  filter(grepl("0/0:[0-9]+", EB7_S118))  %>%  filter(grepl("0/0:[0-9]+", EB8_S119)) %>% 
  filter(grepl("0/0:[0-9]+", EB9_S83))  

fixed_moderate <- moderate_fixed0$chrompos

moderate <- moderate[which(!moderate$chrompos %in% fixed_moderate),]

moderatecounts<-list()
moderatecounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- moderate[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,moderate[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

moderatecounts_df <- do.call(rbind, moderatecounts)

#LOW --------------
low <- read.table("Files/snpEff/low_down10x.txt")
colnames(low) <- c("Chrom","Pos","REF","ALT", samples$V1)
low$chrompos<- paste0(low$Chrom,":",low$Pos)


low_fixed0 <- low %>% filter(grepl("0/0:[0-9]+", EB31_S103_15x)) %>% 
  filter(grepl("0/0:[0-9]+", EB31_S103_5x)) %>% filter(grepl("0/0:[0-9]+", EB33_S105_15x)) %>% 
  filter(grepl("0/0:[0-9]+", EB33_S105_5x))  %>% 
  filter(grepl("0/0:[0-9]+", EB11_S84)) %>%   filter(grepl("0/0:[0-9]+", EB12_S85)) %>% 
  filter(grepl("0/0:[0-9]+", EB13_S86)) %>%  filter(grepl("0/0:[0-9]+", EB14_S87)) %>% 
  filter(grepl("0/0:[0-9]+", EB15_S88))  %>%  filter(grepl("0/0:[0-9]+", EB16_S89)) %>% 
  filter(grepl("0/0:[0-9]+", EB17_S90))  %>%  filter(grepl("0/0:[0-9]+",EB18_S91 )) %>% 
  filter(grepl("0/0:[0-9]+", EB19_S92))  %>%  filter(grepl("0/0:[0-9]+",EB1_S114 )) %>% 
  filter(grepl("0/0:[0-9]+", EB20_S93))  %>%  filter(grepl("0/0:[0-9]+", EB21_S94)) %>% 
  filter(grepl("0/0:[0-9]+", EB22_S95))  %>%  filter(grepl("0/0:[0-9]+", EB24_S96)) %>% 
  filter(grepl("0/0:[0-9]+", EB25_S97))  %>%  filter(grepl("0/0:[0-9]+", EB26_S98)) %>% 
  filter(grepl("0/0:[0-9]+", EB27_S99))  %>%  filter(grepl("0/0:[0-9]+",EB28_S100)) %>% 
  filter(grepl("0/0:[0-9]+", EB29_S101))  %>%  filter(grepl("0/0:[0-9]+", EB2_S115)) %>% 
  filter(grepl("0/0:[0-9]+", EB31_S103)) %>% 
  filter(grepl("0/0:[0-9]+", EB32_S104))  %>%  filter(grepl("0/0:[0-9]+", EB33_S105)) %>% 
  filter(grepl("0/0:[0-9]+", EB34_S106))  %>%  filter(grepl("0/0:[0-9]+", EB35_S107)) %>% 
  filter(grepl("0/0:[0-9]+", EB37_S108))  %>%  filter(grepl("0/0:[0-9]+", EB38_S109)) %>% 
  filter(grepl("0/0:[0-9]+", EB39_S110))  %>%  filter(grepl("0/0:[0-9]+", EB3_S116)) %>% 
  filter(grepl("0/0:[0-9]+", EB40_S111))  %>%  filter(grepl("0/0:[0-9]+", EB41_S112)) %>% 
  filter(grepl("0/0:[0-9]+", EB42_S113))  %>%  filter(grepl("0/0:[0-9]+",EB6_S117 )) %>% 
  filter(grepl("0/0:[0-9]+", EB7_S118))  %>%  filter(grepl("0/0:[0-9]+", EB8_S119)) %>% 
  filter(grepl("0/0:[0-9]+", EB9_S83))  


fixed_low <- low_fixed0$chrompos

low <- low[which(!low$chrompos %in% fixed_low),]

lowcounts<-list()
lowcounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- low[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,low[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

lowcounts_df <- do.call(rbind, lowcounts)


# Repeat count of deleterious without the sites exclusively of modern ----
#only need to filter the modern, historicals are all positions
sites <- gsub("-",":",listfreq3$Variant)

# filter from the global count 
high_shared <- high[which(high$chrompos %in% sites),]
#high_shared<- high

highcounts<-lihighhighcounts<-list()
highcounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- high_shared[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,high_shared[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

highcounts_df_shared <- do.call(rbind, highcounts)


moderate_shared <- moderate[which(moderate$chrompos %in% sites),]
#moderate_shared <- moderate

moderatecounts<-list()
moderatecounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- moderate_shared[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,moderate_shared[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

moderatecounts_df_shared <- do.call(rbind, moderatecounts)

low_shared <- low[which(low$chrompos %in% sites),]
#low_shared <- low

lowcounts<-list()
lowcounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- low_shared[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,low_shared[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

lowcounts_df_shared <- do.call(rbind, lowcounts)

# plotting
lowcounts_df_shared$Class <- "Low"
moderatecounts_df_shared$Class <- "Moderate"
highcounts_df_shared$Class <- "High"


counts_df_shared <- rbind(highcounts_df_shared,moderatecounts_df_shared, lowcounts_df_shared)
counts_df_shared$Total_5 <-counts_df_shared$het_refanc_5 +counts_df_shared$hom_refanc_5

counts_df_low_shared <- counts_df_shared[counts_df_shared$Class=="Low",]

counts_df_shared$HetCountNormLow_5 <- counts_df_shared$het_refanc_5/counts_df_low_shared$Total_5
counts_df_shared$HomCountNormLow_5 <- counts_df_shared$hom_refanc_5/counts_df_low_shared$Total_5
counts_df_shared$TotalNormLow_5 <- counts_df_shared$Total_5/counts_df_low_shared$Total_5


counts_down_shared <- counts_df_shared[,c("SampleID", "Class",
                                          "het_refanc_5","hom_refanc_5","Total_5",
                                          "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5")]
counts_down_shared$Dataset <- "Downsampled"

all_data <- read.csv2("Files/snpEff/Counts_totalCoverage.csv") # Add info on historical samples
all_data <- all_data[,c("SampleID", "Class",
                        "het_refanc_5","hom_refanc_5","Total_5",
                        "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5", "Dataset")]

historical_samples <- c("Grus27863","Grus71441","MCZ-220028", "MCZ-301263","MCZ-35730",
                        "MCZ-42511","MCZ-46873" )

all_data <- all_data[which(all_data$SampleID%in%historical_samples),] # keep only the historical samples in the whole dataset


total_counts_shared <- rbind.data.frame(all_data,counts_down_shared)
samples_cov <- c("EB31_S103_15x","EB31_S103_5x", "EB31_S103",
                 "EB33_S105_15x","EB33_S105_5x","EB33_S105")
total_counts_down_shared <-total_counts_shared[-which(total_counts_shared$SampleID%in%samples_cov[-c(3,6)]),]


total_counts_df3_shared <- total_counts_down_shared[,c("SampleID", "Class",
                                                       "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5","Dataset")]


counts_df_gather_shared <- gather(total_counts_df3_shared, State, Count, -SampleID,  -Class, -Dataset,)

countsnorm_df_metadata_shared <- merge(counts_df_gather_shared, metadata[,c("SampleID","Category")], by="SampleID") 
countsnorm_df_metadata_shared$Class <- factor(countsnorm_df_metadata_shared$Class, levels=c("Low","Moderate","High"), ordered = T)

countsnorm_df_metadata_shared$Type <- countsnorm_df_metadata_shared$Category
countsnorm_df_metadata_shared$Type <- gsub("Founders_wildborn","Modern",countsnorm_df_metadata_shared$Type)
countsnorm_df_metadata_shared$Type <- gsub("Late_captive","Modern",countsnorm_df_metadata_shared$Type)
countsnorm_df_metadata_shared$Type <- gsub("Early_captive","Modern",countsnorm_df_metadata_shared$Type)
countsnorm_df_metadata_shared$Type <- gsub("Wild","Modern",countsnorm_df_metadata_shared$Type)

countsnorm_df_metadata_shared$State <- gsub("HetCountNormLow_5","Heterozygous Load",countsnorm_df_metadata_shared$State)
countsnorm_df_metadata_shared$State <- gsub("HomCountNormLow_5","Homozygous Load",countsnorm_df_metadata_shared$State)
countsnorm_df_metadata_shared$State <- gsub("TotalNormLow_5","Total Load",countsnorm_df_metadata_shared$State)


g_hist <- ggplot(countsnorm_df_metadata_shared[countsnorm_df_metadata_shared$Class!="Low",], 
                 aes(Type, Count, color=Type))
g_hist <- g_hist + geom_boxplot(outlier.shape = NA) + geom_jitter()+theme_classic() + facet_grid(Class~State, scales="free_y",space="free_x") + xlab("") + 
  theme(legend.position = "none") + 
  stat_compare_means() + ylab("Normalized Allele Counts")+
  scale_color_manual(values= c("#6d466bff","#EB9486" ))
g_hist




# Final Figure ----
ggarrange(g_hist, ggarrange(freq_dist, rxy_plot, nrow=2, labels=c("","C")) ,nrow=1, labels=c("A","B"),widths = c(2,1))

ggsave(paste0("Manuscript/MainFigures/Figure2",xxx,".pdf"), height = 6, width = 10)


# Counts with only shared variants NOT INCLUDED------
listfreq_high <- read.csv2("Files/Frequency_derived_high.csv", sep=",")
listfreq_moderate <- read.csv2("Files/Frequency_derived_moderate.csv", sep=",")
listfreq_low <- read.csv2("Files/Frequency_derived_low.csv", sep=",")

listfreq_high_only <- listfreq_high[listfreq_high$FreqHist!=0 & listfreq_high$FreqModern!=0,]$Variant
listfreq_moderate_only <- listfreq_moderate[listfreq_moderate$FreqHist!=0 & listfreq_moderate$FreqModern!=0,]$Variant
listfreq_low_only <- listfreq_low[listfreq_low$FreqHist!=0 & listfreq_low$FreqModern!=0,]$Variant

sites_both <- c(listfreq_high_only, listfreq_moderate_only, listfreq_low_only)
sites_both <- gsub("-",":",sites_both)

# filter from the global count 
high_shared <- high[which(high$chrompos %in% sites_both),]

highcounts<-list()
highcounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- high_shared[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,high_shared[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

highcounts_df_shared <- do.call(rbind, highcounts)

moderate_shared <- moderate[which(moderate$chrompos %in% sites_both),]

moderatecounts<-list()
moderatecounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- moderate_shared[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,moderate_shared[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

moderatecounts_df_shared <- do.call(rbind, moderatecounts)

low_shared <- low[which(low$chrompos %in% sites_both),]

lowcounts<-list()
lowcounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- low_shared[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,low_shared[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

lowcounts_df_shared <- do.call(rbind, lowcounts)

lowcounts_df_shared$Class <- "Low"
moderatecounts_df_shared$Class <- "Moderate"
highcounts_df_shared$Class <- "High"


counts_df_shared <- rbind(highcounts_df_shared,moderatecounts_df_shared, lowcounts_df_shared)
counts_df_shared$Total_5 <-counts_df_shared$het_refanc_5 +counts_df_shared$hom_refanc_5

counts_df_low_shared <- counts_df_shared[counts_df_shared$Class=="Low",]

counts_df_shared$HetCountNormLow_5 <- counts_df_shared$het_refanc_5/counts_df_low_shared$Total_5
counts_df_shared$HomCountNormLow_5 <- counts_df_shared$hom_refanc_5/counts_df_low_shared$Total_5
counts_df_shared$TotalNormLow_5 <- counts_df_shared$Total_5/counts_df_low_shared$Total_5


counts_down_shared <- counts_df_shared[,c("SampleID", "Class",
                                          "het_refanc_5","hom_refanc_5","Total_5",
                                          "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5")]
counts_down_shared$Dataset <- "Downsampled"

#add historical samples counts, now they need to be subset with the sites 
samples<- read.table("Files/snpEff/Samples")
high <- read.table("Files/snpEff/high.txt")
colnames(high) <- c("Chrom","Pos","REF","ALT", samples$V1)
high <- high[,c(1:11,49:52)]

high$chrompos<- paste0(high$Chrom,":",high$Pos)

high_fixed0 <- high %>% filter(grepl("0/0:[0-9]+", Grus27863)) %>% 
  filter(grepl("0/0:[0-9]+", Grus71441)) %>% filter(grepl("0/0:[0-9]+", `MCZ-220028`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-301263`)) %>% filter(grepl("0/0:[0-9]+", `MCZ-35730`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-42511`)) %>%  filter(grepl("0/0:[0-9]+", `MCZ-46873`))

fixed <- high_fixed0$chrompos

high <- high[which(!high$chrompos %in% fixed),]
high_shared <- high[which(high$chrompos %in% sites_both),] # filter sites

moderate <- read.table("Files/snpEff/moderate.txt")
colnames(moderate) <- c("Chrom","Pos","REF","ALT", samples$V1)
moderate <- moderate[,c(1:11,49:52)]

moderate$chrompos<- paste0(moderate$Chrom,":",moderate$Pos)


moderate_fixed0 <- moderate %>% filter(grepl("0/0:[0-9]+", Grus27863)) %>% 
  filter(grepl("0/0:[0-9]+", Grus71441)) %>% filter(grepl("0/0:[0-9]+", `MCZ-220028`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-301263`)) %>% filter(grepl("0/0:[0-9]+", `MCZ-35730`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-42511`)) %>%  filter(grepl("0/0:[0-9]+", `MCZ-46873`))

fixed_moderate <- moderate_fixed0$chrompos
moderate <- moderate[which(!moderate$chrompos %in% fixed_moderate),]
moderate_shared <- moderate[which(moderate$chrompos %in% sites_both),]


low <- read.table("Files/snpEff/low.txt")
colnames(low) <- c("Chrom","Pos","REF","ALT", samples$V1)
low <- low[,c(1:11,49:52)]
low$chrompos<- paste0(low$Chrom,":",low$Pos)

low_fixed0 <- low %>% filter(grepl("0/0:[0-9]+", Grus27863)) %>% 
  filter(grepl("0/0:[0-9]+", Grus71441)) %>% filter(grepl("0/0:[0-9]+", `MCZ-220028`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-301263`)) %>% filter(grepl("0/0:[0-9]+", `MCZ-35730`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-42511`)) %>%  filter(grepl("0/0:[0-9]+", `MCZ-46873`)) 
fixed_low <- low_fixed0$chrompos
low <- low[which(!low$chrompos %in% fixed_low),]
low_shared <- low[which(low$chrompos %in% sites_both),]

# extrat info
highcounts_sharedboth<-list()
highcounts_sharedboth <- lapply(1:7, function(i) {
  j <- i+4
  print(i)
  df_high <- high_shared[,c(j)]
  
  sample_name <- colnames(high_shared)[j]
  
  df_high_2 <- as.data.frame(do.call(rbind,strsplit(df_high, ":")))
  colnames(df_high_2) <- c("Genotype","DP")
  df_high_2$DP <- as.numeric(df_high_2$DP)
  
  df_3_high <- cbind.data.frame(df_high_2,high_shared[,c(1,2,12,13,14,15)])
  df_3_high$Gnigricollis <- gsub(":1","", df_3_high$Gnigricollis)
  df_3_high$birdAnc314 <- gsub(":1","", df_3_high$birdAnc314)
  df_3_high$birdAnc315 <- gsub(":1","", df_3_high$birdAnc315)
  df_3_high$birdAnc316 <- gsub(":1","", df_3_high$birdAnc316)
  
  
  #Filter by DP 5
  df_5 <- df_3_high[df_3_high$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&
                              df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=sample_name,
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})
highcounts_df_sharedboth <- do.call(rbind, highcounts_sharedboth)

moderatecounts_sharedboth<-list()
moderatecounts_sharedboth <- lapply(1:7, function(i) {
  j <- i+4
  print(i)
  df_moderate <- moderate_shared[,c(j)]
  sample_name <- colnames(moderate_shared)[j]
  
  df_moderate_2 <- as.data.frame(do.call(rbind,strsplit(df_moderate, ":")))
  colnames(df_moderate_2) <- c("Genotype","DP")
  df_moderate_2$DP <- as.numeric(df_moderate_2$DP)
  
  df_3_mod <- cbind.data.frame(df_moderate_2,moderate_shared[,c(1,2,12,13,14,15)])
  df_3_mod$Gnigricollis <- gsub(":1","", df_3_mod$Gnigricollis)
  df_3_mod$birdAnc314 <- gsub(":1","", df_3_mod$birdAnc314)
  df_3_mod$birdAnc315 <- gsub(":1","", df_3_mod$birdAnc315)
  df_3_mod$birdAnc316 <- gsub(":1","", df_3_mod$birdAnc316)
  
  
  #Filter by DP 5
  df_5 <- df_3_mod[df_3_mod$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&
                              df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=sample_name,
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})
moderatecounts_df_sharedboth <- do.call(rbind, moderatecounts_sharedboth)

lowcounts_sharedboth<-list()
lowcounts_sharedboth <- lapply(1:7, function(i) {
  j <- i+4
  print(i)
  df_low <- low_shared[,c(j)]
  
  sample_name <- colnames(low_shared)[j]
  
  df_low_2 <- as.data.frame(do.call(rbind,strsplit(df_low, ":")))
  colnames(df_low_2) <- c("Genotype","DP")
  df_low_2$DP <- as.numeric(df_low_2$DP)
  
  df_3_low <- cbind.data.frame(df_low_2,low_shared[,c(1,2,12,13,14,15)])
  df_3_low$Gnigricollis <- gsub(":1","", df_3_low$Gnigricollis)
  df_3_low$birdAnc314 <- gsub(":1","", df_3_low$birdAnc314)
  df_3_low$birdAnc315 <- gsub(":1","", df_3_low$birdAnc315)
  df_3_low$birdAnc316 <- gsub(":1","", df_3_low$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3_low[df_3_low$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&
                              df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=sample_name,
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})
lowcounts_df_sharedboth <- do.call(rbind, lowcounts_sharedboth)

lowcounts_df_sharedboth$Class <- "Low"
moderatecounts_df_sharedboth$Class <- "Moderate"
highcounts_df_sharedboth$Class <- "High"

counts_df_sharedboth <- rbind(highcounts_df_sharedboth,moderatecounts_df_sharedboth, lowcounts_df_sharedboth)
counts_df_sharedboth$Total_5 <-counts_df_sharedboth$het_refanc_5 +counts_df_sharedboth$hom_refanc_5

counts_df_low_sharedboth <- counts_df_sharedboth[counts_df_sharedboth$Class=="Low",]

counts_df_sharedboth$HetCountNormLow_5 <- counts_df_sharedboth$het_refanc_5/counts_df_low_sharedboth$Total_5
counts_df_sharedboth$HomCountNormLow_5 <- counts_df_sharedboth$hom_refanc_5/counts_df_low_sharedboth$Total_5
counts_df_sharedboth$TotalNormLow_5 <- counts_df_sharedboth$Total_5/counts_df_low_sharedboth$Total_5


counts_down_sharedbothhist <- counts_df_sharedboth[,c("SampleID", "Class",
                                                      "het_refanc_5","hom_refanc_5","Total_5",
                                                      "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5")]
counts_down_sharedbothhist$Dataset <- "Alldata"


# add everything together
total_counts_shared <- rbind.data.frame(counts_down_sharedbothhist,counts_down_shared)

total_counts_down_shared <-total_counts_shared[-which(total_counts_shared$SampleID%in%samples_cov[-c(3,6)]),]


total_counts_df3_shared <- total_counts_down_shared[,c("SampleID", "Class",
                                                       "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5","Dataset")]


counts_df_gather_shared <- gather(total_counts_df3_shared, State, Count, -SampleID,  -Class, -Dataset,)

countsnorm_df_metadata_shared <- merge(counts_df_gather_shared, metadata[,c("SampleID","Category")], by="SampleID") 
countsnorm_df_metadata_shared$Class <- factor(countsnorm_df_metadata_shared$Class, levels=c("Low","Moderate","High"), ordered = T)

countsnorm_df_metadata_shared$Type <- countsnorm_df_metadata_shared$Category
countsnorm_df_metadata_shared$Type <- gsub("Founders_wildborn","Modern",countsnorm_df_metadata_shared$Type)
countsnorm_df_metadata_shared$Type <- gsub("Late_captive","Modern",countsnorm_df_metadata_shared$Type)
countsnorm_df_metadata_shared$Type <- gsub("Early_captive","Modern",countsnorm_df_metadata_shared$Type)
countsnorm_df_metadata_shared$Type <- gsub("Wild","Modern",countsnorm_df_metadata_shared$Type)

countsnorm_df_metadata_shared$State <- gsub("HetCountNormLow_5","Heterozygous Load",countsnorm_df_metadata_shared$State)
countsnorm_df_metadata_shared$State <- gsub("HomCountNormLow_5","Homozygous Load",countsnorm_df_metadata_shared$State)
countsnorm_df_metadata_shared$State <- gsub("TotalNormLow_5","Total Load",countsnorm_df_metadata_shared$State)


g_hist <- ggplot(countsnorm_df_metadata_shared[countsnorm_df_metadata_shared$Class!="Low",], 
                 aes(Type, Count, color=Type))
g_hist <- g_hist + geom_boxplot(outlier.shape = NA) + geom_jitter()+theme_classic() + facet_grid(Class~State, scales="free_y",space="free_x") + xlab("") + 
  theme(legend.position = "none") + 
  stat_compare_means() + ylab("Normalized Allele Counts")+
  scale_color_manual(values= c("#6d466bff","#EB9486" ))
g_hist



# Final Figure 
ggarrange(g_hist, ggarrange(freq_dist, rxy_plot, nrow=2, labels=c("","C")) ,nrow=1, labels=c("A","B"),widths = c(2,1))

ggsave(paste0("Manuscript/MainFigures/Figure2",xxx,".pdf"), height = 6, width = 10)


# Keeping sites non 0 frequency in modern ------
#Change in Frequency of deleterious variation through time-----

listfreq_high <- read.csv2("Files/Frequency_derived_high.csv", sep=",")
listfreq_moderate <- read.csv2("Files/Frequency_derived_moderate.csv", sep=",")
listfreq_low <- read.csv2("Files/Frequency_derived_low.csv", sep=",")

listfreq_high$Type <- "High"
listfreq_moderate$Type <- "Moderate"
listfreq_low$Type <- "Low"


# only sites with modern freq >0
listfreq3 <- rbind(listfreq_high,listfreq_moderate,listfreq_low)
listfreq3<- listfreq3[which(listfreq3$FreqModern>0),]  # only sites modern >0 freq

listfreq3$Type <- factor(listfreq3$Type, levels=c("Low","Moderate","High"), ordered = T)
listfreq3$FreqHist <- as.numeric(listfreq3$FreqHist)
listfreq3$FreqModern <- as.numeric(listfreq3$FreqModern)

listfreq3$Difference <- listfreq3$FreqModern - listfreq3$FreqHist

# Plot distribution
freq_dist <- ggplot(listfreq3, aes(Difference, color=Type))
freq_dist <- freq_dist + geom_density() + theme_classic() + 
  scale_color_manual(values=c("#4f7cac","#7dde92","#696047"))+ xlab("Delta Frequency (Modern - Historical)")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")+ theme(legend.position = "top", legend.title = element_blank())

freq_dist

#percentage with higher frequency in modern
length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="High",]$Variant)/length(listfreq3[listfreq3$Type=="High",]$Variant)
length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="Moderate",]$Variant)/length(listfreq3[listfreq3$Type=="Moderate",]$Variant)
length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="Low",]$Variant)/length(listfreq3[listfreq3$Type=="Low",]$Variant)


mean(listfreq3[listfreq3$Type=="High",]$Difference)
sd(listfreq3[listfreq3$Type=="High",]$Difference)
mean(listfreq3[listfreq3$Type=="Moderate",]$Difference)
mean(listfreq3[listfreq3$Type=="Low",]$Difference)

# Rxy ------
#read files
Frequencies_all <- listfreq3

##Rxy
types <- c("High","Moderate", "Low")

Rxy_frame<-data.frame(matrix(ncol = 4, nrow = 0))
Rxy_jn<-data.frame(matrix(ncol = 3, nrow = 0))

y<-Frequencies_all$FreqHist*(1-Frequencies_all$FreqModern)
x<-Frequencies_all$FreqModern*(1-Frequencies_all$FreqHist)
xdata_total<-cbind.data.frame(x,y, Type=Frequencies_all$Type)

for (j in types){
  print(j)
  xdata <- xdata_total[xdata_total$Type==j,]
  
  # size <- floor(length(xdata$x)*0.01) # take 1% of the sites for the jacknife
  size <- ceiling(length(xdata$x)*0.01) # take 1% of the sites for the jacknife ONLY FOR THE SHARED FILES
  rxy<-function(a,xdata,b){
    sum(xdata[seq(a,a+(size-1),1),b])
  } 
  rx_ry<-data.frame()
  for(i in c(1:100)) {
    rx_ry[i,1]<-rxy(size*(i-1)+1,xdata,1)
    rx_ry[i,2]<-rxy(size*(i-1)+1,xdata,2)
  } 
  
  rat<-function(c,rx_ry){sum(rx_ry[c,1])/sum(rx_ry[c,2])}
  jn_h<-jackknife(1:100,rat,rx_ry) #do jackknife for dataframe rx_ry with the function of rat for 100 times
  Rxy_jn<-rbind(Rxy_jn, cbind(jn_h$jack.values, rep(j, 100)))
  Rxy_frame<-rbind(Rxy_frame, cbind(j, sum(xdata$x)/sum(xdata$y), jn_h$jack.se))
} 

colnames(Rxy_frame)<-c("Type","Rxy", "SE")
Rxy_frame$Rxy <- as.numeric(Rxy_frame$Rxy)
Rxy_frame$SE <- as.numeric(Rxy_frame$SE)

# Normalize by low 
Rxy_frame_low <- Rxy_frame[Rxy_frame$Type=="Low",]
Rxy_frame_notlow <- Rxy_frame[Rxy_frame$Type!="Low",]

Rxy_frame_normalized<- data.frame(Type=c("High","Moderate"),
                                  Rxy_norm=as.numeric(Rxy_frame_notlow$Rxy)/as.numeric(Rxy_frame_low$Rxy),
                                  SE_norm=as.numeric(Rxy_frame_notlow$S)/as.numeric(Rxy_frame_low$Rxy))

cut.off = 1
Rxy_frame_normalized$Rxy_norm_2 <- Rxy_frame_normalized$Rxy_norm -1

rxy_plot <- ggplot(Rxy_frame_normalized, aes(Type, Rxy_norm_2, fill=Type))
rxy_plot <- rxy_plot + geom_col(position = "dodge") + theme_classic()+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) +scale_fill_manual(values=c("#696047","#7dde92")) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + ylab("Rxy (Modern/Historical)") + theme(legend.position = "none")
rxy_plot

# counts -----
listfreq_high <- read.csv2("Files/Frequency_derived_high.csv", sep=",")
listfreq_moderate <- read.csv2("Files/Frequency_derived_moderate.csv", sep=",")
listfreq_low <- read.csv2("Files/Frequency_derived_low.csv", sep=",")

listfreq_high_only <- listfreq_high[listfreq_high$FreqModern>0,]$Variant
listfreq_moderate_only <- listfreq_moderate[listfreq_moderate$FreqModern>0,]$Variant
listfreq_low_only <- listfreq_low[listfreq_low$FreqModern>0,]$Variant

sites_both <- c(listfreq_high_only, listfreq_moderate_only, listfreq_low_only)
sites_both <- gsub("-",":",sites_both)

# filter from the global count 
high_shared <- high[which(high$chrompos %in% sites_both),]

highcounts<-list()
highcounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- high_shared[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,high_shared[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

highcounts_df_shared <- do.call(rbind, highcounts)

moderate_shared <- moderate[which(moderate$chrompos %in% sites_both),]

moderatecounts<-list()
moderatecounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- moderate_shared[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,moderate_shared[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

moderatecounts_df_shared <- do.call(rbind, moderatecounts)

low_shared <- low[which(low$chrompos %in% sites_both),]

lowcounts<-list()
lowcounts <- lapply(1:41, function(i) {
  j <- i+4
  df <- low_shared[,c(j)]
  df_2 <- as.data.frame(do.call(rbind,strsplit(df, ":")))
  colnames(df_2) <- c("Genotype","DP")
  df_2$DP <- as.numeric(df_2$DP)
  
  df_3 <- cbind.data.frame(df_2,low_shared[,c(46,47,48,49)])
  df_3$Gnigricollis <- gsub(":1","", df_3$Gnigricollis)
  df_3$birdAnc314 <- gsub(":1","", df_3$birdAnc314)
  df_3$birdAnc315 <- gsub(":1","", df_3$birdAnc315)
  df_3$birdAnc316 <- gsub(":1","", df_3$birdAnc316)
  
  
  #Filter by DP 5
  df_5 <- df_3[df_3$DP >=5,]
  
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=samples$V1[i],
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})

lowcounts_df_shared <- do.call(rbind, lowcounts)

lowcounts_df_shared$Class <- "Low"
moderatecounts_df_shared$Class <- "Moderate"
highcounts_df_shared$Class <- "High"


counts_df_shared <- rbind(highcounts_df_shared,moderatecounts_df_shared, lowcounts_df_shared)
counts_df_shared$Total_5 <-counts_df_shared$het_refanc_5 +counts_df_shared$hom_refanc_5

counts_df_low_shared <- counts_df_shared[counts_df_shared$Class=="Low",]

counts_df_shared$HetCountNormLow_5 <- counts_df_shared$het_refanc_5/counts_df_low_shared$Total_5
counts_df_shared$HomCountNormLow_5 <- counts_df_shared$hom_refanc_5/counts_df_low_shared$Total_5
counts_df_shared$TotalNormLow_5 <- counts_df_shared$Total_5/counts_df_low_shared$Total_5


counts_down_shared <- counts_df_shared[,c("SampleID", "Class",
                                          "het_refanc_5","hom_refanc_5","Total_5",
                                          "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5")]
counts_down_shared$Dataset <- "Downsampled"

#add historical samples counts, now they need to be subset with the sites 
samples<- read.table("Files/snpEff/Samples")
high <- read.table("Files/snpEff/high.txt")
colnames(high) <- c("Chrom","Pos","REF","ALT", samples$V1)
high <- high[,c(1:11,49:52)]

high$chrompos<- paste0(high$Chrom,":",high$Pos)

high_fixed0 <- high %>% filter(grepl("0/0:[0-9]+", Grus27863)) %>% 
  filter(grepl("0/0:[0-9]+", Grus71441)) %>% filter(grepl("0/0:[0-9]+", `MCZ-220028`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-301263`)) %>% filter(grepl("0/0:[0-9]+", `MCZ-35730`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-42511`)) %>%  filter(grepl("0/0:[0-9]+", `MCZ-46873`))

fixed <- high_fixed0$chrompos

high <- high[which(!high$chrompos %in% fixed),]
high_shared <- high[which(high$chrompos %in% sites_both),] # filter sites

moderate <- read.table("Files/snpEff/moderate.txt")
colnames(moderate) <- c("Chrom","Pos","REF","ALT", samples$V1)
moderate <- moderate[,c(1:11,49:52)]

moderate$chrompos<- paste0(moderate$Chrom,":",moderate$Pos)


moderate_fixed0 <- moderate %>% filter(grepl("0/0:[0-9]+", Grus27863)) %>% 
  filter(grepl("0/0:[0-9]+", Grus71441)) %>% filter(grepl("0/0:[0-9]+", `MCZ-220028`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-301263`)) %>% filter(grepl("0/0:[0-9]+", `MCZ-35730`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-42511`)) %>%  filter(grepl("0/0:[0-9]+", `MCZ-46873`))

fixed_moderate <- moderate_fixed0$chrompos
moderate <- moderate[which(!moderate$chrompos %in% fixed_moderate),]
moderate_shared <- moderate[which(moderate$chrompos %in% sites_both),]


low <- read.table("Files/snpEff/low.txt")
colnames(low) <- c("Chrom","Pos","REF","ALT", samples$V1)
low <- low[,c(1:11,49:52)]
low$chrompos<- paste0(low$Chrom,":",low$Pos)

low_fixed0 <- low %>% filter(grepl("0/0:[0-9]+", Grus27863)) %>% 
  filter(grepl("0/0:[0-9]+", Grus71441)) %>% filter(grepl("0/0:[0-9]+", `MCZ-220028`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-301263`)) %>% filter(grepl("0/0:[0-9]+", `MCZ-35730`)) %>% 
  filter(grepl("0/0:[0-9]+", `MCZ-42511`)) %>%  filter(grepl("0/0:[0-9]+", `MCZ-46873`)) 
fixed_low <- low_fixed0$chrompos
low <- low[which(!low$chrompos %in% fixed_low),]
low_shared <- low[which(low$chrompos %in% sites_both),]

# extrat info
highcounts_sharedboth<-list()
highcounts_sharedboth <- lapply(1:7, function(i) {
  j <- i+4
  print(i)
  df_high <- high_shared[,c(j)]
  
  sample_name <- colnames(high_shared)[j]
  
  df_high_2 <- as.data.frame(do.call(rbind,strsplit(df_high, ":")))
  colnames(df_high_2) <- c("Genotype","DP")
  df_high_2$DP <- as.numeric(df_high_2$DP)
  
  df_3_high <- cbind.data.frame(df_high_2,high_shared[,c(1,2,12,13,14,15)])
  df_3_high$Gnigricollis <- gsub(":1","", df_3_high$Gnigricollis)
  df_3_high$birdAnc314 <- gsub(":1","", df_3_high$birdAnc314)
  df_3_high$birdAnc315 <- gsub(":1","", df_3_high$birdAnc315)
  df_3_high$birdAnc316 <- gsub(":1","", df_3_high$birdAnc316)
  
  
  #Filter by DP 5
  df_5 <- df_3_high[df_3_high$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&
                              df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=sample_name,
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})
highcounts_df_sharedboth <- do.call(rbind, highcounts_sharedboth)

moderatecounts_sharedboth<-list()
moderatecounts_sharedboth <- lapply(1:7, function(i) {
  j <- i+4
  print(i)
  df_moderate <- moderate_shared[,c(j)]
  sample_name <- colnames(moderate_shared)[j]
  
  df_moderate_2 <- as.data.frame(do.call(rbind,strsplit(df_moderate, ":")))
  colnames(df_moderate_2) <- c("Genotype","DP")
  df_moderate_2$DP <- as.numeric(df_moderate_2$DP)
  
  df_3_mod <- cbind.data.frame(df_moderate_2,moderate_shared[,c(1,2,12,13,14,15)])
  df_3_mod$Gnigricollis <- gsub(":1","", df_3_mod$Gnigricollis)
  df_3_mod$birdAnc314 <- gsub(":1","", df_3_mod$birdAnc314)
  df_3_mod$birdAnc315 <- gsub(":1","", df_3_mod$birdAnc315)
  df_3_mod$birdAnc316 <- gsub(":1","", df_3_mod$birdAnc316)
  
  
  #Filter by DP 5
  df_5 <- df_3_mod[df_3_mod$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&
                              df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=sample_name,
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})
moderatecounts_df_sharedboth <- do.call(rbind, moderatecounts_sharedboth)

lowcounts_sharedboth<-list()
lowcounts_sharedboth <- lapply(1:7, function(i) {
  j <- i+4
  print(i)
  df_low <- low_shared[,c(j)]
  
  sample_name <- colnames(low_shared)[j]
  
  df_low_2 <- as.data.frame(do.call(rbind,strsplit(df_low, ":")))
  colnames(df_low_2) <- c("Genotype","DP")
  df_low_2$DP <- as.numeric(df_low_2$DP)
  
  df_3_low <- cbind.data.frame(df_low_2,low_shared[,c(1,2,12,13,14,15)])
  df_3_low$Gnigricollis <- gsub(":1","", df_3_low$Gnigricollis)
  df_3_low$birdAnc314 <- gsub(":1","", df_3_low$birdAnc314)
  df_3_low$birdAnc315 <- gsub(":1","", df_3_low$birdAnc315)
  df_3_low$birdAnc316 <- gsub(":1","", df_3_low$birdAnc316)
  
  #Filter by DP 5
  df_5 <- df_3_low[df_3_low$DP >=5,]
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&
                              df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  
  data.frame(SampleID=sample_name,
             het_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="0/1",]$Genotype),
             hom_refanc_5=length(df_refanc_5[df_refanc_5$Genotype=="1/1",]$Genotype)*2)
})
lowcounts_df_sharedboth <- do.call(rbind, lowcounts_sharedboth)

lowcounts_df_sharedboth$Class <- "Low"
moderatecounts_df_sharedboth$Class <- "Moderate"
highcounts_df_sharedboth$Class <- "High"

counts_df_sharedboth <- rbind(highcounts_df_sharedboth,moderatecounts_df_sharedboth, lowcounts_df_sharedboth)
counts_df_sharedboth$Total_5 <-counts_df_sharedboth$het_refanc_5 +counts_df_sharedboth$hom_refanc_5

counts_df_low_sharedboth <- counts_df_sharedboth[counts_df_sharedboth$Class=="Low",]

counts_df_sharedboth$HetCountNormLow_5 <- counts_df_sharedboth$het_refanc_5/counts_df_low_sharedboth$Total_5
counts_df_sharedboth$HomCountNormLow_5 <- counts_df_sharedboth$hom_refanc_5/counts_df_low_sharedboth$Total_5
counts_df_sharedboth$TotalNormLow_5 <- counts_df_sharedboth$Total_5/counts_df_low_sharedboth$Total_5


counts_down_sharedbothhist <- counts_df_sharedboth[,c("SampleID", "Class",
                                                      "het_refanc_5","hom_refanc_5","Total_5",
                                                      "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5")]
counts_down_sharedbothhist$Dataset <- "Alldata"


# add everything together
total_counts_shared <- rbind.data.frame(counts_down_sharedbothhist,counts_down_shared)

total_counts_down_shared <-total_counts_shared[-which(total_counts_shared$SampleID%in%samples_cov[-c(3,6)]),]


total_counts_df3_shared <- total_counts_down_shared[,c("SampleID", "Class",
                                                       "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5","Dataset")]


counts_df_gather_shared <- gather(total_counts_df3_shared, State, Count, -SampleID,  -Class, -Dataset,)

countsnorm_df_metadata_shared <- merge(counts_df_gather_shared, metadata[,c("SampleID","Category")], by="SampleID") 
countsnorm_df_metadata_shared$Class <- factor(countsnorm_df_metadata_shared$Class, levels=c("Low","Moderate","High"), ordered = T)

countsnorm_df_metadata_shared$Type <- countsnorm_df_metadata_shared$Category
countsnorm_df_metadata_shared$Type <- gsub("Founders_wildborn","Modern",countsnorm_df_metadata_shared$Type)
countsnorm_df_metadata_shared$Type <- gsub("Late_captive","Modern",countsnorm_df_metadata_shared$Type)
countsnorm_df_metadata_shared$Type <- gsub("Early_captive","Modern",countsnorm_df_metadata_shared$Type)
countsnorm_df_metadata_shared$Type <- gsub("Wild","Modern",countsnorm_df_metadata_shared$Type)

countsnorm_df_metadata_shared$State <- gsub("HetCountNormLow_5","Heterozygous Load",countsnorm_df_metadata_shared$State)
countsnorm_df_metadata_shared$State <- gsub("HomCountNormLow_5","Homozygous Load",countsnorm_df_metadata_shared$State)
countsnorm_df_metadata_shared$State <- gsub("TotalNormLow_5","Total Load",countsnorm_df_metadata_shared$State)


g_hist <- ggplot(countsnorm_df_metadata_shared[countsnorm_df_metadata_shared$Class!="Low",], 
                 aes(Type, Count, color=Type))
g_hist <- g_hist + geom_boxplot(outlier.shape = NA) + geom_jitter()+theme_classic() + facet_grid(Class~State, scales="free_y",space="free_x") + xlab("") + 
  theme(legend.position = "none") + 
  stat_compare_means() + ylab("Normalized Allele Counts")+
  scale_color_manual(values= c("#6d466bff","#EB9486" ))
g_hist



# Final Figure ----
ggarrange(g_hist, ggarrange(freq_dist, rxy_plot, nrow=2, labels=c("","C")) ,nrow=1, labels=c("A","B"),widths = c(2,1))

ggsave(paste0("Manuscript/MainFigures/Figure2","modernmore0",".pdf"), height = 6, width = 10)
