# Change of frequency from historical to modern of high, moderate and low load
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)


setwd("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/")

# Read data from historical samples and format equally as in ROHs_load_down10.R
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
# extrat info
historical_list<-list()
historical_list <- lapply(1:7, function(i) {
  j <- i+4
  print(i)
  df_high <- high[,c(j)]
  df_moderate <- moderate[,c(j)]
  df_low <- low[,c(j)]
  
  sample_name <- colnames(high)[j]
  
  df_high_2 <- as.data.frame(do.call(rbind,strsplit(df_high, ":")))
  colnames(df_high_2) <- c("Genotype","DP")
  df_high_2$DP <- as.numeric(df_high_2$DP)
  
  df_moderate_2 <- as.data.frame(do.call(rbind,strsplit(df_moderate, ":")))
  colnames(df_moderate_2) <- c("Genotype","DP")
  df_moderate_2$DP <- as.numeric(df_moderate_2$DP)
  
  df_low_2 <- as.data.frame(do.call(rbind,strsplit(df_low, ":")))
  colnames(df_low_2) <- c("Genotype","DP")
  df_low_2$DP <- as.numeric(df_low_2$DP)
  
  
  df_3_high <- cbind.data.frame(df_high_2,high[,c(1,2,12,13,14,15)])
  df_3_high$Gnigricollis <- gsub(":1","", df_3_high$Gnigricollis)
  df_3_high$birdAnc314 <- gsub(":1","", df_3_high$birdAnc314)
  df_3_high$birdAnc315 <- gsub(":1","", df_3_high$birdAnc315)
  df_3_high$birdAnc316 <- gsub(":1","", df_3_high$birdAnc316)
  
  df_3_mod <- cbind.data.frame(df_moderate_2,moderate[,c(1,2,12,13,14,15)])
  df_3_mod$Gnigricollis <- gsub(":1","", df_3_mod$Gnigricollis)
  df_3_mod$birdAnc314 <- gsub(":1","", df_3_mod$birdAnc314)
  df_3_mod$birdAnc315 <- gsub(":1","", df_3_mod$birdAnc315)
  df_3_mod$birdAnc316 <- gsub(":1","", df_3_mod$birdAnc316)
  
  df_3_low <- cbind.data.frame(df_low_2,low[,c(1,2,12,13,14,15)])
  df_3_low$Gnigricollis <- gsub(":1","", df_3_low$Gnigricollis)
  df_3_low$birdAnc314 <- gsub(":1","", df_3_low$birdAnc314)
  df_3_low$birdAnc315 <- gsub(":1","", df_3_low$birdAnc315)
  df_3_low$birdAnc316 <- gsub(":1","", df_3_low$birdAnc316)
  
  
  #Filter by DP 5
  df_5_high <- df_3_high[df_3_high$DP >=5,]
  df_5_high$Class <- "High"
  df_5_moderate <- df_3_mod[df_3_mod$DP >=5,]
  df_5_moderate$Class <- "Moderate"
  df_5_low <- df_3_low[df_3_low$DP >=5,]
  df_5_low$Class <- "Low"
  
  df_5 <- rbind.data.frame(df_5_high,df_5_moderate, df_5_low)
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  cbind.data.frame(df_refanc_5[,c(1:4,9)], SampleID=sample_name)
})

historical_genotypes <- do.call(rbind,historical_list)

# Load data from downsampled 10x only modern 
down_10x_genotypes <- read.csv2("Files/Genotype_load_modern_10x.csv")


# Calculate frequency of variants in historical and modern ----
historical_genotypes$Type <- "Historical"
down_10x_genotypes$Type <- "Modern"


samples_cov <- c("EB31_S103_15x","EB31_S103_5x", "EB31_S103",
                 "EB33_S105_15x","EB33_S105_5x","EB33_S105")

down_10x_genotypes <-down_10x_genotypes[-which(down_10x_genotypes$SampleID%in%samples_cov),]


total_genotypes <- rbind.data.frame(historical_genotypes, down_10x_genotypes)
total_genotypes$Variant <- paste0(total_genotypes$Chrom,"-",total_genotypes$Pos)

# All categories and not 0/0
derived_genotypes_high <- total_genotypes[total_genotypes$Genotype!="0/0"&total_genotypes$Class=="High",]
unique_variants_high <- unique(derived_genotypes_high$Variant)

derived_genotypes_moderate <- total_genotypes[total_genotypes$Genotype!="0/0"&total_genotypes$Class=="Moderate",]
unique_variants_moderate <- unique(derived_genotypes_moderate$Variant)

derived_genotypes_low <- total_genotypes[total_genotypes$Genotype!="0/0"&total_genotypes$Class=="Low",]
unique_variants_low <- unique(derived_genotypes_low$Variant)

#loop over variants

#high
listfreq_high <- do.call(rbind,lapply(1:length(unique_variants_high), function(i) {
  variant <- unique_variants_high[i]
  hist <- derived_genotypes_high[derived_genotypes_high$Type=="Historical"&derived_genotypes_high$Variant==variant,]
  f_hist<-(length(hist[hist$Genotype =="0/1",]$Genotype) + length(hist[hist$Genotype =="1/1",]$Genotype)*2)/14
  
  modern <- derived_genotypes_high[derived_genotypes_high$Type=="Modern"&derived_genotypes_high$Variant==variant,]
  f_modern<-(length(modern[modern$Genotype =="0/1",]$Genotype) + length(modern[modern$Genotype =="1/1",]$Genotype)*2)/70
  cbind.data.frame(Variant=variant,FreqHist=f_hist,FreqModern=f_modern)
}))

#moderate
listfreq_moderate <- do.call(rbind,lapply(1:length(unique_variants_moderate), function(i) {
  variant <- unique_variants_moderate[i]
  hist <- derived_genotypes_moderate[derived_genotypes_moderate$Type=="Historical"&derived_genotypes_moderate$Variant==variant,]
  f_hist<-(length(hist[hist$Genotype =="0/1",]$Genotype) + length(hist[hist$Genotype =="1/1",]$Genotype)*2)/14
  
  modern <- derived_genotypes_moderate[derived_genotypes_moderate$Type=="Modern"&derived_genotypes_moderate$Variant==variant,]
  f_modern<-(length(modern[modern$Genotype =="0/1",]$Genotype) + length(modern[modern$Genotype =="1/1",]$Genotype)*2)/70
  cbind.data.frame(Variant=variant,FreqHist=f_hist,FreqModern=f_modern)
}))

#low
listfreq_low <- do.call(rbind,lapply(1:length(unique_variants_low), function(i) {
  variant <- unique_variants_low[i]
  hist <- derived_genotypes_low[derived_genotypes_low$Type=="Historical"&derived_genotypes_low$Variant==variant,]
  f_hist<-(length(hist[hist$Genotype =="0/1",]$Genotype) + length(hist[hist$Genotype =="1/1",]$Genotype)*2)/14
  
  modern <- derived_genotypes_low[derived_genotypes_low$Type=="Modern"&derived_genotypes_low$Variant==variant,]
  f_modern<-(length(modern[modern$Genotype =="0/1",]$Genotype) + length(modern[modern$Genotype =="1/1",]$Genotype)*2)/70
  cbind.data.frame(Variant=variant,FreqHist=f_hist,FreqModern=f_modern)
}))

# Save frequencies
write.csv(listfreq_high, "Files/Frequency_derived_high.csv", quote=F, row.names = F)
write.csv(listfreq_moderate, "Files/Frequency_derived_moderate.csv", quote=F, row.names = F)
write.csv(listfreq_low, "Files/Frequency_derived_low.csv", quote=F, row.names = F)

#high but now by category ----

# All categories and not 0/0 only for modern
founders_samples <- c("EB25_S97","EB28_S100","EB30_S102","EB31_S103","EB32_S104","EB35_S107")
captive_samples <-c("EB24_S96","EB26_S98","EB27_S99","EB29_S101","EB33_S105","EB34_S106")
captive_late_samples <- c("EB37_S108","EB38_S109","EB39_S110","EB40_S111","EB41_S112","EB42_S113")

wild_samples <- c("EB11_S84","EB12_S85","EB13_S86","EB14_S87","EB15_S88","EB16_S89",
                  "EB17_S90","EB18_S91","EB19_S92","EB1_S114","EB20_S93","EB21_S94",
                  "EB22_S95","EB2_S115","EB3_S116","EB6_S117","EB7_S118","EB8_S119","EB9_S83")


derived_genotypes_high_modern <- total_genotypes[total_genotypes$Type!="Historical"&total_genotypes$Genotype!="0/0"&total_genotypes$Class=="High",]
unique_variants_high_modern <- unique(derived_genotypes_high_modern$Variant)

derived_genotypes_moderate_modern <- total_genotypes[total_genotypes$Type!="Historical"&total_genotypes$Genotype!="0/0"&total_genotypes$Class=="Moderate",]
unique_variants_moderate_modern <- unique(derived_genotypes_moderate_modern$Variant)

derived_genotypes_low_modern <- total_genotypes[total_genotypes$Type!="Historical"&total_genotypes$Genotype!="0/0"&total_genotypes$Class=="Low",]
unique_variants_low_modern <- unique(derived_genotypes_low_modern$Variant)

# with hist
derived_genotypes_high_modern <- total_genotypes[total_genotypes$Genotype!="0/0"&total_genotypes$Class=="High",]
unique_variants_high_modern <- unique(derived_genotypes_high_modern$Variant)

derived_genotypes_moderate_modern <- total_genotypes[total_genotypes$Genotype!="0/0"&total_genotypes$Class=="Moderate",]
unique_variants_moderate_modern <- unique(derived_genotypes_moderate_modern$Variant)

derived_genotypes_low_modern <- total_genotypes[total_genotypes$Genotype!="0/0"&total_genotypes$Class=="Low",]
unique_variants_low_modern <- unique(derived_genotypes_low_modern$Variant)


listfreq_high_cat <- do.call(rbind,lapply(1:length(unique_variants_high_modern), function(i) {
  variant <- unique_variants_high_modern[i]
  hist <- derived_genotypes_high_modern[derived_genotypes_high_modern$Type=="Historical"&derived_genotypes_high_modern$Variant==variant,]
  f_hist<-(length(hist[hist$Genotype =="0/1",]$Genotype) + length(hist[hist$Genotype =="1/1",]$Genotype)*2)/14
  
  founders <- derived_genotypes_high_modern[which(derived_genotypes_high_modern$SampleID%in%founders_samples&
                                                    derived_genotypes_high_modern$Variant==variant),]
  f_founders<-(length(founders[founders$Genotype =="0/1",]$Genotype) + length(founders[founders$Genotype =="1/1",]$Genotype)*2)/12
 
  captive1 <- derived_genotypes_high_modern[which(derived_genotypes_high_modern$SampleID%in%captive_samples&
                                                    derived_genotypes_high_modern$Variant==variant),]
  f_captive1<-(length(captive1[captive1$Genotype =="0/1",]$Genotype) + length(captive1[captive1$Genotype =="1/1",]$Genotype)*2)/12
  
  captive2 <- derived_genotypes_high_modern[which(derived_genotypes_high_modern$SampleID%in%captive_late_samples&
                                                    derived_genotypes_high_modern$Variant==variant),]
  f_captive2<-(length(captive2[captive2$Genotype =="0/1",]$Genotype) + length(captive2[captive2$Genotype =="1/1",]$Genotype)*2)/12
  
  wild <- derived_genotypes_high_modern[which(derived_genotypes_high_modern$SampleID%in%wild_samples&
                                                derived_genotypes_high_modern$Variant==variant),]
  f_wild<-(length(wild[wild$Genotype =="0/1",]$Genotype) + length(wild[wild$Genotype =="1/1",]$Genotype)*2)/38

   cbind.data.frame(Variant=variant,FreqHist=f_hist,FreqFounders=f_founders,FreqEarlyCaptive=f_captive1,
                    FreqLateCaptive=f_captive2,FreqWild=f_wild)
}))

#moderate but now by category 
listfreq_moderate_cat <- do.call(rbind,lapply(1:length(unique_variants_moderate_modern), function(i) {
  variant <- unique_variants_moderate_modern[i]
  hist <- derived_genotypes_moderate_modern[derived_genotypes_moderate_modern$Type=="Historical"&derived_genotypes_moderate_modern$Variant==variant,]
  f_hist<-(length(hist[hist$Genotype =="0/1",]$Genotype) + length(hist[hist$Genotype =="1/1",]$Genotype)*2)/14
 
   founders <- derived_genotypes_moderate_modern[which(derived_genotypes_moderate_modern$SampleID%in%founders_samples&
                                                         derived_genotypes_moderate_modern$Variant==variant),]
  f_founders<-(length(founders[founders$Genotype =="0/1",]$Genotype) + length(founders[founders$Genotype =="1/1",]$Genotype)*2)/12
  
  captive1 <- derived_genotypes_moderate_modern[which(derived_genotypes_moderate_modern$SampleID%in%captive_samples&
                                                        derived_genotypes_moderate_modern$Variant==variant),]
  f_captive1<-(length(captive1[captive1$Genotype =="0/1",]$Genotype) + length(captive1[captive1$Genotype =="1/1",]$Genotype)*2)/12
  
  captive2 <- derived_genotypes_moderate_modern[which(derived_genotypes_moderate_modern$SampleID%in%captive_late_samples&
                                                        derived_genotypes_moderate_modern$Variant==variant),]
  f_captive2<-(length(captive2[captive2$Genotype =="0/1",]$Genotype) + length(captive2[captive2$Genotype =="1/1",]$Genotype)*2)/12
  
  wild <- derived_genotypes_moderate_modern[which(derived_genotypes_moderate_modern$SampleID%in%wild_samples&
                                                    derived_genotypes_moderate_modern$Variant==variant),]
  f_wild<-(length(wild[wild$Genotype =="0/1",]$Genotype) + length(wild[wild$Genotype =="1/1",]$Genotype)*2)/38
  
  cbind.data.frame(Variant=variant,FreqHist=f_hist,FreqFounders=f_founders,FreqEarlyCaptive=f_captive1,
                   FreqLateCaptive=f_captive2,FreqWild=f_wild)
}))

#low but now by category 
listfreq_low_cat <- do.call(rbind,lapply(1:length(unique_variants_low_modern), function(i) {
  variant <- unique_variants_low_modern[i]
  hist <- derived_genotypes_low_modern[derived_genotypes_low_modern$Type=="Historical"&derived_genotypes_low_modern$Variant==variant,]
  f_hist<-(length(hist[hist$Genotype =="0/1",]$Genotype) + length(hist[hist$Genotype =="1/1",]$Genotype)*2)/14
 
   founders <- derived_genotypes_low_modern[which(derived_genotypes_low_modern$SampleID%in%founders_samples&
                                                    derived_genotypes_low_modern$Variant==variant),]
  f_founders<-(length(founders[founders$Genotype =="0/1",]$Genotype) + length(founders[founders$Genotype =="1/1",]$Genotype)*2)/12
  
  captive1 <- derived_genotypes_low_modern[which(derived_genotypes_low_modern$SampleID%in%captive_samples&
                                                   derived_genotypes_low_modern$Variant==variant),]
  f_captive1<-(length(captive1[captive1$Genotype =="0/1",]$Genotype) + length(captive1[captive1$Genotype =="1/1",]$Genotype)*2)/12
  
  captive2 <- derived_genotypes_low_modern[which(derived_genotypes_low_modern$SampleID%in%captive_late_samples&
                                                   derived_genotypes_low_modern$Variant==variant),]
  f_captive2<-(length(captive2[captive2$Genotype =="0/1",]$Genotype) + length(captive2[captive2$Genotype =="1/1",]$Genotype)*2)/12
  
  wild <- derived_genotypes_low_modern[which(derived_genotypes_low_modern$SampleID%in%wild_samples&
                                               derived_genotypes_low_modern$Variant==variant),]
  f_wild<-(length(wild[wild$Genotype =="0/1",]$Genotype) + length(wild[wild$Genotype =="1/1",]$Genotype)*2)/38
  
  cbind.data.frame(Variant=variant,FreqHist=f_hist,FreqFounders=f_founders,FreqEarlyCaptive=f_captive1,
                   FreqLateCaptive=f_captive2,FreqWild=f_wild)
}))


# Save frequencies
#write.csv(listfreq_high_cat, "Files/Frequency_derived_high_captivewild.csv", quote=F, row.names = F)
#write.csv(listfreq_moderate_cat, "Files/Frequency_derived_moderate_captivewild.csv", quote=F, row.names = F)
#write.csv(listfreq_low_cat, "Files/Frequency_derived_low_captivewild.csv", quote=F, row.names = F)

write.csv(listfreq_high_cat, "Files/Frequency_derived_high_captivewild_hist.csv", quote=F, row.names = F)
write.csv(listfreq_moderate_cat, "Files/Frequency_derived_moderate_captivewild_hist.csv", quote=F, row.names = F)
write.csv(listfreq_low_cat, "Files/Frequency_derived_low_captivewild_hist.csv", quote=F, row.names = F)



#remove singletons----
listfreq_high2 <- listfreq_high[-which(listfreq_high$FreqHist==1/14 & listfreq_high$FreqModern==0 ),]
listfreq_high2 <- listfreq_high2[-which(listfreq_high2$FreqModern==1/70 & listfreq_high2$FreqHist==0 ),]

listfreq_moderate2 <- listfreq_moderate[-which(listfreq_moderate$FreqHist==1/14 & listfreq_moderate$FreqModern==0 ),]
listfreq_moderate2 <- listfreq_moderate2[-which(listfreq_moderate2$FreqModern==1/70 & listfreq_moderate2$FreqHist==0 ),]

listfreq_low2 <- listfreq_low[-which(listfreq_low$FreqHist==1/14 & listfreq_low$FreqModern==0 ),]
listfreq_low2 <- listfreq_low2[-which(listfreq_low2$FreqModern==1/70 & listfreq_low2$FreqHist==0 ),]


#unite all
listfreq_high2$Type <- "High"
listfreq_moderate2$Type <- "Moderate"
listfreq_low2$Type <- "Low"

listfreq <- rbind.data.frame(listfreq_high2,listfreq_moderate2, listfreq_low2)

listfreq$Difference <- listfreq$FreqModern- listfreq$FreqHist


# Plot distribution

g <- ggplot(listfreq, aes(Difference, color=Type))
g <- g + geom_density() + theme_classic() +
  geom_vline(xintercept = 0)#+ theme(axis.text.x = element_blank())

g
g <- ggplot(listfreq, aes(Type,Difference, color=Type))
g <- g + geom_violin() + theme_classic() 
#geom_vline(xintercept = 0)#+ theme(axis.text.x = element_blank())

g

#remove sites exclusively in modern----
listfreq_high3 <- listfreq_high[-which(listfreq_high$FreqHist==0 & listfreq_high$FreqModern>0 ),]
listfreq_moderate3 <- listfreq_moderate[-which(listfreq_moderate$FreqHist==0 & listfreq_moderate$FreqModern>0 ),]
listfreq_low3 <- listfreq_low[-which(listfreq_low$FreqHist==0 & listfreq_low$FreqModern>0 ),]


#unite all
listfreq_high3$Type <- "High"
listfreq_moderate3$Type <- "Moderate"
listfreq_low3$Type <- "Low"

listfreq3 <- rbind.data.frame(listfreq_high3,listfreq_moderate3, listfreq_low3)
listfreq3$Difference <- listfreq3$FreqModern- listfreq3$FreqHist
listfreq3$Type <- factor(listfreq3$Type, levels=c("Low","Moderate","High"), ordered = T)

write.csv2(listfreq3,"Files/Change_frequency.csv", quote=F, row.names=F)


# Plot distribution
g <- ggplot(listfreq3, aes(Difference, color=Type))
g <- g + geom_density() + theme_classic() + 
  scale_color_manual(values=c("#4f7cac","#7dde92","#696047"))+ xlab("Freq Modern - Freq Historical")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")#+ theme(axis.text.x = element_blank())

g

ggsave("Plots/Density_frequency.pdf", height = 6, width = 12)


g <- ggplot(listfreq3, aes(Difference, fill=Type))
g <- g + geom_histogram() + theme_classic() + facet_wrap(.~Type, scales = "free_y")+
  scale_fill_manual(values=c("#4f7cac","#7dde92","#696047"))+ xlab("Freq Hist - Freq Modern")+
  theme(legend.position = "none")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")#+ theme(axis.text.x = element_blank())

g
ggsave("Plots/Histogram_frequency.pdf", height = 3, width = 10)


my_comparisons <- list( c("High", "Low"), c("High", "Moderate"), c("Low", "Moderate"))

g1 <- ggplot(listfreq3, aes(Type,Difference, color=Type))
g1 <- g1 +geom_boxplot(outlier.shape = NA ) + geom_violin() + theme_classic() +stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  #geom_vline(xintercept = 0)#+ theme(axis.text.x = element_blank())

g1


g1 <- ggplot(listfreq3[listfreq3$Difference>0,], aes(Type,Difference, color=Type))
g1 <- g1 +geom_boxplot(outlier.shape = NA ) + geom_violin() + theme_classic() +stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
#geom_vline(xintercept = 0)#+ theme(axis.text.x = element_blank())

g1

g <- ggplot(listfreq3, aes(Difference, color=Type))
g <- g + geom_density() + theme_classic() + xlim(c(0,1))+
  scale_color_manual(values=c("#4f7cac","#7dde92","#696047"))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")#+ theme(axis.text.x = element_blank())

g


length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="High",]$Variant)/length(listfreq3[listfreq3$Type=="High",]$Variant)
length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="Moderate",]$Variant)/length(listfreq3[listfreq3$Type=="Moderate",]$Variant)
length(listfreq3[listfreq3$Difference>0 & listfreq3$Type=="Low",]$Variant)/length(listfreq3[listfreq3$Type=="Low",]$Variant)





# exclusively shared sites NOT USED
listfreq_high3 <- listfreq_high[which(listfreq_high$FreqHist>0 & listfreq_high$FreqModern>0 ),]
listfreq_moderate3 <- listfreq_moderate[which(listfreq_moderate$FreqHist>0 & listfreq_moderate$FreqModern>0 ),]
listfreq_low3 <- listfreq_low[which(listfreq_low$FreqHist>0 & listfreq_low$FreqModern>0 ),]

write.csv(listfreq_high3, "Files/Frequency_derived_high_nosingletons.csv", quote=F, row.names = F)
write.csv(listfreq_moderate3, "Files/Frequency_derived_moderate_nosingletons.csv", quote=F, row.names = F)
write.csv(listfreq_low3, "Files/Frequency_derived_low_nosingletons.csv", quote=F, row.names = F)


#unite all
listfreq_high3$Type <- "High"
listfreq_moderate3$Type <- "Moderate"
listfreq_low3$Type <- "Low"

listfreq3 <- rbind.data.frame(listfreq_high3,listfreq_moderate3, listfreq_low3)
listfreq3$Difference <- listfreq3$FreqModern- listfreq3$FreqHist

listfreq3$Type <- factor(listfreq3$Type, levels=c("Low","Moderate","High"), ordered = T)

write.csv2(listfreq3,"Files/Change_frequency_nosingletons.csv", quote=F, row.names=F)

# Plot distribution
g <- ggplot(listfreq3, aes(Difference, color=Type))
g <- g + geom_density() + theme_classic() + 
  scale_color_manual(values=c("#4f7cac","#7dde92","#696047"))+ xlab("Freq Modern - Freq Historical")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")#+ theme(axis.text.x = element_blank())

g

ggsave("Plots/Density_frequency_nosingletons.pdf", height = 6, width = 12)





# Fixed in modern compared to historical, remove sites exclusivly in modern ----
length(listfreq_high[listfreq_high$FreqHist>0&listfreq_high$FreqModern>=0.7,]$Variant)/
  length(listfreq_high[-which(listfreq_high$FreqHist==0 & listfreq_high$FreqModern>0),]$Variant)
length(listfreq_moderate[listfreq_moderate$FreqHist>0&listfreq_moderate$FreqModern>=0.9,]$Variant)/
  length(listfreq_moderate[-which(listfreq_moderate$FreqHist==0 & listfreq_moderate$FreqModern>0),]$Variant)

listfreq_low[listfreq_low$FreqModern>=0.9,]


listfreq_low[listfreq_low$FreqHist==1,]
