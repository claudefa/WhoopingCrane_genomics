# Figure 2 Genetic Load Downsampled to 10 x 
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

setwd("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/")
metadata <- read_excel("Metadata_crane.xlsx", sheet = 1)


samples<- read.table("Files/snpEff/Samples_down")

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

# PLOTTING ----------
# using all ancestral groups
lowcounts_df$Class <- "Low"
moderatecounts_df$Class <- "Moderate"
highcounts_df$Class <- "High"


counts_df <- rbind(highcounts_df,moderatecounts_df, lowcounts_df)
counts_df$Total_5 <-counts_df$het_refanc_5 +counts_df$hom_refanc_5

counts_df_low <- counts_df[counts_df$Class=="Low",]

counts_df$HetCountNormLow_5 <- counts_df$het_refanc_5/counts_df_low$Total_5
counts_df$HomCountNormLow_5 <- counts_df$hom_refanc_5/counts_df_low$Total_5
counts_df$TotalNormLow_5 <- counts_df$Total_5/counts_df_low$Total_5

## plot counts ----
head(counts_df)
counts_df2 <- counts_df[,c("SampleID", "het_refanc_5","hom_refanc_5","Total_5","Class")]
counts_df_gather <- gather(counts_df2, State, Count, -SampleID,  -Class)
counts_df_metadata <- merge(counts_df_gather, metadata, by="SampleID") 

counts_df_metadata$Category<- factor(counts_df_metadata$Category, 
                                     levels=c("Historical","Founders_wildborn","Early_captive","Late_captive","Wild"), 
                                     ordered = T)

counts_df_metadata$Class <- factor(counts_df_metadata$Class, levels=c("Low","Moderate","High"), ordered = T)


counts_df_metadata$State <- gsub("het_refanc_5","Heterozygous Load",counts_df_metadata$State)
counts_df_metadata$State <- gsub("hom_refanc_5","Homozygous Load",counts_df_metadata$State)
counts_df_metadata$State <- gsub("Total_5","Total Load",counts_df_metadata$State)

# counts this one
my_comparisons <- list( c("Founders_wildborn", "Early_captive"), 
                        c("Founders_wildborn", "Late_captive"), 
                        c("Founders_wildborn", "Wild"), 
                        c("Early_captive", "Late_captive"),
                        c("Early_captive", "Wild"),
                        c("Late_captive", "Wild") 
)
g_modern <- ggplot(counts_df_metadata, aes(Category, Count, color=Category ))
g_modern <- g_modern + theme_classic() + facet_grid(Class~State, scales="free") + xlab("")+
  geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitterdodge())+ 
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none") + 
  ylab("Derived Allele Count")+
  #ggtitle("Allele Count") + 
  stat_compare_means(label.x=2)+  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  scale_color_manual(values=c("#8C2717","#E67260","#F1B1A7","#e3e359" ))
g_modern
ggsave("Manuscript/Plots_supplementary/FigureS9_GeneticLoadCount_DP5_down10.pdf", height = 10, width = 10)


## plot counts normalized by low mutations total ----
counts_df3 <- counts_df[,c("SampleID", "Class",
                           "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5")]
counts_df_gather <- gather(counts_df3, State, Count, -SampleID,  -Class)
countsnorm_df_metadata <- merge(counts_df_gather, metadata, by="SampleID") 

countsnorm_df_metadata$Category<- factor(countsnorm_df_metadata$Category, 
                                         levels=c("Historical","Founders_wildborn","Early_captive","Late_captive","Wild"), 
                                         ordered = T)
countsnorm_df_metadata$Class <- factor(countsnorm_df_metadata$Class, levels=c("Low","Moderate","High"), ordered = T)

g <- ggplot(countsnorm_df_metadata[countsnorm_df_metadata$Class!="Low",], 
            aes(State, Count, color=Category))
g <- g + geom_boxplot() + theme_classic() + facet_grid(Class~., scales="free_y") + xlab("") + 
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  ggtitle("All Sites Normalized by total low count") +
  stat_compare_means(comparisons = my_comparisons)+
  scale_color_manual(values=c("#8C2717","#E67260","#F1B1A7","#e3e359" ))+
  stat_compare_means()
  
g




counts_down <- counts_df[,c("SampleID", "Class",
                           "het_refanc_5","hom_refanc_5","Total_5",
                           "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5")]
counts_down$Dataset <- "Downsampled"


# 1st comparison with only 2 samples at different coverage----
samples_cov <- c("EB31_S103_15x","EB31_S103_5x", "EB31_S103",
                 "EB33_S105_15x","EB33_S105_5x","EB33_S105")

total_counts_cov <-counts_down[which(counts_down$SampleID%in%samples_cov),]
total_counts_cov$Coverage <- c(15,5,10,15,5,10,15,5,10,15,5,10,15,5,10,15,5,10)
total_counts_cov$Samples_Final <- total_counts_cov$SampleID
total_counts_cov$SampleID <- gsub("_5x","",total_counts_cov$SampleID) 
total_counts_cov$SampleID <- gsub("_15x","",total_counts_cov$SampleID) 

#Clean names
colnames(total_counts_cov)<- c("SampleID","Class","HeterozygousLoad", "HomozygousLoad","TotalLoad","HetCountNormLow_5",
                               "HomCountNormLow_5","TotalNormLow_5","Dataset", "Coverage","Samples_Final"  )


total_counts_df3 <- total_counts_cov[,c("SampleID", "Class",
                           "HeterozygousLoad","HomozygousLoad","TotalLoad","Dataset",
                           "Coverage","Samples_Final")]

counts_df_gather <- gather(total_counts_df3, State, Count, -SampleID,  -Class, -Dataset,-Coverage,-Samples_Final)

countsnorm_df_metadata <- merge(counts_df_gather, metadata[,c("SampleID","Category")], by="SampleID") 

countsnorm_df_metadata$Category<- factor(countsnorm_df_metadata$Category, 
                                         levels=c("Historical","Founders_wildborn","Early_captive","Late_captive","Wild"), 
                                         ordered = T)
countsnorm_df_metadata$Class <- factor(countsnorm_df_metadata$Class, levels=c("Low","Moderate","High"), ordered = T)
countsnorm_df_metadata$Coverage <- as.character(countsnorm_df_metadata$Coverage)


countsnorm_df_metadata$Coverage <- factor(countsnorm_df_metadata$Coverage, levels=c("5","10","15","Total"), ordered = T)

g <- ggplot(countsnorm_df_metadata, aes(Coverage, Count, color=State, shape=SampleID, group=interaction(SampleID,State)))
g <- g + geom_point()+ geom_line() + theme_classic() + facet_wrap(.~Class, scales="free_y") + xlab("") +
 theme(legend.position = "top")
g

### normalization total low
total_counts_df3 <- total_counts_cov[,c("SampleID", "Class",
                                        "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5","Dataset",
                                        "Coverage","Samples_Final")]

counts_df_gather <- gather(total_counts_df3, State, Count, -SampleID,  -Class, -Dataset,-Coverage,-Samples_Final)

countsnorm_df_metadata <- merge(counts_df_gather, metadata[,c("SampleID","Category")], by="SampleID") 

countsnorm_df_metadata$Category<- factor(countsnorm_df_metadata$Category, 
                                         levels=c("Historical","Founders_wildborn","Early_captive","Late_captive","Wild"), 
                                         ordered = T)
countsnorm_df_metadata$Class <- factor(countsnorm_df_metadata$Class, levels=c("Low","Moderate","High"), ordered = T)
countsnorm_df_metadata$Coverage <- as.character(countsnorm_df_metadata$Coverage)



countsnorm_df_metadata$Coverage <- factor(countsnorm_df_metadata$Coverage, levels=c("5","10","15","Total"), ordered = T)

g1 <- ggplot(countsnorm_df_metadata, aes(Coverage, Count, color=State, shape=SampleID, group=interaction(SampleID,State)))
g1 <- g1 + geom_point()+ geom_line() + theme_classic() + facet_wrap(.~Class, scales="free_y") + xlab("Coverage") +
  ylab("Count/Count_Low")+ theme(legend.position = "none")
g1

ggarrange(g,g1, nrow=2,align="v", heights = c(1.2,1), labels=c("A","B"))
ggsave("Manuscript/Plots_supplementary/FigureS4_Downsampling_load_5_10_15_normalization.pdf", height = 6, width = 10)

#With historical and modern------
all_data <- read.csv2("Files/snpEff/Counts_totalCoverage.csv") # Add info on historical samples
all_data <- all_data[,c("SampleID", "Class",
                        "het_refanc_5","hom_refanc_5","Total_5",
                        "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5", "Dataset")]

historical_samples <- c("Grus27863","Grus71441","MCZ-220028", "MCZ-301263","MCZ-35730",
                        "MCZ-42511","MCZ-46873" )

all_data <- all_data[which(all_data$SampleID%in%historical_samples),] # keep only the historical samples in the whole dataset
total_counts <- rbind.data.frame(all_data,counts_down)

total_counts_down <-total_counts[-which(total_counts$SampleID%in%samples_cov[-c(3,6)]),]


total_counts_df3 <- total_counts_down[,c("SampleID", "Class",
                                        "HetCountNormLow_5","HomCountNormLow_5","TotalNormLow_5","Dataset")]


counts_df_gather <- gather(total_counts_df3, State, Count, -SampleID,  -Class, -Dataset,)

countsnorm_df_metadata <- merge(counts_df_gather, metadata[,c("SampleID","Category")], by="SampleID") 
countsnorm_df_metadata$Class <- factor(countsnorm_df_metadata$Class, levels=c("Low","Moderate","High"), ordered = T)

countsnorm_df_metadata$Type <- countsnorm_df_metadata$Category
countsnorm_df_metadata$Type <- gsub("Founders_wildborn","Modern",countsnorm_df_metadata$Type)
countsnorm_df_metadata$Type <- gsub("Late_captive","Modern",countsnorm_df_metadata$Type)
countsnorm_df_metadata$Type <- gsub("Early_captive","Modern",countsnorm_df_metadata$Type)
countsnorm_df_metadata$Type <- gsub("Wild","Modern",countsnorm_df_metadata$Type)

countsnorm_df_metadata$State <- gsub("HetCountNormLow_5","Heterozygous Load",countsnorm_df_metadata$State)
countsnorm_df_metadata$State <- gsub("HomCountNormLow_5","Homozygous Load",countsnorm_df_metadata$State)
countsnorm_df_metadata$State <- gsub("TotalNormLow_5","Total Load",countsnorm_df_metadata$State)


g_hist <- ggplot(countsnorm_df_metadata[countsnorm_df_metadata$Class!="Low",], 
            aes(Type, Count, color=Type))
g_hist <- g_hist + geom_boxplot(outlier.shape = NA) + geom_jitter()+theme_classic() + facet_grid(Class~State, scales="free_y",space="free_x") + xlab("") + 
   theme(legend.position = "none") + 
  stat_compare_means() + ylab("Normalized Allele Counts")+
  scale_color_manual(values= c("#6d466bff","#EB9486" ))
g_hist

#Change in Frequency of deleterious variation through time. 

listfreq3 <- read.csv2("Files/Change_frequency.csv")

listfreq3$Type <- factor(listfreq3$Type, levels=c("Low","Moderate","High"), ordered = T)

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


# Final Figure ----
ggarrange(g_hist,freq_dist, nrow=1, labels=c("A","B"),widths = c(2,1))

ggsave("Manuscript/MainFigures/Figure2.pdf", height = 5, width = 10)

