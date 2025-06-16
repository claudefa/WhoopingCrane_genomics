# New estimates of Load after revision ---
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(readxl)
library(bootstrap)

setwd("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/")

metadata <- read_excel("Manuscript/Review/Fontsere_etal_whooping_crane_Supp_tables_R1.xlsx", sheet = 2)
###

# Polarization success-----
# get how many sites we lose because of that for each category
#high
outgroup<- c("Gnigricollis","birdAnc314","birdAnc315","birdAnc316")
samples<- read.table("Files/snpEff/Samples_down")

high <- read.table("Files/snpEff/high_down10x.txt")
moderate <- read.table("Files/snpEff/moderate_down10x.txt")
low <- read.table("Files/snpEff/low_down10x.txt")

all_load <- rbind.data.frame(high, moderate,low)

colnames(all_load) <- c("Chrom","Pos","REF","ALT", samples$V1)
all_load$chrompos<- paste0(all_load$Chrom,"-",all_load$Pos)

all_load_outgroup <- all_load[,c(46:50),]
all_load_outgroup$Gnigricollis <- gsub("0/0:[0-9]+","0",all_load_outgroup$Gnigricollis)
all_load_outgroup$birdAnc314 <- gsub("0/0:[0-9]+","0",all_load_outgroup$birdAnc314)
all_load_outgroup$birdAnc315 <- gsub("0/0:[0-9]+","0",all_load_outgroup$birdAnc315)
all_load_outgroup$birdAnc316 <- gsub("0/0:[0-9]+","0",all_load_outgroup$birdAnc316)

all_load_outgroup$Gnigricollis <- gsub("0/1:[0-9]+","1",all_load_outgroup$Gnigricollis)
all_load_outgroup$birdAnc314 <- gsub("0/1:[0-9]+","1",all_load_outgroup$birdAnc314)
all_load_outgroup$birdAnc315 <- gsub("0/1:[0-9]+","1",all_load_outgroup$birdAnc315)
all_load_outgroup$birdAnc316 <- gsub("0/1:[0-9]+","1",all_load_outgroup$birdAnc316)

all_load_outgroup$Gnigricollis <- gsub("1/1:[0-9]+","2",all_load_outgroup$Gnigricollis)
all_load_outgroup$birdAnc314 <- gsub("1/1:[0-9]+","2",all_load_outgroup$birdAnc314)
all_load_outgroup$birdAnc315 <- gsub("1/1:[0-9]+","2",all_load_outgroup$birdAnc315)
all_load_outgroup$birdAnc316 <- gsub("1/1:[0-9]+","2",all_load_outgroup$birdAnc316)

all_load_outgroup$Gnigricollis <- gsub("2/[0-9]+:[0-9]+","3",all_load_outgroup$Gnigricollis)
all_load_outgroup$birdAnc314 <- gsub("2/[0-9]+:[0-9]+","3",all_load_outgroup$birdAnc314)
all_load_outgroup$birdAnc315 <- gsub("2/[0-9]+:[0-9]+","3",all_load_outgroup$birdAnc315)
all_load_outgroup$birdAnc316 <- gsub("2/[0-9]+:[0-9]+","3",all_load_outgroup$birdAnc316)

all_load_outgroup$Gnigricollis <- gsub("[0-9]+/2:[0-9]+","3",all_load_outgroup$Gnigricollis)
all_load_outgroup$birdAnc314 <- gsub("[0-9]+/2:[0-9]+","3",all_load_outgroup$birdAnc314)
all_load_outgroup$birdAnc315 <- gsub("[0-9]+/2:[0-9]+","3",all_load_outgroup$birdAnc315)
all_load_outgroup$birdAnc316 <- gsub("[0-9]+/2:[0-9]+","3",all_load_outgroup$birdAnc316)

all_load_outgroup$Gnigricollis <- gsub("[0-9]+/[0-9]+:[0-9]+","3",all_load_outgroup$Gnigricollis)
all_load_outgroup$birdAnc314 <- gsub("[0-9]+/[0-9]+:[0-9]+","3",all_load_outgroup$birdAnc314)
all_load_outgroup$birdAnc315 <- gsub("[0-9]+/[0-9]+:[0-9]+","3",all_load_outgroup$birdAnc315)
all_load_outgroup$birdAnc316 <- gsub("[0-9]+/[0-9]+:[0-9]+","3",all_load_outgroup$birdAnc316)

x <- length(all_load_outgroup$chrompos)
shared_all <-length(all_load_outgroup[which(all_load_outgroup$Gnigricollis=="0"&all_load_outgroup$birdAnc314=="0"&
                    all_load_outgroup$birdAnc315=="0"&all_load_outgroup$birdAnc316=="0"),]$chrompos)

shared_gn_45 <-length(all_load_outgroup[which(all_load_outgroup$Gnigricollis=="0"&all_load_outgroup$birdAnc314=="0"&
                                 all_load_outgroup$birdAnc315=="0"),]$chrompos)

shared_gn_4 <- length(all_load_outgroup[which(all_load_outgroup$Gnigricollis=="0"&all_load_outgroup$birdAnc314=="0"),]$chrompos)

shared_4_5 <-length(all_load_outgroup[which(all_load_outgroup$birdAnc314=="0"&all_load_outgroup$birdAnc315=="0"),]$chrompos)
shared_4_5_6 <- length(all_load_outgroup[which(all_load_outgroup$birdAnc314=="0"&
                                 all_load_outgroup$birdAnc315=="0"&all_load_outgroup$birdAnc316=="0"),]$chrompos)


shared_out<- data.frame(Total=x, Shared_all=shared_all, Shared_Gn_45=shared_gn_45,
           Shared_Gn_4=shared_gn_4, Shared_45=shared_4_5, Shared_456=shared_4_5_6)

shared_out_gather <- gather(shared_out, Group, Sites, -Total)
shared_out_gather$Perc <- shared_out_gather$Sites/shared_out_gather$Total



###



# Obtain counts for historical ------
# after reviewers comments limit DP max 2times the highest coverage = 16
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
  
  
  #Filter by DP 5 and max DP 16
  df_5_high <- df_3_high[df_3_high$DP >=5& df_3_high$DP<=16,]
  df_5_high$Class <- "High"
  df_5_moderate <- df_3_mod[df_3_mod$DP >=5 &df_3_mod$DP<=16,]
  df_5_moderate$Class <- "Moderate"
  df_5_low <- df_3_low[df_3_low$DP >=5 &df_3_low$DP<=16,]
  df_5_low$Class <- "Low"
  
  df_5 <- rbind.data.frame(df_5_high,df_5_moderate, df_5_low)
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  cbind.data.frame(df_refanc_5[,c(1:4,9)], SampleID=sample_name)
})

historical_genotypes <- do.call(rbind,historical_list)

# write output historical ----

#write.csv2(historical_genotypes, "Files/historical_load_genotypes_R1.csv", quote = F, row.names = F)

# Obtain counts for modern (downsampled to 10x) ------
# after reviewers comments limit DP max 2times the coverage = 20
#HIGH
samples<- read.table("Files/snpEff/Samples_down")

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
  filter(grepl("0/0:[0-9]+", EB30_S102)) %>% filter(grepl("0/0:[0-9]+", EB31_S103)) %>% 
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

#MODERATE 
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
  filter(grepl("0/0:[0-9]+", EB30_S102)) %>% filter(grepl("0/0:[0-9]+", EB31_S103)) %>% 
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

#LOW
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
  filter(grepl("0/0:[0-9]+", EB30_S102)) %>% filter(grepl("0/0:[0-9]+", EB31_S103)) %>% 
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

# extract info genotype 
modern_list<-list()
modern_list <- lapply(1:41, function(i) {
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
  
  
  df_3_high <- cbind.data.frame(df_high_2,high[,c(1,2,46,47,48,49)])
  df_3_high$Gnigricollis <- gsub(":1","", df_3_high$Gnigricollis)
  df_3_high$birdAnc314 <- gsub(":1","", df_3_high$birdAnc314)
  df_3_high$birdAnc315 <- gsub(":1","", df_3_high$birdAnc315)
  df_3_high$birdAnc316 <- gsub(":1","", df_3_high$birdAnc316)
  
  df_3_mod <- cbind.data.frame(df_moderate_2,moderate[,c(1,2,46,47,48,49)])
  df_3_mod$Gnigricollis <- gsub(":1","", df_3_mod$Gnigricollis)
  df_3_mod$birdAnc314 <- gsub(":1","", df_3_mod$birdAnc314)
  df_3_mod$birdAnc315 <- gsub(":1","", df_3_mod$birdAnc315)
  df_3_mod$birdAnc316 <- gsub(":1","", df_3_mod$birdAnc316)
  
  df_3_low <- cbind.data.frame(df_low_2,low[,c(1,2,46,47,48,49)])
  df_3_low$Gnigricollis <- gsub(":1","", df_3_low$Gnigricollis)
  df_3_low$birdAnc314 <- gsub(":1","", df_3_low$birdAnc314)
  df_3_low$birdAnc315 <- gsub(":1","", df_3_low$birdAnc315)
  df_3_low$birdAnc316 <- gsub(":1","", df_3_low$birdAnc316)
  
  
  #Filter by DP 5 and max DP 16
  df_5_high <- df_3_high[df_3_high$DP >=5& df_3_high$DP<=20,]
  df_5_high$Class <- "High"
  df_5_moderate <- df_3_mod[df_3_mod$DP >=5 &df_3_mod$DP<=20,]
  df_5_moderate$Class <- "Moderate"
  df_5_low <- df_3_low[df_3_low$DP >=5 &df_3_low$DP<=20,]
  df_5_low$Class <- "Low"
  
  df_5 <- rbind.data.frame(df_5_high,df_5_moderate, df_5_low)
  
  #Polarize alleles
  df_refanc_5 <- df_5[which(df_5$Gnigricollis == "0/0"&df_5$birdAnc314=="0/0"&df_5$birdAnc315=="0/0"&df_5$birdAnc316=="0/0"),]
  cbind.data.frame(df_refanc_5[,c(1:4,9)], SampleID=sample_name)
})

modern_genotypes <- do.call(rbind,modern_list)

# write output historical ----
#write.csv2(modern_genotypes, "Files/modern_load_genotypes_R1.csv", quote = F, row.names = F)


# Calculate frequencies -----
historical_genotypes <- read.csv2("Files/historical_load_genotypes_R1.csv")
down_10x_genotypes <- read.csv2("Files/modern_load_genotypes_R1.csv")

historical_genotypes$Type <- "Historical"
down_10x_genotypes$Type <- "Modern"


samples_cov <- c("EB31_S103_15x","EB31_S103_5x", 
                 "EB33_S105_15x","EB33_S105_5x")

#only_wild <- c("EB1_S114","EB11_S84","EB12_S85","EB13_S86","EB14_S87","EB15_S88",
#               "EB16_S89","EB17_S90","EB18_S91","EB19_S92","EB2_S115","EB20_S93",
#               "EB21_S94","EB22_S95","EB3_S116","EB6_S117","EB7_S118","EB8_S119","EB9_S83")

down_10x_genotypes <-down_10x_genotypes[-which(down_10x_genotypes$SampleID%in%samples_cov),]
#down_10x_genotypes <-down_10x_genotypes[which(down_10x_genotypes$SampleID%in%only_wild),]


total_genotypes <- rbind.data.frame(historical_genotypes, down_10x_genotypes)
total_genotypes$Variant <- paste0(total_genotypes$Chrom,"-",total_genotypes$Pos)

# extract derived per class
derived_genotypes_high <- total_genotypes[total_genotypes$Genotype!="0/0"&total_genotypes$Class=="High",]
unique_variants_high <- unique(derived_genotypes_high$Variant)

derived_genotypes_moderate <- total_genotypes[total_genotypes$Genotype!="0/0"&total_genotypes$Class=="Moderate",]
unique_variants_moderate <- unique(derived_genotypes_moderate$Variant)

derived_genotypes_low <- total_genotypes[total_genotypes$Genotype!="0/0"&total_genotypes$Class=="Low",]
unique_variants_low <- unique(derived_genotypes_low$Variant)

#loop over variants to extract frequency of variants in each group. 

#high
listfreq_high <- do.call(rbind,lapply(1:length(unique_variants_high), function(i) {
  variant <- unique_variants_high[i]
  hist <- derived_genotypes_high[derived_genotypes_high$Type=="Historical"&derived_genotypes_high$Variant==variant,]
  f_hist<-(length(hist[hist$Genotype =="0/1",]$Genotype) + length(hist[hist$Genotype =="1/1",]$Genotype)*2)/14
  
  modern <- derived_genotypes_high[derived_genotypes_high$Type=="Modern"&derived_genotypes_high$Variant==variant,]
  f_modern<-(length(modern[modern$Genotype =="0/1",]$Genotype) + length(modern[modern$Genotype =="1/1",]$Genotype)*2)/74
  cbind.data.frame(Variant=variant,FreqHist=f_hist,FreqModern=f_modern)
}))

#moderate
listfreq_moderate <- do.call(rbind,lapply(1:length(unique_variants_moderate), function(i) {
  variant <- unique_variants_moderate[i]
  hist <- derived_genotypes_moderate[derived_genotypes_moderate$Type=="Historical"&derived_genotypes_moderate$Variant==variant,]
  f_hist<-(length(hist[hist$Genotype =="0/1",]$Genotype) + length(hist[hist$Genotype =="1/1",]$Genotype)*2)/14
  
  modern <- derived_genotypes_moderate[derived_genotypes_moderate$Type=="Modern"&derived_genotypes_moderate$Variant==variant,]
  f_modern<-(length(modern[modern$Genotype =="0/1",]$Genotype) + length(modern[modern$Genotype =="1/1",]$Genotype)*2)/74
  cbind.data.frame(Variant=variant,FreqHist=f_hist,FreqModern=f_modern)
}))

#low
listfreq_low <- do.call(rbind,lapply(1:length(unique_variants_low), function(i) {
  variant <- unique_variants_low[i]
  hist <- derived_genotypes_low[derived_genotypes_low$Type=="Historical"&derived_genotypes_low$Variant==variant,]
  f_hist<-(length(hist[hist$Genotype =="0/1",]$Genotype) + length(hist[hist$Genotype =="1/1",]$Genotype)*2)/14
  
  modern <- derived_genotypes_low[derived_genotypes_low$Type=="Modern"&derived_genotypes_low$Variant==variant,]
  f_modern<-(length(modern[modern$Genotype =="0/1",]$Genotype) + length(modern[modern$Genotype =="1/1",]$Genotype)*2)/74
  cbind.data.frame(Variant=variant,FreqHist=f_hist,FreqModern=f_modern)
}))

# Save frequencies
listfreq_high$Type <- "High"
listfreq_moderate$Type <- "Moderate"
listfreq_low$Type <- "Low"

write.csv(listfreq_high, "Files/Frequency_derived_high_R1.csv", quote=F, row.names = F)
write.csv(listfreq_moderate, "Files/Frequency_derived_moderate_R1.csv", quote=F, row.names = F)
write.csv(listfreq_low, "Files/Frequency_derived_low_R1.csv", quote=F, row.names = F)

#write.csv(listfreq_high, "Files/Frequency_derived_high_R1_onlywild.csv", quote=F, row.names = F)
#write.csv(listfreq_moderate, "Files/Frequency_derived_moderate_R1_onlywild.csv", quote=F, row.names = F)
#write.csv(listfreq_low, "Files/Frequency_derived_low_R1_onlywild.csv", quote=F, row.names = F)


frequency_all <- rbind.data.frame(listfreq_high,listfreq_moderate,listfreq_low)
frequency_all$Difference <- frequency_all$FreqModern - frequency_all$FreqHist
frequency_all$Type <- factor(frequency_all$Type, levels=c("Low","Moderate","High"), ordered = T)

# No filter ----
freq_dist_all <- ggplot(frequency_all, aes(Difference, color=Type))
freq_dist_all <- freq_dist_all + geom_density() + theme_classic() + 
  scale_color_manual(values=c("#4f7cac","#7dde92","#696047"))+ xlab("Delta Frequency (Modern - Historical)")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")+ theme(legend.position = "top", legend.title = element_blank())

freq_dist_all

x<- length(frequency_all[frequency_all$Difference>0 & frequency_all$Type=="High",]$Variant)
xxx<- length(frequency_all[frequency_all$Difference<0 & frequency_all$Type=="High",]$Variant)
xx <-length(frequency_all[frequency_all$Type=="High",]$Variant)

y<-length(frequency_all[frequency_all$Difference>0 & frequency_all$Type=="Moderate",]$Variant)
yyy<-length(frequency_all[frequency_all$Difference<0 & frequency_all$Type=="Moderate",]$Variant)
yy <-length(frequency_all[frequency_all$Type=="Moderate",]$Variant)

z<-length(frequency_all[frequency_all$Difference>0 & frequency_all$Type=="Low",]$Variant)
zzz<-length(frequency_all[frequency_all$Difference<0 & frequency_all$Type=="Low",]$Variant)
zz<-length(frequency_all[frequency_all$Type=="Low",]$Variant)


df <- rbind.data.frame(c("High",x,xxx,xx),c("Moderate",y,yyy,yy),c("Low",z,zzz,zz))
colnames(df) <- c("Type","VariantsUP","VariantsDown","TotalVariants")


mean(frequency_all[frequency_all$Type=="High",]$Difference)
mean(frequency_all[frequency_all$Type=="Moderate",]$Difference)
mean(frequency_all[frequency_all$Type=="Low",]$Difference)

# Filter : remove variation exclusively present in modern --That is the one we used previously ----
listfreq_3 <- frequency_all[-which(frequency_all$FreqHist==0 & frequency_all$FreqModern>0 ),]
#write.csv2(listfreq_3,"Files/Change_frequency_filterExclMod.csv", quote=F, row.names=F)

freq_dist_3 <- ggplot(listfreq_3, aes(Difference, color=Type))
freq_dist_3 <- freq_dist_3 + geom_density() + theme_classic() + 
  scale_color_manual(values=c("#4f7cac","#7dde92","#696047"))+ xlab("Delta Frequency (Modern - Historical)")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")+ theme(legend.position = "top", legend.title = element_blank())

freq_dist_3

x<- length(listfreq_3[listfreq_3$Difference>0 & listfreq_3$Type=="High",]$Variant)
xxx<- length(listfreq_3[listfreq_3$Difference<0 & listfreq_3$Type=="High",]$Variant)
xx <-length(listfreq_3[listfreq_3$Type=="High",]$Variant)

y<-length(listfreq_3[listfreq_3$Difference>0 & listfreq_3$Type=="Moderate",]$Variant)
yyy<-length(listfreq_3[listfreq_3$Difference<0 & listfreq_3$Type=="Moderate",]$Variant)
yy <-length(listfreq_3[listfreq_3$Type=="Moderate",]$Variant)

z<-length(listfreq_3[listfreq_3$Difference>0 & listfreq_3$Type=="Low",]$Variant)
zzz<-length(listfreq_3[listfreq_3$Difference<0 & listfreq_3$Type=="Low",]$Variant)
zz<-length(listfreq_3[listfreq_3$Type=="Low",]$Variant)


df_3 <- rbind.data.frame(c("High",x,xxx,xx),c("Moderate",y,yyy,yy),c("Low",z,zzz,zz))
colnames(df_3) <- c("Type","VariantsUP","VariantsDown","TotalVariants")


# Filter: remove variation exclusively present in historical -----
listfreq_4 <- frequency_all[-which(frequency_all$FreqHist>0 & frequency_all$FreqModern==0 ),]
#write.csv2(listfreq_4,"Files/Change_frequency_filterExclHist.csv", quote=F, row.names=F)

freq_dist_4 <- ggplot(listfreq_4, aes(Difference, color=Type))
freq_dist_4 <- freq_dist_4 + geom_density() + theme_classic() + 
  scale_color_manual(values=c("#4f7cac","#7dde92","#696047"))+ xlab("Delta Frequency (Modern - Historical)")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")+ theme(legend.position = "top", legend.title = element_blank())

freq_dist_4

x<- length(listfreq_4[listfreq_4$Difference>0 & listfreq_4$Type=="High",]$Variant)
xxx<- length(listfreq_4[listfreq_4$Difference<0 & listfreq_4$Type=="High",]$Variant)
xx <-length(listfreq_4[listfreq_4$Type=="High",]$Variant)

y<-length(listfreq_4[listfreq_4$Difference>0 & listfreq_4$Type=="Moderate",]$Variant)
yyy<-length(listfreq_4[listfreq_4$Difference<0 & listfreq_4$Type=="Moderate",]$Variant)
yy <-length(listfreq_4[listfreq_4$Type=="Moderate",]$Variant)

z<-length(listfreq_4[listfreq_4$Difference>0 & listfreq_4$Type=="Low",]$Variant)
zzz<-length(listfreq_4[listfreq_4$Difference<0 & listfreq_4$Type=="Low",]$Variant)
zz<-length(listfreq_4[listfreq_4$Type=="Low",]$Variant)


df_4 <- rbind.data.frame(c("High",x,xxx,xx),c("Moderate",y,yyy,yy),c("Low",z,zzz,zz))
colnames(df_4) <- c("Type","VariantsUP","VariantsDown","TotalVariants")

# Filter: only shared variants -----
listfreq_5<- frequency_all[-which(frequency_all$FreqHist==0 | frequency_all$FreqModern==0 ),]
#write.csv2(listfreq_5,"Files/Change_frequency_filterOnlyShared.csv", quote=F, row.names=F)

freq_dist_5<- ggplot(listfreq_5, aes(Difference, color=Type))
freq_dist_5 <- freq_dist_5 + geom_density() + theme_classic() + 
  scale_color_manual(values=c("#4f7cac","#7dde92","#696047"))+ xlab("Delta Frequency (Modern - Historical)")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")+ theme(legend.position = "top", legend.title = element_blank())

freq_dist_5

x<- length(listfreq_5[listfreq_5$Difference>0 & listfreq_5$Type=="High",]$Variant)
xxx<- length(listfreq_5[listfreq_5$Difference<0 & listfreq_5$Type=="High",]$Variant)
xx <-length(listfreq_5[listfreq_5$Type=="High",]$Variant)

y<-length(listfreq_5[listfreq_5$Difference>0 & listfreq_5$Type=="Moderate",]$Variant)
yyy<-length(listfreq_5[listfreq_5$Difference<0 & listfreq_5$Type=="Moderate",]$Variant)
yy <-length(listfreq_5[listfreq_5$Type=="Moderate",]$Variant)

z<-length(listfreq_5[listfreq_5$Difference>0 & listfreq_5$Type=="Low",]$Variant)
zzz<-length(listfreq_5[listfreq_5$Difference<0 & listfreq_5$Type=="Low",]$Variant)
zz<-length(listfreq_5[listfreq_5$Type=="Low",]$Variant)


df_5 <- rbind.data.frame(c("High",x,xxx,xx),c("Moderate",y,yyy,yy),c("Low",z,zzz,zz))
colnames(df_5) <- c("Type","VariantsUP","VariantsDown","TotalVariants")

# RXY -----
types <- c("High","Moderate", "Low")

# All sites ------
Rxy_frame<-data.frame(matrix(ncol = 4, nrow = 0))
Rxy_jn<-data.frame(matrix(ncol = 3, nrow = 0))

y<-frequency_all$FreqHist*(1-frequency_all$FreqModern)
x<-frequency_all$FreqModern*(1-frequency_all$FreqHist)
xdata_total<-cbind.data.frame(x,y, Type=frequency_all$Type)

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

rxy_plot_all <- ggplot(Rxy_frame_normalized, aes(Type, Rxy_norm_2, fill=Type))
rxy_plot_all <- rxy_plot_all + geom_col(position = "dodge") + theme_classic()+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) +scale_fill_manual(values=c("#696047","#7dde92")) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + ylab("Rxy (Modern/Historical)") + theme(legend.position = "none")
rxy_plot_all

# Remove Exclusively Modern sites ------
Rxy_frame<-data.frame(matrix(ncol = 4, nrow = 0))
Rxy_jn<-data.frame(matrix(ncol = 3, nrow = 0))

y<-listfreq_3$FreqHist*(1-listfreq_3$FreqModern)
x<-listfreq_3$FreqModern*(1-listfreq_3$FreqHist)
xdata_total<-cbind.data.frame(x,y, Type=listfreq_3$Type)

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

rxy_plot_RmExclMod <- ggplot(Rxy_frame_normalized, aes(Type, Rxy_norm_2, fill=Type))
rxy_plot_RmExclMod <- rxy_plot_RmExclMod + geom_col(position = "dodge") + theme_classic()+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) +scale_fill_manual(values=c("#696047","#7dde92")) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + ylab("Rxy (Modern/Historical)") + theme(legend.position = "none")
rxy_plot_RmExclMod

# Remove Exclusively Historical sites ------
Rxy_frame<-data.frame(matrix(ncol = 4, nrow = 0))
Rxy_jn<-data.frame(matrix(ncol = 3, nrow = 0))

y<-listfreq_4$FreqHist*(1-listfreq_4$FreqModern)
x<-listfreq_4$FreqModern*(1-listfreq_4$FreqHist)
xdata_total<-cbind.data.frame(x,y, Type=listfreq_4$Type)

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

rxy_plot_RmExclHist <- ggplot(Rxy_frame_normalized, aes(Type, Rxy_norm_2, fill=Type))
rxy_plot_RmExclHist <- rxy_plot_RmExclHist + geom_col(position = "dodge") + theme_classic()+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) +scale_fill_manual(values=c("#696047","#7dde92")) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + ylab("Rxy (Modern/Historical)") + theme(legend.position = "none")
rxy_plot_RmExclHist

# Keep only shared sites ------
Rxy_frame<-data.frame(matrix(ncol = 4, nrow = 0))
Rxy_jn<-data.frame(matrix(ncol = 3, nrow = 0))

y<-listfreq_5$FreqHist*(1-listfreq_5$FreqModern)
x<-listfreq_5$FreqModern*(1-listfreq_5$FreqHist)
xdata_total<-cbind.data.frame(x,y, Type=listfreq_5$Type)

for (j in types){
  print(j)
  xdata <- xdata_total[xdata_total$Type==j,]
  
  size <- floor(length(xdata$x)*0.02) # take 1% of the sites for the jacknife
  size <- ceiling(length(xdata$x)*0.02) # take 1% of the sites for the jacknife ONLY FOR THE SHARED FILES
  rxy<-function(a,xdata,b){
    sum(xdata[seq(a,a+(size-1),1),b])
  } 
  rx_ry<-data.frame()
  for(i in c(1:100)) {
    rx_ry[i,1]<-rxy(size*(i-1)+1,xdata,1)
    rx_ry[i,2]<-rxy(size*(i-1)+1,xdata,2)
  } 
  
  rat<-function(c,rx_ry){sum(rx_ry[c,1])/sum(rx_ry[c,2])}
  jn_h<-jackknife(1:73,rat,rx_ry) #do jackknife for dataframe rx_ry with the function of rat for 100 times
  Rxy_jn<-rbind(Rxy_jn, cbind(jn_h$jack.values, rep(j, 73)))
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

rxy_plot_Shared <- ggplot(Rxy_frame_normalized, aes(Type, Rxy_norm_2, fill=Type))
rxy_plot_Shared <- rxy_plot_Shared + geom_col(position = "dodge") + theme_classic()+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) +scale_fill_manual(values=c("#696047","#7dde92")) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + ylab("Rxy (Modern/Historical)") + theme(legend.position = "none")
rxy_plot_Shared


# COUNTS ----
samples_total <- read.table("Files/snpEff/Samples")
samples_total <- samples_total$V1[1:44]
#samples_total<- c(samples_total[1:7],only_wild)

#extract subset of sites 
sites_rmExclMod <- listfreq_3$Variant
sites_rmExclHist <- listfreq_4$Variant
sites_Shared <- listfreq_5$Variant

#extract counts
df_counts <- data.frame()
for (i in samples_total) {
for (j in types){
  type_df <- total_genotypes[total_genotypes$Class==j&total_genotypes$SampleID==i,]
  type_df_3 <- total_genotypes[which(total_genotypes$Class==j&total_genotypes$SampleID==i & total_genotypes$Variant%in%sites_rmExclMod),]
  type_df_4 <- total_genotypes[which(total_genotypes$Class==j&total_genotypes$SampleID==i & total_genotypes$Variant%in%sites_rmExclHist),]
  type_df_5 <- total_genotypes[which(total_genotypes$Class==j&total_genotypes$SampleID==i & total_genotypes$Variant%in%sites_Shared),]
  
  df_counts <- rbind(df_counts,data.frame(SampleID=i,
                              Class=j,
                              het_refanc_All=length(type_df[type_df$Genotype=="0/1",]$Genotype),
                              hom_refanc_All=length(type_df[type_df$Genotype=="1/1",]$Genotype)*2,
                              het_refanc_3=length(type_df_3[type_df_3$Genotype=="0/1",]$Genotype),
                              hom_refanc_3=length(type_df_3[type_df_3$Genotype=="1/1",]$Genotype)*2,
                              het_refanc_4=length(type_df_4[type_df_4$Genotype=="0/1",]$Genotype),
                              hom_refanc_4=length(type_df_4[type_df_4$Genotype=="1/1",]$Genotype)*2,
                              het_refanc_5=length(type_df_5[type_df_5$Genotype=="0/1",]$Genotype),
                              hom_refanc_5=length(type_df_5[type_df_5$Genotype=="1/1",]$Genotype)*2))
}}

# Gather dataset to perform normalization
df_counts$Total_All <-df_counts$het_refanc_All + df_counts$hom_refanc_All
df_counts$Total_3 <-df_counts$het_refanc_3 + df_counts$hom_refanc_3
df_counts$Total_4 <-df_counts$het_refanc_4 + df_counts$hom_refanc_4
df_counts$Total_5 <-df_counts$het_refanc_5 + df_counts$hom_refanc_5

#write.csv2(df_counts,"Files/Counts_allfilter.csv", quote=F, row.names=F)

df_counts_low <- df_counts[df_counts$Class=="Low",]

df_counts_gather <- gather(df_counts, Load, Count, -SampleID, -Class)
df_counts_gather$Sites <- df_counts_gather$Load
df_counts_gather$Sites <-gsub("het_refanc_","", df_counts_gather$Sites)
df_counts_gather$Sites <-gsub("hom_refanc_","", df_counts_gather$Sites)
df_counts_gather$Sites <-gsub("Total_","", df_counts_gather$Sites)

df_counts_gather$Load <-gsub("_refanc_All","", df_counts_gather$Load)
df_counts_gather$Load <-gsub("_refanc_3","", df_counts_gather$Load)
df_counts_gather$Load <-gsub("_refanc_4","", df_counts_gather$Load)
df_counts_gather$Load <-gsub("_refanc_5","", df_counts_gather$Load)
df_counts_gather$Load <-gsub("_3","", df_counts_gather$Load)
df_counts_gather$Load <-gsub("_4","", df_counts_gather$Load)
df_counts_gather$Load <-gsub("_5","", df_counts_gather$Load)

df_counts_gather$LoadNorm <- "NA"

for (i in samples_total) {
  df_counts_gather[df_counts_gather$Sites=="All"&df_counts_gather$SampleID==i,]$LoadNorm <- 
    df_counts_gather[df_counts_gather$Sites=="All"&df_counts_gather$SampleID==i,]$Count/df_counts_low[df_counts_low$SampleID==i,]$Total_All
  df_counts_gather[df_counts_gather$Sites=="3"&df_counts_gather$SampleID==i,]$LoadNorm <- 
    df_counts_gather[df_counts_gather$Sites=="3"&df_counts_gather$SampleID==i,]$Count/df_counts_low[df_counts_low$SampleID==i,]$Total_3
  df_counts_gather[df_counts_gather$Sites=="4"&df_counts_gather$SampleID==i,]$LoadNorm <- 
    df_counts_gather[df_counts_gather$Sites=="4"&df_counts_gather$SampleID==i,]$Count/df_counts_low[df_counts_low$SampleID==i,]$Total_4
  df_counts_gather[df_counts_gather$Sites=="5"&df_counts_gather$SampleID==i,]$LoadNorm <- 
    df_counts_gather[df_counts_gather$Sites=="5"&df_counts_gather$SampleID==i,]$Count/df_counts_low[df_counts_low$SampleID==i,]$Total_5
  
}


# add metadata
df_counts_gather_metadata <- merge(df_counts_gather, metadata[,c("SampleID","Category")], by="SampleID") 
df_counts_gather_metadata$Class <- factor(df_counts_gather_metadata$Class, levels=c("Low","Moderate","High"), ordered = T)

df_counts_gather_metadata$Type <- df_counts_gather_metadata$Category
df_counts_gather_metadata$Type <- gsub("Founders_wildborn","Modern",df_counts_gather_metadata$Type)
df_counts_gather_metadata$Type <- gsub("Late_captive","Modern",df_counts_gather_metadata$Type)
df_counts_gather_metadata$Type <- gsub("Early_captive","Modern",df_counts_gather_metadata$Type)
df_counts_gather_metadata$Type <- gsub("Wild","Modern",df_counts_gather_metadata$Type)

df_counts_gather_metadata$LoadNorm <- as.numeric(df_counts_gather_metadata$LoadNorm)


df_counts_gather_metadata$Load <- gsub("het","Heterozygous Count",df_counts_gather_metadata$Load)
df_counts_gather_metadata$Load <- gsub("hom","Homozygous Count",df_counts_gather_metadata$Load)
df_counts_gather_metadata$Load <- gsub("Total_All","Total Count",df_counts_gather_metadata$Load)

df_counts_gather_metadata_all <- df_counts_gather_metadata[df_counts_gather_metadata$Sites=="All",]

g_hist <- ggplot(df_counts_gather_metadata_all[df_counts_gather_metadata_all$Class!="Low",], 
                 aes(Type, LoadNorm, color=Type))
g_hist <- g_hist + geom_boxplot(outlier.shape = NA) + geom_jitter()+theme_classic() + facet_grid(Class~Load, scales="free_y",space="free_x") + xlab("") + 
  theme(legend.position = "none") + 
  stat_compare_means() + ylab("Normalized Allele Counts")+
  scale_color_manual(values= c("#6d466bff","#c99e00" ))
g_hist


df_counts_gather_metadata_3 <- df_counts_gather_metadata[df_counts_gather_metadata$Sites=="3",]

g_hist_mod <- ggplot(df_counts_gather_metadata_3[df_counts_gather_metadata_3$Class!="Low",], 
                 aes(Type, LoadNorm, color=Type))
g_hist_mod <- g_hist_mod + geom_boxplot(outlier.shape = NA) + geom_jitter()+theme_classic() + facet_grid(Class~Load, scales="free_y",space="free_x") + xlab("") + 
  theme(legend.position = "none") + 
  stat_compare_means() + ylab("Normalized Allele Counts")+
  scale_color_manual(values= c("#6d466bff","#c99e00" ))
g_hist_mod



df_counts_gather_metadata_4 <- df_counts_gather_metadata[df_counts_gather_metadata$Sites=="4",]

g_hist_hist <- ggplot(df_counts_gather_metadata_4[df_counts_gather_metadata_4$Class!="Low",], 
                     aes(Type, LoadNorm, color=Type))
g_hist_hist <- g_hist_hist + geom_boxplot(outlier.shape = NA) + geom_jitter()+theme_classic() + facet_grid(Class~Load, scales="free_y",space="free_x") + xlab("") + 
  theme(legend.position = "none") + 
  stat_compare_means() + ylab("Normalized Allele Counts")+
  scale_color_manual(values= c("#6d466bff","#c99e00" ))
g_hist_hist


df_counts_gather_metadata_5 <- df_counts_gather_metadata[df_counts_gather_metadata$Sites=="5",]

g_hist_shared <- ggplot(df_counts_gather_metadata_5[df_counts_gather_metadata_5$Class!="Low",], 
                      aes(Type, LoadNorm, color=Type))
g_hist_shared <- g_hist_shared + geom_boxplot(outlier.shape = NA) + geom_jitter()+theme_classic() + facet_grid(Class~Load, scales="free_y",space="free_x") + xlab("") + 
  theme(legend.position = "none") + 
  stat_compare_means() + ylab("Normalized Allele Counts")+
  scale_color_manual(values= c("#6d466bff","#c99e00" ))
g_hist_shared



# FINAL FIGURE ----

ggarrange(g_hist, ggarrange(freq_dist_all, rxy_plot_all, nrow=2, labels=c("","C")) ,nrow=1, labels=c("A","B"),widths = c(2,1))
#ggsave(paste0("Plots/Genetic_load_sites_filter/GenLoad_allSites_onlywild.pdf"), height = 6, width = 10)
ggsave(paste0("Plots/Genetic_load_sites_filter/GenLoad_allSites.pdf"), height = 6, width = 10)

ggarrange(g_hist_mod, ggarrange(freq_dist_3, rxy_plot_RmExclMod, nrow=2, labels=c("","C")) ,nrow=1, labels=c("A","B"),widths = c(2,1))
ggsave(paste0("Plots/Genetic_load_sites_filter/GenLoad_RmExclMod.pdf"), height = 6, width = 10)
#ggsave(paste0("Plots/Genetic_load_sites_filter/GenLoad_RmExclMod_onlywild.pdf"), height = 6, width = 10)

ggarrange(g_hist_hist, ggarrange(freq_dist_4, rxy_plot_RmExclHist, nrow=2, labels=c("","C")) ,nrow=1, labels=c("A","B"),widths = c(2,1))
ggsave(paste0("Plots/Genetic_load_sites_filter/GenLoad_RmExclHist.pdf"), height = 6, width = 10)

ggarrange(g_hist_shared, ggarrange(freq_dist_5, rxy_plot_Shared, nrow=2, labels=c("","C")) ,nrow=1, labels=c("A","B"),widths = c(2,1))
#ggsave(paste0("Plots/Genetic_load_sites_filter/GenLoad_Shared_onlywild.pdf"), height = 6, width = 10)
ggsave(paste0("Plots/Genetic_load_sites_filter/GenLoad_Shared.pdf"), height = 6, width = 10)

