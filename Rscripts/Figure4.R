# Figure 4 Effect of captive breeding ---
library("ggplot2")
library("tidyr")
library("ggrepel")
library("readxl")
library("ggpubr")

setwd("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/")
metadata <- read_excel("Metadata_crane.xlsx", sheet = 1)
colors_category_2 <- c("#8C2717","#E67260","#F1B1A7","#e3e359" )
colors_category <- c("#EB9486","#e3e359")

metadata_nonhist <- metadata[metadata$Type!="Historical",]


#Figure 4A Heterozygosity -----
metadata_nonhist$Category <- factor(metadata_nonhist$Category, 
                                       levels=c("Founders_wildborn","Early_captive" ,"Late_captive","Wild"),
                                       ordered = T)


p <- ggplot(metadata_nonhist, aes(Category, `Heterozygosity (all Pos only Modern)`,color=Category))
p <- p+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ 
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("")+
  theme(legend.position = "none") + scale_color_manual(values=colors_category_2) 
p

my_comparisons <- list( c("Founders_wildborn", "Early_captive"), 
                        c("Founders_wildborn", "Late_captive"), 
                        c("Founders_wildborn", "Wild"), 
                        c("Early_captive", "Late_captive"),
                        c("Early_captive", "Wild"),
                        c("Late_captive", "Wild") 
                        )

# downsampled to 10x only min10x

het <- ggplot(metadata_nonhist[metadata_nonhist$Coverage>=10,], aes(Category, `Heterozygosity (all Pos only Modern 10x)`,color=Category))
het <- het+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ 
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("")+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y =0.00105, label.x=2.2)+
  theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1)) +
  scale_color_manual(values=colors_category_2) 
het

mean(metadata_nonhist[metadata_nonhist$Coverage>=10&metadata_nonhist$Category=="Wild",]$`Heterozygosity (all Pos only Modern 10x)`)
mean(metadata_nonhist[metadata_nonhist$Coverage>=10&metadata_nonhist$Category=="Late_captive",]$`Heterozygosity (all Pos only Modern 10x)`)


# supplementary heterozygosity and coverage
summary(lm(metadata[metadata$Category!="Historical",]$`Heterozygosity (all Pos only Modern)`~metadata[metadata$Category!="Historical",]$Coverage))

p1 <- ggplot()
p1 <- p1+ geom_point(data=metadata[metadata$Category!="Historical",],aes(`Heterozygosity (all Pos only Modern)`, Coverage, color=Type)) +
  theme_classic() + xlab("Heterozygosity (bp-1)")+ ylab("Coverage")+
  geom_smooth(data=metadata[metadata$Category!="Historical",],
              aes(`Heterozygosity (all Pos only Modern)`, Coverage), method = "lm")+
  theme(legend.position = c(0.13,0.8), legend.title = element_blank()) + scale_color_manual(values=colors_category) 
p1

summary(lm(metadata[metadata$Category!="Historical",]$`Heterozygosity (all Pos only Modern 10x)`~metadata[metadata$Category!="Historical",]$Coverage))

p2 <- ggplot()
p2 <- p2+ geom_point(data=metadata[metadata$Category!="Historical",],aes(`Heterozygosity (all Pos only Modern 10x)`, Coverage, color=Type)) +
  theme_classic() + xlab("Heterozygosity (bp-1)")+ ylab("Coverage")+
  geom_smooth(data=metadata[metadata$Category!="Historical",],
              aes(`Heterozygosity (all Pos only Modern 10x)`, Coverage), method = "lm")+
  theme(legend.position = c(0.13,0.8), legend.title = element_blank()) + scale_color_manual(values=colors_category) 
p2

ggarrange(p1,p2, nrow=1)
ggsave("Manuscript/Plots_supplementary/FigureS5_covhet_10.pdf", height = 4,width = 11)

# Figure 4B Inbreeding ----
roh_perc_df <- read.csv2("Files/Froh_total.csv")
roh_perc_df_gather <- gather(roh_perc_df, key, value, -sample, -Median,-Average,Average2, -count,-sumMb,-maxMb)
roh_perc_df_gather$key <- gsub("Intermediate2","2.5-5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate","1-2.5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("ExtraLong",">10Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long","5-10Mb", roh_perc_df_gather$key )
order <- c("Total","Average2","1-2.5Mb","2.5-5Mb", "5-10Mb",">10Mb")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$SampleID <- roh_perc_df_gather$sample

roh_perc_df_gather2 <- merge(roh_perc_df_gather, metadata, by="SampleID")

roh_perc_df_gather4 <- roh_perc_df_gather2[which(roh_perc_df_gather2$key=="Total"),]
order_category <- c("Historical","Founders_wildborn","Early_captive",
                    "Late_captive", "Wild")
roh_perc_df_gather4$Category <- factor(roh_perc_df_gather4$Category,
                                       levels=order_category, ordered=T)

froh <- ggplot(roh_perc_df_gather4[roh_perc_df_gather4$Category!="Historical"
                                  &roh_perc_df_gather4$Coverage>10,], aes(Category, value*100, color=Category))
froh <- froh+ geom_boxplot(outlier.shape = NA) +  geom_jitter()+ 
  theme_classic() + ylab("FROH (%)")+ xlab("")+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y =26, label.x=2.3)+
  theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1)) + scale_color_manual(values=colors_category_2) 
froh


roh_perc_df_gather4$Category2 <- roh_perc_df_gather4$Category

roh_perc_df_gather4$Category2 <- gsub("Early_captive", "Captive_born",roh_perc_df_gather4$Category2)
roh_perc_df_gather4$Category2 <- gsub("Late_captive", "Captive_born",roh_perc_df_gather4$Category2)
roh_perc_df_gather4$Category2 <- factor(roh_perc_df_gather4$Category2,
                               levels=c("Founders_wildborn","Captive_born","Wild"), ordered=T)


my_comparisons2 <- list( c("Founders_wildborn", "Captive_born"), 
                         c("Founders_wildborn", "Wild"), 
                         c("Captive_born", "Wild"))


froh2 <- ggplot(roh_perc_df_gather4[roh_perc_df_gather4$Category!="Historical"
                                   &roh_perc_df_gather4$Coverage>10,], aes(Category2, value*100, color=Category2))
froh2 <- froh2+ geom_boxplot(outlier.shape = NA) +  geom_jitter()+ 
  theme_classic() + ylab("FROH (%)")+ xlab("")+
  stat_compare_means(comparisons = my_comparisons2)+
  stat_compare_means(label.y =22, label.x=2)+
  theme(legend.position = "none") + scale_color_manual(values=colors_category_2[c(1,3,4)]) 
froh2



# Average Froh per category and total 
mean(roh_perc_df_gather4$value)*100
mean(roh_perc_df_gather4[roh_perc_df_gather4$Category=="Wild",]$value)*100
mean(roh_perc_df_gather4[roh_perc_df_gather4$Category=="Founders_wildborn",]$value)*100
mean(roh_perc_df_gather4[roh_perc_df_gather4$Category=="Early_captive",]$value)*100
mean(roh_perc_df_gather4[roh_perc_df_gather4$Category=="Late_captive",]$value)*100



#To supplementary 

# Comparison Froh all dataset vs downsampled 10x
roh_final_total_10x <- roh_perc_df
roh_final_total_10x <- read.csv2("Files/Froh_total_10x.csv")

roh_final_comparison <- merge(roh_final_total_10x, roh_perc_df, by="sample")


summary(lm(roh_final_comparison$Total.x~roh_final_comparison$Total.y))

g <- ggplot(roh_final_comparison, aes(Total.x, Total.y))
g <- g +geom_point() + theme_classic() + xlab("FROH - Downsampled to 10x")+
  ylab("FROH - Total Dataset") +geom_smooth(method = "lm")
g

ggsave("Manuscript/Plots_supplementary/FigureS6_corFroh_10x.pdf", height = 4,width = 6)



roh_perc_df_gather3 <- roh_perc_df_gather2[which(roh_perc_df_gather2$key!="Total"&
                                                   roh_perc_df_gather2$key!="Average2" ),]
order_category <- c("Historical","Founders_wildborn","Early_captive",
                    "Late_captive", "Wild")
roh_perc_df_gather3$Category <- factor(roh_perc_df_gather3$Category,
                                       levels=order_category, ordered=T)

pa <- ggplot(roh_perc_df_gather3[roh_perc_df_gather3$Type!="Historical",], 
             aes(sample,value,fill=key))
pa <- pa + geom_col()+ theme_classic() +  xlab("")+ facet_grid(.~Category,  scales="free")+
  ylab("FROH") + 
  theme(legend.position = "right", axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=rev(c("#957583","#356255","#AAC0AF","#E8EEE9")), name="RoH Size") 
pa

#ggsave("Manuscript/Plots_supplementary/FigureS7_barplot_modern_size.pdf", width = 14, height = 4)


p <- ggplot(roh_perc_df_gather3[roh_perc_df_gather3$Category!="Historical"&roh_perc_df_gather3$Coverage>10,], aes(Category, value*100, color=Category))
p <- p+ geom_boxplot(outlier.shape = NA) +  geom_jitter()+ 
  facet_wrap(key~.,scales="free_y", nrow=1) +
  theme_classic() + ylab("FROH (%)")+ xlab("")+
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1)) + scale_color_manual(values=colors_category_2) + 
  #ggtitle("FROH 1Mb window size (>10x)")+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.x=1.5)

p

ggsave("Manuscript/Plots_supplementary/FigureS8_boxplot_sizes.pdf", width = 14, height = 7)

# Figure 4 C Size 10 ----
roh_final_rohan <- read.csv2("Files/Roh_persample_metadata.csv")
finalV1 <- roh_perc_df$sample
roh_perc_2 <- lapply(1:length(finalV1), function(i) {
  df2 <- roh_final_rohan[roh_final_rohan$Sample == finalV1[i],]
  data.frame(sample=finalV1[i],
             Small=sum(df2[df2$ROH_LENGTH<10000000,]$ROH_LENGTH)/1112400596,  ## do this but for category
             Long=sum(df2[df2$ROH_LENGTH>=10000000,]$ROH_LENGTH)/1112400596) })

roh_perc_df_2 <- do.call(rbind, roh_perc_2)

roh_perc_df_gather <- gather(roh_perc_df_2, key, value, -sample)
roh_perc_df_gather$key <- gsub("Small","<10Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long",">10Mb", roh_perc_df_gather$key )
order <- c("<10Mb",">10Mb")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$SampleID <- roh_perc_df_gather$sample

roh_perc_2 <- merge(roh_perc_df_gather, metadata, by="SampleID")
roh_perc_2 <- roh_perc_2[roh_perc_2$Type!="Historical",]

roh_perc_2$Category <- factor(roh_perc_2$Category,
                              levels=order_category, ordered=T)

p <- ggplot(roh_perc_2[roh_perc_2$Coverage>10,], aes(Category, value*100, color=Category))
p <- p+ geom_boxplot(outlier.shape = NA) +  geom_jitter()+ 
  facet_grid(.~key) +
  theme_classic() + ylab("FROH (%)")+ xlab("")+
  theme(legend.position = "right") + scale_color_manual(values=colors_category_2) + 
  ggtitle("FROH 1Mb window size")+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y =22, label.x=2)
p

roh_perc_2$Category2 <- roh_perc_2$Category

roh_perc_2$Category2 <- gsub("Early_captive", "Captive_born",roh_perc_2$Category2)
roh_perc_2$Category2 <- gsub("Late_captive", "Captive_born",roh_perc_2$Category2)
roh_perc_2$Category2 <- factor(roh_perc_2$Category2,
                              levels=c("Founders_wildborn","Captive_born","Wild"), ordered=T)


my_comparisons2 <- list( c("Founders_wildborn", "Captive_born"), 
                        c("Founders_wildborn", "Wild"), 
                        c("Captive_born", "Wild"))

froh3 <- ggplot(roh_perc_2[roh_perc_2$Coverage>10,], aes(Category2, value*100, color=Category2))
froh3 <- froh3+ geom_boxplot(outlier.shape = NA) +  geom_jitter()+ 
  facet_grid(.~key) +
  theme_classic() + ylab("FROH (%)")+ xlab("")+
  theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1)) + scale_color_manual(values=colors_category_2[c(1,3,4)]) + 
  stat_compare_means(comparisons = my_comparisons2)+
  stat_compare_means(label.y =10, label.x=1)
froh3

# averages
mean(roh_perc_2[roh_perc_2$key==">10Mb",]$value)*100
mean(roh_perc_2[roh_perc_2$Category=="Wild"&roh_perc_2$key==">10Mb",]$value)*100
mean(roh_perc_2[roh_perc_2$Category=="Founders_wildborn"&roh_perc_2$key==">10Mb",]$value)*100
mean(roh_perc_2[roh_perc_2$Category=="Early_captive"&roh_perc_2$key==">10Mb",]$value)*100
mean(roh_perc_2[roh_perc_2$Category=="Late_captive"&roh_perc_2$key==">10Mb",]$value)*100



# Figure 4D generations ago -----

species <-c('Red_neckedCrane','quail','CollFly_female','CollFly_male','zebraFinch','HeHo_male',
            'HeHo_female','HouseSp')
r=c(3.42,1.9,2.28,3.56,3.18,1.86,1.71,4.82)
list_gen <- list()
for (i in 1:length(r)){  
  recomb_spcies <- as.character(species[i])
  recomb_rate <- r[i]
  Tgen <- (100/((roh_final_rohan$ROH_LENGTH/1000000)*recomb_rate))/2
  list_gen[[i]] <- Tgen
}

roh_final_rohan$RedNeckedCrane <- list_gen[[1]]
roh_final_rohan$CollFly_female <- list_gen[[3]]

# summary ROH TGen CollFly_female NOT DONE ------
roh_perc <- list()

roh_perc <- lapply(1:length(finalV1), function(i) {
  df2 <- roh_final_rohan[roh_final_rohan$Sample == finalV1[i],]
  data.frame(sample=finalV1[i],
             Small2=sum(df2[df2$CollFly_female>=0 & df2$CollFly_female<2,]$ROH_LENGTH)/1112400596,
             Small=sum(df2[df2$CollFly_female>=2 & df2$CollFly_female<4,]$ROH_LENGTH)/1112400596, 
             Intermediate=sum(df2[df2$CollFly_female>=4 & df2$CollFly_female<8,]$ROH_LENGTH)/1112400596, 
             Intermediate2=sum(df2[df2$CollFly_female>=8 & df2$CollFly_female<16,]$ROH_LENGTH)/1112400596, 
             Long=sum(df2[df2$CollFly_female>=16 & df2$CollFly_female<32,]$ROH_LENGTH)/1112400596, 
             ExtraLong=sum(df2[df2$CollFly_female>=32,]$ROH_LENGTH)/1112400596,
             Average=mean(df2$CollFly_female)
  )
})
roh_perc_df <- do.call(rbind, roh_perc)

roh_perc_df_gather <- gather(roh_perc_df, key, value, -sample,-Average)
roh_perc_df_gather$key <- gsub("Small2","0-2 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Small","2-4 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate2","8-16 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate","4-8 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("ExtraLong",">32 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long","16-32 gen", roh_perc_df_gather$key )

order <- c("2 gen","4 gen","8 gen","16 gen","32 gen",">32 gen")
order <- c(">32 gen","16-32 gen","8-16 gen","4-8 gen","2-4 gen","0-2 gen")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$SampleID <- roh_perc_df_gather$sample

roh_perc_df_gather2 <- merge(roh_perc_df_gather, metadata, by="SampleID")
roh_perc_df_gather2 <- roh_perc_df_gather2[roh_perc_df_gather2$Category!="Historical",]

roh_perc_df_gather2$Category <- factor(roh_perc_df_gather2$Category,
                              levels=order_category, ordered=T)

pa <- ggplot(roh_perc_df_gather2, 
             aes(sample,value*100,fill=key))
pa <- pa + geom_col()+ theme_classic() +  xlab("")+ facet_grid(.~Category,  scales="free")+
  ylab("% of Genome in RoH") + 
  ggtitle("CollFly_female recomb rate")+
  theme(legend.position = "right", axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=rev(c("#8D5B4C","#CDACA2","#587D71","#ABC4BC","#F6AE2D","#FBDA9D")), name="RoH Size") 
pa


pa <- ggplot(roh_perc_df_gather2[roh_perc_df_gather2$Coverage>10&
                                   roh_perc_df_gather2$key!=">32 gen",], 
             
             aes(Category,value*100,color=Category))
pa <- pa + geom_boxplot()+ geom_jitter()+ theme_classic() +  xlab("")+ facet_wrap(.~key,  scales="free", nrow=1)+
  ylab("% of Genome in RoH") + 
  ggtitle("CollFly_female recomb rate")+
  theme(legend.position = "none", axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_color_manual(values=colors_category_2)+
  stat_compare_means(label.x = 1.5)+
  stat_compare_means(comparisons = my_comparisons)

pa

# only to 5 gen and larger than 5
roh_perc <- list()

roh_perc <- lapply(1:length(finalV1), function(i) {
  df2 <- roh_final_rohan[roh_final_rohan$Sample == finalV1[i],]
  data.frame(sample=finalV1[i],
             Small=sum(df2[df2$CollFly_female>=0 & df2$CollFly_female<5,]$ROH_LENGTH)/1112400596,
             Long=sum(df2[df2$CollFly_female>=5,]$ROH_LENGTH)/1112400596, 
             Average=mean(df2$CollFly_female)
  )})
roh_perc_df <- do.call(rbind, roh_perc)

roh_perc_df_gather <- gather(roh_perc_df, key, value, -sample,-Average)
roh_perc_df_gather$key <- gsub("Small","0-5 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long",">=5 gen", roh_perc_df_gather$key )

order <- c("0-5 gen",">=5 gen")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$SampleID <- roh_perc_df_gather$sample

roh_perc_df_gather2 <- merge(roh_perc_df_gather, metadata, by="SampleID")
roh_perc_df_gather2 <- roh_perc_df_gather2[roh_perc_df_gather2$Category!="Historical",]

roh_perc_df_gather2$Category <- factor(roh_perc_df_gather2$Category,
                                       levels=order_category, ordered=T)

pa <- ggplot(roh_perc_df_gather2[which(roh_perc_df_gather2$key%in%c("0-5 gen",">=5 gen")&roh_perc_df_gather2$Coverage>10),], 
             aes(Category,value*100,color=Category))
pa <- pa + geom_boxplot()+ geom_jitter(height=0)+ theme_classic() +  xlab("")+ facet_wrap(.~key,  scales="free", nrow=1)+
  ylab("% of Genome in RoH") + 
  ggtitle("Red Necked Crane recomb rate")+
  theme(legend.position = "none", axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_color_manual(values=colors_category_2)+
  stat_compare_means(label.x = 1.5, label.y=c(21,8))+
  stat_compare_means(comparisons = my_comparisons)

pa


# summary ROH TGen Black Necked Crane ------
roh_perc <- list()
roh_perc <- lapply(1:length(finalV1), function(i) {
  df2 <- roh_final_rohan[roh_final_rohan$Sample == finalV1[i],]
  data.frame(sample=finalV1[i],
             Small2=sum(df2[df2$RedNeckedCrane>=0 & df2$RedNeckedCrane<2,]$ROH_LENGTH)/1112400596,
             Small=sum(df2[df2$RedNeckedCrane>=2 & df2$RedNeckedCrane<4,]$ROH_LENGTH)/1112400596, 
             Intermediate=sum(df2[df2$RedNeckedCrane>=4 & df2$RedNeckedCrane<8,]$ROH_LENGTH)/1112400596, 
             Intermediate2=sum(df2[df2$RedNeckedCrane>=8 & df2$RedNeckedCrane<16,]$ROH_LENGTH)/1112400596, 
         #    Long=sum(df2[df2$RedNeckedCrane>=16 & df2$RedNeckedCrane<32,]$ROH_LENGTH)/1112400596, 
          #   ExtraLong=sum(df2[df2$RedNeckedCrane>=32,]$ROH_LENGTH)/1112400596,
             Average=mean(df2$RedNeckedCrane)
  )
})
roh_perc_df <- do.call(rbind, roh_perc)

roh_perc_df_gather <- gather(roh_perc_df, key, value, -sample,-Average)
roh_perc_df_gather$key <- gsub("Small2","0-2 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Small","2-4 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate2","8-16 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate","4-8 gen", roh_perc_df_gather$key )
#roh_perc_df_gather$key <- gsub("ExtraLong",">32 gen", roh_perc_df_gather$key )
#roh_perc_df_gather$key <- gsub("Long","16-32 gen", roh_perc_df_gather$key )

order <- c("8-16 gen","4-8 gen","2-4 gen","0-2 gen")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$SampleID <- roh_perc_df_gather$sample

roh_perc_df_gather2 <- merge(roh_perc_df_gather, metadata, by="SampleID")
roh_perc_df_gather2 <- roh_perc_df_gather2[roh_perc_df_gather2$Category!="Historical",]

roh_perc_df_gather2$Category <- factor(roh_perc_df_gather2$Category,
                                       levels=order_category, ordered=T)

pa <- ggplot(roh_perc_df_gather2, 
             aes(sample,value*100,fill=key))
pa <- pa + geom_col()+ theme_classic() +  xlab("")+ facet_grid(.~Category,  scales="free")+
  ylab("% of Genome in RoH") + 
  theme(legend.position = "right", legend.title = element_blank(),axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=rev(c("#8D5B4C","#CDACA2","#587D71","#ABC4BC","#F6AE2D","#FBDA9D")), name="RoH Size") 
pa
ggsave("Manuscript/Plots_supplementary/FigureS7_barplot_timing_blacknecked.pdf",  width = 13, height = 4)




pa_rn <- ggplot(roh_perc_df_gather2[which(roh_perc_df_gather2$key%in%c("8-16 gen","4-8 gen","2-4 gen","0-2 gen")&
                                         roh_perc_df_gather2$Coverage>10),], 
             aes(Category,value*100,color=Category))
pa_rn <- pa_rn + geom_boxplot()+ geom_jitter(height=0)+ theme_classic() +  xlab("")+ facet_wrap(.~key,  scales="free", nrow=1)+
  ylab("% of Genome in RoH") + 
 # ggtitle("Red Necked Crane recomb rate")+
  theme(legend.position = "none", axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_color_manual(values=colors_category_2)+
  stat_compare_means(label.x = 1.5)+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")#method = "t.test")

pa_rn


# Summarize it less than 5 generations (65 years) and further away) NOT USED-----
roh_perc <- list()

roh_perc <- lapply(1:length(finalV1), function(i) {
  df2 <- roh_final_rohan[roh_final_rohan$Sample == finalV1[i],]
  data.frame(sample=finalV1[i],
             Small=sum(df2[df2$RedNeckedCrane>=0 & df2$RedNeckedCrane<4,]$ROH_LENGTH)/1112400596,
             Long=sum(df2[df2$RedNeckedCrane>=4,]$ROH_LENGTH)/1112400596, 
             Average=mean(df2$RedNeckedCrane)
  )})
roh_perc_df <- do.call(rbind, roh_perc)

roh_perc_df_gather <- gather(roh_perc_df, key, value, -sample,-Average)
roh_perc_df_gather$key <- gsub("Small","0-4 gen", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long",">=4 gen", roh_perc_df_gather$key )

order <- c("0-4 gen",">=4 gen")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$SampleID <- roh_perc_df_gather$sample

roh_perc_df_gather2 <- merge(roh_perc_df_gather, metadata, by="SampleID")
roh_perc_df_gather2 <- roh_perc_df_gather2[roh_perc_df_gather2$Category!="Historical",]

roh_perc_df_gather2$Category <- factor(roh_perc_df_gather2$Category,
                                       levels=order_category, ordered=T)

pa <- ggplot(roh_perc_df_gather2[which(roh_perc_df_gather2$Coverage>10),], 
             aes(Category,value*100,color=Category))
pa <- pa + geom_boxplot()+ geom_jitter(height=0)+ theme_classic() +  xlab("")+ facet_wrap(.~key,  scales="free", nrow=1)+
  ylab("% of Genome in RoH") + 
  ggtitle("Red Necked Crane recomb rate")+
  theme(legend.position = "none", axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_color_manual(values=colors_category_2)+
  stat_compare_means(label.x = 1.5)+
  stat_compare_means(comparisons = my_comparisons)

pa

roh_perc_df_gather2$Category2 <- roh_perc_df_gather2$Category
roh_perc_df_gather2$Category2 <- gsub("Early_captive", "Captive_born",roh_perc_df_gather2$Category2)
roh_perc_df_gather2$Category2 <- gsub("Late_captive", "Captive_born",roh_perc_df_gather2$Category2)
roh_perc_df_gather2$Category2 <- factor(roh_perc_df_gather2$Category2,
                               levels=c("Founders_wildborn","Captive_born","Wild"), ordered=T)

pa <- ggplot(roh_perc_df_gather2[which(roh_perc_df_gather2$Coverage>10),], 
             aes(Category2,value*100,color=Category2))
pa <- pa + geom_boxplot()+ geom_jitter(height=0)+ theme_classic() +  xlab("")+ facet_wrap(.~key,  scales="free", nrow=1)+
  ylab("% of Genome in RoH") + 
  ggtitle("Red Necked Crane recomb rate")+
  theme(legend.position = "none", axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_color_manual(values=colors_category_2)+
  stat_compare_means(label.x = 1.5)+
  stat_compare_means(comparisons = my_comparisons2)

pa

# Time, generations and size of roh distribution ----
library(data.table)
library(patchwork)

GenTime=13
L=seq(0,20,by = 0.01)
r=c(3.42,2.28)
X=data.table(expand_grid(L=L,r=r))
X$species=NA
X$species[X$r==3.42]='RedNecked_crane'
X$species[X$r==2.28]='CollFly_female'

# formula is L = 100/2t cM (from Thompson 2013, Genetics 194, 301-326)
X$Tgen=(100/(X$L*X$r))/2
X$Tyears=X$Tgen*GenTime
X=X[order(species,L)]

ggplot(X,aes(L,Tyears,col=paste(species,r))) + geom_line() + geom_point()
p1=ggplot(X[L> 0.5 & L < 1],aes(L,Tyears,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH=0.5-1 MB")
p2=ggplot(X[L>1 & L < 2.5],aes(L,Tyears,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH=1-2.5 MB")
p3=ggplot(X[L>2.5 & L < 5],aes(L,Tyears,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH=2.5-5 MB")
p4=ggplot(X[L>5 & L < 10],aes(L,Tyears,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH=5-10 MB")
p5=ggplot(X[L>= 10],aes(L,Tyears,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH>=10 MB")
(p1+p2+p3+p4+p5) +  plot_layout(guides = 'collect')

p1=ggplot(X[L> 0.5 & L < 1],aes(L,Tgen,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH=0.5-1 MB")
p2=ggplot(X[L>1 & L < 2.5],aes(L,Tgen,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH=1-2.5 MB")
p3=ggplot(X[L>2.5 & L < 5],aes(L,Tgen,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH=2.5-5 MB")
p4=ggplot(X[L>5 & L < 10],aes(L,Tgen,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH=5-10 MB")
p5=ggplot(X[L>=10 ],aes(L,Tgen,col=paste(species,r))) + geom_line() + geom_point() + ggtitle("ROH>=10 MB")

(p1+p2+p3+p4+p5) +  plot_layout(guides = 'collect')




# Figure 4E -Adding load in rohs ----
high_list_df <- read.csv2("Files/Genotype_load_modern_10x.csv")
roh_sample <- read_excel("Files/Roh_persample_metadata.xlsx", sheet = 1)
total_hom <-read.table("Files/snpEff/Total_countsHOM.txt") # it is from downsampled at 10x
colnames(total_hom) <- c("SampleID","TotalHOMALT","TotalHOMDERIVED","HOMDERV_ROH")
total_hom$HOMDERV_outROH <- total_hom$TotalHOMDERIVED - total_hom$HOMDERV_ROH

samples_modern <- unique(roh_sample$SampleID)

samples_modern <- samples_modern[-22]

counts_inside_outside <- do.call(rbind,lapply(1:length(samples_modern), function(i){
  load_high <- high_list_df[high_list_df$SampleID==samples_modern[i]&high_list_df$Class=="High",]
  load_moderate <- high_list_df[high_list_df$SampleID==samples_modern[i]&high_list_df$Class=="Moderate",]
  load_low <- high_list_df[high_list_df$SampleID==samples_modern[i]&high_list_df$Class=="Low",]
  
  load_HOM_high <- cbind.data.frame(load_high[load_high$Genotype=="1/1",]$Chrom,
                                    load_high[load_high$Genotype=="1/1",]$Pos-1, load_high[load_high$Genotype=="1/1",]$Pos,
                                    "High")
  load_HOM_moderate <- cbind.data.frame(load_moderate[load_moderate$Genotype=="1/1",]$Chrom,
                                        load_moderate[load_moderate$Genotype=="1/1",]$Pos-1, load_moderate[load_moderate$Genotype=="1/1",]$Pos,
                                        "Moderate")
  load_HOM_low <- cbind.data.frame(load_low[load_low$Genotype=="1/1",]$Chrom,
                                   load_low[load_low$Genotype=="1/1",]$Pos-1, load_low[load_low$Genotype=="1/1",]$Pos,
                                   "Low")
  
  colnames(load_HOM_high) <- c("Scaffold","BEGIN","END","Class")
  colnames(load_HOM_moderate) <- c("Scaffold","BEGIN","END","Class")
  colnames(load_HOM_low) <- c("Scaffold","BEGIN","END","Class")
  
  roh <- roh_sample[roh_sample$SampleID==samples_modern[i],c(2:4)]
  
  load_inRoh_HOM_high <- bedtoolsr::bt.intersect(a = load_HOM_high, b=roh,wa=TRUE) 
  load_inRoh_HOM_moderate <- bedtoolsr::bt.intersect(a = load_HOM_moderate, b=roh,wa=TRUE)
  load_inRoh_HOM_low <- bedtoolsr::bt.intersect(a = load_HOM_low, b=roh,wa=TRUE)
  
  #add load outside ROH
  load_outsideRoh_HOM_high <-  bedtoolsr::bt.subtract(a=load_HOM_high, b=load_inRoh_HOM_high)
  load_outsideRoh_HOM_moderate <-  bedtoolsr::bt.subtract(a=load_HOM_moderate, b=load_inRoh_HOM_moderate)
  load_outsideRoh_HOM_low <-  bedtoolsr::bt.subtract(a=load_HOM_low, b=load_inRoh_HOM_low)
  
  #TOTAL HOM outside and inside roh
  total_hom_sample <- total_hom[total_hom$SampleID==samples_modern[i],]
  hom_roh <- total_hom_sample$HOMDERV_ROH *2
  hom_outroh <- total_hom_sample$HOMDERV_outROH*2
  
  data.frame(SampleID=samples_modern[i], 
             HOM_roh_High=length(load_inRoh_HOM_high$V1)*2,
             HOM_outsideroh_High=length(load_outsideRoh_HOM_high$V1)*2,
             HOM_total_High=length(load_HOM_high$Scaffold)*2,
             
             Prop_HOM_roh_High=(length(load_inRoh_HOM_high$V1)*2)/hom_roh,
             Prop_HOM_outsideroh_High=(length(load_outsideRoh_HOM_high$V1)*2)/hom_outroh, # does this make sense
             
             Prop_HOM_roh_High_LOW=((length(load_inRoh_HOM_high$V1)*2)/length(load_inRoh_HOM_low$V1)*2)/hom_roh,
             Prop_HOM_outsideroh_High_LOW=((length(load_outsideRoh_HOM_high$V1)*2)/length(load_outsideRoh_HOM_low$V1)*2)/hom_outroh, 
             
             HOM_roh_Moderate=length(load_inRoh_HOM_moderate$V1)*2,
             HOM_outsideroh_Moderate=length(load_outsideRoh_HOM_moderate$V1)*2,
             HOM_total_Moderate=length(load_HOM_moderate$Scaffold)*2,
             
             Prop_HOM_roh_Moderate=(length(load_inRoh_HOM_moderate$V1)*2)/hom_roh,
             Prop_HOM_outsideroh_Moderate=(length(load_outsideRoh_HOM_moderate$V1)*2)/hom_outroh, # does this make sense
             Prop_HOM_roh_Moderate_LOW=((length(load_inRoh_HOM_moderate$V1)*2)/length(load_inRoh_HOM_low$V1)*2)/hom_roh,
             Prop_HOM_outsideroh_Moderate_LOW=((length(load_outsideRoh_HOM_moderate$V1)*2)/length(load_outsideRoh_HOM_low$V1)*2)/hom_outroh, 
             
             HOM_roh_Low=length(load_inRoh_HOM_low$V1)*2,
             HOM_outsideroh_Low=length(load_outsideRoh_HOM_low$V1)*2,
             HOM_total_Low=length(load_HOM_low$Scaffold)*2,
             Prop_HOM_roh_Low=(length(load_inRoh_HOM_low$V1)*2)/hom_roh,
             Prop_HOM_outsideroh_Low=(length(load_outsideRoh_HOM_low$V1)*2)/hom_outroh # does this make sense
  )
}))


#  Proportion of load inside and outside roh - corrected by totalHOM (not genome size in roh) 
counts_inside_outside_metadata <- merge(counts_inside_outside, metadata[,c("SampleID","Category","Type","Coverage","FROH 1Mb (All only Coverage >4)")], by="SampleID")
order_category <- c("Founders_wildborn","Early_captive","Late_captive","Wild")
counts_inside_outside_metadata$Category <- factor(counts_inside_outside_metadata$Category, levels=order_category,ordered=TRUE) 



counts_inside_outside_metadata_gather <- gather(counts_inside_outside_metadata, key, value, -Category, -Coverage, -SampleID, -Type,
                                                -HOM_roh_High, -HOM_outsideroh_High, -HOM_total_High, -HOM_roh_Moderate,  -Prop_HOM_roh_High_LOW,
                                                -Prop_HOM_outsideroh_High_LOW,-Prop_HOM_roh_Moderate_LOW,-Prop_HOM_outsideroh_Moderate_LOW,
                                                -HOM_outsideroh_Moderate, -HOM_total_Moderate,
                                                -HOM_roh_Low,-HOM_total_Low,-HOM_outsideroh_Low,-`FROH 1Mb (All only Coverage >4)`)

counts_inside_outside_metadata_gather$LoadType <- unlist(lapply(1:length(counts_inside_outside_metadata_gather$key), function(i) 
  strsplit(counts_inside_outside_metadata_gather$key, "_")[[i]][4]))

counts_inside_outside_metadata_gather$Group <- unlist(lapply(1:length(counts_inside_outside_metadata_gather$key), function(i) 
  strsplit(counts_inside_outside_metadata_gather$key, "_")[[i]][3]))

order_type <- c("Low","Moderate","High")
counts_inside_outside_metadata_gather$LoadType <- factor(counts_inside_outside_metadata_gather$LoadType, levels=order_type,ordered=TRUE) 

counts_inside_outside_metadata_gather$Group <- gsub("outsideroh","Outside ROH",counts_inside_outside_metadata_gather$Group)
counts_inside_outside_metadata_gather$Group <- gsub("roh","Inside ROH",counts_inside_outside_metadata_gather$Group)

g_loadroh <- ggplot(counts_inside_outside_metadata_gather, aes(Group,value, fill=Group))
g_loadroh <- g_loadroh + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size=1, aes(group=Group))+ theme_classic() +
  xlab("") + facet_grid(LoadType~Category, scales="free_y")+ ylab("Normalized Derived Allele Count")+
  scale_fill_manual(values=c("#f5c542","#76b7f5"))+
  stat_compare_means() + theme(axis.text.x = element_text(hjust=1, angle=45), legend.position = "none")
g_loadroh


# supplementary 
g <- ggplot(counts_inside_outside_metadata_gather[which(counts_inside_outside_metadata_gather$Category%in%c("Late_captive","Wild")),],
            aes(Category,value, fill=Category))
g <- g + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size=1, aes(group=Group))+ theme_classic() +
  xlab("") + facet_grid(LoadType~Group, scales="free_y")+ ylab("Normalized Derived Allele Count")+
  scale_fill_manual(values=c("#f5c542","#76b7f5"))+
  stat_compare_means() + theme(axis.text.x = element_text(hjust=1, angle=45), legend.position = "none")
g

ggsave("Manuscript/Plots_supplementary/FigureS12_load_latewild.pdf", height = 8, width = 10)


# Final figure 4 -----
ggarrange(ggarrange(het, froh,froh3, nrow=1, labels=c("A","B","C")),
          pa_rn, g_loadroh,nrow=3, labels=c("","D","E"), heights = c(1,1,2))
ggsave("Manuscript/MainFigures/Figure4.pdf", height = 15, width = 13)


# Relatedness wild supplementary -----
library(corrplot)
library(Matrix)

rel_vcf <- read.table("Files/relatedness/WC_modern_allPos_rel.ml", header = TRUE)
samples <- data.frame()
Names <- read.table("Files/PCA/Samples_modern")$V1
Names_order <- read.table("Files/PCA/Samples_modern")$V1

rel_vcf$a <- rel_vcf$a +1
rel_vcf$b <- rel_vcf$b +1

Names_df <- data.frame(SampleID=Names)
Names_metadata <- merge(Names_df, metadata[,c("SampleID","Category","StudBook Number")], by="SampleID")

Names_metadata <- cbind.data.frame(Names_metadata, a=c(1:37), b=c(1:37))
Names_metadata$IdFinal <- paste0(Names_metadata$SampleID,"_", Names_metadata$Category)


mat <- sparseMatrix(i=rel_vcf$a, j=rel_vcf$b, x=rel_vcf$theta, symmetric = TRUE)
diag(mat) <- 0.5
mat <- as.matrix(mat)

#colnames(mat) <- Names
#row.names(mat) <- Names

colnames(mat) <- Names_metadata$IdFinal
row.names(mat) <- Names_metadata$IdFinal
x <- mat

pdf("Plots/Modern_ngsrelate2.pdf", height=11, width = 11)
corrplot(x,
         tl.col = "black", type="lower",
         tl.cex = 0.95,  
         tl.srt = 45,# method = "cir",
         pch.cex = 1.0, method="square",is.corr = FALSE, )

dev.off()


# average pairwise relatedness within and between groups
rel_vcf2 <- merge(rel_vcf, Names_metadata[,c("a","SampleID","Category","IdFinal","StudBook Number")], by="a")
rel_vcf2 <- merge(rel_vcf2, Names_metadata[,c("b","SampleID","Category","IdFinal","StudBook Number")], by="b")

head(rel_vcf2)
# add the other comparison
rel_vcf2_sub <- rel_vcf2[,c(16,34,35,38,39)]
#write.csv2(rel_vcf2_sub,"Files/relatedness/Relatedness.csv", quote = F, row.names = F) # For supplementary table S3
rel_vcf2_sub2 <- rel_vcf2_sub
colnames(rel_vcf2_sub2) <- c( "theta","SampleID.y", "Category.y", "SampleID.x" ,"Category.x")

rel_vcf2_sub_all <- rbind.data.frame(rel_vcf2_sub,rel_vcf2_sub2)

# average relatedness 
mean(rel_vcf2_sub$theta)
rel_vcf2_sub$theta
min(rel_vcf2_sub$theta)

# Same 0.5 (>0.375)
rel_vcf2_sub[rel_vcf2_sub$theta >=0.375,]

# First degree 0.25 (0.1875-0.375)
rel_vcf2_sub[rel_vcf2_sub$theta >=0.1874&rel_vcf2_sub$theta<0.375,]

# second degree 0.125 (0.09375-0.1875)
rel_vcf2[rel_vcf2$theta >=0.09375&rel_vcf2$theta<0.1875,]

# 3rd degree 0.0625 (0.046875-0.09375)
# forth or more distant <0.046875 [ 0.03125 ]
rel_vcf2[rel_vcf2$theta <0.09375&rel_vcf2$theta>0.001,]


# plot
categories <- c("Founders_wildborn","Early_captive","Late_captive","Wild")

mean_rel <- list()
mean_rel2 <- list()

for(i in 1:length(categories)){
  cat1 <- categories[i]
  for(j in 1:length(categories)){
    cat2 <- categories[j]
    mean_group <- mean(rel_vcf2_sub_all[which(rel_vcf2_sub_all$Category.x==cat1 &rel_vcf2_sub_all$Category.y==cat2),]$theta)
    mean_rel2[[j]] <- data.frame(Category1=cat1, Category2=cat2, MeanRel=mean_group)
  }
  mean_rel[[i]] <- do.call(rbind, mean_rel2)
}
mean_rel_df <- do.call(rbind, mean_rel)


mean_rel_df$Category1 <- factor(mean_rel_df$Category1, levels=categories, ordered = T)
mean_rel_df$Category2 <- factor(mean_rel_df$Category2, levels=categories, ordered = T)

g <- ggplot(mean_rel_df, aes(Category1, Category2, fill=MeanRel))
g <- g + geom_tile() + theme_classic() + scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0)) + theme(axis.line = element_blank())+
  xlab("") + ylab("") +   scale_fill_gradient(low = "#FFFFCC", high = "#8b0000") +
  guides(fill = guide_colourbar(title = "Mean Kinship Coefficient"))+
  geom_text(aes(label = format(round(MeanRel, 4), nsmall =4)), color = "black", size = 4) 
g
ggsave("Manuscript/Plots_supplementary/FigureS7_relatednessgroup.pdf", width = 7.5, height = 4)
