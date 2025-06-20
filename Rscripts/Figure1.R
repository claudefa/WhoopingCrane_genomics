# Figure 1 - Whooping cranes. 
library("ggplot2")
library("tidyr")
library("readxl")
library("ggpubr")
library("ggrepel")
library("rgdal") 
library('maps')
library('mapdata')
library('maptools')
library('png')

setwd("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/")
colors_category <- c("#6d466bff","#e3e359","#EB9486" )
colors_category_2 <- c("#8C2717","#E67260","#F1B1A7","#e3e359" )

metadata <- read_excel("Manuscript/Review/Fontsere_etal_whooping_crane_Supp_tables_R1.xlsx", sheet = 2)
metadata$Type <- factor(metadata$Type, levels = c("Historical","Wild","Captive"), ordered = T)

## Figure 1 A Map -------
map <-readOGR("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/Files/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp")
map2 <-readOGR("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/Files/ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes.shp")
lakes <-readOGR("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/Files/ne_50m_lakes/")
rivers <-readOGR("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/Files/ne_110m_rivers_lake_centerlines/ne_110m_rivers_lake_centerlines.shp")
ocean<-readOGR("~/Documents/OneDrive - University of Copenhagen/CottonTop_Tamarins/Files/ne_10m_ocean/ne_10m_ocean.shp")

# WC
area.1 <- readOGR("Files/redlist_species_data_2af8698a-5d2f-4446-9934-73eecef6dc14/data_0.shp")
area.1.points <- fortify(area.1)
area.all <- rbind.data.frame(area.1.points)

# Plain Map 
k <-ggplot() + geom_polygon(data=ocean,aes(x = long,y = lat,group = group),fill = "#e6f4ff",size=0.3) + 
  geom_path(data=ocean,aes(x = long,y = lat,group = group),colour = "#e6f4ff",size=0.3)

k <- k +  geom_polygon(data = map, aes(x = long,y = lat, group = group),fill = "#fffdf5",size=0.3) +
  geom_path(data = map,aes(x = long,y = lat,group = group),colour = "grey50",size=0.4) + 
  geom_path(data = map2,aes(x = long,y = lat,group = group),colour = "grey50",size=0.2) +
  geom_polygon(data=lakes,aes(x = long,y = lat,group = group),fill = "#e6f4ff",size=0.3)+
  geom_path(data=rivers,aes(x = long,y = lat,group = group),colour = "#e6f4ff",size=0.7)+
  labs(x = "Longitude", y = "Latitude") + coord_fixed(xlim = c(-130, -74), ylim = c(20, 65), ratio=1) + 
  theme(panel.background = element_rect(fill = '#F3F1F1'),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

k

# Subspecies Map
area.all$species <- area.all$group
colors_habitat=rep(c("#7a9e70"),  times=c(15))

x <- k
x <- x + geom_polygon(data = area.all,aes(x = long,y = lat, fill = group, color=group), color=NA, size=0.01, alpha = 0.4)  +
  ylab("") +xlab("")+
  scale_fill_manual(values=colors_habitat) +   
  geom_path(data=rivers,aes(x = long,y = lat,group = group),colour = "#F3F1F1",size=0.4)
x

Site_coordinates <- read_excel("Files/historicalSites.xlsx")

LABELS <- Site_coordinates$Site
LON <- Site_coordinates$Longitud
LAT <- Site_coordinates$Latitude
SIZE <- Site_coordinates$Latitude
TYPE <- Site_coordinates$Type

coord.SITES <- cbind.data.frame(LABELS, LON, LAT,SIZE,TYPE)
coord.SITES$LABELS <- as.factor(coord.SITES$LABELS)
coord.SITES$LAT <- as.numeric(coord.SITES$LAT)
coord.SITES$LON <- as.numeric(coord.SITES$LON)
coord.SITES$SIZE <- as.numeric(coord.SITES$SIZE)
coord.SITES$TYPE <- as.factor(coord.SITES$TYPE)

colors <- c("#6d466bff","#e3e359","#eb9486ff" )
coord.SITES$TYPE <- factor(coord.SITES$TYPE, levels = c("Historical","Wild","Captive"), ordered = T)
s <- x
s <- s +
  geom_point(data=coord.SITES,aes(LON, LAT,color=TYPE, group=TYPE),size=4) + ylab("") +xlab("")+
  geom_text_repel(data=coord.SITES,aes(LON, LAT, label=LABELS),color="gray2", size=4, max.overlaps = Inf, box.padding = 0.7)+
  scale_color_manual(values=colors)+
  scale_size(range = c(0, 10))+
  guides(fill="none", shape="none")+
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA, size=0.8), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.box.background = element_blank(),
        legend.background = element_blank(),legend.title =element_blank(),
        legend.key=element_rect(fill="white"), legend.text = element_text(size=10),
        legend.position = c(0.15,0.2),
        legend.box = "horizontal")
s

## Figure 1 B Heterozygosity many Birds ----
metadata <- read_excel("Metadata_listBirds.xlsx", sheet = 1)
fins <- list.files("Files/Het/AllBirds/",pattern=".ml",full.names=TRUE)
samples <- sub('\\_het.ml$', '', basename(fins))

df_list <- lapply(fins,
                  FUN = function(files) {
                    scan(files)
                  })
df_list2 <- list()
for (i in 1:length(samples)){
  
  df_list2[[i]] <- cbind.data.frame(Sample=samples[i], Het=(df_list[[i]][2]/sum(df_list[[i]])))
}

het_df <- do.call(rbind,df_list2)
het_df_metadata <- merge(het_df, metadata, by="Sample")
het_df_metadata$Sample <- gsub("_"," ", het_df_metadata$Sample)

het_df_metadata$IUCN <- factor(het_df_metadata$IUCN, levels = c( "Least Concern","Near Threatened",       
                                                                 "Vulnerable" ,"Endangered" ,                  
                                                                 "Critically Endangered"), ordered = T )

#Add values for sample with highest coverage historical (MCZ-220028 , 7x coverage =0.00102421148792956*2.89) and average modern samples at 10x (0.000782119)
#ts/tv ratio 2.89

wc_dataset <- data.frame(Sample=c("Grus_americana_historical_8x","Grus_americana_modern_10x"),
                         Het=c(0.00102421148792956*2.89,0.000782119),
                         Order=c("Gruiformes","Gruiformes"),
                         IUCN=c("Endangered","Endangered"))
het_df_metadata2 <- rbind.data.frame(het_df_metadata,wc_dataset)

het_df_metadata2 <- het_df_metadata2[with(het_df_metadata2,order(het_df_metadata2$Het)), ] ## Sorting
samples_order <- het_df_metadata2$Sample
het_df_metadata2$Sample <- factor(het_df_metadata2$Sample, levels=samples_order, ordered=TRUE)

het_df_metadata2$Type <- c(rep("Others", times=6),"Sp","Sp",rep("Others", times=15), "Sp",rep("Others", times=12))

g <- ggplot(het_df_metadata2, aes(Sample, Het, fill=IUCN, shape=Type))
g <- g+ theme_classic()+ 
  geom_segment(aes(x=Sample, xend=Sample, y=0, yend=Het, color=Type), linetype = "dotted", alpha=0.8) +
  scale_color_manual(values=c("grey2","#ff333a"))+
  geom_point(size=4,color="grey2")+ #facet_grid(.~Order, scales="free_x")+
  scale_shape_manual(values=c(21,23))+  guides(fill = guide_legend(override.aes = list(shape = 21), nrow=2, byrow=T), shape="none", color="none")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.007))+
  theme(legend.title = element_blank(), legend.position = "top",panel.grid.major.y=element_line(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(margin=margin(t = 0, r = 20, b = 0, l = 0),
                                   angle=90, hjust=1,vjust=0.5,
                                   face=c("italic"),
                                   color=c("black","black","black","black","black",
                                           "black","#ff333a","#ff333a","black","black",
                                           "black","black","black","black","black",
                                           "black","black","black","black","black",
                                           "black","black","black","#ff333a","black",
                                           "black","black","black","black","black",
                                           "black","black"))) +
  xlab("") + ylab("Heterozygosity (bp-1)")+
  scale_fill_manual(values = c("#384E77","#8BBEB2","#e6f9af","#F7B801","#FF5A5F")) #+  ggtitle("Heterozygosity")
g


## Figure 1 C Heterozygosity only transversions historical, captive and wild----
colors_category <- c("#6d466bff","#e3e359","#EB9486" )
colors_category_2 <- c("#8C2717","#E67260","#F1B1A7","#e3e359" )

metadata <- read_excel("Manuscript/Review/Fontsere_etal_whooping_crane_Supp_tables_R1.xlsx", sheet = 2)
metadata$Type <- factor(metadata$Type, levels = c("Historical","Wild","Captive"), ordered = T)


my_comparisons <- list( c("Historical", "Wild"), 
                        c("Historical", "Captive"), 
                        c("Wild", "Captive")) 



p <- ggplot(metadata, aes(Type, `Heterozygosity (tv)`, color=Type))
p <- p+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ 
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("")+
  scale_color_manual(values = colors_category)+ ylim(0,0.0017)+
  #stat_compare_means(label.x = 1.8, label.y = 0.0015)+
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test")+
  theme(legend.position = "none", legend.title = element_blank()) #+ ggtitle("Heterozygosity, only Transversions")
p


## Figure 1 D FROH historical, captive and wild----

roh_perc_df <- read.csv2("Files/Froh_total.csv")

roh_perc_df_gather <- gather(roh_perc_df, key, value, -sample, -Median,-Average,Average2, -count,-sumMb,-maxMb)
#roh_perc_df_gather$key <- gsub("Small","500Kb-1Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate2","2.5-5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate","1-2.5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("ExtraLong",">10Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long","5-10Mb", roh_perc_df_gather$key )
#order <- c("Total","500Kb-1Mb", "1-2.5Mb","2.5-5Mb", "5-10Mb",">10Mb")
order <- c("Total","Average2","1-2.5Mb","2.5-5Mb", "5-10Mb",">10Mb")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)
roh_perc_df_gather$SampleID <- roh_perc_df_gather$sample

roh_perc_df_gather2 <- merge(roh_perc_df_gather, metadata, by="SampleID")


# total roh
roh_perc_df_gather4 <- roh_perc_df_gather2[which(roh_perc_df_gather2$key=="Total"),]
order_category <- c("Historical","Founders_wildborn","Early_captive",
                    "Late_captive", "Wild")
roh_perc_df_gather4$Category <- factor(roh_perc_df_gather4$Category,
                                       levels=order_category, ordered=T)

roh_perc_df_gather4$Type <- factor(roh_perc_df_gather4$Type, levels = c("Historical","Wild","Captive"), ordered = T)

p2 <- ggplot(roh_perc_df_gather4[roh_perc_df_gather4$Coverage>=4,], aes(Type, value*100, color=Type))
p2 <- p2+ geom_boxplot(outlier.shape = NA) +  geom_jitter()+   guides(color=guide_legend(keywidth = 2, keyheight = 1,
                                                                                         override.aes = list(size = 1, shape = NA)))+
  theme_classic() + ylab("FROH (%)")+ xlab("")+
  theme(legend.position = "none", legend.title = element_blank()) + 
  scale_color_manual(values=colors_category) +
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test")
 # stat_compare_means(label.x = 1.8, label.y = 24)
p2

# Final figure with map and diagram -----
library(png)
s2 <-  readPNG("Manuscript/MainFigures/Figure1b_2.png")
s2 <-  readPNG("Manuscript/MainFigures/Figure1b_R1.png")
s2_g <- ggplot() + background_image(s2) +
  theme(plot.margin = margin(t=3, l=0, r=0, b=2, unit = "cm"))

ggarrange(ggarrange(s,s2_g , ncol=2,nrow=1, labels=c("A","B")),
          ggarrange(g, ggarrange(p,p2,ncol=1, nrow=2, labels=c("D","E")), 
          ncol=2, labels=c("C","")), nrow=2)
  
ggsave("Manuscript/MainFigures/Figure1_option2_R1.pdf", height = 10, width = 10)



# Supplementary figures -----

# Correlation between heterozygosity and coverage-----

summary(lm(metadata[metadata$Category!="Historical",]$`Heterozygosity (tv)`~metadata[metadata$Category!="Historical",]$Coverage))

p <- ggplot()
p <- p+ geom_point(data=metadata[metadata$Category!="Historical",],
                   aes(`Heterozygosity (tv)`, Coverage, color=Type)) +
  theme_classic() + xlab("Heterozygosity (bp-1)")+
  geom_smooth(data=metadata[metadata$Category!="Historical",],aes(`Heterozygosity (tv)`, Coverage), method = "lm")+
  theme(legend.position = c(0.1,0.8), legend.title = element_blank()) + scale_color_manual(values=colors_category[2:3]) + 
  ggtitle("Modern: Heterozygosity vs Coverage")
p

summary(lm(metadata[metadata$Category=="Historical",]$`Heterozygosity (tv)`~metadata[metadata$Category=="Historical",]$Coverage))

p1 <- ggplot()
p1 <- p1+ geom_point(data=metadata[metadata$Category=="Historical",],
                     aes(`Heterozygosity (tv)`, Coverage, color=Type)) +
  theme_classic() + xlab("Heterozygosity (bp-1)")+
  geom_smooth(data=metadata[metadata$Category=="Historical",],aes(`Heterozygosity (tv)`, Coverage), method = "lm")+
  theme(legend.position = c(0.1,0.8), legend.title = element_blank()) + scale_color_manual(values=colors_category) + 
  ggtitle("Historical: Heterozygosity vs Coverage")
p1

ggarrange(p,p1, ncol=2, labels=c("A","B"))
ggsave("Manuscript/Plots_supplementary/FigureS1_covhet_R1.pdf", height = 4,width = 11)


# Heterozygosity downsampled to 4 x------

p <- ggplot(metadata, aes(Type, `Heterozygosity (tv - 4x)`, color=Type))
p <- p+ geom_boxplot(outlier.shape = NA) + geom_jitter()+ 
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("")+
  scale_color_manual(values = colors_category)+ ylim(0,0.001)+
  stat_compare_means(label.x = 1.8, label.y = 0.001)+
  theme(legend.position = c(0.85,0.8), legend.title = element_blank()) #+ ggtitle("Heterozygosity, only Transversions")
p
ggsave("Manuscript/Plots_supplementary/FigureS2_hetdown4x_R1.pdf", height = 4,width = 6)

# Correlation between heterozygosity and coverage after downsampling not included as supplementary figure
summary(lm(metadata[metadata$Category!="Historical",]$`Heterozygosity (tv - 4x)`~metadata[metadata$Category!="Historical",]$CoverageDown4x))

p <- ggplot()
p <- p+ geom_point(data=metadata[metadata$Category!="Historical",],aes(`Heterozygosity (tv - 4x)`, CoverageDown4x, color=Type)) +
  theme_classic() + xlab("Heterozygosity (bp-1)")+ ylab("Coverage (downsampled 4x)")+
  geom_smooth(data=metadata[metadata$Category!="Historical",],
              aes(`Heterozygosity (tv - 4x)`, CoverageDown4x), method = "lm")+
  theme(legend.position = "right") + scale_color_manual(values=colors_category[2:3]) + 
  ggtitle("Heterozygosity vs Coverage - Down 4x Coverage Modern Samples")
p

summary(lm(metadata[metadata$Category=="Historical",]$`Heterozygosity (tv - 4x)`~metadata[metadata$Category=="Historical",]$CoverageDown4x))

p2 <- ggplot()
p2 <- p2+ geom_point(data=metadata[metadata$Category=="Historical",],aes(`Heterozygosity (tv - 4x)`, CoverageDown4x, color=Type)) +
  theme_classic() + xlab("Heterozygosity (bp-1)")+ ylab("Coverage (downsampled 4x)")+
  geom_smooth(data=metadata[metadata$Category=="Historical",],
              aes(`Heterozygosity (tv - 4x)`, CoverageDown4x), method = "lm")+
  theme(legend.position = "right") + scale_color_manual(values=colors_category[1]) + 
  ggtitle("Heterozygosity vs Coverage - Down 4x Coverage Modern Samples")
p2

# After reviewers comment add more context on the downsampling to 4x so there is no correlation between groups
metadata <- read_excel("Manuscript/Review/Het_coverage_year.xlsx")

metadata$Year <- as.numeric(metadata$Year)

summary(lm(metadata[metadata$Type=="Historical",]$`Heterozygosity (tv - 4x)`~
             metadata[metadata$Type=="Historical",]$Year))

p1 <- ggplot()
p1 <- p1+ geom_point(data=metadata[metadata$Type=="Historical",],aes(Year, `Heterozygosity (tv - 4x)`, color=Type)) +
  theme_classic() + ylab("Heterozygosity (bp-1) - Down 4x")+ xlab("Year")+
  geom_smooth(data=metadata[metadata$Type=="Historical",],
              aes(Year, `Heterozygosity (tv - 4x)`, color=Type), method = "lm")+
  theme(legend.position = "none") + scale_color_manual(values=colors_category) + 
  ggtitle("Heterozygosity Downsampled to 4x vs Year")
p1
summary(lm(metadata[metadata$Type=="Historical",]$`Heterozygosity (tv)`~
             metadata[metadata$Type=="Historical",]$Year))
p2 <- ggplot()
p2 <- p2+ geom_point(data=metadata[metadata$Type=="Historical",],aes(Year, `Heterozygosity (tv)`, color=Type)) +
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("Year")+
  geom_smooth(data=metadata[metadata$Type=="Historical",],
              aes(Year, `Heterozygosity (tv)`, color=Type), method = "lm")+
  theme(legend.position = "none") + scale_color_manual(values=colors_category) + 
  ggtitle("Heterozygosity vs Year")
p2

summary(lm(metadata[metadata$Type=="Captive",]$`Heterozygosity (tv - 4x)`~
             metadata[metadata$Type=="Captive",]$Year))
summary(lm(metadata[metadata$Type=="Wild",]$`Heterozygosity (tv - 4x)`~
             metadata[metadata$Type=="Wild",]$Year))
p3 <- ggplot()
p3 <- p3+ geom_point(data=metadata[metadata$Type!="Historical",],aes(Year, `Heterozygosity (tv - 4x)`, color=Type)) +
  theme_classic() + ylab("Heterozygosity (bp-1) - Down 4x")+ xlab("Year")+
  geom_smooth(data=metadata[metadata$Type!="Historical",],
              aes(Year, `Heterozygosity (tv - 4x)`, color=Type), method = "lm")+
  theme(legend.position = "top") + scale_color_manual(values=colors_category[2:3]) + 
  ggtitle("Heterozygosity Downsampled to 4x vs Year")
p3
summary(lm(metadata[metadata$Type=="Captive",]$`Heterozygosity (tv)`~
             metadata[metadata$Type=="Captive",]$Year))
summary(lm(metadata[metadata$Type=="Wild",]$`Heterozygosity (tv)`~
             metadata[metadata$Type=="Wild",]$Year))
p4 <- ggplot()
p4 <- p4+ geom_point(data=metadata[metadata$Type!="Historical",],aes(Year, `Heterozygosity (tv)`, color=Type)) +
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("Year")+
  geom_smooth(data=metadata[metadata$Type!="Historical",],
              aes(Year, `Heterozygosity (tv)`, color=Type), method = "lm")+
  theme(legend.position = "top") + scale_color_manual(values=colors_category[2:3]) + 
  ggtitle("Heterozygosity vs Year")
p4

ggarrange(p1,p3,p2,p4)
ggsave("Manuscript/Plots_supplementary/Revision_figureHetYearDown.pdf", height = 7, width = 10)


p4 <- ggplot()
p4 <- p4+ geom_boxplot(data=metadata,aes(Category, Coverage, color=Type)) +
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("Year")+
  theme(legend.position = "top") + scale_color_manual(values=colors_category) + 
  ggtitle("Heterozygosity vs Year")
p4


# captive and wild with heterozygosity
p4 <- ggplot()
p4 <- p4+ geom_point(data=metadata[metadata$Type!="Historical",],aes(Year, `Heterozygosity (tv)`, color=Category)) +
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("Year")+
  geom_smooth(data=metadata[metadata$Type!="Historical",],
              aes(Year, `Heterozygosity (tv - 4x)`, color=Category), method = "lm")+
  theme(legend.position = "top") + scale_color_manual(values=colors_category_2) + 
  ggtitle("Heterozygosity vs Year")
p4


# Roh Bar plot
roh_perc_df_gather3 <- roh_perc_df_gather2[which(roh_perc_df_gather2$key!="Total"&
                                                   roh_perc_df_gather2$key!="Average2" ),]
order_category <- c("Historical","Founders_wildborn","Early_captive",
                    "Late_captive", "Wild")
roh_perc_df_gather3$Category <- factor(roh_perc_df_gather3$Category,
                                       levels=order_category, ordered=T)
pa <- ggplot(roh_perc_df_gather3[roh_perc_df_gather3$Coverage>=4,], 
             aes(sample,value,fill=key))
pa <- pa + geom_col()+ theme_classic() +  xlab("")+ facet_grid(.~Category,  scales="free")+
  ylab("FROH") + 
  theme(legend.position = "right", axis.text.y = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=rev(c("#957583","#356255","#AAC0AF","#E8EEE9")), name="RoH Size") 
pa

ggsave("Manuscript/Plots_supplementary/FigureS3_barplotROH.pdf", height = 4,width = 13)



# Correlation between heterozygosity and coverage after downsampling to 10x
summary(lm(metadata[metadata$Category!="Historical",]$Coverage~metadata[metadata$Category!="Historical",]$`Heterozygosity (all Pos only Modern 10x)`))

p <- ggplot()
p <- p+ geom_point(data=metadata[metadata$Category!="Historical",],aes(Coverage, `Heterozygosity (all Pos only Modern 10x)`, , color=Type)) +
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("Coverage (pre-downsampled 10x)")+
  geom_smooth(data=metadata[metadata$Category!="Historical",],
              aes(Coverage, `Heterozygosity (all Pos only Modern 10x)`), method = "lm")+
  theme(legend.position = c(0.1,0.8),legend.title = element_blank()) + scale_color_manual(values=colors_category[2:3]) 
p

summary(lm(metadata[metadata$Category!="Historical",]$Coverage~metadata[metadata$Category!="Historical",]$`Heterozygosity (all Pos only Modern)`))

p2 <- ggplot()
p2 <- p2+ geom_point(data=metadata[metadata$Category!="Historical",],aes(Coverage, `Heterozygosity (all Pos only Modern)`, color=Type)) +
  theme_classic() + ylab("Heterozygosity (bp-1)")+ xlab("Coverage")+
  geom_smooth(data=metadata[metadata$Category!="Historical",],aes(Coverage, `Heterozygosity (all Pos only Modern)`), method = "lm")+
  theme(legend.position =  c(0.1,0.8),legend.title = element_blank()) + scale_color_manual(values=colors_category[2:3]) 
p2
ggarrange(p2,p, labels=c("A","B"))

ggsave("Manuscript/Plots_supplementary/FigureS10_hetdown10x_modern_R1.pdf", height = 4,width = 12)


