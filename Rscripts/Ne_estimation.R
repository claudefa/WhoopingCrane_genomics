#Estimation of Ne
library(ggplot2)
library(readxl)

setwd("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/")

# Variance Ne ----
freq_latecaptivewild <- read.table("Files/NeEstimator/WC_latecaptivewild_flw.frq.counts",header=T)
freq_founders_l <- read.table("Files/NeEstimator/WC_founders_flw.frq.counts",header=T)


foundersN=6
lateN=6
wildN=19
latewildN=lateN+wildN

# latecaptive-Wild vs founder
freq_latecaptive_allele <- cbind.data.frame(SNP=paste0(freq_latecaptivewild$CHR,"_",rownames(freq_latecaptivewild)),
                                     p2=freq_latecaptivewild$C1/((latewildN-freq_latecaptivewild$G0)*2))
freq_founders_l_allele <- cbind.data.frame(SNP=paste0(freq_founders_l$CHR,"_",rownames(freq_founders_l)),
                                           p1=freq_founders_l$C1/((foundersN-freq_founders_l$G0)*2))

freq_l_f <- cbind.data.frame(SNP=freq_latecaptive_allele$SNP, 
                             p1=freq_founders_l_allele$p1,
                             p2=freq_latecaptive_allele$p2) 

# Calculate mean allele frequency per SNP 
freq_l_f$p_bar <- (freq_l_f$p1 + freq_l_f$p2) / 2

df <- freq_l_f

# Remove monomorphic or invalid sites 
df <- df[df$p_bar > 0 & df$p_bar < 1, ]

# Compute per-locus F 
df$Fi <- (df$p1 - df$p2)^2 / (df$p_bar * (1 - df$p_bar))

# Compute average F 
Fhat <- mean(df$Fi, na.rm = TRUE)

# Harmonic mean of sample sizes in chromosomes 
n1 <- foundersN 
n2 <- latewildN 


S <- (2 * n1 * 2 * n2) / (2 * n1 + 2 * n2)

# Estimate variance Ne 
t <- 4  
Ne <- t / (2 * (Fhat - (1 / S)))
Ne

Variance_Ne_Late_Founders <- Ne


# Inbreeding Ne ----
metadata <- read_excel("Manuscript/Review/Fontsere_etal_whooping_crane_Supp_tables_R1.xlsx", sheet = 2)

metadata$FROH

Fwild <- mean(metadata[metadata$Category=="Wild",]$FROH)
Flate_captive<- mean(metadata[metadata$Category=="Late_captive",]$FROH)
Ffounder <- mean(metadata[metadata$Category=="Founders_wildborn",]$FROH)
Flate_captive_wild <- mean(metadata[metadata$Category=="Late_captive"|metadata$Category=="Wild" ,]$FROH)


Delta_l_f <- (Flate_captive - Ffounder)/4
Delta_w_f <- (Fwild - Ffounder)/4 # generations estimated by collection year - gen time of 13y
Delta_wl_f <- (Flate_captive_wild - Ffounder)/4 # generations estimated by collection year - gen time of 13y

#only interested in the late - founders and wild-founders
comparisons <- c("l_f","w_f","lw_F")
delta <- c(Delta_l_f,Delta_w_f,Delta_wl_f)
Inbreeding_Ne=unlist(lapply(1:length(delta), function(i){ 1/(2*(delta[i])) }))