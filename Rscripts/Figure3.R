# Figure 2 - Whooping cranes demographic history
library("ggplot2")
library("readxl")
library("scales")
library("psych")
library("ggpubr")

setwd("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/")


# PSMC plotting--------
# Necessary functions to plot
read_psmc <- function(x){
  #x<- "Files/PSMC/Crane1.combined.psmc"
  out <- scan(x, what = "", sep = "\n", quiet = TRUE)
  niters <- as.integer(gsub("^.*n_iterations:|,.*$", "", out[14]))
  parapattern <- gsub("^.*pattern:|,.*$", "", out[13])
  nintervs <- eval(parse(text = parapattern))
  decoding <- as.integer(gsub("^.*is_decoding:", "", out[15]))
  res <- list()
  res$niters <- niters
  res$n <- nintervs
  res <- extract_psmc(out, res, decoding)
  class(res) <- "psmc"
  res
}

extract_psmc <- function(out, res, decoding){
  res$n_free_lambdas <- as.integer(gsub("^.*, n_free_lambdas:", "", out[13]))
  
  ## get the log-likelihoods:
  s <- grep("^LK", out, value = TRUE)
  res$logLik <- as.numeric(gsub("^LK\t", "", s))
  
  ## get Q differences:
  s <- grep("^QD", out, value = TRUE)
  s <- gsub("^QD\t", "", s)
  s <- unlist(strsplit(s, " -> "))
  s <- matrix(as.numeric(s), ncol = 2, byrow = TRUE)
  colnames(s) <- c("before", "after")
  res$EMQ <- s
  
  ## get RI (relative information or KL distance):
  s <- grep("^RI", out, value = TRUE)
  res$RI <- as.numeric(gsub("^RI\t", "", s))
  
  ## get C_pi and n_recomb
  s <- grep("^MM\tC_pi: ", out, value = TRUE)
  s <- gsub("^MM\tC_pi: ", "", s)
  s <- unlist(strsplit(s, ", n_recomb: "))
  s <- matrix(as.numeric(s), ncol = 2, byrow = TRUE)
  res$Cpi <- s[, 1]
  res$Nrecomb <- s[, 2]
  
  ## get theta, rho, and max T:
  s <- grep("^TR", out, value = TRUE)
  s <- gsub("^TR\t", "", s)
  s <- as.numeric(unlist(strsplit(s, "\t")))
  s <- matrix(s, ncol = 2, byrow = TRUE)
  ## append max T:
  maxt <- grep("^MT", out, value = TRUE)
  maxt <- as.numeric(gsub("^MT\t", "", maxt))
  s <- cbind(s, maxt)
  colnames(s) <- c("theta0", "rho0", "maxT")
  res$parameters <- s
  
  s <- grep("^RS", out, value = TRUE)
  s <- gsub("^RS\t", "", s)
  s <- strsplit(s, "\t")
  s <- matrix(as.numeric(unlist(s)), ncol = 6, byrow = TRUE)
  colnames(s) <- c("k", "t_k", "lambda_k", "pi_k", "sum_A_kl", "A_kk")
  iter <- rep(0:res$niters, each = res$n)
  s <- cbind(s, iter = iter)
  res$RS <- s
  
  s <- grep("^TC", out, value = TRUE)
  s <- gsub("^TC\t", "", s)
  s <- strsplit(s, "\t")
  res$TC <- matrix(as.numeric(unlist(s)), ncol = 4, byrow = TRUE)
  
  if (decoding) {
    s <- grep("^DC", out, value = TRUE)
    s <- gsub("^DC\t", "", s)
    s <- strsplit(s, "\t")
    ns <- length(s)
    s <- as.data.frame(matrix(unlist(s), ns, 6, byrow = TRUE),
                       stringsAsFactors = FALSE)
    for (j in 2:6) s[[j]] <- as.numeric(s[[j]])
    names(s) <- c("Chromosome", "begin", "end", "best-k", "t_k+Delta_k", "max-prob")
    res$decoding <- s
  }
  res
}

# Format data for plotting
getXYplot.psmc <- function(x, mutation.rate, g, scaled, bin.size){
  RS <- x$RS[x$RS[, "iter"] == x$niters, ]
  theta0 <- x$parameters[nrow(x$parameters), "theta0"]
  if (scaled) {
    ##xx <- RS[, "t_k"] / (theta0 / bin.size)
    xx <- RS[, "t_k"] / (theta0 * bin.size)
    yy <- theta0 * RS[, "lambda_k"] / bin.size
  } else {
    # denom <- 4 * mutation.rate * g * bin.size #from psmcr --> bad!
    #  N0 <- theta0 / denom
    # xx <- 2 * N0 * RS[, "t_k"]
    
    N0 <- theta0 / (4*mutation.rate * bin.size)
    xx <- 2 * N0 * RS[, "t_k"] *g
    yy <- N0 * RS[, "lambda_k"]
  }
  list(xx = xx, yy = yy)
}
getXYplot_bt.psmc <- function(x, mutation.rate, g, bin.size){
  boostraprep <- 1:101 # 98 with georgettes files it should be 101 but there were 3 faulty files
  RS <- x$RS[x$RS[, "iter"] == x$niters, ]
  RS <- cbind(RS, rep(boostraprep, each=x$n))
  colnames(RS) <- c("k","t_k","lambda_k","pi_k","sum_A_kl","A_kk","iter","bootstrap")
  RS <- as.data.frame(RS)
  theta <- x$parameters
  theta <- cbind(theta, rep(boostraprep, each=x$niters+1))
  colnames(theta) <- c("theta0","rho0","maxT","boostrap")
  theta <- as.data.frame(theta)
  theta0 <- do.call(rbind,lapply(boostraprep, function(i) theta[theta$boostrap==i,][x$niters+1,]))
  
  result_bt <- list()
  for (i in boostraprep) {
    theta0_i <- theta0[theta0$boostrap==i,]$theta0
    RS_i <- RS[RS$bootstrap==i,]
    N0 <- theta0_i / (4*mutation.rate * bin.size)
    xx <- 2 * N0 * RS_i$t_k *g
    yy <- N0 * RS_i$lambda_k
    result_bt[[i]]<- cbind.data.frame(xx = xx, yy = yy, b=as.character(i))
  }
  do.call(rbind,result_bt)
}

# Read files boostratp
psmc1 <- read_psmc("Files/PSMC/Crane1.combined_new.psmc") # georgettes
psmc1 <- read_psmc("Files/PSMC/WC_standard.split.psmc")  
# plot
mut <- 1.45e-8
gen <- 13

df <- as.data.frame(getXYplot_bt.psmc(psmc1, mut, gen, 100))

df2<- df[which(df$b=="1"&df$xx>=5000&df$xx<=1000000),]

ne_summary <- data.frame(Sample="Whooping Crane",
                         harmonic_mean =harmonic.mean(df2$yy),
                         max_NE=max(df2$yy),
                         min_NE=min(df2$yy),
                         delta_NE=max(df2$yy)-min(df2$yy),
                         ratio_NE=max(df2$yy)/min(df2$yy))

ne_summary

g1 <- ggplot() 
g1 <- g1 +  
  #geom_rect(aes(xmin = 20000, ymin =0, xmax = 26000, ymax = 3.5), fill="#abd5ed", alpha=0.2)+
  geom_vline(xintercept = 2527990, linetype="dotted")+
  geom_step(data = df[df$b!="1",],aes(x = (xx)/1000, y = yy, group=b), color="#e8dfc5") +
  geom_step(data = df[df$b=="1",],aes(x = (xx)/1000, y = yy), color="#cf9b00") +
  xlab(paste0("Thousands of Years ago (kya)")) + ylab("Ne")+
  scale_x_log10(limits=c(5,1000),n.breaks = 10,
                labels = scales::comma,expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,35000),n.breaks = 10,labels = scales::comma)+
  theme_classic()
g1
## Stairway Plot ----
stairwayplot <- read.table("Files/StairWayPlot/Folded/WC_all_morebreaks_nosingletons.summary", header = T) # generation time 13y, all indv, no singletons
#stairwayplot <- read.table("Files/StairWayPlot/Folded/WC_wild_morebreaks_nosingletons.summary", header = T) # generation time 13y, wild indv, no singletons

#Plot

g <- ggplot(stairwayplot) 
g <- g +  geom_step(aes(x = year, y = Ne_median), color="#cf9b00") +
  geom_step(aes(x = year, y = Ne_2.5.), color="#e8dfc5") +
  geom_step(aes(x = year, y = Ne_97.5.), color="#e8dfc5") +
  xlab(paste0("Years ago")) + ylab("Ne")+
  scale_x_log10(limits=c(100,5000), n.breaks = 10,
                labels = scales::comma,expand = c(0, 0) ) +
  scale_y_log10(limits=c(30,50000),n.breaks = 10,labels = scales::comma) +
  
  theme_classic()
g
#ggsave("Plots/Stairwayplot_all_morebreaks_nosingletons.pdf", width = 9, height = 9)
##ggsave("Plots/Stairwayplot_wild_morebreaks_nosingletons.pdf", width = 9, height = 9)


## Gone ----
#gone <- read.csv2("Files/Gone/Output_Ne_WC_wild_2.28cMMb", header=T, sep="\t")
gone <- read.csv2("Files/Gone/Output_Ne_WC_modern_3.42", header=T, sep="\t")

gonebs_list <- list.files("Files/Gone/boostrap/",full.names=TRUE) # boostrap taking one sample out
samples <- gsub('Output_Ne_WC_modern_', '', basename(gonebs_list))

df_gonebs_list <- lapply(gonebs_list,
                  FUN = function(files) {
                    read.csv2(files, sep="\t")
                  })
df_list2 <- list()

for (i in 1:length(samples)){
  
  df_list2[[i]] <- cbind.data.frame(Sample=samples[i], df_gonebs_list[[i]])
}

gonebs <- do.call(rbind,df_list2)
head(gonebs)

# Convert to years
gone$Years <- gone$Generation * 13
gone$Geometric_mean <- as.numeric(gone$Geometric_mean)

gonebs$Years <- gonebs$Generation *13
gonebs$Geometric_mean <- as.numeric(gonebs$Geometric_mean)

# limit 100 generations or 1300ya, also add band duringk 19th century 
min19=2020-1800
max=2020-1938

gone$Years2 <- 2020-gone$Years

g2 <- ggplot() 
g2 <- g2 + 
#  geom_rect(aes(xmin = min19, ymin =0, xmax = max, ymax = 11000), fill="#abd5ed", alpha=0.2)+
  geom_rect(aes(xmin = 140, ymin =-5, xmax = 10, ymax = 500), color="red", fill="white",alpha=0.2)+
  geom_line(data=gonebs,aes(x = Years, y = Geometric_mean, group=Sample), color="#e8dfc5") +
  geom_line(data=gone,aes(x = Years, y = Geometric_mean), color="#cf9b00") +
  scale_x_continuous() +
  xlab(paste0("Years ago")) + ylab("Ne")+
  scale_x_continuous(limits=c(0,500),
                     n.breaks = 10,labels = scales::comma,expand=c(0,0))+
  theme_classic()
g2


g3 <- ggplot() 
g3 <- g3 + 
 # geom_rect(aes(xmin = 135, ymin =0, xmax = 11, ymax = 200), color="red", fill="white")+
  geom_line(data=gonebs,aes(x = Years, y = Geometric_mean, group=Sample), color="#e8dfc5") +
  geom_line(data=gone,aes(x = Years, y = Geometric_mean), color="#cf9b00") +
  scale_y_continuous(limits=c(0,200)) +
  xlab(paste0("Years ago")) + ylab("Ne")+
  scale_x_continuous(sec.axis = sec_axis(~ .-2020,),
                     limits=c(10,140),n.breaks = 10,labels = scales::comma,
                     expand=c(0,0))+
  theme_classic() + theme(plot.background = element_blank())
g3



# Correspondance of Ne between groups ----

stairwayplot[stairwayplot$year>99&stairwayplot$year<=100,]$Ne_median
stairwayplot[stairwayplot$year>99&stairwayplot$year<=100,]$Ne_97.5.
stairwayplot[stairwayplot$year>99&stairwayplot$year<=100,]$Ne_2.5.

gone[gone$Years<120,]

stairwayplot[stairwayplot$year>9900&stairwayplot$year<10200,]$Ne_median
stairwayplot[stairwayplot$year>9900&stairwayplot$year<10200,]$Ne_97.5.
stairwayplot[stairwayplot$year>9900&stairwayplot$year<10200,]$Ne_2.5.
df2[df2$xx<33000,]

df2[df2$xx<100000,]

list_min <- list()
for (i in 1:length(samples)){
  df <- df_list2[[i]][df_list2[[i]]$Generation<40,]
  df$Geometric_mean <- as.numeric(df$Geometric_mean)
  list_min[[i]] <- data.frame(NeMin=min(df$Geometric_mean),Generations=df[df$Geometric_mean==min(df$Geometric_mean),]$Generation)
  
}

df_min <- do.call(rbind, list_min)
max(df_min$NeMin)
min(df_min$NeMin)
mean(df_min$NeMin)

# Final Plot

g4 <- g2 + annotation_custom(ggplotGrob(g3),
                              xmin = 220, ymin = 1, ymax=6000,
                              xmax = 500)

ggarrange(ggarrange(g1,g, nrow=1, labels = c("A","B")),g4, nrow=2, labels=c("","C"))
ggsave("Manuscript/MainFigures/Figure3.pdf",width = 10, height = 10)
