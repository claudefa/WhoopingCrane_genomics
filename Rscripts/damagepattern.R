#Plot Mapdamage plots all together for WC

library(ggplot2)

setwd("~/Documents/OneDrive - University of Copenhagen/Whooping_Crane/")

# read 5 prime
Fivep <- list.files("Plots/MapDamage/Frequencies/", pattern="_5pCtoT_freq.txt",full.names=TRUE) # boostrap taking one sample out
samples <- gsub('_5pCtoT_freq.txt', '', basename(Fivep))
type <- rep(c("Modern","Historical"), times=c(37,15))

Fivep_list <- lapply(Fivep,
                         FUN = function(files) {
                           read.csv2(files, sep="\t")
                         })
df_list2 <- list()

for (i in 1:length(samples)){
  df_list2[[i]] <- cbind.data.frame(Sample=samples[i], Fivep_list[[i]], Type=type[i])
}

prime5_df <- do.call(rbind,df_list2)
prime5_df$X5pC.T<-as.numeric(prime5_df$X5pC.T)

g <- ggplot(prime5_df, aes(pos,X5pC.T, color=Type, group=interaction(Type,Sample)))
g <- g+ geom_line() + ylim(c(0,0.1)) + theme_classic() + ylab("Frequency") + xlab("Position")+
  ggtitle("5' - C>T") +  scale_color_manual(values= c("#6d466bff","#EB9486" )) +
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))

g

ggsave("Manuscript/Plots_supplementary/FigureS2_damange.pdf", height = 4,width = 6)

g <- ggplot(prime5_df, aes(pos,X5pC.T, color=Type, group=interaction(Type,Sample)))
g <- g+ geom_line() + ylim(c(0,0.1)) + theme_classic() + ylab("Frequency") + xlab("Position")+
  facet_wrap(.~Sample)+
  ggtitle("5' - C>T") +  scale_color_manual(values= c("#6d466bff","#EB9486" )) +
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8))

g

# read 3 prime
Threep <- list.files("Plots/MapDamage/Frequencies/", pattern="_3pGtoA_freq.txt",full.names=TRUE) # boostrap taking one sample out
samples <- gsub('_3pGtoA_freq.txt', '', basename(Threep))
type <- rep(c("Modern","Historical"), times=c(37,15))

Threep_list <- lapply(Threep,
                     FUN = function(files) {
                       read.csv2(files, sep="\t")
                     })
df_list2 <- list()

for (i in 1:length(samples)){
  df_list2[[i]] <- cbind.data.frame(Sample=samples[i], Threep_list[[i]], Type=type[i])
}

prime3_df <- do.call(rbind,df_list2)
prime3_df$X3pG.A<-as.numeric(prime3_df$X3pG.A)

g <- ggplot(prime3_df, aes(pos,X3pG.A, color=Type, group=interaction(Type,Sample)))
g <- g+ geom_line() + ylim(c(0,0.1)) + theme_classic() + ylab("Frequency") + xlab("Position")+
  ggtitle("3' - G>A")
g
