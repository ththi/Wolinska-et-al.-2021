
#GRAPH
library(Rmisc)
library(car)
library(ggplot2)
library(sciplot)
library(AID)
library(forecast)
library(corrplot)

FW.dir<-""
siliques_stem <- read.table("heat_killed_FW_raw.txt", header=T, sep="\t")



#levels genotype_treatment_plant by genotype
l <- c("sterile","HK","BFO")


siliques_stem$treatment<- factor(siliques_stem$treatment, levels=l)



siliques_stem<-na.omit(siliques_stem)



b<-ggplot(siliques_stem, aes(x=treatment, y=FW_g))+theme_bw()+
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
  #scale_colour_manual(values=as.character(colors$color)) +
  #scale_shape_manual(values=shapes$shape) +
  geom_jitter(position=position_jitter(0.17), size=1.5, alpha=0.7)+
  #coord_cartesian(ylim = c(0, 0.06)) +
  ggtitle("heat kill WT " )
b

ggsave(paste(FW.dir, "fig_S2_A.pdf", sep=""), b)





#STATISTICS


#checking normality

#qqnorm(siliques_stem$FW_g)
#qqline(siliques_stem$FW_g, col="red")
#shapiro.test(siliques_stem$FW_g)



#hist(siliques_stem$FW_g)

#plot(siliques_stem$FW_g)


#Kruskall-Wallis
library("rcompanion")
library("FSA")
library("plyr")
#relabund.phy.sums.WT.pbs3.leaves <- subset(relabund.phy.sums, relabund.phy.sums$WT_pbs3_leaves=="yes")

model3 <- function(data){
  capture.output(kruskal.test(siliques_stem$FW_g ~ siliques_stem$treatment))
}

OT <- dlply(siliques_stem, .(experiment), model3)
cat("\n", file="KW_S2_a.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="KW_S2_a.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(FW_g ~ treatment, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(siliques_stem, .(experiment), model4b)
write.table(ldply(KW_output), file="KW_S2_a.txt", col.names=TRUE, row.names=TRUE, append=TRUE)

#summary(FW_graph$genotype_microbes)





##########################################################################################################

