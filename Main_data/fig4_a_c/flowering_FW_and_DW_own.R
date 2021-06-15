
#GRAPH
library(Rmisc)
library(car)
library(ggplot2)
library(sciplot)
library(AID)
library(forecast)
library(corrplot)


FW_DW <- read.table("flowering_raw_data_full.txt", header=T, sep="\t")
figures.dir <- ""


FW_DW <- subset(FW_DW, FW_DW$biological %in% c("bio2", "bio3") &
                 FW_DW$genotype %in% c("WT","cyp79b2b3"))

FW_rosette <-FW_DW[c(1:nrow(FW_DW)),c(1:9, 23)]
FW_stem <-FW_DW[c(1:nrow(FW_DW)),c(1:9, 24)]
DW_rosette <-FW_DW[c(1:nrow(FW_DW)),c(1:9, 25)]
DW_stem <-FW_DW[c(1:nrow(FW_DW)),c(1:9, 26)]



colors <- data.frame(group=c("cyp79b2b3", "WT"),
                    color=c("blue1", "black"))

FW_rosette$genotype <- factor(FW_rosette$genotype, levels=colors$group)
FW_stem$genotype <- factor(FW_stem$genotype, levels=colors$group)
DW_rosette$genotype <- factor(DW_rosette$genotype, levels=colors$group)
DW_stem$genotype <- factor(DW_stem$genotype, levels=colors$group)

shapes <- data.frame(group=c("bio1","bio2","bio3"),
                     shape=c(15,16,17))
FW_rosette$biological <- factor(FW_rosette$biological, levels=shapes$group)
FW_stem$biological <- factor(FW_stem$biological, levels=shapes$group)
DW_rosette$biological <- factor(DW_rosette$biological, levels=shapes$group)
DW_stem$biological <- factor(DW_stem$biological, levels=shapes$group)



l <- c("WT_sterile","WT_B","WT_F","WT_O","WT_BF","WT_BO","WT_FO","WT_BFO",
	          "cyp79b2b3_sterile","cyp79b2b3_B","cyp79b2b3_F","cyp79b2b3_O","cyp79b2b3_BF","cyp79b2b3_BO", "cyp79b2b3_FO","cyp79b2b3_BFO")




FW_rosette<-na.omit(FW_rosette)
FW_stem<-na.omit(FW_stem)
DW_rosette<-na.omit(DW_rosette)
DW_stem<-na.omit(DW_stem)

FW_rosette$genotype_treatment <- factor(FW_rosette$genotype_treatment, levels=l)
FW_stem$genotype_treatment <- factor(FW_stem$genotype_treatment, levels=l)
DW_rosette$genotype_treatment <- factor(DW_rosette$genotype_treatment, levels=l)
DW_stem$genotype_treatment <- factor(DW_stem$genotype_treatment, levels=l)



b<-ggplot(DW_rosette, aes(x=genotype_treatment, y=DW_rosette,  color=genotype))+theme_bw()+
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
  scale_colour_manual(values=as.character(colors$color)) +
  #scale_shape_manual(values=shapes$shape) +
  geom_jitter(position=position_jitter(0.17), size=1.5, alpha=0.7)+
  #coord_cartesian(ylim = c(0, 100)) +
  ggtitle("Rosette dry weight")+theme(axis.text.x=element_text(angle=90,hjust=1))
b

ggsave(paste(figures.dir, "fig4_A.pdf", sep=""), b)



library("rcompanion")
library("FSA")
library("plyr")
model3 <- function(data){
  capture.output(kruskal.test(DW_rosette$DW_rosette ~ DW_rosette$genotype_treatment))
}

OT <- dlply(DW_rosette, .(experiment), model3)
cat("\n", file="KW_fig4_a.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="KW_fig4_a.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(DW_rosette ~ genotype_treatment, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(DW_rosette, .(experiment), model4b)
write.table(ldply(KW_output), file="KW_fig4_a.txt", col.names=TRUE, row.names=TRUE, append=TRUE)


