
#GRAPH
library(Rmisc)
library(car)
library(ggplot2)
library(sciplot)
library(AID)
library(forecast)
library(corrplot)


siliques_stem <- read.table("flowering_raw_data_full.txt", header=T, sep="\t")
figures.dir <- ""
#subset for FW and DW

siliques_stem <- subset(siliques_stem, siliques_stem$biological %in% c("bio2", "bio3") &
                          siliques_stem$genotype %in% c("WT", "cyp79b2b3"))

siliques <-siliques_stem[c(1:nrow(siliques_stem)),c(1:9, 27)]
stem_cm <-siliques_stem[c(1:nrow(siliques_stem)),c(1:9, 28)]
no_stems <-siliques_stem[c(1:nrow(siliques_stem)),c(1:9, 29)]

colors <- data.frame(group=c("cyp79b2b3", "WT"),
                     color=c("blue1", "black"))

siliques$genotype <- factor(siliques$genotype, levels=colors$group)
stem_cm$genotype <- factor(stem_cm$genotype, levels=colors$group)
no_stems$genotype <- factor(no_stems$genotype, levels=colors$group)


shapes <- data.frame(group=c("bio1","bio2","bio3"),
                     shape=c(15,16,17))
siliques$biological <- factor(siliques$biological, levels=shapes$group)
stem_cm$biological <- factor(stem_cm$biological, levels=shapes$group)
no_stems$biological <- factor(no_stems$biological, levels=shapes$group)


#levels genotype_treatment_plant by genotype
l <- c("WT_sterile","WT_B","WT_F","WT_O","WT_BF","WT_BO","WT_FO","WT_BFO",
       "cyp79b2b3_sterile","cyp79b2b3_B","cyp79b2b3_F","cyp79b2b3_O","cyp79b2b3_BF","cyp79b2b3_BO","cyp79b2b3_FO","cyp79b2b3_BFO")


siliques$genotype_treatment<- factor(siliques$genotype_treatment, levels=l)
stem_cm$genotype_treatment<- factor(stem_cm$genotype_treatment, levels=l)
no_stems$genotype_treatment<- factor(no_stems$genotype_treatment, levels=l)



siliques<-na.omit(siliques)
stem_cm<-na.omit(stem_cm)
no_stems<-na.omit(no_stems)


b<-ggplot(siliques, aes(x=genotype_treatment, y=siliques_number,  color=genotype))+theme_bw()+
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
  scale_colour_manual(values=as.character(colors$color)) +
  #scale_shape_manual(values=shapes$shape) +
  geom_jitter(position=position_jitter(0.17), size=1.5, alpha=0.7)+
  #coord_cartesian(ylim = c(0, 100)) +
  ggtitle("Number of siliques after 9 weeks")+theme(axis.text.x=element_text(angle=90,hjust=1))
b

ggsave(paste(figures.dir, "fig4_C.pdf", sep=""), b)





library("rcompanion")
library("FSA")
library("plyr")
model3 <- function(data){
  capture.output(kruskal.test(siliques$siliques_number ~ siliques$genotype_treatment))
}

OT <- dlply(siliques, .(experiment), model3)
cat("\n", file="KW_fig4_c.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="KW_fig4_c.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(siliques$siliques_number ~ genotype_treatment, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(siliques, .(experiment), model4b)
write.table(ldply(KW_output), file="KW_fig4_c.txt", col.names=TRUE, row.names=TRUE, append=TRUE)

