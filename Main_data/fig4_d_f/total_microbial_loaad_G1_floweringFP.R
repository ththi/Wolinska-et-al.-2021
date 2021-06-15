
figures.dir<-""



bacteria <- read.table( "total_load_bacteria_G1_rem.txt", header=T, sep="\t")
fungi <- read.table("total_load_fungi_G1_rem.txt", header=T, sep="\t")
oomycetes <- read.table("total_load_oomycetes_G1_rem.txt", header=T, sep="\t")

bacteria_noNA <- bacteria[!is.na(bacteria$X16S_UBQ_ref), ]
fungi_noNA <- fungi[!is.na(fungi$ITS1_UBQ_ref), ]
oomycetes_noNA <- oomycetes[!is.na(oomycetes$oITS1_UBQ_ref), ]




colors <- data.frame(group=c("cyp79b2b3", "WT"),
                    color=c("blue1","black"))


colors <- colors[colors$group %in% bacteria_noNA$genotype, ]
bacteria_noNA$genotype<- factor(bacteria_noNA$genotype, levels=colors$group)
colors <- colors[colors$group %in% fungi_noNA$genotype, ]
fungi_noNA$genotype<- factor(fungi_noNA$genotype, levels=colors$group)
colors <- colors[colors$group %in% oomycetes_noNA$genotype, ]
oomycetes_noNA$genotype<- factor(oomycetes_noNA$genotype, levels=colors$group)

#order

l <- c("WT_B","WT_F","WT_O","WT_BF","WT_BO","WT_FO","WT_BFO", "cyp79b2b3_B","cyp79b2b3_F","cyp79b2b3_O","cyp79b2b3_BF","cyp79b2b3_BO","cyp79b2b3_FO","cyp79b2b3_BFO")

bacteria_noNA$genotype_treatment <- factor(bacteria_noNA$genotype_treatment, levels=l)
fungi_noNA$genotype_treatment <- factor(fungi_noNA$genotype_treatment, levels=l)
oomycetes_noNA$genotype_treatment <- factor(oomycetes_noNA$genotype_treatment, levels=l)

bacteria_noNA_bio23 <- subset(bacteria_noNA, bacteria_noNA$biological %in% c("bio2", "bio3"))
fungi_noNA_bio23 <- subset(fungi_noNA, fungi_noNA$biological %in% c("bio2", "bio3"))
oomycetes_noNA_bio23 <- subset(oomycetes_noNA, oomycetes_noNA$biological %in% c("bio2", "bio3"))


#graph
library(ggplot2)
#graph without outliers being marked
p_bacteria <- ggplot(bacteria_noNA_bio23, aes(x=genotype_treatment, y=X16S_UBQ_ref, color=genotype)) +
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent")+
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x="", y="bacteria/plant/ref ratio") +theme_bw()+ggtitle("bacteria abundance in coomparison to plant gDNA floweringFP G1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

p_bacteria
ggsave(paste(figures.dir, "fig4_D.pdf", sep=""), p_bacteria)

p_fungi <- ggplot(fungi_noNA_bio23, aes(x=genotype_treatment, y=ITS1_UBQ_ref, color=genotype)) +
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent")+
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x="", y="fungi/plant/ref ratio") +theme_bw()+ggtitle("fungi abundance in coomparison to plant gDNA floweringFP G1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

p_fungi
ggsave(paste(figures.dir, "fig4_E.pdf", sep=""), p_fungi)

p_oomycetes <- ggplot(oomycetes_noNA_bio23, aes(x=genotype_treatment, y=oITS1_UBQ_ref, color=genotype)) +
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent")+
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x="", y="oomycetes/plant/ref ratio") +theme_bw()+ggtitle("oomycetes abundance in coomparison to plant gDNA floweringFP G1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

p_oomycetes
ggsave(paste(figures.dir, "fig4_F.pdf", sep=""), p_oomycetes)






#############################################################

library("rcompanion")
library("FSA")
library("plyr")
#################################################################
#bacteria
model3 <- function(data){
  capture.output(kruskal.test(data$X16S_UBQ_ref ~ data$genotype_treatment))
}

OT <- dlply(bacteria_noNA_bio23, .(project), model3)
cat("\n", file="Kruskal_wallis_fig4_d.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="Kruskal_wallis_fig4_d.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(X16S_UBQ_ref ~ genotype_treatment, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(bacteria_noNA_bio23, .(project), model4b)
write.table(ldply(KW_output), file="Kruskal_wallis_fig4_d.txt", col.names=TRUE, row.names=TRUE, append=TRUE)
#################################################################
#fungi
model3 <- function(data){
  capture.output(kruskal.test(data$ITS1_UBQ_ref ~ data$genotype_treatment))
}

OT <- dlply(fungi_noNA_bio23, .(project), model3)
cat("\n", file="Kruskal_wallis_fig4_e.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="Kruskal_wallis_fig4_e.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(ITS1_UBQ_ref ~ genotype_treatment, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(fungi_noNA_bio23, .(project), model4b)
write.table(ldply(KW_output), file="Kruskal_wallis_fig4_e.txt", col.names=TRUE, row.names=TRUE, append=TRUE)
#################################################################
#oomycetes
model3 <- function(data){
  capture.output(kruskal.test(data$oITS1_UBQ_ref ~ data$genotype_treatment))
}

OT <- dlply(oomycetes_noNA_bio23, .(project), model3)
cat("\n", file="Kruskal_wallis_fig4_f.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="Kruskal_wallis_fig4_f.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(oITS1_UBQ_ref ~ genotype_treatment, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(oomycetes_noNA_bio23, .(project), model4b)
write.table(ldply(KW_output), file="Kruskal_wallis_fig4_f.txt", col.names=TRUE, row.names=TRUE, append=TRUE)






