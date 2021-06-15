
#GRAPH
library(Rmisc)
library(car)
library(ggplot2)
library(sciplot)
library(AID)
library(forecast)
library(corrplot)


FW_graph <- read.table("FW_D1_and_A6_cyp_WT_with_zero_rem.txt", header=T, sep="\t")



FW_graph <- subset(FW_graph, FW_graph$genotype %in% c("cyp79","WT"))

colors <- data.frame(group=c("cyp79","WT"),
                     color=c("blue1","black"))

FW_graph$genotype <- factor(FW_graph$genotype, levels=colors$group)

#levels BFO/no
l <- c("WT_BFO","WT_B","WT_F","WT_O","WT_BF","WT_BO","WT_FO",
       "cyp79_BFO","cyp79_B","cyp79_F","cyp79_O","cyp79_BF","cyp79_BO","cyp79_FO")

FW_graph$genotype_treatment<- factor(FW_graph$genotype_treatment, levels=l)

FW_graph <- droplevels(FW_graph)


#order for genotype_microbes

FW_graph<-na.omit(FW_graph)
summary(FW_graph$genotype_treatment)


b<-ggplot(FW_graph, aes(x=genotype_treatment, y=FW,  color=genotype))+theme_bw()+
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
  scale_colour_manual(values=as.character(colors$color)) +
  #scale_shape_manual(values=shapes$shape) +
  geom_jitter(position=position_jitter(0.17), size=1.5, alpha=0.7)+theme(axis.text.x=element_text(angle=90,hjust=1))
  ggtitle("FW FlowPot fractionation project D1")
b

ggsave("fig_s9_a.pdf", b)

library("rcompanion")
library("FSA")
library("plyr")
#relabund.phy.sums.WT.pbs3.leaves <- subset(relabund.phy.sums, relabund.phy.sums$WT_pbs3_leaves=="yes")

model3 <- function(data){
  capture.output(kruskal.test(data$FW ~ data$genotype_treatment))
}

OT <- dlply(FW_graph, .(project), model3)
cat("\n", file="kw_s9_a.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="kw_s9_a.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(FW ~ genotype_treatment, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(FW_graph, .(project), model4b)
write.table(ldply(KW_output), file="kw_s9_a.txt", col.names=TRUE, row.names=TRUE, append=TRUE)

