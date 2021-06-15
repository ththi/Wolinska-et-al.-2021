
#GRAPH
library(Rmisc)
library(car)
library(ggplot2)
library(sciplot)
library(AID)
library(forecast)
library(corrplot)


FW_graph <- read.table("FW_project_A5_all_clean_order_for_baplot_rem.txt", header=T, sep="\t")




average<-aggregate(FW_graph$FW_g~FW_graph$genotype_treatment, FUN=mean)


colnames(average) <- c("genotype_microbes", "FW")

WT_sterile<-subset(average, average$genotype_microbes=="WT_sterile")$FW
WT_BFO<-subset(average, average$genotype_microbes=="WT_BFO")$FW

cyp79b2b3_sterile<-subset(average, average$genotype_microbes=="cyp79b2_b3_sterile")$FW
cyp71a12a13_sterile<-subset(average, average$genotype_microbes=="cyp71a12_a13_sterile")$FW
pen2_pad3_sterile<-subset(average, average$genotype_microbes=="pen2_pad3_sterile")$FW
pen2_sterile<-subset(average, average$genotype_microbes=="pen2_sterile")$FW
pad3_sterile<-subset(average, average$genotype_microbes=="pad3_sterile")$FW
mybtriple_sterile<-subset(average, average$genotype_microbes=="myb_triple_sterile")$FW
pykbglu_sterile<-subset(average, average$genotype_microbes=="pyk_bglu_sterile")$FW
cyp71a27_sterile<-subset(average, average$genotype_microbes=="cyp71a27_sterile")$FW

#get WT_BFO average - WT_sterile average
WT_difference <- WT_BFO - WT_sterile

#subset for only microbes
only_microbes<-subset(FW_graph, FW_graph$treatment=="BFO")

test2 <-only_microbes
test2<-droplevels(test2)
test2<-na.omit(test2)



#calculating relative FW change for each mutant, normalized to WT
#WT
for (i in 1:length(test2$FW_g)){ 
  if (test2$genotype_treatment[i]=='WT_BFO') {
    test2$relative_FW[i] <- c(((test2$FW_g[i])-WT_sterile)/WT_difference)
  }
}

#cyp79b2b3
for (i in 1:length(test2$FW_g)){ 
  if (test2$genotype_treatment[i]=='cyp79b2_b3_BFO') {
    test2$relative_FW[i] <- c(((test2$FW_g[i])-cyp79b2b3_sterile)/WT_difference)
  }
}

#cyp71a12/a13
for (i in 1:length(test2$FW_g)){ 
  if (test2$genotype_treatment[i]=='cyp71a12_a13_BFO') {
    test2$relative_FW[i] <- c(((test2$FW_g[i])-cyp71a12a13_sterile)/WT_difference)
  }
}

#pen2/pad3
for (i in 1:length(test2$FW_g)){ 
  if (test2$genotype_treatment[i]=='pen2_pad3_BFO') {
    test2$relative_FW[i] <- c(((test2$FW_g[i])-pen2_pad3_sterile)/WT_difference)
  }
}

#pen2
for (i in 1:length(test2$FW_g)){ 
  if (test2$genotype_treatment[i]=='pen2_BFO') {
    test2$relative_FW[i] <- c(((test2$FW_g[i])-pen2_sterile)/WT_difference)
  }
}

#pad3
for (i in 1:length(test2$FW_g)){ 
  if (test2$genotype_treatment[i]=='pad3_BFO') {
    test2$relative_FW[i] <- c(((test2$FW_g[i])-pad3_sterile)/WT_difference)
  }
}

#myb_triple
for (i in 1:length(test2$FW_g)){ 
  if (test2$genotype_treatment[i]=='myb_triple_BFO') {
    test2$relative_FW[i] <- c(((test2$FW_g[i])-mybtriple_sterile)/WT_difference)
  }
}

#pyk/bglu
for (i in 1:length(test2$FW_g)){ 
  if (test2$genotype_treatment[i]=='pyk_bglu_BFO') {
    test2$relative_FW[i] <- c(((test2$FW_g[i])-pykbglu_sterile)/WT_difference)
  }
}

#cyp71a27
for (i in 1:length(test2$FW_g)){ 
  if (test2$genotype_treatment[i]=='cyp71a27_BFO') {
    test2$relative_FW[i] <- c(((test2$FW_g[i])-cyp71a27_sterile)/WT_difference)
  }
}




test2<-droplevels(test2)


test3<-na.omit(test2)







sum_test<-test3



sum_test <- droplevels(sum_test)

no_pencyp_quin<-subset(sum_test, sum_test$genotype %in% c("cyp71a12_a13", "cyp71a27", "cyp79b2_b3", "myb_triple", "pad3",
                                                     "pen2", "pen2_pad3", "pyk_bglu", "WT"))

no_pencyp_quin <-subset(no_pencyp_quin, no_pencyp_quin$relative_FW <= 9)

b<-ggplot(no_pencyp_quin, aes(x=genotype_barplot, y=relative_FW))+theme_bw()+
  geom_boxplot(alpha=1, outlier.size=0.5, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
  geom_jitter(position=position_jitter(0.17), size=0.5, alpha=0.7)+
  scale_y_continuous(breaks=c(-1,0,1,3,7,10)) +
  #scale_colour_manual(values=as.character(colors$color)) +
  #coord_cartesian(ylim = c(-2, 9)) +
  ggtitle("Relative growth promotion")+theme(axis.text.x=element_text(angle=90,hjust=1))


ggsave( "fig_S7_c.pdf", b)

use_tab<-subset(FW_graph, FW_graph$treatment %in% c("sterile"))

b<-ggplot(use_tab, aes(x=genotype_barplot, y=FW_g))+theme_bw()+
  geom_boxplot(alpha=1, outlier.size=0.5, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
  geom_jitter(position=position_jitter(0.17), size=0.5, alpha=0.7)+
  scale_y_continuous(breaks=c(0,0.02,0.06,0.07)) +
  #scale_colour_manual(values=as.character(colors$color)) +
  #coord_cartesian(ylim = c(-2, 9)) +
  ggtitle("Shoot fresh weight")+theme(axis.text.x=element_text(angle=90,hjust=1))


ggsave( "fig_S7_b.pdf", b)


#STATISTICS

above_9 <-subset(no_pencyp_quin, no_pencyp_quin$relative_FW > 9)
#write.table(sum_test, file="outliers_above_9_no_pencyp_quin_file.txt", col.names=TRUE, row.names=TRUE, append=FALSE)

#checking normality


#qqnorm(no_pencyp_quin$relative_FW)
#qqline(no_pencyp_quin$relative_FW, col="red")
#shapiro.test(no_pencyp_quin$relative_FW)

#hist(FW_graph$FW)

#plot(FW_graph$FW)


#Kruskall-Wallis
library("rcompanion")
library("FSA")
library("plyr")
#relabund.phy.sums.WT.pbs3.leaves <- subset(relabund.phy.sums, relabund.phy.sums$WT_pbs3_leaves=="yes")

model3 <- function(data){
  capture.output(kruskal.test(data$relative_FW ~ data$genotype_barplot))
}

OT <- dlply(no_pencyp_quin, .(experiment), model3)
cat("\n", file="kw_s7_c.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="kw_s7_c.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(relative_FW ~ genotype_barplot, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(no_pencyp_quin, .(experiment), model4b)
write.table(ldply(KW_output), file="kw_s7_c.txt", col.names=TRUE, row.names=TRUE, append=TRUE)






