
#GRAPH
library(Rmisc)

library(car) #### geht nicht mit version ### on dell 13 -> error 
library(ggplot2)
library(sciplot)
library(AID) #### fehler bei installation
library(forecast) #### fehler bei installation
library(corrplot)


FW_graph <- read.table("FW_raw_two_rounds_of_FlowPot_new_column_REM.txt", header=T, sep="\t")
figures.dir<-""

#summary(FW_graph)
FW_graph$biological <- factor(FW_graph$biological)

average<-aggregate(FW_graph$FW~FW_graph$genotype_microbes, FUN=mean)



WT_sterile<-subset(average, average$genotype_microbes=="1aWT_no")$FW
WT2_sterile<-subset(average, average$genotype_microbes=="1aWT_no_2")$FW
bb_sterile<-subset(average, average$genotype_microbes=="bak1bkk1_no")$FW
bbc_sterile<-subset(average, average$genotype_microbes=="bbc_no_2")$FW
efc_sterile<-subset(average, average$genotype_microbes=="efrfls2cerk1_no")$FW
lyk5_sterile<-subset(average, average$genotype_microbes=="lyk5_no_2")$FW
cyp79_sterile<-subset(average, average$genotype_microbes=="cyp79b2b3_no")$FW
cyp792_sterile<-subset(average, average$genotype_microbes=="cyp79_no_2")$FW
wrky3340_sterile<-subset(average, average$genotype_microbes=="wrky3340_no")$FW
wrky33_sterile<-subset(average, average$genotype_microbes=="wrky33_no_2")$FW
deps_sterile<-subset(average, average$genotype_microbes=="deps_no")$FW
pad4_sterile<-subset(average, average$genotype_microbes=="pad4_no")$FW
BRI_sterile<-subset(average, average$genotype_microbes=="35SBRI_no_2")$FW
apex1_sterile<-subset(average, average$genotype_microbes=="apex1_no_2")$FW
apex2_sterile<-subset(average, average$genotype_microbes=="apex2_no_2")$FW
apex3_sterile<-subset(average, average$genotype_microbes=="apex3_no_2")$FW
rar1_sterile<-subset(average, average$genotype_microbes=="rar1_no")$FW


only_no_microbes<-subset(FW_graph, FW_graph$microbes=="no")

test2 <-only_no_microbes
test2<-droplevels(test2)
test2<-na.omit(test2)



test2<-droplevels(test2)


test3<-na.omit(test2)





sum_test<-test3

library(plyr)
sum_test <- transform(sum_test,
                      genotype_microbes=revalue(genotype_microbes,c("1aWT_no_2"="1aWT_no")))

sum_test <- transform(sum_test,
                      genotype_microbes=revalue(genotype_microbes,c("cyp79_no_2"="cyp79b2b3_no")))

sum_test <- subset(sum_test, sum_test$genotype %in% c("1aWT", "35SBRI", "apex1", "apex2", "apex3", "bak1bkk1", "bbc",  
                                                      "cyp79", "cyp79b2b3", "deps", "efrfls2cerk1", 
                                                      "lyk5", "pad4",  "rar1", "wrky33", "wrky3340"))


#order<- c("1aWT_no","bak1bkk1_no","bbc_no_2","efrfls2cerk1_no","lyk5_no_2","apex1_no_2","apex2_no_2","apex3_no_2",
#          "wrky3340_no","wrky33_no_2",
#          "deps_no","pad4_no","cyp79b2b3_no","35SBRI_no_2","rar1_no")

#sum_test$genotype_microbes <- factor(sum_test$genotype_microbes, levels=order)

sum_test[,"genotype_microbes"]=gsub("apex1_no_2","hub1_no_2",sum_test[,"genotype_microbes"])
sum_test[,"genotype_microbes"]=gsub("apex2_no_2","apex_no_2",sum_test[,"genotype_microbes"])
sum_test[,"genotype_microbes"]=gsub("apex3_no_2","hub2_no_2",sum_test[,"genotype_microbes"])

order<- c("1aWT_no","bak1bkk1_no","bbc_no_2","efrfls2cerk1_no","lyk5_no_2","hub1_no_2","apex_no_2","hub2_no_2",
          "wrky3340_no","wrky33_no_2",
          "deps_no","pad4_no","cyp79b2b3_no","35SBRI_no_2","rar1_no")

sum_test$genotype_microbes <- factor(sum_test$genotype_microbes, levels=order)


sum_test <- droplevels(sum_test)
b<-ggplot(sum_test, aes(x=genotype_microbes, y=FW))+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
   geom_jitter(position=position_jitter(0.17), size=1.5, alpha=0.7)+
  #scale_colour_manual(values=as.character(colors$color)) +
#  coord_cartesian(ylim = c(-2, 5)) +
  scale_y_continuous(breaks=c(-1,0,1,5,10)) +
  ggtitle("germ free")
b
ggsave(paste(figures.dir, "fig_S2_B.pdf", sep=""), b)





#STATISTICS


#checking normality


#qqnorm(sum_test$relative_FW)
#qqline(sum_test$relative_FW, col="red")
#shapiro.test(sum_test$relative_FW)

#hist(FW_graph$FW)

#plot(FW_graph$FW)


#Kruskall-Wallis
library("rcompanion")
library("FSA")
library("plyr")
#relabund.phy.sums.WT.pbs3.leaves <- subset(relabund.phy.sums, relabund.phy.sums$WT_pbs3_leaves=="yes")

model3 <- function(data){
  capture.output(kruskal.test(data$FW ~ data$genotype_microbes))
}

OT <- dlply(sum_test, .(microbes), model3)
cat("\n", file="Kruskal_wallis_figS2_b.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="Kruskal_wallis_figS2_b.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(FW ~ genotype_microbes, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(sum_test, .(microbes), model4b)
write.table(ldply(KW_output), file="Kruskal_wallis_figS2_b.txt", col.names=TRUE, row.names=TRUE, append=TRUE)

