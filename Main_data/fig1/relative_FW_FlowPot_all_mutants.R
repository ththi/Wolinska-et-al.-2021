
#GRAPH
library(Rmisc)

library(car)  
library(ggplot2)
library(sciplot)
library(AID) 
library(forecast) 
library(corrplot)


FW_graph <- read.table("fig1_fw_table.txt", header=T, sep="\t")
figures.dir<-""

summary(FW_graph)
FW_graph$biological <- factor(FW_graph$biological)

average<-aggregate(FW_graph$FW~FW_graph$genotype_microbes, FUN=mean)

write.table(average, file="average_FW_both_FlowPots.txt", col.names=TRUE, row.names=TRUE, append=FALSE)

###########
#getting single averages for each sterile genotype
summary(average$FW_graph$genotype_microbes)
colnames(average) <- c("genotype_microbes", "FW")

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


only_microbes<-subset(FW_graph, FW_graph$microbes=="BF")

test2 <-only_microbes
test2<-droplevels(test2)
test2<-na.omit(test2)




#calculating relative FW change for each mutant, normalized to WT
#WT
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='1aWT_BF') {
    test2$relative_FW[i] <- c(((test2$FW[i])-WT_sterile)/0.020766190)
  }
}

#WT 2
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='1aWT_BF_2') {
    test2$relative_FW[i] <- c(((test2$FW[i])-WT2_sterile)/0.012104365)
  }
}

#bb
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='bak1bkk1_BF') {
    test2$relative_FW[i] <- c(((test2$FW[i])-bb_sterile)/0.020766190)
  }
}

#bbc
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='bbc_BF_2') {
    test2$relative_FW[i] <- c(((test2$FW[i])-bbc_sterile)/0.012104365)
  }
}

#efc
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='efrfls2cerk1_BF') {
    test2$relative_FW[i] <- c(((test2$FW[i])-efc_sterile)/0.020766190)
  }
}

#lyk5
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='lyk5_BF_2') {
    test2$relative_FW[i] <- c(((test2$FW[i])-lyk5_sterile)/0.012104365)
  }
}

#cyp79
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='cyp79b2b3_BF') {
    test2$relative_FW[i] <- c(((test2$FW[i])-cyp79_sterile)/0.020766190)
  }
}

#cyp79-2
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='cyp79_BF_2') {
    test2$relative_FW[i] <- c(((test2$FW[i])-cyp792_sterile)/0.012104365)
  }
}


#wrky3340
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='wrky3340_BF') {
    test2$relative_FW[i] <- c(((test2$FW[i])-wrky3340_sterile)/0.020766190)
  }
}

#wrky33
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='wrky33_BF_2') {
    test2$relative_FW[i] <- c(((test2$FW[i])-wrky33_sterile)/0.012104365)
  }
}

#deps
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='deps_BF') {
    test2$relative_FW[i] <- c(((test2$FW[i])-deps_sterile)/0.020766190)
  }
}

#pad4
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='pad4_BF') {
    test2$relative_FW[i] <- c(((test2$FW[i])-pad4_sterile)/0.020766190)
  }
}


#35SBRI
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='35SBRI_BF_2') {
    test2$relative_FW[i] <- c(((test2$FW[i])-BRI_sterile)/0.012104365)
  }
}

#apex1
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='apex1_BF_2') {
    test2$relative_FW[i] <- c(((test2$FW[i])-apex1_sterile)/0.012104365)
  }
}

#apex2
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='apex2_BF_2') {
    test2$relative_FW[i] <- c(((test2$FW[i])-apex2_sterile)/0.012104365)
  }
}

#apex3
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='apex3_BF_2') {
    test2$relative_FW[i] <- c(((test2$FW[i])-apex3_sterile)/0.012104365)
  }
}

#rar1
for (i in 1:length(test2$FW)){ 
  if (test2$genotype_microbes[i]=='rar1_BF') {
    test2$relative_FW[i] <- c(((test2$FW[i])-rar1_sterile)/0.020766190)
  }
}



test2<-droplevels(test2)


test3<-na.omit(test2)





sum_test<-test3

library(plyr)
sum_test <- transform(sum_test,
                      genotype_microbes=revalue(genotype_microbes,c("1aWT_BF_2"="1aWT_BF")))

sum_test <- transform(sum_test,
                      genotype_microbes=revalue(genotype_microbes,c("cyp79_BF_2"="cyp79b2b3_BF")))

sum_test <- subset(sum_test, sum_test$genotype %in% c("1aWT", "35SBRI", "apex1", "apex2", "apex3", "bak1bkk1", "bbc",  
                                                      "cyp79", "cyp79b2b3", "deps", "efrfls2cerk1", 
                                                      "lyk5", "pad4",  "rar1", "wrky33", "wrky3340"))


#order<- c("1aWT_BF","bak1bkk1_BF","bbc_BF_2","efrfls2cerk1_BF","lyk5_BF_2","apex1_BF_2","apex2_BF_2","apex3_BF_2",
#          "wrky3340_BF","wrky33_BF_2",
#          "deps_BF","pad4_BF","cyp79b2b3_BF","35SBRI_BF_2","rar1_BF")

#sum_test$genotype_microbes <- factor(sum_test$genotype_microbes, levels=order)

#sum_test <- droplevels(sum_test)


##### change genotype names to latest version

sum_test[,"genotype_microbes"]=gsub("apex1_BF_2","hub1_BF_2",sum_test[,"genotype_microbes"])
sum_test[,"genotype_microbes"]=gsub("apex2_BF_2","apex_BF_2",sum_test[,"genotype_microbes"])
sum_test[,"genotype_microbes"]=gsub("apex3_BF_2","hub2_BF_2",sum_test[,"genotype_microbes"])

order<- c("1aWT_BF","bak1bkk1_BF","bbc_BF_2","efrfls2cerk1_BF","lyk5_BF_2","hub1_BF_2","apex_BF_2","hub2_BF_2",
         "wrky3340_BF","wrky33_BF_2",
          "deps_BF","pad4_BF","cyp79b2b3_BF","35SBRI_BF_2","rar1_BF")

sum_test$genotype_microbes <- factor(sum_test$genotype_microbes, levels=order)

#sum_test <- droplevels(sum_test)


b<-ggplot(sum_test, aes(x=genotype_microbes, y=relative_FW))+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
   geom_jitter(position=position_jitter(0.17), size=1.5, alpha=0.7)+
  #scale_colour_manual(values=as.character(colors$color)) +
#  coord_cartesian(ylim = c(-2, 5)) +
  scale_y_continuous(breaks=c(-1,0,1,5,10)) +
  ggtitle("Relative FW FlowPot both experiments")
b
ggsave(paste(figures.dir, "fig_1_C.pdf", sep=""), b)


xx=grep("a1aWT_no|a1aWT_BF",FW_graph[,"genotype_microbes_graph"])
new_fw=FW_graph[xx,]

c<-ggplot(new_fw, aes(x=genotype_microbes_graph, y=FW))+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
   geom_jitter(position=position_jitter(0.17), size=1.5, alpha=0.7)+
  #scale_colour_manual(values=as.character(colors$color)) +
#  coord_cartesian(ylim = c(-2, 5)) +
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15)) +
  ggtitle("FW")
  c
ggsave(paste(figures.dir, "fig_1_B.pdf", sep=""), c)


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
  capture.output(kruskal.test(data$FW ~ data$genotype_microbes_graph))
}

OT <- dlply(new_fw, .(microbes), model3)
cat("\n", file="Kruskal_wallis_fig1_b.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="Kruskal_wallis_fig1_b.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(FW ~ genotype_microbes_graph, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(new_fw, .(microbes), model4b)
write.table(ldply(KW_output), file="Kruskal_wallis_fig1_b.txt", col.names=TRUE, row.names=TRUE, append=TRUE)

########
########
model3 <- function(data){
  capture.output(kruskal.test(data$relative_FW ~ data$genotype_microbes))
}

OT <- dlply(sum_test, .(microbes), model3)
cat("\n", file="Kruskal_wallis_fig1_c.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="Kruskal_wallis_fig1_c.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(relative_FW ~ genotype_microbes, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(sum_test, .(microbes), model4b)
write.table(ldply(KW_output), file="Kruskal_wallis_fig1_c.txt", col.names=TRUE, row.names=TRUE, append=TRUE)







