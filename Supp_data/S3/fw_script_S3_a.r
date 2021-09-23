

#GRAPH FW WEEK 5 AND 8 TOGETHER
library(Rmisc)
library(car)
library(ggplot2)
library(sciplot)
library(AID)
library(forecast)
library(corrplot)


FW_graph <- read.table("FW_both_experiments_S3.txt", header=T, sep="\t")


summary(FW_graph)
FW_graph$biological <- factor(FW_graph$biological)
FW_graph$time_point <- factor(FW_graph$time_point)


#FW_graph_T1<-subset(FW_graph,FW_graph$time_point=="1")
#FW_graph_T1<-FW_graph_T1[c(1:nrow(FW_graph_T1)),c(1:6)]
#FW_graph_T1$genotype=factor(FW_graph_T1$genotype,levels(FW_graph_T1$Genotype)[c(10,2,5,)])

#colors time point
colors <- data.frame(group=c("tp1","tp2"),
                     color=c("grey","grey"))


#levels genotype
l <- c("WT","WT2","bak1bkk1","bak1bkk1cerk1","efrfls2cerk1","lyk5","cyp79b2b3","quin","wrky3340","wrky33","wrky40","deps","pad4","bri301",
       "35SBRI","apex1","apex2","apex3","rar1","kai2")

#FW_graph$genotype <- factor(FW_graph$genotype, levels=colors$group)
FW_graph$time_point <- factor(FW_graph$time_point, levels=colors$group)

FW_graph$genotype <- factor(FW_graph$genotype, levels=l)
#FW_graph$biological <- factor(FW_graph$biological, levels=colors$group)


FW_graph<-na.omit(FW_graph)




#making relative FW in comparison to WT


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
#relative FW but different = data point / average WT per bio replicate

test2 <-FW_graph

for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio1'& test2$time_point[i]=="tp1" & test2$experiment[i]=="B1") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.01206667)
  }
}

#Bio 1 time point 2 old
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio1'& test2$time_point[i]=="tp2" & test2$experiment[i]=="B1") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.09581250)
  }
}

#Bio 2 time point 1 old
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio2'& test2$time_point[i]=="tp1" & test2$experiment[i]=="B1") {
    test2$relative_FW[i] <- c(log2((test2$FW[i])/0.04335000))
  }
}

#Bio 2 time point 2 old
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio2'& test2$time_point[i]=="tp2" & test2$experiment[i]=="B1") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.09673684)
  }
}

#Bio 3 time point 1 old
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio3'& test2$time_point[i]=="tp1"  & test2$experiment[i]=="B1") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.01445000)
  }
}

#Bio 3 time point 2 old
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio3'& test2$time_point[i]=="tp2" & test2$experiment[i]=="B1") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.06388889)
  }
}




#Bio 1 time point 1 new
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio1'& test2$time_point[i]=="tp1" & test2$experiment[i]=="B2") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.069105)
  }
}

#Bio 1 time point 2 new
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio1'& test2$time_point[i]=="tp2" & test2$experiment[i]=="B2") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.199115)
  }
}

#Bio 2 time point 1 new
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio2'& test2$time_point[i]=="tp1" & test2$experiment[i]=="B2") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.02649)
  }
}

#Bio 2 time point 2 new
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio2'& test2$time_point[i]=="tp2" & test2$experiment[i]=="B2") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.04615)
  }
}

#Bio 3 time point 1 new
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio3'& test2$time_point[i]=="tp1" & test2$experiment[i]=="B2") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.0217)
  }
}

#Bio 3 time point 2 new
for (i in 1:length(test2$FW)){ 
  if (test2$biological[i]=='bio3'& test2$time_point[i]=="tp2" & test2$experiment[i]=="B2") {
    test2$relative_FW[i] <- c((test2$FW[i])/0.0513)
  }
}


#checking if the if function worked (1st trial with substration was fine)
test2$check<-c(-(test2$FW) + test2$relative_FW)
unique(test2$check)

relFW_week5<-subset(test2,test2$time_point=="tp1")

#combining WT names
relFW_week5[,"genotype"]=gsub("WT2","WT",relFW_week5[,"genotype"])

relFW_week5_thesis <- subset(relFW_week5, relFW_week5$genotype %in% c("WT","bak1bkk1","bak1bkk1cerk1","efrfls2cerk1","lyk5","apex1","apex2","apex3","wrky3340","wrky33","wrky40",
                                                                      "deps","pad4","cyp79b2b3","35SBRI","bri301","rar1"))
#new order
l <- c("WT","bak1bkk1","bak1bkk1cerk1","efrfls2cerk1","lyk5","apex1","apex2","apex3","wrky3340","wrky33","wrky40",
       "deps","pad4","cyp79b2b3","35SBRI","bri301","rar1")

relFW_week5_thesis$genotype <- factor(relFW_week5_thesis$genotype, levels=l)

####

relFW_week5_thesis[,"genotype"]=gsub("apex1","hub1",relFW_week5_thesis[,"genotype"])
relFW_week5_thesis[,"genotype"]=gsub("apex2","apex",relFW_week5_thesis[,"genotype"])
relFW_week5_thesis[,"genotype"]=gsub("apex3","hub2",relFW_week5_thesis[,"genotype"])

l <- c("WT","bak1bkk1","bak1bkk1cerk1","efrfls2cerk1","lyk5","hub1","apex","hub2","wrky3340","wrky33","wrky40",
       "deps","pad4","cyp79b2b3","35SBRI","bri301","rar1")

relFW_week5_thesis$genotype <- factor(relFW_week5_thesis$genotype, levels=l)

#relative FW week 5 
b<-ggplot(relFW_week5_thesis, aes(x=genotype, y=relative_FW, color=time_point))+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent", outlier.colour="black", outlier.shape=16)+
  geom_jitter(position=position_jitter(0.17), size=1.5, alpha=0.7)+
  scale_colour_manual(values=as.character(colors$color)) +
  ggtitle("Relative FW ")+ ylim(-2, 3.5)
b

ggsave( "relativeFW_fig_s3_a.pdf", b)

#stat testing

library("rcompanion")
library("FSA")
library("plyr")


########
########

kw_result=kruskal.test(relFW_week5$relative_FW ~ relFW_week5$genotype)

OT1=dunnTest(relFW_week5[,"relative_FW"]~relFW_week5[,"genotype"],method="bh")
OT1 <- OT1$res
OT1a <- cldList(comparison = OT1$Comparison, p.value = OT1$P.adj, threshold = 0.05)








