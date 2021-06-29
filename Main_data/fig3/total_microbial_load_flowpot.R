
figure.dir<-""

#########################################################################################



#qPCR_test <- read.table("bacteria_total_abundance_FlowPot_nosoil_2aWT_REM.txt", header=T, sep="\t")
#qPCR_test <- read.table("total_fungal_abundance_FlowPot_corrected_nosoil_2aWT_REM.txt", header=T, sep="\t")
qPCR_test <- read.table("oomycetes_total_abundance_FlowPot_no_soil_2aWT_REM.txt", header=T, sep="\t")

###

qPCR_BFO = subset(qPCR_test, 
                  qPCR_test$microbes %in% c("BFO")) 
qPCR_BFO<-droplevels(qPCR_BFO)


#for correlation plot to get a mean/median

#for fungi
#qPCR_BFO[,"genotype"]=gsub("WT_new","WT",qPCR_BFO[,"genotype"])
#qPCR_BFO[,"genotype"]=gsub("cyp79b2b3_new","cyp79b2b3",qPCR_BFO[,"genotype"])

#qPCR_BFO[,"genotype2"]=gsub("WT_new","WT",qPCR_BFO[,"genotype2"])
#qPCR_BFO[,"genotype2"]=gsub("cyp79b2b3_new","cyp79b2b3",qPCR_BFO[,"genotype2"])

#average<-aggregate(qPCR_BFO$ITS_UBQ_ref~qPCR_BFO$genotype, FUN=mean)

#write.table(average, file="U:/PhD/cyp79_manuscript/FIGURES/additional_files/Figure_R3/files_for_nwe_scatterplot/fungi_load_mean.txt", col.names=TRUE, row.names=TRUE, append=FALSE)


colors <- data.frame(group=c( "2aWT","35SBRI","apex1","apex2","apex3","bak1bkk1","bak1bkk1cerk1","bri301", "cyp79b2b3", "deps", "efrfls2cerk1",
                              "lyk5", "pad4",  "rar1","wrky33", "wrky3340"),
                     color=c("black","peachpuff1","thistle","thistle2","thistle4","red2","red4","peachpuff3","blue1","orange3","firebrick3",
                             "darkred","darkgoldenrod1", "green2","purple","purple3"))
colors <- colors[colors$group %in% qPCR_BFO$genotype, ]
qPCR_BFO$genotype<- factor(qPCR_BFO$genotype, levels=colors$group)

qPCR_BFO<-na.omit(qPCR_BFO)

qPCR_BFO <- subset(qPCR_BFO, qPCR_BFO$genotype %in% c("2aWT","35SBRI", "apex1", "apex2", "apex3" ,"bak1bkk1", "bak1bkk1cerk1", "bri301", "cyp79b2b3", "deps", "efrfls2cerk1",
                                                      "lyk5", "pad4",  "rar1", "wrky33", "wrky3340"))


### change names from old version

qPCR_BFO[,"genotype2"]=gsub("apex1","hub1",qPCR_BFO[,"genotype2"])
qPCR_BFO[,"genotype2"]=gsub("apex2","apex",qPCR_BFO[,"genotype2"])
qPCR_BFO[,"genotype2"]=gsub("apex3","hub2",qPCR_BFO[,"genotype2"])
													


l <- c("2aWT","bak1bkk1","bak1bkk1cerk1","efrfls2cerk1","lyk5","hub1","apex", "hub2", "wrky3340","wrky33",
       "deps","pad4","cyp79b2b3","35SBRI","bri301","rar1")
qPCR_BFO$genotype2 <- factor(qPCR_BFO$genotype2, levels=l)

library(ggplot2)
#graph without outliers being marked
p_gen_shannon <- ggplot(qPCR_BFO, aes(x=genotype2, y=X16S_UBQ_ref, color=genotype)) +
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent")+
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x="", y="microbe/plant/ref ratio") +theme_bw()+ggtitle(" Oomycete - relative abundance")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")+theme(axis.text.x=element_text(angle=90,hjust=1))

p_gen_shannon
ggsave(paste(figure.dir, "fig_3_C.pdf", sep=""), p_gen_shannon)

#total_load_flowpot_bacteria_no_kai2quin

###############################################################################################################
###############################################################################################################
#STATISTICS
###############################################################################################################
###############################################################################################################

#qqnorm(qPCR_BFO$X16S_UBQ_ref)
#qqline(qPCR_BFO$X16S_UBQ_ref, col="red")
#shapiro.test(qPCR_BFO$X16S_UBQ_ref)
#############################################################

library("rcompanion")
library("FSA")
library("plyr")

model3 <- function(data){
  capture.output(kruskal.test(data$X16S_UBQ_ref ~ data$genotype))
}

OT <- dlply(qPCR_BFO, .(compartment), model3)
cat("\n", file="Kruskal_Wallis_fig3_c.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="Kruskal_Wallis_fig3_c.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(X16S_UBQ_ref ~ genotype, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(qPCR_BFO, .(compartment), model4b)
write.table(ldply(KW_output), file="Kruskal_Wallis_fig3_c.txt", col.names=TRUE, row.names=TRUE, append=TRUE)



