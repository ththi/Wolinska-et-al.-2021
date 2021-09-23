#plotting the total fungal biomass relative to root biomass




qPCR_test <- read.table( "A5_oomycetes_total_load_rem.txt", header=T, sep="\t")

###

qPCR_BFO = subset(qPCR_test, 
                  qPCR_test$compartment %in% c("root")) 
qPCR_BFO<-droplevels(qPCR_BFO)



qPCR_BFO<-na.omit(qPCR_BFO)

l <- c("WT","cyp79b2_b3",
       "cyp71a12_a13","pen2_pad3","pen2","pad3","myb_triple","pyk_bglu",
        "cyp71a27")
qPCR_BFO$genotype <- factor(qPCR_BFO$genotype, levels=l)


no_pencyp_quin <- subset(qPCR_BFO, qPCR_BFO$genotype %in% c("WT", "cyp79b2_b3", "cyp71a12_a13", 
                                                       "pen2_pad3", "pen2", "pad3", "myb_triple", "pyk_bglu", "cyp71a27"))



library(ggplot2)
#graph without outliers being marked
p_gen_shannon <- ggplot(no_pencyp_quin, aes(x=genotype, y=X16S_UBQ_ref)) +
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent")+
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7) +
  #scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  #coord_cartesian(ylim = c(0, 40)) +
  labs(x="", y="oomycete/plant/ref ratio") +theme_bw()+ggtitle("oomycete realtive abundance")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")+theme(axis.text.x=element_text(angle=90,hjust=1))

p_gen_shannon
ggsave("s10_f.pdf", p_gen_shannon)



qPCR_BFO[,"genotype"]=gsub("WT","1aWT",qPCR_BFO[,"genotype"])
qPCR_BFO$genotype <- as.factor(qPCR_BFO$genotype)

no_pencyp <- subset(qPCR_BFO, qPCR_BFO$genotype %in% c("1aWT", "cyp79b2_b3", "quintuple", "cyp71a12_a13", 
                                                       "pen2_pad3", "pen2", "pad3", "myb_triple", "pyk_bglu", "cyp71a27"))

no_pencyp_quin <- subset(qPCR_BFO, qPCR_BFO$genotype %in% c("1aWT", "cyp79b2_b3", "cyp71a12_a13", 
                                                            "pen2_pad3", "pen2", "pad3", "myb_triple", "pyk_bglu", "cyp71a27"))


#qqnorm(no_pencyp_quin$oITS_UBQ_ref)
#qqline(no_pencyp_quin$oITS_UBQ_ref, col="red")
#shapiro.test(no_pencyp_quin$oITS_UBQ_ref)
#############################################################

library("rcompanion")
library("FSA")
library("plyr")

model3 <- function(data){
  capture.output(kruskal.test(data$X16S_UBQ_ref ~ data$genotype))
}

OT <- dlply(no_pencyp_quin, .(compartment), model3)
cat("\n", file="KW_s10_f.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="KW_s10_f.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(X16S_UBQ_ref ~ genotype, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(no_pencyp_quin, .(compartment), model4b)
write.table(ldply(KW_output), file="KW_s10_f.txt", col.names=TRUE, row.names=TRUE, append=TRUE)



