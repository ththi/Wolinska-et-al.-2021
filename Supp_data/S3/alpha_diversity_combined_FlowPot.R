#
# originally by Ruben Garrido-Oter
#

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")


library(ggplot2, quietly=T, warn.conflicts=F)
library(scales, quietly=T, warn.conflicts=F)
library(grid, quietly=T, warn.conflicts=F)
library(vegan, quietly=T, warn.conflicts=F)



# load data

figures.dir<-""

design <- read.table("../../Main_data/fig2/bac_data/bac_design.txt", header=T, sep="\t")

shannon <- read.table("../../Main_data/fig2/bac_data/bac_chao.txt", sep="\t", header=T, check.names=F)


 
#otu_table = subset(otu_table[, as.vector(design$SampleID)])


design[,"genotype"]=gsub("WT_new","WT",design[,"genotype"])
design[,"genotype"]=gsub("cyp79b2b3_new","cyp79b2b3",design[,"genotype"])
design[,"genotype"]=gsub("initial_new","initial",design[,"genotype"])
design[,"genotype"]=gsub("soil_new","soil",design[,"genotype"])
design[,"genotype_microbes"]=gsub("WT_new_BFO","WT_BFO",design[,"genotype_microbes"])
design[,"genotype_microbes"]=gsub("cyp79b2b3_new_BFO","cyp79b2b3_BFO",design[,"genotype_microbes"])
design[,"genotype_microbes"]=gsub("soil_new_BFO","soil_BFO",design[,"genotype_microbes"])
design[,"genotype_microbes"]=gsub("initial_new","initial",design[,"genotype_microbes"])

#subset the data
design = subset(design, 
                design$genotype %in% c("initial","soil","WT","bak1bkk1","bak1bkk1cerk1", "efrfls2cerk1", "lyk5","apex1","apex2", "apex3", "wrky3340",
                                       "wrky33", "deps","pad4","cyp79b2b3","35SBRI","bri301",  "rar1"))

### alpha diversity
colors <- data.frame(group=c("35SBRI_BFO","apex1_BFO","apex2_BFO","apex3_BFO","bak1bkk1_BFO","bak1bkk1cerk1_BFO","bri301_BFO", "cyp79b2b3_BFO",
                              "deps_BFO", "efrfls2cerk1_BFO",
                             "initial", "kai2_BFO","lyk5_BFO", "pad4_BFO", "quin_BFO", "rar1_BFO","soil_BFO","wrky33_BFO",
                             "wrky3340_BFO", "WT_BFO"),
                     color=c("peachpuff1","thistle","thistle2","thistle4","red2","red4","peachpuff3","blue1","orange3","firebrick3","violetred",
                             "yellow2","darkred","darkgoldenrod1", "deepskyblue1","green2","tan4","purple","purple3", "black"))
colors <- colors[colors$group %in% design$genotype_microbes, ]




index_shannon <- cbind(shannon[, 2], design[match(shannon[, 1], rownames(design)), ])

colnames(index_shannon)[1] <- "value"
index_shannon<-na.omit(index_shannon)


### renaming of genotypes from old version

index_shannon[,"genotype"]=gsub("apex1","hub1",index_shannon[,"genotype"])
index_shannon[,"genotype"]=gsub("apex2","apex",index_shannon[,"genotype"])
index_shannon[,"genotype"]=gsub("apex3","hub2",index_shannon[,"genotype"])


l <- c("initial","soil","WT","bak1bkk1","bak1bkk1cerk1", "efrfls2cerk1", "lyk5","hub1","apex", "hub2", "wrky3340",
       "wrky33", "deps","pad4","cyp79b2b3","35SBRI","bri301",  "rar1")

index_shannon$genotype <- factor(index_shannon$genotype, levels=l)




#colors <- colors[match(l, colors$group), ]
library(ggplot2)
#graph without outliers being marked
p_gen_shannon <- ggplot(index_shannon, aes(x=genotype, y=value, color=genotype_microbes)) +
  geom_boxplot(alpha=1, outlier.size=1, width=0.4, fill="transparent")+
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x="", y="chao index") +theme_bw()+ggtitle("bacteria ")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + ylim(0, 80)+theme(axis.text.x=element_text(angle=90,hjust=1))

p_gen_shannon
ggsave(paste(figures.dir, "fig_S3_x.pdf", sep=""), p_gen_shannon)






library("rcompanion")
library("FSA")
library("plyr")

model3 <- function(data){
  capture.output(kruskal.test(data$value ~ data$genotype))
}

OT <- dlply(index_shannon, .(library_type), model3)
cat("\n", file="kw_S3_x.txt", sep="\n", append=TRUE)
write.table(ldply(OT), file="kw_S3_x.txt", col.names=TRUE, row.names=FALSE, append=TRUE)

model4b <- function(Data){
  OT2 <- dunnTest(value ~ genotype, data=Data, method="bh")
  OT2 <- OT2$res
  OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)
  OT2a
}

KW_output <- dlply(index_shannon, .(library_type), model4b)
write.table(ldply(KW_output), file="kw_S3_x.txt", col.names=TRUE, row.names=TRUE, append=TRUE)




