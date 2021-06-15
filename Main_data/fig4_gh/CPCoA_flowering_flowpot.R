
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")
source("cpcoa_funct_Ruben.R")

############# UNCONSTRAINED PCoA PLOT ##################

library(ggplot2, quietly=T, warn.conflicts=F)


# load data

design <- read.table("fun_data/design_fun.txt", header=T, sep="\t")

otu_table <- read.table("fun_data/count_tab_fun.txt", sep="\t", header=T, check.names=F)




#subset the data #include root_soil samples (for dead plants)
design = subset(design, 
               design$biological %in% c("bio2","bio3") & 
                 design$compartment %in% c("root","root_soil") &
                 design$genotype %in% c("WT", "cyp79b2b3"))
otu_table = subset(otu_table[, as.vector(design$Sample)])

#only samples above 1000 reads
threshold <- 1000
idx <- colSums(otu_table) > threshold 
otu_table <- otu_table[,idx ]



#To normalize (total sum), you can use:
otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))
# and to do unconstrained PCoA


# re-order data matrices
otu_table<-otu_table_norm
idx <- design$Sample %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$Sample, colnames(otu_table))
otu_table <- otu_table[, idx]

library(vegan)
bray_curtis <- vegdist(t(otu_table), method="bray")


##
colors <- data.frame(group=c("cyp79b2b3", "WT"),
          				   	color=c("blue1", "black"))
colors <- colors[colors$group %in% design$genotype, ]

#################unconstr shapes 

#shapes bacteria
shapes <- data.frame(group=c("B","BF","BFO","BO","F","FO","O","initial"),
                     shape=c(15,16,18,17,0,5,1,8))
shapes <- shapes[shapes$group %in% design$treatment, ]


#shapes <- data.frame(group=c("bio1","bio2", "bio3","bio4"),
#                     shape=c(15,18,16,17))
#shapes <- shapes[shapes$group %in% design$biological, ]
#
#shapes <- data.frame(group=c("high","low","medium","not_apply", "unknwon"),
#                     shape=c(16,17,15,0,3))
#shapes <- shapes[shapes$group %in% design$growth_promotion_three, ]

# PCoA Bray-Curtis

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$Sample), ])

#points colors

points$genotype <- factor(points$genotype, levels=colors$group)

#points shapes

points$treatment <- factor(points$treatment, levels=shapes$group)


#############
#points$Treatment <- factor(points$Treatment, levels=shapes$group)
#points$Tissue <- factor(points$Tissue, levels=colors$group)
###########
# plot PCo 1 and 2
pointsnoNA<-na.omit(points)
library(ggplot2)
unconstr_p <- ggplot(pointsnoNA, aes(x=x, y=y, color=genotype, shape=treatment)) +
  geom_point(alpha=.7, size=5) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste("Unconstrained PCoA flowering flowpot"))+  theme_bw() + 
  theme(plot.title = element_text(size = 25), axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        legend.text=element_text(size=10),legend.title=element_text(size=12),
        legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
unconstr_p

figures.dir <-""
ggsave(paste(figures.dir, "pcoa_floweringFP.pdf", sep=""), unconstr_p)




