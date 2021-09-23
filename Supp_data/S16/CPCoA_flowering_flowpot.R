


options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")
source("cpcoa_funct_Ruben.R")

############# UNCONSTRAINED PCoA PLOT ##################

library(ggplot2, quietly=T, warn.conflicts=F)




# load data

design <- read.table("../fig4_gh/bac_data/design_bac.txt", header=T, sep="\t")

otu_table <- read.table("../fig4_gh/bac_data/count_tab_bac.txt", sep="\t", header=T, check.names=F)




#only samples above 1000 reads
threshold <- 1000
idx <- colSums(otu_table) > threshold 
otu_table <- otu_table[,idx ]



#subset design file
design = subset(design, 
                design$biological %in% c("bio2","bio3") &
                  design$genotype %in% c("WT") &
                  design$compartment %in% c("root","root_soil"))
design<-droplevels(design)



#To normalize (total sum), you can use:
otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))


# re-order data matrices
otu_table<-otu_table_norm
idx <- design$Sample %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$Sample, colnames(otu_table))
otu_table <- otu_table[, idx]

#bray_curtis <- vegdist(t(otu_table), method="bray")

#################constr colors 
###NEW COLORS


colors <- data.frame(group=c("soil","endophytes"),
                     color=c("tan4", "black"))
colors <- colors[colors$group %in% design$compartment, ]




#shapes bacteria
shapes <- data.frame(group=c("B","BF","BFO","BO","F","FO","O","initial"),
                     shape=c(15,16,18,17,0,5,1,8))
shapes <- shapes[shapes$group %in% design$treatment, ]










sqrt_transform <- T

d <- design


capscale.gen <- capscale(t(otu_table) ~ treatment + Condition(library+biological), data=d, add=F, sqrt.dist=sqrt_transform, distance="bray")



# ANOVA-like permutation analysis

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen <- variability_table(capscale.gen)

eig <- capscale.gen$CCA$eig

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

# extract the weighted average (sample) scores

points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$Sample), ])

##########################
# WHOLE DATASET ############################
##########################

p_gen <- ggplot(points, aes(x=x, y=y, color=genotype, shape=treatment)) +
  geom_point(alpha=.7, size=4) +
  #scale_colour_manual(values=as.character(colors$color)) +
  #scale_shape_manual(values=shapes$shape)+
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(" genotype x treatment FlowerFP fungi ",format(100 * variance, digits=3), " % of variance; p=",
                format(p.val, digits=2),
                sep=""))+  theme_bw() + theme(plot.title = element_text(size = 25), axis.title=element_text(size=18),
                                              axis.text=element_text(size=14),
                                              legend.text=element_text(size=18),legend.title=element_text(size=19)
                                              #legend.position = "none"
                                              ,panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p_gen
ggsave( "fig_s_16.pdf", p_gen)


