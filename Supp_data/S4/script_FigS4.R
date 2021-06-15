#####Loading libraries#####
library(phyloseq)
library(dplyr)
library(car)
library(ggplot2)
##################
#####Bacteria#####
##################

#####Loading datasets, filtering and renaming#####
biomass=read.table("FW_microbial_load.txt", h=T)
biomass$genotype[biomass$genotype=="bak1bkk1"]<-"bak1/bkk1"
biomass$genotype[biomass$genotype=="bak1bkk1cerk1"]<-"bak1/bkk1/cerk1"

design=read.table("design_copy_bacteria_for_CPCoA.txt", h=T, sep="\t", row.names = 1)
design$genotype[design$genotype=="WT_new"]<-"WT"
design$genotype[design$genotype=="cyp79b2b3_new"]<-"cyp79b2b3"
design$genotype[design$genotype=="bak1bkk1"]<-"bak1/bkk1"
design$genotype[design$genotype=="bak1bkk1cerk1"]<-"bak1/bkk1/cerk1"
design=design[which(design$genotype%in%biomass$genotype),]
design=design[-which(design$genotype=="bri301"),]
design$genotype<-droplevels(as.factor(design$genotype))

otu=read.csv("otu_tab_filter_98_bacteria.txt", sep="")
otu_filtered=otu[,which(colnames(otu)%in% rownames(design))]

#####Making the PCOA anf extracting the coordinates on PCOA1 and PCOA2#####
samples = sample_data(design)
OTU = otu_table(otu_filtered, taxa_are_rows = TRUE)
GP1 <- phyloseq(OTU,samples)
GP1_norm  = transform_sample_counts(GP1, function(x) x / sum(x) )
GP1_filter=filter_taxa(GP1_norm, function(x) mean(x) > 0.001, TRUE)
GP1_filter_merged_mutants=merge_samples(GP1_norm,"genotype")
GP1_filter_merged_mutants= transform_sample_counts(GP1_filter_merged_mutants, function(x) x / sum(x) )
GP.ord <- ordinate(GP1_filter_merged_mutants, "MDS", "bray")

p1 = plot_ordination(GP1_filter_merged_mutants, GP.ord, type="samples",axes = c(1,2))
p1

x=GP.ord$vectors
x=as.matrix(x)
x=x[,1]
x=as.data.frame(x)
x$names<-rownames(x)
y=GP.ord$vectors
y=as.matrix(y)
y=y[,2]
y=as.data.frame(y)
y$names<-rownames(y)

colnames(biomass)[1]<-"names"
tab=inner_join(biomass, x, by="names")
tab=inner_join(tab, y, by="names")
colnames(tab)[2]<-"growth_promotion"
colnames(tab)[1]<-"genotype"

####Regression between Growth promotion and PCOA1####
mod=lm(growth_promotion~x, data=tab)
qqnorm(residuals(mod))
qqline(residuals(mod))
shapiro.test(residuals(mod))
Anova(mod)
summary(mod)
p1 <- ggplot(tab, aes(x=x, y= growth_promotion)) +labs(x="PCOA1", y="Mean Relative FW")+
  geom_point(aes(colour=genotype)) +
  scale_colour_manual(values=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1","green2","purple","purple3", "black"))+
  geom_smooth(method=lm, se=TRUE,  color="black")+
  ggtitle("p-value=0.7299, R2=-0.07231 ")
p1
####Regression between Growth promotion and PCOA2####
mod=lm(growth_promotion~y, data=tab)
qqnorm(residuals(mod))
qqline(residuals(mod))
shapiro.test(residuals(mod))
Anova(mod)
summary(mod)
p2 <- ggplot(tab, aes(x=y, y= growth_promotion)) +labs(x="PCOA2", y="Mean Relative FW")+
  geom_point(aes(colour=genotype)) +
  scale_colour_manual(values=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1","green2","purple","purple3", "black"))+
  geom_smooth(method=lm, se=TRUE,  color="black")+
  ggtitle("p-value=0.09106, R2=0.1427")
p2


###############
#####FUNGI#####
###############

#####Loading datasets, filtering and renaming#####
biomass=read.table("FW_microbial_load.txt", h=T)
biomass$genotype[biomass$genotype=="bak1bkk1"]<-"bak1/bkk1"
biomass$genotype[biomass$genotype=="bak1bkk1cerk1"]<-"bak1/bkk1/cerk1"

design=read.table("design_copy_fungi_for_CPCoA.txt", h=T, sep="\t", row.names = 1)
design$genotype[design$genotype=="WT_new"]<-"WT"
design$genotype[design$genotype=="cyp79b2b3_new"]<-"cyp79b2b3"
design$genotype[design$genotype=="bak1bkk1"]<-"bak1/bkk1"
design$genotype[design$genotype=="bak1bkk1cerk1"]<-"bak1/bkk1/cerk1"
design=design[which(design$genotype%in%biomass$genotype),]
design=design[-which(design$genotype=="bri301"),]
design$genotype<-droplevels(as.factor(design$genotype))

otu=read.csv("otu_tab_filter_98_fungi.txt", sep="")
otu_filtered=otu[,which(colnames(otu)%in% rownames(design))]

#####Making the PCOA anf extracting the coordinates on PCOA1 and PCOA2#####
samples = sample_data(design)
OTU = otu_table(otu_filtered, taxa_are_rows = TRUE)
GP1 <- phyloseq(OTU,samples)
GP1_norm  = transform_sample_counts(GP1, function(x) x / sum(x) )
GP1_filter=filter_taxa(GP1_norm, function(x) mean(x) > 0.001, TRUE)
GP1_filter_merged_mutants=merge_samples(GP1_norm,"genotype")
GP1_filter_merged_mutants= transform_sample_counts(GP1_filter_merged_mutants, function(x) x / sum(x) )
GP.ord <- ordinate(GP1_filter_merged_mutants, "MDS", "bray")

p1 = plot_ordination(GP1_filter_merged_mutants, GP.ord, type="samples",axes = c(1,2))
p1

x=GP.ord$vectors
x=as.matrix(x)
x=x[,1]
x=as.data.frame(x)
x$names<-rownames(x)
y=GP.ord$vectors
y=as.matrix(y)
y=y[,2]
y=as.data.frame(y)
y$names<-rownames(y)

colnames(biomass)[1]<-"names"
tab=inner_join(biomass, x, by="names")
tab=inner_join(tab, y, by="names")
colnames(tab)[2]<-"growth_promotion"
colnames(tab)[1]<-"genotype"

####Regression between Growth promotion and PCOA1####
mod=lm(growth_promotion~x, data=tab)
qqnorm(residuals(mod))
qqline(residuals(mod))
shapiro.test(residuals(mod))
Anova(mod)
summary(mod)
p1 <- ggplot(tab, aes(x=x, y= growth_promotion)) +labs(x="PCOA1", y="Mean Relative FW")+
  geom_point(aes(colour=genotype)) +
  scale_colour_manual(values=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1","green2","purple","purple3", "black"))+
  geom_smooth(method=lm, se=TRUE,  color="black")+
  ggtitle("p-value=0.894, R2=-0.0754 ")
p1
####Regression between Growth promotion and PCOA2####
mod=lm(growth_promotion~y, data=tab)
qqnorm(residuals(mod))
qqline(residuals(mod))
shapiro.test(residuals(mod))
Anova(mod)
summary(mod)
p2 <- ggplot(tab, aes(x=y, y= growth_promotion)) +labs(x="PCOA2", y="Mean Relative FW")+
  geom_point(aes(colour=genotype)) +
  scale_colour_manual(values=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1","green2","purple","purple3", "black"))+
  geom_smooth(method=lm, se=TRUE,  color="black")+
  ggtitle("p-value=0.1694, R2=0.07387")
p2

###################
#####OOMYCETES#####
###################

#####Loading datasets, filtering and renaming#####
biomass=read.table("FW_microbial_load.txt", h=T)
biomass$genotype[biomass$genotype=="bak1bkk1"]<-"bak1/bkk1"
biomass$genotype[biomass$genotype=="bak1bkk1cerk1"]<-"bak1/bkk1/cerk1"

design=read.table("design_copy_oomycetes_for_CPCoA.txt", h=T, sep="\t", row.names = 1)
design$genotype[design$genotype=="WT_new"]<-"WT"
design$genotype[design$genotype=="cyp79b2b3_new"]<-"cyp79b2b3"
design$genotype[design$genotype=="bak1bkk1"]<-"bak1/bkk1"
design$genotype[design$genotype=="bak1bkk1cerk1"]<-"bak1/bkk1/cerk1"
design=design[which(design$genotype%in%biomass$genotype),]
design=design[-which(design$genotype=="bri301"),]
design$genotype<-droplevels(as.factor(design$genotype))

otu=read.csv("otu_tab_filter_98_oomycetes.txt", sep="")
otu_filtered=otu[,which(colnames(otu)%in% rownames(design))]

#####Making the PCOA anf extracting the coordinates on PCOA1 and PCOA2#####
samples = sample_data(design)
OTU = otu_table(otu_filtered, taxa_are_rows = TRUE)
GP1 <- phyloseq(OTU,samples)
GP1_norm  = transform_sample_counts(GP1, function(x) x / sum(x) )
GP1_filter=filter_taxa(GP1_norm, function(x) mean(x) > 0.001, TRUE)
GP1_filter_merged_mutants=merge_samples(GP1_norm,"genotype")
GP1_filter_merged_mutants= transform_sample_counts(GP1_filter_merged_mutants, function(x) x / sum(x) )
GP.ord <- ordinate(GP1_filter_merged_mutants, "MDS", "bray")

p1 = plot_ordination(GP1_filter_merged_mutants, GP.ord, type="samples",axes = c(1,2))
p1

x=GP.ord$vectors
x=as.matrix(x)
x=x[,1]
x=as.data.frame(x)
x$names<-rownames(x)
y=GP.ord$vectors
y=as.matrix(y)
y=y[,2]
y=as.data.frame(y)
y$names<-rownames(y)

colnames(biomass)[1]<-"names"
tab=inner_join(biomass, x, by="names")
tab=inner_join(tab, y, by="names")
colnames(tab)[2]<-"growth_promotion"
colnames(tab)[1]<-"genotype"

####Regression between Growth promotion and PCOA1####
mod=lm(growth_promotion~x, data=tab)
qqnorm(residuals(mod))
qqline(residuals(mod))
shapiro.test(residuals(mod))
Anova(mod)
summary(mod)
p1 <- ggplot(tab, aes(x=x, y= growth_promotion)) +labs(x="PCOA1", y="Mean Relative FW")+
  geom_point(aes(colour=genotype)) +
  scale_colour_manual(values=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1","green2","purple","purple3", "black"))+
  geom_smooth(method=lm, se=TRUE,  color="black")+
  ggtitle("p-value=0.3404, R2=-0.001476  ")
p1
####Regression between Growth promotion and PCOA2####
mod=lm(growth_promotion~y, data=tab)
qqnorm(residuals(mod))
qqline(residuals(mod))
shapiro.test(residuals(mod))
Anova(mod)
summary(mod)
p2 <- ggplot(tab, aes(x=y, y= growth_promotion)) +labs(x="PCOA2", y="Mean Relative FW")+
  geom_point(aes(colour=genotype)) +
  scale_colour_manual(values=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1","green2","purple","purple3", "black"))+
  geom_smooth(method=lm, se=TRUE,  color="black")+
  ggtitle("p-value=0.944, R2=-0.0765")
p2
