
library("vegan")
source("cpcoa.func.R")


#otu_tab <- read.table("bac_data/bac_count_tab_rem.txt", header=T)
#otu_tab <- read.table("fun_data/fun_count_tab_rem.txt", header=T)
otu_tab <- read.table("oo_data/oo_count_tab_rem.txt", header=T)


#note, header of first column should be removed (e.g. SampleID)
#bc_nam <- read.table("bac_data/bac_design.txt", header=T)
#bc_nam <- read.table("fun_data/fun_design.txt", header=T)
bc_nam <- read.table("oo_data/oo_design.txt", header=T)


# remove outlier  apex2 sample - only Bacteria # 
#bc_nam=bc_nam[-32,]

# change names for "_new" samples
bc_nam[,"genotype"]=gsub("WT_new","WT",bc_nam[,"genotype"])
bc_nam[,"genotype"]=gsub("cyp79b2b3_new","cyp79b2b3",bc_nam[,"genotype"])

# change genotype names from old version
bc_nam[,"genotype"]=gsub("apex1","hub1",bc_nam[,"genotype"])
bc_nam[,"genotype"]=gsub("apex2","apex",bc_nam[,"genotype"])
bc_nam[,"genotype"]=gsub("apex3","hub2",bc_nam[,"genotype"])

# use design, choose samples /colors

gen_list=c("35SBRI","hub1","apex","hub2","bak1bkk1","bak1bkk1cerk1","bri301", "cyp79b2b3", "deps", "efrfls2cerk1","lyk5", "pad4", "rar1","wrky33","wrky3340", "WT")
col_list=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","peachpuff3","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1", "green2","purple","purple3", "black")


bc_nam = subset(bc_nam,bc_nam$microbes %in% c("BFO") & bc_nam$compartment %in% c("endophytes") & bc_nam$genotype %in% gen_list)


match_row=match(rownames(bc_nam),colnames(otu_tab))

#oomycetes, some NAs are present in match_row and was giving error from next line on without na.omit
match_row<-na.omit(match_row)

otu_fin=otu_tab[,match_row]
use_row=intersect(colnames(otu_fin),rownames(bc_nam))

# for oomycetes
bc_nam=bc_nam[use_row,]

otu_ra=sweep(otu_fin,2,colSums(otu_fin),"/")


cs=capscale(t(otu_ra[,use_row])~bc_nam[use_row,"genotype"]+Condition(bc_nam[use_row,"experiment"])+Condition(bc_nam[use_row,"biological"])+Condition(bc_nam[use_row,"library"]),add=F,sqrt.dist=T,distance="bray") ### otu tab 




term=anova.cca(cs)
p.val <- term[1, 4]
eig=cs$CCA$eig
var1=format(100 * eig[1] / sum(eig))
var2=format(100 * eig[2] / sum(eig))

var_tab=variability_table(cs)
variance=var_tab["constrained", "proportion"]
expl_var=format(100 * variance, digits=3)

par(mar=c(6,6,6,10))
plot(cs$CCA$wa,xlab=var1,ylab=var2,main=paste(expl_var," pval=",p.val),col="white",cex=0.5)


#for(i in 1:length(gen_list)){
  
  ### use these if individual points should be plotted	
  
  #aa_1=which(bc_nam[,"genotype"]==gen_list[i])
  #points(cs$CCA$wa[aa_1,],col=col_list[i],pch=19)
  
  #}



x_max=max(cs$CCA$wa[,1])+0.1
x_min=min(cs$CCA$wa[,1])-0.1

y_max=max(cs$CCA$wa[,2])+0.1
y_min=min(cs$CCA$wa[,2])-0.1

legend(x_max,y_max,gen_list,xpd=T,fill=col_list)


###for min_max


for(i in 1:length(gen_list)){
  
  use_row=which(bc_nam[,"genotype"]==gen_list[i])
  
  x_max=max(cs$CCA$wa[use_row,1])
  x_min=min(cs$CCA$wa[use_row,1])
  x_mean=mean(cs$CCA$wa[use_row,1])
  
  y_max=max(cs$CCA$wa[use_row,2])
  y_min=min(cs$CCA$wa[use_row,2])
  y_mean=mean(cs$CCA$wa[use_row,2])
  
  # plot x
  segments(x_min,y_mean,x_max,y_mean,col=col_list[i])
  
  # plot y
  segments(x_mean,y_min,x_mean,y_max,col=col_list[i])
  
}

dev.print(pdf,"fig2_f.pdf",useDingbats=F)





