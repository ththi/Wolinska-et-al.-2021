
library("vegan")
source("cpcoa.func.R")


otu_tab <- read.table("oo_data/oo_count_tab_rem.txt", header=T)

#note, header of first column should be removed (e.g. SampleID)
bc_nam <- read.table("oo_data/oo_design.txt", header=T)



# change names for "_new" samples
bc_nam[,"genotype"]=gsub("WT_new","WT",bc_nam[,"genotype"])
bc_nam[,"genotype"]=gsub("bulk_new","bulk",bc_nam[,"genotype"])



# use design, choose samples /colors
gen_list=c("35SBRI","apex1","apex2","apex3","bak1bkk1","bak1bkk1cerk1","bri301", "cyp79b2b3", "deps", "efrfls2cerk1","lyk5", "pad4", "rar1","wrky33","wrky3340", "WT")
col_list=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","peachpuff3","blue1","yellow2","deeppink1","darkgoldenrod1", "deepskyblue1","green2","purple","purple2","purple3", "black")


bc_nam = subset(bc_nam,bc_nam$microbes %in% c("BFO") & bc_nam$compartment %in% c("endophytes") & bc_nam$genotype %in% gen_list)


match_row=match(rownames(bc_nam),colnames(otu_tab))

#oomycetes, some NAs are present in match_row and was giving error from next line on without na.omit
match_row<-na.omit(match_row)

otu_fin=otu_tab[,match_row]
use_row=intersect(colnames(otu_fin),rownames(bc_nam))

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



for(i in 1:length(gen_list)){
  
  aa_1=which(bc_nam[use_row,"genotype"]==gen_list[i])	
  
  if(gen_list[i]=="WT"){
    
    points(cs$CCA$wa[aa_1,],col=col_list[i],pch=19)
    next
  }
  
  ### do analysis only for current genotype & WT
  
  bc_nam_sub = subset(bc_nam,bc_nam$microbes %in% c("BFO") & bc_nam$compartment %in% c("endophytes") & bc_nam$genotype %in% gen_list[c(i,16)] )
  
  use_row=intersect(colnames(otu_fin),rownames(bc_nam_sub))
  
  cs_x=capscale(t(otu_ra[,use_row])~bc_nam_sub[use_row,"genotype"]+Condition(bc_nam_sub[use_row,"experiment"])+Condition(bc_nam_sub[use_row,"biological"])+Condition(bc_nam_sub[use_row,"library"]),add=F,sqrt.dist=T,distance="bray")
  
  term=anova.cca(cs_x)
  p.val <- term[1, 4]
  var_tab=variability_table(cs_x)
  variance=var_tab["constrained", "proportion"]
  expl_var=format(100 * variance, digits=3)
  
  print(paste(expl_var," pval=",p.val," ",gen_list[i]))
  
  if(p.val <=0.05){
    points(cs$CCA$wa[aa_1,],col=col_list[i],pch=19)
  }else{
    ### print non significant ones
    #points(cs$CCA$wa[aa_1,],col="grey",pch=19)
  }
  
}




