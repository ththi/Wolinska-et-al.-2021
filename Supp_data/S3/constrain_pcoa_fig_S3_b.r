


library("vegan")
source("cpcoa.func.R")

otu_tab <- read.table("otu_table_oo_ra.txt", header=T)

bc_nam <- read.table("design_greenh_oo.txt", header=T)


# change names for "_new" samples
bc_nam[,"genotype"]=gsub("WT_new","WT",bc_nam[,"genotype"])



# use design, choose samples /colors

gen_list=c("35SBRI","apex1","apex2","apex3","bb","bbc","bri301", "cyp79", "deps", "efc","lyk5", "pad4", "rar1","wrky33","wrky40", "wrky3340", "WT")
col_list=c("peachpuff1","thistle","thistle3","thistle4","red2","red4","peachpuff3","blue1","orange3","firebrick3","deeppink1","darkgoldenrod1", "green2","purple","purple2","purple3", "black")





# only tp1 is used !
bc_nam = subset(bc_nam,bc_nam$time_point %in% c("tp1") & bc_nam$compartment %in% c("endophytes") & bc_nam$genotype %in% gen_list)


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
  
  ### use these if individual points should be plotted	
  
  #aa_1=which(bc_nam[use_row,"genotype"]==gen_list[i])
  #points(cs$CCA$wa[aa_1,],col=col_list[i],pch=19)
  
}



x_max=max(cs$CCA$wa[,1])+0.1
x_min=min(cs$CCA$wa[,1])-0.1

y_max=max(cs$CCA$wa[,2])+0.1
y_min=min(cs$CCA$wa[,2])-0.1

legend(x_max,y_max,gen_list,xpd=T,fill=col_list)


###for min_max

use_row_old=use_row

for(i in 1:length(gen_list)){
  
  use_row=which(bc_nam[use_row_old,"genotype"]==gen_list[i])
  
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

###for standard deviation

#for(i in 1:length(gen_list)){

#	use_row=grep(gen_list[i],rownames(cs$CCA$wa))

#	x_sd=sd(cs$CCA$wa[use_row,1])
#	y_sd=sd(cs$CCA$wa[use_row,2])

#	x_mean=mean(cs$CCA$wa[use_row,1])
#	y_mean=mean(cs$CCA$wa[use_row,2])

# plot x
#segments(x_mean-x_sd,y_mean,x_mean+x_sd,y_mean,col=col_list[i],lty=3)

# plot y
#segments(x_mean,y_mean-y_sd,x_mean,y_sd+y_mean,col=col_list[i],lty=3)
#}





