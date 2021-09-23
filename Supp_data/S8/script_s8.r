library("FSA")
library("rcompanion")

tab=read.table("tab_with_aux_S8.txt",header=T)

dev.new(width=9.095022, height=6.597285)
par(mfrow=c(3,4),cex.main=0.75)
col_list=c("black","black","black","#2156a6","#2156a6","#2156a6")

for(i in 4:15){
	
	boxplot(tab[,i]~tab[,3],las=2,xlab="",main=colnames(tab)[i],border=col_list,at=c(2,3,1,5,6,4))
	stripchart(tab[,i]~tab[,3],add=T,vertical=T,method="jitter",pch=19,cex=0.5,col=col_list,at=c(2,3,1,5,6,4))
	aa=kruskal.test(tab[,i]~tab[,3])
	
	write("\n",file="stat_file.txt",append=TRUE)
	write(paste(colnames(tab)[i],aa$p.value,sep=" "),file="stat_file.txt",append=TRUE)
	
	
	
	if(aa$p.value <=0.05){
		OT1=dunnTest(tab[,i]~tab[,3],method="bh")
		OT1 <- OT1$res
		OT1a <- cldList(comparison = OT1$Comparison, p.value = OT1$P.adj, threshold = 0.05)
		write(as.matrix(OT1a[c(3,1,2,6,4,5),1:2]),file="stat_file.txt",append=TRUE)
	}	
}



