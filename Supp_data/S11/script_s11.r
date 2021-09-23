#library("FSA")
#library("rcompanion")

tab=read.table("table_S_11.txt",header=T)

dev.new(width=11.45022, height=7.597285)
par(mfrow=c(3,4),cex.main=0.75)

tab$WTgenotype=gsub("_"," ",tab$WTgenotype)

for(i in 2:10){
	
	
	aa=kruskal.test(tab[,i]~tab[,1])
	
	write("\n",file="stat_file.txt",append=TRUE)
	write(paste(colnames(tab)[i],aa$p.value,sep=" "),file="stat_file.txt",append=TRUE)
	
	
	
	if(aa$p.value <=0.05){
		OT1=dunnTest(tab[,i]~tab[,1],method="bh")
		OT1 <- OT1$res
		OT1a <- cldList(comparison = OT1$Comparison, p.value = OT1$P.adj, threshold = 0.05)
		write(as.matrix(OT1a[c(2,1,5,3,9,8,7,6,10,4),1:2]),file="stat_file.txt",append=TRUE)
	}	
	
	stat_vec=OT1a[c(2,1,5,3,9,8,7,6,10,4),2]
	title_stat=paste(colnames(tab)[i],stat_vec,sep="\n")
	
	boxplot(tab[,i]~tab[,1],las=2,xlab="",main=paste("xxxxxx",paste(stat_vec,collapse="    "),sep="\n"),border="black",at=c(2,1,4,10,3,8,7,6,5,9))
	stripchart(tab[,i]~tab[,1],add=T,vertical=T,method="jitter",pch=19,cex=0.5,col="black",at=c(2,1,4,10,3,8,7,6,5,9))
	
	#boxplot(tab[,i]~tab[,1],las=2,xlab="",main=paste("xxxxxx",paste(stat_vec,collapse="    "),sep="\n"),border="black")
	#stripchart(tab[,i]~tab[,1],add=T,vertical=T,method="jitter",pch=19,cex=0.5,col="black")
	
	
}



