
tab1=read.table("tab_1.txt",header=T)

tab2=read.table("tab_2.txt",header=T)

dev.new(width=9.43,height=5.53)
par(mfrow=c(1,2),mar=c(7,4,4,4))

col_list=c("black","black","black","black","#2156A6","#2156A6","#2156A6","#2156A6")

boxplot(tab1[,"FW_g"]~tab1$treat+tab1$gen,las=2,at=c(2,4,3,1,6,8,7,5),ylim=c(0,0.21),border=col_list,ylab="Shoot FW [g]",xlab="",main="Fungal SynCom 1")
stripchart(tab1[,"FW_g"]~tab1$treat+tab1$gen,add=T,vertical=T,method="jitter",pch=19,cex=0.5,col=col_list,at=c(2,4,3,1,6,8,7,5),ylim=c(0,0.21))

boxplot(tab2[,"FW_g"]~tab2$treat+tab2$gen,las=2,at=c(2,4,3,1,6,8,7,5),ylim=c(0,0.21),border=col_list,ylab="Shoot FW [g]",xlab="",main="Fungal SynCom 2")
stripchart(tab2[,"FW_g"]~tab2$treat+tab2$gen,add=T,vertical=T,method="jitter",pch=19,cex=0.5,col=col_list,at=c(2,4,3,1,6,8,7,5),ylim=c(0,0.21))

dev.print(pdf,"fig_s_15.pdf",useDingbats=F)

### stat tests

new_tab=cbind(tab1,paste(tab1[,"treat"],tab1[,"gen"]))
new_tab2=cbind(tab2,paste(tab2[,"treat"],tab2[,"gen"]))

new_tab[,8]=gsub("-","_",new_tab[,8])
new_tab2[,8]=gsub("-","_",new_tab2[,8])

kruskal.test(tab1[,"FW_g"]~new_tab[,8])

kruskal.test(tab2[,"FW_g"]~new_tab2[,8])

library("FSA")
library("rcompanion")

### stat results will be stored in OT1 and OT2

OT1=dunnTest(new_tab[,"FW_g"]~new_tab[,8],method="bh")
OT1 <- OT1$res
OT1a <- cldList(comparison = OT1$Comparison, p.value = OT1$P.adj, threshold = 0.05)


OT2=dunnTest(new_tab2[,"FW_g"]~new_tab2[,8],method="bh")
OT2 <- OT2$res
OT2a <- cldList(comparison = OT2$Comparison, p.value = OT2$P.adj, threshold = 0.05)