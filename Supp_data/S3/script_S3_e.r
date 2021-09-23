


tab=read.csv("week_5_ratio_S3_e.txt",header=T)

dev.new(width=7.520362 , height=6.000000)

xx=by(tab[,4],tab[,5],median)
zz=order(xx)
aa=c(1:17)
new_ord=match(aa,zz)
boxplot(tab[,4]~tab[,5],las=2,log="y",at=new_ord,xlab="",ylab="Relative fungal/plant/ref ratio",main="Greenhouse (CAS) fungal load")


#stat testing

kruskal.test(tab[,4]~tab[,5])

### no more testing, since NS

