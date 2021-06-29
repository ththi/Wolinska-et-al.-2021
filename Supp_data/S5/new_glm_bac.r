library("gplots")
library("edgeR")

bc_nam <- read.table("../../Main_data/fig2/bac_data/bac_design.txt", header=T)

otu_tab <- read.table("../../Main_data/fig2/bac_data/bac_count_tab_rem.txt", header=T)

otu_ra_orig=sweep(otu_tab,2,colSums(otu_tab),"/")

### for orig figure 2000 reads
keep_row=rownames(otu_tab[rowSums(otu_tab)>=1000,])
otu_tab=otu_tab[keep_row,]

otu_ra=otu_ra_orig[keep_row,]

# remove strange apex2 sample
bc_nam=bc_nam[-32,]

# change names for "_new" samples
bc_nam[,"genotype"]=gsub("WT_new","WT",bc_nam[,"genotype"])
bc_nam[,"genotype"]=gsub("cyp79b2b3_new","cyp79b2b3",bc_nam[,"genotype"])


tab_ra=sweep(otu_tab,2,colSums(otu_tab),"/")


# use design, choose samples /colors

# order is crucial for sample selection in loop !!!

gen_list=c("35SBRI","apex1","apex2","apex3","bak1bkk1","bak1bkk1cerk1","bri301", "cyp79b2b3", "deps", "efrfls2cerk1","lyk5", "pad4", "rar1","wrky33", "wrky3340","WT")



bc_nam = subset(bc_nam,bc_nam$microbes %in% c("BFO") & bc_nam$compartment %in% c("endophytes") & bc_nam$genotype %in% gen_list)




gen_list=gen_list[-16]

#dev.new(width=6.87500, height=11.16988)

sig_tab_glm=matrix(data=0,nrow=nrow(otu_tab),ncol=15,dimnames=list(c(rownames(otu_tab)),gen_list) )
fc_tab_glm=matrix(data=0,nrow=nrow(otu_tab),ncol=15,dimnames=list(c(rownames(otu_tab)),gen_list) )



#for(i in 1:length(gen_list)){
for(i in 1:15){

	#aa=unique(bc_nam[bc_nam[,"genotype"]==gen_list[i],"experiment_biological"])
	#print(paste(gen_list[i]))
	#print(paste(aa))

	use_row=which( bc_nam[,"genotype"]==gen_list[i] | bc_nam[,"genotype"]=="WT" )
	bc_nam_use=bc_nam[use_row,]


	if(i==8){

		col_use_a1=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A1_bio1",])
		col_use_a2=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A1_bio2",])
		col_use_a3=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A1_bio4",])
		col_use_a4=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A4_bio1",])
		col_use_a5=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A4_bio2",])
		col_use_a6=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A4_bio4",])

		all_col=c(col_use_a1,col_use_a2,col_use_a3,col_use_a4,col_use_a5,col_use_a6)
		
	}

	if(i==5||i==9||i==10||i==12||i==13||i==15){

		col_use_a1=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A1_bio1",])
		col_use_a2=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A1_bio2",])
		col_use_a3=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A1_bio4",])

		all_col=c(col_use_a1,col_use_a2,col_use_a3)
				
	}

	if(i==1||i==2||i==3||i==4||i==6||i==7||i==11||i==14||i==16){

		col_use_a1=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A4_bio1",])
		col_use_a2=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A4_bio2",])
		col_use_a3=rownames(bc_nam_use[bc_nam_use[,"experiment_biological"]=="A4_bio4",])
		
		all_col=c(col_use_a1,col_use_a2,col_use_a3)	
		
	
		
	}

	sub_otu_tab=otu_tab[,all_col]
	
	wt_sam=grep("WT",colnames(sub_otu_tab))
	mut_sam=grep("WT",colnames(sub_otu_tab),invert=T)
	
	sam_list=c(rep(1,length(wt_sam)),rep(2,length(mut_sam)))
	
	sub_otu_tab2=sub_otu_tab[,c(wt_sam,mut_sam)]
	
	sub_bc=bc_nam_use[c(all_col[wt_sam],all_col[mut_sam]),]
	
	new_obj=DGEList(counts=sub_otu_tab2,lib.size=colSums(sub_otu_tab2),group=sam_list,remove.zeros=T)
	new_obj <- calcNormFactors(new_obj)
	
	#exp_num=factor(substring(colnames(sub_otu_tab2),4,6))
	exp_num=factor(sub_bc[,"experiment_biological"])
	cond_num=factor(sam_list)
	
	design=model.matrix(~exp_num+cond_num)
	rownames(design)=colnames(sub_otu_tab2)
	
	new_obj_disp=estimateGLMCommonDisp(new_obj,design,verbose=F)
	new_obj_disp=estimateGLMTrendedDisp(new_obj_disp,design)
	new_obj_disp=estimateGLMTagwiseDisp(new_obj_disp,design)

	fit=glmFit(new_obj_disp,design)
	lrt=glmLRT(fit)
		
	sig_ones=rownames(lrt$table[lrt$table$PValue <=0.05,])
	sig_tab_glm[sig_ones,i]=1
	fc_tab_glm[rownames(lrt$table),i]=lrt$table[,1]


}

bb=sig_tab_glm[,]
bb[bb==1]="*"
bb[bb==0]=NA

#### renaming apex samples
colnames(fc_tab_glm)=gsub("apex1","hub1",colnames(fc_tab_glm))
colnames(fc_tab_glm)=gsub("apex2","apex",colnames(fc_tab_glm))
colnames(fc_tab_glm)=gsub("apex3","hub2",colnames(fc_tab_glm))


dev.new(height=7.7,width=3.5)
use_col=colorRampPalette(c("dark blue","blue","blue","white","red","red","brown"))(100)
heatmap.2(as.matrix(fc_tab_glm[,c(5,6,10,11,2,3,4,15,14,9,12,8,1,7,13)]),col=use_col,tracecol=F,Colv=F,cellnote=bb[,c(5,6,10,11,2,3,4,15,14,9,12,8,1,7,13)],notecol="black",margin=c(8,6),cexCol=0.75,cexRow=0.5)
#heatmap.2(as.matrix(fc_tab_glm[,]),col=use_col,tracecol=F,Colv=T,cellnote=bb[,],notecol="black",margin=c(8,6),cexCol=0.75,cexRow=0.5)

dev.print(pdf,"fig_S5_a.pdf",useDingbats=F)


