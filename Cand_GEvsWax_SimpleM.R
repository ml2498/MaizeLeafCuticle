### gene expression and cuticular profile: 250K all candidates
###########################################
##### First run Part 1 to get the gene list around GWAS hits
##########################################################
gene_DE.all<-read.table("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\RNAseq_ReadDepth\\RNAseq_ReadDepth.txt",
                        header=T,sep="\t",stringsAsFactors=FALSE)
gene_DE.all<-read.table("/Users/Meng/Google Drive/MLC_AZ_2017/fromPengfei/RNAseq_ReadDepth/RNAseq_ReadDepth.txt",
                        header=T,sep="\t",stringsAsFactors=FALSE)
gene_DE.all<-gene_DE.all[,-seq(2,43,by=2)]
gene_DE.all<-gene_DE.all[,order(colnames(gene_DE.all))]
rownames(gene_DE.all)<-gene_DE.all[,1]
gene_DE.all<-gene_DE.all[,-1]

##
cuticle<-read.csv("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\RNAseq_ReadDepth\\Cuticle_PQ_new.csv",
                  header=T,stringsAsFactors=FALSE)
cuticle<-read.csv("/Users/Meng/Google Drive/MLC_AZ_2017/fromPengfei/RNAseq_ReadDepth/Cuticle_PQ_new.csv",
                  header=T,stringsAsFactors=FALSE)
cuticle<-cuticle[order(cuticle$Sample),]

##
all_genes<-read.table("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\wrapper_cand\\all_100_cand_genes_250kb.txt",header=T,sep="\t")
######################################
# wax components were correlated, so the tests were not independent
# apply modified simpleM method to obtain the number of effective tests
######################################
# 1. correlation matrix
# 2. eigen values of the correlation matrix (PCA)
# 3. M(eff), number of eigen vectors that explained >99.5% total variation
# 4. p_adj=0.05/M(eff)
##########################################
library(rrBLUP)
wax_cor<-cor(cuticle[,-1], method = "pearson", use = "complete.obs")
eigenPCA <- eigen(wax_cor)
eigen_val<-eigenPCA$values
accum_var=0
i=1

while (accum_var<0.995*sum(eigen_val)){
  accum_var=accum_var+eigen_val[i]
  M_eff=i ## number of effective tests across all 41 wax components
  i=i+1
}


## parse out the names of wax components
temp <- file("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\RNAseq_ReadDepth\\Cuticle_PQ_new.txt","r")
chem_names <- readLines(temp,n=1)
close(temp)
chem_names <-strsplit(chem_names,"\t")
chem_names<-matrix(unlist(chem_names),ncol=1,byrow=T)
chem_names<-chem_names[-1,]

###
#candidates<-read.delim("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\GAPIT_cand_250K_v3AtRice_0830_1.txt",header=T,sep="\t")
#candidates2<-c("GRMZM2G032528","GRMZM6G133375","GRMZM5G892627","AC148152.3_FG008","AC210976.4_FG017",
#               "GRMZM2G034153","GRMZM2G070304","GRMZM2G413835") # 1:3 are genes on chr7, the rest are orthologs within the 250kb window

#############################
# candidate gene - wax correlation
##############################
#library(qvalue)
all_genes<-read.table("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\wrapper_cand\\all_100_cand_genes_250kb.txt",header=T,sep="\t")
POS<-unique(all_genes$Position)


nm_wax<-ncol(cuticle)
nm_gene<-nrow(all_genes)

####################################
## consider all the test together, not regionally
################################################
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\RNAseq_ReadDepth")

gene_DE<-gene_DE.all[which(rownames(gene_DE.all) %in% all_genes$v3_gene_model),]

gene_name<-vector()
wax_name<-vector()
r2<-vector()
r<-vector()
pvalue<-vector()
Chr=vector()
Position=vector()
FDR=vector()

nm_tested_gene<-100-length(which(rowSums(gene_DE !=0)==0)) #74 genes with non-zero expression were tested
p_adj<-0.05/(nm_tested_gene*M_eff) #5.63e-05

pdf("gene_DE_vs_cuticle_comp_wrapper_simpleM_02252019.pdf",width=4.5,height=4.5)
for (i in 1:nm_gene){
   gene_DE_1<-t(gene_DE[i,])
   chr=all_genes$Chromosome[which(all_genes$v3_gene_model==rownames(gene_DE)[i])]
   pos=all_genes$Position[which(all_genes$v3_gene_model==rownames(gene_DE)[i])]
   if (any(gene_DE_1!=0)){
     for (j in 2:nm_wax){

       wax<-cuticle[,j]
       p.1<-cor.test(wax,gene_DE_1[,1])$p.value
       pvalue<-c(pvalue,p.1)

       r.1<-cor(wax,gene_DE_1[,1],use="complete.obs")
       r<-c(r,r.1)
       r2.1<-r.1^2
       r2<-c(r2,r2.1)
       gene_name<-c(gene_name,rownames(gene_DE)[i])
       wax_name<-c(wax_name,chem_names[j-1])
       Chr=c(Chr,chr)
       Position=c(Position,pos)



       if (p.1<p_adj){
          plot(wax~gene_DE_1[,1],pch=19,col=rep(c(1:6,"darkorange"),3), xlab="Normalized gene expression level",
               ylab=chem_names[j-1],
               main=paste(rownames(gene_DE)[i]," vs ",chem_names[(j-1)],sep="")
               )
          abline(fit <- lm(wax ~ gene_DE_1[,1]), col='black')
          R2<-format(summary(fit)$adj.r.squared, digits=3)
          legend("top",bty="n", legend=bquote(italic(R)^{2}~"="~.(R2)))
          legend("right",legend=paste("section",2:8,sep=" "),pch=19,col=c(1:6,"darkorange"),cex = 0.6)
          
         }

     }
   }
 }
 dev.off()
 
 cor_res<-as.data.frame(cbind(gene_name,wax_name,Chr,Position,r,r2,pvalue))
 cor_res$pvalue<-as.numeric(as.character(cor_res$pvalue))
 cor_res$gene_name<-as.character(cor_res$gene_name)
 all_genes$v3_gene_model<-as.character(all_genes$v3_gene_model)
# setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\RNAseq_ReadDepth")
# write.table(cor_res,"AllCand_vs_Wax_cor_12212018.txt",col.names=T,row.names=T,sep="\t",quote=F)

###### Heatmap of correlations
rsquared<-cor_res
rsquared$r2<-as.numeric(as.character(rsquared$r2))
mean(rsquared$r2,na.rm=T)
sd(rsquared$r2,na.rm=T)


rsquared.1<-matrix(rsquared$r,byrow=T,ncol=41) # not all genes has non-zero expression
colnames(rsquared.1)<-rsquared$wax_name[1:41]
genelist<-unique(rsquared$gene_name)
rownames(rsquared.1)<-rsquared$gene_name[match(genelist,rsquared$gene_name)]

pmatrix<-matrix(rsquared$pvalue,byrow=T,ncol=41)

starmatrix<-matrix(nrow=nrow(rsquared.1),ncol=ncol(rsquared.1))
class(rsquared.1) <- "numeric"
starmatrix[which(pmatrix<p_adj)]<-"*"
#colors = seq(from=min(range(rsquared.1)), to=max(range(rsquared.1)), length.out=30)

collab<-chem_names
rowlab<-rownames(rsquared.1)
for (i in 1:length(rowlab)){
  if (any("*" %in% starmatrix[i,])){
    rowlab[i]<-rowlab[i]
  } else {rowlab[i]<-NA}
}

library(RColorBrewer)
library(gplots)
#pdf("corr_74gene_41wax.pdf",width=12, height=10)
pdf("corr_100gene_41wax_simpleM_02252019.pdf",width=13, height=12)
heatmap.2(rsquared.1,
          col=colorRampPalette(brewer.pal(9, "RdBu"))(50),
          #col=redblue(75),
          cellnote = starmatrix,key=TRUE,trace="none",
          dendrogram="none",cexRow=0.9,cexCol=0.8,margins=c(14,10),
          Rowv=TRUE,Colv=TRUE,
          labRow=rowlab,labCol=collab,
          keysize=0.8, key.par = list(cex=0.5)
)
dev.off()

#####
qqnorm(cor_res$pvalue, pch = 1, frame = FALSE)



