########## pre-process the gff file on server
cd MaizeLeafCuticle/Meng_GAPIT/GeneAnnotation/
awk '$3 == "gene"' Zea_mays.AGPv3.31.gff3 > Zea_mays.AGPv3.31.short.gff3
#############################################
# download the short form to 
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\v3_gff_annotation")
setwd("/Users/menglin/Google Drive/MLC_AZ_2017/gene_study/v3_gff_annotation")

v3_gene_gff<-read.delim("Zea_mays.AGPv3.31.short.gff3", header=FALSE, stringsAsFactors=FALSE)

#### parse out gene names #########
v3_gene_gff$V9<-as.character(v3_gene_gff$V9)
gene_name<-strsplit(v3_gene_gff[,9],";") #the v4 gene names are in the last anotation column and need be pick out
gene_name<-sapply(gene_name, "[", 1)  #take the 2nd element of each component in a list
gene_name<- substr(gene_name,9,nchar(gene_name))
v3_gene_gff<-cbind(v3_gene_gff[,1:5],gene_name)
colnames(v3_gene_gff)<-c("chr","source","type","start","end","v3_gene_model")

######################################

## gapit hits
gapit_hits<-read.table("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\sum_GP_SNP_forCandidate_wrapper.txt",header=T,sep="\t")
#gapit_hits<-read.table("/Users/Meng/Google Drive/MLC_AZ_2017/gene_study/sum_GP_SNP_forCandidate_1.txt",header=T,sep="\t")
gapit_hits<-read.table("/Users/menglin/Google Drive/MLC_AZ_2017/gene_study/sum_GP_SNP_forCandidate_wrapper.txt",header=T,sep="\t")

## gene annotation v3 
gene_anno<-read.delim("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\v3_gff_annotation\\ZmB73_5a_gene_descriptors.txt",header=T,sep="\t")
gene_anno<-read.delim("/Users/menglin/Google Drive/MLC_AZ_2017/gene_study/v3_gff_annotation/ZmB73_5a_gene_descriptors.txt",header=T,sep="\t")

gene_anno<-gene_anno[,c(1,4)]

## SNP positions on v3
#v3_pos<-read.table("/Users/Meng/Google Drive/Genotype/liftover2018/v2_v3/GBSvmatch_v2_v3_pos_20181025.txt",header=T,sep="\t",stringsAsFactors=F)
#v3_pos<-read.table("C:\\Users\\ml2498\\Google Drive\\Genotype\\liftover2018\\v2_v3\\GBSvmatch_v2_v3_pos_20181025.txt",header=T,sep="\t",stringsAsFactors=F)

# incorperate gapit hits with v3 coordinates
keeplist<-unique(gapit_hits$SNP)
gapit_hits_uniq<-gapit_hits[match(keeplist,gapit_hits$SNP),]
#colnames(v3_pos)[1]<-"SNP"
#gapit_hits_uniq<-merge(v3_pos,gapit_hits_uniq,by="SNP",all=F)
#gapit_hits_uniq<-gapit_hits_uniq[,-(2:3)]

## find genes names in v3 and v4 within +/- 250 kb window from gapit hits
gapit_hits_uniq$left250<-gapit_hits_uniq$Position-250000
gapit_hits_uniq$right250<-gapit_hits_uniq$Position+250000

all_genes<-matrix(ncol=15,nrow=0)
for (i in 1:nrow(gapit_hits_uniq)){
  one_hit<-gapit_hits_uniq[i,]
  sub_genes<-v3_gene_gff[which(v3_gene_gff$chr==one_hit$Chromosome & v3_gene_gff$start<one_hit$right250 & v3_gene_gff$end>one_hit$left250),c(1,4:6)]
  nm_subgenes<-nrow(sub_genes)
  hit_info<-one_hit[rep(1, each=nm_subgenes),]
  sub_genes<-cbind(hit_info,sub_genes)
  all_genes<-rbind(all_genes,sub_genes)
}

## remove duplicated genes
length(unique(all_genes$v3_gene_model))
keeplist2<-unique(all_genes$v3_gene_model)
all_genes<-all_genes[match(keeplist2,all_genes$v3_gene_model),]
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\wrapper_cand")
write.table(all_genes,"all_100_cand_genes_250kb.txt",col.names=T,row.names=F,quote=F,sep="\t")

setwd("/Users/Meng/Google Drive/MLC_AZ_2017/gene_study/wrapper_cand")
write.table(all_genes,"all_100_cand_genes_250kb.txt",col.names=T,row.names=F,quote=F,sep="\t")

########################################
###################################################################################
## find gene functions of the above genes (v3 annotation)
#
colnames(gene_anno)[1]<-"v3_gene_model"
all_genes$v3_gene_model<-as.character(all_genes$v3_gene_model)

all_genes_v3ANNO<-merge(all_genes,gene_anno,by="v3_gene_model",all.y=F,all.x=T)

setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\wrapper_cand")
write.table(all_genes_v3ANNO,"GAPIT_cand_250K_v3anno_1218.txt",row.names=F,col.names=T,sep="\t",quote=F)

setwd("/Users/Meng/Google Drive/MLC_AZ_2017/gene_study/wrapper_cand")
write.table(all_genes_v3TrANNO,"GAPIT_cand_250K_v3anno_1218.txt",row.names=F,col.names=T,sep="\t",quote=F)

###################################################################################

## find gene functions of the above genes: Isabel orthologs
source.all<-read.table("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\from_Isabel\\LipGene_function_maize_AT_3files.txt",header=T,sep="\t")
source.all<-read.table("/Users/menglin/Google Drive/MLC_AZ_2017/gene_study/from_Isabel/LipGene_function_maize_AT_3files.txt",header=T,sep="\t")
colnames(source.all)[2]<-"v3_gene_model"
all_genes_v3ORTH<-merge(all_genes,source.all,by="v3_gene_model",all=F)
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\wrapper_cand")
write.table(all_genes_v3ORTH,"GAPIT_cand_250K_wrapper_v3ORTH_1218.txt",row.names=F,col.names=T,sep="\t",quote=F)

setwd("/Users/Meng/Google Drive/MLC_AZ_2017/gene_study/wrapper_cand")
write.table(all_genes_v3ORTH,"GAPIT_cand_250K_wrapper_v3ORTH_1218.txt",row.names=F,col.names=T,sep="\t",quote=F)

##################################################################################3
## find gene functions of the above genes: Pengfei's v3 gene annotation
gene_anno_pq<-read.delim("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\combined_maize_orth_annotation_info_3.30.txt",header=T,sep="\t",comment.char="")
gene_anno_pq<-read.delim("/Users/menglin/Google Drive/MLC_AZ_2017/fromPengfei/combined_maize_orth_annotation_info_3.30.txt",header=T,sep="\t",comment.char="")

#R_A_annotation<-Anno_database_maize_at_rice # This annotation only has v4 maize gene ID that correspond to AT and Rice

colnames(gene_anno_pq)[1]<-"v3_gene_model"
all_genes_v3ANNO2<-merge(all_genes_v3ANNO,gene_anno_pq,by="v3_gene_model",all.x=T,all.y=F)

###### add a column to indicate if the gene is in the lipid biosynthesis pathway in at
all_genes_v3ANNO2$LipidPathway<-NA
for (i in 1:nrow(all_genes_v3ANNO2)){
  if (all_genes_v3ANNO2$v3_gene_model[i] %in% all_genes_v3ORTH$v3_gene_model){
    all_genes_v3ANNO2$LipidPathway[i]<-"Y"
  }
}

####### Distance ############
all_genes_v3ANNO2$Dist<-NA
for (i in 1:nrow(all_genes_v3ANNO2)){
  if (all_genes_v3ANNO2$Position[i]>all_genes_v3ANNO2$start[i]&all_genes_v3ANNO2$Position[i]<all_genes_v3ANNO2$end[i]){
    all_genes_v3ANNO2$Dist[i]<-0
  }else {
    all_genes_v3ANNO2$Dist[i]<-min(c(abs(all_genes_v3ANNO2$Position[i]-all_genes_v3ANNO2$start[i]),abs(all_genes_v3ANNO2$Position[i]-all_genes_v3ANNO2$end[i])))
  }
  
}

###### r2 ################
install.packages("lubridate")
library(lubridate)

geno.all<-read.table("C:\\Users\\ml2498\\Google Drive\\Genotype\\merge4\\Individual_Set\\GBS_454_286K_CR06_AGPv3.hmp.txt",header=T,sep="\t")
geno.all<-read.table("/Users/menglin/Google Drive/Genotype/merge4/Individual_Set/GBS_454_286K_CR06_AGPv3.hmp.txt",header=T,sep="\t")

numericGeno<-function(geno){
  geno<-t(geno)
  colnames(geno)<-geno[1,];geno<-as.matrix(geno[-1,])
  
  alleles<-c("A","T","G","C")
  
  for (j in 1:ncol(geno)){
    
    freq<-table(geno[,j])
    freq<-freq[which(names(freq) %in% alleles)]
    
    freq<-freq[order(freq,decreasing=T)]
    major<-names(freq[1])
    minor<-names(freq[2])
    geno[,j]<-as.character(geno[,j])
    geno[which(geno[,j]==major),j]<-0
    geno[which(geno[,j]==minor),j]<-2
    geno[which(geno[,j]!=0&geno[,j]!=2),j]<-NA
  }
  return(geno)
}

both_end<-c("start","end")
all_genes_v3ANNO2$LD_r2_max<-NA
all_genes_v3ANNO2$LD_r2_median<-NA
all_genes_v3ANNO2$LD_r2_mean<-NA
all_genes_v3ANNO2$SNPs_in_gene<-NA

for (i in 1:nrow(all_genes_v3ANNO2)){
  SNP<-as.character(all_genes_v3ANNO2$SNP[i])
  snp_chr<-all_genes_v3ANNO2$Chromosome[i]
  snp_pos<-all_genes_v3ANNO2$Position[i]
  
  gene_start<-all_genes_v3ANNO2$start[i]
  gene_end<-all_genes_v3ANNO2$end[i]
  
  Interval<-seq(from=min(c(gene_start,gene_end)),to=max(c(gene_start,gene_end)),by=1)
  snps_in_gene<-geno.all[which(geno.all$chrom==snp_chr&geno.all$pos %in% Interval),-(2:11)]
  
  if (all_genes_v3ANNO2$Dist[i]==0){
    all_genes_v3ANNO2$LD_r2_max[i]<-1
    all_genes_v3ANNO2$LD_r2_median[i]<-1
    all_genes_v3ANNO2$LD_r2_mean[i]<-1
  }else{
    if (nrow(snps_in_gene)>0){
      geno<-snps_in_gene
      all_genes_v3ANNO2$SNPs_in_gene[i]<-"Y"
    } else {
      pick_end<-both_end[which.min(c(abs(snp_pos-gene_start),abs(snp_pos-gene_end)))]
      Interval2<-seq(from=min(c(all_genes_v3ANNO2[i,pick_end],snp_pos)),to=max(c(all_genes_v3ANNO2[i,pick_end],snp_pos)),by=1)
      geno.pre<-geno.all[which(geno.all$chrom==snp_chr,geno.all$pos %in% Interval2),-(2:11)]
      if (snp_pos>gene_start){
        geno<-geno.pre[1:10,]
      }else {
        geno<-tail(geno.pre,10)
      }
    }
  
  geno<-numericGeno(geno)
  ## geno.1 is the genotype of top SNP
  geno.1<-geno.all[c(which(geno.all$rs==SNP)),-(2:11)]
  geno.1<-numericGeno(geno.1)

  
  LD<-matrix(nrow=ncol(geno),ncol=2)
  for (j in 1:ncol(geno)){
    
    geno.2<-geno[,j] 
    geno.11<-geno.1[which(!is.na(geno.1)&!is.na(geno.2))]
    geno.22<-geno.2[which(!is.na(geno.1)&!is.na(geno.2))]
    #temp<-as.data.frame(cbind(as.numeric(geno.11),as.numeric(geno.22)))
    #temp$comb<-paste(temp[,1],temp[,2],sep="_")
    
    len<-length(geno.11)
    p1<-length(which(geno.11==0))/len;p2<-1-p1
    q1<-length(which(geno.22==0))/len;q2<-1-q1
    p1q1<-length(which(geno.11==0&geno.22==0))/len # 0 is major allele, so there has to be individuals have 00 genotype
    r2<-round((p1q1-p1*q1)^2/(p1*q1*p2*q2),3)
    
    #LD[j,1]<-colnames(geno)[j]
    LD[j,2]<-r2
  }
  all_genes_v3ANNO2$LD_r2_max[i]<-round(max(LD[,2]),3)
  all_genes_v3ANNO2$LD_r2_median[i]<-round(median(LD[,2]),3)
  all_genes_v3ANNO2$LD_r2_mean[i]<-round(mean(LD[,2]),3)
  }
  
}


##########################

setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\wrapper_cand")
setwd("/Users/menglin/Google Drive/MLC_AZ_2017/gene_study/wrapper_cand")
#write.table(all_genes_v3GeneANNO,"FCPU_cand_250K_v3AtRice_1029.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(all_genes_v3ANNO2,"GAPIT_cand_250K_wrapper_v3ZeaAtRice_1218.txt",row.names=F,col.names=T,sep="\t",quote=F)
############################################################################

## parse out the genes have not been searched on NCBI
ToBeBlast<-all_genes_v3ANNO2[which(all_genes_v3ANNO2$Chromosome %in% c(1,10)&all_genes_v3ANNO2$Position %in% c(194700493,19364955,144113824)),c(1,10:12)]
keeplist3<-unique(ToBeBlast[,1])
ToBeBlast<-ToBeBlast[match(keeplist3,ToBeBlast[,1]),]
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\gene_study\\BLAST_NCBI")
write.table(ToBeBlast,"genes_tobeBlast_set2_12192018.txt",col.names=T,row.names=F,sep="\t",quote=F)



keeplist3<-unique(all_genes_v3GeneANNO$v3_gene_model)
all_genes_v3GeneANNO_uniq<-all_genes_v3GeneANNO[match(keeplist3,all_genes_v3GeneANNO$v3_gene_model),]
write.table(all_genes_v3GeneANNO_uniq,"GAPIT_cand_250K_wrapper_v3AtRice_1124_uniq.txt",row.names=F,col.names=T,sep="\t",quote=F)

###############################
########################################

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


cuticle<-read.csv("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\RNAseq_ReadDepth\\Cuticle_PQ_new.csv",
                  header=T,stringsAsFactors=FALSE)
cuticle<-read.csv("/Users/Meng/Google Drive/MLC_AZ_2017/fromPengfei/RNAseq_ReadDepth/Cuticle_PQ_new.csv",
                  header=T,stringsAsFactors=FALSE)
cuticle<-cuticle[order(cuticle$Sample),]

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
library(qvalue)
POS<-unique(all_genes$Position)


nm_wax<-ncol(cuticle)
#rsquared<-data.frame(name=character(),r2=numeric(),stringsAsFactors=FALSE)
gene_name<-vector()
wax_name<-vector()
r2<-vector()
r<-vector()
pvalue<-vector()
Chr=vector()
Position=vector()
FDR=vector()


for (pos in POS){
  all_genes_pos<-all_genes[which(all_genes$Position==pos),]
  gene_DE<-gene_DE.all[which(rownames(gene_DE.all) %in% all_genes_pos$v3_gene_model),] ## testing 67 genes; 110 genes in wrapper's BLUP
  nm_gene<-nrow(gene_DE)
  p.hit<-vector()
  chr=unique(all_genes_pos$Chromosome)
  #rsquared<-matrix(nrow=nm_gene,ncol=(nm_wax-1))
  #colnames(rsquared)<-colnames(cuticle)[-1]
  #rownames(rsquared)<-rownames(gene_DE)
  for (i in 1:nm_gene){
    gene_DE_1<-t(gene_DE[i,])
    if (any(gene_DE_1!=0)){
      for (j in 2:nm_wax){
        wax<-cuticle[,j]
        p.1<-cor.test(wax,gene_DE_1[,1])$p.value
        pvalue<-c(pvalue,p.1)
        p.hit<-c(p.hit,p.1)
        
        r.1<-cor(wax,gene_DE_1[,1],use="complete.obs")
        r<-c(r,r.1)
        r2.1<-r.1^2
        r2<-c(r2,r2.1)
        gene_name<-c(gene_name,rownames(gene_DE)[i])
        wax_name<-c(wax_name,chem_names[j-1])
        Chr=c(Chr,chr)
        Position=c(Position,pos)
        #rsquared_1<-c(name,r2)
        #rsquared<-rbind(rsquared,rsquared_1)
      }
    }
  }
  qobj <- qvalue(p = p.hit)
  qvalues <- qobj$qvalues
  FDR<-c(FDR,qvalues)
}
cor_res<-as.data.frame(cbind(gene_name,wax_name,Chr,Position,r,r2,pvalue,FDR))
cor_res$pvalue<-as.numeric(as.character(cor_res$pvalue))
cor_res$FDR<-as.numeric(as.character(cor_res$FDR))
cor_res$gene_name<-as.character(cor_res$gene_name)
all_genes$v3_gene_model<-as.character(all_genes$v3_gene_model)
# cor_res$chr<-NA
# cor_res$SNP<-NA
# 
# for (i in 1:nrow(cor_res)){
#   cor_res$chr[i]<-all_genes$Chromosome[which(all_genes$v3_gene_model==cor_res$gene_name[i])]
#   cor_res$SNP[i]<-as.character(all_genes$SNP[which(all_genes$v3_gene_model==cor_res$gene_name[i])])
# }
#cor_res$SNP<-as.factor(cor_res$SNP)
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\RNAseq_ReadDepth")
write.table(cor_res,"AllCand_vs_Wax_cor_12212018.txt",col.names=T,row.names=T,sep="\t",quote=F)


#setwd("/Users/Meng/Google Drive/MLC_AZ_2017/fromPengfei/RNAseq_ReadDepth")
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\RNAseq_ReadDepth")
#pdf("gene_DE_vs_cuticle_comp_10312018.pdf",width=4.5,height=4.5)
#pdf("gene_DE_vs_cuticle_comp_wrapper_11262018.pdf",width=4.5,height=4.5)
pdf("gene_DE_vs_cuticle_comp_wrapper_fdr_12212018.pdf",width=4.5,height=4.5)

BelowFDR<-cor_res[which(cor_res$FDR<0.05),]
nm_entry<-nrow(BelowFDR)
for (i in 1:nm_entry){
  chem<-as.character(BelowFDR$wax_name[i])
  gene<-BelowFDR$gene_name[i]
  chem_index<-which(chem_names==chem)+1
  gene_index<-which(rownames(gene_DE.all)==gene)
  
  gene_DE_1<-t(gene_DE.all[gene_index,])
  wax<-cuticle[,chem_index]

  plot(wax~gene_DE_1[,1],pch=19,col=rep(c(1:6,"darkorange"),3), xlab="Normalized gene expression level",
       ylab=chem, main=paste(gene," vs ",chem,sep="")
        )
  abline(fit <- lm(wax ~ gene_DE_1[,1]), col='black')
  R2<-format(summary(fit)$adj.r.squared, digits=3)
  legend("top",bty="n", legend=bquote(italic(R)^{2}~"="~.(R2)))
  legend("right",legend=paste("section",2:8,sep=" "),pch=19,col=c(1:6,"darkorange"),cex = 0.6)
    
}
dev.off()



# for (i in 1:nm_gene){
#   gene_DE_1<-t(gene_DE[i,])
#   if (any(gene_DE_1!=0)){
#     for (j in 2:nm_wax){
#       wax<-cuticle[,j]
#       r.1<-cor(wax,gene_DE_1[,1],use="complete.obs")
#       r<-c(r,r.1)
#       r2.1<-r.1^2
#       r2<-c(r2,r2.1)
#       gene_name<-c(gene_name,rownames(gene_DE)[i])
#       wax_name<-c(wax_name,colnames(cuticle)[j]) 
#       # if (r2.1>0.7){
#       #   plot(wax~gene_DE_1[,1],pch=19,col=rep(c(1:6,"darkorange"),3), xlab="Normalized gene expression level", 
#       #        ylab=chem_names[j-1],
#       #        main=paste(rownames(gene_DE)[i]," vs ",chem_names[(j-1)],sep="")
#       #        )
#       #   abline(fit <- lm(wax ~ gene_DE_1[,1]), col='black')
#       #   R2<-format(summary(fit)$adj.r.squared, digits=3)
#       #   legend("top",bty="n", legend=bquote(italic(R)^{2}~"="~.(R2)))
#       #   legend("right",legend=paste("section",2:8,sep=" "),pch=19,col=c(1:6,"darkorange"),cex = 0.6)
#       #   }
#       #rsquared_1<-c(name,r2)
#       #rsquared<-rbind(rsquared,rsquared_1)
#     }
#   }
# }
# dev.off()
# rsquared<-as.data.frame(cbind(gene_name,wax_name,r,r2))
# write.table(rsquared,"R2_cor_110gene42wax.txt",col.names=T,row.names=F,quote=F,sep="\t")


###### Heatmap of correlations
rsquared<-cor_res
rsquared$r2<-as.numeric(as.character(rsquared$r2))
mean(rsquared$r2,na.rm=T)
sd(rsquared$r2,na.rm=T)


rsquared.1<-matrix(rsquared$r,byrow=T,ncol=41) # not all genes has non-zero expression
colnames(rsquared.1)<-rsquared$wax_name[1:41]
genelist<-unique(rsquared$gene_name)
rownames(rsquared.1)<-rsquared$gene_name[match(genelist,rsquared$gene_name)]

FDRmatrix<-matrix(rsquared$FDR,byrow=T,ncol=41)

starmatrix<-matrix(nrow=nrow(rsquared.1),ncol=ncol(rsquared.1))
class(rsquared.1) <- "numeric"
starmatrix[which(FDRmatrix<0.05)]<-"*"
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
pdf("corr_100gene_41wax_12212018.pdf",width=12, height=10)
heatmap.2(rsquared.1,
          col=colorRampPalette(brewer.pal(9, "RdBu"))(50),
          #col=redblue(75),
          cellnote = starmatrix,key=TRUE,Rowv=TRUE,Colv=TRUE,trace="none",
          dendrogram="none",cexRow=0.5,cexCol=0.5,margins=c(9,10),
          labRow=rowlab,labCol=collab,
          keysize=0.8, key.par = list(cex=0.5)
)
dev.off()





###
rsquared<-read.table("/Users/Meng/Google Drive/MLC_AZ_2017/fromPengfei/RNAseq_ReadDepth/R2_cor_110gene42wax.txt",header=T,sep="\t")
rsquared$r2<-as.numeric(as.character(rsquared$r2))
out_std<-rsquared[which(rsquared$r2>0.7),]
length(unique(out_std$gene_name)) # 15 unique genes
write.table(out_std,"R2_cor_76genexwax.txt",col.names=T,row.names=F,quote=F,sep="\t")

## for a single graph
genename<-"GRMZM2G150166"
cuticlename<-"Hydroxy.factty.acid.16.0.16.OH"
cuticlename<-"Total.cutin"
cuticlename<-"Alkane.C25"

setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\fromPengfei\\RNAseq_ReadDepth")
setwd("/Users/Meng/Google Drive/MLC_AZ_2017/fromPengfei/RNAseq_ReadDepth")

pdf(paste(genename,"_vs_",cuticlename,"_10312018.pdf",sep=""),width=5,height=5)
par(mar=c(5.1,5.1,4.1,2.1))
gene_DE_1<-t(gene_DE[which(rownames(gene_DE)==genename),])
wax<-cuticle[,which(colnames(cuticle)==cuticlename)]
r<-round(cor(wax,gene_DE_1[,1],use="complete.obs"),3)
plot(wax~gene_DE_1[,1],pch=19,col=rep(c(1:6,"darkorange"),3), xlab="Normalized gene expression of CAP", 
     #ylab=expression(paste("Hydroxy factty acid 16:0 16-OH (ug/dm"^"2",")",sep="")),
     ylab=expression(paste("Total cutin (ug/dm"^"2",")",sep="")),
     #ylab=expression(paste("Alkane C25 (ug/dm"^"2",")",sep="")),
     #main=paste("CAP vs Hydroxy factty acid 16:0 16-OH",sep="")
     main=paste("CAP vs Total cutin",sep="")
     #main=paste("CAP vs Alkane C25",sep="")
)
abline(fit <- lm(wax ~ gene_DE_1[,1]), col='black')
R2<-format(summary(fit)$adj.r.squared, digits=3)
legend("top",bty="n", legend=bquote(italic(R)^{2}~"="~.(R2)))
legend("right",legend=paste("section",2:8,sep=" "),pch=19,col=c(1:6,"darkorange"),cex = 0.6)

dev.off()

