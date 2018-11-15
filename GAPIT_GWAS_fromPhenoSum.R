######### Combine GBS and RNA-seq after imputed seperately ######
#geno_G<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set/GBS_454_286K_CR06_AGPv3.hmp.txt",
#                     header=T,sep="\t",comment.char="")
#geno_R<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set/RNA_fastPHASE_451_419K_CR06_maf05_AGPv3.hmp.txt",
#                     header=T,sep="\t",comment.char="")
info_G<-geno_G[,1:11];geno_G<-geno_G[,-(1:11)]
info_R<-geno_R[,1:11];geno_R<-geno_R[,-(1:11)]
info<-rbind(info_G,info_R)

nm_Gsnp<-nrow(geno_G)
nm_Rsnp<-nrow(geno_R)

common_taxa<-intersect(colnames(geno_G),colnames(geno_R))
G_taxa<-colnames(geno_G)[-which(colnames(geno_G) %in% common_taxa)]
R_taxa<-colnames(geno_R)[-which(colnames(geno_R) %in% common_taxa)]

common_genoG<-geno_G[,which(colnames(geno_G) %in% common_taxa)]
common_genoG<-common_genoG[,order(colnames(common_genoG))]
common_genoR<-geno_R[,which(colnames(geno_R) %in% common_taxa)]
common_genoR<-common_genoR[,order(colnames(common_genoR))]
uniq_G<-geno_G[,G_taxa]
uniq_R<-as.data.frame(geno_R[,R_taxa])

common_geno<-rbind(common_genoG,common_genoR)
common_geno<-as.matrix(common_geno)
uniq_geno<-matrix(nrow=nrow(info),ncol=(ncol(uniq_G)+ncol(uniq_R)))
uniq_geno[1:nrow(uniq_G),1:ncol(uniq_G)]<-as.matrix(uniq_G)
uniq_geno[(nrow(uniq_G)+1):nrow(uniq_geno),(ncol(uniq_G)+1):ncol(uniq_geno)]<-as.matrix(uniq_R)
colnames(uniq_geno)<-c(G_taxa,R_taxa)

for (i in 1:ncol(uniq_geno)){
  uniq_geno[is.na(uniq_geno[,i]),i]<-"N"
}

geno.all<-cbind(common_geno,uniq_geno)
geno.all<-cbind(info,geno.all)
geno.all<-geno.all[order(geno.all[,3],geno.all[,4],geno.all[,1]),]
geno.all$pos<-format(geno.all$pos, scientific = FALSE)

setwd("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set")
write.table(geno.all,"GBS_RNA_bind_aftImpute_AGPv3.hmp.txt",col.names=T,row.names=F,quote=F,sep="\t")

###### subset unique positions ###########
###### or it cannot open by TASSEL #######
###### will use this set with unique positions for MLMM ######
geno.all$pos<-trimws(geno.all$pos, which = c("both", "left", "right"))
geno.all$chr_pos<-paste(geno.all$chrom,geno.all$pos,sep="_")
uniq_pos<-geno.all$chr_pos[!duplicated(geno.all$chr_pos)]

geno.all.1<-geno.all[match(uniq_pos,geno.all$chr_pos),] #only keep the GBS SNPs with unique postions

geno.all.1$pos<-as.numeric(as.character(geno.all.1$pos))
geno.all.1<-geno.all.1[order(geno.all.1[,3],geno.all.1[,4]),]
geno.all.1<-geno.all.1[,-ncol(geno.all.1)]
geno.all.1$pos<-format(geno.all.1$pos, scientific = FALSE)
geno.all.1$pos<-trimws(geno.all.1$pos, which = c("both", "left", "right"))

setwd("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set")
write.table(geno.all.1[1:2,],"GBS_RNA_bind_aftImpute_AGPv3_uniqpos_test.hmp.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(geno.all.1,"GBS_RNA_bind_aftImpute_AGPv3_uniqpos.hmp.txt",col.names=T,row.names=F,quote=F,sep="\t")


## for checking genotype string lengths #####
for (i in 12:ncol(geno.all)){
#for (i in 1:14){
  geno<-as.character(geno.all[,i])
  abnormal<-which(nchar(geno)>1)
  print(abnormal)
}

for (i in c(1,14,15)){
  for (j in 1:100){
    curr_geno<-as.character(geno.all[j,i])
    if (nchar(curr_geno)>1){
      print(paste(i,j,sep=","))
    }
  }
}
#############################################

######## Only Server can handle merged data #####################
setwd("/workdir/ml2498/MaizeLeafCuticle/Meng_GAPIT")
#pheno.all<-read.table("CE_FT_alllines_04072018.txt",header=T,sep="\t")
pheno.all<-read.table("CE_FT_alllines_11092018_wrapper.txt",header=T,sep="\t")
pheno.all<-read.table("CE_FT_alllines_09242018.txt",header=T,sep="\t")
pheno.all<-read.table("CE_FT_alllines_incldGXE_07302018.txt",header=T,sep="\t")
pheno.all<-read.table("CE_FT_alllines_incldGXEint_08222018.txt",header=T,sep="\t")
pheno.all<-read.table("STRUCTURE_K=3.txt",header=T,sep="\t")
#pheno.all<-read.table("asreml_ar1xar1_blup_summary.txt",header=T,sep="\t")
#### remove sweet/pop corn########
swt_pop<-read.table("SwtPopList.txt")
pheno.all<-pheno.all[-which(pheno.all$MLC_STANDARD %in% swt_pop[,1]),]
##################################

#geno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/merge4/cmb_MLC_435_463K_AGPv4.hmp.txt",
#                     header=F,sep="\t",comment.char="")
#geno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set/GBS_454_258K_CR06_AGPv4.hmp.txt",
#                     header=F,sep="\t",comment.char="")
#geno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set/RNA_451_348K_CR06_AGPv4.hmp.txt",
#                     header=F,sep="\t",comment.char="")
#geno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set/RNA_fastPHASE_451_396K_CR07_AGPv4.hmp.txt",
#                    header=F,sep="\t",comment.char="")
#geno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set/RNA_Beagle_452_388K_CR07_AGPv4.hmp.txt",
#                     header=F,sep="\t",comment.char="")
geno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set/GBS_454_286K_CR06_AGPv3.hmp.txt",
                     header=F,sep="\t",comment.char="")
#geno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set/RNA_fastPHASE_451_419K_CR06_maf05_AGPv3.hmp.txt",
#                     header=F,sep="\t",comment.char="")
#geno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/GBSdata_WiDiv_Swtcorn/Individual_Set/GBS_RNA_bind_aftImpute_AGPv3.hmp.txt",
#                     header=F,sep="\t",comment.char="")

# imputed from Hapmap3
geno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/Beagle/Guillaume/MLC_GBSSNP459_877K_v3crossmap_chr10.hmp.txt",
                     header=F,sep="\t",comment.char="")


myG<-geno.all

trait<-c(2,16:20)# 09242018 CE rate, just AZ, SD and all4
trait<-c(1,16:20,7:11)
trait<-c(2,12:15,18:19) # 12:15 for single env, 18:19 for 16 combined and 17 combined
trait<-c(2,14:22)# AZ, SD combined and all4
trait<-c(2,18:19,22,20:21,14:17)
trait<-c(2,14:17)# individual exp.
trait<-c(2,3:6,8:11,13)
#trait<-c(1,26,28) # for OLS slope and Gibbs slope in FW regresion,without intercept
#trait<-c(1,30,32) # for OLS slope and Gibbs slope in FW regresion, with intercept
myY<-pheno.all[,trait]
#geno.all<-matrix(unlist(geno.all), nrow=dim(geno.all)[1], byrow=F)

#taxa<-intersect(t(as.character(geno.all[1,])),as.character(pheno.all[,2]))

#myG<-geno.all[,c(1:11,which(geno.all[1,] %in% taxa))]
#myY<-pheno.all[which(pheno.all$MLC_STANDARD %in% taxa),c(2,trait)]

#myY<-myY[order(myY[,1]),]
#myG.sub<-myG[,12:dim(myG)[2]]
#myG.sub<-myG.sub[,order(myG.sub[1,])]
#myG<-cbind(myG[,1:11],myG.sub)

#setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\GAPIT")
setwd("/workdir/ml2498/MaizeLeafCuticle/Meng_GAPIT")

#write.table(myG,"merged_temp.hmp.txt",col.names=F,row.names=F,quote=FALSE,sep='\t')
#myG<-read.table("merged_temp.hmp.txt",header=F,sep="\t")

#myKI<-read.table("centeredIBS_merge4_prnLD02_noSwtPop.txt")
#myKI<-read.table("centeredIBS_merge4_prnLD02.txt")
myKI<-read.table("centeredIBS_GBS_Final_CR06_AGPv4_LD02.txt")
myKI<-read.table("centeredIBS_RNA_Final_CR07_AGPv4_LD02.txt") # RNA

myKI<-read.table("centeredIBS_GBS_Final_CR06_AGPv4_LD02_noSWTPOP.txt")

source("http://www.bioconductor.org/biocLite.R") 
biocLite("multtest")
install.packages("gplots") 
install.packages("LDheatmap") 
install.packages("genetics")
install.packages("EMMREML") 
install.packages("scatterplot3d") #The downloaded link at: http://cran.r-project.org/package=scatterplot3d

library(multtest) 
library(gplots) 
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R 
library("scatterplot3d")
#source("http://zzlab.net/GAPIT/gapit_functions.txt")
#source("/Users/Meng/Google Drive/MLC_AZ_2017/GAPIT/gapit_functions.R")
source("gapit_functions.R")
source("http://zzlab.net/GAPIT/emma.txt")


#setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\GAPIT\\GBSrelax_All3_DWadj_untr_CEless1100")
setwd("/workdir/ml2498/MaizeLeafCuticle/Meng_GAPIT/CE_v3_final_09242018")
setwd("/workdir/ml2498/MaizeLeafCuticle/Meng_GAPIT/CE_v3_final_09242018/Imp_hmp_chr4_31733810")
setwd("/workdir/ml2498/MaizeLeafCuticle/Meng_GAPIT/CE_v3_final_09242018/HMP3_GWAS/all4/chr10")
setwd("/ml2498/MaizeLeafCuticle/Meng_GAPIT/CE_v3_wrapper_11092018")

# Merged_FTall3_K
nm_ind<-nrow(myY)

myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  #PCA.total=20,
  KI=myKI,
  #kinship.algorithm="Zhang",
  group.from=nm_ind, # 353 for (lsa adjusted ce;RNA), 359 for (lsa adjusted ce;GBS)
  group.to=nm_ind,
  group.by=1,
  Major.allele.zero=T,
  Model.selection=TRUE
) #automatically include best number of K groups

##### multiple genotype files
setwd("/workdir/ml2498/MaizeLeafCuticle/Meng_GAPIT/CE_v3_final_09242018/HMP3_GWAS/all4/chr10")

nm_ind<-nrow(myY)
myGAPIT <- GAPIT(
  Y=myY,
  file.G="MLC_GBSSNP459_877K_v3crossmap_chr",
  file.Ext.G="hmp.txt",
  file.from=1,
  file.to=10,
  file.total=10,
  file.path="/workdir/ml2498/MaizeLeafCuticle/Hapmap3/Beagle/Guillaume/",
  
  KI=myKI,
  #kinship.algorithm="Zhang",
  group.from=nm_ind, # 353 for (lsa adjusted ce;RNA), 359 for (lsa adjusted ce;GBS)
  group.to=nm_ind,
  group.by=1,
  Major.allele.zero=T,
  Model.selection=TRUE
  
)




## GLM
myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  #PCA.total=20,
  KI=myKI,
  #kinship.algorithm="Zhang",
  group.from=0, # 353 for (lsa adjusted ce;RNA), 359 for (lsa adjusted ce;GBS)
  group.to=0,
  #group.by=1,
  Major.allele.zero=T,
  #Model.selection=TRUE
) #automatically include best number of K groups

##
setwd("/home/meng/MaizeLeafCuticle/Meng_GAPIT/CEall4_201804_PC20")
myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=20,
  #kinship.algorithm="Zhang",
  group.from=0, # 353 for (lsa adjusted ce;RNA), 359 for (lsa adjusted ce;GBS)
  group.to=0,
  #group.by=1,
  Major.allele.zero=T
  #Model.selection=TRUE
) #automatically include best number of K groups

###############################################
# With covariates
###############################################
# ...
# myG<-read.table(...)
myCV<-t(myG[which(myG[,1]=="rna5_92354212"),-(1:11)])
myCV<-as.data.frame(as.numeric(as.factor(myCV)))
myCV<-cbind(as.character(myY[,1]),myCV)
colnames(myCV)<-c("Taxa","rna5_92354212")
setwd("/home/meng/MaizeLeafCuticle/Meng_GAPIT/Merged_all4_logCE_K_chr5cov")
nm_ind<-nrow(myY)

myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  CV=myCV,
  #kinship.algorithm="Zhang",
  group.from=nm_ind, # 353 for (lsa adjusted ce;RNA), 359 for (lsa adjusted ce;GBS)
  group.to=nm_ind,
  group.by=1,
  Major.allele.zero=T,
  Model.selection=TRUE
)
############################################################################

## Manhattan plot & QQ plot ############3
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\GAPIT\\CE_v3_final_09242018")
library(qqman)

Files<-c("16","17","ALL4","AZ","SD")
for (f in Files){
  file<-paste("GAPIT.MLM.ce_",f,"_untr.GWAS.Results.csv",sep="")
  gemma<-read.csv(file,header=T)
  
  if (length(gemma$P.value[which(gemma$FDR_Adjusted_P.values<0.1)])>0){
    
    threshold<-(max(gemma$P.value[which(gemma$FDR_Adjusted_P.values<0.1)])+min(gemma$P.value[which(gemma$FDR_Adjusted_P.values>0.1)]))/2
    
    pdf(paste("CE_",f,"_GBS_man.pdf",sep=""),width=8,height=4)
    manhattan(gemma, chr = "Chromosome", bp = "Position", p = "P.value", snp = "SNP",
              col = c("navy", "darkorange1"), chrlabs = NULL,
              suggestiveline = FALSE,genomewideline =-log10(threshold), main=NULL)
    dev.off()
  }
}
###############################################################
###############################################################
# HapMap3 GWAS
###############################################################
################################################################
chr=10
setwd("/workdir/ml2498/MaizeLeafCuticle/Meng_GAPIT")
#pheno.all<-read.table("CE_FT_alllines_04072018.txt",header=T,sep="\t")
pheno.all<-read.table("CE_FT_alllines_09242018.txt",header=T,sep="\t")

#geno.all<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/Beagle/Guillaume/MLC_GBSSNP459_877K_v3crossmap_chr",chr,".hmp.txt",sep=""),
#                     header=F,sep="\t",comment.char="")


myG<-geno.all

#trait<-c(2,12:20)# 09242018 CE rate, just AZ, SD and all4
trait<-c(2,12:19) # 12:15 for single env, 18:19 for 16 combined and 17 combined
myY<-pheno.all[,trait]
myKI<-read.table("centeredIBS_GBS_Final_CR06_AGPv4_LD02.txt")

# source("http://www.bioconductor.org/biocLite.R") 
# biocLite("multtest")
# install.packages("gplots") 
# install.packages("LDheatmap") 
# install.packages("genetics")
# install.packages("EMMREML") 
# install.packages("scatterplot3d") #The downloaded link at: http://cran.r-project.org/package=scatterplot3d

library(multtest) 
library(gplots) 
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R 
library("scatterplot3d")
#source("http://zzlab.net/GAPIT/gapit_functions.txt")
#source("/Users/Meng/Google Drive/MLC_AZ_2017/GAPIT/gapit_functions.R")
source("gapit_functions.R")
source("http://zzlab.net/GAPIT/emma.txt")


#setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\GAPIT\\GBSrelax_All3_DWadj_untr_CEless1100")
setwd(paste("/workdir/ml2498/MaizeLeafCuticle/Meng_GAPIT/CE_v3_final_09242018/HMP3_GWAS/chr",chr,sep=""))


# Merged_FTall3_K
nm_ind<-nrow(myY)

myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  #PCA.total=20,
  KI=myKI,
  #kinship.algorithm="Zhang",
  group.from=nm_ind, # 353 for (lsa adjusted ce;RNA), 359 for (lsa adjusted ce;GBS)
  group.to=nm_ind,
  group.by=1,
  Major.allele.zero=T,
  Model.selection=TRUE
) #automatically include best number of K groups
