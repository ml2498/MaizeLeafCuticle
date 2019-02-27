library(RColorBrewer)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(compiler)
library(scatterplot3d)
library(EMMREML)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

## IF GAPIT resutls
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\GAPIT\\CE_v3_wrapper_11092018")
peak_snp=read.table('summary_CE_GAPIT_G_wrapper_forPlot_11132018.txt',header=T,sep="\t")

setwd("/Users/Meng/Google Drive/MLC_AZ_2017/GAPIT/CE_v3_wrapper_11092018")
peak_snp=read.table('summary_CE_GAPIT_G_wrapper_forPlot_11132018.txt',header=T,sep="\t")

sample_file=read.csv('GAPIT.MLM.ce_All_untr.GWAS.Results.csv')
softw<-"GAPIT"
## IF FCPU resutls
setwd("/Users/Meng/Google Drive/MLC_AZ_2017/FarmCPU/CErate_GBS_09242018")
peak_snp=read.table('summary_CE_FCPU_Gfinal_forPlot_10282018.txt',header=T,sep="\t")
sample_file=read.csv('FarmCPU.ce_ALL4_untr.GWAS.Results.csv')
softw<-"FCPU"



#get trait names
#traits=unique(peak_snp$exp)
traits=c("ce_AZ_untr","ce_SD_untr","ce_All_untr")

######################################
# for AZ, SD, All only, no Y16 or Y17
######################################
peak_snp=peak_snp[which(peak_snp$exp %in% traits),]

#traits_short=gsub("_untr","",as.character(traits))
#traits_short=c("MA","SD","Y16","Y17","AllEnv")
traits_short=c("MA","SD","AllEnv")
y.lim <- ceiling(max(-log10(peak_snp$P.value)))
# y.lim=15


GI.MP=sample_file[,c("Chromosome","Position","P.value")]
GI.MP=matrix(as.numeric(as.matrix(GI.MP) ) ,nrow(GI.MP),ncol(GI.MP))
#Remove all SNPs that do not have a choromosome, bp position and p value(NA)
GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
GI.MP <- GI.MP[!is.na(GI.MP[,3]),]
#Retain SNPs that have P values between 0 and 1 (not na etc)
GI.MP <- GI.MP[GI.MP[,3]>0,]
GI.MP <- GI.MP[GI.MP[,3]<=1,]
#Remove chr 0 and 99
GI.MP <- GI.MP[GI.MP[,1]!=0,]
numMarker=nrow(GI.MP)
#Replace P the -log10 of the P-values
GI.MP[,3] <-  -log10(GI.MP[,3])
chm.to.analyze <- unique(GI.MP[,1])
chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
numCHR= length(chm.to.analyze)
nchr=length(chm.to.analyze)
#
GI.MP <- GI.MP[order(GI.MP[,2]),]
GI.MP <- GI.MP[order(GI.MP[,1]),]
ticks=NULL
lastbase=0
lastbase_chr=c()
#change base position to accumulatives (ticks)
for (i in chm.to.analyze){
  index=(GI.MP[,1]==i)
  ticks <- c(ticks, lastbase+mean(GI.MP[index,2]))
  GI.MP[index,2]=GI.MP[index,2]+lastbase
  lastbase_chr=c(lastbase_chr,lastbase)
  lastbase=max(GI.MP[index,2])
}
x0 <- as.numeric(GI.MP[,2])
y0 <- as.numeric(GI.MP[,3])
z0 <- as.numeric(GI.MP[,1])
position=order(y0,decreasing = TRUE)
index0=GAPIT.Pruning(y0[position])
index=position[index0]


##format peak snps file
peak_snp$Position_cum=lastbase_chr[peak_snp$Chromosome]+peak_snp$Position
peak_snp$logp=-log10(peak_snp$P.value)

{
  x=x0[index]
  y=y0[index]
  z=z0[index]
  
  pdf(paste("GWAS.",softw,".Manhattan.Genomewise_3comb_5Genes.pdf" ,sep = ""), width = 12,height=5)
  #pdf(paste("GWAS.",softw,".Manhattan.Genomewise_all_0.1_nolegends.pdf" ,sep = ""), width = 12,height=5)
   
     #jpeg(paste("GWAS.GAPIT.Manhattan.Genomewise_all_0.1.jpeg" ,sep = ""), width = 1200,height=500)
  
  par(mar = c(5,5,5,1)) #original: par(mar = c(10,8,6,1))
  par(oma=c(0,0,0,0))
  plot(y~x,xlab="",ylab="",ylim=c(0,y.lim+1),xlim=c(0,lastbase), cex.lab=3,pch=16,type="n",axes=FALSE,
       #main=paste('CE rate (',softw,')',sep="")
       main=paste('CE rate')
       )
  
  
  #rectangles for chromosome
  for(a in 1:floor(nchr/2)){
    rect(xleft=lastbase_chr[2*a-1],xright=lastbase_chr[2*a],ybottom=0,ytop=y.lim+2,col=rgb(239,247,250,max=255),border=NA)
  }
  #add last rect
  rect(xleft=lastbase_chr[11],xright=lastbase,ybottom=0,ytop=y.lim+2,col=rgb(239,247,250,max=255),border=NA)
  #add horizontal line
  #abline(h=5,lwd=1.2,col='gainsboro')
  #abline(h=10,lwd=1.2,col='gainsboro')
  #abline(h=15,lwd=1.2,col='gainsboro')
  #abline(h=20,lwd=1.2,col='gainsboro')
  
  #Set axises
  axis(1, at=ticks,cex.axis=1,labels=chm.to.analyze,tick=F,las=1)
  axis(2, at=c(0,5,10,15,20),cex.axis=1,labels=c(0,5,10,15,20),tick=F,las=1)
  mtext(side = 1, text = "Chromosome", line = 3,cex=1.2)
  palette=brewer.pal(n =12, name = "Paired")[c(2,4,6,8,10)] #remove gray
  palette[3]='darkorchid1';palette[2]='gold1'
  #plot significant snps
  for (b in 1:nrow(peak_snp)){
    points(x=peak_snp$Position_cum[b],y=peak_snp$logp[b],pch=16,cex=1.5,col=palette[which(traits==peak_snp$exp[b])])
  }
  
  #plot non-significant snps
  non_sig=GI.MP[which(GI.MP[,3]<min(peak_snp$logp)),]
  points(non_sig[,3]~non_sig[,2],pch=16,col='gray')
  abline(h=min(peak_snp$logp),col="red",lty="dashed", lwd=0.5)
  
  mtext(side = 2, line = 2.5, text=expression(paste(-log[10](italic(p)))),cex=1)
  
  for (c in 1:length(lastbase_chr)){
    segments(x0=lastbase_chr[c],x1=lastbase_chr[c],y0=-0.5,y1=y.lim+2,col=rgb(239,247,250,max=255),border=NA)}
  
  legend('topright',legend=traits_short,col=palette,pch=16,horiz = T, cex=1,bty='n')
  
  if (softw=="GAPIT"){
  text(x=270417003+lastbase_chr[1],y=(-log10(4.401134e-07)+1),labels='CAP',cex=1) #mute this line if for FCPU
  text(x=270417003+lastbase_chr[1],y=(-log10(4.401134e-07)+1.5),labels='CER7',cex=1) #mute this line if for FCPU
  text(x=30226164+lastbase_chr[4],y=(-log10(2.577895e-08)+1),labels='Vps4 regulator',cex=1)
  text(x=195018498+lastbase_chr[4],y=(-log10(1.194420e-07)+1),labels='ABCG',cex=1)
  text(x=193359341+lastbase_chr[5],y=(-log10(1.988904e-07)+1),labels='Aquaporin',cex=1) #mute this line if for FCPU
  #text(x=194700493+lastbase_chr[1],y=(-log10(9.41e-07)+0.5),labels='LTPG1',cex=1) #mute this line if for FCPU
  #text(x=19364955+lastbase_chr[10],y=(-log10(6.39e-07)+1),labels='HIBCH1',cex=1) #mute this line if for FCPU
  #text(x=231664302+lastbase_chr[2],y=(-log10(1.446443e-06)+1),labels='GPAT',cex=1) #mute this line if for FCPU
  # #abline(v=270417003+lastbase_chr[1],col="blue")
  #abline(v=30226164+lastbase_chr[4],col="red")
  #abline(v=195028437+lastbase_chr[4],col="red")
  #abline(v=193359341+lastbase_chr[5],col="red")
  }
  # if (softw=="FCPU"){
  # ## for FCPU
  # text(x=30226164+lastbase_chr[4],y=(-log10(2.544788e-11)+1),labels='Vps4 regulator',cex=1)
  # text(x=195028437+lastbase_chr[4],y=(-log10(3.422983e-15)+1),labels='ABC G',cex=1)
  # text(x=193359341+lastbase_chr[5],y=(-log10(1.438346e-09)+1),labels='Aquaporin',cex=1) #mute this line if for FCPU
  # text(x=231664302+lastbase_chr[2],y=(-log10(4.567095e-16)+1),labels='GPAT',cex=1) #mute this line if for FCPU
  # }
   
  dev.off()
}
