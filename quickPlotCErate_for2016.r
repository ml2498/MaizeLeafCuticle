setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\16FromRawData\\aggregatedSampleData_old")
setwd("/Users/Meng/Google Drive/MLC_AZ_2017/Pheno/16FromRawData/aggregatedSampleData_old")
raw<-read.csv("MLC_AZ16_AggregatedSampleData_Filtered.csv",header=T,sep="\t")
raw<-read.csv("MLC_SD16_AggregatedSampleData_Filtered.csv",header=T,sep="\t")

#i=1
pdf("check_SD16_rawCErate_plot_11022018.pdf",width=5,height=5)
for (i in 1:dim(raw)[1]){
  #for (i in 1:20){
  x<-t(raw[i,c(5,9,13,17,21,25)])
  x1<-x[!is.na(x)]
  y<-t(raw[i,c(6,10,14,18,22,26)])
  y1<-y[!is.na(x)]
  y1<-(y1-y1[1])/60
  #x<-x-min(x)
  plot(y1,x1,col="blue",cex=1.5,xlab="Time (min)",ylab="Weight (g)",main=paste(as.character(raw[i,1]),as.character(raw[i,2]),sep="_"))
  abline(lm(x1~y1),col="red")
  #slope[i,2]<-(-3600)*coef(lm(y~x))["x"] # (-3600) change the unit from g/s to g/hour
  #slope[i,1]<-paste(as.character(raw[i,1]),as.character(raw[i,2]),sep="_")
}
dev.off()


######### calculate raw CE rate with filtered wet weight #########
raw<-read.csv("MLC_AZ16_AggregatedSampleData_Filtered_Meng.csv",header=T,sep=",")
raw<-read.csv("MLC_SD16_AggregatedSampleData_Filtered_Meng.csv",header=T,sep=",")

slope<-data.frame(Geno=as.character(),Leaf=character(), 
                  RawRate=numeric(),DW=numeric(),stringsAsFactors=FALSE)

for (i in 1:dim(raw)[1]){
  x<-t(raw[i,seq(5,28,4)])
  x1<-x[!is.na(x)]
  y<-t(raw[i,seq(6,28,4)])
  y1<-y[!is.na(x)]
  slope[i,3]<-(-1)*3600*coef(lm(x1~y1))["y1"] # the unit from g/s
  slope[i,1]<-as.character(raw[i,1])
  slope[i,2]<-raw[i,2]
  slope[i,4]<-raw[i,3]
}
slope$CEadj<-slope[,3]/slope[,4]

###########################
CEmean<-aggregate(RawRate~Geno,slope,mean)
CEadj<-aggregate(CEadj~Geno,slope,mean)
CEmean<-cbind(CEmean,CEadj[,2])

colnames(CEmean)<-c("Geno","CEmean","CEadj")
hist(CEmean[,2])
hist(CEmean[,3])
################################
FieldDesign<-read.table("C:\\Meng\\MaizeLeafCuticle_by02012018\\MLC_AZ_2017\\Pheno\\CE_BLUP\\re-calcu_16CErate\\Maricopa_Field_Design_with_Checks_20160422_nk_FINAL.csv",header=T,sep=",")
FieldDesign<-read.table("C:\\Meng\\MaizeLeafCuticle_by02012018\\MLC_AZ_2017\\Pheno\\CE_BLUP\\re-calcu_16CErate\\UCSD16_Field_Design_with_Checks_05042016.txt",header=T,sep="\t")

colnames(FieldDesign)[7]<-"Geno"
forCEBLUP<-merge(CEmean,FieldDesign,by="Geno",all=F)

##

map<-read.table("C:\\Meng\\MaizeLeafCuticle_by02012018\\MLC_AZ_2017\\Pheno\\CE_BLUP\\re-calcu_16CErate\\TaxaMap_sasname_16FieldBook.txt",header=T,sep="\t")

colnames(map)[2]<-"name" # AZ16
colnames(map)[3]<-"name" # SD16

forCEBLUP<-merge(map,forCEBLUP,by="name",all=F)


for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,2]=="MO17") {forCEBLUP$CHECK[i]<-1}
  else if (forCEBLUP[i,2]=="B73") {forCEBLUP$CHECK[i]<-2} 
  else forCEBLUP$CHECK[i]<-99
}
for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,2]=="MO17" | forCEBLUP[i,2]=="B73") {forCEBLUP$IS_EXEXPERIMENTAL[i]<-0} else forCEBLUP$IS_EXEXPERIMENTAL[i]<-1
}
#forCEBLUP$IS_EXEXPERIMENTAL[which(as.character(forCEBLUP[,17]=="1"))] <- 0
#forCEBLUP$IS_EXEXPERIMENTAL[which(as.character(forCEBLUP[,17]=="99"))] <- 1
forCEBLUP1<-forCEBLUP[,c(2,4,6,7,9,10,15,16)]#AZ
forCEBLUP1<-forCEBLUP[,c(2,4,6,7,9,10,14,15)]#SD

forCEBLUP1$Year<-"16"
forCEBLUP1$Loc<-"AZ"
forCEBLUP1$ENV<-1

forCEBLUP1$Loc<-"SD"
forCEBLUP1$ENV<-2

colnames(forCEBLUP1)[5:6]<-c("COL","BLOCK")
forCEBLUP1<-forCEBLUP1[,c(1,2,9,10,11,5:8,3,4)]#re-order columns
write.table(forCEBLUP1,"SD16_CERawAdj_BLUPinput2018.txt",col.names=T,row.names=F,sep="\t",quote=FALSE)
