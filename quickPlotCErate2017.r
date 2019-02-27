raw<-read.table("~/Documents/Meng/MaizeLeafCuticle/MLC_Pilot/MLC_PT17_ManuallyFiltered_Sample_Datapoints_adjPT1.txt",header=T,sep="\t")
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\arhieve\\SD17")
raw<-read.csv("/Users/Meng/Google Drive/MLC_AZ_2017/Pheno/arhieve/MLC_SD17_Unfiltered_Sample_Datapoints.csv",header=T,sep=",")
raw<-read.csv("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\arhieve\\MLC_SD17_Unfiltered_Sample_Datapoints.csv",header=T,sep=",")

#setwd("C:\\Users\\ml2498\\Desktop")
#raw<-read.csv("C:\\Users\\ml2498\\Desktop\\MLC_SD17_Unfiltered_Sample_Datapoints.csv",header=T,sep="\t")

##*** correct one for SD17
raw<-read.csv("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\arhieve\\SD17\\MLC_SD17_Filtered_Sample_Datapoints(Final).csv",header=T,sep=",")


#raw <- raw[with(raw, order(as.Date(raw$Date,format="%m/%d/%Y"),Genotype, Leaf_Number)),]

write.table(raw, file ="MLC_AZ17_Unfiltered_Sample_Datapoints_sorted.txt",quote=FALSE, sep="\t", row.names = FALSE)
raw<-read.csv("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\arhieve\\MLC_SD17_Unfiltered_Sample_Datapoints_sorted_filtered.csv",header=T,sep=",")
raw<-read.csv("/Users/Meng/Google Drive/MLC_AZ_2017/Pheno/arhieve/SD17/MLC_SD17_Unfiltered_Sample_Datapoints.csv",header=T,sep="\t")

#i=1
pdf("check_SD17_rawCErate_plot_sorted.pdf")
for (i in 1:dim(raw)[1]){
#for (i in 1:20){
  x<-t(raw[i,c(10,12,14,16,18)])
  x1<-x[!is.na(x)]
  y<-t(raw[i,c(11,13,15,17,19)])
  y1<-y[!is.na(x)]
  #x<-x-min(x)
  plot(y1,x1,main=paste(as.character(raw[i,1]),as.character(raw[i,2]),sep="_"))
  #slope[i,2]<-(-3600)*coef(lm(y~x))["x"] # (-3600) change the unit from g/s to g/hour
  #slope[i,1]<-paste(as.character(raw[i,1]),as.character(raw[i,2]),sep="_")
}
dev.off()

library(ggplot2)
pdf("check_SD17_rawCErate_plot_sorted_filtered.pdf")
for (i in 1:nrow(raw)){
  x<-t(raw[i,c(10,12,14,16,18)])
  x1<-x[!is.na(x)]
  y<-t(raw[i,c(11,13,15,17,19)])
  y1<-as.numeric(y[!is.na(x)])
  data<-as.data.frame(cbind(x1,y1))
  p<-ggplot(data,aes(x=y1, y=x1)) +
    geom_point(shape=1, size=3) +    # Use hollow circles
    geom_smooth(method=lm,se=FALSE)+ # Add linear regression line
    ggtitle(paste(as.character(raw[i,1]),as.character(raw[i,2]),sep="_"))
  print(p)
}
dev.off()

######### calculate raw CE rate with filtered wet weight #########
raw<-read.csv("/Users/menglin/Google Drive/MLC_AZ_2017/Pheno/arhieve/MLC_AZ17_Filtered_Sample_Datapoints.csv",header=T,sep=",")
raw<-read.csv("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\arhieve\\MLC_AZ17_Filtered_Sample_Datapoints.csv",header=T,sep=",")


raw<-read.csv("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\arhieve\\SD17\\MLC_SD17_Filtered_Sample_Datapoints(Final).csv",header=T,sep=",")
raw<-read.csv("/Users/menglin/Google Drive/MLC_AZ_2017/Pheno/arhieve/SD17/MLC_SD17_Filtered_Sample_Datapoints(Final).csv",header=T,sep=",")



raw <- raw[with(raw, order(as.Date(raw$Date,format="%m/%d/%Y"),Genotype, Leaf_Number)),]
slope<-matrix(,nrow=dim(raw)[1],6)

for (i in 1:dim(raw)[1]){
  x<-t(raw[i,c(10,12,14,16,18)])
  x1<-x[!is.na(x)]
  y<-t(raw[i,c(11,13,15,17,19)])
  y1<-y[!is.na(x)]
  slope[i,2]<-(-1)*3600*coef(lm(x1~y1))["y1"] # the unit from g/h
  slope[i,1]<-paste(as.character(raw[i,1]),as.character(raw[i,2]),sep="_")
  slope[i,3]<-y1[1]
  slope[i,4]<-y1[length(y1)]
  slope[i,5]<-as.character(raw[i,8])
  slope[i,6]<-as.character(raw[i,9])
}
colnames(slope)<-c("Geno_Leaf","RawRate","StartTS","EndTS","Sensor","wwDate")
##############################################
LSA<-read.table("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\LeafImaging\\FinalLeafArea_MichaelM_07062017.txt",header=T,sep="\t")
LSA<-read.table("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\LeafImaging\\SD17_LSA_Meng_Susanne.txt",header=T,sep="\t")

WW_LSA<-merge(slope,LSA,by="Geno_Leaf",all=F) # change all=T to F when analyze
#WW_LSA_DW<-merge(WW_LSA,DW,by="Geno_Leaf",all=T) # Do not have DW data after 6/20
nm_leaf<-dim(WW_LSA)[1]
RawRate_col<-which(colnames(WW_LSA)=="RawRate")
area_col<-which(colnames(WW_LSA)=="area")

WW_LSA$CE_adj<-as.numeric(as.character(WW_LSA[,RawRate_col]))/as.numeric(as.character(WW_LSA[,area_col]))

write.table(WW_LSA,"C:\\Meng\\MaizeLeafCuticle\\MLC_AZ_2017\\Pheno\\arhieve\\for_outlier_temp.txt")

###############################################
split<-strsplit(sub("(.*)_(.+)$", "\\1 \\2", WW_LSA[,1]), ' ') #spilt the string by last "_
split.1<-matrix(unlist(split), nrow=length(split), byrow=T)
colnames(split.1)<-c("Geno","Leaf")
WW_LSA<-cbind(split.1,WW_LSA)


###############################################
# Corrected by leaf Dry Weight
##############################################
DW<-read.csv("/Users/menglin/Google Drive/MLC_AZ_2017/Pheno/arhieve/MLC_AZ17_Master_DryWeights_Processed.csv",header=T,sep="\t")

DW<-read.csv("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\arhieve\\SD17\\MLC_SD17_Master_DryWeights_Processed.csv",header=T,sep="\t")
DW<-read.csv("/Users/menglin/Google Drive/MLC_AZ_2017/Pheno/arhieve/SD17/MLC_SD17_Master_DryWeights_Processed.csv",header=T,sep="\t")

DW$Geno_Leaf<-paste(as.character(DW[,1]),as.character(DW[,2]),sep="_")

DW_LSA<-merge(slope,DW,by="Geno_Leaf",all=F) # change all=T to F when analyze
#WW_LSA_DW<-merge(WW_LSA,DW,by="Geno_Leaf",all=T) # Do not have DW data after 6/20
nm_leaf<-dim(DW_LSA)[1]
RawRate_col<-which(colnames(DW_LSA)=="RawRate")
dw_col<-which(colnames(DW_LSA)=="Dry_Weight")

DW_LSA$CE_adj<-as.numeric(as.character(DW_LSA[,RawRate_col]))/as.numeric(as.character(DW_LSA[,dw_col]))

#write.table(DW_LSA,"C:\\Meng\\MaizeLeafCuticle\\MLC_AZ_2017\\Pheno\\arhieve\\for_outlier_temp.txt")

###############################################
split<-strsplit(sub("(.*)_(.+)$", "\\1 \\2", DW_LSA[,1]), ' ') #spilt the string by last "_
split.1<-matrix(unlist(split), nrow=length(split), byrow=T)
colnames(split.1)<-c("Geno","Leaf")
DW_LSA<-cbind(split.1,DW_LSA)
WW_LSA<-DW_LSA
########### Incorperate RH and Field Design ######
########### Option 1: Use Room RH #########
colnames(WW_LSA)[1]<-"Genotype"

CEmean<-aggregate(as.numeric(as.character(CE_adj))~Genotype+wwDate,WW_LSA,mean)
CEraw<-aggregate(as.numeric(as.character(RawRate))~Genotype+wwDate,WW_LSA,mean)
pos_list<-match(CEmean[,1],WW_LSA$Genotype,nomatch=0)
WWDate<-WW_LSA$wwDate[pos_list]
STARTTS<-aggregate(as.numeric(as.character(StartTS))~Genotype+wwDate,WW_LSA,min)
ENDTS<-aggregate(as.numeric(as.character(EndTS))~Genotype+wwDate,WW_LSA,max)
CEmean<-cbind(CEraw[,c(1,3)],CEmean[,3],STARTTS[,3],ENDTS[,3],WWDate)

colnames(CEmean)<-c("Geno_mdf","rawCE","CEadj.mean","StartTS","EndTS","WWDate")
CEmean$CEadj.mean<-CEmean$CEadj.mean*1e4 #change the unit from g*h-1*mm2(-1) to g*h-1*m2(-1)

hist(as.numeric(as.character(CEmean[,2])),xlab="SD17_rawCErate",main="SD17_rawCErate_Distribution")
hist(as.numeric(as.character(CEmean[,3])),xlab="SD17_adjCErate",main="SD17_adjCErate_Distribution")

############# Common part (for LSA & DW correction) to incorperate Room RH ########################################
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\PB_database\\ProbeDB_FromNick_07062017")
setwd("/Users/menglin/Google Drive/MLC_AZ_2017/PB_database/ProbeDB_FromNick_07062017")
setwd("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\PB_database\\ProbeDB_FromSusanne_08112017")
setwd("/Users/menglin/Google Drive/MLC_AZ_2017/PB_database/ProbeDB_FromSusanne_08112017")

CR300<-read.table("HMP155A_reformatted_0623.txt",header=T,sep=",")
CR300<-read.table("HMP155A_reformatted_0811.txt",header=T,sep=",")
CR300$Timestamp<-as.numeric(as.POSIXlt(CR300$TIMESTAMP))#convert local time to unix time
probe0 = data.frame(matrix(vector(), dim(CR300)[1], 4,dimnames=list(c(), c("Timestamp", "ID","Temperature" ,"RH"))),stringsAsFactors=F)
probe0$Timestamp<-CR300$Timestamp
probe0$ID<-CR300$RECORD
probe0$Temperature<-CR300$AirTC_Avg
probe0$RH<-CR300$RH

probe6<-read.csv("probe6_Balerion.csv",sep="\t",header=T)
probe6_clean<-probe6[which(probe6$RH<110),]
probe0_clean<-probe0[which(probe0$RH<110),]
###for AZ17
time1<-as.numeric(as.POSIXlt("2017-06-18 16:00:00")) #EST 16:00 = AZ 13:00
time2<-as.numeric(as.POSIXlt("2017-06-18 17:00:00"))
nm_taxa<-dim(CEmean)[1]
for (i in 1:nm_taxa){
  t1<-CEmean$StartTS[i]
  t2<-CEmean$EndTS[i]
  if (CEmean$StartTS[i]>time1 & CEmean$StartTS[i]<time2) CEmean$RH[i]<-mean(probe0_clean[which(probe0_clean$Timestamp>t1&probe0_clean$Timestamp<t2),4])
  else CEmean$RH[i]<-mean(probe6_clean[which(probe6_clean$Timestamp>t1&probe6_clean$Timestamp<t2),4])
}
hist(CEmean$RH,main="Distribution of RH AZ17",xlab="Relative Humidity")
## for SD17
time1<-as.numeric(as.POSIXlt("2017-08-13 17:20:00")) #EST 16:00 = AZ 13:00
time2<-as.numeric(as.POSIXlt("2017-08-13 19:00:00"))
nm_taxa<-dim(CEmean)[1]
for (i in 1:nm_taxa){
  t1<-CEmean$StartTS[i]
  t2<-CEmean$EndTS[i]
  if (CEmean$StartTS[i]>time1 & CEmean$StartTS[i]<time2) CEmean$RH[i]<-mean(probe0_clean[which(probe0_clean$Timestamp>time1&probe0_clean$Timestamp<time2),4])
  else CEmean$RH[i]<-mean(probe6_clean[which(probe6_clean$Timestamp>t1&probe6_clean$Timestamp<t2),4])
}
hist(CEmean$RH,main="Distribution of RH SD17",xlab="Relative Humidity")
####
CEmean_noDish<-CEmean # use this line when correct raw CE rate by DW

####### Prepare field augmented analysis (with RH) ##########
check.list<-CEmean[which(grepl("Chk",CEmean[,1],fixed=TRUE)),1]
CEmean_noDish<-CEmean[-which(CEmean[,1] %in% check.list),]


FieldDesign<-read.table("/Users/menglin/Google Drive/MLC_AZ_2017/Pheno/field_design/AZ17_FieldDesign_Barcode_MLC.txt",header=T,sep="\t")
FieldDesign<-read.table("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\field_design\\AZ17_FieldDesign_Barcode_MLC.txt",header=T,sep="\t")


FieldDesign<-read.table("C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\field_design\\SD17_FieldDesign_Barcode_MLC.txt",header=T,sep="\t")
FieldDesign<-read.table("/Users/menglin/Google Drive/MLC_AZ_2017/Pheno/field_design/SD17_FieldDesign_Barcode_MLC.txt",header=T,sep="\t")


colnames(CEmean_noDish)[1]<-"MLC_mf"
forCEBLUP<-merge(CEmean_noDish,FieldDesign,by="MLC_mf",all=F)

colnm<-which(colnames(forCEBLUP)=="MLC_STANDARD")
for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,colnm]=="MO17") {forCEBLUP$CHECK[i]<-1}
  else if (forCEBLUP[i,colnm]=="N28HT") {forCEBLUP$CHECK[i]<-3} 
  else forCEBLUP$CHECK[i]<-99
}
for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,colnm]=="MO17" | forCEBLUP[i,colnm]=="N28HT") {forCEBLUP$IS_EXEXPERIMENTAL[i]<-0} else forCEBLUP$IS_EXEXPERIMENTAL[i]<-1
}
#forCEBLUP$IS_EXEXPERIMENTAL[which(as.character(forCEBLUP[,17]=="1"))] <- 0
#forCEBLUP$IS_EXEXPERIMENTAL[which(as.character(forCEBLUP[,17]=="99"))] <- 1
forCEBLUP1<-forCEBLUP[,c(9,7,12,13,15,16,21,22,2,3)]
#hist(forCEBLUP1[,9])
colnames(forCEBLUP1)<-c("GENOTYPE","RH","Year","LOC","COL","BLOCK","CHECK","IS_EXPERIMENTAL","rawCE","CEadj.mean")
forCEBLUP1$COL1<-forCEBLUP1$COL
## AZ17 only
forCEBLUP1$COL1[which(forCEBLUP1$BLOCK<=13&forCEBLUP1$COL>10)]<-21-forCEBLUP1$COL[which(forCEBLUP1$BLOCK<=13&forCEBLUP1$COL>10)]
forCEBLUP1$COL1[which(forCEBLUP1$BLOCK>13&forCEBLUP1$COL<11)]<-21-forCEBLUP1$COL[which(forCEBLUP1$BLOCK>13&forCEBLUP1$COL<11)]

write.table(forCEBLUP1,"C:\\Users\\ml2498\\Google Drive\\MLC_AZ_2017\\Pheno\\CE_BLUP\\SD17_CEadjLSA_roomRH_BLUPinput0413.txt",col.names=T,row.names=F,sep="\t",quote=FALSE)
write.table(forCEBLUP1,"/Users/menglin/Google Drive/MLC_AZ_2017/Pheno/CE_BLUP/AZ17_CEadjLSA_roomRH_BLUPinput2018.txt",col.names=T,row.names=F,sep="\t",quote=FALSE)






############# prepare Wet Weight augmented analysis #############
#fd.mo17.list<-CEmean[which(grepl("MO17",CEmean[,1],fixed=TRUE)),1]
#fd.n28ht.list<-CEmean[which(grepl("N28HT",CEmean[,1],fixed=TRUE)),1]
#fd.check.list<-c(fd.mo17.list,fd.n28ht.list)
#CEmean_noFdChk<-CEmean[-which(CEmean[,1] %in% fd.check.list),]

FieldDesign<-read.table("C:\\Meng\\MaizeLeafCuticle\\MLC_AZ_2017\\Pheno\\field_design\\AZ17_FieldDesign_Barcode_MLCstd.txt",header=T,sep="\t")
forCEBLUP<-merge(CEmean,FieldDesign,by="Geno_mdf",all.x=T,all.y=F)

forCEBLUP[] <- lapply(forCEBLUP, as.character)
dish.list<-forCEBLUP[which(grepl("Chk",forCEBLUP[,1],fixed=TRUE)),1]

forCEBLUP[which(forCEBLUP[,1] %in% dish.list),c(7,8)]<-"Chk"
forCEBLUP[which(forCEBLUP[,1] %in% dish.list),9]<-"SD"
forCEBLUP[which(forCEBLUP[,1] %in% dish.list),10]<-"17"
forCEBLUP[which(forCEBLUP[,1] %in% dish.list),11]<-"AZ"
forCEBLUP[which(forCEBLUP[,1] %in% dish.list),12]<-"521"
forCEBLUP[which(forCEBLUP[,1] %in% dish.list),13]<-"1"
forCEBLUP[which(forCEBLUP[,1] %in% dish.list),14]<-"1"
forCEBLUP[which(forCEBLUP[,1] %in% dish.list),15]<-"521"
forCEBLUP[which(forCEBLUP[,1] %in% dish.list),16]<-"2"
forCEBLUP[which(forCEBLUP[,1] %in% dish.list),17]<-"521"


for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,7]=="MO17") {forCEBLUP$CHECK1[i]<-1}
  else if (forCEBLUP[i,7]=="N28HT") {forCEBLUP$CHECK1[i]<-2} 
  else forCEBLUP$CHECK1[i]<-99
}

for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,7]=="MO17" | forCEBLUP[i,7]=="N28HT") {forCEBLUP$IS_EXEXPERIMENTAL1[i]<-0} else forCEBLUP$IS_EXEXPERIMENTAL1[i]<-1
}

for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,7]=="Chk") {forCEBLUP$CHECK2[i]<-1}
  else forCEBLUP$CHECK2[i]<-99
}

for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,7]=="Chk") {forCEBLUP$IS_EXEXPERIMENTAL2[i]<-0} else forCEBLUP$IS_EXEXPERIMENTAL2[i]<-1
}

forCEBLUP1<-forCEBLUP[,c(7,6,5,10,11,13,14,18:21,2)]
#hist(forCEBLUP1[,9])
colnames(forCEBLUP1)<-c("GENOTYPE","RH","WWDATE","Year","LOC","COL","BLOCK","CHECK1","IS_EXPERIMENTAL1","CHECK2","IS_EXPERIMENTAL2","CEadj.mean")
write.table(forCEBLUP1,"C:\\Meng\\MaizeLeafCuticle\\MLC_AZ_2017\\Pheno\\CE_BLUP\\AZ17_CEadjLSA_roomRH_BLUPinput_wDish.txt",col.names=T,row.names=F,sep="\t",quote=FALSE)







############## Option 2: Use RH in each quadrant ###############
setwd("C:\\Meng\\MaizeLeafCuticle\\MLC_AZ_2017\\PB_database\\FromNick")

CR300<-read.table("HMP155A_reformatted_0623.txt",header=T,sep=",")
CR300$Timestamp<-as.numeric(as.POSIXlt(CR300$TIMESTAMP))#convert local time to unix time
probe0 = data.frame(matrix(vector(), dim(CR300)[1], 4,dimnames=list(c(), c("Timestamp", "ID","Temperature" ,"RH"))),stringsAsFactors=F)
probe0$Timestamp<-CR300$Timestamp
probe0$ID<-CR300$RECORD
probe0$Temperature<-CR300$AirTC_Avg
probe0$RH<-CR300$RH

probe1<-read.csv("probe1_rpithon.csv",sep="\t",header=T)
probe3<-read.csv("probe3_Smaug.csv",sep="\t",header=T)
probe4<-read.csv("probe4_Scatha.csv",sep="\t",header=T)
probe6<-read.csv("probe6_Balerion.csv",sep="\t",header=T)

probe1_clean<-probe1[which(probe1$RH<110),]
probe3_clean<-probe3[which(probe3$RH<110),]
probe4_clean<-probe4[which(probe4$RH<110),]
probe6_clean<-probe6[which(probe6$RH<110),]
probe0_clean<-probe0[which(probe0$RH<110),]

WW_LSA[] <- lapply(WW_LSA, as.character)
time1<-as.numeric(as.POSIXlt("2017-06-18 16:00:00")) #EST 16:00 = AZ 13:00
time2<-as.numeric(as.POSIXlt("2017-06-18 17:00:00"))
nm_taxa<-dim(WW_LSA)[1]
Diff16<-0.651 # This is calculated in the Calibration.R
Diff36<-(-0.337)
Diff46<-0.494
i=1194
i=38
for (i in 1:nm_taxa){
  t1<-WW_LSA$StartTS[i]
  t2<-WW_LSA$EndTS[i]
    if (as.numeric(WW_LSA$StartTS[i])>time1 & as.numeric(WW_LSA$StartTS[i])<time2 &WW_LSA$Sensor[i]=="1") WW_LSA$RH[i]<-mean(probe0_clean[which(probe0_clean$Timestamp>t1&probe0_clean$Timestamp<t2),4])+Diff16
    else if (as.numeric(WW_LSA$StartTS[i])>time1 & as.numeric(WW_LSA$StartTS[i])<time2 &WW_LSA$Sensor[i]=="3") WW_LSA$RH[i]<-mean(probe0_clean[which(probe0_clean$Timestamp>t1&probe0_clean$Timestamp<t2),4])+Diff36
    else if (as.numeric(WW_LSA$StartTS[i])>time1 & as.numeric(WW_LSA$StartTS[i])<time2 &WW_LSA$Sensor[i]=="4") WW_LSA$RH[i]<-mean(probe0_clean[which(probe0_clean$Timestamp>t1&probe0_clean$Timestamp<t2),4])+Diff46
    else if (as.numeric(WW_LSA$StartTS[i])>time1 & as.numeric(WW_LSA$StartTS[i])<time2 &WW_LSA$Sensor[i]=="6") WW_LSA$RH[i]<-mean(probe0_clean[which(probe0_clean$Timestamp>t1&probe0_clean$Timestamp<t2),4])
    else if ((as.numeric(WW_LSA$StartTS[i])<time1 | as.numeric(WW_LSA$StartTS[i])>time2) & WW_LSA$Sensor[i]=="1") WW_LSA$RH[i]<-mean(probe6_clean[which(probe6_clean$Timestamp>t1&probe6_clean$Timestamp<t2),4])+Diff16
    else if ((as.numeric(WW_LSA$StartTS[i])<time1 | as.numeric(WW_LSA$StartTS[i])>time2) &WW_LSA$Sensor[i]=="3") WW_LSA$RH[i]<-mean(probe6_clean[which(probe6_clean$Timestamp>t1&probe6_clean$Timestamp<t2),4])+Diff36
    else if ((as.numeric(WW_LSA$StartTS[i])<time1 | as.numeric(WW_LSA$StartTS[i])>time2) &WW_LSA$Sensor[i]=="4") WW_LSA$RH[i]<-mean(probe6_clean[which(probe6_clean$Timestamp>t1&probe6_clean$Timestamp<t2),4])+Diff46
    else if ((as.numeric(WW_LSA$StartTS[i])<time1 | as.numeric(WW_LSA$StartTS[i])>time2) &WW_LSA$Sensor[i]=="6") WW_LSA$RH[i]<-mean(probe6_clean[which(probe6_clean$Timestamp>t1&probe6_clean$Timestamp<t2),4])
  }
#which(is.na(WW_LSA$RH))

#################################################

check.list<-WW_LSA[which(grepl("Chk",WW_LSA[,1],fixed=TRUE)),1]
WW_LSA_noDish<-WW_LSA[-which(WW_LSA[,1] %in% check.list),]
colnames(WW_LSA_noDish)[9]<-"Geno_mdf"


FieldDesign<-read.table("C:\\Meng\\MaizeLeafCuticle\\MLC_AZ_2017\\Pheno\\field_design\\AZ17_FieldDesign_Barcode_MLCstd.txt",header=T,sep="\t")
forCEBLUP<-merge(WW_LSA_noDish,FieldDesign,by="Geno_mdf",all=F)

for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,15]=="MO17") {forCEBLUP$CHECK[i]<-1}
  else if (forCEBLUP[i,15]=="N28HT") {forCEBLUP$CHECK[i]<-2} 
   else forCEBLUP$CHECK[i]<-99
}
for (i in 1:dim(forCEBLUP)[1]){
  if (forCEBLUP[i,15]=="MO17" | forCEBLUP[i,15]=="N28HT") {forCEBLUP$IS_EXEXPERIMENTAL[i]<-0} else forCEBLUP$IS_EXEXPERIMENTAL[i]<-1
}
#forCEBLUP$IS_EXEXPERIMENTAL[which(as.character(forCEBLUP[,17]=="1"))] <- 0
#forCEBLUP$IS_EXEXPERIMENTAL[which(as.character(forCEBLUP[,17]=="99"))] <- 1
forCEBLUP1<-forCEBLUP[,c(15,14,18,19,21,22,26,27,13)]
#hist(forCEBLUP1[,9])
colnames(forCEBLUP1)<-c("GENOTYPE","RH","Year","LOC","COL","BLOCK","CHECK","IS_EXPERIMENTAL","CEadj.mean")
write.table(forCEBLUP1,"C:\\Meng\\MaizeLeafCuticle\\MLC_AZ_2017\\Pheno\\CE_BLUP\\AZ17_CEadj_4QRH_BLUPinput.txt",col.names=T,row.names=F,sep="\t",quote=FALSE)

##############################################





############ add RH and Temperature ###########
probe3<-read.csv("~/Documents/Meng/MaizeLeafCuticle/MLC_Pilot/ProbeDB/probe3_Smaug.csv",sep="\t",header=T)
probe7<-read.csv("~/Documents/Meng/MaizeLeafCuticle/MLC_Pilot/ProbeDB/HMP155A_reformatted.csv",sep="\t",header=T)
probe3_clean<-probe3[which(probe3$RH<110&probe3$RH>20),]
probe7_clean<-probe7[which(probe7$RH<110&probe7$RH>20),]
#i=24
for (i in 1:dim(raw)[1]){
  lastpoint<-(max(which(!is.na(raw[i,c(3:16)])))+2)
  firstpoint<-(min(which(!is.na(raw[i,c(3:16)])))+3)
  if (raw$Design=="32"){
    raw$averRH[i]<-mean(probe3_clean[which(probe3_clean$Timestamp>raw[i,firstpoint]&probe3_clean$Timestamp<raw[i,lastpoint]),4])# not adjusted RH
    raw$averTemp[i]<-mean(probe3_clean[which(probe3_clean$Timestamp>raw[i,firstpoint]&probe3_clean$Timestamp<raw[i,lastpoint]),3])# raw T, not adjusted
  }
  else {
    raw$averRH[i]<-mean(probe7_clean[which(probe7_clean$Timestamp>raw[i,firstpoint]&probe7_clean$Timestamp<raw[i,lastpoint]),4])# HMP155A RH
    raw$averTemp[i]<-mean(probe7_clean[which(probe7_clean$Timestamp>raw[i,firstpoint]&probe7_clean$Timestamp<raw[i,lastpoint]),3])
  #raw$averRH1[i]<-mean(probe3_clean[which(probe3_clean$Timestamp>raw[i,4]&probe3_clean$Timestamp<raw[i,6]),5])
  #raw$averRH2[i]<-mean(probe3_clean[which(probe3_clean$Timestamp>raw[i,6]&probe3_clean$Timestamp<raw[i,8]),5])
  #raw$averRH3[i]<-mean(probe3_clean[which(probe3_clean$Timestamp>raw[i,8]&probe3_clean$Timestamp<raw[i,10]),5])
  #raw$averRH4[i]<-mean(probe3_clean[which(probe3_clean$Timestamp>raw[i,10]&probe3_clean$Timestamp<raw[i,12]),5])
  }}
########### 
RH<-raw$averRH # relative humidity in %, e.g. 70
Temp<-raw$averTemp+273.15
SA<- raw$SA/10000 # change the unit of leaf area from mm2 to m2
CE<-as.numeric(slope[,2]) # raw cuticular evaporation rate
deltapsi<-(-310)*log(RH)
raw$conductance<-SA*deltapsi/CE # calculate for each sample and average over genotype. Compare between RH levers
CEadj<-CE/SA
raw<-cbind(raw,CEadj)
##############################
#raw$Geno<-as.character(raw$Geno)
AW30<-apply(raw[which(raw$Geno=="ArmyWorm"&raw$Design=="32"),c(20:23)],2,mean)
AW50<-apply(raw[which(raw$Geno=="ArmyWorm"&raw$Design=="50"),c(20:23)],2,mean)
AW70<-apply(raw[which(raw$Geno=="ArmyWorm"&raw$Design=="75"),c(20:23)],2,mean)
AW50cuff<-apply(raw[which(raw$Geno=="ArmyWorm"&raw$Design=="Cuff_50"),c(20:23)],2,mean)

B85.30<-apply(raw[which(raw$Geno=="B85"&raw$Design=="32"),c(20:23)],2,mean)
B85.50<-apply(raw[which(raw$Geno=="B85"&raw$Design=="50"),c(20:23)],2,mean)
B85.70<-apply(raw[which(raw$Geno=="B85"&raw$Design=="75"),c(20:23)],2,mean)

CK30<-apply(raw[which(raw$Geno=="CheckDish"&raw$Design=="32"),c(20:23)],2,mean)
CK50<-apply(raw[which(raw$Geno=="CheckDish"&raw$Design=="50"),c(20:23)],2,mean)
CK70<-apply(raw[which(raw$Geno=="CheckDish"&raw$Design=="75"),c(20:23)],2,mean)

ND252.30<-apply(raw[which(raw$Geno=="ND252"&raw$Design=="32"),c(20:23)],2,mean)
ND252.50<-apply(raw[which(raw$Geno=="ND252"&raw$Design=="50"),c(20:23)],2,mean)
ND252.70<-apply(raw[which(raw$Geno=="ND252"&raw$Design=="75"),c(20:23)],2,mean)

summary<-rbind(AW30,AW50,AW50cuff,AW70, B85.30,B85.50,B85.70, CK30,CK50,CK70,ND252.30,ND252.50,ND252.70)
summary[,2]<-summary[,2]+273.15
write.table(summary,"17MLC_PT_CEsummary.txt",sep="\t")
###########################
(summary[8,4]-summary[9,4])*(summary[8,1]*log(summary[8,2])-summary[10,1]*log(summary[10,2]))/(summary[8,1]*log(summary[8,2])-summary[9,1]*log(summary[9,2]))
summary[8,4]-summary[10,4]
(summary[8,4]-summary[10,4])*(log(summary[8,2])-log(summary[9,2]))/(log(summary[8,2])-log(summary[10,2]))
summary[8,4]-summary[9,4]
(summary[8,4]-summary[9,4])*(log(summary[9,2])-log(summary[10,2]))/(log(summary[8,2])-log(summary[9,2]))
summary[9,4]-summary[10,4]

(summary[1,4]-summary[2,4])*(summary[1,1]*log(summary[1,2])-summary[4,1]*log(summary[4,2]))/(summary[1,1]*log(summary[1,2])-summary[2,1]*log(summary[2,2]))
summary[1,4]-summary[4,4]
(summary[1,4]-summary[4,4])*(log(summary[1,2])/(summary[1,2])-log(summary[2,2])/summary[2,2])/(log(summary[1,2])/(summary[1,2])-log(summary[4,2])/summary[4,2])
summary[1,4]-summary[2,4]
(summary[1,4]-summary[2,4])*(log(summary[2,2])-log(summary[4,2]))/(log(summary[1,2])-log(summary[2,2]))
(summary[1,4]-summary[2,4])*(log(summary[2,2])/summary[2,2]-log(summary[4,2])/summary[4,2])/(log(summary[1,2])/summary[1,2]-log(summary[2,2])/summary[2,2])

summary[2,4]-summary[4,4]


(summary[5,4]-summary[6,4])*(summary[5,1]*log(summary[5,2])-summary[7,1]*log(summary[7,2]))/(summary[5,1]*log(summary[5,2])-summary[6,1]*log(summary[6,2]))
summary[5,4]-summary[7,4]
(summary[5,4]-summary[7,4])*(log(summary[5,2])-log(summary[6,2]))/(log(summary[5,2])-log(summary[7,2]))
summary[5,4]-summary[6,4]
(summary[5,4]-summary[6,4])*(log(summary[6,2])-log(summary[7,2]))/(log(summary[5,2])-log(summary[6,2]))
summary[6,4]-summary[7,4]


(summary[11,4]-summary[12,4])*(summary[11,1]*log(summary[11,2])-summary[13,1]*log(summary[13,2]))/(summary[11,1]*log(summary[11,2])-summary[12,1]*log(summary[12,2]))
summary[11,4]-summary[13,4]
(summary[11,4]-summary[13,4])*(log(summary[11,2])-log(summary[12,2]))/(log(summary[11,2])-log(summary[13,2]))
summary[11,4]-summary[12,4]
(summary[11,4]-summary[12,4])*(log(summary[12,2])-log(summary[13,2]))/(log(summary[11,2])-log(summary[12,2]))
summary[12,4]-summary[13,4]



(ck30-ck50)*(T30*log(RH30)-T70*log(RH70))/(T30*log(RH30)-T50*log(RH50))
ck30-ck70

(ck30-ck70)*(T30*log(RH30)-T50*log(RH50))/(T30*log(RH30)-T70*log(RH70))
ck30-ck50

(ck30-ck50)*(T50*log(RH50)-T70*log(RH70))/(T30*log(RH30)-T50*log(RH50))
ck50-ck70
############# ANOVA ###################

data<-raw[which(raw$Design != "Cuff_50"),(17:23)]
data$logRH<-log(data$averRH)
setwd("~/Documents/Meng/MaizeLeafCuticle/MLC_Pilot")
write.table(data,"17MLC_pilot_forRH_sasinput_adjPT1.txt",row.names=F,sep="\t")

data2<-raw[which(raw$Design =="Cuff_50" | (raw$Design ==50 & raw$Geno=="ArmyWorm")),(17:23)]
write.table(data2,"17MLC_pilot_forLeafCuff_sasinput.txt",row.names=F)

data.aov1 = aov(data$CEadj ~ data$Geno * data$Design)
summary(data.aov1)
TukeyHSD(data.aov1)
#LSD.test(data.aov1,"data$Geno * data$Design")
#pairwise.t.test(data$CEadj, data$)

data.ancova<-aov(data$CEadj ~ data$averRH + as.factor(data$Geno) + data$averRH*as.factor(data$Geno))
#Anova(data.ancova,type="III")
summary(data.ancova)
data.ancova.pr<-proj(data.ancova)
residuals<-data.ancova.pr[,"Residuals"]#extract residuals from the model and see if it still related to RH (design)
residuals<-data.frame(residuals)
test<-aov(residuals[,1]~data$Design)
summary(test)

data.fit<-lm(data$CEadj ~ data$averRH + as.factor(data$Geno) + data$averRH*as.factor(data$Geno))
summary(data.fit)

data.aov2 = aov(data$conductance ~ data$Geno * data$Design)
summary(data.aov2)
TukeyHSD(data.aov2)
