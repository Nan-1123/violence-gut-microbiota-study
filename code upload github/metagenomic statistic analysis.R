#################### metagenomic analysis
########## vio vs nonvio
### phenodat
library(readxl)
viophe0<-read_excel("../renalaysis_pheno.xlsx")
viophe0<-data.frame(viophe0)
viophe<-viophe0[which(viophe0$Group!="A"),c("Sample","Group")] # 234
table(viophe$Group) # nonvio 113; vio 121
colnames(viophe)<-c("xlid","Group")

date="date" # date for name of result file 

##### metagenomic analysis
### metagenomic data ready
source("K:\\mycode\\microdat_ready.R")
levellist<-c("kingdom","class","family","genus","species","order","phylum")
for (l in levellist) {
  level<-paste(l,"data",sep = "")
  dataready("4.TaxAnnotation/MAT/Relative/",
            level,BFphe,date,".")
}

### group barplot
source("K:\\mycode\\microcomp_barplot.R")
levellist<-c("kingdom","class","family","genus","species","order","phylum")
for (l in levellist[-1]) {
  microcomp_barplot(l,10,
                    paste(l,"data_dat0317.csv",sep = ""),
                    viophe,".","Group",c("vio","nonvio"),date)
}

### alpha diversity
source("K:\\mycode\\alpha_divanalysis.R")
alpha_index<-read.csv("alpha_diversity_indexdata.csv",header = T,stringsAsFactors = F) # different alpha_diversity index
alpha_distcheck(alpha_index,viophe,".","Group",c("vio","nonvio"),date)
alphares_group<-read.csv("Group/alpha divanalysis/alphares_groupdatdate.csv",header = T)
alphalist<-colnames(alphares_group)[2:7]
alpha_boxplot(alphalist,alphares_group,
              ".","Group","wilcox.test",list(c("vio","nonvio"))) # Non-parametric test for Non-normal distribution

### beta diviersity 
library(vegan)
library(reshape2)
source("K:\\mycode\\betadiversity_binary.R")
for (m in levellist) {
  fundif_analysis(m,
                  paste(m,"data_dat0317.csv",sep = ""),
                  viophe,"Group","nonvio","vio",".",1)
}

### LEFSE analysis
library(readxl)

library(vegan)
library(reshape2)
source("micro_lefse.R") 

colormatch<-data.frame(matrix(c("vio","nonvio","#FEB3AE","#70C1B3"),nrow = 2,ncol = 2),stringsAsFactors = F)
colnames(colormatch)<-c("groups","color")

speciesdat<-read.csv("speciesdata_dat0317.csv",header = T)
speciesdat<-speciesdat[,-1]

otu<-paste("otu_",0:(nrow(speciesdat)-1),sep = "")
rownames(speciesdat)<-otu

tax_table0<-read.csv("tax_table.csv",header = T)
rownames(tax_table0)<-otu
tax_table0<-tidy_taxonomy(tax_table0)

dir.create(file.path("./","vio"))
microlefse_analysis(microlevel='all',speciesdat,tax_table0,viophe,"Group",c("nonvio","vio"),"lefse","vio")

###### logistic regression with covariate for violent vs non-violent
# check the p significance of tax in a logistic regression model with 4 covariate
library(readxl)
sigtax<-read.csv("sigtax.csv",header = T) # significance tax from lesfe analysis
mtdat<-read.table("mtdata.txt",header = T,sep = " ") # covariate data
covdat<-mtdat[match(sigtax$sampleid,mtdat$xlid),-c(1,2,13)] 

for (i in 1:(ncol(sigtax)-1)) {
  mbi<-data.frame(sigtax[,c(12,i+1)],covdat[,c(1:2,9:10)],stringsAsFactors = F)
  mbi$y<-as.factor(mbi$y)
  mbi_fomula<-paste("y ~ ",paste(colnames(mbi)[-1],collapse = "+"),sep = "")
  mbi_lr<-glm(mbi_fomula,data = mbi,family = binomial(link = "logit"))
  
  if(i==1){
    sig_res<-data.frame(summary(mbi_lr)[["coefficients"]][2,],stringsAsFactors = F)
    colnames(sig_res)<-colnames(sigtax)[i+1]
  }else{
    mbi_res<-data.frame(summary(mbi_lr)[["coefficients"]][2,],stringsAsFactors = F)
    colnames(mbi_res)<-colnames(sigtax)[i+1]
    sig_res<-data.frame(sig_res,mbi_res,stringsAsFactors = F)
  }  
  
}
write.table(sig_res,"sig_res.txt")

########## personality traits
###### heatmap for spearman correlation between tax and personality traits
taxdatall<-read.csv("alltax_relativeabundance.csv",header = T) # relative abundance data
taxdatall$tax<-gsub(";",".",taxdatall$tax)
# siglist is significance tax from logistic regression analysis and KW test with personality traits
boxdat<-t(taxdatall[match(siglist,taxdatall$tax),match(sigtax$sampleid,colnames(taxdatall))])
colnames(boxdat)<-boxplotlist
boxdat1<-data.frame(boxdat,violent=sigtax$y,stringsAsFactors = F)

psydat0<-read.csv("psydat.csv",header = T) # personality traits score dataset
psydat<-psydat0[match(rownames(boxdat1),psydat0$xlid),]
cordat<-data.frame(tax=colnames(boxdat1)[1:14],
                   taxname=NA,
                   Violent=c("Violent","Violent",rep("Non-violent",12)),
                   ASPD=NA,
                   AvPD=NA,
                   Openness=NA,
                   Neuroticism=NA,stringsAsFactors = F)
aspdid<-psydat$xlid[which(is.na(psydat$PDQ_fsh)==F)]
avpdid<-psydat$xlid[which(is.na(psydat$PDQ_hb)==F)]
openid<-psydat$xlid[which(is.na(psydat$BF_kfx)==F)]
neuroid<-psydat$xlid[which(is.na(psydat$BF_sjz)==F)]

for (i in 1:14) {
  cordat$taxname[i]<-strsplit(cordat$tax[i],split = "\\.")[[1]][2]
  
  cordat$ASPD[i]<-cor(boxdat1[aspdid,i], psydat$PDQ_fsh[match(aspdid,psydat$xlid)], method = "spearman")
  cordat$AvPD[i]<-cor(boxdat1[avpdid,i], psydat$PDQ_hb[match(avpdid,psydat$xlid)], method = "spearman")
  cordat$Openness[i]<-cor(boxdat1[openid,i], psydat$BF_kfx[match(openid,psydat$xlid)], method = "spearman")
  cordat$Neuroticism[i]<-cor(boxdat1[neuroid,i], psydat$BF_sjz[match(neuroid,psydat$xlid)], method = "spearman")
}
write.csv(cordat,"spearmancorr.csv",row.names = F)

########## moderation analysis
# read the phenotype data
PDQphe<-read.csv("PDQphe.csv",header = T)

scaledat<-read.csv("scalesdat0320.csv",header = T)
xlid<-read.csv("xlid.csv",header = T)
BFphe<-data.frame(xlid=xlid$XLID[match(scaledat$ID,xlid$ID)],
                  scaledat[,19:23],stringsAsFactors = F)
colnames(BFphe)[-1]<-paste(colnames(BFphe)[-1],"t",sep = "")

library(readxl)
viophe0<-read_excel("pheno.xlsx")
viophe<-viophe0[which(viophe0$Group!="A"),c("Sample","Group")] # 234

# ready the moderation case
sigtaxres<-read.csv("sigtax_res.csv",header = T)
# 只筛选在logistic回归和lefse中均显著的菌
taxlist<-sigtaxres$sigtaxlist[which(sigtaxres$logitp<=0.05 | sigtaxres$clrlogitp<=0.05)]
sigtaxres_info<-sigtaxres[match(taxlist,sigtaxres$sigtaxlist),grep("_LDA",colnames(sigtaxres))]

moderate_case<-matrix(,nrow = 1,ncol = 2)
colnames(moderate_case)<-c("tax","x")

for (i in 1:nrow(sigtaxres)) {
  wi<-which(is.na(sigtaxres_info[i,])==F)
  if(length(wi)>0){
    mi<-matrix(,ncol = 2,nrow = length(wi))
    mi[,1]<-rep(taxlist[i],length(wi))
    mi[,2]<-colnames(sigtaxres_info)[wi]
    colnames(mi)<-c("tax","x")
    moderate_case<-rbind(moderate_case,mi)
  }
}
moderate_case<-moderate_case[-1,]
length(unique(moderate_case[,1])) #14

moderate_case[,2]<-gsub('_LDA', '', moderate_case[,2])

xtype<-as.matrix(data.frame(strsplit(moderate_case[,2],split = "_")))

moderate_case1<-data.frame(moderate_case,
                           xtype=xtype[1,],stringsAsFactors = F)

write.csv(moderate_case1,"moderate_case.csv",row.names = F)

# logistic interaction model > moderation analysis
taxdatall<-read.csv("alltax_relativeabundance.csv",header = T)
taxdatall$tax<-gsub(";",".",taxdatall$tax)
sigtaxdat<-t(taxdatall[match(taxlist,taxdatall$tax),match(viophe$Sample,colnames(taxdatall))])
colnames(sigtaxdat)<-taxlist

colnames(viophe)<-c("sampleid","y")
sigtaxdat<-data.frame(sigtaxdat,
                      viophe,stringsAsFactors = F)
#install.packages("visreg")
library("visreg")
moderate_case1$interactionp<-NA
for (i in 1:nrow(moderate_case1)) {
  dati<-sigtaxdat[,c("sampleid",moderate_case1$tax[i],"y")]
  colnames(dati)[2]<-"z"
  
  dati$y<-as.factor(dati$y)
  dati$z<-c(scale(dati$z, center=TRUE, scale=FALSE)) # 交互项仅做中心化
  
  if(moderate_case1$xtype[i]=="BF"){
    dati$x<-BFphe[match(dati$sampleid,BFphe$xlid),moderate_case1$x[i]]
    dati<-na.omit(dati)
    
    dati$x<-c(scale(dati$x, center=TRUE, scale=FALSE))
  }else{
    dati$x<-PDQphe[match(dati$sampleid,PDQphe$xlid),moderate_case1$x[i]]
    dati<-na.omit(dati)
  }
  # logit interaction
  fi<-glm(y~x+z+x*z,family = binomial(link = logit),data = dati)
  fi_coef<-data.frame(summary(fi)$coefficients)
  moderate_case1$interactionp[i]<-fi_coef["x:z","Pr...z.."]
  if(fi_coef["x:z","Pr...z.."]<=0.05){
    write.csv(fi_coef,paste(moderate_case1$tax[i],
                            moderate_case1$x[i],"glmres_sig.csv",sep = "_"),row.names = T)
    png(paste(moderate_case1$tax[i],
              moderate_case1$x[i],"siginteractionplot.png",sep = "_"))
    plot(visreg(fi,xvar ="x",by="z",plot =  F,
                breaks=c(mean(dati$z)-sd(dati$z),mean(dati$z),mean(dati$z)+sd(dati$z))),
         overlay = T,band = FALSE,
         xlab=moderate_case1$x[i],
         ylab="logit-vio")
    dev.off()
    print(i)
  }else{
    write.csv(fi_coef,paste(moderate_case1$tax[i],
                            moderate_case1$x[i],"glmres.csv",sep = "_"),row.names = T)
  }
  
}
moderate_case1$viogroup<-sigtaxres$Group[match(moderate_case1$tax,sigtaxres$sigtaxlist)]
write.csv(moderate_case1,"moderate_casesummary.csv",row.names = F)
