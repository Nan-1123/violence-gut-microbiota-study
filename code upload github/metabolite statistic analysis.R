#################### metabolite analysis
library("readxl")
metdat0<-read_excel("met_ALL_data.xlsx") # metabolite clean data
# phenotype data
criphe<-read.table("criminal_phenotype.txt",header = T,encoding = "UTF-8",sep = "\t")
colnames(criphe)[3]<-"xlid"
criphe<-criphe[which(criphe$pheno!="" & criphe$xlid%in%colnames(metdat0)==T),]
table(criphe$pheno)

metdat<-metdat0[,c(1,match(criphe$xlid,colnames(metdat0)))]
write.csv(metdat,"VCmetdat.csv",row.names = F)

# OPLS-DA
install.packages("devtools")

BiocManager::install("RBGL")
BiocManager::install("BiocParallel")
BiocManager::install("edgeR")
BiocManager::install("fgsea")
BiocManager::install("impute")
BiocManager::install("pcaMethods")
BiocManager::install("MSnbase")
BiocManager::install("siggenes")

devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, buildvignettes = FALSE)
library("MetaboAnalystR")

install.packages("ggrepel")
install.packages("iheatmapr")
install.packages("ellipse")

rm(list = ls())
lab<-matrix(c("Label",criphe$pheno[match(colnames(metdat)[-1],criphe$xlid)]),
            nrow = 1,ncol = ncol(metdat))
colnames(lab)<-colnames(metdat)
metdat1<-rbind(lab,metdat)
colnames(metdat1)[1]<-"Sample"
metdat1[1,which(metdat1[1,]=="暴力型")]<-"vio"
metdat1[1,which(metdat1[1,]=="非暴力型")]<-"non"
write.csv(metdat1,"viomet_analysisdat.csv",row.names = F)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "viomet_analysisdat.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)
# The final dataset should contain no more than than 5000 variables for effective computing.
mSet<-FilterVariable(mSet, var.filter="iqr", qc.filter="F", var.cutoff=50)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", transNorm="LogNorm",scaleNorm="MeanCenter", "NULL", ratio=FALSE, ratioNum=20)

# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, "fdr", FALSE)
mSet<-PlotTT(mSet, imgName = "tt_0_", format = "png", dpi = 72, width=NA)

# Perform oPLS-DA analysis
mSet<-OPLSR.Anal(mSet, reg=TRUE)
# Create a 2D oPLS-DA score plot
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "png", dpi=72, width=NA, 1,2,0.95,1,0)
# Create a plot of the model overview
mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "png", dpi=72, width=NA)
# Perform and plot oPLS-DA permutation 
mSet<-OPLSDA.Permut(mSet, 200)
mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "png", dpi=72, width=NA)

##### select metabolite marker with vip >1 | test fdr q<=0.05
vipres<-read.csv("oplsda_vip.csv",header = T)
viplist<-vipres$X[which(vipres$V1>1)] 
ttres<-read.csv("t_test.csv",header = T) 
markerlist<-unique(c(viplist,ttres$X)) 
write.table(as.matrix(markerlist),"markerlist.txt",row.names = F)

length(which(ttres$X%in%viplist==T))
marker_vtsig<-ttres$X[which(ttres$X%in%viplist==T)]
write.table(marker_vtsig,"marker_vtsig.txt",row.names = F)

########## partial correlation between taxa and metabolites
library(readxl)
viophe0<-read_excel("pheno.xlsx")
viophe<-viophe0[which(viophe0$Group!="A"),c("Sample","Group")]

# CLR for metagenomic 
sigtaxclr0<-read.csv("sigtaxa_clrdat.csv",header = T)
sigtaxclr<-sigtaxclr[match(viophe$Sample,sigtaxclr0$xlid),] # 234 samples

# metabolite 50% NA filter & log & auto scale
metdat0<-read.csv("sigmet_dat.csv",header = T) # significant metabolites by VIP socre and q-value of t.test
metdat<-metdat0[,which(colnames(metdat0)%in%viophe$Sample==T)]
metdat<-data.frame(xlid=colnames(metdat),t(metdat),stringsAsFactors = F)
colnames(metdat)[-1]<-metdat0$Index

metdat_log <- log2(metdat[,-1])
metdatlogs <- scale(metdat_log, scale = TRUE, center = TRUE)
metdatlogs_all<-data.frame(xlid=metdat$xlid,metdatlogs,stringsAsFactors = F)

sigtaxclr<-sigtaxclr[match(metdatlogs_all$xlid,sigtaxclr$xlid),] 

# covariate data
mtdat<-read.table("mtdata.txt",header = T,sep = " ")
#adjusting for age, BMI, smoke status, and pre-incarceration drinking behavior
covdat<-mtdat[match(metdatlogs_all$xlid,mtdat$xlid),c(3:4,11:13)]

# partial correlation test
sigtaxlist<-colnames(sigtaxclr)[-1]
metaalllist<-colnames(metdatlogs_all)[-1]

tm_p<-matrix(,nrow = length(metaalllist),ncol = length(sigtaxlist))
tm_r<-matrix(,nrow = length(metaalllist),ncol = length(sigtaxlist))

sigpcor<-matrix(,nrow = 1,ncol = 4)
colnames(sigpcor)<-c("tax","metabo","r","p")

library(ppcor)
for (t in 1:length(sigtaxlist)) {
  for (m in 1:length(metaalllist)) {
    cordat<-data.frame(tax=sigtaxclr[,sigtaxlist[t]],
                       met=metdatlogs_all[,metaalllist[m]],
                       covdat[,-ncol(covdat)],stringsAsFactors = F)
    
    pcor_test <- pcor.test(cordat$met, cordat$tax, 
                           cordat[, 3:ncol(cordat)], method = "spearman") 

    tm_r[m,t]<-pcor_test$estimate
    tm_p[m,t]<-pcor_test$p.value
    
    if(pcor_test$p.value<=0.05){
      sigresi<-matrix(c(sigtaxlist[t],
                        metaalllist[m],
                        pcor_test$estimate,
                        pcor_test$p.value),nrow = 1,ncol = 4)
      colnames(sigresi)<-c("tax","metabo","r","p")
      sigpcor<-rbind(sigpcor,sigresi)
    }
  }
  print(t)
}
sigpcor1<-data.frame(sigpcor[-1,],stringsAsFactors = F)
sigpcor1$r<-as.numeric(sigpcor1$r)
sigpcor1$p<-as.numeric(sigpcor1$p)
write.csv(sigpcor1,"sigpcor_res.csv",row.names = F)

