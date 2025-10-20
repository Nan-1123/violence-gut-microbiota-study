##### microbiome data ready for analysis
print("dataready(path,level,phedat,date,savepath)")
dataready<-function(path,level,phedat,date,savepath){
  dat0<-read.table(paste(path,level,".txt",sep = ""),header = T,sep = "\t")
  dat<-dat0[,c(1,as.numeric(na.omit(match(phedat$xlid,colnames(dat0)))))]
  write.csv(dat,paste(savepath,"/",level,"_dat",date,".csv",sep = ""),row.names = F)
}