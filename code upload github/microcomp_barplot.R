library(reshape2)
library(ggplot2)
library(ggprism)
library(plyr)
library(ggsci)

print("microcomp_barplot(mlevel,n,datapath,phedat,phenoname,pheno,grouplist,date)")
### group microbiome level barplot
microcomp_barplot<-function(mlevel,n,datapath,phedat,phenoname,pheno,grouplist,date){
  sub_dir<-paste(phenoname,"/group barplot",sep = "")
  output_dir<-file.path("./",sub_dir)
  
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }else{
    print(paste(sub_dir,"already exists!",sep = ""))
  }
  
  phedat<-na.omit(phedat[,c("xlid",pheno)])
  
  dat<-read.csv(datapath,header = T)
  groupdat<-data.frame(dat[,mlevel])
  colnames(groupdat)<-mlevel
  
  # calculate relative abundance for each group type(mean of all sample from each group)
  for (t in grouplist) {
    itid<-phedat$xlid[which(phedat[,pheno]==t)]
    itmm<-apply(dat[,which(colnames(dat)%in%itid==T)], 1,mean)
    groupdat<-data.frame(groupdat,itmm,stringsAsFactors = F)
    colnames(groupdat)[ncol(groupdat)]<-paste(pheno,"_",t,sep = "")
  }
  write.csv(groupdat,paste(phenoname,"/group barplot/",mlevel,"_group_relativedat.csv",sep = ""),row.names = F)
  
  # plot the barplot
  grdat<-groupdat
  grdat$asum<-apply(grdat[,-1],1,sum)
  grdat<-grdat[order(grdat$asum,decreasing = T),]
  
  grdat_top<-grdat[-which(grdat[,mlevel]=="Others"),]
  grdat_top<-grdat_top[1:n,]
  grdat_top<-rbind(grdat_top,rep(0,ncol(grdat_top)))
  grdat_top[n+1,1]<-"Others"
  grdat_top[n+1,-1]<-apply(grdat_top[,-1],2,function(x){1-sum(as.numeric(x[-(n+1)]))})
  grdat_top<-grdat_top[,-ncol(grdat_top)]
  
  df<- melt(grdat_top)
  names(df)[1:2] <- c("Taxonomy","sample")
  p1 <- ggplot(df, aes( x = sample,y=100 * value,fill = Taxonomy))+
    geom_col(position = 'stack', width = 0.6)+#geom_bar(position = "stack", stat = "identity", width = 0.6)+ 
    scale_y_continuous(expand = c(0,0))+# 
    labs(x="group",y="Relative Abundance(%)",
         fill="Taxonomy")+
    ggtitle(paste(phenoname," top10 ",mlevel,sep = "")) +
    guides(fill=guide_legend(keywidth = 1, keyheight = 1))
  # +ggsci::scale_fill_npg()
  ggsave(
    filename = paste(phenoname,"/group barplot/",mlevel,"_groupbarplot.png",sep = "")
  )
}