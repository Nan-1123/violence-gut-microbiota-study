library(ggplot2)
library(ggsignif)
library(ggsci)

print("alpha_distcheck(alpha_index,phedat,phenoname,pheno,grouplist,date)")
print("alpha_boxplot<-function(alphalist=alphalist,alphares_group=alphares_group,pheno,phenoname,tmethod(t.test or wilcox.test),mycompare)")
print("mycompare must be a list!")

# 检验正态性+方差一致性
alpha_distcheck<-function(alpha_index,phedat,phenoname,pheno,grouplist,date){
  sub_dir<-paste(phenoname,"/alpha divanalysis",sep = "")
  output_dir<-file.path("./",sub_dir)
  
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }else{
    print(paste(sub_dir,"already exists!",sep = ""))
  }
  
  phedat<-na.omit(phedat[which(phedat[,pheno]%in%grouplist==T & phedat[,"xlid"]%in%alpha_index[,1]==T),c("xlid",pheno)])
  
  alphares_group<-data.frame(alpha_index[match(phedat$xlid,alpha_index[,1]),],phedat[,pheno],stringsAsFactors = F)
  #colnames(alphares_group)[2:4]<-c("chao1","ACE","obs")
  colnames(alphares_group)[ncol(alphares_group)]<-pheno
  write.csv(alphares_group,paste(phenoname,"/alpha divanalysis/alphares_groupdat",date,".csv",sep = ""),row.names = F)
  
  alphalist<-colnames(alphares_group)[2:7]
  
  for (a in alphalist) {
    checkdista<-matrix(,nrow = 1,ncol = 1+length(grouplist))
    colnames(checkdista)<-c("vartest",grouplist)
    rownames(checkdista)<-pheno
    
    y<-pheno
    # 组间方差对齐检验；若p值小于0.05，可以认为方差不齐
    btdat<-na.omit(alphares_group[,c(a,y)])
    vartest<-bartlett.test(btdat[,a],btdat[,y])
    checkdista[1,"vartest"]<-vartest$p.value
    
    typelist<-grouplist
    for(t in typelist){
      # shapiro test for shannon index of each group; 若p小于0.05则不满足正态分布
      itshap<-shapiro.test(alphares_group[which(alphares_group[,y]==t),a])
      checkdista[1,which(colnames(checkdista)==t)]<-itshap$p.value
    }
    colnames(checkdista)<-paste(a,"_",colnames(checkdista),sep = "")
    if(a==alphalist[1]){
      alpha_checkdist<-data.frame(checkdista,stringsAsFactors = F)
    }else{
      alpha_checkdist<-data.frame(alpha_checkdist,checkdista,stringsAsFactors = F)
    }
  }
  write.csv(alpha_checkdist,paste(phenoname,"/alpha divanalysis/alpha_checkdist",date,".csv",sep = ""),row.names = F)
}

# alpha boxplot
alpha_boxplot<-function(alphalist=alphalist,alphares_group=alphares_group,pheno,phenoname,tmethod,mycompare){
  for (a in alphalist) {
    y<-pheno
    alphadat<-na.omit(alphares_group[,c(a,y)])
    colnames(alphadat)<-c("alpha","groups")
    alphadat$groups<-as.character(alphadat$groups)
    groups<-unique(alphadat$groups)
    colorlist<-colormatch$color[match(groups,colormatch$groups)]
    alpha_boxplot<-ggplot(alphadat,aes(x=groups,y=alpha,fill=groups))+
      geom_boxplot()+#箱线图
      labs(title=paste(y,"Alpha diversity"), x="Group", y=a)+# 坐标轴标题啥的
      scale_fill_manual(values = colorlist[1:length(groups)])+#挑了两个喜欢的颜色;对应组别失败
      #scale_color_manual(values =c("#FEB3AE","#FFC24B","#70C1B3"),labels=c("H","N","L"),guide="legend")+ # 无法识别#号对应的颜色
      theme(plot.title=element_text(hjust=0.5), legend.title=element_blank())+ #去除背景网格
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+#去除背景网格
      geom_signif(comparisons = mycompare, test = tmethod)
    #+scale_fill_simpsons() # ggsci 的配色
    alpha_boxplot
    ggsave(
      filename = paste(phenoname,"/alpha divanalysis/",a,"_boxplot.png",sep = "")
    )
  }
}