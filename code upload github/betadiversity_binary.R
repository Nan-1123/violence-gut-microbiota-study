library(reshape2)
library(ggplot2)
library(ggprism)
library(plyr)
library(ggsci)

library(patchwork)
library(tidyverse)
library(vegan)

library(devtools)
library(ggbiplot)
library(scatterplot3d)
library(ade4)

library(readxl)

print("fundif_analysis(microlevel,datapath,phedat,pheno,y1,y2,phenoname,datype)")
# fundif for two class
fundif_analysis<-function(microlevel,datapath,phedat,pheno,y1,y2,phenoname,datype){
  if(datype==1){
    phenoname<-paste(phenoname,"/micro_",sep = "")
  }else if(datype==2){
    phenoname<-paste(phenoname,"/fun_",sep = "")
  }else{
    print("datype must be 1 or 2")
    stop()
  }
  
  analysislist<-c("PCA plot","PCoA plot","NMDS plot")
  for(d in analysislist){
    sub_dir<-paste(phenoname,d,sep = "")
    output_dir<-file.path("./",sub_dir)
    
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }else{
      print(paste(sub_dir,"already exists!",sep = ""))
    }
  }
  
  mdat<-read.csv(datapath,header = T)
  dat<-mdat[,-1]
  rownames(dat)<-mdat[,1]
  
  group<-na.omit(phedat[which(phedat[,pheno]==y1 | phedat[,pheno]==y2),c("xlid",pheno)])
  colnames(group) <- c("samples","group")
  group<-group[which(group$samples%in%colnames(dat)==T),]
  
  dat<-dat[,match(group$samples,colnames(dat))]
  microdat<-data.frame(t(dat))
  micro.sum<-apply(microdat, 2, sum)
  if(length(which(micro.sum==0))==0){
    microdat1<-microdat
  }else{
    microdat1<-microdat[,-which(micro.sum==0)]
  }
  
  ##### PCA
  micro.pca<-prcomp(microdat1,scale. = T)
  # PCA 2d plot
  df_pca<-micro.pca
  df_pcs <-data.frame(df_pca$x, groups = group$group,stringsAsFactors = F)  
  percentage<-round(summary(micro.pca)$importance[2,]*100,digits = 2)
  percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
  p1<-ggplot(df_pcs,aes(x=PC1,y=PC2,color = groups))+ geom_point()+stat_ellipse(level = 0.95, show.legend = F) + 
    xlab(percentage[1]) +
    ylab(percentage[2]) + 
    geom_point(size= 0.5 )+
    labs(title=paste(microlevel,"PCA"), 
         subtitle=" PC1 and PC2 principal components ")
  ggsave(
    filename = paste(phenoname,"PCA plot/",microlevel,"_pca.png",sep = "")
  )
  # PCA 3d plot
  colorlist<-group$group
  colorlist<-replace(colorlist,which(colorlist==y1),'#f94144')
  #colorlist<-replace(colorlist,which(colorlist=="N"),'#f9c74f')
  colorlist<-replace(colorlist,which(colorlist==y2),'#5390d9')
  
  shapes<-group$group
  shapes<-replace(shapes,which(shapes==y1),16)
  #shapes<-replace(shapes,which(shapes=="N"),17)
  shapes<-replace(shapes,which(shapes==y2),18)
  shapes<-as.numeric(shapes)
  
  png( 
    filename = paste(phenoname,"PCA plot/",microlevel,"_pca3d.png",sep = ""), # 文件名称
    width = 480,           # 宽
    height = 480,          # 高
    units = "px",          # 单位
    bg = "white",          # 背景颜色
    res = 72)              # 分辨率
  # 2. 绘图
  scatterplot3d(micro.pca$x[,1:3], 
                color = colorlist, 
                #angle=10,
                cex.symbols=1.5,cex.axis=0.8, 
                pch = shapes,
                main=paste(pheno,microlevel,"PCA 3d plot"),
                xlab = paste("PCA1: ",round(summary(micro.pca)$importance[2,1]*100,digits = 2),"%",sep = ""), 
                ylab = paste("PCA2: ",round(summary(micro.pca)$importance[2,2]*100,digits = 2),"%",sep = ""), 
                zlab = paste("PCA3: ",round(summary(micro.pca)$importance[2,3]*100,digits = 2),"%",sep = ""))
  legend("right",legend = c(y1,y2),
         col = c('#f94144','#5390d9'),
         pch=c(16,17))
  # 3. 关闭画布
  dev.off()
  
  ### PCoA
  distance<-vegdist(microdat1,method = "bray")
  pcoa<-dudi.pco(distance,
                 scannf = F,   # 一种逻辑值，指示是否应该显示特征值条形图
                 nf=3)
  pc12<-pcoa$li[,1:2]
  pc123<-pcoa$li
  pc<-round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
  
  pc12 <- as.data.frame(pc12)
  pc12$samples <- row.names(pc12)
  
  pc123 <- as.data.frame(pc123)
  pc123$samples <- row.names(pc123)
  
  # PCoA 2d
  df <- merge(pc12,group,by="samples")
  #color=c("#1597A5","#FFC24B","#FEB3AE","#8B008B")
  #color=c("#1597A5","#FFC24B","#FEB3AE")
  p1<-ggplot(data=df,aes(x=A1,y=A2,
                         color=group,shape=group))+
    theme_bw()+
    geom_point(size=1.8)+
    stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
    theme(panel.grid = element_blank())+
    geom_vline(xintercept = 0,lty="dashed")+
    geom_hline(yintercept = 0,lty="dashed")+
    #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+
    #guides(color=guide_legend(title=NULL))+
    labs(x=paste0("PCo1 ",pc[1],"%",sep=""),
         y=paste0("PCo2 ",pc[2],"%",sep=""),
         title = paste(pheno,microlevel,"PCoA",sep = " "))+
    #scale_color_manual(values = color) +
    #scale_fill_manual(values = c("#1597A5","#FFC24B","#FEB3AE","#C1FFC1"))+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank())
  ggsave(
    filename = paste(phenoname,"PCoA plot/",microlevel,"_pcoa.png",sep = "")
  )
  
  # PCoA 3d
  df1 <- merge(pc123,group,by="samples")
  colorlist<-df1$group
  colorlist<-replace(colorlist,which(colorlist==y1),'#f94144')
  #colorlist<-replace(colorlist,which(colorlist=="N"),'#f9c74f')
  colorlist<-replace(colorlist,which(colorlist==y2),'#5390d9')
  
  shapes<-df1$group
  shapes<-replace(shapes,which(shapes==y1),16)
  #shapes<-replace(shapes,which(shapes=="N"),17)
  shapes<-replace(shapes,which(shapes==y2),18)
  shapes<-as.numeric(shapes)
  
  png( 
    filename = paste(phenoname,"PCoA plot/",microlevel,"_pcoa3d.png",sep = ""), # 文件名称
    width = 480,           # 宽
    height = 480,          # 高
    units = "px",          # 单位
    bg = "white",          # 背景颜色
    res = 72)              # 分辨率
  # 2. 绘图
  scatterplot3d(df1$A1, df1$A2, df1$A3, 
                color = colorlist, 
                #angle=10,
                cex.symbols=1.5,cex.axis=0.8, 
                pch = shapes,
                main=paste(pheno,microlevel,"PCoA 3d plot"),
                xlab = paste("PCoA1: ",pc[1],"%",sep = ""), 
                ylab = paste("PCoA2: ",pc[2],"%",sep = ""), 
                zlab = paste("PCoA3: ",pc[3],"%",sep = ""))
  legend("right",legend = c(y1,y2),
         col = c('#f94144','#5390d9'),
         pch=c(16,17))
  # 3. 关闭画布
  dev.off()
  
  ##### NMDS
  df_nmds <- metaMDS(distance, k = 2)
  df_nmds_stress <- df_nmds$stress
  write.table(matrix(df_nmds_stress,nrow = 1,ncol = 1),
              paste(phenoname,"NMDS plot/",microlevel,"_nmds_stress.txt",sep = ""),row.names = F)
  #提取作图数据
  df_points <- as.data.frame(df_nmds$points)
  #添加samp1es变量
  df_points$samples <- row.names(df_points)
  #修改列名
  names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
  df_nmds_group <-data.frame(df_points, groups = group$group[match(df_points$samples,group$samples)],stringsAsFactors = F)  
  p2<-ggplot(df_nmds_group,aes(x=NMDS1,y=NMDS2,color = groups))+ geom_point()+stat_ellipse(level = 0.95, show.legend = F) + 
    geom_point(size= 0.1 )+
    labs(title=paste(microlevel,"NMDS"))
  ggsave(
    filename = paste(phenoname,"NMDS plot/",microlevel,"_nmds.png",sep = "")
  )
  
  ##### ANOSIM
  df_anosim <- anosim(distance,df_nmds_group$groups,permutations = 999)
  df_anosim_res<-data.frame(x=df_anosim$class.vec,y=df_anosim$dis.rank)
  p3<-ggplot(df_anosim_res,aes(x=x,y=y))+
    stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
    geom_boxplot(aes(fill=x), 
                 outlier.colour="white",size=0.5)+
    theme(panel.background =element_blank(), 
          axis.line=element_line(),
          legend.position="none",plot.title = element_text(size=14))+
    scale_fill_manual(values=c("#FFC24B",'#f94144','#5390d9'))+ #指定颜色
    ggtitle("Bray-Curtis Anosim")+
    theme_prism(palette = "candy_bright",
                base_fontface = "plain",
                base_family = "serif", 
                base_size = 14,  
                base_line_size = 0.8, 
                axis_text_angle = 45)+
    theme(legend.position = 'none')+
    labs(x = paste("R=",df_anosim$statistic,", ","p=", df_anosim$signif),
         y = "Rank of Distance (Bray_Curtis)")
  ggsave(
    filename = paste(phenoname,"NMDS plot/",microlevel,"_anosim.png",sep = "")
  )
  ##### MRPP
  MRPP <- mrpp(distance,df_nmds_group$groups,permutations = 999)
  MRPPres<-data.frame(microlevel=microlevel,A=MRPP$A,P=MRPP$Pvalue,stringsAsFactors = F)
  write.table(MRPPres,
              paste(phenoname,"NMDS plot/",microlevel,"_MRPP_res.txt",sep = ""),row.names = F)
  
  ##### adonis
  adonis<-adonis2(microdat1~group, data=group, permutations = 999)
  adonisres<-data.frame(as.matrix(adonis),stringsAsFactors = F)
  write.table(adonisres,
              paste(phenoname,"NMDS plot/",microlevel,"_adonis_res.txt",sep = ""),row.names = F)
}

