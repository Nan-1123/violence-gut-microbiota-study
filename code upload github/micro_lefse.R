library(microeco)
library(ggtree)

print("microlefse_analysis(microlevel='all',speciesdat,tax_table0,phedat,y,grouplist,lmethod,phenoname)")
# use relative dat
microlefse_analysis<-function(microlevel="all",speciesdat,tax_table0,phedat,y,grouplist,lmethod,phenoname){
  sub_dir<-paste(phenoname,"/micro_lefse result",sep = "")
  output_dir<-file.path("./",sub_dir)
  
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }else{
    print(paste(sub_dir,"already exists!",sep = ""))
  }
  
  phedat<-phedat[which(is.na(match(phedat[,y],grouplist))==F),]
  # sample table
  sample_table<-na.omit(phedat[which(phedat[,"xlid"]%in%colnames(speciesdat)==T),c("xlid",y)])
  colnames(sample_table)<-c("sampleid","group")
  rownames(sample_table)<-sample_table$sampleid
  # feature table
  feature_table0<-speciesdat[,match(sample_table$sampleid,colnames(speciesdat))]
  feature_table0$sumab<-apply(feature_table0, 1, sum)
  feature_table<-feature_table0[which(feature_table0$sumab!=0),-ncol(feature_table0)]
  # tax table
  tax_table<-tax_table0[match(rownames(feature_table),rownames(tax_table0)),]
  
  dataset<-microtable$new(sample_table = sample_table,
                          otu_table = feature_table,
                          tax_table = tax_table)
  if(microlevel=="all"){
    a<-0.05
    pa<-"none"
  }else if(length(unique(tax_table[,microlevel]))<=10){
    a<-0.05
    pa<-"none"
  }else{
    a<-0.1
    pa<-"fdr"
  }
  
  # lefse 分析 taxa_level default "all"
  lefse<-trans_diff$new(dataset = dataset,
                        method = lmethod,
                        group = "group",
                        taxa_level = microlevel,
                        alpha = a,
                        lefse_subgroup = NULL,
                        p_adjust_method = pa)
  lefse_resdiff<-lefse$res_diff
  write.table(lefse_resdiff,paste(phenoname,"/micro_lefse result/",microlevel,"_",lmethod,"_",pa,a,".txt",sep=""),row.names = F,sep = ",")
  
  # lesfe 差异柱状图
  if(lmethod=="lefse"){
    if(nrow(lefse_resdiff)<30){
      un<-nrow(lefse_resdiff)
    }else{
      un<-30
    }
    png( 
      filename = paste(phenoname,"/micro_lefse result/",microlevel,"_",lmethod,"_bar",pa,a,".png",sep = ""),# 文件名称
      width = 720,           # 宽
      height = 480,          # 高
      units = "px",          # 单位
      bg = "white",          # 背景颜色
      res = 72)              # 分辨率
    # 2. 绘图
    lp<-lefse$plot_diff_bar(use_number=1:un,
                            width=0.8,
                            group_order=grouplist)
    print(lp)
    # 3. 关闭画布
    dev.off()
    
    # lesfe 分类树状图
    png( 
      filename = paste(phenoname,"/micro_lefse result/",microlevel,"_",lmethod,"_tree",pa,a,".png",sep = ""),# 文件名称
      width = 720,           # 宽
      height = 900,          # 高
      units = "px",          # 单位
      bg = "white",          # 背景颜色
      res = 72) 
    lefsetree<-lefse$plot_diff_cladogram(use_taxa_num=200,
                                         use_feature_num = 50,
                                         clade_label_level = 5,
                                         group_order = grouplist)#+ theme(text = element_text(size = 1))
    print(lefsetree)
    dev.off()
  }
}