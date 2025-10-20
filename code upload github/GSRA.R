install.packages("ReporterScore")
library(ReporterScore)

rm(list = ls())
# relative abundance data of functional profile
ko_raw<-read.csv("keggko.csv",header = T)
# y
library(readxl)
viophe0<-read_excel("../renalaysis_pheno.xlsx")
viophe<-viophe0[which(viophe0$Group!="A"),c("Sample","Group")] # 234

# data format ready for reporter_score function
kodf<-ko_raw[,match(viophe$Sample,colnames(ko_raw))]
rownames(kodf)<-ko_raw$KOs

metadata<-data.frame(group=viophe$Group,stringsAsFactors = F)
rownames(metadata)<-viophe$Sample

# GSRA analysis
reporter_result <- reporter_score(
  kodf,
  "group",
  metadata = metadata,
  mode = c("directed", "mixed")[1],
  type = c("pathway", "module")[1],
  p.adjust.method2 = "fdr"
)

write.csv(reporter_result$reporter_s,"reporter_result.csv",row.names = F)
write.csv(reporter_result$ko_stat,"gsra_kostat.csv",row.names = F)

library(dplyr)
library(ggplot2)
plot_report(reporter_result,rs_threshold = c(-5,5),y_text_size=10)+
  labs(title = "violent vs nonviolent")
plot_report(reporter_result,rs_threshold=c(-1.96,1.96),mode = 2,y_text_size=5)+
  labs(title = "violent vs nonviolent")

plot_report_bar(reporter_result, 
                rs_threshold = 5, 
                facet_level = TRUE,y_text_size=10,)+
  theme(axis.text.x=element_text(vjust=1,size=10))+
  labs(title = "violent vs nonviolent")

p<-plot_report(reporter_result,
            rs_threshold=1.96,
            mode = 2,
            y_text_size=5,
            facet_level = TRUE)+
  theme(axis.text.x=element_text(vjust=1,size=5))+
  labs(title = "")
p
ggsave("gsra_t196bubble.png",p,dpi=300)
