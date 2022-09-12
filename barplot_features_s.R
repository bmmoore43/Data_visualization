#data is feature list with associated pos and neg weights
#args = commandArgs(TRUE)
#setwd(args[1])
#nc= args[2]

#read in data
nc = "Slyc_Athrem-RF_Slyc-Athremshared2-class_imp_scaled_top50.txt"
fcq1<-read.table(nc, header=T, sep='\t') 
#write.table(fcq1, file = "Slyc_Athrem-RF_500feat_imp_sclaed_top50.txt", sep="\t")

nd=paste(c(basename(nc),"_barplot.pdf"), collapse ="")
nd
library("ggplot2")
#keep sorted order
fcq1<- transform(fcq1, feature=factor(feature, levels = unique(feature))) #maintains order of dataframe read in

#barplot

pdf(file = nd)

p.b<-ggplot(fcq1,aes(x=feature,y=imp_score_scaled, fill= enrichment)) + ggtitle("Important scores for Final Model SM vs GM")+
  xlab("Feature") + ylab("Importance score scaled") +geom_bar(stat="identity")+
  stat_summary(fun.y=mean, geom="bar", position = "dodge")+
  theme(plot.background = element_blank(),text = element_text(size =9),
  axis.text.x  = element_text(angle=90, vjust=0.4, colour= "black"), axis.text.y  = element_text(vjust=0.75, colour= "black")) + 
  scale_y_continuous() + coord_flip()+
  scale_fill_manual(values = c("turquoise3","indianred2"))
p.b
   
dev.off()
## without fill
pdf(file = nd)

ggplot(fcq1,aes(x=feature,y=importance.score)) + ggtitle("Important scores for SM vs GM")+
  xlab("Feature") + ylab("Importance score scaled") +geom_bar(stat="identity")+
  stat_summary(fun.y=mean, geom="bar", position = "dodge")+
  theme(plot.background = element_blank(),text = element_text(size =9),
        axis.text.x  = element_text(angle=90, vjust=0.4, colour= "black"), axis.text.y  = element_text(vjust=0.75, colour= "black")) + 
  scale_y_continuous() + coord_flip()


dev.off()
