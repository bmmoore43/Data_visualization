#percentage barplot
library("ggplot2")
library("scales")
library("reshape2")

#read in table with feature_presence, feature, gene_type, percent
args = commandArgs(TRUE)
setwd(args[1])
nc= args[2]
nc="datatable_percent_logratio.txt"

p1<-read.csv(nc, header=T, sep=',')
p1<-read.table(nc, header=T, sep='\t', na.strings=c("","NA")) 

nb=paste(c(basename(nc),"_percent.pdf"),collapse='')
nb
p2 <- subset(p1, class == "BM_GM"|class == 'BM_SM')
p2 <- subset(p1, class == 'GM-GM'|class == 'aaa'|class == 'SM-SM')
p3 <- subset(p2, feature == "p450"| feature == "NAD_binding_8"| feature == "UDPGT")

p2<- melt(p2)
#p1$feature <- factor(p1$feature, levels = p1$feature[order(p1$log2ratio, decreasing=TRUE)])
p2 <- transform(p2, variable=factor(variable, levels = unique(variable))) #maintains order of dataframe read in
ord <- c("node_c_0", "node_c_1","node_c_3", "node_c_4","node_c_5", "node_c_6", "node_c_7", "node_c_8", "node_c_10", "node_c_12", "node_c_14", "node_c_17", "node_c_19", "node_c_22", "node_c_24", "node_single")
ord2 <- rev(ord)
ord2
length(p2[,"feature"])
length(p2[,"class"])
length(p2[,"percent_pos"])

p2

p2<- transform(p2, feature=factor(feature, levels = ord)) #maintains order based on ord
p2<- transform(p1, type=factor(type, levels = unique(type))) #maintains order of dataframe read in
#p1 <- melt(p1)
p2
#logratio
data <- p1[,c(1,4,6)]#subset data
#bar plot of subsetted data
plot1<-ggplot(p1,aes(x=model,y=percent.SM, fill=Class))+ scale_fill_manual(values = alpha(c("turquoise3", "blueviolet", "cornflowerblue", "indianred2"))) +
   stat_summary(fun.y=mean,geom="bar") + coord_flip()
plot1
#make pdf file
nb=paste(basename(nc),"_logratio.pdf",collapse='')
pdf(file = nb)
plot1
dev.off()
#percent
#subset percent data
p2<-p1[,c(1,2,3)]
p2
#melt function
p2 <- melt(p1)
p2
#OR subset the specific features you want
p2<-p2[,c(1,2,4)]
p3 <- subset(p2, feature=="c200_stressfc_run8_cl_111"|feature=="c200_stressfc_run8_cl_87"|feature=="h50_dev_average_cl_2")

p2a<- transform(p2a, feature=factor(feature, levels = ord2))
#melt data
p2a <- melt(p2)
#pdf file name
nb=paste(c(basename(nc),"_percent.pdf"),collapse='')

#############
pdf(file = nb)

ggplot(p2a,aes(x = feature, y = value, fill = variable, main = "Percent Duplication Point")) + 
  geom_bar(position = "dodge",stat = "identity") + theme(plot.background = element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size =8.5),
  axis.text.x  = element_text(angle=45, vjust=0.5, colour= "black", size= rel(1.5))) + 
  scale_y_continuous() + xlab("Feature") + ylab("Percent genes at duplication point") +
  coord_flip()
dev.off()

pdf(file = nb)

