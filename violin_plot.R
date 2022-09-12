#ggplot violin plot
#write directory name in setwd()
setwd()
# load ggplot and reshape packages
library(reshape2)
library(ggplot2)
##filename
nc = "Aracyc_SMSMvsGMGMvsAAA.txt_cont-cat-bin.combined_matrix4.0_NAimputed_only.txt"
##read data
p1 = read.table(nc, header=T, sep="\t",row.names=1, na.strings=c("","NA"))
#write data (if needed)
write.table(p1, file = "Aracyc_SMvsGMvsSMGM_20180723_only.txt_continuous_matrix_4.0_mod.txt", sep="\t")
p1 = read.table(nc, header=T, sep="\t")
#p1 <- na.omit(p1)
##subset data
newdata <- subset(p1, Class == 'GM'|Class == 'SM'|Class == 'Slyc-GM'| Class == 'Slyc-SM') 
newdata <- subset(p1, Class.Predicted == 'GM-GM'| Class.Predicted == 'aaa'|Class.Predicted == 'SM-SM')
newdata <- subset(p1, Class == 'Slyc-GM_GM_Slyc-GM'|Class == 'Slyc-GM_GM_Slyc-SM'|Class == 'Slyc-GM_SM_Slyc-GM'|Class == 'Slyc-GM_SM_Slyc-SM'|Class == 'Slyc-SM_GM_Slyc-GM'|Class == 'Slyc-SM_GM_Slyc-SM'|Class == 'Slyc-SM_SM_Slyc-GM'|Class == 'Slyc-SM_SM_Slyc-SM')
unique(p1$class2)
##transform to maintain order 
df<- transform(newdata, Class=factor(Class, levels = c('GM','SM','Slyc-GM','Slyc-SM')))
df<- transform(p1, Class.Predicted=factor(Class.Predicted, levels = c('GM-GM','aaa','SM-SM')))
df<- transform(p1, Class=factor(Class, levels = c('PM','SM-PM','SM')))
df<- transform(newdata, Class=factor(Class, levels = c('Slyc-GM','Slyc-SM'))) #'up_0_dwn_2','up_1_dwn_0','up_mulit_dwn_0',
df<- transform(newdata, Class=factor(Class, levels = c('Slyc-GM_GM_Slyc-GM','Slyc-GM_GM_Slyc-SM','Slyc-GM_SM_Slyc-GM','Slyc-GM_SM_Slyc-SM','Slyc-SM_GM_Slyc-GM','Slyc-SM_GM_Slyc-SM','Slyc-SM_SM_Slyc-GM','Slyc-SM_SM_Slyc-SM')))
#Convert to "long-form"
melted.df2<-melt(df)
#omit NAs
#new.mdf <- na.omit(melted.df2)

#merge class and variable ##not needed mostly
#library(tidyr)
#df = data.frame(new.mdf$Class, new.mdf$variable, new.mdf$value)
#df2 <- unite(df, newcol, c(new.mdf.Class, new.mdf.variable), remove=FALSE)

#make output file name
nb=paste(c(basename(nc),"_ltm-log_violin.pdf"),collapse='')

#pdf for one plot
#new.mdf1 <- subset(melted.df2, variable=='Fold_changes_for_2017_pathosystem.txt.maxSP.txt_GM')
#new.mdf1 <- na.omit(new.mdf1)
pdf(file = nb) #8 classes
p.b1 <-ggplot(df, aes(x=factor(class2), y=Median, fill=factor(class2))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE, alpha=0.5) + 
  scale_fill_manual(values =c("darkcyan","cyan3","dodgerblue2","darkorchid","goldenrod1","mediumvioletred","firebrick1","indianred2"))+
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  ggtitle("SM score")+ theme(axis.text.x  = element_text(angle=45, vjust=0.5, colour= "black", size= rel(1)))
p.b1
dev.off()

nb=paste(c(basename(nc),"_violin2.pdf"),collapse='')

df2 <- df[,c(3,11)]
df2 <- na.omit(df2)
pdf(file = nb) #4 classes
p.b1 <-ggplot(df2, aes(x=factor(class2), y=Median, fill=factor(class2))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE, alpha=0.5) + 
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","firebrick1"))+ #"indianred2"
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  ggtitle("SM score junction")+ theme(axis.text.x  = element_text(angle=45, vjust=0.5, colour= "black", size= rel(1)))
p.b1
dev.off()

pdf(file = nb) #1 class
df= newdata
p.b1 <-ggplot(df, aes(x=factor(Class), y=Median.SM.score, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE, alpha=0.5) + 
  scale_fill_manual(values =c("yellow1"))+
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  ggtitle("SM score")
p.b1
dev.off()
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

new.mdf1 <- subset(melted.df2, variable=='FamilySize')
new.mdf1 <- na.omit(new.mdf1)

p.b1a <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Gene Family Size") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b1a

new.mdf1 <- subset(melted.df2, variable=='hormone_GM_PCC')
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("PCC to GM- Hormone") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b2

new.mdf1 <- subset(melted.df2, variable== 'abiotic.shoot_expr_breadth_uponly')
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("abiotic expression breadth") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b3

new.mdf1 <- subset(melted.df2, variable=='biotic_expr_breadth_uponly')
new.mdf1 <- na.omit(new.mdf1)

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("biotic expression breadth") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b1

new.mdf1 <- subset(melted.df2, variable== "develop_PM_corr")

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Corr to GM develop") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b2

new.mdf1 <- subset(melted.df2, variable== "abiotic_PM_corr")

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Corr to GM abiotic") +
  stat_summary(fun=median, fun.min = median, fun.max = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b3

new.mdf1 <- subset(melted.df2, variable== "maxPCC_biotic")

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("max PCC under biotic") +
  stat_summary(fun=median, fun.min = median, fun.max = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b4

nb=paste(c(basename(nc),"_expr_violin.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1,p.b2, p.b3, p.b4, cols=2)
dev.off()

new.mdf1 <- subset(melted.df2, variable=='aaLength.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("amino acid length") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b1

new.mdf1 <- subset(melted.df2, variable=='ptri_kaks.insig.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("dN/dS to Poplar") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b2

new.mdf1 <- subset(melted.df2, variable=='pi.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("pi- nucleotide diversity") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b3

new.mdf1 <- subset(melted.df2, variable== "gene_family_size_orthomcl_1.5.continuous.MLD")

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("gene family size") +
  stat_summary(fun=median, fun.min = median, fun.max = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b4

nb=paste(c(basename(nc),"_misc_violin.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1,p.b2, p.b3, p.b4, cols=2)
dev.off()

new.mdf1 <- subset(melted.df2, variable== "expression_atgen.breadth_99.continuous.MLD")

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Expr breadth develop") +
  stat_summary(fun=median, fun.min = median, fun.max = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b1

new.mdf1 <- subset(melted.df2, variable== "expression_atgen.unlogged.median.continuous.MLD")

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Median expr level") +
  stat_summary(fun=median, fun.min = median, fun.max = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b2

new.mdf1 <- subset(melted.df2, variable== "functional_likelihood")

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("functional likelihood") +
  stat_summary(fun=median, fun.min = median, fun.max = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b3

new.mdf1 <- subset(melted.df2, variable== "develop_SM_corr")

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class.Predicted), y=value, fill=factor(Class.Predicted))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("PCC to SM under development") +
  stat_summary(fun=median, fun.min = median, fun.max = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","mediumorchid2","indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())#
p.b4

nb=paste(c(basename(nc),"_expr-func_violin.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1,p.b2, p.b3, p.b4, cols=2)
dev.off()

new.mdf1 <- subset(melted.df2, variable=='develop_expr_breadth')
new.mdf1 <- na.omit(new.mdf1)

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Expression breadth- development") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2")) #
p.b1

new.mdf1 <- subset(melted.df2, variable=='develop_SM_corr') 
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("PCC to SM- development") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2")) #"dodgerblue2",
p.b2

new.mdf1 <- subset(melted.df2, variable=='stress_SM_corr_max')
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("PCC to SM- stress") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2"))
p.b3

new.mdf1 <- subset(melted.df2, variable=='hormone_SM_corr')
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Correlation to SM- hormone") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2"))
p.b4

new.mdf1 <- subset(melted.df2, variable=='Ccanephora_maxKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("C.canephora max KaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2"))
p.b5

new.mdf1 <- subset(melted.df2, variable=='Alyrata_medKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("A.lyrata med Ka/Ks") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2"))
p.b6

new.mdf1 <- subset(melted.df2, variable=='Osativa_medKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b7 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("O. sativa med Ka/Ks") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2"))
p.b7

new.mdf1 <- subset(melted.df2, variable=='Spennellii_medKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b8 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("S. pennellii medKaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan2","cyan4","dodgerblue","dodgerblue4","mediumslateblue","mediumorchid2","indianred4","indianred2"))
p.b8

nb=paste(c(basename(nc),"_expr-evo-Ath-Slyc_violin.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1,p.b2, p.b3, p.b4, p.b5, p.b6, p.b7, p.b8, cols=2)
dev.off()



new.mdf1 <- subset(melted.df2, variable=='Osativa_medKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("O. sativa med Ka/Ks") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","indianred2","cyan4","indianred3"))
p.b2

nb=paste(c(basename(nc),"_kaks_violin.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1, p.b2, cols=2)
dev.off()

nb=paste(c(basename(nc),"_kaks_violin.pdf"),collapse='')
new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Coffea_canephora__cross_sp_orthologs_maxKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("C. canephora max KaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b1

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Alyrata_384_v2.1__cross_sp_orthologs_medKaKs') 
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("A. lyrata") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b2

new.mdf1 <- subset(melted.df2, variable=='tomato_paralogs_maxKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("tomato_paralogs_maxKaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b3

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__BrapaFPsc_277_v1.3__cross_sp_orthologs_medKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("B. rapa medKaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b4

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Smoellendorffii_91_v1.0__cross_sp_orthologs_maxKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("S. moellendorffii, max KaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b5

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Ptrichocarpa_210_v3.0__cross_sp_orthologs_maxKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("P. trichocarpa, maxKaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b6

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Ppatens_318_v3.3__cross_sp_orthologs_maxKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b7 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("P.patens, maxKaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b7

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Acoerulea_322_v3.1__cross_sp_orthologs_maxKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b8 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("A. coerulea, maxKaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b8

pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, p.b7, p.b8, cols=2)
dev.off()

nb=paste(c(basename(nc),"_kaks2_violin.pdf"),collapse='')

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Atrichopoda_291_v1.0__cross_sp_orthologs_maxKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("A. trichopoda, maxKaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b1

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Osativa_323_v7.0__cross_sp_orthologs_maxKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("O. sativa max Ka/Ks") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b2

new.mdf1 <- subset(melted.df2, variable=='NeighborParalogsCount') 
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Neighbor paralogs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b3

new.mdf1 <- subset(melted.df2, variable=='NeighborSMCount')
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Neighbors SM") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b4

new.mdf1 <- subset(melted.df2, variable=='FamilySize')
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Family Size, log") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b5

new.mdf1 <- subset(melted.df2, variable=='NeighborGMCount')
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("GM Neighbors") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b6

new.mdf1 <- subset(melted.df2, variable=='NeighborParalogsCount') 
new.mdf1 <- na.omit(new.mdf1)

p.b7 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Neighbor paralogs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b7

new.mdf1 <- subset(melted.df2, variable=='NeighborSMCount')
new.mdf1 <- na.omit(new.mdf1)

p.b8 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Neighbors SM") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b8

nb
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, p.b7, p.b8, cols=2)
dev.off()

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Solanum_pennellii__cross_sp_orthologs_medKaKs')
new.mdf1 <- na.omit(new.mdf1)

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("S. pennellii, medKaKs") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b1

new.mdf1 <- subset(melted.df2, variable=='X2015_circadian.txt_max_x')
new.mdf1 <- na.omit(new.mdf1)

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Circadien Max") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b1

new.mdf1 <- subset(melted.df2, variable=='hormone_Sly_20180124.txt_max') 
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Hormone max FC") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b2

new.mdf1 <- subset(melted.df2, variable=='X2015_time_course.txt_median_x')
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Time course med FC") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b3

new.mdf1 <- subset(melted.df2, variable=='all_combination_Sly_20180124.txt_max')
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("All combinations max FC") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b4

new.mdf1 <- subset(melted.df2, variable=='X2017_yield.txt_max_x')
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Yield max FC") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b5

new.mdf1 <- subset(melted.df2, variable=='mutation_Sly_20180124.txt_max')
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Mutation max FC") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p.b6

nb=paste(c(basename(nc),"_expr_violin.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, cols=2)
dev.off()

new.mdf1 <- subset(melted.df2, variable=='X2017_ltm.txt_MAD_x')
new.mdf1 <- na.omit(new.mdf1)

p.b7 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("LTM MAD") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b7

new.mdf1 <- subset(melted.df2, variable=='X2016_fruit_ripening.txt_MAD_x') 
new.mdf1 <- na.omit(new.mdf1)

p.b8 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Fruit Ripening MAD") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b8

new.mdf1 <- subset(melted.df2, variable=='X2016_inflorescence.txt_median_y')
new.mdf1 <- na.omit(new.mdf1)

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Log Med Inflorescence") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b1

new.mdf1 <- subset(melted.df2, variable=='X2016_inflorescence.txt_max_y')
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Log Max Inflorescence") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b2

new.mdf1 <- subset(melted.df2, variable=='X2015_root.txt_max_y')
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Log Max FC Root") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b3
new.mdf1 <- subset(melted.df2, variable=='X2015_root.txt_median_y')
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=log10(value), fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Log Median FC Root") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b4

new.mdf1 <- subset(melted.df2, variable=='Spearman_Results_Fold_changes_for_2015_time_courseSM_GM.max')
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Spearmans to DA, Time course") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b5

new.mdf1 <- subset(melted.df2, variable=='PCC_Results_Fold_changes_for_2015_time_courseSM_GM.max')
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("PCC to DA, Time course") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b6

nb=paste(c(basename(nc),"_expr2_violin.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, p.b7, p.b8, cols=2)
dev.off()


new.mdf1 <- subset(melted.df2, variable=='corpcor_Results_Fold_changes_for_2015_DC3000_GM.med')
new.mdf1 <- na.omit(new.mdf1)

p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("CorpCor to GM, DC3000") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b1

new.mdf1 <- subset(melted.df2, variable=='PCC_Results_Fold_changes_for_2013_bacterialSM_GM.med') 
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("PCC to DA, Bacterial FC") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b2

new.mdf1 <- subset(melted.df2, variable=='corpcor_Results_Fold_changes_for_2014_ARF_SM.med')
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Med Corpcor to SM, ARF 2014") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b3

new.mdf1 <- subset(melted.df2, variable=='PCC_Results_Fold_changes_for_2015_circadian_SM.max')
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("max PCC to SM, circadien") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b4

new.mdf1 <- subset(melted.df2, variable=='X2013_TYLCV.txt_median_x')
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Median TYLCV") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b5

new.mdf1 <- subset(melted.df2, variable=='X2015_root.txt_median_y')
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + ggtitle("Median Root") +
  stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("cyan3","dodgerblue2","darkorchid","indianred2"))+ 
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b6

nb=paste(c(basename(nc),"expr_violin3.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, cols=2)
dev.off()

## 3 class
nc = "Slyc_combinedclass_SMvsPM_Tcyc-BM.txt-Sl.allgenes.v2.5_generate_matrix.txt"
p1 = read.table(nc, header=T, sep="\t",row.names=1)

#subset data
newdata <- subset(p1, Class == 'SM' | Class == 'PM'| Class == 'SM-PM') 
#transform to maintain order 
df<- transform(newdata, Class=factor(Class, levels = c('SM','PM','SM-PM')))

#Convert to "long-form"
melted.df2<-melt(df, id.vars='Class')
#get variable you want to make graph of
new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Mpolymorpha_320_v3.1__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

#plot
p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("M. polymorpha med Ka/Ks")
p.b1

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Ppatens_318_v3.3__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("P.patens med Ka/Ks")
p.b2

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Smoellendorffii_91_v1.0__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("S. moellendorffii med Ka/Ks")
p.b3

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Pabies_HCMC__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("P. abies med Ka/Ks")
p.b4

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Atrichopoda_291_v1.0__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("A. trichopoda med Ka/Ks")
p.b5

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Osativa_323_v7.0__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("O. sativa med Ka/Ks")
p.b6


nb=paste("SMvsPMvsSM-PM_kaks_Mpoly-Osat_violin.pdf")
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, cols=2)
dev.off()

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Acoerulea_322_v3.1__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

#plot
p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("A. coerulea med Ka/Ks")
p.b1

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Crubella_183_v1.0__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("C. rubella med Ka/Ks")
p.b2

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Athaliana_167_TAIR10__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("A. thaliana med Ka/Ks")
p.b3

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Alyrata_384_v2.1__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("A. lyrata med Ka/Ks")
p.b4

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__BrapaFPsc_277_v1.3__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("B. rapa med Ka/Ks")
p.b5

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Tcacao_233_v1.1__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("T. cacao med Ka/Ks")
p.b6


nb=paste("SMvsPMvsSM-PM_kaks_Acoe-Tcac_violin.pdf")
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, cols=2)
dev.off()

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Nicotiana_tabacum.TN90_AYMY.SS__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

#plot
p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("N. tabacum med Ka/Ks")
p.b1

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Nicotiana_tomen__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("N. tomentosiformis med Ka/Ks")
p.b2

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Solanum_melongena__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("S. melongena med Ka/Ks")
p.b3

new.mdf1 <- subset(melted.df2, variable=='tomato_paralogs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("S. lycopersicum med Ka/Ks")
p.b4

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Solanum_pennellii__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("S. pennellii med Ka/Ks")+scale_y_continuous(limits= c(0,1.5))
p.b5

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Solanum_tuberosum__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("S. tuberosum med Ka/Ks")
p.b6
new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Capsicum_annuum_Pepper.v.1.55__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b7 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("C. annuum med Ka/Ks")
p.b7

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Capsicum_annuum_var.glabr_Chil__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b8 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("C. annuum var.Chil med Ka/Ks")
p.b8


nb=paste("SMvsPMvsSM-PM_kaks_Ntab-Can_violin.pdf")
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, p.b7, p.b8, cols=2)
dev.off()

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Ptrichocarpa_210_v3.0__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

#plot
p.b1 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("P. trichocarpa med Ka/Ks")
p.b1

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Mtruncatula_285_Mt4.0v1__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b2 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("M. truncatula med Ka/Ks")
p.b2

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Vvinifera_145__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b3 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("V. vinifera med Ka/Ks")
p.b3

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Coffea_canephora__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b4 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("C. canephora med Ka/Ks")
p.b4

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Ipomoea_trifida_ITR_r1.0__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b5 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("I. trifida med Ka/Ks")
p.b5

new.mdf1 <- subset(melted.df2, variable=='Slycopersicum__v__Petunia_axillaris__cross_sp_orthologs_medKaKs')
#omit NAs
new.mdf1 <- na.omit(new.mdf1)

p.b6 <-ggplot(new.mdf1, aes(x=factor(Class), y=value, fill=factor(Class))) + 
  geom_violin(stat = "ydensity",position = "dodge", draw_quantiles = TRUE, trim = TRUE,
              scale = "area", show.legend = NA, inherit.aes = TRUE) + stat_summary(fun.y=median, fun.ymin = median, fun.ymax = median, geom="crossbar",width = 0.5, size= 0.3, color="black")+
  scale_fill_manual(values =c("indianred2","cyan3","blueviolet"))+ ggtitle("P. axillaris med Ka/Ks")
p.b6


nb=paste("SMvsPMvsSM-PM_kaks_Ptr-Pax_violin.pdf")
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, cols=2)
dev.off()
