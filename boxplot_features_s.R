library(ggplot2)

args = commandArgs(TRUE)
setwd(args[1])
nc= args[2]
#nb = paste(c(basename(nc),".boxplot.pdf"),collapse='')
nc = "combined_class-GO-aracyc_SMvsPM_allgenes.txt-colocal_SM-PMgenes_5-10.txt"
p1 = read.table(nc, header=T, sep="\t",row.names=1) #will remove NAs later
p1 <- subset(p1, Class == "SM" | Class == "PM")
dim(p1)
colnames(p1)

# Convert to "long-form".
library(reshape2)
library(ggplot2)

##### to drop outliers and rescale- plot ka/ks
df= data.frame(p1$Class, p1$vvin_kaks, p1$alyr_kaks,  p1$ptri_kaks,	p1$slyc_kaks,	p1$osat_kaks,	p1$ppat_kaks,	p1$within_atha_ka_ks)
melted.df2<-melt(df,id.vars='p1.Class')
new.mdf <- na.omit(melted.df2)
nb=paste(c(basename(nc),"_kaks_boxplot.pdf"),collapse='')

# Plot with ggplot2
pdf(file = nb)
p.b <-ggplot(new.mdf,aes(x=p1.Class,y=value)) + 
  geom_boxplot(outlier.shape=NA, aes(fill = factor(p1.Class))) +
  facet_grid(variable~.,scales='free') + scale_y_continuous(limits = c(0, 0.6)) + 
  scale_fill_manual(values=c("cyan3","indianred2")) + coord_flip()
p.b + annotate("text", x = 4, y = 25, label = "Some text")
dev.off()

#instead of pdf, can use ggsave
ggsave(filename = default_name(plot), plot = last_plot(),
       device = default_device(filename), path = NULL,
       scale = 1, width = par("din")[1],
       height = par("din")[2], units = c("in", "cm", "mm"),
       dpi = 300, ...)
##other categories  
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
   
## multiplot expression data

df = data.frame(p1$Class, p1$expression_atgen.unlogged.maximum.insig.continuous.MLD)
newdf <- na.omit(df)
p.b <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.expression_atgen.unlogged.maximum.insig.continuous.MLD), las=1) + 
  geom_boxplot(outlier.shape=NA, aes(fill = p1.Class)) + xlab("Gene Class") + ylab("expression.unlogged.max") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + scale_y_continuous(limits = c(0, 10000)) +
  scale_fill_manual(values = c("cyan3","indianred2")) + ggtitle("Expression Max")

df = data.frame(p1$Class, p1$expression_atgen.unlogged.median.continuous.MLD)
newdf <- na.omit(df)
p.b1 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.expression_atgen.unlogged.median.continuous.MLD), las=1) + 
  geom_boxplot(outlier.shape=NA, aes(fill = p1.Class)) + xlab("Gene Class") + ylab("expression.unlogged.median") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + scale_y_continuous(limits = c(0, 5000)) +
  scale_fill_manual(values = c("cyan3","indianred2")) + ggtitle("Expression Median")

df = data.frame(p1$Class, p1$expression_atgen.unlogged.var_MADoverMed.continuous.MLD)
newdf <- na.omit(df)
p.b2 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.expression_atgen.unlogged.var_MADoverMed.continuous.MLD), las=1) + 
  geom_boxplot(outlier.shape=NA, aes(fill = p1.Class)) + xlab("Gene Class") + ylab("expression.unlogged_MadoverMed") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = c("cyan3","indianred2")) + ggtitle("Expression MAD/Median")

df = data.frame(p1$Class, p1$expression_atgen.breadth_99.continuous.MLD)
newdf <- na.omit(df)
p.b3 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.expression_atgen.breadth_99.continuous.MLD), las=1) + 
  geom_boxplot(outlier.shape=NA, aes(fill = p1.Class)) + xlab("Gene Class") + ylab("expression.breadth") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + scale_y_continuous() +
  scale_fill_manual(values = c("cyan3","indianred2")) + ggtitle("Expression Breadth")

nb=paste(c(basename(nc), "_expression_boxplot.pdf"), collapse='')
pdf(file = nb)
multiplot(p.b, p.b1, p.b2, p.b3, cols=2)
dev.off()

#expression correlation1
df = data.frame(p1$Class, p1$expression_atgen.unlogged.expr_corr.insig.continuous.MLD)
newdf <- na.omit(df)
p.b <-ggplot(newdf, aes(x=factor(p1.Class), 
                         y=p1.expression_atgen.unlogged.expr_corr.insig.continuous.MLD), las=1) + 
  geom_boxplot(outlier.shape=NA, aes(fill = p1.Class)) + xlab("Gene Class") + ylab("expression.unlogged.expr_corr") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + scale_y_continuous() +
  scale_fill_manual(values = c("cyan3","indianred2")) + ggtitle("Expression Corr")

df = data.frame(p1$Class, p1$coexpression_clust_k2000.continuous.MLD)
newdf <- na.omit(df)
p.b1 <-ggplot(newdf, aes(x=factor(p1.Class), 
                         y=p1.coexpression_clust_k2000.continuous.MLD), las=1) + 
  geom_boxplot(outlier.shape=NA, aes(fill = p1.Class))+xlab("Gene Class") + ylab("coexpression_clust_k2000") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + scale_y_continuous() +
  scale_fill_manual(values = c("cyan3","indianred2")) + ggtitle("Coexpression clustering")

nb=paste(c(basename(nc), "_expr.corr_boxplot1.pdf"), collapse='')
pdf(file = nb)
multiplot(p.b, p.b1, cols=2)
dev.off()

#expression correlation2

df = data.frame(p1$Class, p1$abiotic_max.exprs.corr)
newdf <- na.omit(df)
p.b <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.abiotic_max.exprs.corr), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("PCC Abiotic Stress") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))+ ggtitle("Abiotic stress Corr") +
  scale_fill_manual(values = c("cyan3","indianred2")) + scale_y_continuous(limits = c(0, 1))

df = data.frame(p1$Class, p1$biotic__max.exprs.corr)
newdf <- na.omit(df)
p.b1 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.biotic__max.exprs.corr), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("PCC Biotic Stress") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + ggtitle("Biotic stress Corr") +
  scale_fill_manual(values = c("cyan3","indianred2")) + scale_y_continuous(limits = c(0, 1))

df = data.frame(p1$Class, p1$develop_maxPCC)
newdf <- na.omit(df)
p.b2 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.develop_maxPCC), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("PCC Develop") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + ggtitle("Development Corr") +
  scale_fill_manual(values = c("cyan3","indianred2")) + scale_y_continuous(limits = c(0, 1))

df = data.frame(p1$Class, p1$hormone_maxPCC)
newdf <- na.omit(df)
p.b3 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.hormone_maxPCC), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("PCC Hormone") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + ggtitle("Hormone Corr") +
  scale_fill_manual(values = c("cyan3","indianred2")) + scale_y_continuous(limits = c(0, 1))

nb=paste(c(basename(nc), "_expr.corr_boxplot2.pdf"), collapse='')
pdf(file = nb)
multiplot(p.b, p.b1, p.b2, p.b3, cols=2)
dev.off()

#gene interactions
df = data.frame(p1$Class, p1$interactors_AIMC.insig.continuous.MLD)
newdf <- na.omit(df)
 p.b1 <-ggplot(newdf, aes(x=factor(p1.Class), y=p1.interactors_AIMC.insig.continuous.MLD), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("protein-interactors_AIMC") +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_manual(values = c("cyan3","indianred2")) + scale_y_continuous(limits = c(0, 11))

df = data.frame(p1$Class, p1$interactions_AraNet.continuous.MLD)
newdf <- na.omit(df)
p.b2 <-ggplot(newdf, aes(x=factor(p1.Class), y=p1.interactions_AraNet.continuous.MLD), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("gene-interactions_Aranet") +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  scale_fill_manual(values = c("cyan3","indianred2")) + scale_y_continuous()

df= data.frame(p1$Class, p1$aaLength.continuous.MLD)
melted.df2<-melt(df,id.vars='p1.Class')
new.mdf <- na.omit(melted.df2)
p.b3 <-ggplot(new.mdf,aes(x=p1.Class,y=value), las=1) +
  geom_boxplot(outlier.shape=NA, aes(fill = p1.Class)) + 
  scale_y_continuous(limits = c(0, 1000)) + scale_fill_manual(values = c("cyan3","indianred2")) +
  xlab("Class") + ylab("AA length") + ggtitle("AA Length") + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5, colour= "black", size= rel(1)))

df= data.frame(p1$Class, p1$numOfDomains.continuous.MLD)
melted.df2<-melt(df,id.vars='p1.Class')
new.mdf <- na.omit(melted.df2)
p.b4 <-ggplot(new.mdf,aes(x=p1.Class,y=value), las=1) +
  geom_boxplot(outlier.shape=NA, aes(fill = p1.Class)) + scale_y_continuous(limits = c(0, 11))+
  scale_fill_manual(values = c("cyan3","indianred2")) +
  xlab("Class") + ylab("Number of Domains") + ggtitle("Num of Domains") + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5, colour= "black", size= rel(1)))


nb=paste(c(basename(nc),"_interactions_boxplot.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, cols=2)
dev.off()

#co-localization
df = data.frame(p1$Class, p1$num_co.local_genesSM_5)
newdf <- na.omit(df)
p.b1 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.num_co.local_genesSM_5), las=1) + 
  geom_boxplot(aes(fill=p1.Class)) + xlab("Gene Class") + ylab("co-localized genes 5") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
p.b1 <- p.b1 + scale_y_continuous() + scale_fill_manual(values = c("cyan3","indianred2")) +
  xlab("Class") + ylab("Co-localized genes SM 5") + ggtitle("Co-localized Genes 5") + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), text = element_text(size =14),
        axis.text.x  = element_text(angle=45, vjust=0.5, colour= "black", size= rel(1)))

df = data.frame(p1$Class, p1$num_co.local_genesSM_10)
newdf <- na.omit(df)
p.b2 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.num_co.local_genesSM_10), las=1) + 
  geom_boxplot(aes(fill=p1.Class)) + xlab("Gene Class") + ylab("co-localized genes 10") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
p.b2 <- p.b2 + scale_y_continuous()  + scale_fill_manual(values = c("cyan3","indianred2")) + 
  xlab("Class") + ylab("Co-localized genes SM 10") + ggtitle("Co-localized Genes 10") + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), text = element_text(size =14),
        axis.text.x  = element_text(angle=45, vjust=0.5, colour= "black", size= rel(1)))

df = data.frame(p1$Class, p1$num_co.local_genesPM_5)
newdf <- na.omit(df)
p.b3 <-ggplot(newdf, aes(x=factor(p1.Class), 
                         y=p1.num_co.local_genesPM_5), las=1) + 
  geom_boxplot(aes(fill=p1.Class)) + xlab("Gene Class") + ylab("co-localized genes 5") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
p.b3 <- p.b3 + scale_y_continuous() + scale_fill_manual(values = c("cyan3","indianred2")) +
  xlab("Class") + ylab("Co-localized genes PM 5") + ggtitle("Co-localized Genes 5") + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), text = element_text(size =14),
        axis.text.x  = element_text(angle=45, vjust=0.5, colour= "black", size= rel(1)))

df = data.frame(p1$Class, p1$num_co.local_genesPM_10)
newdf <- na.omit(df)
p.b4 <-ggplot(newdf, aes(x=factor(p1.Class), 
                         y=p1.num_co.local_genesPM_10), las=1) + 
  geom_boxplot(aes(fill=p1.Class)) + xlab("Gene Class") + ylab("co-localized genes 10") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
p.b4 <- p.b4 + scale_y_continuous()  + scale_fill_manual(values = c("cyan3","indianred2")) + 
  xlab("Class") + ylab("Co-localized genes PM 10") + ggtitle("Co-localized Genes 10") + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), text = element_text(size =14),
        axis.text.x  = element_text(angle=45, vjust=0.5, colour= "black", size= rel(1)))


nb=paste(c(basename(nc),"_colocal_boxplot.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, cols=2)
dev.off()

#measures of selection

df= data.frame(p1$Class, p1$within_atha_ka_ks)
melted.df2<-melt(df,id.vars='p1.Class')
new.mdf <- na.omit(melted.df2)

# Plot with ggplot2

p.b <-ggplot(new.mdf,aes(x=p1.Class,y=value)) + 
    geom_boxplot(aes(fill = factor(p1.Class))) + #outlier.shape=NA,
    scale_y_continuous(limits = c(0, 0.6)) + 
    scale_fill_manual(values=c("cyan3","indianred2")) + xlab("Gene Class") + 
    ylab("Athal.KaKs")
p.b
df = data.frame(p1$Class, p1$FayWuH.insig.continuous.MLD)
newdf <- na.omit(df)

p.b1 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.FayWuH.insig.continuous.MLD), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("FayWuH.insig") +
  scale_y_continuous()  + scale_fill_manual(values = c("cyan3","indianred2")) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
p.b1

df = data.frame(p1$Class, p1$MK.G.insig.continuous.MLD)
newdf <- na.omit(df)

p.b2 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.MK.G.insig.continuous.MLD), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("MK.G.insig") +
  scale_y_continuous()  + scale_fill_manual(values = c("cyan3","indianred2")) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
p.b2

df = data.frame(p1$Class, p1$pi.continuous.MLD)
newdf <- na.omit(df)
p.b3 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.pi.continuous.MLD), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("pi.continuous") +
  scale_y_continuous()  + scale_fill_manual(values = c("cyan3","indianred2")) + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
p.b3

nb=paste(c(basename(nc),"_selection.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b, p.b1, p.b2, p.b3, cols=2)
dev.off()

#plot ks, gene turnover, family size: put in multiplot

df= data.frame(p1$Class, p1$gene_family_size_orthomcl_1.5.continuous.MLD)
melted.df2<-melt(df,id.vars='p1.Class')
new.mdf <- na.omit(melted.df2)
new.mdf$p1.Class <- factor(new.mdf$p1.Class, levels = new.mdf$p1.Class[order(new.mdf$p1.Class)])
p.b <-ggplot(new.mdf,aes(x=p1.Class,y=value), las=1) +
  geom_boxplot(aes(fill = p1.Class)) +
  scale_y_continuous(limits = c(0, 75)) + scale_fill_manual(values = c("cyan3","indianred2")) +
  xlab("Class") + ylab("Gene Family Size") + ggtitle("Gene Family Size") + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5, colour= "black", size= rel(1)))
p.b

df = data.frame(p1$Class, p1$sequence_similarity_to_paralogs_max_perID.continuous.MLD)
newdf <- na.omit(df)

p.b1 <-ggplot(newdf, aes(x=factor(p1.Class), y=p1.sequence_similarity_to_paralogs_max_perID.continuous.MLD), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("sequence_similarity_to_paralogs_max_perID") +
  ggtitle("Paralog Sequence Similarity") + scale_y_continuous()  + scale_fill_manual(values = c("cyan3","indianred2")) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text.x  = element_text(angle=90, 
       vjust=0.5, colour= "black", size= rel(1)))
p.b1

df= data.frame(p1$Class, p1$atha_paralog_ks)
melted.df2<-melt(df,id.vars='p1.Class')
new.mdf <- na.omit(melted.df2)

p.b2 <-ggplot(new.mdf,aes(x=factor(p1.Class),y=value)) + 
  geom_boxplot(aes(fill = p1.Class)) +
  scale_y_continuous() + xlab("Gene Class") + ylab("Athal paralog Ks") +
  ggtitle("Ks distribution to paralog") + scale_fill_manual(values=c("cyan3","indianred2")) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text.x  = element_text(angle=90, 
                                  vjust=0.5, colour= "black", size= rel(1)))
p.b2

df = data.frame(p1$Class, p1$functional_likelihood)
newdf <- na.omit(df)
p.b3 <-ggplot(newdf, aes(x=factor(p1.Class), 
                        y=p1.functional_likelihood), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("functional likelihood score") +
  ggtitle("Functional likelihood") + scale_y_continuous()  + scale_fill_manual(values = c("cyan3","indianred2")) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text.x  = element_text(angle=90, 
        vjust=0.5, colour= "black", size= rel(1)))
p.b3

nb=paste(c(basename(nc),"_ks-geneturnover_boxplot.pdf"),collapse='')
pdf(file = nb)
multiplot(p.b, p.b1, p.b2, p.b3, cols=2)
dev.off()

df = data.frame(p1$Class, p1$Retention.rate)
newdf <- na.omit(df)
p.b3 <-ggplot(newdf, aes(x=factor(p1.Class), 
                         y=p1.Retention.rate), las=1) + 
  geom_boxplot(aes(fill = p1.Class)) + xlab("Gene Class") + ylab("Retention.rate") +
  ggtitle("Gene Retention") + scale_y_continuous(limits = c(0.7, 1))  + scale_fill_manual(values = c("cyan3","indianred2")) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text.x  = element_text(angle=90, 
                                                                        vjust=0.5, colour= "black", size= rel(1)))
p.b3

nb=paste(c(basename(nc),"_retRate_boxplot.pdf"),collapse='')
pdf(file = nb)
p.b3
dev.off()
