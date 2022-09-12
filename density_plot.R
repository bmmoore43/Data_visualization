#ggplot density plot
#write directory name in setwd()
setwd()
# load ggplot and reshape packages
library(reshape2)
library(ggplot2)
#get dataframe
#nc should equal your file name
nc = "PTAL_mlc_cds_TVMG4_T2.raxml_alt_1.1_dnds_pairs.txt"
p1 = read.table(nc, header=T, sep="\t")#,row.names=1)
p1 <- na.omit(p1)
#subset data
#3-class
newdata <- subset(p1, combined_class == 'SM' | combined_class == 'PM'| combined_class == 'SM-PM') 
#mispredicted genes
newdata <- subset(p1, Class == 'SM.SM' | Class == 'PM.PM' | Class == 'SM.PM' | Class == 'PM.SM') 
#by data type
newdata <- newdata[,c(1,3:4)] #abiotic
newdata <- newdata[,c(1,2,5)] #biotic
newdata <- newdata[,c(1,6,8)] #develop
newdata <- newdata[,c(1,7,9)] #hormone
newdata <- newdata[,c(1,2,3)] #up and down
newdata <- p1[,c(7,12:16)] #dN.dS PTAL monocot
newdata <- p1[,c(2,3)]
sorted_new <-newdata[order(newdata[,3],decreasing = T),]
sorted_new <- sorted_new[1:200,2:3]
#transform to maintain order 
df<- transform(sorted_new, Class=factor(Class, levels = unique(Class)))
#Convert to "long-form"
melted.df2<-melt(newdata, id.vars='type')
#omit NAs
new.mdf <- na.omit(melted.df2)
#merge class and variable
#library(tidyr)
#df = data.frame(new.mdf$Class, new.mdf$variable, new.mdf$value)
#df2 <- unite(df, newcol, c(new.mdf.Class, new.mdf.variable), remove=FALSE)

#new.mdf <- df[,1:3]
#make output file name
nb=paste(c(basename(nc),"_density.pdf"),collapse='')

#pdf for one plot
pdf(file = nb)
p.b1 <-ggplot(new.mdf, aes(x=value, fill =factor(type))) + 
  geom_density(alpha=0.5) + 
  scale_fill_manual(values =c("cyan3","blueviolet", "indianred2"))
p.b1
dev.off()

#pdf for multiple plots or use multiple plot function
pdf(file = nb)
p.b <-ggplot(new.mdf, aes(x=value, fill =factor(Class))) + 
  geom_density(alpha=0.5)+ facet_grid(variable~.) + scale_x_continuous(limits = c(0,0.65)) + #, scales="free_y"
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b
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
# nc = "Col_SH_pinwound_0.25hr_cluster_dwn.fa_fdrN_df_p0.01-DAP_seq_matrix.txt-Dnase_seq_matrix.txt_RF_imp_scaled.txt"
# p1 = read.table(nc, header=T, sep="\t",row.names=1)
# newdata <- p1[,c(1,6)]
# #Convert to "long-form"
# melted.df2<-melt(newdata, id.vars='feat_type')
# #omit NAs
# new.mdf <- na.omit(melted.df2)

#pdfs for multiple plots, then combined with multiplot function
new.mdf <- subset(newdata, clade_gene1 == 'Basal-angiosperm' | clade_gene2 == 'Basal-angiosperm') 
#
p.b1 <-ggplot(new.mdf, aes(x=dN.dS, fill =factor(type3))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Basal-angiosperm") + xlab("dN/dS") +
  scale_fill_manual(values =c("cyan3","blueviolet","darkgoldenrod1"))+
theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b1

new.mdf <- subset(newdata, clade_gene1 == 'Poaceae' | clade_gene2 == 'Poaceae') 
p.b2 <-ggplot(new.mdf, aes(x=dN.dS, fill =factor(type3))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Poaceae") +
  scale_fill_manual(values =c("blueviolet","green","darkgoldenrod1"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b2

new.mdf <- subset(newdata, clade_gene1 == 'Basal-angiosperm' | clade_gene2 == 'Basal-angiosperm') 
new.mdf1 <- subset(new.mdf, type3 == 'PAL'| type3 == 'PTAL') 

p.b3 <-ggplot(new.mdf1, aes(x=dN.dS, fill =factor(type3))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Basal-angiosperm") + xlab("dN/dS") +
  scale_fill_manual(values =c("blueviolet","darkgoldenrod1"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b3

new.mdf <- subset(newdata, clade_gene1 == 'Non-grass_graminid' | clade_gene2 == 'Non-grass_graminid') 
new.mdf1 <- subset(new.mdf, clade_gene1 == 'Poaceae' | clade_gene2 == 'Poaceae') 

new.mdf1 <- subset(new.mdf, type3 == 'PAL'| type3 == 'PTAL') 

p.b4 <-ggplot(new.mdf1, aes(x=dN.dS, fill =factor(type3))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Non-grass_graminid") + xlab("dN/dS") +
  scale_fill_manual(values =c("blueviolet","green","darkgoldenrod1"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b4


pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, cols=2)
dev.off()

