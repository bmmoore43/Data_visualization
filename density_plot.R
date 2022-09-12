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

new.mdf <- subset(newdata, type1 == 'basal_angiosperm')
new.mdf1 <- subset(new.mdf, type2 == "eudicot-canonical"| type2 == 'eudicot-noncanonical'| type2 == 'monocot-canonical'| type2 == 'monocot-noncanonical')
new.mdf2 <- subset(new.mdf, type3 == "TyrA1"| type3 == 'TyrA2'| type3 == 'TyrA3')

#Convert to "long-form"
#melted.df2<-melt(newdata, id.vars='feat_type')
#omit NAs
#new.mdf <- na.omit(melted.df2)
p.b2 <-ggplot(new.mdf1, aes(x=ka.ks, fill =factor(type2))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Basal Angiosperm") +
  scale_fill_manual(values =c("cyan3","blueviolet", "green","darkgoldenrod1"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b2

p.b2 <-ggplot(new.mdf2, aes(x=ka.ks, fill =factor(type3))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Basal Angiosperm") +
  scale_fill_manual(values =c("firebrick1","firebrick4", "dimgrey"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b2

new.mdf <- subset(newdata, type1 == 'basal_angiosperm-noncanonical')
new.mdf1 <- subset(new.mdf, type2 == "eudicot-canonical"| type2 == 'eudicot-noncanonical'| type2 == 'monocot-canonical'| type2 == 'monocot-noncanonical')
new.mdf2 <- subset(new.mdf, type3 == "TyrA1"| type3 == 'TyrA2'| type3 == 'TyrA3')

#Convert to "long-form"
#melted.df2<-melt(newdata, id.vars='feat_type')
#omit NAs
#new.mdf <- na.omit(melted.df2)
p.b3 <-ggplot(new.mdf1, aes(x=ka.ks, fill =factor(type2))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Basal Angiosperm Noncanonical") +
  scale_fill_manual(values =c("cyan3","blueviolet", "green","darkgoldenrod1"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b3

p.b3 <-ggplot(new.mdf2, aes(x=ka.ks, fill =factor(type3))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Basal Angiosperm Noncanonical") +
  scale_fill_manual(values =c("firebrick1","firebrick4", "dimgrey"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b3

new.mdf <- subset(newdata, type1 == 'gymnosperm') 
new.mdf1 <- subset(new.mdf, type2 == "eudicot-canonical"| type2 == 'eudicot-noncanonical'| type2 == 'monocot-canonical'| type2 == 'monocot-noncanonical')
new.mdf2 <- subset(new.mdf, type3 == "TyrA1"| type3 == 'TyrA2'| type3 == 'TyrA3')

#Convert to "long-form"
#melted.df2<-melt(newdata, id.vars='feat_type')
#omit NAs
#new.mdf <- na.omit(melted.df2)
p.b4 <-ggplot(new.mdf1, aes(x=ka.ks, fill =factor(type2))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Gymnosperm") +
  scale_fill_manual(values =c("cyan3","blueviolet", "green","darkgoldenrod1"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b4

p.b4 <-ggplot(new.mdf2, aes(x=ka.ks, fill =factor(type3))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Gymnosperm") +
  scale_fill_manual(values =c("firebrick1","firebrick4", "dimgrey"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b4

new.mdf <- subset(newdata, type.merged == 'dicot_canonical')
#Convert to "long-form"
#melted.df2<-melt(newdata, id.vars='feat_type')
#omit NAs
#new.mdf <- na.omit(melted.df2)
p.b4 <-ggplot(new.mdf, aes(x=ka.ks, fill =factor(brachy_gene))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Dicot Canonical") +
  scale_fill_manual(values =c("cyan3","blueviolet", "indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b4

new.mdf <- subset(newdata, type.merged == 'dicot_non-canonical')
#Convert to "long-form"
#melted.df2<-melt(newdata, id.vars='feat_type')
#omit NAs
#new.mdf <- na.omit(melted.df2)
p.b5 <-ggplot(new.mdf, aes(x=ka.ks, fill =factor(brachy_gene))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Dicot Non-canonical") +
  scale_fill_manual(values =c("cyan3","blueviolet", "indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b5

new.mdf <- subset(newdata, type.merged == 'gymnosperm')

p.b6 <-ggplot(new.mdf, aes(x=ka.ks, fill =factor(brachy_gene))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Gymnosperm") +
  scale_fill_manual(values =c("cyan3","blueviolet", "indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b6

new.mdf <- subset(newdata, type.merged == 'monocot-grass_canonical')

p.b7 <-ggplot(new.mdf, aes(x=ka.ks, fill =factor(brachy_gene))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Monocot-grass canonical") +
  scale_fill_manual(values =c("cyan3","blueviolet", "indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b7

new.mdf <- subset(newdata, type.merged == 'monocot-grass_non-canonical')

p.b8 <-ggplot(new.mdf, aes(x=ka.ks, fill =factor(brachy_gene))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Monocot-grass Non-canonical") +
  scale_fill_manual(values =c("cyan3","blueviolet", "indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b8

new.mdf <- subset(newdata, type.merged == 'monocot-nongrass_canonical')

p.b9 <-ggplot(new.mdf, aes(x=ka.ks, fill =factor(brachy_gene))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Monocot-Nongrass canonical") +
  scale_fill_manual(values =c("cyan3","blueviolet", "indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b9

new.mdf <- subset(newdata, type.merged == 'monocot-nongrass_non-canonical')

p.b10 <-ggplot(new.mdf, aes(x=ka.ks, fill =factor(brachy_gene))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Monocot-Nongrass Noncanonical") +
  scale_fill_manual(values =c("cyan3","blueviolet", "indianred2"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p.b10

nb="KaKs_TyrA_density.pdf"
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, cols=2)
dev.off()
#prediction density plots
new.mdf1 <- subset(melted.df2, variable=='functional_likelihood_mod')
new.mdf1 <- na.omit(new.mdf1)
p.b1 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Functional likelihood") + xlab("FL score") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b1

new.mdf1 <- subset(melted.df2, variable=='pi.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)
p.b2 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Nucleotide diversity") + xlab("nt diversity") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b2


new.mdf1 <- subset(melted.df2, variable=='sequence_similarity_to_paralogs_max_perID.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)
p.b3 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Max Percent ID") +
  xlab("maximum BLAST percent ID") + scale_x_continuous(limits=c(0,100)) +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b3

new.mdf1 <- subset(melted.df2, variable=='gene_family_size_orthomcl_1.5.continuous.mod')
new.mdf1 <- na.omit(new.mdf1)
p.b4 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Gene Family Size") +
  xlab("log number of paralogs") + scale_x_continuous() +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b4
nb= "continous_fl_pi_percentID-genefam_predicts_density.pdf"
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, cols=2)
dev.off()

##kaks
new.mdf1 <- subset(melted.df2, variable=='vvin_kaks.insig.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)
p.b1 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Vitis vinifera dN/dS") + xlab("dN/dS") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b1

new.mdf1 <- subset(melted.df2, variable=='alyr_kaks.continuous.insig.MLD')
new.mdf1 <- na.omit(new.mdf1)
p.b2 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Arabidopsis lyrata") + xlab("dN/dS") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b2


new.mdf1 <- subset(melted.df2, variable=='ptri_kaks.insig.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)
p.b3 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Populus trichocarpa") +
  xlab("dN/dS") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b3

new.mdf1 <- subset(melted.df2, variable=='slyc_kaks.insig.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)
p.b4 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Solanum lycopersicum") +
  xlab("dN/dS") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b4

new.mdf1 <- subset(melted.df2, variable=='osat_kaks.insig.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)
p.b5 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Oryza sativa") +
  xlab("dN/dS") + 
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b5

new.mdf1 <- subset(melted.df2, variable=='within_atha_ka_ks.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)
p.b6 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Arabidopsis thaliana") +
  xlab("dN/dS") + 
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b6

new.mdf1 <- subset(melted.df2, variable=='ppat_kaks.insig.continuous.MLD')
new.mdf1 <- na.omit(new.mdf1)
p.b7 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Physcomitrella patens") +
  xlab("dN/dS") + 
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b7
nb= "kaks_predicts_density.pdf"
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, p.b7, cols=2)
dev.off()

##expression characteristics

##expression breadth
new.mdf1 <- subset(new.mdf, variable=='up')
new.mdf1 <- na.omit(new.mdf1)
p.b1 <-ggplot(new.mdf1, aes(x=value,fill =factor(combined_class))) + 
  geom_density(alpha=0.5, bw = 0.5) + ggtitle("Abiotic up-regulation") + xlab("number of up-conditions") +
  scale_fill_manual(values =c("indianred2","cyan3","darkorchid1"))#+ scale_x_continuous(limits=c(0,5))
p.b1

new.mdf1 <- subset(melted.df2, variable=='hormone_expr_breadth_uponly')
new.mdf1 <- na.omit(new.mdf1)
p.b1 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Hormone up-regulation") + xlab("number of conditions") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b1

new.mdf1 <- subset(melted.df2, variable=='hormone_expr_breadth_downonly')
new.mdf1 <- na.omit(new.mdf1)
p.b2 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Hormone down-regulation") + xlab("number of conditions") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b2

new.mdf1 <- subset(melted.df2, variable=='biotic_expr_breadth_uponly')
new.mdf1 <- na.omit(new.mdf1)
p.b3 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Biotic up-regulation") +
  xlab("number of conditions") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b3

new.mdf1 <- subset(melted.df2, variable=='biotic_expr_breadth_downonly')
new.mdf1 <- na.omit(new.mdf1)
p.b4 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Biotic down-regulation") +
  xlab("number of conditions") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b4

new.mdf1 <- subset(melted.df2, variable=='abiotic.shoot_expr_breadth_uponly_mod1')
new.mdf1 <- na.omit(new.mdf1)
p.b5 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Abiotic shoot up-regulation") +
  xlab("number of conditions") + 
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b5

new.mdf1 <- subset(melted.df2, variable=='abiotic.shoot_expr_breadth_downonly')
new.mdf1 <- na.omit(new.mdf1)
p.b6 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Abiotic shoot down-regulation") +
  xlab("number of conditions") + 
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b6

new.mdf1 <- subset(melted.df2, variable=='abiotic.root_expr_breadth_uponly')
new.mdf1 <- na.omit(new.mdf1)
p.b7 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Abiotic root up-regulation") +
  xlab("number of conditions") + 
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b7

new.mdf1 <- subset(melted.df2, variable=='abiotic.root_expr_breadth_downonly')
new.mdf1 <- na.omit(new.mdf1)
p.b8 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Abiotic root down-regulation") +
  xlab("number of conditions") + 
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b8
nb= "expression_breadth_predicts_density.pdf"
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, p.b7, p.b8, cols=2)
dev.off()

#other expression features
new.mdf1 <- subset(melted.df2, variable=='expression_atgen.logged.maximum')
new.mdf1 <- na.omit(new.mdf1)
p.b <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = "nrd0") + ggtitle("Max Expression") +
  xlab("logged maximum expression intensity") + scale_x_continuous(limits=c(0,5)) +
  scale_fill_manual(values =c("indianred2","cyan3", "chartreuse3", "blueviolet"))
p.b

new.mdf1 <- subset(melted.df2, variable=='expression_atgen.logged.median')
new.mdf1 <- na.omit(new.mdf1)
p.b1 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Median Expression") +
  xlab(" logged median expression intensity") + scale_x_continuous(limits=c(0,5)) +
  scale_fill_manual(values =c("indianred2","cyan3", "chartreuse3", "blueviolet"))
p.b1

new.mdf1 <- subset(melted.df2, variable=='expression_atgen.breadth_99.continuous.mod')
new.mdf1 <- na.omit(new.mdf1)
p.b2 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Expression breadth") +
  xlab("expression breadth") +
  scale_fill_manual(values =c("indianred2","cyan3", "chartreuse3", "blueviolet"))
p.b2

new.mdf1 <- subset(melted.df2, variable=='expression_atgen.logged.var_MADoverMed')
new.mdf1 <- na.omit(new.mdf1)
p.b3 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Expression variation") +
  xlab("expression variation") +
  scale_fill_manual(values =c("indianred2","cyan3", "chartreuse3", "blueviolet"))
p.b3
nb= "continous_expression_predicts_density.pdf"
pdf(file = nb)
multiplot(p.b, p.b1, p.b2, p.b3, cols=2)
dev.off()

new.mdf1 <- subset(melted.df2, variable=='maxPCC_abiotic')
new.mdf1 <- na.omit(new.mdf1)
p.b <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to paralog- abiotic stress") +
  xlab("PCC") +
  scale_fill_manual(values =c("indianred2","cyan3", "chartreuse3", "blueviolet"))
p.b
new.mdf1 <- subset(melted.df2, variable=='maxPCC_biotic')
new.mdf1 <- na.omit(new.mdf1)
p.b1 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to paralog- biotic stress") +
  xlab("PCC") +
  scale_fill_manual(values =c("indianred2","cyan3", "chartreuse3", "blueviolet"))
p.b1
new.mdf1 <- subset(melted.df2, variable=='develop_SM_corr')
new.mdf1 <- na.omit(new.mdf1)
p.b2 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to SM genes-development") +
  xlab("PCC") +
  scale_fill_manual(values =c("indianred2","cyan3", "chartreuse3", "blueviolet"))
p.b2

new.mdf1 <- subset(melted.df2, variable=='abiotic.shoot_expr_breadth_uponly_mod_logged')
new.mdf1 <- na.omit(new.mdf1)
p.b3 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Abiotic-shoot Up-regulation") +
  xlab("number of conditions up-reg") +
  scale_fill_manual(values =c("indianred2","cyan3", "chartreuse3", "blueviolet"))
p.b3
nb = "imp_expression_features_predicted_density.pdf"
pdf(file = nb)
multiplot(p.b, p.b1, p.b2, p.b3, cols=2)
dev.off()

#correlation plots
new.mdf1 <- subset(melted.df2, variable=='develop_PM_corr')
new.mdf1 <- na.omit(new.mdf1)
p.b1 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to PM genes-development") +
  xlab("PCC") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b1

new.mdf1 <- subset(melted.df2, variable=='develop_SM_corr')
new.mdf1 <- na.omit(new.mdf1)
p.b2 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to SM genes-development") +
  xlab("PCC") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b2

new.mdf1 <- subset(melted.df2, variable=='abiotic_PM_corr')
new.mdf1 <- na.omit(new.mdf1)
p.b3 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to PM genes-abiotic") +
  xlab("PCC") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b3

new.mdf1 <- subset(melted.df2, variable=='abiotic_SM_corr')
new.mdf1 <- na.omit(new.mdf1)
p.b4 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to SM genes-abiotic") +
  xlab("PCC") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b4

new.mdf1 <- subset(melted.df2, variable=='biotic_PM_corr')
new.mdf1 <- na.omit(new.mdf1)
p.b5 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to PM genes-biotic") +
  xlab("PCC") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b5

new.mdf1 <- subset(melted.df2, variable=='biotic_SM_corr')
new.mdf1 <- na.omit(new.mdf1)
p.b6 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to SM genes-biotic") +
  xlab("PCC") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b6

new.mdf1 <- subset(melted.df2, variable=='hormone_PM_PCC')
new.mdf1 <- na.omit(new.mdf1)
p.b7 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to PM genes-hormone") +
  xlab("PCC") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b7

new.mdf1 <- subset(melted.df2, variable=='hormone_SM_corr')
new.mdf1 <- na.omit(new.mdf1)
p.b8 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5) + ggtitle("Max PCC to SM genes-hormone") +
  xlab("PCC") +
  scale_fill_manual(values =c("cyan3","darkorchid1", "indianred2","blue"))
p.b8
nb= "expression_corr_predicts_density.pdf"
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, p.b5, p.b6, p.b7, p.b8, cols=2)
dev.off()
##tomato expr breadth
new.mdf1 <- subset(new.mdf, variable=='up')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Abiotic Up-regulation") +
  xlab("number of up-reg conditions") + scale_x_continuous(limits=c(0,7)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b

new.mdf1 <- subset(new.mdf, variable=='down')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b1 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Abiotic Down-regulation") +
  xlab("number of up-reg conditions") + scale_x_continuous(limits=c(0,7)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b1

new.mdf1 <- subset(new.mdf, variable=='up')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b2 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Biotic Up-regulation") +
  xlab("number of up-reg conditions") + scale_x_continuous(limits=c(0,12)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b2

new.mdf1 <- subset(new.mdf, variable=='down')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b3 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Biotic Down-regulation") +
  xlab("number of down-reg conditions") + scale_x_continuous(limits=c(0,12)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b3

new.mdf1 <- subset(new.mdf, variable=='up')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b4 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Plant Dev Up-regulation") +
  xlab("number of up-reg conditions") + scale_x_continuous(limits=c(0,35)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b4

new.mdf1 <- subset(new.mdf, variable=='down')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b5 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Plant-Dev Down-regulation") +
  xlab("number of down-reg conditions") + scale_x_continuous(limits=c(0,35)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b5

new.mdf1 <- subset(new.mdf, variable=='up')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b6 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Fruit Up-regulation") +
  xlab("number of up-reg conditions") + scale_x_continuous(limits=c(0,20)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b6

new.mdf1 <- subset(new.mdf, variable=='down')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b7 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Fruit Down-regulation") +
  xlab("number of down-reg conditions") + scale_x_continuous(limits=c(0,20)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b7

new.mdf1 <- subset(new.mdf, variable=='up')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b8 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Hormone Up-regulation") +
  xlab("number of up-reg conditions") + scale_x_continuous(limits=c(0,10)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b8

new.mdf1 <- subset(new.mdf, variable=='down')
#changed bw for increased bandwidth (to avoid multiple peaks, should be at min 1)
p.b9 <-ggplot(new.mdf1, aes(x=value,fill =factor(Class))) + 
  geom_density(alpha=0.5, bw = 1) + ggtitle("Hormone Down-regulation") +
  xlab("number of down-reg conditions") + scale_x_continuous(limits=c(0,10)) +
  scale_fill_manual(values =c("cyan3","indianred2"))
p.b9

pdf(file = nb)
multiplot(p.b, p.b1, p.b2, p.b3, p.b4, p.b5, cols=2)
dev.off()



pdf(file = nb)
multiplot(p.b6, p.b7, p.b8, p.b9, cols=2)
dev.off()

##exprs corr multiplot

p.b1 <-ggplot(df2, aes(x=new.mdf.value,fill =factor(newcol))) + 
  geom_density(alpha=0.5) + 
  scale_fill_manual(values =c("cyan3","blue","blueviolet", "indianred2"))
p.b1

p.b2 <-ggplot(df2, aes(x=new.mdf.value,fill =factor(newcol))) + 
  geom_density(alpha=0.5) + 
  scale_fill_manual(values =c("cyan3","blue","blueviolet", "indianred2"))
p.b2

p.b3 <-ggplot(df2, aes(x=new.mdf.value,fill =factor(newcol))) + 
  geom_density(alpha=0.5) + 
  scale_fill_manual(values =c("cyan3","blue","blueviolet", "indianred2"))
p.b3

p.b4 <-ggplot(df2, aes(x=new.mdf.value,fill =factor(newcol))) + 
  geom_density(alpha=0.5) + 
  scale_fill_manual(values =c("cyan3","blue","blueviolet", "indianred2"))
p.b4
#multiplot
pdf(file = nb)
multiplot(p.b1, p.b2, p.b3, p.b4, cols=2)
dev.off()

nb= "FL.genefam.breadth_transformed3_predicts_density.pdf"
pdf(file = nb)
multiplot(p.b1, p.b2, p.b4, p.b5, cols=2)
dev.off()
#density plot with a cutoff line
p8 <- ggplot(airquality, aes(x = Ozone)) +
  geom_density(colour = lines, fill = fill,
               size = 1) +
  scale_x_continuous(name = "Mean ozone in\nparts per billion",
                     breaks = seq(0, 200, 25),
                     limits=c(0, 200)) +
  scale_y_continuous(name = "Density") +
  ggtitle("Density plot of mean ozone") +
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 9),
        axis.text.y=element_text(colour="black", size = 9))

fill <- "#4271AE"
line <- "#1F3552"

p8 <- p8 + geom_vline(xintercept = 75, size = 1, colour = "#FF3721",
                      linetype = "dashed")
p8
