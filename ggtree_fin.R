#ggtree
## get ggtree
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
install.packages("ggrepel")
library(ggplot2)
library(ggtree)
library(ggrepel)
library(treeio)
library(ape)
# your gene tree
setwd() # set to your current directory
getwd() # get working directory
## new read in for raxml-ng
nb="C5H_clade_cds.fa.p2n.aln_TrNIG4_T1.raxml.supportFBP" # tree file name
t <- read.tree(file= nb)
# branch label
# get labels from txt file with header and gene\tlabel
nc= "C5H_clade_prot-cds_names_nodups.txt"
p2<-read.table(nc, header=T, sep='\t', row.names=1) 
p2
# group labels into a list- choose green plant or monocot
# green-plant
grp <- list(Amborellales=rownames(subset(p2, clade == 'basal-angiosperm')),
            Monocot=rownames(subset(p2, clade == 'monocot')),
            Eudicot=rownames(subset(p2, clade == 'dicot')),            
            Gymnosperm=rownames(subset(p2, clade == 'gymnosperm')),
            Seed.free.land.plants=rownames(subset(p2, clade == 'basal-nonflower')),
            Green.algae=rownames(subset(p2, clade == 'green-algae')))
#monocot
grp <- list(Amborellales=rownames(subset(p2, clade == 'Basal-angiosperm')),
            Alismatales.and.Acorales=rownames(subset(p2, clade == 'Alismatales-and-Acorales')),
            Other.Poales=rownames(subset(p2, clade == 'other_Poales')),            
            Asparagales=rownames(subset(p2, clade == 'Asparagales')),
            Poaceae=rownames(subset(p2, clade == 'Poaceae')),
            Arecaceae=rownames(subset(p2, clade == 'Arecaceae')),
            Non.grass_graminid=rownames(subset(p2, clade == 'Non-grass_graminid')),
            Musaceae=rownames(subset(p2, clade == 'Musaceae')),
            Dioscoreales.and.Pandanales=rownames(subset(p2, clade == 'Dioscoreales-and-Pandanales')))
print(grp)
length(grp)
## add clade to tree
t <- groupOTU(t,grp)
## get colors for clade
# all plant
cols2 <- c("Amborellales"="#882255","Seed.free.land.plants"="#117733","Eudicot"="#332288",
           "Gymnosperm"="#88CCEE","Monocot"="#DDCC77","Green.algae"="black","0"="grey")
# monocot
cols3 <- c("Amborellales"="#882255","Alismatales.and.Acorales"="#CE2220","Asparagales"="#F4A736","Arecaceae"="#D0B440","Musaceae"="#57A2AC",
           "Other.Poales"="springgreen2","Non.grass_graminid"="cyan2","Poaceae"="blue2","0"="grey",
           "Dioscoreales.and.Pandanales"="#521913")
length(cols3)
## draw tree- change scale_color_manual values to colors you want (ie. cols2 or cols3)
t1 <-ggtree(t,aes(color=group))+
  geom_tiplab(geom="text", offset=0.025, hjust=-.025, size=1.5)+scale_color_manual(values=cols3)+
  theme(legend.position = c(0.9, 0.75)) # "right"
t1
## add bootstrap
t2 <- t1 + geom_nodelab(aes(label=bootstrap)+
                          theme(legend.position = "right"),color ="black", hjust = 1, vjust= -0.25, size=2)
t2
## make pdf
nd=paste(c(basename(nb),"_tree.pdf"),collapse='')
nd
pdf(file = nd, width= 8.5 , height= 11)
t2
dev.off()
## other tree formats
## radial
t3<-ggtree(t,layout="daylight",aes(color=group))+
  scale_color_manual(values=cols3)+ geom_tiplab(geom="text", offset=0.025, hjust=-.05, size=1.5)+
  theme(legend.position = "right")+ geom_label2(aes(label=label, subset = !is.na(as.numeric(label))),color ="black", size=1.0)
t3
## circular
t4<-ggtree(t, layout="circular",aes(color=group))+
  geom_tiplab(geom="text", offset=0.025, hjust=-.05, size=1.5)+
  scale_color_manual(values=cols3)+
  theme(legend.position = "right")+ geom_label2(aes(label=label, subset = !is.na(as.numeric(label))),color ="black", size=1) #& as.numeric(label) > 80)
  #theme(legend.position = "right"),color ="black", hjust = 1.5, vjust= -0.75, size=2)
t4
#t4$label
## can comment out geom_tiplab to get rid of gene labels

## draw tree in pdf
nd=paste(c(basename(nb),"_radtree.pdf"),collapse='') #_notips
nd
pdf(file = nd, width= 8.5 , height= 11)
t3
dev.off()

nd=paste(c(basename(nb),"_cirtree.pdf"),collapse='')
nd
pdf(file = nd, width= 8.5 , height= 11)
t4
dev.off()

