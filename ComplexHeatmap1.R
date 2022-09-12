if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install(c("ComplexHeatmap"))

library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library("reshape2")
setwd() #choose working directory
#Read in data- takes a matrix file (column 1 = genes, header = treatment, values are numbers.. like logFC or 1,0)
nc ="phenylpropanoid_genes_v3.2.txt_grplantFeb1921_count_matrix_mod_noTFs.txt"
##read in matrix file
fcq1 <- read.table(nc, header=T, sep='\t', na.strings=c(""," ","NA"),row.names=1)#,fileEncoding="latin1")
#write table- can write table back
#write.table(fcq2, file = "Slyc_combinedclass_SMvsPMvsSMPM_Tcyc-BM_continuous_matrix_allfeat_20190222.txt.NAimputed.txtgenepair_distance.matrix_SMpaths_mod.txt", sep="\t")
## dataframe for classes
df1 <- fcq1[2:21,1]
df1<- as.data.frame(df1)
df2<- fcq1[1,2:45]
df2<- as.data.frame(df2)
## or read in file to annotate classes (optional)
df <- read.table("annotation.txt", header=T, sep='\t')
## subset data if you want
fcq2<- fcq1[2:21,2:45]
fcq2<- as.numeric(fcq2)

fcq2 <- as.numeric(fcq1[2:353,])#c(1:3,16)]
fcq2 <- fcq1[complete.cases(fcq1),]
fcq2 <- subset(fcq1, GO.category=='biological_process')

## write output if wanted
write.table(fcq2, file = "RESULTS_fin_mod.txt", sep="\t", row.names = TRUE)

## get a matrix
sorted_fcq2 <- fcq2[order(row.names(fcq2)), ]
df1 <- fcq1[,1,drop=FALSE]
df2<- fcq1[,2,drop=FALSE]
df1<- df1[order(row.names(df1)),]
df1<- as.data.frame(df1)
df2<- df2[order(row.names(df2)),]
df2<- as.data.frame(df2)
mat <- as.matrix(sapply(fcq2, as.numeric, rownames=TRUE))  
y <- rownames(fcq2)
y
rownames(mat) <- y
print(fcq2[1,1])
print(mat[1,1])

#transpose if desired
mat_t <- t(mat)

## omit NAs
#mat <- na.omit(mat)

#heatmap with row annotations and kmeans
ha = rowAnnotation(df = df, col = list(df = c("evolutionary" = "dodgerblue2", "duplication" = "springgreen4","expression_development"="springgreen3","expression_hormone"="dodgerblue4","expression_mutant"="blueviolet","expression_stress"="orange"))) #"expression_all" = "gold2","expression_circadian" = "darkgreen",

hm<- Heatmap(mat, name= "TF", km = 5, column_title = "Pathway", 
             row_title = "Features", column_title_side = "top", col = colorRamp2(c(-20, 0, 20), c("blue", "white", "red")),
             cluster_columns = TRUE, cluster_rows= TRUE, show_row_dend = TRUE, column_names_gp = gpar(fontsize = 7), show_row_names = TRUE,row_names_gp = gpar(fontsize = 7))
p.b1 <- hm + ha
nb=paste(c(basename(nc),"_heatmap_k5col.pdf"),collapse='')
pdf(file = nb)
p.b1
dev.off()

### heatmap with annotation examples
print(unlist(df1[1,])) # if annotating column
col_fun = colorRamp2(c(0, 0.6), c("white", "green"))
ha = HeatmapAnnotation(imp = df1[,1], which='row',col = list(imp = col_fun))
column_ha = HeatmapAnnotation(clade = unlist(df2[1,]))#, alg= unlist(df2[1,]))
column_ha = HeatmapAnnotation(clade = unlist(df1[1,])) 
ha2 = rowAnnotation(neg_set = df1[,1])
ha3 = rowAnnotation(promoter= df2[,1])
ha2 = rowAnnotation(TF = df$TF.fam, w0.25hup = row_anno_points(mat2[,1]), w0.5hup = row_anno_points(mat2[,2]), w1hrup= row_anno_points(mat2[,3]),
                    w3hrup= row_anno_points(mat2[,4]),w6hrup= row_anno_points(mat2[,5]),w12hrup= row_anno_points(mat2[,6]),w24hrup= row_anno_points(mat2[,7]),
                    w0.25hrdwn= row_anno_points(mat2[,8]),w0.5hrdwn= row_anno_points(mat2[,9]),w1hrdwn= row_anno_points(mat2[,10]),w24hrdwn= row_anno_points(mat2[,11]))
ha3 = rowAnnotation(TF = df$TF.fam, imp.rank = row_anno_points(mat2[,1:11], pch = 1:11, gp = gpar(col = 1:11)))
ha3 = rowAnnotation(DAP.fam = df$DAPfam, PCCtoDAP = df$PCC, cisBP.fam= df$cisBPfam, PCCtoCisBP= df$PCC.1)
ha_barplot= rowAnnotation(count = anno_barplot(as.numeric(df1[,1]), baseline = 0, gp = gpar(fill = "cyan"))) #right_annotation
column_bar = HeatmapAnnotation(count = anno_barplot(as.numeric(unlist(df1[1,]))),baseline = 0) #top annotation  
##heatmap
hm<- Heatmap(mat, name= "-log qvalue", column_title = "Cluster", 
             row_title = "Pathway", column_title_side = "bottom", col = colorRamp2(c(0, 2), c("white","darkred")),
             cluster_columns = F, cluster_rows= F, show_row_dend = F, column_names_gp = gpar(fontsize = 10), show_column_names = F, 
             show_row_names = T, row_names_gp = gpar(fontsize = 10), top_annotation = column_ha)
hm
hm<- Heatmap(mat_t, name= "Importance score", column_title = "Cluster", 
             row_title = "Kmer", column_title_side = "bottom", col = colorRamp2(c(0, 1), c("white","darkred")),
             cluster_columns = F, cluster_rows= F, show_row_dend = F, column_names_gp = gpar(fontsize = 8), show_column_names = T, 
             show_row_names = T, row_names_gp = gpar(fontsize = 5))#, top_annotation = column_ha, right_annotation = ha_barplot)
hm
# add text in cell- importance
hm<- Heatmap(mat, name= "Scaled importance score", column_title = "Cluster", 
             row_title = "pCRE", column_title_side = "top",row_title_side = "left", col = colorRamp2(c(0.75, 1), c("white","darkred")),
             cluster_columns = F, cluster_rows= F, show_row_dend = F, column_names_gp = gpar(fontsize = 8), show_column_names = T, 
             show_row_names = T, row_names_gp = gpar(fontsize = 8), #right_annotation = c(ha2,ha3), #row_split = df1[,1],
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 8))})
hm
# gene count
hm<- Heatmap(mat, name= "gene count", column_title = "Species", 
             row_title = "enzyme family", column_title_side = "bottom", col = colorRamp2(c(0,30), c("white","darkred")),
             cluster_columns = F, cluster_rows= F, show_row_dend = F, column_names_gp = gpar(fontsize = 8), show_column_names = T, 
             show_row_names = T, row_names_gp = gpar(fontsize = 8), right_annotation = ha_barplot, top_annotation = column_ha, 
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.0f", mat[i, j]), x, y, gp = gpar(fontsize = 6))})
hm
#gene change
hm<- Heatmap(mat, name= "gene change", column_title = "Species", 
             row_title = "enzyme family", column_title_side = "bottom", col = colorRamp2(c(-1,0,100), c("blue","white","darkred")),
             cluster_columns = F, cluster_rows= F, show_row_dend = F, column_names_gp = gpar(fontsize = 8), show_column_names = T, 
             show_row_names = T, row_names_gp = gpar(fontsize = 8),top_annotation = column_ha,
             cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.0f", mat[i, j]), x, y, gp = gpar(fontsize = 6))})
hm
# kmer count
hm<- Heatmap(mat, name= "kmer count", column_title = "Frequency", 
             row_title = "Kmer", column_title_side = "bottom", col = colorRamp2(c(0,12), c("white","darkred")),
             cluster_columns = F, cluster_rows= F, show_row_dend = F, column_names_gp = gpar(fontsize = 8), show_column_names = T, 
             show_row_names = T, row_names_gp = gpar(fontsize = 6), right_annotation = ha,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.0f", mat[i, j]), x, y, gp = gpar(fontsize = 6))})
hm
p.b1 <- hm + ha
p.b2 <- hm + ha2
p.b3 <- hm + ha3
nb=paste(c(basename(nc),"_BP_heatmap.pdf"),collapse='')
nb=paste(c(basename(nc),"_heatmap.pdf"),collapse='')
pdf(file = nb, width= 8.5 , height= 11)
hm
dev.off()
pdf(file = nb, width= 8.5 , height= 8.5)
hm
dev.off()
#heatmap with k=X clusters
library(dendextend)
row_dend = as.dendrogram(hclust(dist(mat)))
row_dend = color_branches(row_dend, k = 9) # `color_branches()` returns a dendrogram object
hm<- Heatmap(mat, cluster_rows = row_dend, name= "1-PCC distance", column_title = "TF family", 
        row_title = "CREs", column_title_side = "top", cluster_columns = TRUE, show_row_dend = TRUE, 
        column_names_gp = gpar(fontsize = 7), show_row_names = FALSE)
#add annotation
ha = rowAnnotation(df = df, col = list(df = c("wound_0.25hr_up" = "lightblue2","wound_0.5hr_up" = "lightblue3", "wound_1hr_up"="dodgerblue1", "wound_3hr_up"="dodgerblue3", "wound_6hr_up"="orchid2","wound_12hr_up"="darkorchid2","wound_24hr_up"="darkorchid4")))#,"expression"= "gold2","co-expression_cluster"= "yellow", "co-expression_correlation" = "orange2"))) #
p.b1 <- hm + ha
nb=paste(c(basename(nc),"_heatmap6.pdf"),collapse='')
pdf(file = nb)
p.b1
dev.off()

#heatmap cluster within group
hm<- Heatmap(mat, name= "1-PCC distance", column_title = "TF family", 
             row_title = "CREs", column_title_side = "top", cluster_columns = TRUE, show_row_dend = TRUE, 
             row_km = 7, column_names_gp = gpar(fontsize = 7), show_row_names = FALSE)

Heatmap(mat, name = "mat", row_km = 2, column_km = 3, show_parent_dend_line = FALSE)

#add annotation
ha = rowAnnotation(df = df, col = list(df = c("wound_0.25hr_up" = "lightblue2","wound_0.5hr_up" = "lightblue3", "wound_1hr_up"="dodgerblue1", "wound_3hr_up"="dodgerblue3", "wound_6hr_up"="orchid2","wound_12hr_up"="darkorchid2","wound_24hr_up"="darkorchid4")))#,"expression"= "gold2","co-expression_cluster"= "yellow", "co-expression_correlation" = "orange2"))) #
p.b1 <- hm + ha
nb=paste(c(basename(nc),"_heatmap6.pdf"),collapse='')
pdf(file = nb)
p.b1
dev.off()

#add barplots
df1 <- subset(df, class == 'SM')
x <- df1$feature
df2 <- data.frame("value" = df1$value) 
colnames(fcq1) <- x
df2 <- t(df2)
ha1 <- subset(df, class=='SM')
ha2 <- subset(df, class=='SM-PM')
ha3 <- subset(df, class=='PM')
ha_barplot= HeatmapAnnotation(barplot1 = anno_barplot(ha1$percent_pos, baseline = 0, gp = gpar(fill = "indianred1")), 
                              barplot2= anno_barplot(ha2$percent_pos, baseline = 0, gp = gpar(fill = "purple")), 
                              barplot3= anno_barplot(ha3$percent_pos, baseline = 0, gp = gpar(fill = "cyan")))
Heatmap(mat, name= "S.lycopersicum Domains", column_title = "Domains", 
        row_title = "Genes", column_title_side = "bottom", col = colorRamp2(c(0, 1), c("grey", "forestgreen")),
        cluster_columns = function(m) hclust(dist(m)), cluster_rows= TRUE, show_row_dend = TRUE, show_column_dend = FALSE, 
        show_column_names = TRUE, show_row_names = FALSE, top_annotation = ha_barplot)


