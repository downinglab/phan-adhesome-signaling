# Figure 2
# By Andrew Phan

rm(list=ls()) #clears variables
cat("\014") #clears console

library(dplyr)
library(Seurat)
library(patchwork)

# Load the combined dataset. In CellRanger /outs/raw_feature_bc_matrix or filtered_feature_bc_matrix
#Contains barcodes.tsv, features.tsv, and matrix.mtx
combined.data <- Read10X(data.dir = "~/DowningLab_Git/AQPhan/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data)
combined <- CreateSeuratObject(counts = combined.data, project = "combined3k", min.cells = 3, min.features = 200)
combined

########################### QC and selecting cells for further analysis ###########################

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#filter cells that have unique feature counts over 2,500 or less than 200
#filter cells that have >5% mitochondrial counts
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

########################### Normalizing the data ###########################

#use "LogNormalize" method by default. normalizes the feature expression measurements for each cell by the total expression,
#multiplies this by a scale factor (10,000 by default), and log-transforms the result
#Normalized values are stored in combined[["RNA"]]@data
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

#combined <- NormalizeData(combined) #same as above line, not as explicit with defaults

########################### Identification of highly variable features (feature selection) ###########################

#calculate a subset of features that exhibit high cell-to-cell variation in the dataset. (2000 features returned is default)
#(i.e, they are highly expressed in some cells, and lowly expressed in others) (https://www.nature.com/articles/nmeth.2645) (https://www.biorxiv.org/content/early/2018/11/02/460147.full.pdf)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(combined), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

########################### Scaling the data ###########################

#apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques, like PCA
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1 (This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate)
#The results of this are stored in combined[["RNA"]]@scale.data
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

#can remove unwanted sources of variation too, see https://satijalab.org/seurat/v3.2/combined3k_tutorial.html


########################### Perform linear dimensional reduction ###########################

#perform PCA on the scaled data. only the previously determined variable features are used as input, default, but can be defined using features argument
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

# Examine and visualize PCA results a few different ways
print(combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(combined, dims = 1:2, reduction = "pca")
DimPlot(combined, reduction = "pca")
DimHeatmap(combined, dims = 1, cells = 500, balanced = TRUE) #allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses
#Setting cells to a number plots the 'extreme' cells on both ends of the spectrum, which dramatically speeds plotting for large datasets
DimHeatmap(combined, dims = 1:15, cells = 500, balanced = TRUE) #displays PCs 1 through 15


########################### Determine the 'dimensionality' of the dataset ###########################

#To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a 'metafeature' that combines information across a correlated feature set. The top principal components therefore represent a robust compression of the dataset.

#randomly permute a subset of the data (1% by default) and rerun PCA, constructing a 'null distribution' of feature scores, and repeat this procedure. Identifies 'significant' PCs as those who have a strong enrichment of low p-value features

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
combined <- JackStraw(combined, num.replicate = 100)
combined <- ScoreJackStraw(combined, dims = 1:20)

#provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line) 
#'Significant' PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line)
JackStrawPlot(combined, dims = 1:15)

#alternatively, can use elbow plot to find inflection point
ElbowPlot(combined)

#NOTES: 
#-encourage users to repeat downstream analyses with a different number of PCs. often do not differ dramatically
#-advise users to err on the higher side when choosing this parameter. For example, performing downstream analyses with only 5 PCs does significantly and adversely affect results

########################### Cluster the cells ###########################

#first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity)
combined <- FindNeighbors(combined, dims = 1:10)
#next apply modularity optimization techniques such as the Louvain algorithm, to iteratively group cells together, with the goal of optimizing the standard modularity function
#resolution or granularity is typically good between 0.4-1.2
combined <- FindClusters(combined, resolution = 0.5)

combined <- FindClusters(combined, resolution = 0.44) #res of 0.43 or 0.44 makes 9 clusters. 0.44 a little cleaner


# Look at cluster IDs of the first 5 cells
head(Idents(combined), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
combined <- RunUMAP(combined, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(combined, reduction = "umap")

#to save plot
#saveRDS(combined, file = "../output/combined_tutorial.rds")


########################### Finding differentially expressed features (cluster biomarkers) ###########################
#find markers that define clusters via differential expression (identifies pos. and neg. markers of a cluster compared to all other cells by default, ident.1)
#min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups

# find all markers of cluster 1
cluster1.markers <- FindMarkers(combined, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(combined, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect)
cluster1.markers <- FindMarkers(combined, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#violin plot of probability distributions of marker expressions
VlnPlot(combined, features = c("MS4A1", "CD79A"))

#you can plot raw counts as well
VlnPlot(combined, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

################ Fibroblasts ##############
#Cluster 3
FeaturePlot(combined, features = c("CAV1", "TFPI2", "COL5A2", "COL3A1", "COL6A2", "CYTL1", "TWIST1", "NUPR1", "CDKN1A", "CDK6", "PDGFRA", "DDR2", "SFRP1"))

#Cluster 0
FeaturePlot(combined, features = c("COL5A2", "COL3A1", "CYTL1", "TWIST1", "NUPR1", "CDKN1A", "CDK6"))
FeaturePlot(combined, features = c("TFPI2"))


#####################

############# Neural ###########

#neural; BMP4: neuronal fate commitment/differentiation. Other BMPs: neuronal or dendrite development
#FeaturePlot(combined, features = c("NOTCH1", "TGFB2", "HES1", "OCLN", "SYT11", "DLX2", "AXL", "BMP4", "BMP5", "BMP6", "BMP7", "NFASC", "NGFR"))
FeaturePlot(combined, features = c("NEUROG3", "NEUROD1", "NEUROD4", "NEUROG1", "POU3F2", "NHLH1", "NHLH2")) #Neural spike
#FeaturePlot(combined, features = c("NEFL", "NEFM", "NEXN")) #transition to neuronal?

FeaturePlot(combined, features = c("CITED2", "CUL3", "HOPX", "SP3", "JUNB"))
##############

########### Pluripotency #############

#CLDN7, CRB3 epithelial markers critical for MET found in iPSC cluster

#Primed Pluripotency Upregulated https://www.nature.com/articles/s41598-018-24051-5
FeaturePlot(combined, features = c("LEFTY1", "LEFTY2", "TGFB1", "INHBA", "FURIN")) #EpiSCs upregulated TGF-beta
FeaturePlot(combined, features = c("PDGFRA")) #LINEAGE SPECIFIC
FeaturePlot(combined, features = c("FGFR1", "DUSP6")) # FGF signaling

#Ground Naive Pluripotency Upregulated https://www.nature.com/articles/s41598-018-24051-5
FeaturePlot(combined, features = c("DPPA3", "DNMT3L", "DPPA5", "ZFP42")) #DPPA5, DNMT3L (https://www.cell.com/cell-reports/pdfExtended/S2211-1247(18)32074-6) DPPA3, ZFP42 (https://www.nature.com/articles/s41598-018-24051-5)

#Naive LIF+Serum https://www.nature.com/articles/s41598-018-24051-5
FeaturePlot(combined, features = c("ZSCAN4", "ID1", "SYCP3", "UTF1"))
FeaturePlot(combined, features = c("CALD1"))

FeaturePlot(combined, features = c("NODAL")) #IN AND OUT, MORE OUTER
FeaturePlot(combined, features = c("TDGF1")) #IN AND OUT
FeaturePlot(combined, features = c("EPCAM")) #IN AND OUT
FeaturePlot(combined, features = c("LIN28A")) #IN AND OUT
FeaturePlot(combined, features = c("DPPA3")) #INNER
FeaturePlot(combined, features = c("UTF1")) #IN AND OUT, POTENTIALLY MORE INNER
FeaturePlot(combined, features = c("LEFTY2")) #OUTER, NOT AT INNERMOST iPSC
FeaturePlot(combined, features = c("DNMT3B")) #IN AND OUT, MORE INNER
FeaturePlot(combined, features = c("NANOG")) #IN AND OUT
FeaturePlot(combined, features = c("DPPA2")) #INNER
FeaturePlot(combined, features = c("LEFTY1")) #OUTER, LIKE LEFTY2
FeaturePlot(combined, features = c("TET1")) #IN AND OUT
FeaturePlot(combined, features = c("CDH1")) #IN AND OUT
FeaturePlot(combined, features = c("FOXH1")) #IN AND OUT
FeaturePlot(combined, features = c("DPPA5")) #INNER
FeaturePlot(combined, features = c("DNMT3L")) #INNER
FeaturePlot(combined, features = c("LIN28B")) #IN AND OUT
FeaturePlot(combined, features = c("TEAD4")) #IN AND OUT
FeaturePlot(combined, features = c("ITGB5")) #OUTER
FeaturePlot(combined, features = c("DNMT3A")) #IN AND OUT
FeaturePlot(combined, features = c("SALL4")) #IN AND OUT
#################

############# Neural Crest and SMCs ##############
#Neuralepithelium/Neural Crest Markers - Bottom right of Cluster 1 https://www.nature.com/articles/cr201211
FeaturePlot(combined, features = c("HES1", "OCLN", "CDH1", "NES", "SOX9", "NGFR", "HOXA5", "HOXD9"))
combined2 = AddModuleScore(combined, features = list(c("HES1", "OCLN", "NGFR", "NES", "HOXA5", "HOXD9")), name = "Neural.Crest")
FeaturePlot(combined2, features = c("Neural.Crest1"))
#Markers for CNS, should not be present in NCC https://www.nature.com/articles/srep19727
FeaturePlot(combined, features = c("HES5", "PAX6", "DACH1", "SOX1"))

#NCC markers from pluripotency to differentiation and migration
FeaturePlot(combined2, features = c("CCND1", "E2F1")) #proliferation markers in NCCs
FeaturePlot(combined2, features = c("POU5F1", "NANOG")) #medially located pluripotency NCC marker https://www.nature.com/articles/s41467-017-01561-w
#FeaturePlot(combined2, features = c("SOX2", "MSI1")) #laterally located pluripotency neural stem cell marker https://www.nature.com/articles/s41467-017-01561-w
FeaturePlot(combined2, features = c("KRT19")) #migratory NCC marker https://www.nature.com/articles/s41467-017-01561-w
FeaturePlot(combined, features = c("COL2A1", "TFAP2A")) #NCC migration/differentiation. Does not overlap with pluripotent NCCs https://www.nature.com/articles/s41467-017-01561-w

FeaturePlot(combined, features = c("AKT1", "AKT2", "ARPP19", "C1QTNF2", "CD36", "DGAT2", "DYRK2", "EPM2AIP1", "FOXO1", "GPLD1", "GPT2", "HIF1A", "HMGB1", "IGF1", "INSR")) # ARPP19 LEFT SIDE C1QTNF2 TOP CD36 FIBROBLASTS DGAT2 BOTTOM RIGHT GPLD1 BOTTLENECK/NCCS GPT2 iPSCs+NCCs HMGB1 proliferation+iPSC IGF1 MIGRATORY NCCs INSR iPSCs IRS1 SMCs
FeaturePlot(combined, features = c("IRS1")) # https://www.nature.com/articles/s41467-017-01561-w

FeaturePlot(combined, features = c("ACTA2")) 


#Markers for Trunk NCCs https://www.nature.com/articles/srep19727
FeaturePlot(combined, features = c("HOXA5", "HOXD9"))

#NCC Border Specifier Genes https://www.sciencedirect.com/science/article/pii/S0012160610002988
FeaturePlot(combined, features = c("MSX1", "ZIC1")) 

FeaturePlot(combined, features = c("COL2A1")) #NCC migration/differentiation. Does not overlap with pluripotent NCCs https://www.nature.com/articles/s41467-017-01561-w


#Smooth Muscle Cell Markers Top left of Cluster 1
FeaturePlot(combined, features = c("ACTA2", "TAGLN", "MYL9", "ID3", "ACTC1"))
combined2 = AddModuleScore(combined, features = list(c("ACTA2", "TAGLN", "MYL9", "ID3", "ACTC1")), name = "SMC")
FeaturePlot(combined2, features = c("SMC1"))

####################

########### Epithelial Cell #########
FeaturePlot(combined2, features = c("EPCAM", "KRT8", "KRT18", "FBP1", "OCLN")) #KRT8/18 https://www.cell.com/cell-reports/pdfExtended/S2211-1247(14)00705-0
FeaturePlot(combined2, features = c("EED")) #NT5E THY1 ENG 
###################

########## Trophoblast########
FeaturePlot(combined2, features = c("EED", "ID1", "ID2", "ETS2"))
FeaturePlot(combined2, features = c("CCNB1"))

###################

############# ZGA ###########

FeaturePlot(combined2, features = c("CCNB1"))
FeaturePlot(combined2, features = c("CCNB1"))

###############

#generates an expression heatmap for given cells and features. ex: top 20 markers
top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(combined, features = top10$gene) + NoLegend()

########################### Assigning cell type identity to clusters ###########################

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet", "lol")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#to save
#saveRDS(combined, file = "../output/combined3k_final.rds")

metadata <- combined@meta.data

metadata %>% 
  ggplot(aes(color=seurat_clusters, x=nCount_RNA, fill= seurat_clusters)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)


combined <- RenameIdents(object = combined, `0` = "Fibr.", `1` = "NCCs", `2` = "Prolif.Troph.", `3` = "cIMHC Neg.", `4` = "Epi.", `5` = "iPSCs",  `6` = "Troph.", `7` = "Neur.Epi.", `8` = "SMCs", `9` = "Neur.")

#function to create split violin plots with ggplot2
#https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "width", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

########################### QC and selecting cells for further analysis ###########################

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#filter cells that have unique feature counts over 2,500 or less than 200
#filter cells that have >5% mitochondrial counts
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

########################### Normalizing the data ###########################

#use "LogNormalize" method by default. normalizes the feature expression measurements for each cell by the total expression,
#multiplies this by a scale factor (10,000 by default), and log-transforms the result
#Normalized values are stored in combined[["RNA"]]@data
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

#combined <- NormalizeData(combined) #same as above line, not as explicit with defaults

########################### Identification of highly variable features (feature selection) ###########################

#calculate a subset of features that exhibit high cell-to-cell variation in the dataset. (2000 features returned is default)
#(i.e, they are highly expressed in some cells, and lowly expressed in others) (https://www.nature.com/articles/nmeth.2645) (https://www.biorxiv.org/content/early/2018/11/02/460147.full.pdf)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(combined), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

########################### Scaling the data ###########################

#apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques, like PCA
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1 (This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate)
#The results of this are stored in combined[["RNA"]]@scale.data
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

#can remove unwanted sources of variation too, see https://satijalab.org/seurat/v3.2/combined3k_tutorial.html


########################### Perform linear dimensional reduction ###########################

#perform PCA on the scaled data. only the previously determined variable features are used as input, default, but can be defined using features argument
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

# Examine and visualize PCA results a few different ways
print(combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(combined, dims = 1:2, reduction = "pca")
DimPlot(combined, reduction = "pca")
DimHeatmap(combined, dims = 1, cells = 500, balanced = TRUE) #allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses
#Setting cells to a number plots the 'extreme' cells on both ends of the spectrum, which dramatically speeds plotting for large datasets
DimHeatmap(combined, dims = 1:15, cells = 500, balanced = TRUE) #displays PCs 1 through 15


########################### Determine the 'dimensionality' of the dataset ###########################

#To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a 'metafeature' that combines information across a correlated feature set. The top principal components therefore represent a robust compression of the dataset.

#randomly permute a subset of the data (1% by default) and rerun PCA, constructing a 'null distribution' of feature scores, and repeat this procedure. Identifies 'significant' PCs as those who have a strong enrichment of low p-value features

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
combined <- JackStraw(combined, num.replicate = 100)
combined <- ScoreJackStraw(combined, dims = 1:20)

#provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line)
#'Significant' PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line)
JackStrawPlot(combined, dims = 1:15)

#alternatively, can use elbow plot to find inflection point
ElbowPlot(combined)

#NOTES:
#-encourage users to repeat downstream analyses with a different number of PCs. often do not differ dramatically
#-advise users to err on the higher side when choosing this parameter. For example, performing downstream analyses with only 5 PCs does significantly and adversely affect results

########################### Cluster the cells ###########################

#first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity)
combined <- FindNeighbors(combined, dims = 1:20)
#next apply modularity optimization techniques such as the Louvain algorithm, to iteratively group cells together, with the goal of optimizing the standard modularity function
#resolution or granularity is typically good between 0.4-1.2
combined <- FindClusters(combined, resolution = 0.44) #res of 0.43 or 0.44 makes 9 clusters. 0.44 a little cleaner

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(combined, reduction = "umap")

md <- combined@meta.data %>% as.data.table
md[, .N, by = c("Condition")]

#to save plot
#saveRDS(combined, file = "../output/combined_tutorial.rds")


########################### Finding differentially expressed features (cluster biomarkers) ###########################
#find markers that define clusters via differential expression (identifies pos. and neg. markers of a cluster compared to all other cells by default, ident.1)
#min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups

# find markers for every cluster compared to all remaining cells, report only the positive ones
plan("multiprocess", workers = 24)
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
x <- DoHeatmap(combined, features = clusterMarkers) + scale_fill_gradientn(colors = c(rev(brewer.pal(n = 11, name = "RdBu")))) + NoLegend()
x <- DoMultiBarHeatmap(combined, features = top10$gene, additional.group.by = c("Time", "Condition")) + scale_fill_gradientn(colors = c(rev(brewer.pal(n = 11, name = "RdBu"))))
show(x)
#fwrite(combined.markers, file = "/home/aqphan/DowningLab_Git/AQPhan/scRNA-seq/allClusterMarkers.txt", sep = "\t")

#SMC NCC NEURAL iPSC EPI TROPH FIBR NEUROEPI
clusterMarkers = c("ACTA2", "TAGLN", "SOX9", "NEUROD1", "CDH1", "NANOG", "KRT8", "EPCAM", "CENPF", "BIRC5", "CAV1", "PAK1", "CRABP1", "CRYAB")

#ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect)
cluster1.markers <- FindMarkers(combined, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#violin plot of probability distributions of marker expressions
VlnPlot(combined, features = c("MS4A1", "CD79A"))
VlnPlot(combined, features = c("SHROOM3")) + NoLegend()


#you can plot raw counts as well
VlnPlot(combined, features = c("SHROOM3"), slot = "counts", log = TRUE)

################ Fibroblasts ##############
#Cluster 3
FeaturePlot(combined, features = c("CAV1", "TFPI2", "COL5A2", "COL3A1", "COL6A2", "CYTL1", "TWIST1", "NUPR1", "CDKN1A", "CDK6", "PDGFRA", "DDR2", "SFRP1", "COL3A1")) #COL3A1 in fibroblasts https://www.sciencedirect.com/science/article/pii/S0021925819753605?via%3Dihub

#Cluster 0
FeaturePlot(combined, features = c("COL5A2", "COL3A1", "CYTL1", "TWIST1", "NUPR1", "CDKN1A", "CDK6", "COL3A1"))
FeaturePlot(combined, features = c("S100A4"))


#####################


#UMAP sorted by timepoint of both conditions together
#Change alpha of points: https://github.com/satijalab/seurat/issues/2835
del = DimPlot(combined, reduction = "umap", shuffle = T, seed = 1, group.by = "Time", pt.size = .1, label = F)
del[[1]]$layers[[1]]$aes_params$alpha = .6
del

#UMAP sorted by condition together
del = DimPlot(combined, reduction = "umap", group.by = "Condition", seed = 1, shuffle = T, pt.size = .1, label = F)
del[[1]]$layers[[1]]$aes_params$alpha = .4
del

#UMAP split by condition. could downsample so S3 looks similar to LacZ
del = DimPlot(combined, reduction = "umap", split.by = "Condition", pt.size = .1, label = F)
del

############# Neural ###########

#neural; BMP4: neuronal fate commitment/differentiation. Other BMPs: neuronal or dendrite development
#FeaturePlot(combined, features = c("NOTCH1", "TGFB2", "HES1", "OCLN", "SYT11", "DLX2", "AXL", "BMP4", "BMP5", "BMP6", "BMP7", "NFASC", "NGFR"))
FeaturePlot(combined, features = c("NEUROG3", "NEUROD1", "NEUROD4", "NEUROG1", "POU3F2", "NHLH1", "NHLH2")) #Neural spike
#FeaturePlot(combined, features = c("NEFL", "NEFM", "NEXN")) #transition to neuronal?

FeaturePlot(combined, features = c("CITED2", "CUL3", "HOPX", "SP3", "JUNB"))
##############

########### Pluripotency #############

#CLDN7, CRB3 epithelial markers critical for MET found in iPSC cluster

#Primed Pluripotency Upregulated https://www.nature.com/articles/s41598-018-24051-5
FeaturePlot(combined, features = c("LEFTY1", "LEFTY2", "TGFB1", "FURIN")) #EpiSCs upregulated TGF-beta
#FeaturePlot(combined, features = c("PDGFRA")) #LINEAGE SPECIFIC
FeaturePlot(combined, features = c("FGFR1", "DUSP6")) # FGF signaling

#Ground Naive Pluripotency Upregulated https://www.nature.com/articles/s41598-018-24051-5
FeaturePlot(combined, features = c("DPPA3", "ZFP42")) #DPPA3, ZFP42 (https://www.nature.com/articles/s41598-018-24051-5)

#Naive LIF+Serum https://www.nature.com/articles/s41598-018-24051-5
FeaturePlot(combined, features = c("ZSCAN4", "ID1", "SYCP3", "UTF1")) #SYCP3 ZGA gene https://www.mdpi.com/1422-0067/21/21/8170

#Naive https://www.cell.com/cell-reports/pdfExtended/S2211-1247(18)32074-6
FeaturePlot(combined, features = c("KLF17", "TBX3", "KLF5", "DPPA5", "DNMT3L", "DPPA3", "KHDC1L", "OLAH", "FAM151A", "TRIM60", "HORMAD1", "KHDC3L"))
#Transition troph./NCCs into iPSCs: KLF5, TBX3
#ZGA Tip: KLF17, KHDC1L, FAM151A, TRIM60, ALPP

#Primed https://www.cell.com/cell-reports/pdfExtended/S2211-1247(18)32074-6
FeaturePlot(combined2, features = c("ZIC2", "SFRP2", "PTPRZ1", "DUSP6")) #many are in NCC/Troph. transition populations

#Intermediate https://www.cell.com/cell-reports/pdfExtended/S2211-1247(18)32074-6
FeaturePlot(combined2, features = c("ABCG2", "CLDN4"))

#PGC colonization  https://www.pnas.org/content/116/51/25677
FeaturePlot(combined2, features = c("PDLIM1", "BAIAP2L1", "TDRD12", "SYCP3", "MAEL"))

#PGCs and iPSCs https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0004013
FeaturePlot(combined2, features = c("TBXT", "FGF8", "SOX17")) #PGC markers in ES cells

FeaturePlot(combined2, features = c("FUT4")) #PGC markers in ES cells

FeaturePlot(combined2, features = c("DND1", "PIWIL2", "TEX14", "LEFTY1")) #TEX14 tip DND1, PIWIL2 naive iPSC cluster

FeaturePlot(combined, features = c("NODAL")) #IN AND OUT, MORE OUTER
FeaturePlot(combined, features = c("TDGF1")) #IN AND OUT
FeaturePlot(combined, features = c("EPCAM")) #IN AND OUT
FeaturePlot(combined, features = c("LIN28A")) #IN AND OUT
FeaturePlot(combined, features = c("DPPA3")) #INNER
FeaturePlot(combined, features = c("UTF1")) #IN AND OUT, POTENTIALLY MORE INNER
FeaturePlot(combined, features = c("LEFTY2")) #OUTER, NOT AT INNERMOST iPSC
FeaturePlot(combined, features = c("DNMT3B")) #IN AND OUT, MORE INNER
FeaturePlot(combined, features = c("NANOG")) #IN AND OUT
FeaturePlot(combined, features = c("DPPA2")) #INNER
FeaturePlot(combined, features = c("LEFTY1")) #OUTER, LIKE LEFTY2
FeaturePlot(combined, features = c("TET1")) #IN AND OUT
FeaturePlot(combined, features = c("CDH1")) #IN AND OUT
FeaturePlot(combined, features = c("FOXH1")) #IN AND OUT
FeaturePlot(combined, features = c("DPPA5")) #INNER
FeaturePlot(combined, features = c("DNMT3L")) #INNER
FeaturePlot(combined, features = c("LIN28B")) #IN AND OUT
FeaturePlot(combined, features = c("TEAD4")) #IN AND OUT
FeaturePlot(combined, features = c("ITGB5")) #OUTER
FeaturePlot(combined, features = c("DNMT3A")) #IN AND OUT
FeaturePlot(combined, features = c("SALL4")) #IN AND OUT

combined = AddModuleScore(combined, features = list(c("NODAL", "TDGF1", "EPCAM", "LIN28A", "UTF1", "LEFTY2", "LEFTY1", "CDH1", "LIN28B", "DNMT3A", "SALL4")), name = "recreate.pluripotency")

VlnPlot(combined, features = "Pluripotency", idents = c(1,8,5), split.by = "Condition")

#################

############# Neural Crest and SMCs ##############
#Neuralepithelium/Neural Crest Markers - Bottom right of Cluster 1 https://www.nature.com/articles/cr201211
FeaturePlot(combined, features = c("HES1", "OCLN", "CDH1", "NES", "SOX9", "NGFR", "HOXA5", "HOXD9"))
combined = AddModuleScore(combined, features = list(c("HES1", "OCLN", "NGFR", "NES", "HOXA5", "HOXD9")), name = "Neural.Crest")
FeaturePlot(combined, features = c("Neural.Crest1"))
#Markers for CNS, should not be present in NCC https://www.nature.com/articles/srep19727
FeaturePlot(combined, features = c("HES5", "PAX6", "DACH1", "SOX1"))

#NCC markers from pluripotency to differentiation and migration
FeaturePlot(combined2, features = c("CCND1", "E2F1")) #proliferation markers in NCCs
FeaturePlot(combined2, features = c("POU5F1", "NANOG")) #medially located pluripotency NCC marker https://www.nature.com/articles/s41467-017-01561-w
#FeaturePlot(combined2, features = c("SOX2", "MSI1")) #laterally located pluripotency neural stem cell marker https://www.nature.com/articles/s41467-017-01561-w
FeaturePlot(combined, features = c("KRT19")) #migratory NCC marker https://www.nature.com/articles/s41467-017-01561-w
FeaturePlot(combined, features = c("COL2A1", "TFAP2A")) #NCC migration/differentiation. Does not overlap with pluripotent NCCs https://www.nature.com/articles/s41467-017-01561-w

FeaturePlot(combined, features = c("AKT1", "AKT2", "ARPP19", "C1QTNF2", "CD36", "DGAT2", "DYRK2", "EPM2AIP1", "FOXO1", "GPLD1", "GPT2", "HIF1A", "HMGB1", "IGF1", "INSR")) # ARPP19 LEFT SIDE C1QTNF2 TOP CD36 FIBROBLASTS DGAT2 BOTTOM RIGHT GPLD1 BOTTLENECK/NCCS GPT2 iPSCs+NCCs HMGB1 proliferation+iPSC IGF1 MIGRATORY NCCs INSR iPSCs IRS1 SMCs
FeaturePlot(combined, features = c("IRS1")) # https://www.nature.com/articles/s41467-017-01561-w

FeaturePlot(combined, features = c("TFAP2A"))

#Markers for Trunk NCCs https://www.nature.com/articles/srep19727
FeaturePlot(combined, features = c("HOXA5", "HOXD9"))

#NCC Border Specifier Genes https://www.sciencedirect.com/science/article/pii/S0012160610002988
FeaturePlot(combined, features = c("MSX1", "ZIC1"))

FeaturePlot(combined, features = c("COL2A1")) #NCC migration/differentiation. Does not overlap with pluripotent NCCs https://www.nature.com/articles/s41467-017-01561-w

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5811199/
#POTENCY
FeaturePlot(combined, features = c("SOX5"))
FeaturePlot(combined, features = c("MYCN")) #oncogene
FeaturePlot(combined, features = c("CDH1"))
FeaturePlot(combined, features = c("GEMIN2")) #inhibit BMP signaling required for EMT

combined = AddModuleScore(combined, features = list(c("SOX5", "MYCN", "CDH1", "GEMIN2")), name = "NCC_POTENCY")
FeaturePlot(combined, features = c("NCC_POTENCY1"), split.by = "Condition")
VlnPlot(combined, features = c("NCC_POTENCY1"), split.by = "Condition", idents = 1)

#DIFFERENTIATION/EMT
FeaturePlot(combined, features = c("BMP4"), split.by = "Condition") #chick trunk NCCs, higher in LacZ
FeaturePlot(combined, features = c("SOX9"))  #higher in LacZ
FeaturePlot(combined, features = c("SNAI2"), split.by = "Condition") #not clear which is more
FeaturePlot(combined, features = c("CDH2"), split.by = "Condition") #seems like more in LacZ, more highly expressed
FeaturePlot(combined, features = c("CDH11"), split.by = "Condition") #not a big difference https://onlinelibrary.wiley.com/doi/10.1002/dvg.23028

combined = AddModuleScore(combined, features = list(c("BMP4", "SOX9", "SNAI2", "CDH2", "CDH11")), name = "EMT")
FeaturePlot(combined, features = c("EMT1"), split.by = "Condition")

FeaturePlot(combined, features = c("NRP1")) #cell-environment communication - Neuropilin receptors interact with SEMA ligands
FeaturePlot(combined, features = c("SEMA3A"))
FeaturePlot(combined, features = c("NRP2"), split.by = "Condition") #higher in LacZ in NCCs

FeaturePlot(combined, features = c("EPHA4"), split.by = "Condition")
FeaturePlot(combined, features = c("EPHB1"), split.by = "Condition")
FeaturePlot(combined, features = c("EFNB2"), split.by = "Condition")

FeaturePlot(combined, features = c("EPHB3"), split.by = "Condition")
FeaturePlot(combined, features = c("EFNB1"), split.by = "Condition")

FeaturePlot(combined, features = c("VCAN"), split.by = "Condition") #versican, ECM for cancer and NCC migration

FeaturePlot(combined, features = c("WNT5A")) #Canonical WNT signaling induction and specification of NCCs, and differentiation, CSCs - self-renewal and migration
#FeaturePlot(combined, features = c("WNT3A"))
FeaturePlot(combined, features = c("FZD2")) #WNT receptor for migration and proliferation in cancer

#Smooth Muscle Cell Markers Top left of Cluster 1
FeaturePlot(combined, features = c("ACTA2", "TAGLN", "MYL9", "ID3", "ACTC1"))
combined = AddModuleScore(combined, features = list(c("ACTA2", "TAGLN", "MYL9", "ID3", "ACTC1")), name = "SMC")
FeaturePlot(combined, features = c("SMC1"))

####################

########### Epithelial Cell #########
FeaturePlot(combined, features = c("EPCAM", "KRT8", "KRT18", "FBP1", "OCLN")) #KRT8/18 https://www.cell.com/cell-reports/pdfExtended/S2211-1247(14)00705-0

FeaturePlot(combined, features = c("CLDN6", "CLDN7")) #transition from NCCs and trophoblasts to iPSCs? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5853091/
FeaturePlot(combined, features = c("CDH11"))

###################

########## Trophoblast########
FeaturePlot(combined2, features = c("EED", "ID1", "ID2", "ETS2"))
FeaturePlot(combined2, features = c("CCNB1"))
FeaturePlot(combined2, features = c("TOP2A", "UBE2C")) #tumorigenesis/drug resistance. link to iPSCs?
#FeaturePlot(combined2, features = c("CCNA1"))


###################

############# ZGA ###########

FeaturePlot(combined2, features = c("TRIM24", "YAP1", "NFYA")) #Major wave https://www.cell.com/developmental-cell/comments/S1534-5807(17)30602-0 NFYA #https://www.cell.com/cell/pdf/S0092-8674(16)30653-5.pdf
FeaturePlot(combined2, features = c("POU5F1", "NANOG")) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6558659/
FeaturePlot(combined2, features = c("CNOT6L"))
combined2 = AddModuleScore(combined2, features = list(c("TRIM24", "YAP1", "NFYA", "POU5F1", "NANOG")), name = "Major.ZGA")
FeaturePlot(combined2, features = c("Major.ZGA1"))

FeaturePlot(combined2, features = c("TEAD2", "TEAD3", "TEAD4")) #TEAD2-4 seem to colocalize on left side or with YAP1
FeaturePlot(combined2, features = c("CNOT6", "CHAF1A"))

#YAP1 targets and co-activators



FeaturePlot(combined2, features = c("ZSCAN4", "KDM4E", "DUXB", "ARGFX", "LEUTX", "DPRX", "TPRX1", "CCNA1", "TRIM43", "PRAMEF12", "RFPL4B", "TRIM51", "KLF17", "PRAMEF2", "SLC34A2")) # in tip of iPSCs https://www.nature.com/articles/ng.3844
FeaturePlot(combined2, features = c("TOP2A"))
#ZGA Tip: KLF17, KHDC1L, FAM151A, TRIM60, ALPP

FeaturePlot(combined2, features = c("GATA6", "SOX17", "GATA4", "OTX2")) #primitive endoderm markers, gives rise to yolk sac https://cob.silverchair-cdn.com/cob/content_public/journal/dev/145/21/10.1242_dev.167833/7/dev167833.pdf?Expires=1624495286&Signature=UijJwuzueCcjD4waX5ae5MDr7ZnA-WKiJu2-QPSE~bSjWWTWqB6A21asBCZkdCopRndVVdZN8fZwgRiM7-zp9jC~CEu6gyROOB~G-S26w56DnXw066GXTsPqTaEBWdtPZ6p4FwrUcvfEQYNiUEyGP4MzOmKvI9uKCHd-Zy38j8qyRg2~lwxO6vOBcrQ1mif3mT6nYBfu75~6UulBnoqafFTfslKNK9PENly5XWbGD8xhqkG74eXqyoBPaylk4IEGKN5bK3Z4BXmn02hl6RGWJp4IuC7wv3zNVlM36upU0zGGF24JSzWrUVK~p1xDy0Yro5O4hxxCm4trCgf2SvaULA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA

combined2 = AddModuleScore(combined, features = list(c("POU5F1", "NANOG", "ZSCAN4", "KDM4E", "DUXB", "ARGFX", "LEUTX", "DPRX", "TPRX1", "CCNA1", "TRIM43", "PRAMEF12", "RFPL4B", "TRIM51", "KLF17", "PRAMEF2", "SLC34A2")), name = "ZGA")
combined2 = AddModuleScore(combined, features = list(c("ZSCAN4", "KDM4E", "DUXB", "ARGFX", "LEUTX", "DPRX", "TPRX1", "CCNA1", "TRIM43", "PRAMEF12", "RFPL4B", "TRIM51", "KLF17", "PRAMEF2", "SLC34A2")), name = "minor.ZGA")
FeaturePlot(combined2, features = c("minor.ZGA1"))

FeaturePlot(combined2, features = c("HHLA1")) #endogenous retrovirus activation
#MHC class I reactivation in response to HERVs? Part of ZGA? In tumors/tumorigenesis?
#tumor signatures: metastatic, stemness/static, benign

#3rd ZGA wave (mid-preimplantation gene activation) https://www.sciencedirect.com/science/article/pii/S0070215316301016?via%3Dihub#bb0035
FeaturePlot(combined2, features = c("TRIM43"))

###############

####### Notch ##########
# FeaturePlot(combined, features = c("POU5F1")) #similar to OCT4 expression

#Receptors
FeaturePlot(combined, features = c("NOTCH1"))
#FeaturePlot(combined, features = c("NOTCH2"))
FeaturePlot(combined, features = c("NOTCH3"))
#FeaturePlot(combined, features = c("NOTCH4"))

#Ligands
FeaturePlot(combined, features = c("JAG1"))
FeaturePlot(combined, features = c("JAG2"))
FeaturePlot(combined, features = c("DLL1")) #neural
#FeaturePlot(combined, features = c("DLL3")) #neural
#FeaturePlot(combined, features = c("DLL4"))

#combined = AddModuleScore(combined, features = list(c("NOTCH1", "NOTCH3", "JAG1", "JAG2", "DLL1")), name = "NOTCH.SIGNALING")

combined = AddModuleScore(combined, features = list(c("NOTCH1", "NOTCH3", "JAG1", "JAG2", "DLL1", "DLL3", "NOTCH2")), name = "NOTCH.SIGNALING")

#Nodal Signaling
FeaturePlot(combined, features = c("NODAL"))
FeaturePlot(combined, features = c("ACVR2B"))
#FeaturePlot(combined, features = c("ACVRL1")) #NCC

#Nodal TFs
FeaturePlot(combined, features = c("TP53"))
FeaturePlot(combined, features = c("FOXH1"))

FeaturePlot(combined, features = c("FURIN"))

#nodal antagonists
FeaturePlot(combined, features = c("LEFTY1"))
FeaturePlot(combined, features = c("LEFTY2"))
FeaturePlot(combined, features = c("BMP7"))

FeaturePlot(combined, features = c("DACT2")) #expressed in neural crest

FeaturePlot(combined, features = c("TGIF1"))
FeaturePlot(combined, features = c("TGIF2"))

FeaturePlot(combined, features = c("PITX2")) #left-right symmetry

FeaturePlot(combined, features = c("TDGF1")) #cripto-1 gene that regulates stem cell properties and cancer progression https://www.intechopen.com/chapters/19271



combined = AddModuleScore(combined, features = list(c("ACVR2B", "TP53", "FOXH1", "LEFTY2")), name = "NODAL.SIGNALING")
FeaturePlot(combined, features = c("NODAL.SIGNALING1"))
VlnPlot(combined, features = c("NOTCH.SIGNALING1")) + NoLegend()


##########

####### Planar Cell Polarity Pathway ##########

#Non-Canonical WNT/PCP Involved in NCC Migration/EMT/Differentiation https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5811199/
FeaturePlot(combined, features = c("FZD7")) #receptor for WNT11 ligand. Required for CIL migration
FeaturePlot(combined, features = c("WNT11"))
FeaturePlot(combined, features = c("DVL1")) #also required for CIL migration
FeaturePlot(combined, features = c("PRICKLE1"), split.by = "Condition")

#PCP Pathway
FeaturePlot(combined, features = c("WNT5A"))
FeaturePlot(combined, features = c("VANGL2"))
FeaturePlot(combined, features = c("FZD1"))
FeaturePlot(combined, features = c("FZD3"))
#FeaturePlot(combined, features = c("FZD6"))

FeaturePlot(combined, features = c("DVL1"))
FeaturePlot(combined, features = c("DVL2"))
#FeaturePlot(combined, features = c("DVL3"))

FeaturePlot(combined, features = c("CELSR1"))
#FeaturePlot(combined, features = c("CELSR2"))
FeaturePlot(combined, features = c("CELSR3"))

FeaturePlot(combined, features = c("PRICKLE1"))
FeaturePlot(combined, features = c("PRICKLE2"))
#FeaturePlot(combined, features = c("PRICKLE3"))

combined = AddModuleScore(combined, features = list(c("WNT5A", "VANGL2", "FZD3", "DVL2", "CELSR1", "PRICKLE1", "FZD7", "WNT11", "DVL1")), name = "PCP.PATHWAY")
FeaturePlot(combined, features = c("PCP.PATHWAY1"))

VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "Condition", idents = c(1,8), split.plot = T)

#############

####### EPH Signaling #######
FeaturePlot(combined, features = c("CDH1")) #regulates EPHA2 expression

FeaturePlot(combined, features = c("EPHA1"))
FeaturePlot(combined, features = c("EPHA2"))
#FeaturePlot(combined, features = c("EPHA3"))
FeaturePlot(combined, features = c("EPHA4"))
#FeaturePlot(combined, features = c("EPHA5"))
#FeaturePlot(combined, features = c("EPHA6"))
#FeaturePlot(combined, features = c("EPHA7"))
#FeaturePlot(combined, features = c("EPHA8"))
#FeaturePlot(combined, features = c("EPHA10"))

#FeaturePlot(combined, features = c("EPHB3")) #NCCs
#FeaturePlot(combined, features = c("EPHB4")) #NCCs/iPSCs

#FeaturePlot(combined, features = c("EFNA1"))
#FeaturePlot(combined, features = c("EFNA2"))
#FeaturePlot(combined, features = c("EFNA3"))
FeaturePlot(combined, features = c("EFNA4"))
FeaturePlot(combined, features = c("EFNA5"))

#EPHA effectors
#https://www.tandfonline.com/doi/full/10.1080/13543784.2020.1762566
FeaturePlot(combined, features = c("CXCL12"))
FeaturePlot(combined, features = c("CASP3"))

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2731651/
FeaturePlot(combined, features = c("CLDN4"))
FeaturePlot(combined, features = c("ADAM10"))


#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3438283/
FeaturePlot(combined, features = c("OCLN"))

FeaturePlot(combined, features = c("SRC"), split.by = "Condition")

#combined = AddModuleScore(combined, features = list(c("CDH1", "EPHA1", "EPHA2", "EPHA4", "EFNA4", "EFNA5", "CXCL12", "OCLN")), name = "EPH.PATHWAY") #OLD DO NOT USE
combined = AddModuleScore(combined, features = list(c("EPHA1", "EPHA2", "EPHA4", "EFNA4", "EFNA5", "EPHA8", "EPHB2")), name = "EPH.PATHWAY")

VlnPlot(combined, features = c("EPH.PATHWAY1"), split.by = "Condition", pt.size = 0, cols = c("#152852", "#fd5e53")) + NoLegend()


##########

##### FN1 Signaling ########

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2132721/
FeaturePlot(combined, features = c("FN1"))
FeaturePlot(combined, features = c("ITGA5"))
#FeaturePlot(combined, features = c("ITGB1"))
FeaturePlot(combined, features = c("PTK2"))
#FeaturePlot(combined, features = c("TNS1"))
FeaturePlot(combined, features = c("PXN"))
#FeaturePlot(combined, features = c("RBL2"))

#TGFB https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415016/
FeaturePlot(combined, features = c("TGFB1"))
FeaturePlot(combined, features = c("TGFBR1"))

#Lung Cancer and FN https://www.nature.com/articles/6605154
FeaturePlot(combined, features = c("MMP9"))
FeaturePlot(combined, features = c("CAPN2"))

#combined = AddModuleScore(combined, features = list(c("FN1", "ITGA5", "PTK2", "PXN", "TGFB1", "MMP9", "CAPN2")), name = "FN1.SIGNALING") #OLD DO NOT USE

combined = AddModuleScore(combined, features = list(c("FN1", "CD44", "SDC1", "SDC4", "ITGA3", "ITGA4", "ITGA5", "ITGAV", "ITGB1", "ITGB8")), name = "FN1.SIGNALING") #from CellChat
VlnPlot(combined, features = c("FN1.SIGNALING1")) + NoLegend()

#########

#md <- combined@meta.data %>% as.data.table #can see metadata such as number of cells
#md[, .N, by = c("Sample")]

#generates an expression heatmap for given cells and features. ex: top 10 markers
top10 <- combined2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(combined2, features = top10$gene) + NoLegend()

VlnPlot(combined, features = "NOTCH2", split.by = "Condition", group.by = c("Time"), split.plot = T)

########################### Assigning cell type identity to clusters ###########################

new.cluster.ids <- c("Class I MHC Expressing", "Neural Crest Cells", "Highly Proliferating Trophoblasts", "Fibroblasts", "Epithelial", "iPSCs", "Trophoblasts", "Neuroepithelial", "Smooth Muscle Cells", "Neural")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(object = combined, `0` = "Fibr.", `1` = "NCCs", `2` = "Prolif.Troph.", `3` = "cIMHC Neg.", `4` = "Epi.", `5` = "iPSCs",  `6` = "Troph.", `7` = "Neur.Epi.", `8` = "SMCs", `9` = "Neur.")
#combined3 <- RenameIdents(object = combined2, `0` = "zero", `1` = "one", `2` = "two", `3` = "Class I MHC+", `4` = "four", `5` = "five",  `6` = "six", `7` = "seven", `8` = "eight", `9` = "nine")
#combined2 <- RenameIdents(object = combined2, '0' = "Class I MHC+", '1' = "Neural Crest", '2' = "Prolif. Trophoblasts", '3' = "Fibroblasts", '4' = "Epithelial", '5' = "iPSCs",  '6' = "Trophoblasts", '7' = "Neuroepithelial", '8' = "Smooth Muscle", '9' = "Neural")
#View(as.data.frame(Idents(combined2)))

DimPlot(combined, reduction = "umap", label = F, pt.size = 0.5) #+ NoLegend()

#to save
#saveRDS(combined, file = "../output/combined3k_final.rds")

metadata <- combined@meta.data

metadata %>%
  ggplot(aes(color=seurat_clusters, x=nCount_RNA, fill= seurat_clusters)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

#subset just cells from D12
D12combined = subset(x = combined, subset = Time == "D12")
FeaturePlot(combined, features = c("ESAM"), split.by = "Condition") #IN AND OUT
FeaturePlot(combined, features = c("TWIST1"))

############ Signatures For Each Cluster ###########
#according to markers unique to each cluster
combined = AddModuleScore(combined, features = list(c("ACTA2", "TAGLN", "MYL9", "ID3")), name = "SMC.SIG")
combined = AddModuleScore(combined, features = list(c("TFAP2A", "SOX9", "HES1", "HOXA5", "NES", "NGFR")), name = "NCC.SIG")
combined = AddModuleScore(combined, features = list(c("NEUROG3", "NEUROD1", "NEUROG1", "POU3F2", "NHLH1")), name = "NEURAL.SIG")
combined = AddModuleScore(combined, features = list(c("POU5F1", "NANOG", "CDH1", "DNMT3A", "UTF1", "LIN28A", "NODAL", "LEFTY2")), name = "iPSC.SIG")
combined = AddModuleScore(combined, features = list(c("KRT8", "KRT18", "FBP1", "EPCAM")), name = "EPI.SIG")
combined = AddModuleScore(combined, features = list(c("UBE2C", "CENPF", "BIRC5", "TOP2A", "EED", "CCNB1", "ID1", "ETS2")), name = "TROPH.SIG")
combined = AddModuleScore(combined, features = list(c("CAV1", "COL5A2", "COL3A1", "TWIST1")), name = "FIBR.SIG")
combined = AddModuleScore(combined, features = list(c("CRABP1", "CRYAB", "NBEAL1", "S100A6", "PMEPA1")), name = "NEUREPI.SIG")


################### Figures for manuscript ###################
#UMAP plots of timepoint in each condition
s3seurat = subset(x = combined, subset = Condition == "Shroom3") #separate conditions
lacZseurat = subset(x = combined, subset = Condition == "LacZ")

#cols_time_iPSCs = c("#005fba", "#0095c4", "#34bbe4", "#81d3c7")#, "#bde1d7") #blues
#cols_time_iPSCs = c("#006774", "#00a6a5", "#00c4c6", "#9dc1b6") #greens
cols_time_iPSCs = c("#006774", "#00a6a5", "#017abc", "#9dc1b6") #green to blue, winner

cols_iPSCs = viridis(n = 11, option = "C")
cols_iPSCs = c(cols_iPSCs[1], cols_iPSCs[9], cols_iPSCs[3], cols_iPSCs[4], cols_iPSCs[5], cols_iPSCs[10], cols_iPSCs[7], cols_iPSCs[8], cols_iPSCs[2], cols_iPSCs[5]) #iPSCs are yellow
#cols_iPSCs = c(cols_iPSCs[1], cols_iPSCs[2], cols_iPSCs[3], cols_iPSCs[4], cols_iPSCs[6], cols_iPSCs[5], cols_iPSCs[7], cols_iPSCs[8], cols_iPSCs[9], cols_iPSCs[10]) #harmonious

#UMAP by time
DimPlot(s3seurat, reduction = "umap", label = F, pt.size = 0.5, group.by = "Time", shuffle = T, cols = cols_time_iPSCs) #+ NoLegend()
DimPlot(lacZseurat, reduction = "umap", label = F, pt.size = 0.5, group.by = "Time", shuffle = T, cols = cols_time_iPSCs) #+ NoLegend()
#export as 6inx6in and scale to 47% for fig

#UMAP of clusters
DimPlot(s3seurat, reduction = "umap", label = F, pt.size = 0.5, shuffle = T, cols = cols_iPSCs) #+ NoLegend()
DimPlot(lacZseurat, reduction = "umap", label = F, pt.size = 0.5, shuffle = T, cols = cols_iPSCs) #+ NoLegend()
#export as 6inx6in and scale to 45% for fig


#Cluster Marker Heatmap
combined2 = combined

#have to rename because DoMultiHeatmap automatically reorders by alphabetical order
combined2 <- RenameIdents(object = combined, `0` = "aCIMHC+", `1` = "cNCC", `2` = "eProlif.Troph.", `3` = "bFibro.", `4` = "fEpi.", `5` = "giPSCs",  `6` = "dTroph.", `7` = "hNeuroepi.", `8` = "iSMC", `9` = "jNeuro.")
##########reorder to correct order (optional) ###########
combined2@active.ident <- relevel(combined2@active.ident, 10) #9 0 1 2 3 4 5 6 7 8
combined2@active.ident <- relevel(combined2@active.ident, 10) #8 9 0 1 2 3 4 5 6 7
combined2@active.ident <- relevel(combined2@active.ident, 10) #7 8 9 0 1 2 3 4 5 6
combined2@active.ident <- relevel(combined2@active.ident, 9) #5 7 8 9 0 1 2 3 4 6
combined2@active.ident <- relevel(combined2@active.ident, 9) #4 5 7 8 9 0 1 2 3 6
combined2@active.ident <- relevel(combined2@active.ident, 8) #2 5 7 8 9 0 1 3 6
combined2@active.ident <- relevel(combined2@active.ident, 10) #6 2 5 7 8 9 0 1 3
combined2@active.ident <- relevel(combined2@active.ident, 9) #1 6 2 5 7 8 9 0 3
combined2@active.ident <- relevel(combined2@active.ident, 10) #3 1 6 2 5 7 8 9 0
combined2@active.ident <- relevel(combined2@active.ident, 10) #0 3 1 6 2 5 7 8 9
############

combined2 = AddMetaData(combined2, col.name = "ColClusters", metadata = combined2@active.ident)

#define cluster markers per cell type
clusterMarkers = c("COL5A2", "PAK1", #FIBR
                   "NES", "SOX9", ##NCC
                   "CENPF", "TOP2A", #TROPH
                   "KRT8", "KRT18", #EPI
                   "NANOG", "CDH1", ##iPSCs
                   "CRYAB", "NBEAL1", #NEUROEPI
                   "TAGLN", "ACTA2", #SMCs
                   "NEUROG1", "NEUROD1") #NEUR

x <- DoMultiBarHeatmap(combined2, features = clusterMarkers, additional.group.by = c("Time", "Condition"), disp.min = 0, disp.max = 2) + scale_fill_viridis(option = "C") #adjust contrast with disp.min/max
show(x)


#Barplot of cell type proportions per timepoint
s3seurat = subset(x = combined, subset = Condition == "Shroom3") #separate conditions
lacZseurat = subset(x = combined, subset = Condition == "LacZ")

#seuratobj = lacZseurat #NOTE: run each condition separately
seuratobj = s3seurat

#Find the % of population in each cluster for each timepoint in a given condition
breakdown<-as.data.frame(table(seuratobj@meta.data$seurat_clusters, seuratobj@meta.data$Time)) #table functions counts occurances based on given crossfactors
breakdown$Perc = breakdown$Freq #find percent of total cells per condition in each cluster

for(i in length(unique(breakdown$Var2))){
  day = unique(breakdown$Var2)[i]
  rows = breakdown$Var2==day
  total = sum(breakdown$Freq[rows])
  breakdown$Perc[rows] = breakdown$Freq[rows]/total*100
}

#freq_lacZ = breakdown
freq_s3 = breakdown

#plot percentage stacked barplots of each condition over time
ggplot(freq_lacZ, aes(fill=Var1, y=Freq, x=Var2)) +
  geom_bar(position="fill", stat="identity") + theme_classic() + scale_fill_manual(values = cols_iPSCs)

ggplot(freq_s3, aes(fill=Var1, y=Freq, x=Var2)) +
  geom_bar(position="fill", stat="identity") + theme_classic() + scale_fill_manual(values = cols_iPSCs)

#################
#TODO: REDO VIOLIN PLOTS WITH MEAN
#################

#SHROOM3 Expression
VlnPlot(combined, features = c("SHROOM3"), split.by = "Condition", idents = c(1,8), split.plot = T, pt.size = 0) + geom_hline(yintercept = c(0.5521554, 0.3754184, 0.6533828, 0.4168875), color = c("blue", "blue", "red", "red"))
test1 = VlnPlot(combined, features = c("SHROOM3"), split.by = "Condition", split.plot = T)
test1 = test1$data
id = 0
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==1,1], y = test1[test1$split=="Shroom3" & test1$ident==1,1]) #p-value <2.2e-16
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==8,1], y = test1[test1$split=="Shroom3" & test1$ident==8,1]) #p-value <2.2e-16
mean(test1[test1$split=="LacZ" & test1$ident==1,1]) #mean LacZ NCCs = 0.5521554
mean(test1[test1$split=="Shroom3" & test1$ident==1,1]) #mean Shroom3 NCCs = 0.3754184
mean(test1[test1$split=="LacZ" & test1$ident==8,1]) #mean LacZ SMCs = 0.6533828
mean(test1[test1$split=="Shroom3" & test1$ident==8,1]) #mean Shroom3 SMCs = 0.4168875

#all groups SHROOM3 expression
VlnPlot(combined, features = c("SHROOM3"), split.by = "Condition", split.plot = T, pt.size = 0, cols = c("#152852", "#fd5e53")) +
  stat_summary(fun.y = mean, geom='point', size = 25, shape = 95) + stat_compare_means(comparisons = c("LacZ", "Shroom3"), label = "p.signif")
test1 = VlnPlot(combined, features = c("SHROOM3"), split.by = "Condition", split.plot = F, pt.size = 0)
test1 = test1$data

id = 9 #cluster number
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==id,1], y = test1[test1$split=="Shroom3" & test1$ident==id,1]) #p-value 0.3563
#0 CIMHC+: p-value = 5.171e-11
#1 NCCS: p-value < 2.2e-16
#2 Prolif.Troph.:p-value < 2.2e-16
#3 Fibr.: p-value = 0.0007575
#4 Epi.: p-value < 2.2e-16
#5 iPSCs: p-value < 2.2e-16
#6 Troph.: p-value = 1.314e-13
#7 Neur.Epi.: p-value = 0.9266
#8 SMCs: p-value < 2.2e-16
#9 Neur.: p-value = 0.08272
mean(test1[test1$split=="LacZ" & test1$ident==id,1]) #mean LacZ NCCs = 0.06109879
mean(test1[test1$split=="Shroom3" & test1$ident==id,1]) #mean LacZ NCCs = 0.06109879

#SHROOM3 expression LacZ vs S3
VlnPlot(combined, features = c("SHROOM3"), split.by = "Condition", group.by = "Condition", split.plot = T, pt.size = 0, cols = c("#152852", "#fd5e53")) + stat_summary(fun.y = mean, geom='point', size = 25, shape = 95) + stat_compare_means(comparisons = c("LacZ", "Shroom3"), label = "p.signif")
test1 = VlnPlot(combined, features = c("SHROOM3"), split.by = "Condition", group.by = "Condition", split.plot = F, pt.size = 0)
test1 = test1$data

wilcox.test(x = test1[test1$split=="LacZ",1], y = test1[test1$split=="Shroom3",1]) #p-value < 2.2e-16

#NCC Exploration
#PCP Pathway Expression
VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "Condition", idents = c(1,8), split.plot = T, pt.size = 0) + geom_hline(yintercept = c(0.06109879, 0.0240372, 0.08410178, 0.03275648), color = c("blue", "blue", "red", "red"))
test1 = VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "Condition", idents = c(1,8), split.plot = T)
test1 = test1$data
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==1,1], y = test1[test1$split=="Shroom3" & test1$ident==1,1]) #p-value <2.2e-16
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==8,1], y = test1[test1$split=="Shroom3" & test1$ident==8,1]) #p-value <2.2e-16
mean(test1[test1$split=="LacZ" & test1$ident==1,1]) #mean LacZ NCCs = 0.06109879
mean(test1[test1$split=="Shroom3" & test1$ident==1,1]) #mean Shroom3 NCCs = 0.0240372
mean(test1[test1$split=="LacZ" & test1$ident==8,1]) #mean LacZ SMCs = 0.08410178
mean(test1[test1$split=="Shroom3" & test1$ident==8,1]) #mean Shroom3 SMCs = 0.03275648

#EMT
VlnPlot(combined, features = c("EMT1"), split.by = "Condition", idents = c(1,8), split.plot = T, pt.size = 0) + geom_hline(yintercept = c(0.1531841, 0.08013826, 0.5152968, 0.3851453), color = c("blue", "blue", "red", "red"))
test1 = VlnPlot(combined, features = c("EMT1"), split.by = "Condition", split.plot = T)
test1 = test1$data
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==1,1], y = test1[test1$split=="Shroom3" & test1$ident==1,1]) #p-value <2.2e-16
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==8,1], y = test1[test1$split=="Shroom3" & test1$ident==8,1]) #p-value <2.2e-16

mean(test1[test1$split=="LacZ" & test1$ident==1,1]) #mean LacZ NCCs = 0.1531841
mean(test1[test1$split=="Shroom3" & test1$ident==1,1]) #mean Shroom3 NCCs = 0.08013826
mean(test1[test1$split=="LacZ" & test1$ident==8,1]) #mean LacZ SMCs = 0.5152968
mean(test1[test1$split=="Shroom3" & test1$ident==8,1]) #mean Shroom3 SMCs = 0.3851453

#Pluripotency
VlnPlot(combined, features = c("Pluripotency"), split.by = "Condition", split.plot = T, pt.size = 0, idents = c("NCCs","SMCs","iPSCs")) + geom_hline(yintercept = c(0.1209884, 0.1179736, 0.003666251, 0.02245533, 0.2489648, 0.2559908), color = c("blue", "blue", "red", "red", "green", "green"))
test1 = VlnPlot(combined, features = c("Pluripotency"), split.by = "Condition", split.plot = T)
test1 = test1$data
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==1,1], y = test1[test1$split=="Shroom3" & test1$ident==1,1]) #p-value 0.3563
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==8,1], y = test1[test1$split=="Shroom3" & test1$ident==8,1]) #p-value 4.051e-16
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==5,1], y = test1[test1$split=="Shroom3" & test1$ident==5,1]) #p-value 0.00004412

mean(test1[test1$split=="LacZ" & test1$ident==1,1]) #mean LacZ NCCs = 0.1209884
mean(test1[test1$split=="Shroom3" & test1$ident==1,1]) #mean Shroom3 NCCs = 0.1179736
mean(test1[test1$split=="LacZ" & test1$ident==8,1]) #mean LacZ SMCs = 0.003666251
mean(test1[test1$split=="Shroom3" & test1$ident==8,1]) #mean Shroom3 SMCs = 0.02245533
mean(test1[test1$split=="LacZ" & test1$ident==5,1]) #mean LacZ iPSCs = 0.2489648
mean(test1[test1$split=="Shroom3" & test1$ident==5,1]) #mean Shroom3 iPSCs = 0.2559908

#LvS pluripotency
VlnPlot(combined, features = c("Pluripotency"), split.by = "Condition", group.by = "Condition", split.plot = T, pt.size = 0) + geom_hline(yintercept = c(0.07791321, 0.09381119), color = c("red", "blue"))
test1 = VlnPlot(combined, features = c("Pluripotency"), split.by = "Condition", group.by = "Condition", split.plot = T, pt.size = 0)
test1 = test1$data
wilcox.test(x = test1[test1$split=="LacZ",1], y = test1[test1$split=="Shroom3",1]) #p-value < 2.2e-16

median(test1[test1$split=="LacZ",1]) #0.0956584 mean 0.07791321 median
median(test1[test1$split=="Shroom3",1]) #0.1038581 mean 0.09381119 median

#EPHA
VlnPlot(combined, features = c("EPH.PATHWAY1"), split.by = "Condition", idents = c("NCCs", "SMCs", "iPSCs"), split.plot = T, pt.size = 0, cols = c("#152852", "#fd5e53")) +
  stat_summary(fun.y = mean, geom='point', size = 25, shape = 95) + stat_compare_means(comparisons = c("LacZ", "Shroom3"), label = "p.signif")
test1 = VlnPlot(combined, features = c("EPH.PATHWAY1"), split.by = "Condition", idents = c("NCCs", "SMCs", "iPSCs"), split.plot = T, pt.size = 0, cols = c("#152852", "#fd5e53"))
test1 = test1$data
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==1,1], y = test1[test1$split=="Shroom3" & test1$ident==1,1]) #NCCs: p-value 0.0001336
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==5,1], y = test1[test1$split=="Shroom3" & test1$ident==5,1]) #iPSCs: p-value < 2.2e-16
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident=="SMCs",1], y = test1[test1$split=="Shroom3" & test1$ident=="SMCs",1]) #SMCs: p-value = 6.077e-07

mean(test1[test1$split=="LacZ" & test1$ident=="SMCs",1]) #mean LacZ NCCs = 0.02485683
mean(test1[test1$split=="Shroom3" & test1$ident=="SMCs",1]) #mean Shroom3 NCCs = 0.03116016
mean(test1[test1$split=="LacZ" & test1$ident==5,1]) #mean LacZ iPSCs = -0.02155388
mean(test1[test1$split=="Shroom3" & test1$ident==5,1]) #mean Shroom3 iPSCs = 0.0031577

#FN1
VlnPlot(combined, features = c("FN1.SIGNALING1"), split.by = "Condition", idents = c("NCCs", "SMCs", "iPSCs"), split.plot = T, pt.size = 0, cols = c("#152852", "#fd5e53")) +
  stat_summary(fun.y = mean, geom='point', size = 25, shape = 95) + stat_compare_means(comparisons = c("LacZ", "Shroom3"), label = "p.signif")
test1 = VlnPlot(combined, features = c("FN1.SIGNALING1"), split.by = "Condition", idents = c("NCCs", "SMCs", "iPSCs"), split.plot = T, pt.size = 0, cols = c("#152852", "#fd5e53"))
test1 = test1$data
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==1,1], y = test1[test1$split=="Shroom3" & test1$ident==1,1]) #NCCs: p-value < 2.2e-16
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==8,1], y = test1[test1$split=="Shroom3" & test1$ident==8,1]) #SMCs: p-value < 2.2e-16
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident=="iPSCs",1], y = test1[test1$split=="Shroom3" & test1$ident=="iPSCs",1]) #iPSCs: p-value < 2.2e-16

mean(test1[test1$split=="LacZ" & test1$ident==1,1]) #mean LacZ NCCs = 0.1364116
mean(test1[test1$split=="Shroom3" & test1$ident==1,1]) #mean Shroom3 NCCs = 0.07695954
mean(test1[test1$split=="LacZ" & test1$ident==8,1]) #mean LacZ SMCs = 0.4678382
mean(test1[test1$split=="Shroom3" & test1$ident==8,1]) #mean Shroom3 SMCs = 0.3378847

#NOTCH
VlnPlot(combined, features = c("NOTCH.SIGNALING1"), split.by = "Condition", idents = c("NCCs", "SMCs", "iPSCs"), split.plot = T, pt.size = 0, cols = c("#152852", "#fd5e53")) +
  stat_summary(fun.y = mean, geom='point', size = 25, shape = 95) + stat_compare_means(comparisons = c("LacZ", "Shroom3"), label = "p.signif")
test1 = VlnPlot(combined, features = c("NOTCH.SIGNALING1"), split.by = "Condition", idents = c("NCCs", "SMCs", "iPSCs"), split.plot = T, pt.size = 0, cols = c("#152852", "#fd5e53"))
test1 = test1$data
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==1,1], y = test1[test1$split=="Shroom3" & test1$ident==1,1]) #NCCs: p-value 0.001491
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident==5,1], y = test1[test1$split=="Shroom3" & test1$ident==5,1]) #iPSCs: p-value < 2.2e-16
wilcox.test(x = test1[test1$split=="LacZ" & test1$ident=="SMCs",1], y = test1[test1$split=="Shroom3" & test1$ident=="SMCs",1]) #SMCs: p-value = 1.068e-06

mean(test1[test1$split=="LacZ" & test1$ident==1,1]) #mean LacZ NCCs = 0.03129338
mean(test1[test1$split=="Shroom3" & test1$ident==1,1]) #mean Shroom3 NCCs = 0.01782587
mean(test1[test1$split=="LacZ" & test1$ident==5,1]) #mean LacZ SMCs = -0.0747546
mean(test1[test1$split=="Shroom3" & test1$ident==5,1]) #mean Shroom3 SMCs = -0.04494451

#Plot Gene Signatures Over Pseudotime NOTE: AFTER YOU RUN MONOCLE3 AND ADD PSEUDOTIME VALUES USING AddMetaData()
#separate conditions
s3seurat = subset(x = combined, subset = Condition == "Shroom3")
lacZseurat = subset(x = combined, subset = Condition == "LacZ")

seurat_obj = lacZseurat

df = as.data.frame(seurat_obj@meta.data$Condition)
colnames(df) = "Condition"
df$Time = seurat_obj@meta.data$Pseudotime
df$Cluster = seurat_obj@active.ident
df$NOTCH.SIG = seurat_obj@meta.data$NOTCH.SIGNALING1
df$FN.SIG = seurat_obj@meta.data$FN1.SIGNALING1
df$EPHA.SIG = seurat_obj@meta.data$EPH.PATHWAY1
df$PCP.SIG = seurat_obj@meta.data$PCP.PATHWAY1

# df$iPSC.SIG = seurat_obj@meta.data$iPSC.SIG1
# df$SMC.SIG = seurat_obj@meta.data$SMC.SIG1
# df$NCC.SIG = seurat_obj@meta.data$NCC.SIG1
# df$Neural.SIG = seurat_obj@meta.data$NEURAL.SIG1
# df$EPI.SIG = seurat_obj@meta.data$EPI.SIG1
# df$TROPH.SIG = seurat_obj@meta.data$TROPH.SIG1
# df$FIBR.SIG = seurat_obj@meta.data$FIBR.SIG1
# df$NEUREPI.SIG = seurat_obj@meta.data$NEUREPI.SIG1

df <- df[order(df$Time),] #order based on pseudotime

seurat_obj = s3seurat

df2 = as.data.frame(seurat_obj@meta.data$Condition)
colnames(df2) = "Condition"
df2$Time = seurat_obj@meta.data$Pseudotime
df2$Cluster = seurat_obj@active.ident
df2$NOTCH.SIG = seurat_obj@meta.data$NOTCH.SIGNALING1
df2$FN.SIG = seurat_obj@meta.data$FN1.SIGNALING1
df2$EPHA.SIG = seurat_obj@meta.data$EPH.PATHWAY1
df2$PCP.SIG = seurat_obj@meta.data$PCP.PATHWAY1

df2 <- df2[order(df2$Time),] #order based on pseudotime

dfcomb <- rbind(df, df2) #combine the conditions

ggplot(dfcomb, aes(x = Time, y = PCP.SIG, col = Condition)) + geom_split_violin(scale = "area")  #width looks nicer than area (default)

ggplot(dfcomb, aes(x = Time, y = PCP.SIG, col = Cluster)) + geom_violin(scale = "width") + facet_grid(Condition ~ Time)
#Notch: D12 S3 SMC much lower than LacZ

wilcox.test(df$NOTCH.SIG[df$Time=="D12"], df2$NOTCH.SIG[df2$Time=="D12"]) # p-value < 2.2e-16 S3<LacZ
wilcox.test(df$FN.SIG[df$Time=="D12"], df2$FN.SIG[df2$Time=="D12"]) # p-value = 0.0326 8 S3<LacZ
wilcox.test(df$EPHA.SIG[df$Time=="D12"], df2$EPHA.SIG[df2$Time=="D12"]) # p-value < 2.2e-16

t.test(unlist(subset(df, Time=="D12" & (Cluster==8 ), select = c(NOTCH.SIG))), unlist(subset(df2, Time=="D12" & (Cluster==8 ), select = c(NOTCH.SIG))))

test = wilcox.test(unlist(subset(df, Time=="D12" & (Cluster==8 ), select = c(NOTCH.SIG))), unlist(subset(df2, Time=="D12" & (Cluster==8 ), select = c(NOTCH.SIG))))

mean(unlist(subset(df2, Time=="D12" & (Cluster==8 ), select = c(NOTCH.SIG))))
test = df$NOTCH.SIG[df$Time=="D12" & df$Cluster==8]

#plot a single gene signature score over pseudotime with a loess
#Input: dataframe with pseudotime and gene signature values; the column name for the x axis; the column name for the y axis; the number of breaks for the bins on the x axis [NOTE: TOO MANY BREAKS WILL PREVENT LOESS LINE FROM PLOTTING CORRECTLY]; a loess span value [Default=0.2]
#Output: A plot
plot_GeneSigScore_Time <- function(df, x = "Time", y, span=0.2, ylim = NA){

  dfagg=aggregate(df, by=list(df[[x]]), mean) #average value per break

  lw1 <- loess(df[[y]] ~ df[[x]], data=df, span = span)

  if(!is.na(ylim)) #if ylim given, plot using limits
    plot(dfagg[[y]] ~ Group.1, data=dfagg,pch=19,cex=0.1, ylim = ylim)
  else
    plot(dfagg[[y]] ~ Group.1, data=dfagg,pch=19,cex=0.1)
  # pl <- ggplot(data = dfagg, aes(x=dfagg[[y]],y=dfagg[[y]])) + geom_point()
  # print(pl)
  j <- order(df[[x]])
  lines(df[[x]][j],lw1$fitted[j],col="red",lwd=3)
}

# breakdown$perc.diff = breakdown$perc #find differences between percent compositions between conditions
# for(r in 1:(nrow(breakdown)/2)){
#   breakdown$perc.diff[r] = breakdown$perc[r]-breakdown$perc[r+9]
#   breakdown$perc.diff[r+9] = -1*breakdown$perc.diff[r]
# }

plot_GeneSigScore_Time(df = df[!is.infinite(df$Time),], y = "PCP.SIG", span = 0.1)
par(new=T)
plot_GeneSigScore_Time(df = df2[!is.infinite(df2$Time),], y = "PCP.SIG", span = 0.1)

#lacZseurat = AddModuleScore(lacZseurat, features = list(c("WNT5A", "VANGL2", "FZD3", "DVL2", "CELSR1", "PRICKLE1", "FZD7", "WNT11", "DVL1")), name = "PCP.PATHWAY") #add PCP module scores to seurat objects

#CV Mean Analyses
combined <- RenameIdents(object = combined, `0` = "CIMHC+", `1` = "NCCs", `2` = "Prolif.Troph.", `3` = "Fibr.", `4` = "Epi.", `5` = "iPSCs",  `6` = "Troph.", `7` = "Neur.Epi.", `8` = "SMCs", `9` = "Neur.")

#function that for each gene, calculates: mean expression, CV^2 (used in Desai et al. Science paper), and Fano Factor
calcCV <- function(df) {
  mean = rowMeans(df)
  sd = apply(df, 1, sd) #find sd
  cv = (sd/mean)^2 #CV^2 = (mean/sd)^2
  fano = (sd^2)/mean #fano factor = sd^2/mean
  ans = cbind(data.frame(mean),data.frame(cv), data.frame(fano))
  ans = log10(ans) #take log of each element
  rownames(ans) = rownames(df) #retain rownames (gene names)
  return(ans)
}

#Identify potential gene signatures
sigs = list(pluri = c("RHOX5","TDGF1","UTF1","MKRN1","DPPA5A","UPP1","CHCHD10","KLF2","TRAP1A","MYLPF","1700013H16RIK","AA467197","DHX16","MT2","UBE2A","KHDC3","PYCARD","HSP90AA1","PRRC1","HAT1","CALCOCO2","IMPA2","SAA3","OOEP","BNIP3","MT1","ASNS","ALDOA","TDH","GJB3","RBPMS2","PRPS1","FAM25C","EIF2S2","CENPM","NANOG","NDUFA4L2","SYCE2","GM13251","TAF7","NUDT4","COX5A","SOD2","S100A13","FKBP6","RHOX9","GDF3","2700094K13RIK","FMR1NB","HMGN2","UBALD2","LACTB2","FOLR1","GM7325","AGTRAP","SPP1","HELLS","DPPA4","GABARAPL2","RHOX6","RHOX1","CDC5L","TEX19.1","TRIM28","ATP5G1","SOX2","JAM2","FKBP3","COX7B","ASH2L","DUT","DTYMK","GPX4","EIF4EBP1","MORC1","FABP3","ZFP428","AQP3","GRHPR","HIGD1A","RPP25","RBPMS","MMP3","APOBEC3","SPC24","XLR3A","REC114","MTF2","SNRPN","GM13580","GMNN","CHMP4C","HSF2BP","POLR2E","BLVRB","LDHB","APOC1","SYNGR1","BEX1","NR2C2AP"))#,
            #pcp = c("WNT5A","VANGL2","FZD3","DVL2","CELSR1","PRICKLE1","FZD7","WNT11","DVL1"),
            #iPSC = c("POU5F1","NANOG","CDH1","DNMT3A","UTF1","LIN28A","NODAL","LEFTY2"),
            #FN = c("FN1","CD44","SDC1","SDC4","ITGA3","ITGA4","ITGA5","ITGAV","ITGB1","ITGB8"),
            #EMT = c("SNAI2","BMP4","SOX9","CDH2","CDH11"))

#extract the normalized expression data
expr_all = as.data.frame(combined@assays$RNA@data)
expr_all$gene_sig = 0 #set default gene sig to 0

for(i in 1:length(sigs[])){ #go through each gene signature and identify genes in the signature
  expr_all$gene_sig[match(sigs[[i]], rownames(expr_all))] = names(sigs)[i]
}

#extract cell type metadata and filter for only cells from cell type
ctype = as.data.frame(combined@active.ident)
ctype$id = rownames(ctype)
celltype_ids = list(CIMHC=ctype[ctype[,1]=="CIMHC+",], NCCs=ctype[ctype[,1]=="NCCs",], ProlifTroph=ctype[ctype[,1]=="Prolif.Troph.",], Fibr=ctype[ctype[,1]=="Fibr.",], Epi=ctype[ctype[,1]=="Epi.",], iPSCs=ctype[ctype[,1]=="iPSCs",], Troph=ctype[ctype[,1]=="Troph.",], NeurEpi=ctype[ctype[,1]=="Neur.Epi.",], SMCs=ctype[ctype[,1]=="SMCs",], Neur=ctype[ctype[,1]=="Neur.",])
#CIMHC=ctype[ctype[,1]=="CIMHC+",]

#choose cell type to subset data by (optional)
#CIMHC NCCs ProlifTroph Fibr Epi iPSCs Troph NeurEpi SMCs Neur
cellsub = celltype_ids[["Neur"]] #choose cell type to subset by
expr = expr_all[,match(cellsub$id, colnames(expr_all))] #subset the cells for CV and mean calculations
expr$gene_sig = expr_all$gene_sig

#calculate CV and mean of normalized expression for each gene in the LacZ condition
lacZ_norm = grepl("_LacZ",colnames(expr))
lacZ_norm = expr[,lacZ_norm]
lacZ_cvmean = calcCV(lacZ_norm)
lacZ_cvmean$gene = rownames(expr)
lacZ_cvmean$cond = "LacZ"
lacZ_cvmean$gene_sig = expr$gene_sig

#calculate CV and mean of normalized expression for each gene in the SHROOM3 condition
S3_norm = grepl("_Shroom3",colnames(expr))
S3_norm = expr[,S3_norm]
S3_cvmean = calcCV(S3_norm)
S3_cvmean$gene = rownames(expr)
S3_cvmean$cond = "SHROOM3"
S3_cvmean$gene_sig = expr$gene_sig

#merge the two conditions to one dataframe
comb_cvmean = rbind(lacZ_cvmean, S3_cvmean)

#create plots of CV vs mean
cvL = ggplot(lacZ_cvmean, aes(x=mean, y=cv)) + geom_point() + theme_classic()
cvS = ggplot(S3_cvmean, aes(x=mean, y=cv)) + geom_point() + theme_classic()
cv_comb = ggplot(comb_cvmean[sample(nrow(comb_cvmean)),], aes(x=mean, y=cv, col=cond)) + geom_point(alpha=0.75) + theme_classic() + ylim(-3.25,3.6) + xlim(-4.5,1) #shuffles rows to randomize order of points
#cv_comb

#look at CV vs mean for only subsets of genes
#sub = pcp
sub = "pluri"
#sub = FN
#sub = EMT

cv_sub = grepl(paste(sub, collapse = "|"), comb_cvmean$gene)
cv_sub = comb_cvmean[cv_sub,]
cv_pl = ggplot(cv_sub[sample(nrow(cv_sub)),], aes(x=mean, y=cv, col=cond)) + geom_point(alpha=0.75) + theme_classic()
#cv_pl

#plot the difference of mean and CV of the same genes
dif_cv <- merge(lacZ_cvmean, S3_cvmean, by=0, all=TRUE) #if you want to see comparisons of conditions for each gene
#te_pluri = te[grepl(paste(pluri, collapse = "|"), te$gene.x),] #only look at genes from a specified signature
dif_cv = dif_cv[complete.cases(dif_cv), ] #keep only rows without NA

#dif_cv$dMean = (dif_cv$mean.y - dif_cv$mean.x)/dif_cv$mean.x
dif_cv$dMean = (dif_cv$mean.y - dif_cv$mean.x)#/dif_cv$mean.x
dif_cv$dCV = dif_cv$cv.y - dif_cv$cv.x
dif_cv$dFano = dif_cv$fano.y - dif_cv$fano.x


dif_cv = dif_cv[order(dif_cv$gene_sig.x, decreasing = F),]

#plot all genes, colored by gene signature
lmFit_pluri = lm(dif_cv[dif_cv$gene_sig.x=="pluri",15] ~ dif_cv[dif_cv$gene_sig.x=="pluri",14]) #find linear regression of pluripotency genes
lmFit = lm(dif_cv[,15] ~ dif_cv[,14]) #find linear regression of all genes

ggplot(dif_cv, aes(x=dMean, y=dCV, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#D3D3D3", "#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_abline(intercept = coefficients(lmFit)[1], slope = coefficients(lmFit)[2]) +
  geom_abline(intercept = coefficients(lmFit_pluri)[1], slope = coefficients(lmFit_pluri)[2], col = "red") +
  ylim(-1.5, 1.5) + xlim(-1.6, 1.6) +
  theme(legend.position='none')

lmFit
lmFit_pluri

#plot only genes in gene signatures
ggplot(dif_cv[dif_cv$gene_sig.x!=0,], aes(x=dMean, y=dCV, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")
ggplot(dif_cv[dif_cv$gene_sig.x=="pluri",], aes(x=dMean, y=dCV, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")

#separate genes by change in CV and mean (increase or decrease in each) in S3 vs LacZ
cv_gene_clusters = list(mposcpos = subset(dif_cv, dMean > 0 & dCV > 0), mnegcpos = subset(dif_cv, dMean < 0 & dCV > 0), mposcneg = subset(dif_cv, dMean > 0 & dCV < 0), mnegcneg = subset(dif_cv, dMean < 0 & dCV < 0))

# library(xlsx)
# write.xlsx(cv_gene_clusters[[1]], file="/home/data/aqphan/RStudio/DowningLab_Git/cv_gene_clusters.xlsx", sheetName="mposcpos", row.names=FALSE)
# write.xlsx(cv_gene_clusters[[2]], file="/home/data/aqphan/RStudio/DowningLab_Git/cv_gene_clusters.xlsx", sheetName="mnegcpos", append=TRUE, row.names=FALSE)
# write.xlsx(cv_gene_clusters[[3]], file="/home/data/aqphan/RStudio/DowningLab_Git/cv_gene_clusters.xlsx", sheetName="mposcneg", append=TRUE, row.names=FALSE)
# fwrite(cv_gene_clusters[[4]], file="/home/data/aqphan/RStudio/DowningLab_Git/cv_gene_clusters2.csv", row.names = F, col.names = T, showProgress = T)

#plot mean LvS and CV LvS
#blue line is regression, red line is if a 1:1 relationship (slope=1)
ggplot(dif_cv, aes(x=mean.x, y=mean.y)) + geom_point(alpha=0.25) + theme_classic() + geom_segment(aes(x=-5, y=-5, xend=1, yend=1, colour="red")) + stat_smooth(method = "lm")
ggplot(dif_cv, aes(x=cv.x, y=cv.y)) + geom_point(alpha=0.25) + theme_classic() + geom_segment(aes(x=-2.5, y=-2.5, xend=5, yend=5, colour="red")) + stat_smooth(method = "lm")

#plot fano factor LvS
ggplot(dif_cv, aes(x=fano.x, y=fano.y, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#000000", "#D3D3D3", "#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + geom_segment(aes(x=-2, y=-2, xend=.75, yend=.75, colour="#000000"), linetype="dashed") #+ stat_smooth(method = "lm") #+ geom_density_2d(bins=20, adjust=8)

fano_gene_clusters = list(pos = subset(dif_cv, dFano > 0), neg = subset(dif_cv, dFano < 0))

# write.xlsx(fano_gene_clusters[[1]], file="/home/data/aqphan/RStudio/DowningLab_Git/fano_gene_clusters.xlsx", sheetName="pos", row.names=FALSE)
# write.xlsx(fano_gene_clusters[[2]], file="/home/data/aqphan/RStudio/DowningLab_Git/fano_gene_clusters.xlsx", sheetName="neg", append=TRUE, row.names=FALSE)

####################### DEG and GSEA ###################

library(msigdbr)
library(fgsea)
library(ggpubr)

#tutorial: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_05_dge.html#Gene_Set_Enrichment_Analysis_(GSEA)

#identify differentially expressed genes of cell types between LacZ and S3 KD
nccMarkers = FindMarkers(combined, ident.1 = "Shroom3", ident.2 = "LacZ", group.by = 'Condition', subset.ident = 'NCCs')
iMarkers = FindMarkers(combined, ident.1 = "Shroom3", ident.2 = "LacZ", group.by = 'Condition', subset.ident = 'iPSCs')

#Create a gene rank based on the gene expression fold change
markers = iMarkers #nccMarkers
gene_rank <- setNames(markers$avg_log2FC, casefold(rownames(markers), upper = T))

# Download gene sets
msigdbgmt <- msigdbr::msigdbr("Homo sapiens")
msigdbgmt <- as.data.frame(msigdbgmt)

# List available gene sets
unique(msigdbgmt$gs_subcat)

#subset gene sets of interest
sets = c("GO:BP")#, "CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS")
sets = c("CP:WIKIPATHWAYS")

msigdbgmt_subset <- msigdbgmt[msigdbgmt$gs_subcat %in% sets, ]
gmt <- lapply(unique(msigdbgmt_subset$gs_name), function(x) {
  msigdbgmt_subset[msigdbgmt_subset$gs_name == x, "gene_symbol"]
})

names(gmt) <- unique(paste0(msigdbgmt_subset$gs_name, "_", msigdbgmt_subset$gs_exact_source))

# Perform enrichment analysis
fgseaRes <- fgsea(pathways = gmt, stats = gene_rank, minSize = 15)
fgseaRes <- fgseaRes[order(fgseaRes$padj, decreasing = T), ]


plotEnrichment(gmt[["GOBP_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT_GO:0007417"]],
               gene_rank) + labs(title="Regulation of Binding")

#GOBP_BIOLOGICAL_ADHESION_GO:0022610, GOBP_NEURON_DIFFERENTIATION_GO:0030182, GOBP_CELL_PROJECTION_ORGANIZATION_GO:0030030
#GOBP_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT_GO:0007417, GOBP_NEUROGENESIS_GO:0022008


top10_UP <- fgseaRes$pathway[1:2]
dev.off()
plotGseaTable(gmt, gene_rank, fgseaRes, gseaParam = 0.5)
