# Figure 6 Seurat Analysis
# By Andrew Phan

library(dplyr)
library(Seurat)
library(patchwork)
library(future)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(DoMultiBarHeatmap)

########################### QC and selecting cells for further analysis ###########################

### Load in data from GSE241274 using Read10X function from Seurat and store in object named "combined"

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
#filter cells that have >15% mitochondrial counts
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)

########################### Normalizing the data ###########################

#use "LogNormalize" method by default. normalizes the feature expression measurements for each cell by the total expression,
#multiplies this by a scale factor (10,000 by default), and log-transforms the result
#Normalized values are stored in combined[["RNA"]]@data
plan("multiprocess", workers = 4)
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
plot1 + plot2 + theme(legend.position="none") #remove legend if "Viewport has zero dimension(s)"

########################### Scaling the data ###########################

#apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques, like PCA
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1 (This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate)
#The results of this are stored in head(combined[["RNA"]]@scale.data)
all.genes <- rownames(combined)
plan("multiprocess", workers = 1)
combined <- ScaleData(combined, features = all.genes)

#can remove unwanted sources of variation too, see https://satijalab.org/seurat/v3.2/combined3k_tutorial.html


########################### Perform linear dimensional reduction ###########################

#perform PCA on the scaled data. only the previously determined variable features are used as input, default, but can be defined using features argument
#run PCA on unscaled data
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

combined <- RunPCA(combined, features = c("POU5F1", "CD44", "EPCAM", "ALDH1A3"))
combined <- RunPCA(combined, features = c("AGR2", "GPX2", "SLC12A2", "LCN2", "ALDH3A1", "ALDH1A1", "CD44", "EPCAM", "ABCG2", "POU5F1"))


# Examine and visualize PCA results a few different ways
print(combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(combined, dims = 1:2, reduction = "pca")
DimPlot(combined, reduction = "pca")
DimHeatmap(combined, dims = 1, cells = 500, balanced = TRUE) #allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses
#Setting cells to a number plots the 'extreme' cells on both ends of the spectrum, which dramatically speeds plotting for large datasets
DimHeatmap(combined, dims = 1:3, cells = 500, balanced = TRUE) #displays PCs 1 through 15

########################### Determine the 'dimensionality' of the dataset ###########################

#To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a 'metafeature' that combines information across a correlated feature set. The top principal components therefore represent a robust compression of the dataset.

#randomly permute a subset of the data (1% by default) and rerun PCA, constructing a 'null distribution' of feature scores, and repeat this procedure. Identifies 'significant' PCs as those who have a strong enrichment of low p-value features

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
#seems to not be effective for this dataset
combined <- JackStraw(combined, num.replicate = 100)
combined <- ScoreJackStraw(combined, dims = 1:20)

#provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line)
#'Significant' PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line)
#JackStrawPlot(combined, dims = 1:20)

#alternatively, can use elbow plot to find inflection point
ElbowPlot(combined)
#elbow at 6

#NOTES:
#-encourage users to repeat downstream analyses with a different number of PCs. often do not differ dramatically
#-advise users to err on the higher side when choosing this parameter. For example, performing downstream analyses with only 5 PCs does significantly and adversely affect results

########################### Cluster the cells ###########################

#first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity)
combined <- FindNeighbors(combined, dims = 1:2)
#next apply modularity optimization techniques such as the Louvain algorithm, to iteratively group cells together, with the goal of optimizing the standard modularity function
#resolution or granularity is typically good between 0.4-1.2
combined <- FindClusters(combined, resolution = 0.2) #res = 0.2 for 4 clusters. res of 0.62 is original 9 clusters from Suoqin. 0.9 = 13 clusters that look more defined in middle area
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(combined, reduction = "umap")

del = subset(combined, idents=c(3))

combined <- FindNeighbors(combined, dims = 1:2, annoy.metric = "euclidean")
combined <- FindClusters(combined, graph.name = 'RNA_snn', resolution = .1)
combined <- RunUMAP(combined, graph = 'RNA_snn', n.neighbors = 30)
DimPlot(combined, reduction = "umap", group.by = "orig.ident")

FeaturePlot(combined, features = c("CAV1"), split.by = "orig.ident") #
VlnPlot(combined, features = c("CAV1"), split.by = "orig.ident") #



#to save plot
#saveRDS(combined, file = "combined_tutorial.rds")

########################### Finding differentially expressed features (cluster biomarkers) ###########################
#find markers that define clusters via differential expression (identifies pos. and neg. markers of a cluster compared to all other cells by default, ident.1)
#min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups

# find markers for every cluster compared to all remaining cells, report only the positive ones
plan("multiprocess", workers = 4)
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10

x <- DoHeatmap(combined, features = top10$gene) + scale_fill_gradientn(colors = c(rev(brewer.pal(n = 11, name = "RdBu"))))
x <- DoMultiBarHeatmap(combined, features = A549markerGenes, additional.group.by = c("conditions"), disp.max = 1.5, disp.min = 0) + scale_fill_viridis(option = "E")
show(x)

A549markerGenes = c("TOP2A", "CENPF", #cluster 1
                    "ALDH3A1", "GPX2", #cluster 2
                    "PHLDA2", "TFPI", #cluster 3
                    "UQCRQ", "S100A10") #cluster 4

#fwrite(combined.markers, file = "A549_allClusterMarkers_4clust.txt", sep = "\t")

#ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect)
cluster1.markers <- FindMarkers(combined, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#violin plot of probability distributions of marker expressions
VlnPlot(combined, features = c("MS4A1", "CD79A"))

#you can plot raw counts as well
VlnPlot(combined, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

####### EMT ###########
#EMT markers from https://www.nature.com/articles/s41467-020-16066-2
FeaturePlot(combined, features = c("CDH2"))
FeaturePlot(combined, features = c("CDH1"))
#FeaturePlot(combined, features = c("SNAI1")) #not very helpful
FeaturePlot(combined, features = c("PMEPA1"))
FeaturePlot(combined, features = c("TGM2"))
FeaturePlot(combined, features = c("CD44"))
FeaturePlot(combined, features = c("FN1"))
FeaturePlot(combined, features = c("TGFB1"))
#FeaturePlot(combined, features = c("TMSB10"))
FeaturePlot(combined, features = c("TGM2"))
#FeaturePlot(combined, features = c("ANXA2"))
#FeaturePlot(combined, features = c("PHLDA1"))
#FeaturePlot(combined, features = c("AREG"))
FeaturePlot(combined, features = c("B4GALT1"))
#FeaturePlot(combined, features = c("RPL22L1"))
#FeaturePlot(combined, features = c("TUBA1C"))
FeaturePlot(combined, features = c("TPM1"))
#FeaturePlot(combined, features = c("HLA-C"))
FeaturePlot(combined, features = c("LAMC2"))
#FeaturePlot(combined, features = c("IL32"))
#FeaturePlot(combined, features = c("ACTB"))
#FeaturePlot(combined, features = c("GADD45A"))
#FeaturePlot(combined, features = c("MAP1B"))
#FeaturePlot(combined, features = c("TUBA4A"))
#FeaturePlot(combined, features = c("KLF6"))
#FeaturePlot(combined, features = c("LAMB3"))
#FeaturePlot(combined, features = c("CALM2"))
#FeaturePlot(combined, features = c("ITGA2"))
#FeaturePlot(combined, features = c("SH3BGRL3"))
#FeaturePlot(combined, features = c("FHL2"))
FeaturePlot(combined, features = c("CDKN1A"))
#FeaturePlot(combined, features = c("ITGAV"))
FeaturePlot(combined, features = c("TGFBI"))
#FeaturePlot(combined, features = c("ARF4"))
#FeaturePlot(combined, features = c("SPHK1"))
FeaturePlot(combined, features = c("MT2A"))
#FeaturePlot(combined, features = c("IER3"))
#FeaturePlot(combined, features = c("S100A13"))
#FeaturePlot(combined, features = c("B2M"))
#FeaturePlot(combined, features = c("LGALS1"))
#FeaturePlot(combined, features = c("CD63"))
#FeaturePlot(combined, features = c("SDCBP"))
#FeaturePlot(combined, features = c("CD59"))
#FeaturePlot(combined, features = c("MYL12A"))
#FeaturePlot(combined, features = c("S100A11"))
FeaturePlot(combined, features = c("HLA-A"))
#FeaturePlot(combined, features = c("OCIAD2"))
#FeaturePlot(combined, features = c("HSPA5"))
#FeaturePlot(combined, features = c("HSP90B1"))
#FeaturePlot(combined, features = c("PDIA3"))
FeaturePlot(combined, features = c("PLEK2"))
#FeaturePlot(combined, features = c("PLAUR"))
#FeaturePlot(combined, features = c("PLAU"))
#FeaturePlot(combined, features = c("SMOX"))
FeaturePlot(combined, features = c("DHRS7"))
#FeaturePlot(combined, features = c("NT5E"))
#FeaturePlot(combined, features = c("TIMP1"))
#FeaturePlot(combined, features = c("ANGPTL4"))
FeaturePlot(combined, features = c("HMGA2"))
#FeaturePlot(combined, features = c("STEAP1"))
FeaturePlot(combined, features = c("IL11"))
FeaturePlot(combined, features = c("SMTN"))
FeaturePlot(combined, features = c("P4HA2"))
FeaturePlot(combined, features = c("SERPINE2"))
#FeaturePlot(combined, features = c("TP53I3"))
FeaturePlot(combined, features = c("PTHLH"))
#FeaturePlot(combined, features = c("DLC1"))
#FeaturePlot(combined, features = c("TUBA1A"))
#FeaturePlot(combined, features = c("UBC"))
#FeaturePlot(combined, features = c("GAL"))
#FeaturePlot(combined, features = c("FABPS"))
#FeaturePlot(combined, features = c("HN1"))
#FeaturePlot(combined, features = c("TFPI2"))
#FeaturePlot(combined, features = c("DUSP4"))
#FeaturePlot(combined, features = c("ASNS"))
#FeaturePlot(combined, features = c("GARS"))
#FeaturePlot(combined, features = c("CTSL"))

#https://ar.iiarjournals.org/content/38/7/3797
FeaturePlot(combined, features = c("SNAI2")) #
FeaturePlot(combined, features = c("MMP2")) #



combined = AddModuleScore(combined, features = list(c("CDH2", "CDH1", "PMEPA1", "TGM2", "CD44", "FN1", "TGFB1", "TGM2", "B4GALT1", "TPM1", "LAMC2", "CDKN1A", "TGFBI", "MT2A", "HLA-A", "PLEK2", "DHRS7", "HMGA2", "IL11", "SMTN", "P4HA2", "SERPINE2", "PTHLH")), name = "Up.EMT")
FeaturePlot(combined, features = c("Up.EMT1"), cols = heat.colors(100))


#EMT Downregulated
FeaturePlot(combined, features = c("ELF3"))
FeaturePlot(combined, features = c("KRT8"))
FeaturePlot(combined, features = c("KRT18"))
FeaturePlot(combined, features = c("AGR2"))
FeaturePlot(combined, features = c("ATP1B1"))
FeaturePlot(combined, features = c("EPHX1"))
FeaturePlot(combined, features = c("KRT19"))
FeaturePlot(combined, features = c("PLA2G16"))
#FeaturePlot(combined, features = c("HCFC1R1")) #not very helpful
FeaturePlot(combined, features = c("DHCR24"))
FeaturePlot(combined, features = c("CDKN3"))
FeaturePlot(combined, features = c("GINS2"))
#FeaturePlot(combined, features = c("ACAT2")) #not very helpful
FeaturePlot(combined, features = c("HIST1H4C"))

combined = AddModuleScore(combined, features = list(c("ELF3", "KRT8", "KRT18", "AGR2", "ATP1B1", "EPHX1", "KRT19", "PLA2G16", "DHCR24", "CDKN3", "GINS2", "HIST1H4C")), name = "Down.EMT")
FeaturePlot(combined, features = c("Down.EMT1"), cols = heat.colors(100))

#Consistent EMT list
combined = AddModuleScore(combined, features = list(c("BMP4", "SOX9", "SNAI2", "CDH2", "CDH11")), name = "EMT")
FeaturePlot(combined, features = c("EMT1"), split.by = "orig.ident")

#########

#### Cellular Aldehyde Metabolic Process #######
#in top 100 of 4 clusters
FeaturePlot(combined, features = c("ALDH3A1"))
FeaturePlot(combined, features = c("TKT"))
FeaturePlot(combined, features = c("ALDH3B1"))
FeaturePlot(combined, features = c("IDH1"))
FeaturePlot(combined, features = c("ALDH1A1"))

FeaturePlot(combined, features = c("NOTCH3")) #https://cancerres.aacrjournals.org/content/70/23/9937.short


combined = AddModuleScore(combined, features = list(c("ALDH3A1", "TKT", "ALDH3B1", "IDH1", "ALDH1A1")), name = "ALDH.up")
FeaturePlot(combined, features = c("ALDH.up1"), cols = heat.colors(100))

VlnPlot(combined, features = c("ALDH.up1"), split.by = "orig.ident", group.by = "orig.ident", split.plot = T, pt.size = 0)#, idents = 1)

VlnPlot(combined, features = c("CSRP1"), split.by = "orig.ident",  split.plot = T, pt.size = 0)#, idents = 1)


#######

########### Pluripotency #############
FeaturePlot(combined, features = c("TET1")) #Higher in CSRP1
FeaturePlot(combined, features = c("CDH1"))#, split.by = "orig.ident") #Higher in CSRP1
FeaturePlot(combined, features = c("ID1")) #Higher in CSRP1

FeaturePlot(combined, features = c("LIN28B")) #Higher in CSRP1
FeaturePlot(combined, features = c("TEAD4")) #Higher in CSRP1
FeaturePlot(combined, features = c("ITGB5")) #Higher in CSRP1 - Invasion, EMT, cancer stem cell biology. enhances invasiveness and proliferation in breast cancer https://www.nature.com/articles/onc2012320
FeaturePlot(combined, features = c("DNMT3A")) #Some expression
FeaturePlot(combined, features = c("SALL4")) #Some expression

FeaturePlot(combined, features = c("CTNNB1")) #Some expression
FeaturePlot(combined, features = c("CD44")) #Higher in CSRP1
#FeaturePlot(combined, features = c("CD24")) #Higher in CSRP1

FeaturePlot(combined, features = c("ALCAM")) #CD166 - Higher in LacZ

#################

###### CSC Markers ########
FeaturePlot(combined, features = c("AGR2")) #
FeaturePlot(combined, features = c("GPX2"), split.by = "orig.ident") #
FeaturePlot(combined, features = c("SLC12A2")) #
FeaturePlot(combined, features = c("LCN2")) #
FeaturePlot(combined, features = c("ALDH3A1")) #
FeaturePlot(combined, features = c("ALDH1A1")) #More ubiquitous

FeaturePlot(combined, features = c("ABCG2")) #
#FeaturePlot(combined, features = c("ABCC1")) #MRP1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5897668/
FeaturePlot(combined, features = c("THY1")) #CD90
FeaturePlot(combined, features = c("CD44")) #
#FeaturePlot(combined, features = c("PROM1")) #CD133

FeaturePlot(combined, features = c("NOTCH1")) #NOTCH1-3 in middle. 3 is specific to the bottom half
FeaturePlot(combined, features = c("DLL4")) #JAG1 (and to some extent JAG2) and DLL3 in middle

#https://www.mdpi.com/2072-6694/13/12/2996/pdf
FeaturePlot(combined, features = c("EPCAM"))


FeaturePlot(combined, features = c("SPP1")) #OPN
FeaturePlot(combined, features = c("ABCG2"))
FeaturePlot(combined, features = c("CCNB1"))
FeaturePlot(combined, features = c("TM4SF4"))
FeaturePlot(combined, features = c("TCF4"))

#POTENTIALLY ADD TO FIGURE OR SUPPLEMENTARY
#https://www.sciencedirect.com/science/article/pii/S0092867418303581?via%3Dihub
FeaturePlot(combined, features = c("TP53")) #correlated with higher stemness and smoking in LUAD
FeaturePlot(combined, features = c("FOXM1")) #correlated with higher LUAD stemness; higher in CSRP1
FeaturePlot(combined, features = c("CCNB1")) #correlated with higher LUAD stemness; slightly higher in CSRP1
FeaturePlot(combined, features = c("ANXA1")) #correlated with lower LUAD stemness; lower in CSRP1

#combined = AddModuleScore(combined, features = list(c("AGR2", "GPX2", "SLC12A2", "LCN2", "ALDH3A1", "ALDH1A1")), name = "CSC") #DID NOT USE
FeaturePlot(combined, features = c("CSC1"), cols = heat.colors(100))

FeaturePlot(combined, features = c("SOX2")) #




########

####### Metastasis ########
FeaturePlot(combined, features = c("LCN2")) #
FeaturePlot(combined, features = c("S100A4")) #
FeaturePlot(combined, features = c("ANG")) #
FeaturePlot(combined, features = c("AKR1C2")) #
FeaturePlot(combined, features = c("CEACAM6")) #
FeaturePlot(combined, features = c("TFF1")) #
FeaturePlot(combined, features = c("ANXA13")) #

###########


####### Fibrinogen ########
FeaturePlot(combined, features = c("FGL1")) #
FeaturePlot(combined, features = c("FGB")) #
FeaturePlot(combined, features = c("FGG")) #

#########


########### Epithelial Cell #########
FeaturePlot(combined, features = c("EPCAM", "KRT8", "KRT18")) #KRT8/18 https://www.cell.com/cell-reports/pdfExtended/S2211-1247(14)00705-0

###################

########## Trophoblast########
FeaturePlot(combined, features = c("TOP2A", "UBE2C")) #tumorigenesis/drug resistance. link to iPSCs?

###################

####### Notch ##########
FeaturePlot(combined, features = c("JAG1"))
# FeaturePlot(combined, features = c("NODAL"))
# FeaturePlot(combined, features = c("LEFTY1"))
# FeaturePlot(combined, features = c("LEFTY2"))
# FeaturePlot(combined, features = c("ROCK1"))
# FeaturePlot(combined, features = c("ROCK2"))
FeaturePlot(combined, features = c("SHROOM3"))

############

####### Planar Cell Polarity Pathway ##########

#Non-Canonical WNT/PCP Involved in NCC Migration/EMT/Differentiation https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5811199/
FeaturePlot(combined, features = c("FZD7"), split.by = "orig.ident") #receptor for WNT11 ligand. Required for CIL migration
#FeaturePlot(combined, features = c("WNT11"))
FeaturePlot(combined, features = c("DVL1")) #also required for CIL migration
FeaturePlot(combined, features = c("DVL2"), split.by = "orig.ident") #interacts with CSRP1
FeaturePlot(combined, features = c("PRICKLE1"))

#PCP Pathway
FeaturePlot(combined, features = c("ANKRD6"))#, split.by = "orig.ident") #diversin gene, interacts with CSRP1 in zebrafish
#ANKRD10
#FeaturePlot(combined, features = c("FZD1"))
FeaturePlot(combined, features = c("FZD3"), split.by = "orig.ident")
FeaturePlot(combined, features = c("FZD6"), split.by = "orig.ident")

FeaturePlot(combined, features = c("DVL1"), split.by = "orig.ident")
FeaturePlot(combined, features = c("DVL2"), split.by = "orig.ident")
FeaturePlot(combined, features = c("DVL3"), split.by = "orig.ident")

FeaturePlot(combined, features = c("CELSR1"), split.by = "orig.ident")
#FeaturePlot(combined, features = c("CELSR2"), split.by = "orig.ident")
#FeaturePlot(combined, features = c("CELSR3"), split.by = "orig.ident")

FeaturePlot(combined, features = c("PRICKLE1"), split.by = "orig.ident")
FeaturePlot(combined, features = c("PRICKLE2"), split.by = "orig.ident") #higher in LacZ
FeaturePlot(combined, features = c("PRICKLE3"), split.by = "orig.ident")

combined = AddModuleScore(combined, features = list(c("DVL1", "DVL2", "DVL3", "PRICKLE1", "PRICKLE2", "PRICKLE3", "FZD3", "FZD6", "CELSR1", "DAAM1")), name = "PCP.PATHWAY_A549")
FeaturePlot(combined, features = c("PCP.PATHWAY_A5491"))
VlnPlot(combined, features = c("PCP.PATHWAY_A5491"), split.by = "orig.ident",group.by = "orig.ident", split.plot = T, pt.size = 0)

combined = AddModuleScore(combined, features = list(c("WNT5A", "VANGL2", "FZD3", "DVL2", "CELSR1", "PRICKLE1", "FZD7", "WNT11", "DVL1")), name = "PCP.PATHWAY")
VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "orig.ident",group.by = "orig.ident", split.plot = T, pt.size = 0)


########## MET/EMT ##########
#https://www.nature.com/articles/s41598-019-43021-z

#EMT
FeaturePlot(combined, features = c("ZEB1"))
FeaturePlot(combined, features = c("ZEB2"))
FeaturePlot(combined, features = c("TWIST1"))
FeaturePlot(combined, features = c("TWIST2"))
FeaturePlot(combined, features = c("SNAI1"))
FeaturePlot(combined, features = c("SNAI2"))
FeaturePlot(combined, features = c("VIM"))

VlnPlot(combined, features = c("ZEB1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("ZEB2"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("TWIST1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("TWIST2"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("SNAI1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("SNAI2"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("VIM"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("COL6A1"), split.by = "orig.ident", pt.size = 0)#, idents = 1) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7046610/#MOESM1

#VIM significantly correlated in TCGA LUAD cohort according to EMTome.org. Spearman Correlation >0.65 (PTRF not in data)
VlnPlot(combined, features = c("EMP3"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("AXL"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("GLIPR1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("GNAI2"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
#see EMTome paper: https://www.nature.com/articles/s41416-020-01178-9

combined = AddModuleScore(combined, features = list(c("VIM", "EMP3", "AXL", "GLIPR1", "GNAI2")), name = "EMT.LUAD")
VlnPlot(combined, features = c("EMT.LUAD1"), split.by = "orig.ident", group.by = "orig.ident", pt.size = 0)#, idents = 1)


#MET
FeaturePlot(combined, features = c("GATA3"))
FeaturePlot(combined, features = c("HNF1A"))
FeaturePlot(combined, features = c("TP63"))
FeaturePlot(combined, features = c("FOXA3"))
FeaturePlot(combined, features = c("YBX2"))
FeaturePlot(combined, features = c("KLF4"))
FeaturePlot(combined, features = c("GRHL3"))
FeaturePlot(combined, features = c("EHF"))
FeaturePlot(combined, features = c("ANKRD22"))
FeaturePlot(combined, features = c("FOXA1"))
FeaturePlot(combined, features = c("GRHL1"))
FeaturePlot(combined, features = c("ZNF165"))
#FeaturePlot(combined, features = c("GRHL2"))
FeaturePlot(combined, features = c("ELF3"))
FeaturePlot(combined, features = c("IRF6"))
FeaturePlot(combined, features = c("OVOL2"))
FeaturePlot(combined, features = c("ID1"))
FeaturePlot(combined, features = c("CDH1"))

VlnPlot(combined, features = c("GATA3"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("HNF1A"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("TP63"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("FOXA3"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("YBX2"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("KLF4"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("GRHL3"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("EHF"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("ANKRD22"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("FOXA1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("GRHL1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("ZNF165"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("GRHL2"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("ELF3"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("IRF6"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("OVOL2"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("ID1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)

VlnPlot(combined, features = c("CDH1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
#CDH1 significant correlation (Spearman correlation coefficient > 0.5) in LUAD with RNA expression with FDR(BH) < 0.05 (ANKRD56 not expressed)
VlnPlot(combined, features = c("TMEM30B"), split.by = "orig.ident", pt.size = 0)#, idents = 1)  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7046610/#MOESM1
VlnPlot(combined, features = c("CGN"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("PRSS8"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("CDS1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("ESRP1"), split.by = "orig.ident", pt.size = 0)#, idents = 1)
VlnPlot(combined, features = c("ESRP2"), split.by = "orig.ident", pt.size = 0)#, idents = 1)

combined = AddModuleScore(combined, features = list(c("CDH1", "TMEM30B", "CGN", "CDS1", "ESRP1", "ESRP2")), name = "MET.try")
VlnPlot(combined, features = c("MET.try1"), split.by = "orig.ident", group.by = "orig.ident", split.plot = T)#, idents = 1)

#JUN signaling tied to CSRP1 expression in zebrafish: https://www.pnas.org/content/104/27/11274
FeaturePlot(combined, features = c("DVL1"), split.by = "orig.ident")
VlnPlot(combined, features = c("DVL3"), split.by = "orig.ident", group.by = "orig.ident", split.plot = T, pt.size = 0)

#############

#generates an expression heatmap for given cells and features. ex: top 10 markers
top10 <- combined2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(combined2, features = top10$gene) + NoLegend()

########################### Assigning cell type identity to clusters ###########################

new.cluster.ids <- c("Class I MHC Expressing", "Neural Crest Cells", "Highly Proliferating Trophoblasts", "Fibroblasts", "Epithelial", "iPSCs", "Trophoblasts", "Neuroepithelial", "Smooth Muscle Cells", "Neural")
names(new.cluster.ids) <- levels(combined2)
combined3 <- RenameIdents(object = combined2, `0` = "CIMHC+", `1` = "NCCs", `2` = "Prolif.Troph.", `3` = "Fibr.", `4` = "Epi.", `5` = "iPSCs",  `6` = "Troph.", `7` = "Neur.Epi.", `8` = "SMCs", `9` = "Neur.")
#combined3 <- RenameIdents(object = combined2, `0` = "zero", `1` = "one", `2` = "two", `3` = "Class I MHC+", `4` = "four", `5` = "five",  `6` = "six", `7` = "seven", `8` = "eight", `9` = "nine")
#combined2 <- RenameIdents(object = combined2, '0' = "Class I MHC+", '1' = "Neural Crest", '2' = "Prolif. Trophoblasts", '3' = "Fibroblasts", '4' = "Epithelial", '5' = "iPSCs",  '6' = "Trophoblasts", '7' = "Neuroepithelial", '8' = "Smooth Muscle", '9' = "Neural")
#View(as.data.frame(Idents(combined2)))

DimPlot(combined3, reduction = "umap", label = F, pt.size = 0.5) #+ NoLegend()

#to save
#saveRDS(combined, file = "combined3k_final.rds")

metadata <- combined@meta.data

metadata %>%
  ggplot(aes(color=seurat_clusters, x=nCount_RNA, fill= seurat_clusters)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

############# Figures ################

#UMAP of clusters
cols_A549 = viridis(n = 4, option = "E")
cols_A549 = c(cols_A549[1], cols_A549[4], cols_A549[3], cols_A549[2]) #rearrange colors to make CSC yellow
DimPlot(combined, reduction = "umap", cols = cols_A549)

#separate conditions
lacZ_ident = unique(combined$orig.ident)[1]
CSRP1_ident = unique(combined$orig.ident)[2]
lacZseurat = subset(x = combined, subset = orig.ident == lacZ_ident)
CSRP1seurat = subset(x = combined, subset = orig.ident == CSRP1_ident)

DimPlot(lacZseurat, reduction = "umap", cols = cols_A549)
DimPlot(CSRP1seurat, reduction = "umap", cols = cols_A549)
#export as 6x6in

#Heatmap of marker genes per cluster
x <- DoHeatmap(combined, features = top10$gene) + scale_fill_gradientn(colors = c(rev(brewer.pal(n = 11, name = "RdBu"))))
x <- DoMultiBarHeatmap(combined, features = top10$gene, disp.min = 0, disp.max = 2) + scale_fill_viridis(option = "E") #adjust contrast with disp.min/max
show(x)

#violin plots of marker genes (orig.ident = CSRP1 vs. LacZ KD conditions)
VlnPlot(combined, features = c("ALDH1A3"), split.by = "orig.ident", group.by = "orig.ident")#, idents = 1)
VlnPlot(combined, features = c("POU5F1"), split.by = "orig.ident", group.by = "orig.ident")#, idents = 1) #similar if only group 1, but much fewer cells in lacZ
VlnPlot(combined, features = c("CD44"), split.by = "orig.ident", group.by = "orig.ident")#, idents = 1)
VlnPlot(combined, features = c("EPCAM"), split.by = "orig.ident", group.by = "orig.ident")#, idents = 1)
#all show higher expression in CSRP1 KD
VlnPlot(combined, features = c("ALDH1A3", "POU5F1", "CD44", "EPCAM", "ABCG2", "NOTCH3"), split.plot = T, split.by = "orig.ident", group.by = "orig.ident", idents = 2)#, idents = 1)

#Cancer Stem Cell Signature
cscMarkerGenes = c("ALDH1A3", "POU5F1", "CD44", "EPCAM", "ABCG2", "NOTCH3")
combined = AddModuleScore(combined, features = list(cscMarkerGenes), name = "CSC")
FeaturePlot(combined, features = c("CSC1"))

VlnPlot(combined, features = c("CSC_TRY1"), split.by = "orig.ident", group.by = "orig.ident", pt.size = 0, split.plot = T) + geom_hline(yintercept = c(-0.1441363, -0.0986022), col = c("blue", "red"))
test1 = VlnPlot(combined, features = c("CSC1"), split.by = "orig.ident", group.by = "orig.ident", pt.size = 0, split.plot = T)
test1 = test1$data
lacZ_ident = unique(combined$orig.ident)[1]
CSRP1_ident = unique(combined$orig.ident)[2]
wilcox.test(x = test1[test1$split== lacZ_ident,1], y = test1[test1$split==CSRP1_ident,1]) #p-value <2.2e-16
median(test1[test1$split==lacZ_ident,1]) #median LacZ 0 = -0.1441363
median(test1[test1$split==CSRP1_ident,1]) #median CSRP1 0 = -0.0986022

VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "orig.ident") #, idents = 1) #can look at cells from cluster 1
#cluster 0 has much less of an increase in CSC sig than cluster 1 for LvC

##################
VlnPlot(combined, features = c("CELSR1"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("VANGL1"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("VANGL2"), split.by = "orig.ident", split.plot = T)

VlnPlot(combined, features = c("PRICKLE1"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("PRICKLE2"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("PRICKLE3"), split.by = "orig.ident", split.plot = T)

VlnPlot(combined, features = c("DVL3"), split.by = "orig.ident", split.plot = T)

VlnPlot(combined, features = c("ANKRD6"), split.by = "orig.ident", split.plot = T)

VlnPlot(combined, features = c("WNT5A"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("WNT5B"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("WNT11"), split.by = "orig.ident", split.plot = T)

VlnPlot(combined, features = c("FZD1"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("FZD2"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("FZD3"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("FZD4"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("FZD5"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("FZD6"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("FZD7"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("FZD8"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("FZD9"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("FZD10"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("MCAM"), split.by = "orig.ident", split.plot = T)

VlnPlot(combined, features = c("ROR1"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("ROR2"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("PTK7"), split.by = "orig.ident", split.plot = T)

VlnPlot(combined, features = c("MAGI3"), split.by = "orig.ident", split.plot = T)

#https://www.spandidos-publications.com/or/14/6/1583
VlnPlot(combined, features = c("DAAM1"), split.by = "orig.ident", split.plot = T, pt.size = 0, group.by = "orig.ident") #higher in LacZ. JNK signaling through MAP3K and MAP2K4/7
VlnPlot(combined, features = c("DAAM2"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("MAP2K4"), split.by = "orig.ident", group.by = "orig.ident", split.plot = T, pt.size = 0)
VlnPlot(combined, features = c("MAP2K7"), split.by = "orig.ident", group.by = "orig.ident", split.plot = T, pt.size = 0)


VlnPlot(combined, features = c("ANKRD6"), split.by = "orig.ident", group.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("NKD1"), split.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("NKD2"), split.by = "orig.ident", split.plot = T)

VlnPlot(combined, features = c("WNT5A"), split.by = "orig.ident", group.by = "orig.ident", split.plot = T)
VlnPlot(combined, features = c("WNT5B"), split.by = "orig.ident", group.by = "orig.ident", split.plot = T)

################

#RidgePlot(combined, features = c(cscMarkerGenes), group.by = "orig.ident", idents = 1) #ridgeplots do not look as nice as violin plots

#saveRDS(object = combined, file = "A549_combined_reclustered.rds")

DotPlot(combined, features = c("ALDH1A3", "POU5F1", "CD44", "EPCAM", "ABCG2", "NOTCH3"), split.by = "orig.ident", dot.scale = 20, idents = c(0,1), scale.max = 30, col.min = 0)
DotPlot(lacZseurat, features = c("ALDH1A3", "POU5F1", "CD44", "EPCAM", "ABCG2", "NOTCH3"), dot.scale = 20, idents = c(1,2), scale.max = 30, col.min = 0)
DotPlot(CSRP1seurat, features = c("ALDH1A3", "POU5F1", "CD44", "EPCAM", "ABCG2", "NOTCH3"), dot.scale = 20, idents = c(1,2), scale.max = 30 ,col.min = 0) #c("lightgrey", "#cb181d")

DotPlot(CSRP1seurat, features = c("ALDH1A3", "POU5F1", "CD44", "EPCAM", "ABCG2", "NOTCH3"), dot.scale = 20, idents = c(1,2), scale.max = 30) +
  scale_color_viridis(option = "E")
  #scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")

#CSRP1 Expression
VlnPlot(combined, features = c("CSRP1"), split.by = "orig.ident", split.plot = T, pt.size = 0, idents = c(0,1)) + geom_hline(yintercept = c(0.7894767, 2.616364e-06, 0.6463445, 2.169778e-06), color = c("blue", "blue", "red", "red"))
test1 = VlnPlot(combined, features = c("CSRP1"), split.by = "orig.ident", idents = c(0,1), split.plot = T)
test1 = test1$data
lacZ_ident = unique(combined$orig.ident)[1]
CSRP1_ident = unique(combined$orig.ident)[2]
wilcox.test(x = test1[test1$split== lacZ_ident & test1$ident==0,1], y = test1[test1$split==CSRP1_ident & test1$ident==0,1]) #p-value <2.2e-16
wilcox.test(x = test1[test1$split==lacZ_ident & test1$ident==1,1], y = test1[test1$split==CSRP1_ident & test1$ident==1,1]) #p-value <2.2e-16
median(test1[test1$split==lacZ_ident & test1$ident==0,1]) #median LacZ 0 = 0.7894767
median(test1[test1$split==CSRP1_ident & test1$ident==0,1]) #median CSRP1 0 = 2.616364e-06
median(test1[test1$split==lacZ_ident & test1$ident==1,1]) #median LacZ 1 = 0.6463445
median(test1[test1$split==CSRP1_ident & test1$ident==1,1]) #median CSRP1 1 = 2.169778e-06

#PCP Pathway Expression
VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "orig.ident", split.plot = T, pt.size = 0, idents = c(0,1))# + geom_hline(yintercept = c(0.5521554, 0.3754184, 0.6533828, 0.4168875), color = c("blue", "blue", "red", "red"))
test1 = VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "orig.ident", idents = c(0,1), split.plot = T)
test1 = test1$data
lacZ_ident = unique(combined$orig.ident)[1]
CSRP1_ident = unique(combined$orig.ident)[2]
wilcox.test(x = test1[test1$split== lacZ_ident & test1$ident==0,1], y = test1[test1$split==CSRP1_ident & test1$ident==0,1]) #p-value <2.2e-16
wilcox.test(x = test1[test1$split==lacZ_ident & test1$ident==1,1], y = test1[test1$split==CSRP1_ident & test1$ident==1,1]) #p-value 2.69e-05
median(test1[test1$split==lacZ_ident & test1$ident==0,1]) #median LacZ 0 = 0.5521554
median(test1[test1$split==CSRP1_ident & test1$ident==0,1]) #median CSRP1 0 = 0.3754184
median(test1[test1$split==lacZ_ident & test1$ident==1,1]) #median LacZ 1 = 0.6533828
median(test1[test1$split==CSRP1_ident & test1$ident==1,1]) #median CSRP1 1 = 0.4168875

#test PCP LvC
#separate clusters 0 and 1
test1 = VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "orig.ident", split.plot = T, idents = c(0,1))
test1 = test1$data
wilcox.test(x = test1[test1$split==lacZ_ident,1], y = test1[test1$split==CSRP1_ident,1]) #p-value < 2.2e-16
median(test1[test1$split==lacZ_ident,1]) #median LacZ 1 = -0.01344884
median(test1[test1$split==CSRP1_ident,1]) #median CSRP1 1 = -0.004000061
VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "orig.ident",group.by = "orig.ident", split.plot = T, pt.size = 0) + geom_hline(yintercept = c(-0.01344884,-0.004000061), color = c("blue", "red"))
lacZ_ident = unique(combined$orig.ident)[1]
CSRP1_ident = unique(combined$orig.ident)[2]
median(test1[test1$split==lacZ_ident & test1$ident==0,1]) #median LacZ 0 = 0.5521554
median(test1[test1$split==CSRP1_ident & test1$ident==0,1]) #median CSRP1 0 = 0.3754184
median(test1[test1$split==lacZ_ident & test1$ident==1,1]) #median LacZ 1 = 0.6533828
median(test1[test1$split==CSRP1_ident & test1$ident==1,1]) #median CSRP1 1 = 0.4168875

#separate by condition
VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "orig.ident",group.by = "orig.ident", split.plot = T, pt.size = 0)
test1 = VlnPlot(combined, features = c("PCP.PATHWAY1"), split.by = "orig.ident",group.by = "orig.ident", split.plot = T)
test1 = test1$data
lacZ_ident = unique(combined$orig.ident)[1]
CSRP1_ident = unique(combined$orig.ident)[2]
median(test1[test1$split==lacZ_ident,1]) #median LacZ 0 = 0.5521554
median(test1[test1$split==CSRP1_ident,1]) #median LacZ 0 = 0.5521554

#EMT
# VlnPlot(combined, features = c("EMT1"), split.by = "orig.ident", group.by = "orig.ident", pt.size = 0)# + geom_hline(yintercept = c(0.5521554, 0.3754184, 0.6533828, 0.4168875), color = c("blue", "blue", "red", "red"))
# VlnPlot(combined, features = c("EMT1"), split.by = "orig.ident", split.plot = T, pt.size = 0)# + geom_hline(yintercept = c(0.5521554, 0.3754184, 0.6533828, 0.4168875), color = c("blue", "blue", "red", "red"))
# test1 = VlnPlot(combined, features = c("EMT1"), split.by = "orig.ident", idents = c(0,1), split.plot = T)
# test1 = test1$data
# lacZ_ident = unique(combined$orig.ident)[1]
# CSRP1_ident = unique(combined$orig.ident)[2]
# wilcox.test(x = test1[test1$split== lacZ_ident & test1$ident==0,1], y = test1[test1$split==CSRP1_ident & test1$ident==0,1]) #p-value <2.2e-16
# wilcox.test(x = test1[test1$split==lacZ_ident & test1$ident==1,1], y = test1[test1$split==CSRP1_ident & test1$ident==1,1]) #p-value <2.2e-16
# median(test1[test1$split==lacZ_ident & test1$ident==0,1]) #median LacZ 0 = 0.5521554
# median(test1[test1$split==CSRP1_ident & test1$ident==0,1]) #median CSRP1 0 = 0.3754184
# median(test1[test1$split==lacZ_ident & test1$ident==1,1]) #median LacZ 1 = 0.6533828
# median(test1[test1$split==CSRP1_ident & test1$ident==1,1]) #median CSRP1 1 = 0.4168875

#test EMT LvC
test1 = VlnPlot(combined, features = c("EMT.try1"), split.by = "orig.ident",group.by = "orig.ident", split.plot = T)
test1 = test1$data
lacZ_ident = unique(combined$orig.ident)[1]
CSRP1_ident = unique(combined$orig.ident)[2]
wilcox.test(x = test1[test1$split==lacZ_ident,1], y = test1[test1$split==CSRP1_ident,1]) #p-value < 2.2e-16
median(test1[test1$split==lacZ_ident,1]) #median LacZ 1 = 0.2672915
median(test1[test1$split==CSRP1_ident,1]) #median CSRP1 1 = 0.1614985
VlnPlot(combined, features = c("EMT.try1"), split.by = "orig.ident",group.by = "orig.ident", split.plot = T, pt.size = 0) + geom_hline(yintercept = c(0.2672915,0.1614985), color = c("blue", "red"))

#test MET LvC
test1 = VlnPlot(combined, features = c("MET.try1"), split.by = "orig.ident",group.by = "orig.ident", split.plot = T)
test1 = test1$data
lacZ_ident = unique(combined$orig.ident)[1]
CSRP1_ident = unique(combined$orig.ident)[2]
wilcox.test(x = test1[test1$split==lacZ_ident,1], y = test1[test1$split==CSRP1_ident,1]) #p-value < 2.2e-16
mean(test1[test1$split==lacZ_ident,1]) #median LacZ 1 = -0.004446068
mean(test1[test1$split==CSRP1_ident,1]) #median CSRP1 1 = 0.0002877844
VlnPlot(combined, features = c("MET.try1"), split.by = "orig.ident",group.by = "orig.ident", split.plot = T, pt.size = 0) + geom_hline(yintercept = c(-0.004446068,0.0002877844), color = c("blue", "red"))
#NOTE: MEAN FOR MET, NOT MEDIAN


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
sigs = list(CSC = c("ALDH1A3", "POU5F1", "CD44", "EPCAM", "ABCG2", "NOTCH3"),
            MET = c("CDH1", "TMEM30B", "CGN", "CDS1", "ESRP1", "ESRP2"),
            EMT = c("VIM", "EMP3", "AXL", "GLIPR1", "GNAI2"))

#extract the normalized expression data
expr_all = as.data.frame(combined@assays$RNA@data)
expr_all$gene_sig = 0 #set default gene sig to 0

for(i in 1:length(sigs[])){ #go through each gene signature and identify genes in the signature
  expr_all$gene_sig[match(sigs[[i]], rownames(expr_all))] = names(sigs)[i]
}

#extract cell type metadata and filter for only cells from cell type
ctype = as.data.frame(combined@active.ident)
ctype$id = rownames(ctype)
celltype_ids = list(zero=ctype[ctype[,1]==0,], one=ctype[ctype[,1]==1,], two=ctype[ctype[,1]==2,], three=ctype[ctype[,1]==3,])
#CIMHC=ctype[ctype[,1]=="CIMHC+",]

#choose cell type to subset data by (optional)
#zero one two three
cellsub = celltype_ids[["one"]] #choose cell type to subset by
expr = expr_all[,match(cellsub$id, colnames(expr_all))] #subset the cells for CV and mean calculations
expr$gene_sig = expr_all$gene_sig

#expr = expr_all #for all cell types

#calculate CV and mean of normalized expression for each gene in the LacZ condition
lacZ_norm = grepl("_LacZ",colnames(expr))
lacZ_norm = expr[,lacZ_norm]
lacZ_cvmean = calcCV(lacZ_norm)
lacZ_cvmean$gene = rownames(expr)
lacZ_cvmean$cond = "LacZ"
lacZ_cvmean$gene_sig = expr$gene_sig

#calculate CV and mean of normalized expression for each gene in the SHROOM3 condition
CSRP1_norm = grepl("_CSRP1",colnames(expr))
CSRP1_norm = expr[,CSRP1_norm]
CSRP1_cvmean = calcCV(CSRP1_norm)
CSRP1_cvmean$gene = rownames(expr)
CSRP1_cvmean$cond = "CSRP1"
CSRP1_cvmean$gene_sig = expr$gene_sig

#merge the two conditions to one dataframe
comb_cvmean = rbind(lacZ_cvmean, CSRP1_cvmean)

#create plots of CV vs mean
cvL = ggplot(lacZ_cvmean, aes(x=mean, y=cv)) + geom_point() + theme_classic()
cvC = ggplot(CSRP1_cvmean, aes(x=mean, y=cv)) + geom_point() + theme_classic()
cv_comb = ggplot(comb_cvmean[sample(nrow(comb_cvmean)),], aes(x=mean, y=cv, col=cond)) + geom_point(alpha=0.75) + theme_classic() #shuffles rows to randomize order of points
cv_comb

#look at CV vs mean for only subsets of genes
sub = "EMT"
sub = "MET"
sub = "CSC"

cv_sub = grepl(paste(sub, collapse = "|"), comb_cvmean$gene)
cv_sub = comb_cvmean[cv_sub,]
cv_pl = ggplot(cv_sub[sample(nrow(cv_sub)),], aes(x=mean, y=cv, col=cond)) + geom_point(alpha=0.75) + theme_classic()
cv_pl

#plot the difference of mean and CV of the same genes
dif_cv <- merge(lacZ_cvmean, CSRP1_cvmean, by=0, all=TRUE) #if you want to see comparisons of conditions for each gene
#te_pluri = te[grepl(paste(pluri, collapse = "|"), te$gene.x),] #only look at genes from a specified signature
dif_cv = dif_cv[complete.cases(dif_cv), ] #keep only rows without NA

#dif_cv$dMean = (dif_cv$mean.y - dif_cv$mean.x)/dif_cv$mean.x
dif_cv$dMean = (dif_cv$mean.y - dif_cv$mean.x)#/dif_cv$mean.x
dif_cv$dCV = dif_cv$cv.y - dif_cv$cv.x
dif_cv$dFano = dif_cv$fano.y - dif_cv$fano.x


dif_cv = dif_cv[order(dif_cv$gene_sig.x, decreasing = F),]

#plot all genes, colored by gene signature
ggplot(dif_cv, aes(x=dMean, y=dCV, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#D3D3D3", "#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + geom_abline(intercept = 0, slope = -1)
#plot only genes in gene signatures
ggplot(dif_cv[dif_cv$gene_sig.x!=0,], aes(x=dMean, y=dCV, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")
ggplot(dif_cv[dif_cv$gene_sig.x=="pluri",], aes(x=dMean, y=dCV, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")

#separate genes by change in CV and mean (increase or decrease in each) in CSRP1 vs LacZ
cv_gene_clusters = list(mposcpos = subset(dif_cv, dMean > 0 & dCV > 0), mnegcpos = subset(dif_cv, dMean < 0 & dCV > 0), mposcneg = subset(dif_cv, dMean > 0 & dCV < 0), mnegcneg = subset(dif_cv, dMean < 0 & dCV < 0))

#plot mean LvS and CV LvS
#blue line is regression, red line is if a 1:1 relationship (slope=1)
ggplot(dif_cv, aes(x=mean.x, y=mean.y)) + geom_point(alpha=0.25) + theme_classic() + geom_segment(aes(x=-5, y=-5, xend=1, yend=1, colour="red")) + stat_smooth(method = "lm")
ggplot(dif_cv, aes(x=cv.x, y=cv.y)) + geom_point(alpha=0.25) + theme_classic() + geom_segment(aes(x=-2.5, y=-2.5, xend=5, yend=5, colour="red")) + stat_smooth(method = "lm")

#plot fano factor LvS
ggplot(dif_cv, aes(x=fano.x, y=fano.y, col=gene_sig.x)) + geom_point() + theme_classic() + scale_color_manual(values=c("#000000", "#D3D3D3", "#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) + geom_segment(aes(x=-2, y=-2, xend=.75, yend=.75, colour="#000000"), linetype="dashed") #+ stat_smooth(method = "lm") #+ geom_density_2d(bins=20, adjust=8)

fano_gene_clusters = list(pos = subset(dif_cv, dFano > 0), neg = subset(dif_cv, dFano < 0))
