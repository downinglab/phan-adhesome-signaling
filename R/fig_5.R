# Figure 5 Monocle3 Pseudotime Trajectory Analysis
# 10/5/21
# By Andrew Phan

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(viridis)
library(dplyr)
#/opt/R/4.2.0/lib/R/library

rm(list=ls()) #clears variables
cat("\014") #clears console

#README: To plot a module score over pseudotime, need to load in cds objects containing pseudotime info for L and S conditions, respectively.
#Next, load in Seurat objects (can recreate by splitting combined in Seurat analysis R code). Add the pseudotime metadata (and/or desired module score) to split Seurat objects.
#Finally run the script at the end of the Seurat analysis R code.

#import Seurat data
# old combined:
#load("~/DowningLab_Git/AQPhan/scRNA-seq/results_integration_iPSC_D6D9D12D15_LacZ&Shroom3KD.RData")
# combined = AddModuleScore(combined, features = list(c("ACTA2", "TAGLN", "MYL9", "ID3")), name = "SMC.SIG")
# combined = AddModuleScore(combined, features = list(c("TFAP2A", "SOX9", "HES1", "HOXA5", "NES", "NGFR")), name = "NCC.SIG")
# combined = AddModuleScore(combined, features = list(c("NEUROG3", "NEUROD1", "NEUROG1", "POU3F2", "NHLH1")), name = "NEURAL.SIG")
# combined = AddModuleScore(combined, features = list(c("POU5F1", "NANOG", "CDH1", "DNMT3A", "UTF1", "LIN28A", "NODAL", "LEFTY2")), name = "iPSC.SIG")
# combined = AddModuleScore(combined, features = list(c("KRT8", "KRT18", "FBP1", "EPCAM")), name = "EPI.SIG")
# combined = AddModuleScore(combined, features = list(c("UBE2C", "CENPF", "BIRC5", "TOP2A", "EED", "CCNB1", "ID1", "ETS2")), name = "TROPH.SIG")
# combined = AddModuleScore(combined, features = list(c("CAV1", "COL5A2", "COL3A1", "TWIST1")), name = "FIBR.SIG")
# combined = AddModuleScore(combined, features = list(c("CRABP1", "CRYAB", "NBEAL1", "S100A6", "PMEPA1")), name = "NEUREPI.SIG")
#
# combined <- FindClusters(combined, resolution = 0.44) #res of 0.43 or 0.44 makes 9 clusters. 0.44 a little cleaner

#plot a single gene signature score over pseudotime with a loess
#Input: dataframe with pseudotime and gene signature values; the column name for the x axis; the column name for the y axis; the number of breaks for the bins on the x axis [NOTE: TOO MANY BREAKS WILL PREVENT LOESS LINE FROM PLOTTING CORRECTLY]; a loess span value [Default=0.2]
#Output: A plot
plot_GeneSigScore <- function(df, x = "Pseudotime", y, breaks=50, span=0.2, ylim = NA){
  dfcut = df
  dfcut[[x]] <- cut(dfcut[[x]], breaks = breaks) #find breaks
  dfagg=aggregate(dfcut, by=list(dfcut[[x]]), mean) #average value per break

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

load("/home/data/aqphan/RStudio/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/Seurat/combined_10clust+signatures.rdata") #load seurat object

#add PCP and EMT module scores
combined = AddModuleScore(combined, features = list(c("WNT5A", "VANGL2", "FZD3", "DVL2", "CELSR1", "PRICKLE1", "FZD7", "WNT11", "DVL1")), name = "PCP.PATHWAY")
combined = AddModuleScore(combined, features = list(c("BMP4", "SOX9", "SNAI2", "CDH2", "CDH11")), name = "EMT")


DimPlot(combined, reduction = "umap")
#separate conditions
s3seurat = subset(x = combined, subset = Condition == "Shroom3")
lacZseurat = subset(x = combined, subset = Condition == "LacZ")

#add finalized clusters as metadata to each object
s3seurat <- AddMetaData(
  object = s3seurat,
  metadata = s3seurat@active.ident,
  col.name = "active.clusters"
)

lacZseurat <- AddMetaData(
  object = lacZseurat,
  metadata = lacZseurat@active.ident,
  col.name = "active.clusters"
)


seurat_obj = lacZseurat
seurat_obj = s3seurat

cds <- as.cell_data_set(seurat_obj)

#cds <- preprocess_cds(cds, num_dim = 50)
#cds <- reduce_dimension(cds)

cds <- cluster_cells(cds)
#Retrieve the clustering from the Seurat object
#cds@clusters$UMAP$clusters <- seurat_obj@meta.data$seurat_clusters
#names(cds@clusters$UMAP$clusters) <- Cells(seurat_obj) #rename the cells
cds@clusters$UMAP$clusters <- seurat_obj@active.ident

plot_cells(cds, color_cells_by = "cluster")

#cds <- learn_graph(cds, learn_graph_control = list(ncenter = 550))
cds <- learn_graph(cds)

plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster")

plot_cells(cds, #see timepoint taken to help determine root
           color_cells_by = "Time",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

cds <- order_cells(cds)

plot_cells(cds, #see results of pseudotime
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

############ a helper function to identify the root principal points: ##############
get_earliest_principal_node <- function(cds, time_bin="D6"){
  cell_ids <- which(colData(cds)[, "Time"] == time_bin)

  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

#################

traj_lacZ_TROPH <- choose_graph_segments(cds) #subset cells from a trajectory

cds_sub = traj_s3_NCC

cds_sub <- preprocess_cds(cds_sub, num_dim = 50)

cds_sub <- reduce_dimension(cds_sub, reduction_method = "UMAP")

cds_sub <- cluster_cells(cds_sub)
cds_sub@clusters$UMAP$clusters <- s3seurat@active.ident

plot_cells(cds_sub, group_cells_by = "cluster")

cds_sub <- learn_graph(cds_sub)

plot_cells(cds_sub,
           color_cells_by = "Time",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

cds_sub <- order_cells(cds_sub)

plot_cells(cds_sub, color_cells_by = "pseudotime")


cds_sub2 = cds_sub
cds_sub2@clusters$UMAP$clusters <- s3seurat@active.ident

plot_cells(cds_sub2, group_cells_by = "cluster")


cds_subset <- choose_cells(cds)
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=64)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < .1))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.1)

agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module,
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

###############
############ scRNA-seq for Monocle3 example: https://github.com/mamtagiri/SingleCellSequencinginR/blob/master/monoclewithSeuratobject.R #################
s3seurat = subset(x = combined, subset = Condition == "Shroom3")
lacZseurat = subset(x = combined, subset = Condition == "LacZ")

seuratobj = s3seurat

data <- as(as.matrix(seuratobj@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = seuratobj@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds2 <- newCellDataSet(data,
                       phenoData  = pd,
                       featureData = fd,expressionFamily = negbinomial.size())

#cds2@normalized_data_projection <- Embeddings(object = combined, reduction = "pca")

cds2 <- estimateSizeFactors(cds2)
#cds2 <- estimateDispersions(cds2)

cds2 <- reduceDimension(cds2,reduction_method = "DDRTree", pseudo_expr = 1)
cds2 <- orderCells(cds2)
plot_cell_trajectory(lacZCDS, color_by = "Pseudotime")
#can see phenotype data here: lacZCDS@phenoData@data$

disp_table <- dispersionTable(cds2)
clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds2 <- setOrderingFilter(cds2, clustering_genes$gene_id)

cds2 <- reduceDimension(cds2, num_dim = 40, reduction_method = 'tSNE')
cds2 <- clusterCells(cds2, method = "louvain")
plot_cell_clusters(cds2, cell_size = 0.5) +
  theme(legend.position = "none") +
  labs(x = "tSNE1", y = "tSNE2")

disp_table <- dispersionTable(cds2)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.5 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id

cds2 <- setOrderingFilter(cds2, ordering_genes = ordering_genes)
plot_ordering_genes(cds2)

laczCDS <- clusterCells(lacZCDS)

plot_pc_variance_explained(cds2, return_all = F)

cds2 <- clusterCells(cds2, num_clusters = 10, frequency_thresh = 0.1, cell_type_hierarchy = cth)
plot_cells(cds, color_cells_by = "partition")

cds2 <- reduceDimension(cds2, max_components = 2, num_dim = 4,
                        reduction_method = 'tSNE', verbose = TRUE)

########### plot gene signatures over pseudotime ##############

#extract pseudotime values and signature scores. change the object name, "object", "metadata"
s3seurat <- AddMetaData(
  object = s3seurat,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

lacZseurat <- AddMetaData(
  object = lacZseurat,
  metadata = lacZcds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

seurat_obj = lacZseurat

df = as.data.frame(seurat_obj@meta.data$seurat_clusters)
colnames(df) = "seurat.clusters"
df$Pseudotime = pseudotime(lacZcds)

df$PCP.SIG = seurat_obj@meta.data$PCP.PATHWAY1
df$EMT = seurat_obj@meta.data$EMT1
df$Pluripotency = seurat_obj@meta.data$Pluripotency

df$iPSC.SIG = seurat_obj@meta.data$iPSC.SIG1
df$SMC.SIG = seurat_obj@meta.data$SMC.SIG1
df$NCC.SIG = seurat_obj@meta.data$NCC.SIG1
df$Neural.SIG = seurat_obj@meta.data$NEURAL.SIG1
df$EPI.SIG = seurat_obj@meta.data$EPI.SIG1
df$TROPH.SIG = seurat_obj@meta.data$TROPH.SIG1
df$FIBR.SIG = seurat_obj@meta.data$FIBR.SIG1
df$NEUREPI.SIG = seurat_obj@meta.data$NEUREPI.SIG1

df <- df[order(df$Pseudotime),] #order based on pseudotime

matplot(x=df$Pseudotime, y=df[,3:5], type = c("b")) #matplot can plot multiple columns in one graph and add lines

seurat_obj = s3seurat

df2 = as.data.frame(seurat_obj@meta.data$seurat_clusters)
colnames(df) = "seurat.clusters"
df2$Pseudotime = pseudotime(s3cds)

df2$PCP.SIG = seurat_obj@meta.data$PCP.PATHWAY1
df2$EMT = seurat_obj@meta.data$EMT1
df2$Pluripotency = seurat_obj@meta.data$Pluripotency

df2$iPSC.SIG = seurat_obj@meta.data$iPSC.SIG1
df2$SMC.SIG = seurat_obj@meta.data$SMC.SIG1
df2$NCC.SIG = seurat_obj@meta.data$NCC.SIG1
df2$Neural.SIG = seurat_obj@meta.data$NEURAL.SIG1
df2$EPI.SIG = seurat_obj@meta.data$EPI.SIG1
df2$TROPH.SIG = seurat_obj@meta.data$TROPH.SIG1
df2$FIBR.SIG = seurat_obj@meta.data$FIBR.SIG1
df2$NEUREPI.SIG = seurat_obj@meta.data$NEUREPI.SIG1

df2 <- df2[order(df2$Pseudotime),] #order based on pseudotime

############ find z-score = (value-mean)/SD ###############
#optional, not used
dfzscore = df
for(x in 3:ncol(dfzscore)){ #for every signature score (column), calculate the z score
  mean = mean(df[,x])
  sd = sd(df[,x])
  dfzscore[,x] = (dfzscore[,x]-mean)/sd
}
##############

df = df[!is.infinite(df$Pseudotime),] #remove infinite pseudotime cells
df2 = df2[!is.infinite(df2$Pseudotime),] #remove infinite pseudotime cells

sig = "PCP"

par(mfrow=c(1,2))    # set the plotting area into a 1x2 array
plot_GeneSigScore(df, y=sig, breaks = 50, span = 0.15, ylim = c(-0.05, 0.09))
plot_GeneSigScore(df2, y=sig, breaks = 50, span = 0.15, ylim = c(-0.05, 0.09))

#FIGURES:
par(mfrow=c(1,2))    # set the plotting area into a 1x2 array
plot_GeneSigScore(df_lacZ, y=sig, breaks = 50, span = 0.4, ylim = c(0.05, 0.45)) #lacZ CSC.sig zscore
plot_GeneSigScore(df_CSRP1, y=sig, breaks = 50, span = 0.4, ylim = c(0.05, 0.45)) #lacZ CSC.sig zscore
#save as 10x20in pdf

#next need to bin and aggregate, then plot values and plot lines
dfcut = df
dfcut$Pseudotime <- cut(dfcut$Pseudotime, breaks = 50)
dfagg=aggregate(dfcut, by=list(dfcut$Pseudotime), mean)

#matplot(dfagg[,4:ncol(dfagg)], x=dfagg$Group.1, type = "b", pch = 15:19)

lw1 <- loess(Neural.SIG ~ Pseudotime, data=df, span = span)
plot(Neural.SIG ~ Pseudotime, data=df,pch=19,cex=0.1)
j <- order(df$Pseudotime)
lines(df$Pseudotime[j],lw1$fitted[j],col="red",lwd=3)

df = df[!is.infinite(df$Pseudotime),]
dfzscore = dfzscore[!is.infinite(dfzscore$Pseudotime),]

plot_GeneSigScore(lacZdf, y="FIBR.SIG", breaks = 40, ylim = c(-0.7, .6))
plot_GeneSigScore(df, y="FIBR.SIG", breaks = 40, ylim = c(-0.7, .6))

plot(df$Pseudotime, df$iPSC.SIG)
plot(lacZdf$Pseudotime, lacZdf$iPSC.SIG)


ggplot(data=df, aes(x=Pseudotime, y=iPSC.Sig, col=seurat.clusters)) +
  geom_point()

lw1 <- loess(df[,3] ~ df[,2], data=df, span = span)
plot(iPSC.Sig ~ Group.1, data=dfagg,pch=19,cex=0.1)
j <- order(df$Pseudotime)
lines(df$Pseudotime[j],lw1$fitted[j],col="red",lwd=3)














######## ATTEMPT TO RUN DIFFERENTIAL GENE EXPRESSION ANALYSES IN MONOCLE3; DOES NOT CARRY OVER ALL NECESSARY DATA WITH SEURATWRAPPERS #############
#try adding pseudotime to column in seurat. then remake cds and try running
lacZseurat <- AddMetaData(
  object = lacZseurat,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

s3seurat <- AddMetaData(
  object = s3seurat,
  metadata = s3cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

cds <- as.cell_data_set(s3seurat)

#will correct missing data from seurat wrappers as.cell_data_set
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(s3seurat[["RNA"]])


#cds@colData$Pluripotency

gene_fits <- fit_models(cds, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
time_terms <- fit_coefs %>% filter(term == "pseudotime")


#Find differentially expressed genes according to pseudotime
gene_fits <- fit_models(cds, model_formula_str = "~Pluripotency")
gene_fitss3 <- fit_models(s3cds, model_formula_str = "~pseudotime")

#extract a table of coefficients from each model
fit_coefs <- coefficient_table(gene_fits)
fit_coefss3 <- coefficient_table(gene_fitss3)

#extract only pseudotime models
time_terms <- fit_coefs %>% filter(term == "pseudotime")
time_termss3 <- fit_coefss3 %>% filter(term == "pseudotime")

#filter out genes significantly differing over pseudotime
time_terms %>% filter (q_value < 0.05) %>% select(gene_id, term, q_value, estimate)


test <- evaluate_fits(gene_fits)

#try graph_test
colData(s3cds)$Cluster = s3cds@clusters$UMAP$clusters

ncctraj_cds = s3cds[,colData(s3cds)$Cluster %in% c(1,8,5,4)]
plot_cells(ncctraj_cds, color_cells_by="partition")

pr_graph_test_res <- graph_test(ncctraj_cds, neighbor_graph="knn", cores=64)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(ncctraj_cds[pr_deg_ids,], resolution=1e-2)

cell_group_df <- tibble::tibble(cell=row.names(colData(ncctraj_cds)),
                                cell_group=partitions(s3cds)[colnames(ncctraj_cds)])
agg_mat <- aggregate_gene_expression(ncctraj_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)

"POU5F1" %in% row.names(rowData(ncctraj_cds))

ncctraj_cds <- estimate_size_factors(ncctraj_cds)
ncctraj_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(s3seurat[["RNA"]])



plot_cells(ncctraj_cds,
           genes=gene_module_df,
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)

plot_cells(ncctraj_cds,
           genes=gene_module_df %>% filter(module %in% c(8)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)

##################################

#https://satijalab.org/signac/articles/monocle.html
#how to add pseudotime metadata to seurat object
bone <- AddMetaData(
  object = bone,
  metadata = erythroid.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Erythroid"
)

FeaturePlot(bone, c("Erythroid", "Lymphoid"), pt.size = 0.1) & scale_color_viridis_c()

# s3cds = cds
# s3df = df
# s3dfzscore = dfzscore
# save(s3cds, s3df, s3dfzscore, file = "~/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/Monocle3/Saved_Objs/s3_monocle+dfs.RDATA")
#
# lacZcds = cds
# lacZdf = df
# lacZdfzscore = dfzscore
# save(lacZcds, lacZdf, lacZdfzscore, file = "~/DowningLab_Git/AQPhan/scRNA-seq/iPSC_Reprog/Monocle3/Saved_Objs/lacZ_monocle+dfs.RDATA")

