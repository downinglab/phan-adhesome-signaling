# Figure 6 Monocle3 Pseudotime Trajectory Analysis
# 10/5/21
# By Andrew Phan

library(Seurat)
library(monocle3)
library(SeuratWrappers)

rm(list=ls()) #clears variables
cat("\014") #clears console

### Load A549 Seurat data


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

DimPlot(combined, reduction = "umap")
#separate conditions
CSRP1seurat = subset(x = combined, subset = orig.ident == "A549_D8_CSRP1_KD")
lacZseurat = subset(x = combined, subset = orig.ident == "A549_D8_LacZ")

seurat_obj = lacZseurat
seurat_obj = CSRP1seurat

cds <- as.cell_data_set(seurat_obj)

cds <- cluster_cells(cds)
#Retrieve the clustering from the Seurat object
#cds@clusters$UMAP$clusters <- seurat_obj@meta.data$seurat_clusters
#names(cds@clusters$UMAP$clusters) <- Cells(seurat_obj) #rename the cells
cds@clusters$UMAP$clusters <- seurat_obj@active.ident

plot_cells(cds, color_cells_by = "cluster")

#cds <- learn_graph(cds, learn_graph_control = list(ncenter = 550))
cds <- learn_graph(cds)

plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster")

cds <- order_cells(cds)

plot_cells(cds, #see results of pseudotime
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

############ a helper function to identify the root principal points: ##############
get_earliest_principal_node <- function(cds){ #finds vertex with most cells belonging to it given condition
  cell_ids <- which(colData(cds)[, "orig.ident"] == "A549_D8_CSRP1_KD")

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
           color_cells_by = "pseudotime",
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
cds_sub@clusters$UMAP$clusters <- CSRP1seurat@active.ident

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
cds_sub2@clusters$UMAP$clusters <- CSRP1seurat@active.ident

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

########### plot gene signatures over pseudotime ##############

#extract pseudotime values and signature scores

CSRP1seurat <- AddMetaData(
  object = CSRP1seurat,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

lacZseurat <- AddMetaData(
  object = lacZseurat,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

#need to update by running this line
seurat_obj = lacZseurat

#df = as.data.frame(seurat_obj@meta.data$orig.ident)
#colnames(df) = "Condition"
#df$Pseudotime = seurat_obj@meta.data$Pseudotime

#reorder and extract pseudotime from df object
df = df_lacZ[order(as.numeric(row.names(df_lacZ))), 1:2]

df$clust = seurat_obj@meta.data$Clusters
df$CSC.SIG = seurat_obj@meta.data$CSC1
df$EMT = seurat_obj@meta.data$EMT.LUAD1
df$PCP.SIG = seurat_obj@meta.data$PCP.PATHWAY_A5491


df <- df[order(df$Pseudotime),] #order based on pseudotime

seurat_obj = CSRP1seurat

#df2 = as.data.frame(seurat_obj@meta.data$orig.ident)
#colnames(df2) = "Condition"
#df2$Pseudotime = seurat_obj@meta.data$Pseudotime

#reorder and extract pseudotime from df object
df2 = df_CSRP1[ order(as.numeric(row.names(df_CSRP1))), 1:2]

df2$clust = seurat_obj@meta.data$Clusters
df2$CSC.SIG = seurat_obj@meta.data$CSC1
df2$EMT = seurat_obj@meta.data$EMT.LUAD1
df2$PCP.SIG = seurat_obj@meta.data$PCP.PATHWAY_A5491

df2 <- df2[order(df2$Pseudotime),] #order based on pseudotime

matplot(x=df$Pseudotime, y=df[,3:ncol(df)], type = c("b")) #matplot can plot multiple columns in one graph and add lines

########### find z-score = (value-mean)/SD ##############
#did not use, optional
dfzscore = df
for(x in 3:ncol(dfzscore)){ #for every signature score (column), calculate the z score
  mean = mean(df[,x])
  sd = sd(df[,x])
  dfzscore[,x] = (dfzscore[,x]-mean)/sd
}

dfzscore = dfzscore[!is.infinite(dfzscore$Pseudotime),]

plot_GeneSigScore(df, y="PCP", breaks = 15, span = 0.4) #CSRP1 CSC.sig zscore
##########

df = df[!is.infinite(df$Pseudotime),] #remove infinite pseudotime cells
df2 = df2[!is.infinite(df2$Pseudotime),] #remove infinite pseudotime cells

sig = "PCP.SIG"
plot_GeneSigScore(df, y=sig, breaks = 10, span = 0.2, ylim = c(-0.02, 0.025))
plot_GeneSigScore(df2, y=sig, breaks = 10, span = 0.2, ylim = c(-0.02, 0.025))


#FIGURES:
par(mfrow=c(1,2))    # set the plotting area into a 1x2 array
plot_GeneSigScore(df_lacZ, y=sig, breaks = 15, span = 0.4, ylim = c(0.05, 0.45)) #lacZ CSC.sig zscore
plot_GeneSigScore(df_CSRP1, y=sig, breaks = 15, span = 0.4, ylim = c(0.05, 0.45)) #lacZ CSC.sig zscore
#save as 10x20in pdf

##################################

FeaturePlot(lacZseurat, c("CSC1", "EMT.LUAD1"), pt.size = 0.1) & scale_color_viridis_c()

x="Pseudotime"
breaks = 20

dfcut = df
dfcut[[x]] <- cut(dfcut[[x]], breaks = breaks) #find breaks
dfagg=aggregate(dfcut, by=list(dfcut[[x]]), mean) #average value per break
plot(dfagg[["PCP.SIG"]] ~ Group.1, data=dfagg,pch=19,cex=0.1)

lw1 <- loess(formula = PCP.SIG ~ Pseudotime, data = df)
plot(PCP.SIG ~ Pseudotime, data=df,pch=19,cex=0.1)
j <- order(df$Pseudotime)
lines(df$Pseudotime[j],lw1$fitted[j],col="red",lwd=3)

lw1 <- loess(formula = PCP.SIG ~ Pseudotime, data = df2)
plot(PCP.SIG ~ Pseudotime, data=df2,pch=19,cex=0.1)
j <- order(df2$Pseudotime)
lines(df2$Pseudotime[j],lw1$fitted[j],col="red",lwd=3)

