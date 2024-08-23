# Figure 1 iGraph Network
# By Andrew Phan

library(igraph)

ipdf=read.delim("adhesome_bulk_seq.txt", sep="\t", header=T)
rownames(ipdf) = ipdf$Name

genelist = ipdf[,c(1,20)] #get table of Ensembl Gene IDs to HGNC symbols

ipdf = ipdf[,c(-1,-20)] #get rid of Ensembl Gene IDs and HGNC symbols

#ipdf is the input dataframe. At this stage, it should only be a matrix of values, the different columns being the conditions, and different rows being the genes.

################# write.table(ipdf,"Log_ipdf.txt", sep="\t")

cm=data.frame(NULL)
#cm is the correlation matrix

df = ipdf[,c(1,4,7,10,13,16)] #NV
#df = ipdf[,c(2,5,8,11,14,17)] #LacZ
#df = ipdf[,c(3,6,9,12,15,18)] #S3

df = as.matrix(df)
thresholdforcorrelation=0.9

for(i in 1:nrow(df))
{
  for(j in 1:nrow(df))
  {
    if(cov(df[i,],df[j,],method="pearson")>thresholdforcorrelation & i!=j)
    {
      cm[i,j]=1
    }
    else if(cov(df[i,],df[j,],method="pearson")<(-thresholdforcorrelation) & i!=j)
    {
      cm[i,j]=-1
    }
    else
    {
      cm[i,j]=0
    }
  }
}

rownames(cm)=genelist[,1]
colnames(cm)=genelist[,1]

#cm2 is the same as cm, except that all negative ones are converted to ones 
cm2=cm
for(i in 1:nrow(cm))
{
  for(j in 1:nrow(cm))
  {
    if(cm2[i,j]==-1)
    {
      cm2[i,j]=1
    }
  }
}

rownames(cm2)=genelist[,1]
colnames(cm2)=genelist[,1]
net=graph_from_adjacency_matrix(as.matrix(cm2))
#network has now been created

V(net)$name = unlist(list(strsplit(gsub(" ", "", toString(genelist[,1])),",")))
E(net)$arrow.mode = 0
#nodes labeled be gene name

#deleting the isolated/disconnected nodes
Isolated = which(degree(net)==0)
net = delete.vertices(net, Isolated)

#labeling the edges based on correlation color
for(i in 1:sum(cm2))   #note sum(cm2) will be equal to number of edges, that can also be gotten by gsize() of net
{
  a=rownames(as.matrix(tail_of(net, E(net)[i])))
  b=rownames(as.matrix(head_of(net, E(net)[i])))
  E(net)[i]$edge.color = ifelse(cm[a,b] ==1 , "green","red")
}

#labeling the nodes based on gene identity
rownames(genelist)=genelist[,1]
for(i in 1:gorder(net))  #gorder() of net gives the number of vertices/nodes, 
{
  a=rownames(as.matrix(V(net)[i]))
  V(net)[i]$vertex.color = ifelse(toString(genelist[a,2]) == "Inflammatory_panel" , "red","gray")
}

#labeling the nodes based on gene identity
genecluster=read.delim("List of gene clusters and category from H3Ac analysis.txt", sep="\t", header=T)
rownames(genecluster)=genecluster[,1]
for(i in 1:gorder(net))  #gorder() of net gives the number of vertices/nodes, 
{
  a=rownames(as.matrix(V(net)[i]))
  V(net)[i]$vertex.frame.color = switch(toString(genecluster[a,3]), "Cluster1"="darkorange", "Cluster2"="blueviolet", "NA" = "black")
}


seed=40
dev.new(width=50,height=50,noRStudioGD = TRUE)
plot.igraph(net, vertex.size=2, edge.width=0.2, vertex.label.cex=0.3, rescale=T, edge.color=E(net)$edge.color, arrow.mode=0, vertex.color=V(net)$vertex.color  , vertex.frame.color=V(net)$vertex.frame.color)

noderanking=sort(rowSums(cm2), decreasing=T)

#write.table(noderanking, "noderanking.txt", sep="\t")

#write.table(as_long_data_frame(net), "net.txt", sep="\t")
