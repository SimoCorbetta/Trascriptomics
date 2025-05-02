library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
rownames(sm)
old_names<-rownames(sm)
new_row_names <- sub("_ENSMUSG\\d+.*", "", rownames(sm))

# Assign the new row names to the data frame
rownames(sm) <- new_row_names


p <- CreateSeuratObject(counts =sm, project = "name_of_the_experiment", min.cells =
                                           3, min.features = 200)
#data vengono da esperimento su topi
head(colnames(p)) # ogni colonna rappresenta una cellula
head(rownames(p)) # ogni riga rappresenta un gene
head(p@meta.data) # ottengo un df in cui ogni riga rappresenta una cellula, per essa calcoliamo la library size(ncount_RNA) e il numero di geni espressi nella cellula
p@meta.data$orig.ident<-NULL
grep("^mt-",rownames(p@assays$RNA@counts),value = TRUE)#becco 17 geni mitocondriali
rib <- grep("^Rp[ls]",rownames(p),value = TRUE) # becco 253 proteine ribosomiali
# Calculates the percentage of counts originating from Mitocondrial RNA
p[["percent.mt"]] <- PercentageFeatureSet(p, pattern = "^mt-")
#And from Ribosomial RNA
PercentageFeatureSet(p, features = rib) -> p[["percent.rbp"]]
############################################
# Visualize QC metrics as violin plots - also adding the RPL genes
x11()
VlnPlot(p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4,group.by = NULL)


#without dots:
VlnPlot(p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4, pt.size=0)
x11()
FeatureScatter(p, feature1 = "nCount_RNA", feature2 = "percent.mt")
x11()
plot<-FeatureScatter(p, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot + geom_hline(yintercept = c(200, 2500), linetype = "dashed")
plot
# Display the plot
x11()
FeatureScatter(p, feature1 = "nCount_RNA", feature2 = "percent.rbp")
x11()
FeatureScatter(p, feature1 = "nFeature_RNA", feature2 = "percent.mt")
x11()
FeatureScatter(p, feature1 = "nFeature_RNA", feature2 = "percent.rbp")
x11()
FeatureScatter(p, feature1 = "percent.mt", feature2 = "percent.rbp")
# filtering the dataframe colums(the cells)
#The thresholds are on the number of genes detected (between 200 and 2500), and on the MT DNA/RNA (5%).
p_sub <- subset(p, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)# removed 1095 cellule
p_norm <- NormalizeData(p_sub, normalization.method = "LogNormalize", scale.factor = 10000)
p_norm@assays
head(p_norm@assays$RNA@counts)
#normalized counts are here
head(p_norm@assays$RNA@data)
#we can take a look to the genes that have the highest mean expression across our cells
apply(p_norm@assays$RNA@data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
head(gene.expression, n=50)
#expression of malat1 vs housekeeping gene
x11()
VlnPlot(p_norm, features = c("Malat1","Gapdh"))
#Seurat contains a pre-computed list of cell cycle specific genes
cc.genes.updated.2019
#That can be used to “guess” which CC phase each cell is in
CellCycleScoring(p_norm, s.features = tolower(cc.genes.updated.2019$s.genes),
                 g2m.features = tolower(cc.genes.updated.2019$g2m.genes), set.ident = TRUE) -> p_norm
# finding the principal components
# choosing the 2k genes with highest variance
p_norm <- FindVariableFeatures(p_norm, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(p_norm), 10)
# plot variable features with and without labels
x11()
plot1 <- VariableFeaturePlot(p_norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
# The idea is to shift the expression of each gene, 
#so that the mean expression across cells is 0 and the variance across cells is 1.
all.genes <- rownames(p_norm)
p_norm <- ScaleData(p_norm, features = all.genes)
#PCA
p_norm <- RunPCA(p_norm, features = VariableFeatures(object = p_norm))
# Examine and visualize PCA results a few different ways
print(p_norm[["pca"]], dims = 1:2, nfeatures = 5)


x11()
VizDimLoadings(p_norm, dims = 1:2, reduction = "pca")
#And the projection of the cells in the first two principal components
x11()
DimPlot(p_norm, reduction = "pca") # cells colored based on predicted cell cycle status
# elbow plot needed for understanding the number of PC
x11()
#ElbowPlot(p_norm, ndims=30) # 10 o 11 sembra essere il numero ottimale
ElbowPlot(p_norm, ndims = 30)


#Seurat first constructs a kNN graph based on the euclidean distance in PCA space
p_norm <- FindNeighbors(p_norm, dims = 1:10) # 10 PC
# finding clusters
p_norm <- FindClusters(p_norm, resolution = 0.4)
cluster_assignments <- Idents(p_norm)
# Count the number of objects in each cluster
S <- table(cluster_assignments)
# Seurat authors find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells.
#Optimal resolution often increases for larger datasets.
# Look at cluster IDs of the first 5 cells
head(Idents(p_norm), 5)
#We can plot them in the space of the first two PCA components
x11()
DimPlot(p_norm, reduction = "pca")
# plotting for human visualization using TISNE
x11()
p_norm <- RunTSNE(p_norm, dims=1:10)
DimPlot(p_norm, reduction = "tsne")
#plotting for visualization using UMAP
p_norm <- RunUMAP(p_norm, dims = 1:10)
x11()
DimPlot(p_norm, reduction = "umap")
#We can also check whether some of the critical quality parameters influenced the clustering we got
x11()
VlnPlot(p_norm,features="nCount_RNA")
VlnPlot(p_norm,features="nFeature_RNA")
VlnPlot(p_norm,features="percent.mt")
VlnPlot(p_norm,features="percent.rbp")
# effect of cell cycle on clustering
x11()
p_norm@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")
# marker genes of each clusters(genes higher expressed in that cluster with respect of all others)
#we return only genes "over expressed", found in at least 25% of the cells, and with a logFC threshold of at least 0.25
p_norm.markers <- FindAllMarkers(p_norm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#tabella con 5 gene markers per sample
print(p_norm.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC),n=75) #they are sorted by logFC 

features<-c("Ptgds","Hapln2","Aldoc","Klk6","Cldn5","Tmem141","Ctss","Meg3","Cspg5","Rgs5","Ly6c1","Fyn",
            "Mbp","Fos","Myl9" ) # first gene of the list for each cluster

#we can plot gene marker expression with a heatmap:
x11()
FeaturePlot(p_norm, features = features)
#n single cells grouped by cluster
x11()
p_norm.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(p_norm, features = top10$gene) + NoLegend()

#more gene markers for problematic clusters identified with heatmap
#gene markers for cluster 0
cluster0.markers <- FindMarkers(p_norm, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster0.markers <- cluster0.markers[order(-cluster0.markers$avg_log2FC),]
head(cluster0.markers, n = 10)
#gene markers for cluster 1
cluster1.markers <- FindMarkers(p_norm, ident.1 = 1, min.pct = 0.25, test.use = "wilcox")
cluster1.markers <- cluster1.markers[order(-cluster1.markers$avg_log2FC),]
head(cluster1.markers, n = 10)
#gene markers for cluster 2
cluster2.markers <- FindMarkers(p_norm, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
cluster2.markers <- cluster2.markers[order(-cluster2.markers$avg_log2FC),]
head(cluster2.markers, n = 10)

#gene markers for cluster 3
cluster3.markers <- FindMarkers(p_norm, ident.1 = 3, min.pct = 0.25, test.use = "wilcox")
cluster3.markers <- cluster3.markers[order(-cluster3.markers$avg_log2FC),]
head(cluster3.markers, n = 10)
#gene markers for cluster 4
cluster4.markers <- FindMarkers(p_norm, ident.1 = 4, min.pct = 0.25, test.use = "wilcox")
cluster4.markers <- cluster4.markers[order(-cluster4.markers$avg_log2FC),]
head(cluster4.markers, n = 10)
#gene markers cluster5
cluster5.markers <- FindMarkers(p_norm, ident.1 = 5, min.pct = 0.25, test.use = "wilcox")
cluster5.markers <- cluster5.markers[order(-cluster5.markers$avg_log2FC),]
head(cluster5.markers, n = 10)
#gene markers for cluster7
cluster7.markers <- FindMarkers(p_norm, ident.1 = 7, min.pct = 0.25, test.use = "wilcox")
cluster7.markers <- cluster7.markers[order(-cluster7.markers$avg_log2FC),]
head(cluster7.markers, n = 10)
#gene markers cluster8
cluster8.markers <- FindMarkers(p_norm, ident.1 = 8, min.pct = 0.25, test.use = "wilcox")
cluster8.markers <- cluster8.markers[order(-cluster8.markers$avg_log2FC),]
head(cluster8.markers, n = 15)
#gene markers for cluster 9
cluster9.markers <- FindMarkers(p_norm, ident.1 = 9, min.pct = 0.25, test.use = "wilcox")
cluster9.markers <- cluster9.markers[order(-cluster9.markers$avg_log2FC),]
head(cluster9.markers, n = 15)
#gene markers for cluster 10
cluster10.markers <- FindMarkers(p_norm, ident.1 = 10, min.pct = 0.25, test.use = "wilcox")
cluster10.markers <- cluster10.markers[order(-cluster10.markers$avg_log2FC),]
head(cluster10.markers, n = 10)
#gene markers for cluster 11
cluster11.markers <- FindMarkers(p_norm, ident.1 = 11, min.pct = 0.25, test.use = "wilcox")
cluster11.markers <- cluster11.markers[order(-cluster11.markers$avg_log2FC),]
head(cluster11.markers, n = 20)
#gene markers cluster 12
cluster12.markers <- FindMarkers(p_norm, ident.1 = 12, min.pct = 0.25, test.use = "wilcox")
cluster12.markers <- cluster12.markers[order(-cluster12.markers$avg_log2FC),]
head(cluster12.markers, n = 15)


#gene markers cluster 13
cluster13.markers <- FindMarkers(p_norm, ident.1 = 13, min.pct = 0.25, test.use = "wilcox")
cluster13.markers <- cluster13.markers[order(-cluster13.markers$avg_log2FC),]
head(cluster13.markers, n = 15)
#gene markers cluster 14
cluster14.markers <- FindMarkers(p_norm, ident.1 = 14, min.pct = 0.25, test.use = "wilcox")
cluster14.markers <- cluster14.markers[order(-cluster14.markers$avg_log2FC),]
head(cluster14.markers, n = 15)

#finding marker genes for distinguish between two clusters
cluster01.markers <- FindMarkers(p_norm, ident.1 = 0, ident.2 = 1, min.pct = 0.25, test.use = "wilcox")
cluster01.markers <- cluster01.markers[order(-cluster01.markers$avg_log2FC),]
head(cluster01.markers, n = 10) # genes whose expression is higher in cluster0 to respect to cluster1
#genes whose expression is higher in cluster 3 with respect to cluster1
cluster31.markers <- FindMarkers(p_norm, ident.1 = 3, ident.2 = 1, min.pct = 0.25, test.use = "wilcox")
cluster31.markers <- cluster31.markers[order(-cluster31.markers$avg_log2FC),]
head(cluster31.markers, n = 10)
#genes whose expression is higher in cluster4 with respect to cluster9
cluster49.markers <- FindMarkers(p_norm, ident.1 = 4, ident.2 = 9, min.pct = 0.25, test.use = "wilcox")
cluster49.markers <- cluster49.markers[order(-cluster49.markers$avg_log2FC),]
head(cluster49.markers, n = 10)
#genes whose expression is higher in cluster4 with respect to cluster10
cluster410_9.markers <- FindMarkers(p_norm, ident.1 = c(4,10), ident.2 = 9, min.pct = 0.25, test.use = "wilcox")
cluster410_9.markers <- cluster410_9.markers[order(-cluster410_9.markers$avg_log2FC),]
head(cluster410_9.markers, n = 10)
#genes whose expression is higher in cluster 9 with respect to cluster10
cluster910.markers <- FindMarkers(p_norm, ident.1 = 9, ident.2 = 10, min.pct = 0.25, test.use = "wilcox")
cluster910.markers <- cluster910.markers[order(-cluster910.markers$avg_log2FC),]
head(cluster910.markers, n = 10)
#labeling my clusters
new.cluster.ids <- c("MFOL","MOL" ,"Astrocytes","MOL","endotelial cells","MFOL",
                     "microglia","neurons","Oligodendrocytes precursors","endotelial cells","endotelial cells",
                     "Oligodendrocytes precursors","MFOL","Ependymal cells","smooth muscle cells")
names(new.cluster.ids) <- levels(p_norm)
p_norm <- RenameIdents(p_norm, new.cluster.ids)
x11()
DimPlot(p_norm, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
x11()
VlnPlot(p_norm, features = "Ctss")+ggtitle("Microglia marker")


























