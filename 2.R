library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
options(warn=-1)
pbmc.data <- Read10X(data.dir="I:/MAC/BOAO/BOAO_GCB1_2_filtered_Rawdata/")
pbmc.data <- as.data.frame(pbmc.data)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,Ighd <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB1_filtered", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F)
ElbowPlot(pbmc,ndims = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:16, verbose=F)
pbmc <- FindClusters(pbmc, resolution = 0.6, verbose=F)
pbmc<-RunUMAP(pbmc,dims = 1:16)
pbmc <- RunTSNE(pbmc, dims = 1:16)
pbmc<-subset(pbmc,idents = c("0","1","2","3","4","5","6","7"))
new.cluster.ids <- c("0", "1", "2", "3", "4", "5", 
                     "6", "7")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
new.cluster.ids <- c("DZ", "LZ", "LZ", "DZ", "DZ", "LZ", 
                     "Myc_1", "Myc_2")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$seurat_clusters<-Idents(pbmc)
library(monocle)
data<-as.matrix(pbmc@assays$RNA@counts)
pd <-pbmc@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,expressionFamily=negbinomial())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

pbmc.marker<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0.25,min.pct = 0.1)
ordering_genes <- row.names (subset(pbmc.marker,avg_logFC >0.3 & p_val_adj < 1e-140))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds, method = 'DDRTree',max_components = 2)
monocle_cds <-orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters",cell_size = 0.75)
monocle_expressed_genes <-  row.names(subset(fData(monocle_cds),
                                             num_cells_expressed >= 10))
monocle_filtered <- monocle_cds[monocle_expressed_genes,]
my_genes <- row.names(subset(fData(monocle_filtered),
                             gene_short_name %in% c("Myc", "Irf4", "Kdm6b")))
cds_subset <- monocle_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")