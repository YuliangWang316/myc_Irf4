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

disp_table <- dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
#plot_ordering_genes(monocle_cds)
#plot_pc_variance_explained(monocle_cds, return_all = F)
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, num_dim = 16,
                        reduction_method = 'tSNE', verbose = T)
monocle_cds <- clusterCells(monocle_cds, num_clusters = 4)
#pData(monocle_cds)$Cluster <- pData(monocle_cds)$seurat_cluster
#plot_cell_clusters(monocle_cds, 1, 2, color = "CellType",
#                   markers = c("Myc", "Irf4"))
diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~Cluster",cores = 20) 
ordering_genes <- row.names (subset(diff_test_res, qval < 1E-180))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds, method = 'DDRTree',max_components = 2)
monocle_cds <-orderCells(monocle_cds)

plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters",cell_size = 0.75)

blast_genes <- row.names(subset(fData(monocle_cds),
                                gene_short_name %in% c("Myc", "Irf4")))
plot_genes_jitter(monocle_cds[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)
monocle_expressed_genes <-  row.names(subset(fData(monocle_cds),
                                          num_cells_expressed >1))
monocle_filtered <- monocle_cds[monocle_expressed_genes,]
monocle_genes <- row.names(subset(fData(monocle_filtered),
                             gene_short_name %in% c("Myc", "Irf4")))
cds_subset <- monocle_filtered[monocle_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "State")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")
