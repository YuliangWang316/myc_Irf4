library(patchwork)
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
options(warn=-1)
set.seed(1)
library(Hmisc)

GCB1.data<-Read10X("D:/GCB1_2.2.filtered_feature_bc_matrix/")
GCB2.data<-Read10X("D:/GCB2_2.2.filtered_feature_bc_matrix/")

GCB1.data<-as.data.frame(GCB1.data)
GCB2.data<-as.data.frame(GCB2.data)

# GCB1.metadata<-read.table("d:/GCB1metadata.txt",sep = "\t",header = TRUE,row.names = 1)
# GCB2.metadata<-read.table("d:/GCB2metadata.txt",sep = "\t",header = TRUE,row.names = 1)
# 
# GCB1.data<-GCB1.data[,rownames(GCB1.metadata)]
# GCB2.data<-GCB2.data[,rownames(GCB2.metadata)]

for (i in 1:length(colnames(GCB1.data))) {
  colnames(GCB1.data)[i] <- paste(colnames(GCB1.data)[i],"GCB1",i,sep = "-")  
}

for (i in 1:length(colnames(GCB2.data))) {
  colnames(GCB2.data)[i] <- paste(colnames(GCB2.data)[i],"GCB2",i,sep = "-")  
}

# rownames(GCB1.metadata)<-colnames(GCB1.data)
# rownames(GCB2.metadata)<-colnames(GCB2.data)
# 
# GCB1.metadata$group<-rep("GCB1",length(colnames(GCB1.data)))
# GCB2.metadata$group<-rep("GCB2",length(colnames(GCB2.data)))

# pbmc.metadata<-rbind(GCB1.metadata,GCB2.metadata)

pbmc.data<-cbind(GCB1.data,GCB2.data)
pbmc.data <- as.data.frame(pbmc.data)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,Ighd <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB",min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 30)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F,npcs = 100)
ElbowPlot(pbmc,ndims = 100)

library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ref.data<-read.table("d:/CJJ immunity/GSE188617_counts_info.tsv/counts_info.tsv",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)

ref.metadata<-read.table("d:/CJJ immunity/GSE188617_counts_info.tsv/vdj_info_new.txt",sep = "\t",header = TRUE,row.names = 1)
ref.metadata$ID<-rownames(ref.metadata)

ref.data<-ref.data[,rownames(ref.metadata)]
library(Hmisc)
rownames(ref.data)<-capitalize(tolower(rownames(ref.data)))
ref <- CreateSeuratObject(counts = ref.data, project = "ref3k",meta.data = ref.metadata,min.cells = 3)

ds_seurat<-SummarizedExperiment(assays=list(counts=ref@assays$RNA@counts),colData = ref.metadata[,c("ID","Cluster")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc@assays$RNA@data, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Cluster)

Idents(pbmc)<-pred.hesc@listData$labels
pbmc@meta.data$Cluster<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$Cluster)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:14)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:14)
DimPlot(pbmc,reduction = "ldaumap")
DimPlot(pbmc,reduction = "ldatsne")
VlnPlot(pbmc,features = c("Irf4","Gls","Kdm6b","Ly75"),sort = TRUE)
DotPlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"))
remove(ds_seurat,GCB1.data,GCB2.data,pbmc.data,pbmc.data_t,pred.hesc,ref.data,ref.data2,ref.metadata,ref.metadata2,all.genes,i,new.cluster.ids)
gc()

pbmc2<-readRDS("e:/CJJGSE109732.rds")
pbmc2.data<-pbmc2@assays$RNA@counts
pbmc2.metadata<-pbmc2@meta.data[,c(5,22,23,25)]
remove(pbmc2)
gc()
saveRDS(pbmc,"e:/CJJGSE188617.rds")
remove(pbmc)
gc()
Idents(ref)<-ref$Cluster
ref_new<-subset(ref,idents = c("DZ-4","DZ-3","INT-4","LZ-2","INT-1","DZ-1","LZ-1","INT-2","DZ-2","PreM","Undetermined","INT-3","FCRL2/3","LZ-3"))
pbmc2.metadata$mouse_ID<-rownames(pbmc2.metadata)
pbmc2.metadata$Clustername<-pbmc2.metadata$newcluster
rownames(pbmc2.data)<-toupper(rownames(pbmc2.data)) ##big upper   
ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc2.data),colData = pbmc2.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = ref_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
ref@meta.data$newcluster<-rep("PC",length(rownames(ref@meta.data)))
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(ref@meta.data))) {
    if(rownames(ref@meta.data)[j] ==  rownames(pred.hesc)[i]){
      ref@meta.data$newcluster[j]<-pred.hesc@listData$labels[i]
    }
  }
}
remove(i,j,ref_new,pred.hesc,ds_seurat)
gc()
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genes <- rownames(ref)
ref <- ScaleData(ref, features = all.genes, verbose=F)

ref<-RunLDA(ref,labels = ref$newcluster)
Idents(ref)<-ref$newcluster
ref<-RunUMAP(ref,reduction = "lda",reduction.name = "ldaumap",dims = 1:10)
ref<-RunTSNE(ref,reduction = "lda",reduction.name = "ldatsne",dims = 1:10)
Idents(ref)<-ref$newcluster
DimPlot(ref,reduction = "ldaumap")
DimPlot(ref,reduction = "ldatsne")
remove(all.genes)
new.cluster.ids <- c("LZ","prePB","PC","DZ","LZ","LZ","DZ","DZ","LZ","LZ","LZ")
names(new.cluster.ids) <- levels(ref)
ref <- RenameIdents(ref, new.cluster.ids)
ref@meta.data$newnewcluster<-Idents(ref)
remove(new.cluster.ids)
remove(pbmc2.data,pbmc2.metadata,ref.data,ref.metadata)
gc()

saveRDS(ref,"e:/GSE188617smallformonocle.rds")
saveRDS(ref,"e:/GSE188617bigformonocle.rds")
ref<-RunLDA(ref,labels = ref$newnewcluster,reduction.name = "lda2")
ref<-RunUMAP(ref,reduction = "lda2",reduction.name = "lda2umap",dims = 1:3)
ref<-RunTSNE(ref,reduction = "lda2",reduction.name = "lda2tsne",dims = 1:3)

pbmc_new_new<-ref
remove(ref)
gc()
PC<-subset(pbmc_new_new,idents = c("PC"))
prePB<-subset(pbmc_new_new,idents = c("prePB"))
LZ<-subset(pbmc_new_new,idents = c("LZ"))
DZ<-subset(pbmc_new_new,idents = c("DZ"))

PC_new<-PC[,sample(1:ncol(PC),100)]
prePB_new<-prePB[,sample(1:ncol(prePB),100)]
LZ_new<-LZ[,sample(1:ncol(LZ),100)]
DZ_new<-DZ[,sample(1:ncol(DZ),100)]
pbmc_new <- merge(PC_new, y = c(prePB_new, LZ_new,DZ_new), add.cell.ids = c("PC", "prePB", "LZ","DZ"), project = "Total")
pbmc_new <- NormalizeData(pbmc_new)
pbmc_new <- FindVariableFeatures(pbmc_new, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_new)
pbmc_new <- ScaleData(pbmc_new, features = all.genes)
pbmc_new <- RunPCA(pbmc_new, features = VariableFeatures(object = pbmc_new))
ElbowPlot(pbmc_new)
pbmc_new <- FindNeighbors(pbmc_new, dims = 1:20)
pbmc_new <- FindClusters(pbmc_new, resolution = 1.2)
pbmc_new <- RunUMAP(pbmc_new, dims = 1:20)
pbmc_new <- RunTSNE(pbmc_new, dims = 1:20)
remove(all.genes,DZ,DZ_new,LZ,LZ_new,prePB,prePB_new,PC,PC_new)
gc()
DimPlot(pbmc_new, reduction = "umap")
DimPlot(pbmc_new, reduction = "umap",split.by = "group")
DimPlot(pbmc_new, reduction = "tsne")
DimPlot(pbmc_new, reduction = "tsne",split.by = "group")
DotPlot(pbmc_new,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"))
VlnPlot(pbmc_new,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"),pt.size = 0,sort = TRUE)
DimPlot(pbmc_new,reduction = "pca")
DimPlot(pbmc_new,reduction = "pca",split.by = "group")
pbmc_new <- RunLDA(pbmc_new, labels = pbmc_new$newnewcluster)
pbmc_new <- RunUMAP(pbmc_new ,reduction = "lda",dims = 1:3,reduction.name = "lda_umap")
pbmc_new <- RunTSNE(pbmc_new, reduction = "lda",dims = 1:3,reduction.name = "lda_tsne" )
Idents(pbmc_new)<-pbmc_new$newnewcluster
DimPlot(pbmc_new,reduction = "lda")
DimPlot(pbmc_new,reduction = "lda_umap")
DimPlot(pbmc_new,reduction = "lda_tsne")
pbmc_new_new<-subset(pbmc_new,idents =c("DN","DP","PC"))
remove(pbmc)
gc()
gc()
library(monocle)
trace('project2MST',edit = T,where = asNamespace("monocle"))
data<-as.sparse(pbmc_new@assays$RNA@counts)
pd <-pbmc_new@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,
                              expressionFamily = VGAM::negbinomial.size())
remove(pbmc_new_new,data,fd,fData,pd,pbmc_new)
gc()
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 
# 
# 
pbmc.marker_1<-FindAllMarkers(pbmc_new_new,only.pos = TRUE,min.pct = 0)
diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~newnewcluster",cores = 20)

pbmcmarkers_new<-pbmc.marker_1[which(pbmc.marker_1$p_val_adj < 0.05 ),]
pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
# ordering_genes <- diff_test_res$gene_short_name
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-30))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2)

monocle_cds <-orderCells(monocle_cds )
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "newnewcluster",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "newnewcluster",cell_size = 0.75,)+ facet_wrap(~newnewcluster, nrow = 1)
plotdf=pData(monocle_cds)
library(ggridges)
mycolor<-c("#619CFF","#00BA38","#F8766D","Black")
ggplot(plotdf, aes(x=Pseudotime,y=newnewcluster,fill=newnewcluster))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )+scale_fill_manual(values = mycolor)
saveRDS(monocle_cds,"d:/CJJ_HUMANpublic_pseudotime/GSE188617/monocle2_sampling.rds")
