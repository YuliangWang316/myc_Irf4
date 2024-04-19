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

pbmc<-readRDS("e:/GSE237987_final.arid1a.msobj.pseudotime.rds")
pbmc.data<-pbmc@assays$RNA@counts
pbmc.metadata<-pbmc@meta.data
remove(pbmc)
gc()
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "ref3k",meta.data = pbmc.metadata,min.cells = 3)
remove(pbmc.data,pbmc.metadata)
gc()

pbmc2<-readRDS("e:/CJJGSE109732.rds")
Idents(pbmc)<-pbmc$Genotype
pbmc<-subset(pbmc,idents = "WT")
gc()
pbmc.metadata<-pbmc@meta.data[,c(12,26,32,42)]
pbmc.metadata$mouse_ID<-rownames(pbmc.metadata)
pbmc.metadata$Clustername<-pbmc.metadata$predicted.id
library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc@assays$RNA@counts),colData = pbmc.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc2@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc2)<-pred.hesc@listData$labels
pbmc2@meta.data$melnickCluster1<-Idents(pbmc2)
pbmc2<-RunLDA(pbmc2,labels = pbmc2$melnickCluster1,reduction.name = "lda5")
pbmc2<-RunUMAP(pbmc2,reduction = "lda5",reduction.name = "lda5umap",dims = 1:10)
pbmc2<-RunTSNE(pbmc2,reduction = "lda5",reduction.name = "lda5tsne",dims = 1:10)
Idents(pbmc2)<-pbmc2$melnickCluster1
DimPlot(pbmc2,reduction = "lda5umap")
DimPlot(pbmc2,reduction = "lda5tsne")


pbmc.metadata$Clustername<-pbmc.metadata$celltype
library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc@assays$RNA@counts),colData = pbmc.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc2@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc2)<-pred.hesc@listData$labels
pbmc2@meta.data$melnickCluster2<-Idents(pbmc2)
pbmc2<-RunLDA(pbmc2,labels = pbmc2$melnickCluster2,reduction.name = "lda6")
pbmc2<-RunUMAP(pbmc2,reduction = "lda6",reduction.name = "lda6umap",dims = 1:5)
pbmc2<-RunTSNE(pbmc2,reduction = "lda6",reduction.name = "lda6tsne",dims = 1:5)
Idents(pbmc2)<-pbmc2$melnickCluster2
DimPlot(pbmc2,reduction = "lda6umap")
DimPlot(pbmc2,reduction = "lda6tsne")

pbmc.metadata$Clustername<-pbmc.metadata$fine.celltype
library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc@assays$RNA@counts),colData = pbmc.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc2@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc2)<-pred.hesc@listData$labels
pbmc2@meta.data$melnickCluster3<-Idents(pbmc2)
pbmc2<-RunLDA(pbmc2,labels = pbmc2$melnickCluster3,reduction.name = "lda7")
pbmc2<-RunUMAP(pbmc2,reduction = "lda7",reduction.name = "lda7umap",dims = 1:10)
pbmc2<-RunTSNE(pbmc2,reduction = "lda7",reduction.name = "lda7tsne",dims = 1:10)
Idents(pbmc2)<-pbmc2$melnickCluster3
DimPlot(pbmc2,reduction = "lda7umap")
DimPlot(pbmc2,reduction = "lda7tsne")

pbmc.metadata$Clustername<-pbmc.metadata$broad.celltype
library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc@assays$RNA@counts),colData = pbmc.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc2@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc2)<-pred.hesc@listData$labels
pbmc2@meta.data$melnickCluster4<-Idents(pbmc2)
pbmc2<-RunLDA(pbmc2,labels = pbmc2$melnickCluster4,reduction.name = "lda8")
pbmc2<-RunUMAP(pbmc2,reduction = "lda8",reduction.name = "lda8umap",dims = 1:5)
pbmc2<-RunTSNE(pbmc2,reduction = "lda8",reduction.name = "lda8tsne",dims = 1:5)
Idents(pbmc2)<-pbmc2$melnickCluster4
DimPlot(pbmc2,reduction = "lda8umap")
DimPlot(pbmc2,reduction = "lda8tsne")
remove(ds_seurat,pbmc,pbmc.metadata,pred.hesc)
gc()
saveRDS(pbmc2,"e:/CJJGSE237987.rds")

pbmc2.data<-pbmc2@assays$RNA@counts
pbmc2.metadata<-pbmc2@meta.data[,c(5,22,23,25)]
remove(pbmc2)
gc()
pbmc2.metadata$mouse_ID<-rownames(pbmc2.metadata)
pbmc2.metadata$Clustername<-pbmc2.metadata$newcluster
#
Idents(pbmc)<-pbmc$predicted.id
pbmc_new<-subset(pbmc,idents = c("CB_S_G2M","Prememory_Naive","Centrocyte","Transitioning_Sphase","CB_G2M",                                            
                                 "CB_Rec_Sphase","Centroblast","Transitioning_CB_CC","CC_Rec",                
                                 "Recycling","Prememory_Memory"))

ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc2.data),colData = pbmc2.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
pbmc@meta.data$melnick1<-rep("PC",length(rownames(pbmc@meta.data)))
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(rownames(pbmc@meta.data)[j] ==  rownames(pred.hesc)[i]){
      pbmc@meta.data$melnick1[j]<-pred.hesc@listData$labels[i]
    }
  }
}
remove(i,j,pbmc_new,pred.hesc,ds_seurat)
gc()

Idents(pbmc)<-pbmc$celltype
pbmc_new<-subset(pbmc,idents = c("Centroblast","Prememory_Naive","Centrocyte","Transitioning_CB_CC","Recycling",                                            
                                 "Prememory_Memory"))

ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc2.data),colData = pbmc2.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
pbmc@meta.data$melnick2<-rep("PC",length(rownames(pbmc@meta.data)))
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(rownames(pbmc@meta.data)[j] ==  rownames(pred.hesc)[i]){
      pbmc@meta.data$melnick2[j]<-pred.hesc@listData$labels[i]
    }
  }
}
remove(i,j,pbmc_new,pred.hesc,ds_seurat)
gc()

Idents(pbmc)<-pbmc$fine.celltype
pbmc_new<-subset(pbmc,idents = c("Prememory_Naive","Prememory_Memory","Centroblast","CB_G2M","CB_S_G2M",                                            
                                 "Transitioning_Sphase","Transitioning_CB_CC","Centrocyte","CC_Rec",                
                                 "Recycling","CB_Rec_Sphase"))

ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc2.data),colData = pbmc2.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
pbmc@meta.data$melnick3<-rep("PC",length(rownames(pbmc@meta.data)))
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(rownames(pbmc@meta.data)[j] ==  rownames(pred.hesc)[i]){
      pbmc@meta.data$melnick3[j]<-pred.hesc@listData$labels[i]
    }
  }
}
remove(i,j,pbmc_new,pred.hesc,ds_seurat)
gc()

Idents(pbmc)<-pbmc$broad.celltype
pbmc_new<-subset(pbmc,idents = c("Centroblast","PrePB","Centrocyte","Transitioning_CB_CC","Prememory",                                            
                                 "Recycling"))

ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc2.data),colData = pbmc2.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
pbmc@meta.data$melnick4<-rep("PC",length(rownames(pbmc@meta.data)))
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(rownames(pbmc@meta.data)[j] ==  rownames(pred.hesc)[i]){
      pbmc@meta.data$melnick4[j]<-pred.hesc@listData$labels[i]
    }
  }
}
remove(i,j,pbmc_new,pred.hesc,ds_seurat)
gc()


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)

remove(all.genes,pbmc2.data,pbmc2.metadata)
gc()
pbmc<-RunLDA(pbmc,labels = pbmc$melnick1,reduction.name = "lda")
pbmc<-RunLDA(pbmc,labels = pbmc$melnick2,reduction.name = "lda2")
pbmc<-RunLDA(pbmc,labels = pbmc$melnick3,reduction.name = "lda3")
pbmc<-RunLDA(pbmc,labels = pbmc$melnick4,reduction.name = "lda4")

pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:10)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:10)
pbmc<-RunUMAP(pbmc,reduction = "lda2",reduction.name = "lda2umap",dims = 1:10)
pbmc<-RunTSNE(pbmc,reduction = "lda2",reduction.name = "lda2tsne",dims = 1:10)
pbmc<-RunUMAP(pbmc,reduction = "lda3",reduction.name = "lda3umap",dims = 1:10)
pbmc<-RunTSNE(pbmc,reduction = "lda3",reduction.name = "lda3tsne",dims = 1:10)
pbmc<-RunUMAP(pbmc,reduction = "lda4",reduction.name = "lda4umap",dims = 1:10)
pbmc<-RunTSNE(pbmc,reduction = "lda4",reduction.name = "lda4tsne",dims = 1:10)
gc()
saveRDS(pbmc,"E:/CJJGSE237987_new.rds")
# Idents(pbmc)<-pbmc$melnick2
# DimPlot(pbmc,reduction = "ldaumap")
# DimPlot(pbmc,reduction = "ldatsne")
# remove(pbmc2.metadata,pbmc2.data)
# remove(all.genes)
# remove(pbmc.metadata)
# gc()
# new.cluster.ids <- c("PC","prePB","LZ","LZ","DZ","LZ","LZ","DZ","LZ","DZ")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# pbmc@meta.data$newnewcluster<-Idents(pbmc)
# remove(new.cluster.ids)
# 
# gc()
# 
# pbmc<-RunLDA(pbmc,labels = pbmc$newnewcluster,reduction.name = "lda3")
# pbmc<-RunUMAP(pbmc,reduction = "lda3",reduction.name = "lda3umap",dims = 1:3)
# pbmc<-RunTSNE(pbmc,reduction = "lda3",reduction.name = "lda3tsne",dims = 1:3)
# 
# 
# pbmc_new_new<-pbmc
# remove(pbmc)
# gc()
# PC<-subset(pbmc_new_new,idents = c("PC"))
# prePB<-subset(pbmc_new_new,idents = c("prePB"))
# LZ<-subset(pbmc_new_new,idents = c("LZ"))
# DZ<-subset(pbmc_new_new,idents = c("DZ"))
# 
# PC_new<-PC[,sample(1:ncol(PC),100)]
# prePB_new<-prePB[,sample(1:ncol(prePB),100)]
# LZ_new<-LZ[,sample(1:ncol(LZ),100)]
# DZ_new<-DZ[,sample(1:ncol(DZ),100)]
# pbmc_new <- merge(PC_new, y = c(prePB_new, LZ_new,DZ_new), add.cell.ids = c("PC", "prePB", "LZ","DZ"), project = "Total")
# pbmc_new <- NormalizeData(pbmc_new)
# pbmc_new <- FindVariableFeatures(pbmc_new, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(pbmc_new)
# pbmc_new <- ScaleData(pbmc_new, features = all.genes)
# pbmc_new <- RunPCA(pbmc_new, features = VariableFeatures(object = pbmc_new))
# ElbowPlot(pbmc_new)
# pbmc_new <- FindNeighbors(pbmc_new, dims = 1:20)
# pbmc_new <- FindClusters(pbmc_new, resolution = 1.2)
# pbmc_new <- RunUMAP(pbmc_new, dims = 1:20)
# pbmc_new <- RunTSNE(pbmc_new, dims = 1:20)
# remove(all.genes,DZ,DZ_new,LZ,LZ_new,prePB,prePB_new,PC,PC_new)
# gc()
# DimPlot(pbmc_new, reduction = "umap")
# DimPlot(pbmc_new, reduction = "umap",split.by = "group")
# DimPlot(pbmc_new, reduction = "tsne")
# DimPlot(pbmc_new, reduction = "tsne",split.by = "group")
# DotPlot(pbmc_new,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"))
# VlnPlot(pbmc_new,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"),pt.size = 0,sort = TRUE)
# DimPlot(pbmc_new,reduction = "pca")
# DimPlot(pbmc_new,reduction = "pca",split.by = "group")
# pbmc_new <- RunLDA(pbmc_new, labels = pbmc_new$newnewcluster)
# pbmc_new <- RunUMAP(pbmc_new ,reduction = "lda",dims = 1:3,reduction.name = "lda_umap")
# pbmc_new <- RunTSNE(pbmc_new, reduction = "lda",dims = 1:3,reduction.name = "lda_tsne" )
# Idents(pbmc_new)<-pbmc_new$newnewcluster
# DimPlot(pbmc_new,reduction = "lda")
# DimPlot(pbmc_new,reduction = "lda_umap")
# DimPlot(pbmc_new,reduction = "lda_tsne")
# pbmc_new_new<-subset(pbmc_new,idents =c("DN","DP","PC"))
# remove(pbmc)
# gc()
# library(monocle)
# trace('project2MST',edit = T,where = asNamespace("monocle"))
# data<-as.sparse(pbmc_new@assays$RNA@counts)
# pd <-pbmc_new@meta.data
# pd <- new('AnnotatedDataFrame', data = pd)  
# fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
# fd <- new('AnnotatedDataFrame', data = fData)
# monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,
#                               expressionFamily = VGAM::negbinomial.size())
# monocle_cds <- estimateSizeFactors(monocle_cds)
# monocle_cds <- estimateDispersions(monocle_cds)
# monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# remove(data,fd,fData,pbmc_new,pbmc_new_new,pd)
# gc()
# # 
# # 
# # 
# pbmc.marker_1<-FindAllMarkers(pbmc_new_new,only.pos = TRUE,min.pct = 0)
# diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~newnewcluster",cores = 20)
# 
# pbmcmarkers_new<-pbmc.marker_1[which(pbmc.marker_1$p_val_adj < 0.05 ),]
# pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
# remove(data,fd,fData,pd)
# gc()
# # ordering_genes <- diff_test_res$gene_short_name
# ordering_genes <- row.names (subset(diff_test_res, qval < 1e-30))
# monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
# monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2)
# 
# monocle_cds <-orderCells(monocle_cds )
# plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
# 
# 
# plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
# 
# plot_cell_trajectory(monocle_cds, color_by = "newnewcluster",cell_size = 0.75)
# plot_cell_trajectory(monocle_cds, color_by = "newnewcluster",cell_size = 0.75,)+ facet_wrap(~newnewcluster, nrow = 1)
# plotdf=pData(monocle_cds)
# library(ggridges)
# mycolor<-c("#619CFF","#00BA38","#F8766D","Black")
# ggplot(plotdf, aes(x=Pseudotime,y=newnewcluster,fill=newnewcluster))+
#   geom_density_ridges(scale=1) +
#   geom_vline(xintercept = c(5,10),linetype=2)+
#   scale_y_discrete("")+
#   theme_minimal()+
#   theme(
#     panel.grid = element_blank()
#   )+scale_fill_manual(values = mycolor)
# saveRDS(monocle_cds,"d:/CJJ_HUMANpublic_pseudotime/EMTAB9005/Monocle2_Sampling_NEW.rds")