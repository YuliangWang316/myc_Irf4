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
ref.data<-read.table("d:/CJJ immunity/GSE139891_RAW/GSE139891.txt",sep = "\t",header = TRUE,check.names = FALSE)
ref.data<-ref.data[!duplicated(ref.data$Gene),]
ref.data<-na.omit(ref.data)
rownames(ref.data)<-ref.data$Gene
ref.data<-ref.data[,-1]
ref.metadata<-read.table("d:/CJJ immunity/GSE139891_RAW/GSE139891metadata.txt",sep = "\t",header = TRUE)
rownames(ref.metadata)<-ref.metadata$ID

ref.data<-ref.data[,rownames(ref.metadata)]
library(Hmisc)
rownames(ref.data)<-capitalize(tolower(rownames(ref.data)))
ref <- CreateSeuratObject(counts = ref.data, project = "ref3k",meta.data = ref.metadata,min.cells = 3)
colnames(ref.metadata)<-c("mouse_ID","Clustername")
ds_seurat<-SummarizedExperiment(assays=list(counts=ref@assays$RNA@counts),colData = ref.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc)<-pred.hesc@listData$labels
pbmc@meta.data$Cluster<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$Cluster)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:12)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:12)
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
saveRDS(pbmc,"e:/CJJGSE139891.rds")
remove(pbmc)
gc()
