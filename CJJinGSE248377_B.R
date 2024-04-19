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
remove(all.genes,i,GCB1.data,GCB2.data,pbmc.data,pbmc.data_t)
gc()
library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ref.data<-Read10X("d:/CJJ immunity/GSE248377_RAW/GSM7099104_GEX_Henry_2_Processed/filtered_feature_bc_matrix/")
ref1.data<-Read10X("d:/CJJ immunity/GSE248377_RAW/GSM7099103_Ab_Henry_2_Processed/filtered_feature_bc_matrix/")
ref2.data<-Read10X("d:/CJJ immunity/GSE248377_RAW/GSM7099101_GEX_Henry_1_Processed/filtered_feature_bc_matrix/")
ref3.data<-Read10X("d:/CJJ immunity/GSE248377_RAW/GSM7099100_Ab_Henry_1_Processed/filtered_feature_bc_matrix/")
ref4.data<-Read10X("d:/CJJ immunity/Sample1_MD4_CSP9_Processed/Transcriptome_and_Antibody_Capture/sample_filtered_feature_bc_matrix/")
ref5.data<-Read10X("d:/CJJ immunity/Sample2_MD4_CSP27_Processed/Transcriptome_and_Antibody_Capture/sample_filtered_feature_bc_matrix/")
ref6.data<-Read10X("d:/CJJ immunity/Sample3_WT_CSP9_Processed/Transcriptome_and_Antibody_Capture/sample_filtered_feature_bc_matrix/")
ref7.data<-Read10X("d:/CJJ immunity/Sample4_WT_CSP27_Processed/Transcriptome_and_Antibody_Capture/sample_filtered_feature_bc_matrix/")

ref1.data<-ref1.data$`Gene Expression`
ref3.data<-ref3.data$`Gene Expression`
ref4.data<-ref4.data$`Gene Expression`
ref5.data<-ref5.data$`Gene Expression`
ref6.data<-ref6.data$`Gene Expression`
ref7.data<-ref7.data$`Gene Expression`

ref.metadata<-read.table("e:/ScienceDirect_files_18Apr2024_01-31-16.072/refmetadata.txt",sep = "\t",header = TRUE)
ref1.metadata<-read.table("e:/ScienceDirect_files_18Apr2024_01-31-16.072/ref1metadata.txt",sep = "\t",header = TRUE,row.names = 1)
ref2.metadata<-read.table("e:/ScienceDirect_files_18Apr2024_01-31-16.072/ref2metadata.txt",sep = "\t",header = TRUE)
ref3.metadata<-read.table("e:/ScienceDirect_files_18Apr2024_01-31-16.072/ref3metadata.txt",sep = "\t",header = TRUE,row.names = 1)
ref4.metadata<-read.table("e:/ScienceDirect_files_18Apr2024_01-31-16.072/ref4metadata.txt",sep = "\t",header = TRUE,row.names = 1)
ref5.metadata<-read.table("e:/ScienceDirect_files_18Apr2024_01-31-16.072/ref5metadata.txt",sep = "\t",header = TRUE,row.names = 1)
ref6.metadata<-read.table("e:/ScienceDirect_files_18Apr2024_01-31-16.072/ref6metadata.txt",sep = "\t",header = TRUE,row.names = 1)
ref7.metadata<-read.table("e:/ScienceDirect_files_18Apr2024_01-31-16.072/ref7metadata.txt",sep = "\t",header = TRUE,row.names = 1)

remove(ref.metadata,ref2.metadata)
remove(ref.data,ref2.data)
gc()

ref1.data<-ref1.data[,rownames(ref1.metadata)]
ref3.data<-ref3.data[,rownames(ref3.metadata)]
ref4.data<-ref4.data[,rownames(ref4.metadata)]
a<-intersect(colnames(ref5.data),rownames(ref5.metadata))
ref5.data<-ref5.data[,a]
ref5.metadata<-ref5.metadata[a,]
ref5.data<-ref5.data[,rownames(ref5.metadata)]
a<-intersect(colnames(ref6.data),rownames(ref6.metadata))
ref6.data<-ref6.data[,a]
ref6.metadata<-ref6.metadata[a,]
ref6.data<-ref6.data[,rownames(ref6.metadata)]
a<-intersect(colnames(ref7.data),rownames(ref7.metadata))
ref7.data<-ref7.data[,a]
ref7.metadata<-ref7.metadata[a,]
ref7.data<-ref7.data[,rownames(ref7.metadata)]
remove(a)
gc()

ref1<-CreateSeuratObject(counts = ref1.data,meta.data = ref1.metadata,min.cells = 3,min.features = 200)
ref3<-CreateSeuratObject(counts = ref3.data,meta.data = ref3.metadata,min.cells = 3,min.features = 200)
ref4<-CreateSeuratObject(counts = ref4.data,meta.data = ref4.metadata,min.cells = 3,min.features = 200)
ref5<-CreateSeuratObject(counts = ref5.data,meta.data = ref5.metadata,min.cells = 3,min.features = 200)
ref6<-CreateSeuratObject(counts = ref6.data,meta.data = ref6.metadata,min.cells = 3,min.features = 200)
ref7<-CreateSeuratObject(counts = ref7.data,meta.data = ref7.metadata,min.cells = 3,min.features = 200)

ref<-merge(ref4, y = c(ref5,ref6,ref7), add.cell.ids = c("ref4","ref5","ref6","ref7"), project = "Total")
remove(ref1.metadata,ref3.metadata,ref4.metadata,ref5.metadata,ref6.metadata,ref7.metadata)
remove(ref1.data,ref3.data,ref4.data,ref5.data,ref6.data,ref7.data)
remove(ref1,ref3,ref4,ref5,ref6,ref7)
gc()
Idents(ref)<-ref$Clustername
ref_new<-subset(ref,idents = c("DZd","DZp","LZ3","LZ2","LZ1","PC"))
ref_new_new<-subset(ref_new,idents = c("DZd","DZp","LZ3","LZ2","LZ1"))

ref_new.metadata<-ref_new@meta.data
ref_new.metadata$barcodes<-rownames(ref_new.metadata)
ds_seurat<-SummarizedExperiment(assays=list(counts=ref_new@assays$RNA@counts),colData = ref_new.metadata[,c("barcodes","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc)<-pred.hesc@listData$labels
pbmc@meta.data$Cluster<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$Cluster)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:5)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:5)
DimPlot(pbmc,reduction = "ldaumap")
DimPlot(pbmc,reduction = "ldatsne")
remove(ds_seurat,pred.hesc,ref_new.metadata)
gc()
saveRDS(pbmc,"e:/CJJGSE248377_B.rds")
VlnPlot(pbmc,features = c("Irf4","Gls","Kdm6b","Ly75"),sort = TRUE)
DotPlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"))
remove(ref,ref_new,ref_new_new)
gc()
saveRDS(pbmc,file = "e:/CJJGSE248377_1.rds")

Glutamine_1<-read.table("d:/GSE60927.txt",sep = "\t")
Glutamine_1_new<-Glutamine_1[1:length(rownames(Glutamine_1)),]
Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmc))
Glutamine1<-FetchData(pbmc,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
remove(i)
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine1<-mean(Glutamine1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
remove(j)
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_Glutamine1)
Sortdatalo<-filter(Sortdata,S < average_Glutamine1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((Glutamine_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$plasmacell<-Glutamine1$Average-control$Average
remove(p,control,Data,Data_T,Glutamine_1,Glutamine1,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total,average_Glutamine1,Glutamine_1_new,higene,logene,k)

FeaturePlot(pbmc,features =c("plasmacell"),order = TRUE) + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))

pbmc2<-readRDS("e:/CJJGSE109732.rds")
pbmc2.data<-pbmc2@assays$RNA@counts
pbmc2.metadata<-pbmc2@meta.data[,c(5,22,23,25)]
remove(pbmc2)
gc()
saveRDS(pbmc,"e:/CJJGSE248377.rds")
remove(pbmc)
gc()

Idents(ref_new)<-ref_new$Clustername
ref_new_new<-subset(ref_new,idents = c("DZd","DZp" ,"LZ3", "LZ2", "LZ1"))
pbmc2.metadata$mouse_ID<-rownames(pbmc2.metadata)
pbmc2.metadata$Clustername<-pbmc2.metadata$newcluster

ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc2.data),colData = pbmc2.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = ref_new_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
ref_new@meta.data$newcluster<-rep("PC",length(rownames(ref_new@meta.data)))
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(ref_new@meta.data))) {
    if(rownames(ref_new@meta.data)[j] ==  rownames(pred.hesc)[i]){
      ref_new@meta.data$newcluster[j]<-pred.hesc@listData$labels[i]
    }
  }
}
remove(i,j,ref_new_new,pred.hesc,ds_seurat)
remove(ref,pbmc2.metadata,pbmc2.data)
gc()
ref_new <- FindVariableFeatures(ref_new, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genes <- rownames(ref_new)
ref_new <- ScaleData(ref_new, features = all.genes, verbose=F)

ref_new<-RunLDA(ref_new,labels = ref_new$newcluster)
ref_new<-RunUMAP(ref_new,reduction = "lda",reduction.name = "ldaumap",dims = 1:9)
ref_new<-RunTSNE(ref_new,reduction = "lda",reduction.name = "ldatsne",dims = 1:9)
Idents(ref_new)<-ref_new$newcluster
DimPlot(ref_new,reduction = "ldaumap")
DimPlot(ref_new,reduction = "ldatsne")
remove(all.genes)
new.cluster.ids <- c("LZ","DZ","LZ","prePB","DZ","DZ","PC","LZ","LZ","LZ")
names(new.cluster.ids) <- levels(ref_new)
ref_new <- RenameIdents(ref_new, new.cluster.ids)
ref_new@meta.data$newnewcluster<-Idents(ref_new)
remove(new.cluster.ids)
saveRDS(ref_new,"e:/CJJGSE248377_B_new.rds")
pbmc_new_new<-ref_new
library(monocle)
trace('project2MST',edit = T,where = asNamespace("monocle"))
data<-as.sparse(pbmc_new_new@assays$RNA@counts)
pd <-pbmc_new_new@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,
                              expressionFamily = VGAM::negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 
# 
# 
pbmc.marker_1<-FindAllMarkers(pbmc_new_new,only.pos = TRUE,min.pct = 0)
diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~newnewcluster",cores = 20)
remove(data,fd,fData,pbmc_new_new,pd,ref_new)
gc()
pbmcmarkers_new<-pbmc.marker_1[which(pbmc.marker_1$p_val_adj < 0.05 ),]
pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
# ordering_genes <- diff_test_res$gene_short_name
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-100))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2)
diff_test_res<-cbind(rownames(diff_test_res),diff_test_res)
colnames(diff_test_res)[1]<-"Gene_Symbol"
write.table(diff_test_res,"diff_Test_res_CJJGSE248377_B.txt",sep = "\t",row.names = FALSE,quote = FALSE)

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
saveRDS(monocle_cds,"e:/CJJGSE243877_B_monocle.rds")
