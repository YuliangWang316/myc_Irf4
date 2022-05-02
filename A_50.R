library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
#seurat 3.2.3
options(warn=-1)
set.seed(1)
pbmc.data <- Read10X(data.dir="/Users/yuliangwang/Desktop/BOAO_GCB1_2_filtered_Rawdata/")
pbmc.data <- as.data.frame(pbmc.data)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,Ighd <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB1_filtered", min.cells = 3, min.features = 200)
rm(pbmc.data,pbmc.data_t)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 50)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F)
ElbowPlot(pbmc,ndims = 50)

mtor_1<-read.table("/Users/yuliangwang/Desktop/PENG_LEUCINE_DEPRIVATION_DN_new.txt",sep = "\t",header = TRUE)
mtor_1_new<-mtor_1[2:length(rownames(mtor_1)),]
mtor_1_new<-tolower(mtor_1_new)
library(Hmisc)
mtor_1_new<-capitalize(mtor_1_new)
mtor1<-FetchData(pbmc,vars = mtor_1_new)
for (i in 1:length(rownames(mtor1))) {
  mtor1$Average[i]<-mean(as.numeric(mtor1[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_mtor1<-mean(mtor1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_mtor1)
Sortdatalo<-filter(Sortdata,S < average_mtor1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((mtor_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((mtor_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$Mtor_1<-mtor1$Average-control$Average

mtor_2<-read.table("/Users/yuliangwang/Desktop/PENG_RAPAMYCIN_RESPONSE_DN_new.txt",sep = "\t",header = TRUE)
mtor_2_new<-mtor_2[2:length(rownames(mtor_2)),]
mtor_2_new<-tolower(mtor_2_new)
library(Hmisc)
mtor_2_new<-capitalize(mtor_2_new)
mtor2<-FetchData(pbmc,vars = mtor_2_new)
for (i in 1:length(rownames(mtor2))) {
  mtor2$Average[i]<-mean(as.numeric(mtor2[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_mtor2<-mean(mtor2$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_mtor2)
Sortdatalo<-filter(Sortdata,S < average_mtor2)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((mtor_2_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((mtor_2_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$Mtor_2<-mtor2$Average-control$Average

mtor_3<-read.table("/Users/yuliangwang/Desktop/MTORC1_hallmarks.txt",sep = "\t",header = TRUE)
mtor_3_new<-mtor_3[2:length(rownames(mtor_3)),]
mtor_3_new<-tolower(mtor_3_new)
library(Hmisc)
mtor_3_new<-capitalize(mtor_3_new)
mtor3<-FetchData(pbmc,vars = mtor_3_new)
for (i in 1:length(rownames(mtor3))) {
  mtor3$Average[i]<-mean(as.numeric(mtor3[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_mtor3<-mean(mtor3$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_mtor3)
Sortdatalo<-filter(Sortdata,S < average_mtor3)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((mtor_3_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((mtor_3_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$Mtor_3<-mtor3$Average-control$Average

Glutamine_1<-read.table("/Users/yuliangwang/Desktop/GO_BP_Glutamine.txt",sep = "\t",header = TRUE)
Glutamine_1_new<-Glutamine_1[2:length(rownames(Glutamine_1)),]
Glutamine_1_new<-tolower(Glutamine_1_new)
library(Hmisc)
Glutamine_1_new<-capitalize(Glutamine_1_new)
Glutamine1<-FetchData(pbmc,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine1<-mean(Glutamine1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
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

pbmc@meta.data$Glutamine_1<-Glutamine1$Average-control$Average

Glutamine_2<-read.table("/Users/yuliangwang/Desktop/Reactome_glutamine_metabolism_new.txt",sep = "\t",header = TRUE)
Glutamine_2_new<-Glutamine_2[2:length(rownames(Glutamine_2)),]
Glutamine_2_new<-tolower(Glutamine_2_new)
library(Hmisc)
Glutamine_2_new<-capitalize(Glutamine_2_new)
Glutamine2<-FetchData(pbmc,vars = Glutamine_2_new)
for (i in 1:length(rownames(Glutamine2))) {
  Glutamine2$Average[i]<-mean(as.numeric(Glutamine2[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine2<-mean(Glutamine2$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_Glutamine2)
Sortdatalo<-filter(Sortdata,S < average_Glutamine2)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_2_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((Glutamine_2_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$Glutamine_2<-Glutamine2$Average-control$Average

remove(control,Data,Data_T,Glutamine_1,Glutamine1,Glutamine_1_new,mtor_1,mtor_1_new,mtor_2,mtor_2_new,mtor_3,mtor_3_new,mtor1,mtor2,mtor3)
remove(Total,Sortdata,Sortdatahi,Sortdatalo,rownamesDataT,all.genes,average_Glutamine1,average_mtor1,average_mtor2,average_mtor3,higene,logene,i,j,k,p)
remove(pbmcs)
pbmc <- FindNeighbors(pbmc, dims = 1:8, verbose=F)
pbmc <- FindClusters(pbmc, resolution = 1, verbose=F)
VlnPlot(pbmc,features = "Myc",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Irf4",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Kdm6b",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Ly75",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Mtor_3",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Glutamine_1",pt.size = 0,sort = TRUE)
pbmc<-RunUMAP(pbmc,dims = 1:8)
DimPlot(pbmc,reduction = "umap",label = TRUE)
pbmc<-RunTSNE(pbmc,dims = 1:8)
DimPlot(pbmc,reduction = "tsne",label = TRUE)
FeaturePlot(pbmc,features = "nCount_RNA")
FeaturePlot(pbmc,features = "nFeature_RNA")
FeaturePlot(pbmc,features = "percent.mt")
FeaturePlot(pbmc,features = "Myc",reduction = "umap",label = TRUE)
FeaturePlot(pbmc,features = "Kdm6b",reduction = "umap",label = TRUE)
FeaturePlot(pbmc,features = "Myc",reduction = "tsne",label = TRUE)
FeaturePlot(pbmc,features = "Kdm6b",reduction = "tsne",label = TRUE)
FeaturePlot(pbmc,features = "Ly75",reduction = "umap",label = TRUE)
FeaturePlot(pbmc,features = "Ly75",reduction = "tsne",label = TRUE)
FeaturePlot(pbmc,features = "Irf4",reduction = "umap",label = TRUE)
FeaturePlot(pbmc,features = "Irf4",reduction = "tsne",label = TRUE)
