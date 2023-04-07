library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#seurat 4.2.0
options(warn=-1)
set.seed(1)

pbmc.data <- read.table("c:/Users/xjmik/Downloads/GSE188617_counts_info.tsv/counts_info.tsv",sep = "\t",header = TRUE,row.names = 1)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,IGHD <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB1_filtered", min.cells = 3, min.features = 200)
rm(pbmc.data,pbmc.data_t,pbmc.data1,pbmc.data2)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


pbmcn <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmcn <- NormalizeData(pbmcn)
pbmcn <- FindVariableFeatures(pbmcn, selection.method = "vst", nfeatures = 2000, verbose=F)

all.genes <- rownames(pbmcn)
pbmcn <- ScaleData(pbmcn, features = all.genes, verbose=F)
pbmcn <- RunPCA(pbmcn, features = VariableFeatures(object = pbmcn), verbose=F)

mtor_1<-read.table("c:/Users/xjmik/Desktop/BOAO_geneset/PENG_LEUCINE_DEPRIVATION_DN_new.txt",sep = "\t",header = TRUE)
mtor_1_new<-mtor_1[2:length(rownames(mtor_1)),]
mtor_1_new<-intersect(mtor_1_new,rownames(pbmcn))
mtor1<-FetchData(pbmcn,vars = mtor_1_new)
for (i in 1:length(rownames(mtor1))) {
  mtor1$Average[i]<-mean(as.numeric(mtor1[i,]))
}

Total<-FetchData(pbmcn,vars = rownames(pbmcn))
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
control<-FetchData(pbmcn,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmcn@meta.data$Mtor_1<-mtor1$Average-control$Average

mtor_2<-read.table("c:/Users/xjmik/Desktop/BOAO_geneset/PENG_RAPAMYCIN_RESPONSE_DN_new.txt",sep = "\t",header = TRUE)
mtor_2_new<-mtor_2[2:length(rownames(mtor_2)),]
mtor_2_new<-intersect(mtor_2_new,rownames(pbmcn))
mtor2<-FetchData(pbmcn,vars = mtor_2_new)
for (i in 1:length(rownames(mtor2))) {
  mtor2$Average[i]<-mean(as.numeric(mtor2[i,]))
}

Total<-FetchData(pbmcn,vars = rownames(pbmcn))
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
control<-FetchData(pbmcn,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmcn@meta.data$Mtor_2<-mtor2$Average-control$Average

Glutamine_1<-read.table("c:/Users/xjmik/Desktop/BOAO_geneset/GO_BP_Glutamine.txt",sep = "\t",header = TRUE)
Glutamine_1_new<-Glutamine_1[2:length(rownames(Glutamine_1)),]
Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmcn))
Glutamine1<-FetchData(pbmcn,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmcn,vars = rownames(pbmcn))
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
control<-FetchData(pbmcn,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmcn@meta.data$Glutamine_1<-Glutamine1$Average-control$Average



mtor_3<-read.table("c:/Users/xjmik/Desktop/BOAO_geneset/MTORC1_hallmarks.txt",sep = "\t",header = TRUE)
mtor_3_new<-mtor_3[2:length(rownames(mtor_3)),]
mtor_3_new<-intersect(mtor_3_new,rownames(pbmcn))
mtor3<-FetchData(pbmcn,vars = mtor_3_new)
for (i in 1:length(rownames(mtor3))) {
  mtor3$Average[i]<-mean(as.numeric(mtor3[i,]))
}

Total<-FetchData(pbmcn,vars = rownames(pbmc))
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
control<-FetchData(pbmcn,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmcn@meta.data$Mtor_3<-mtor3$Average-control$Average
remove(control,Data,Data_T,mtor_1,mtor_2,mtor_3,mtor1,mtor2,mtor3,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total)
remove(average_mtor1,average_mtor2,average_mtor3,higene,i,j,k,logene,mtor_1_new,mtor_2_new,mtor_3_new,p)
remove(Glutamine_1,Glutamine_1_new,Glutamine_1_new,all.genes,average_Glutamine1)
remove(Glutamine1)
pbmcK <- FindNeighbors(pbmcn, dims = 1:3, verbose=F)
pbmcK <- FindClusters(pbmcK, resolution = 0.6, verbose=F)
DotPlot(pbmcK, features = c("MYC","IRF4","KDM6B","LY75","Mtor_1","Mtor_2","Mtor_3","Glutamine_1")) + RotatedAxis()
VlnPlot(pbmcK,features = c("MYC","IRF4","KDM6B","LY75","Mtor_1","Mtor_2","Mtor_3","Glutamine_1"),sort = TRUE,pt.size = 0)
