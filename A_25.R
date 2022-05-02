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
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
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

S<-data.frame(0,0)
colnames(S)<-c("i","j")

for (i in 1:50) {
  for (j in seq(0.1,3,by=0.1)) {
    
    pbmcs <- FindNeighbors(pbmc, dims = 1:i, verbose=F)
    pbmcs <- FindClusters(pbmcs, resolution = j, verbose=F)
    
    Myc<-FetchData(pbmcs,vars = "Myc")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    Myc$ident<-idents
    noise<-rnorm(n=length(x=Myc[,"Myc"]))/1e+05
    Myc[,"Myc"]<-Myc[,"Myc"]+noise
    Myc$ident<-factor(Myc$ident,levels = names(x=rev(x=sort(x=tapply(X=Myc[,"Myc"],INDEX = Myc$ident,FUN = mean),decreasing = FALSE))))
    
    Irf4<-FetchData(pbmcs,vars = "Irf4")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    Irf4$ident<-idents
    noise<-rnorm(n=length(x=Irf4[,"Irf4"]))/1e+05
    Irf4[,"Irf4"]<-Irf4[,"Irf4"]+noise
    Irf4$ident<-factor(Irf4$ident,levels = names(x=rev(x=sort(x=tapply(X=Irf4[,"Irf4"],INDEX = Irf4$ident,FUN = mean),decreasing = FALSE))))
    
    Kdm6b<-FetchData(pbmcs,vars = "Kdm6b")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    Kdm6b$ident<-idents
    noise<-rnorm(n=length(x=Kdm6b[,"Kdm6b"]))/1e+05
    Kdm6b[,"Kdm6b"]<-Kdm6b[,"Kdm6b"]+noise
    Kdm6b$ident<-factor(Kdm6b$ident,levels = names(x=rev(x=sort(x=tapply(X=Kdm6b[,"Kdm6b"],INDEX = Kdm6b$ident,FUN = mean),decreasing = FALSE))))
    
    Ly75<-FetchData(pbmcs,vars = "Ly75")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    Ly75$ident<-idents
    noise<-rnorm(n=length(x=Ly75[,"Ly75"]))/1e+05
    Ly75[,"Ly75"]<-Ly75[,"Ly75"]+noise
    Ly75$ident<-factor(Ly75$ident,levels = names(x=rev(x=sort(x=tapply(X=Ly75[,"Ly75"],INDEX = Ly75$ident,FUN = mean),decreasing = FALSE))))
    
    MTOR_1<-FetchData(pbmcs,vars = "Mtor_1")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    MTOR_1$ident<-idents
    noise<-rnorm(n=length(x=MTOR_1[,"Mtor_1"]))/1e+05
    MTOR_1[,"Mtor_1"]<-MTOR_1[,"Mtor_1"]+noise
    MTOR_1$ident<-factor(MTOR_1$ident,levels = names(x=rev(x=sort(x=tapply(X=MTOR_1[,"Mtor_1"],INDEX = MTOR_1$ident,FUN = mean),decreasing = FALSE))))
    
    GLUTAMINE_1<-FetchData(pbmcs,vars = "Glutamine_1")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    GLUTAMINE_1$ident<-idents
    noise<-rnorm(n=length(x=GLUTAMINE_1[,"Glutamine_1"]))/1e+05
    GLUTAMINE_1[,"Glutamine_1"]<-GLUTAMINE_1[,"Glutamine_1"]+noise
    GLUTAMINE_1$ident<-factor(GLUTAMINE_1$ident,levels = names(x=rev(x=sort(x=tapply(X=GLUTAMINE_1[,"Glutamine_1"],INDEX = GLUTAMINE_1$ident,FUN = mean),decreasing = FALSE))))
    
    MTOR_2<-FetchData(pbmcs,vars = "Mtor_2")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    MTOR_2$ident<-idents
    noise<-rnorm(n=length(x=MTOR_2[,"Mtor_2"]))/1e+05
    MTOR_2[,"Mtor_2"]<-MTOR_2[,"Mtor_2"]+noise
    MTOR_2$ident<-factor(MTOR_2$ident,levels = names(x=rev(x=sort(x=tapply(X=MTOR_2[,"Mtor_2"],INDEX = MTOR_2$ident,FUN = mean),decreasing = FALSE))))
    
    MTOR_3<-FetchData(pbmcs,vars = "Mtor_3")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    MTOR_3$ident<-idents
    noise<-rnorm(n=length(x=MTOR_3[,"Mtor_3"]))/1e+05
    MTOR_3[,"Mtor_3"]<-MTOR_3[,"Mtor_3"]+noise
    MTOR_3$ident<-factor(MTOR_3$ident,levels = names(x=rev(x=sort(x=tapply(X=MTOR_3[,"Mtor_3"],INDEX = MTOR_3$ident,FUN = mean),decreasing = FALSE))))
    
    
    
    if (levels(Kdm6b$ident)[1] == levels(Irf4$ident)[1]){
      if (levels(Kdm6b$ident)[1] == levels(Ly75$ident)[1]){
        if(levels(Kdm6b$ident)[1] == levels(Myc$ident)[1]){
          
          if (levels(Myc$ident)[1] == levels(MTOR_1$ident)[1]){
            if (levels(Myc$ident)[1] == levels(GLUTAMINE_1$ident)[1]){
              
              P<-data.frame(i,j)
              colnames(P)<-c("i","j")
              S<-rbind(S,P)
            }
          }
        }
      }
    }
  }
}
remove(control,Data,Data_T,Glutamine_1,Glutamine1,Glutamine_1_new,mtor_1,mtor_1_new,mtor_2,mtor_2_new,mtor_3,mtor_3_new,mtor1,mtor2,mtor3)
remove(Total,Sortdata,Sortdatahi,Sortdatalo,rownamesDataT,all.genes,average_Glutamine1,average_mtor1,average_mtor2,average_mtor3,higene,logene,i,j,k,p)
remove(pbmcs)
pbmc <- FindNeighbors(pbmc, dims = 1:8, verbose=F)
pbmc <- FindClusters(pbmc, resolution = 0.3, verbose=F)
VlnPlot(pbmc,features = "Myc",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Irf4",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Kdm6b",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Ly75",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Mtor_2",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Glutamine_1",pt.size = 0,sort = TRUE)
pbmc<-RunUMAP(pbmc,dims = 1:6)
DimPlot(pbmc,reduction = "umap",label = TRUE)
pbmc<-RunTSNE(pbmc,dims = 1:6)
DimPlot(pbmc,reduction = "tsne",label = TRUE)
FeaturePlot(pbmc,features = "nCount_RNA")
FeaturePlot(pbmc,features = "nFeature_RNA")
FeaturePlot(pbmc,features = "percent.mt")
FeaturePlot(pbmc,features = "Myc",reduction = "umap",label = TRUE)
FeaturePlot(pbmc,features = "Kdm6b",reduction = "umap",label = TRUE)
FeaturePlot(pbmc,features = "Myc",reduction = "tsne",label = TRUE)
FeaturePlot(pbmc,features = "Kdm6b",reduction = "tsne",label = TRUE)

DZ<-read.table("/Users/yuliangwang/Desktop/DZ_NI_new.txt",sep = "\t",header = TRUE)
LZ<-read.table("/Users/yuliangwang/Desktop/LZ_NI_new.txt",sep = "\t",header = TRUE)


DZ_new<-DZ[2:length(rownames(DZ)),]
remove(DZ)
DZ<-FetchData(pbmc,vars = DZ_new)
for (i in 1:length(rownames(DZ))) {
  DZ$Average[i]<-mean(as.numeric(DZ[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_DZ<-mean(DZ$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_DZ)
Sortdatalo<-filter(Sortdata,S < average_DZ)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((DZ_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((DZ_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$DZ<-DZ$Average-control$Average

LZ_new<-LZ[2:length(rownames(LZ)),]
remove(LZ)
LZ<-FetchData(pbmc,vars = LZ_new)
for (i in 1:length(rownames(LZ))) {
  LZ$Average[i]<-mean(as.numeric(LZ[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_LZ<-mean(LZ$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_LZ)
Sortdatalo<-filter(Sortdata,S < average_LZ)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((LZ_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((LZ_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$LZ<-LZ$Average-control$Average
remove(control,Data,Data_T,DZ,LZ,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total,average_DZ,average_LZ,DZ_new,LZ_new,higene,logene,i,j,k,p)

DZLZ<-subset(pbmc,idents = c("0","1","2","3","4","5","6","7","8","10","12","13"))

VlnPlot(DZLZ,features = "LZ",pt.size = 0,sort = TRUE)
VlnPlot(DZLZ,features = "DZ",pt.size = 0,sort = TRUE)
VlnPlot(DZLZ,features = c("Cxcr4","Cd86"),pt.size = 0,sort = TRUE)
FeaturePlot(DZLZ,features = "DZ",label = TRUE)
FeaturePlot(DZLZ,features = "LZ",label = TRUE)
FeaturePlot(DZLZ,features = "Cxcr4",label = TRUE)
FeaturePlot(DZLZ,features = "Cd86",label = TRUE)
FeaturePlot(DZLZ,features = "Cd83",label = TRUE)
FeaturePlot(DZLZ,features = "Mki67",label = TRUE)
FeaturePlot(DZLZ,features = "Aicda",label = TRUE)
FeaturePlot(DZLZ,features = "Hmmr",label = TRUE)
VlnPlot(DZLZ,features = c("Hmmr","Cd83"),pt.size = 0,sort = TRUE)
DZ<-subset(DZLZ,idents = c("1","2","5","6","8","10","12","13"))
LZ<-subset(DZLZ,idents = c("0","3","7","4"))
new.cluster.ids <- c("LZ", "DZ", "DZ", "LZ", "LZ", "DZ",
                     "DZ", "LZ", "DZ","Myc_2","DZ","Myc_1","DZ","DZ")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
Myc_1<-subset(pbmc,idents = "Myc_1")
Myc_2<-subset(pbmc,idents = "Myc_2")
remove(DZ,DZLZ,LZ,Myc_1,Myc_2)
remove(new.cluster.ids)
VlnPlot(pbmc,features  = "Mtor_1",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features  = "Glutamine_1",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features  = "Mtor_2",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features  = "Mtor_3",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features  = "Glutamine_2",pt.size = 0,sort = TRUE)
