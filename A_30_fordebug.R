library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
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
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 30)
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

plasma<-read.table("/Users/yuliangwang/Desktop/plasma_vs_naive_IRF4_new.tsv",sep = "\t",header = TRUE)
plasma_new<-plasma[2:length(rownames(plasma)),]
plasma_new<-tolower(plasma_new)
library(Hmisc)
plasma_new<-capitalize(plasma_new)
remove(plasma)
plasma<-FetchData(pbmc,vars = plasma_new)
for (i in 1:length(rownames(plasma))) {
  plasma$Average[i]<-mean(as.numeric(plasma[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_plasma<-mean(plasma$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_plasma)
Sortdatalo<-filter(Sortdata,S < average_plasma)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((plasma_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((plasma_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$plasma<-plasma$Average-control$Average

KEGGcellcycle<-read.table("/Users/yuliangwang/Desktop/Cellcycle.txt",sep = "\t",header = TRUE)
KEGGcellcycle_new<-KEGGcellcycle[2:length(rownames(KEGGcellcycle)),]
KEGGcellcycle_new<-tolower(KEGGcellcycle_new)
library(Hmisc)
KEGGcellcycle_new<-capitalize(KEGGcellcycle_new)

cellcycle<-FetchData(pbmc,vars = KEGGcellcycle_new)
for (i in 1:length(rownames(cellcycle))) {
  cellcycle$Average[i]<-mean(as.numeric(cellcycle[i,]))
}

Total<-FetchData(pbmc,vars = rownames(pbmc))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_cellcycle<-mean(cellcycle$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_cellcycle)
Sortdatalo<-filter(Sortdata,S < average_cellcycle)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((KEGGcellcycle_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((KEGGcellcycle_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc@meta.data$KEGGcellcycle<-cellcycle$Average-control$Average
remove(control,Data,Data_T,Glutamine_1,Glutamine1,Glutamine_1_new,mtor_1,mtor_1_new,mtor_2,mtor_2_new,mtor_3,mtor_3_new,mtor1,mtor2,mtor3)
remove(Total,Sortdata,Sortdatahi,Sortdatalo,rownamesDataT,all.genes,average_Glutamine1,average_mtor1,average_mtor2,average_mtor3,higene,logene,i,j,k,p)
remove(pbmcs,cellcycle,KEGGcellcycle,plasma,average_cellcycle,average_plasma,KEGGcellcycle_new,plasma_new)

pbmc <- FindNeighbors(pbmc, dims = 1:8, verbose=F)
pbmc <- FindClusters(pbmc, resolution = 0.6, verbose=F)
#VlnPlot(pbmc,features = "plasma",pt.size = 0,sort = TRUE)
#VlnPlot(pbmc,features = "KEGGcellcycle",pt.size = 0,sort = TRUE)
VlnPlot(pbmc,features = "Myc",pt.size = 0)
library(scales)
show_col(hue_pal()(11)) 
VlnPlot(pbmc,features = "Irf4",pt.size = 0)
VlnPlot(pbmc,features = "Kdm6b",pt.size = 0)
VlnPlot(pbmc,features = "Ly75",pt.size = 0)
VlnPlot(pbmc,features = "Mtor_1",pt.size = 0)+geom_boxplot(width=.1,col="black",fill="white")
VlnPlot(pbmc,features = "Mtor_2",pt.size = 0)+geom_boxplot(width=.1,col="black",fill="white")
VlnPlot(pbmc,features = "Mtor_3",pt.size = 0)+geom_boxplot(width=.1,col="black",fill="white")
VlnPlot(pbmc,features = "Glutamine_1",pt.size = 0)+geom_boxplot(width=.1,col="black",fill="white")
pbmc<-RunUMAP(pbmc,dims = 1:8)
DimPlot(pbmc,reduction = "umap",label = TRUE)
DimPlot(pbmc,reduction = "umap")
#pbmc<-RunTSNE(pbmc,dims = 1:8)
#DimPlot(pbmc,reduction = "tsne",label = TRUE)
MTOR_1<-FetchData(pbmc,vars = "Mtor_1")
idents<-Idents(pbmc)[colnames(pbmc)]
MTOR_1$ident<-idents
noise<-rnorm(n=length(x=MTOR_1[,"Mtor_1"]))/1e+05
MTOR_1[,"Mtor_1"]<-MTOR_1[,"Mtor_1"]+noise
MTOR_1$ident<-factor(MTOR_1$ident,levels = names(x=rev(x=sort(x=tapply(X=MTOR_1[,"Mtor_1"],INDEX = MTOR_1$ident,FUN = mean),decreasing = FALSE))))
write.table(MTOR_1,file = "/Users/yuliangwang/Desktop/untitled folder/MTOR1.txt",sep = "\t")
remove(MTOR_1,idents,noise)
GLUTAMINE_1<-FetchData(pbmc,vars = "Glutamine_1")
idents<-Idents(pbmc)[colnames(pbmc)]
GLUTAMINE_1$ident<-idents
noise<-rnorm(n=length(x=GLUTAMINE_1[,"Glutamine_1"]))/1e+05
GLUTAMINE_1[,"Glutamine_1"]<-GLUTAMINE_1[,"Glutamine_1"]+noise
GLUTAMINE_1$ident<-factor(GLUTAMINE_1$ident,levels = names(x=rev(x=sort(x=tapply(X=GLUTAMINE_1[,"Glutamine_1"],INDEX = GLUTAMINE_1$ident,FUN = mean),decreasing = FALSE))))
write.table(GLUTAMINE_1,file = "/Users/yuliangwang/Desktop/untitled folder/GLUTAMINE1.txt",sep = "\t")
remove(GLUTAMINE_1,idents,noise)
MTOR_2<-FetchData(pbmc,vars = "Mtor_2")
idents<-Idents(pbmc)[colnames(pbmc)]
MTOR_2$ident<-idents
noise<-rnorm(n=length(x=MTOR_2[,"Mtor_2"]))/1e+05
MTOR_2[,"Mtor_2"]<-MTOR_2[,"Mtor_2"]+noise
MTOR_2$ident<-factor(MTOR_2$ident,levels = names(x=rev(x=sort(x=tapply(X=MTOR_2[,"Mtor_2"],INDEX = MTOR_2$ident,FUN = mean),decreasing = FALSE))))
write.table(MTOR_2,file = "/Users/yuliangwang/Desktop/untitled folder/MTOR2.txt",sep = "\t")
remove(MTOR_2,idents,noise)
#FeaturePlot(pbmc,features = "nCount_RNA")
#FeaturePlot(pbmc,features = "nFeature_RNA")
#FeaturePlot(pbmc,features = "percent.mt")
FeaturePlot(pbmc,features = "Myc",reduction = "umap")
#FeaturePlot(pbmc,features = "Myc",reduction = "umap") + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "Kdm6b",reduction = "umap")

#FeaturePlot(pbmc,features = "Myc",reduction = "tsne",label = TRUE)
#FeaturePlot(pbmc,features = "Kdm6b",reduction = "tsne",label = TRUE)
FeaturePlot(pbmc,features = "Ly75",reduction = "umap")
#FeaturePlot(pbmc,features = "Ly75",reduction = "tsne",label = TRUE)
FeaturePlot(pbmc,features = "Irf4",reduction = "umap")
#FeaturePlot(pbmc,features = "Irf4",reduction = "tsne",label = TRUE)
FeaturePlot(pbmc,features = "Mtor_1",reduction = "umap",label = FALSE)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "Mtor_1",reduction = "umap",label = FALSE,min.cutoff = 0)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "Mtor_2",reduction = "umap",label = FALSE)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "Mtor_2",reduction = "umap",label = FALSE,min.cutoff = 0)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "Mtor_3",reduction = "umap",label = FALSE)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "Mtor_3",reduction = "umap",label = FALSE,min.cutoff = 0)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "Glutamine_1",reduction = "umap",label = FALSE)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "Glutamine_1",reduction = "umap",label = FALSE,min.cutoff = 0)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))

FeaturePlot(pbmc,features = "plasma",label = TRUE)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "KEGGcellcycle",label = TRUE)

Cluster10<-subset(pbmc,idents = "10")
Cluster8<-subset(pbmc,idents = "8")
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

DZLZ<-subset(pbmc,idents = c("0","1","2","3","4","5","6","7","9"))

VlnPlot(DZLZ,features = "LZ",pt.size = 0,sort = TRUE)
VlnPlot(DZLZ,features = "DZ",pt.size = 0,sort = TRUE)
VlnPlot(DZLZ,features = c("Cxcr4"),pt.size = 0,sort = TRUE)
VlnPlot(DZLZ,features = c("Cd86"),pt.size = 0,sort = TRUE)
VlnPlot(DZLZ,features = c("Cd83"),pt.size = 0,sort = TRUE)
VlnPlot(DZLZ,features = c("Mki67"),pt.size = 0,sort = TRUE)

FeaturePlot(DZLZ,features = "DZ")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(DZLZ,features = "LZ")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(DZLZ,features = "Cxcr4")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(DZLZ,features = "Cd86")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(DZLZ,features = "Cd83")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(DZLZ,features = "Mki67")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(DZLZ,features = "Aicda",label = TRUE)
FeaturePlot(DZLZ,features = "Hmmr")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
s.genes<-tolower(s.genes)
library(Hmisc)
s.genes<-capitalize(s.genes)
g2m.genes<-tolower(g2m.genes)
library(Hmisc)
g2m.genes<-capitalize(g2m.genes)
s.genes_new<-s.genes[-which(s.genes=="Mlf1ip")]
g2m.genes_new<-g2m.genes[-which(g2m.genes=="Fam64a")]
pbmc<-CellCycleScoring(pbmc,s.features = s.genes_new,g2m.features = g2m.genes_new)
pbmc<-AddModuleScore(pbmc,features = s.genes_new,name = "sgenes")

DZLZ<-CellCycleScoring(DZLZ,s.features = s.genes,g2m.features = g2m.genes)
FeaturePlot(DZLZ,features = c("S.Score","G2M.Score"))
DZLZ2<-CellCycleScoring(DZLZ,s.features = s.genes,g2m.features = g2m.genes,set.ident = TRUE)
DZLZ2<-RunPCA(DZLZ2,features = c(s.genes,g2m.genes))
DimPlot(DZLZ2,reduction = "umap")
DZ<-subset(DZLZ,idents = c("1","2","4","5","6","9"))
LZ<-subset(DZLZ,idents = c("0","3","7"))
new.cluster.ids <- c("LZ1", "DZ1", "DZ2", "LZ2", "DZ3", "DZ4",
                     "DZ5", "LZ3", "Myc_2","DZ6","Myc_1")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
Myc_1<-subset(pbmc,idents = "Myc_1")
Myc_2<-subset(pbmc,idents = "Myc_2")
remove(DZ,DZLZ,LZ)
remove(new.cluster.ids,g2m.genes,s.genes)
remove(DZLZ2)
VlnPlot(pbmc,features = "Myc",pt.size = 0)
pbmc.copy<-pbmc
my_levels <- c("DZ1","DZ2","DZ3","DZ4","DZ5","DZ6","LZ1","LZ2","LZ3","Myc_1","Myc_2")
Idents(pbmc.copy) <- factor(Idents(pbmc.copy), levels= my_levels)

VlnPlot(pbmc.copy,features = "Irf4",pt.size = 0)+NoLegend()+scale_fill_manual(values = c("red","red","red","red","red","red","red","red","red","#FF63B6","red"))
VlnPlot(pbmc.copy,features = "Kdm6b",pt.size = 0)+NoLegend()+scale_fill_manual(values = c("red","red","red","red","red","red","red","red","red","#FF63B6","red"))
VlnPlot(pbmc.copy,features = "Ly75",pt.size = 0)+NoLegend()+scale_fill_manual(values = c("red","red","red","red","red","red","red","red","red","#FF63B6","red"))
VlnPlot(pbmc.copy,features = "Mtor_1",pt.size = 0)+geom_boxplot(width=.1,col="black",fill="white")+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Mtor_2",pt.size = 0)+geom_boxplot(width=.1,col="black",fill="white")+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Mtor_3",pt.size = 0)+geom_boxplot(width=.1,col="black",fill="white")+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Kdm6a",pt.size = 0)+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Ezh2",pt.size = 0)+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Glutamine_1",pt.size = 0)+geom_boxplot(width=.1,col="black",fill="white")+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Cd40",pt.size = 0)+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Cd274",pt.size = 0)+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Cd81",pt.size = 0)+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Cd82",pt.size = 0)+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Cd83",pt.size = 0)+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
VlnPlot(pbmc.copy,features = "Cd86",pt.size = 0)+scale_fill_manual(values = c("#DB8E00","#AEA200","#00BD5C","#00C1A7","#00BADE","#EF67EB","#F8766D","#64B200","#00A6FF","#FF63B6","#B385FF"))
FeaturePlot(pbmc,features = "Ezh2")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "Kdm6a")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))

Myc1Myc2<-subset(pbmc,idents = c("Myc_1","Myc_2"))

VlnPlot(Myc1Myc2,features = "plasma",pt.size = 0,sort = TRUE)+geom_boxplot(width=.1,col="black",fill="white")+NoLegend()+scale_fill_manual(values = c("#FF63B6","#B385FF"))
Myc1Myc2.copy<-Myc1Myc2
my_levels2<-c("Myc_1","Myc_2")
Idents(Myc1Myc2.copy) <- factor(Idents(Myc1Myc2.copy), levels= my_levels2)
VlnPlot(Myc1Myc2.copy,features = "KEGGcellcycle",pt.size = 0)+geom_boxplot(width=.1,col="black",fill="white")+NoLegend()+scale_fill_manual(values = c("#FF63B6","#B385FF"))
FeaturePlot(Myc1Myc2,features = "plasma",reduction = "umap",label = TRUE)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(Myc1Myc2,features = "KEGGcellcycle",reduction = "umap",label = TRUE)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
Myc_1<-subset(pbmc,idents = "10")
Myc_2<-subset(pbmc,idents = "8")
remove(Myc1Myc2)
myc1<-as.data.frame(t(FetchData(Myc_1,vars = rownames(Myc_1))))
myc2<-as.data.frame(t(FetchData(Myc_2,vars = rownames(Myc_2))))
myc1s<-as.data.frame(Myc_1@assays[["RNA"]]@scale.data)
myc2s<-as.data.frame(Myc_2@assays[["RNA"]]@scale.data)
write.table(myc1,file = "/Users/yuliangwang/Desktop/myc1.txt",sep = "\t")
write.table(myc2,file = "/Users/yuliangwang/Desktop/myc2.txt",sep = "\t")
write.table(myc1s,file = "/Users/yuliangwang/Desktop/myc1s.txt",sep = "\t")
write.table(myc2s,file = "/Users/yuliangwang/Desktop/myc2s.txt",sep = "\t")
pbmc.markers<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
write.table(pbmc.markers,file = "/Users/yuliangwang/Desktop/untitled folder/pbmc.markers.txt",sep = "\t")
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
pbmc.avg<-AverageExpression(pbmc)
write.table(pbmc.avg$RNA,file = "/Users/yuliangwang/Desktop/untitled folder/pbmc.avg.txt",sep = "\t")
WXMlist<-read.table("/Users/yuliangwang/Desktop/WXMlist.txt",header = TRUE,sep = "\t")
WXMlist_new<-WXMlist[2:length(rownames(WXMlist)),]
WXMlist_new<-tolower(WXMlist_new)
library(Hmisc)
WXMlist_new<-capitalize(WXMlist_new)
DoHeatmap(pbmc, features = WXMlist_new) 
pbmc.copy<-pbmc
my_levels <- c("DZ1","DZ2","DZ3","DZ4","DZ5","DZ6","LZ1","LZ2","LZ3","Myc_1","Myc_2")
Idents(pbmc.copy) <- factor(Idents(pbmc.copy), levels= my_levels)
DoHeatmap(pbmc.copy, features = WXMlist_new) 
myc1_myc2<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "Myc_2",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc1_myc2,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_myc2.txt",sep = "\t")
pbmc.markers_bimod<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,test.use = "bimod")
pbmc.markers_roc<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,test.use = "roc")
pbmc.markers_t<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,test.use = "t")
pbmc.markers_poisson<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,test.use = "poisson")
pbmc.markers_negbinom<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,test.use = "negbinom")
pbmc.markers_LR<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,test.use = "LR")
pbmc.markers_MAST<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,test.use = "MAST")
pbmc.markers_DESeq2<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,test.use = "DESeq2")
myc_1_DZ1<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "DZ1",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc_1_DZ1,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_DZ1.txt",sep = "\t")

myc_1_LZ1<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "LZ1",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc_1_LZ1,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_LZ1.txt",sep = "\t")
myc_1_LZ2<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "LZ2",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc_1_LZ2,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_LZ2.txt",sep = "\t")
myc_1_LZ3<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "LZ3",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc_1_LZ3,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_LZ3.txt",sep = "\t")
myc_1_DZ2<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "DZ2",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc_1_DZ2,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_DZ2.txt",sep = "\t")
myc_1_DZ3<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "DZ3",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc_1_DZ3,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_DZ3.txt",sep = "\t")
myc_1_DZ4<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "DZ4",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc_1_DZ4,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_DZ4.txt",sep = "\t")
myc_1_DZ5<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "DZ5",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc_1_DZ5,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_DZ5.txt",sep = "\t")
myc_1_DZ6<-FindMarkers(pbmc,ident.1 = "Myc_1",ident.2 = "DZ6",only.pos = TRUE,min.pct = 0,logfc.threshold = 0)
write.table(myc_1_DZ6,file = "/Users/yuliangwang/Desktop/untitled folder/my1_vs_DZ6.txt",sep = "\t")
VlnPlot(pbmc.copy,features = "Kdm6a",pt.size = 0)
plasma<-FetchData(pbmc,vars = "plasma")
idents<-Idents(pbmc)[colnames(pbmc)]
plasma$ident<-idents
noise<-rnorm(n=length(x=plasma[,"plasma"]))/1e+05
plasma[,"plasma"]<-plasma[,"plasma"]+noise
plasma$ident<-factor(plasma$ident,levels = names(x=rev(x=sort(x=tapply(X=plasma[,"plasma"],INDEX = plasma$ident,FUN = mean),decreasing = FALSE))))
write.table(plasma,file = "/Users/yuliangwang/Desktop/untitled folder/plasma.txt",sep = "\t")
KEGGcellcyle<-FetchData(pbmc,vars = "KEGGcellcycle")
idents<-Idents(pbmc)[colnames(pbmc)]
KEGGcellcyle$ident<-idents
noise<-rnorm(n=length(x=KEGGcellcyle[,"KEGGcellcycle"]))/1e+05
KEGGcellcyle[,"KEGGcellcycle"]<-KEGGcellcyle[,"KEGGcellcycle"]+noise
KEGGcellcyle$ident<-factor(KEGGcellcyle$ident,levels = names(x=rev(x=sort(x=tapply(X=KEGGcellcyle[,"KEGGcellcycle"],INDEX = KEGGcellcyle$ident,FUN = mean),decreasing = FALSE))))
write.table(KEGGcellcyle,file = "/Users/yuliangwang/Desktop/untitled folder/KEGGcellcycle.txt",sep = "\t")
pbmc@meta.data$Clusters<-Idents(pbmc)

new.cluster.ids <- c("LZ", "DZ", "DZ", "LZ", "DZ", "DZ",
                     "DZ", "LZ", "Myc_2","DZ","Myc_1")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$clusters<-Idents(pbmc)
Idents(pbmc)<-pbmc@meta.data$clusters
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
pbmc.marker<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~clusters")
ordering_genes<-c("Myc","Irf4","Kdm6b","Ly75","Xbp1")
ordering_genes <- row.names (subset(pbmc.marker,avg_logFC >0.2 & p_val_adj < 1e-200))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds, method = 'DDRTree',max_components = 2)
monocle_cds <-orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "Clusters",cell_size = 0.75)
remove(data,diff_test_res,fd,fData,pbmc.marker,pd,new.cluster.ids)

a.genes<-c("Irf4","Kdm6b")
b.genes<-c("Myc","Ly75")
pbmcs<-AddModuleScore(pbmc,features = a.genes)
pbmcs<-CellCycleScoring(pbmc,s.features = a.genes,g2m.features = b.genes)
