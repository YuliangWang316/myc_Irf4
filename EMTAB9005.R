library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#seurat 4.2.0
options(warn=-1)
set.seed(1)

pbmc1<-readRDS("c:/Users/xjmik/Downloads/E-MTAB-9005/SEURAT_OBJECTS/CompleteIntegrated_scRNA_SeuratObject.rds")
pbmc2<-readRDS("c:/Users/xjmik/Downloads/E-MTAB-9005/SEURAT_OBJECTS/HumanTonsil_BCells_scRNA_SeuratObject.rds")
pbmc3<-readRDS("c:/Users/xjmik/Downloads/E-MTAB-9005/SEURAT_OBJECTS/HumanTonsil_MemoryBCells_scRNA_SeuratObject.rds")

pbmc.data <- read.table("c:/Users/xjmik/Downloads/GSE188617_counts_info.tsv/counts_info.tsv",sep = "\t",header = TRUE,row.names = 1)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,IGHD <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB1_filtered", min.cells = 3, min.features = 200)
rm(pbmc.data,pbmc.data_t,pbmc.data1,pbmc.data2)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 30)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F)
ElbowPlot(pbmc,ndims = 50)

mtor_1<-read.table("c:/Users/xjmik/Desktop/BOAO_geneset/PENG_LEUCINE_DEPRIVATION_DN_new.txt",sep = "\t",header = TRUE)
mtor_1_new<-mtor_1[2:length(rownames(mtor_1)),]
#mtor_1_new<-tolower(mtor_1_new)
#library(Hmisc)
#mtor_1_new<-capitalize(mtor_1_new)
mtor_1_new<-intersect(mtor_1_new,rownames(pbmc))
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

mtor_2<-read.table("c:/Users/xjmik/Desktop/BOAO_geneset/PENG_RAPAMYCIN_RESPONSE_DN_new.txt",sep = "\t",header = TRUE)
mtor_2_new<-mtor_2[2:length(rownames(mtor_2)),]
#mtor_2_new<-tolower(mtor_2_new)
#library(Hmisc)
#mtor_2_new<-capitalize(mtor_2_new)
mtor_2_new<-intersect(mtor_2_new,rownames(pbmc))
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

mtor_3<-read.table("c:/Users/xjmik/Desktop/BOAO_geneset/MTORC1_hallmarks.txt",sep = "\t",header = TRUE)
mtor_3_new<-mtor_3[2:length(rownames(mtor_3)),]
#mtor_3_new<-tolower(mtor_3_new)
#library(Hmisc)
#mtor_3_new<-capitalize(mtor_3_new)
mtor_3_new<-intersect(mtor_3_new,rownames(pbmc))
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

Glutamine_1<-read.table("c:/Users/xjmik/Desktop/BOAO_geneset/GO_BP_Glutamine.txt",sep = "\t",header = TRUE)
Glutamine_1_new<-Glutamine_1[2:length(rownames(Glutamine_1)),]
#Glutamine_1_new<-tolower(Glutamine_1_new)
#library(Hmisc)
#Glutamine_1_new<-capitalize(Glutamine_1_new)
Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmc))
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
    
    MYC<-FetchData(pbmcs,vars = "MYC")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    MYC$ident<-idents
    noise<-rnorm(n=length(x=MYC[,"MYC"]))/1e+05
    MYC[,"MYC"]<-MYC[,"MYC"]+noise
    MYC$ident<-factor(MYC$ident,levels = names(x=rev(x=sort(x=tapply(X=MYC[,"MYC"],INDEX = MYC$ident,FUN = mean),decreasing = FALSE))))
    
    IRF4<-FetchData(pbmcs,vars = "IRF4")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    IRF4$ident<-idents
    noise<-rnorm(n=length(x=IRF4[,"IRF4"]))/1e+05
    IRF4[,"IRF4"]<-IRF4[,"IRF4"]+noise
    IRF4$ident<-factor(IRF4$ident,levels = names(x=rev(x=sort(x=tapply(X=IRF4[,"IRF4"],INDEX = IRF4$ident,FUN = mean),decreasing = FALSE))))
    
    KDM6B<-FetchData(pbmcs,vars = "KDM6B")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    KDM6B$ident<-idents
    noise<-rnorm(n=length(x=KDM6B[,"KDM6B"]))/1e+05
    KDM6B[,"KDM6B"]<-KDM6B[,"KDM6B"]+noise
    KDM6B$ident<-factor(KDM6B$ident,levels = names(x=rev(x=sort(x=tapply(X=KDM6B[,"KDM6B"],INDEX = KDM6B$ident,FUN = mean),decreasing = FALSE))))
    
    LY75<-FetchData(pbmcs,vars = "LY75")
    idents<-Idents(pbmcs)[colnames(pbmcs)]
    LY75$ident<-idents
    noise<-rnorm(n=length(x=LY75[,"LY75"]))/1e+05
    LY75[,"LY75"]<-LY75[,"LY75"]+noise
    LY75$ident<-factor(LY75$ident,levels = names(x=rev(x=sort(x=tapply(X=LY75[,"LY75"],INDEX = LY75$ident,FUN = mean),decreasing = FALSE))))
    
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
    
    
    
    if (levels(KDM6B$ident)[1] == levels(IRF4$ident)[1]){
      if (levels(KDM6B$ident)[1] == levels(LY75$ident)[1]){
        if(levels(KDM6B$ident)[1] == levels(MYC$ident)[1]){
          
          if (levels(MYC$ident)[1] == levels(MTOR_1$ident)[1]){
            if (levels(MYC$ident)[1] == levels(GLUTAMINE_1$ident)[1]){
              
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
