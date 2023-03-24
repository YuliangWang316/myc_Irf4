library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#seurat 4.2.0
options(warn=-1)
set.seed(1)
pbmc.data1 <- read.table("c:/Users/xjmik/Downloads/GSE139891_RAW/GSM4148370_GC1_umi_grch38.txt/RK10001_counts.txt",sep = "\t",header = TRUE)
pbmc.data1 <- pbmc.data1[!duplicated(pbmc.data1$X),]
rownames(pbmc.data1)<-pbmc.data1[,1]
pbmc.data1<-pbmc.data1[,-1]
pbmc.data2 <- read.table("c:/Users/xjmik/Downloads/GSE139891_RAW/GSM4148371_GC2_umi_grch38.txt/RK10004_counts.txt",sep = "\t",header = TRUE)
pbmc.data2 <- pbmc.data2[!duplicated(pbmc.data2$X),]
rownames(pbmc.data2)<-pbmc.data2[,1]
pbmc.data2<-pbmc.data2[,-1]
pbmc.data3 <- read.table("c:/Users/xjmik/Downloads/GSE139891_RAW/GSM4148372_DZ1_umi_grch38.txt/RK10002_counts.txt",sep = "\t",header = TRUE)
pbmc.data3 <- pbmc.data3[!duplicated(pbmc.data3$X),]
rownames(pbmc.data3)<-pbmc.data3[,1]
pbmc.data3<-pbmc.data3[,-1]
pbmc.data4 <- read.table("c:/Users/xjmik/Downloads/GSE139891_RAW/GSM4148373_DZ2_umi_grch38.txt/RK10005_counts.txt",sep = "\t",header = TRUE)
pbmc.data4 <- pbmc.data4[!duplicated(pbmc.data4$X),]
rownames(pbmc.data4)<-pbmc.data4[,1]
pbmc.data4<-pbmc.data4[,-1]
pbmc.data5 <- read.table("c:/Users/xjmik/Downloads/GSE139891_RAW/GSM4148374_LZ1_umi_grch38.txt/RK10003_counts.txt",sep = "\t",header = TRUE)
pbmc.data5 <- pbmc.data5[!duplicated(pbmc.data5$X),]
rownames(pbmc.data5)<-pbmc.data5[,1]
pbmc.data5<-pbmc.data5[,-1]
pbmc.data6 <- read.table("c:/Users/xjmik/Downloads/GSE139891_RAW/GSM4148375_LZ2_umi_grch38.txt/RK10006_counts.txt",sep = "\t",header = TRUE)
pbmc.data6 <- pbmc.data6[!duplicated(pbmc.data6$X),]
rownames(pbmc.data6)<-pbmc.data6[,1]
pbmc.data6<-pbmc.data6[,-1]
pbmc.data7 <- read.table("c:/Users/xjmik/Downloads/GSE139891_RAW/GSM4560814_GC3a_umi.txt/counts_no_bad.txt",sep = "\t",header = TRUE)
pbmc.data7 <- pbmc.data7[!duplicated(pbmc.data7$X),]
rownames(pbmc.data7)<-pbmc.data7[,1]
pbmc.data7<-pbmc.data7[,-1]
pbmc.data8 <- read.table("c:/Users/xjmik/Downloads/GSE139891_RAW/GSM4560815_GC3b_umi.txt/counts_no_bad.txt",sep = "\t",header = TRUE)
pbmc.data8 <- pbmc.data8[!duplicated(pbmc.data8$X),]
rownames(pbmc.data8)<-pbmc.data8[,1]
pbmc.data8<-pbmc.data8[,-1]

for (i in 1:length(colnames(pbmc.data1))) {
  colnames(pbmc.data1)[i] <- paste(colnames(pbmc.data1)[i],"A",i,sep = "-")  
}

for (i in 1:length(colnames(pbmc.data2))) {
  colnames(pbmc.data2)[i] <- paste(colnames(pbmc.data2)[i],"B",i,sep = "-")  
}
for (i in 1:length(colnames(pbmc.data3))) {
  colnames(pbmc.data3)[i] <- paste(colnames(pbmc.data3)[i],"C",i,sep = "-")  
}
for (i in 1:length(colnames(pbmc.data4))) {
  colnames(pbmc.data4)[i] <- paste(colnames(pbmc.data4)[i],"D",i,sep = "-")  
}
for (i in 1:length(colnames(pbmc.data5))) {
  colnames(pbmc.data5)[i] <- paste(colnames(pbmc.data5)[i],"E",i,sep = "-")  
}
for (i in 1:length(colnames(pbmc.data6))) {
  colnames(pbmc.data6)[i] <- paste(colnames(pbmc.data6)[i],"F",i,sep = "-")  
}
for (i in 1:length(colnames(pbmc.data7))) {
  colnames(pbmc.data7)[i] <- paste(colnames(pbmc.data7)[i],"G",i,sep = "-")  
}
for (i in 1:length(colnames(pbmc.data8))) {
  colnames(pbmc.data8)[i] <- paste(colnames(pbmc.data8)[i],"H",i,sep = "-")  
}
a<-intersect(rownames(pbmc.data1),rownames(pbmc.data2))
b<-intersect(rownames(pbmc.data3),rownames(pbmc.data4))
c<-intersect(rownames(pbmc.data5),rownames(pbmc.data6))
d<-intersect(rownames(pbmc.data7),rownames(pbmc.data8))
a<-intersect(a,b)
remove(b)
c<-intersect(c,d)
remove(d)
a<-intersect(a,c)
remove(c,i)
pbmc.data1<-pbmc.data1[a,]
pbmc.data2<-pbmc.data2[a,]
pbmc.data3<-pbmc.data3[a,]
pbmc.data4<-pbmc.data4[a,]
pbmc.data5<-pbmc.data5[a,]
pbmc.data6<-pbmc.data6[a,]
pbmc.data7<-pbmc.data7[a,]
pbmc.data8<-pbmc.data8[a,]
remove(a)
pbmc.data<-cbind(pbmc.data2,pbmc.data1,pbmc.data3,pbmc.data4,pbmc.data5,pbmc.data6,pbmc.data7,pbmc.data8)
remove(pbmc.data2,pbmc.data1,pbmc.data3,pbmc.data4,pbmc.data5,pbmc.data6,pbmc.data7,pbmc.data8)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,IGHD <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB1_filtered", min.cells = 3, min.features = 200)
rm(pbmc.data,pbmc.data_t,pbmc.data1,pbmc.data2)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


S<-data.frame(0,0,0,0)
colnames(S)<-c("i","j","h","k")
for (h in c(3000,2500)) {
  for (k in c(15,10,5)) {
    pbmcs <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < h & percent.mt < k)
    pbmcs <- NormalizeData(pbmcs)
    pbmcs <- FindVariableFeatures(pbmcs, selection.method = "vst", nfeatures = 2000, verbose=F)
    
    all.genes <- rownames(pbmcs)
    pbmcs <- ScaleData(pbmcs, features = all.genes, verbose=F)
    pbmcs <- RunPCA(pbmcs, features = VariableFeatures(object = pbmcs), verbose=F)

    for (i in 1:50) {
      for (j in seq(0.1,6,by=0.1)) {
        
        pbmcs <- FindNeighbors(pbmcs, dims = 1:i, verbose=F)
        pbmcs <- FindClusters(pbmcs, resolution = j, verbose=F)
        
        MYC<-FetchData(pbmcs,vars = "MYC")
        idents<-Idents(pbmcs)[colnames(pbmcs)]
        MYC$ident<-idents
        noise<-rnorm(n=length(x=MYC[,"MYC"]))/1e+05
        MYC[,"MYC"]<-MYC[,"MYC"]+noise
        MYC$ident<-factor(MYC$ident,levels = names(x=rev(x=sort(x=tapply(X=MYC[,"MYC"],INDEX = MYC$ident,FUN = mean),decreasing = FALSE))))
        
        LY75<-FetchData(pbmcs,vars = "LY75")
        idents<-Idents(pbmcs)[colnames(pbmcs)]
        LY75$ident<-idents
        noise<-rnorm(n=length(x=LY75[,"LY75"]))/1e+05
        LY75[,"LY75"]<-LY75[,"LY75"]+noise
        LY75$ident<-factor(LY75$ident,levels = names(x=rev(x=sort(x=tapply(X=LY75[,"LY75"],INDEX = LY75$ident,FUN = mean),decreasing = FALSE))))
        
        KDM6B<-FetchData(pbmcs,vars = "KDM6B")
        idents<-Idents(pbmcs)[colnames(pbmcs)]
        KDM6B$ident<-idents
        noise<-rnorm(n=length(x=KDM6B[,"KDM6B"]))/1e+05
        KDM6B[,"KDM6B"]<-KDM6B[,"KDM6B"]+noise
        KDM6B$ident<-factor(KDM6B$ident,levels = names(x=rev(x=sort(x=tapply(X=KDM6B[,"KDM6B"],INDEX = KDM6B$ident,FUN = mean),decreasing = FALSE))))
        
        PRDM1<-FetchData(pbmcs,vars = "PRDM1")
        idents<-Idents(pbmcs)[colnames(pbmcs)]
        PRDM1$ident<-idents
        noise<-rnorm(n=length(x=PRDM1[,"PRDM1"]))/1e+05
        PRDM1[,"PRDM1"]<-PRDM1[,"PRDM1"]+noise
        PRDM1$ident<-factor(PRDM1$ident,levels = names(x=rev(x=sort(x=tapply(X=PRDM1[,"PRDM1"],INDEX = PRDM1$ident,FUN = mean),decreasing = FALSE))))
        
        
        if (levels(MYC$ident)[1] == levels(LY75$ident)[1]){
          if (levels(MYC$ident)[1] == levels(PRDM1$ident)[1]){
            if(levels(MYC$ident)[1] == levels(KDM6B$ident)[1]){
              
              
              
              P<-data.frame(i,j,h,k)
              colnames(P)<-c("i","j","h","k")
              S<-rbind(S,P)
              
            }
          }
        }
      }
    }
  }
  
  
}
remove(KDM6B,LY75,MYC,P,S,PRDM1,all.genes,pbmcs,h,i,idents,j,k,noise)
pbmcn <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
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
pbmcK <- FindNeighbors(pbmcn, dims = 1:3, verbose=F)
pbmcK <- FindClusters(pbmcK, resolution = 0.3, verbose=F)
DotPlot(pbmcK, features = c("MYC","PRDM1","KDM6B","LY75","Mtor_1","Mtor_2","Mtor_3","IRF4")) + RotatedAxis()
VlnPlot(pbmcK,features = c("MYC","PRDM1","KDM6B","LY75","Mtor_1","Mtor_2","Mtor_3","IRF4"),sort = TRUE,pt.size = 0)

remove(pbmcK,pbmcn,all.genes)
