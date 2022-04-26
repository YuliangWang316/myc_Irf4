library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#seurat 4.1.0
options(warn=-1)
set.seed(1)
pbmc1.data <- Read10X(data.dir="D:/BOAO single-cell sequencing/GCB1/outs/filtered_feature_bc_matrix/")
pbmc2.data <- Read10X(data.dir="D:/BOAO single-cell sequencing/GCB2/outs/filtered_feature_bc_matrix/")
pbmc1.data <- as.data.frame(pbmc1.data)
pbmc2.data <- as.data.frame(pbmc2.data)
for (i in 1:length(colnames(pbmc1.data))) {
  colnames(pbmc1.data)[i] <- paste(colnames(pbmc1.data)[i],"pbmc1",i,sep = "-")  
}

for (i in 1:length(colnames(pbmc2.data))) {
  colnames(pbmc2.data)[i] <- paste(colnames(pbmc2.data)[i],"pbmc2",i,sep = "-")  
}

pbmc.data <-cbind(pbmc1.data,pbmc2.data)
pbmc.data <- as.data.frame(pbmc.data)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,Ighd <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB_filtered_noIgd", min.cells = 3, min.features = 200)
rm(pbmc.data,pbmc.data_t,pbmc1.data,pbmc2.data)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
counts<-quantile(pbmc@meta.data$nCount_RNA,c(0.025,0.975))
Features<-quantile(pbmc@meta.data$nFeature_RNA,c(0.025,0.975))
percentmt<-quantile(pbmc@meta.data$percent.mt,c(0.025,0.975))
pbmc <- subset(pbmc, subset = nFeature_RNA > 754 & nFeature_RNA < 5169 & percent.mt < 41.25 & nCount_RNA < 32574 & nCount_RNA >1656 & percent.mt > 4.13)
pbmc1 <- NormalizeData(pbmc)
mtor_1<-read.table("I:/MAC/Memory_B cells/PENG_LEUCINE_DEPRIVATION_DN_new.txt",sep = "\t",header = TRUE)
mtor_1_new<-mtor_1[2:length(rownames(mtor_1)),]
mtor_1_new<-tolower(mtor_1_new)
library(Hmisc)
mtor_1_new<-capitalize(mtor_1_new)
mtor1<-FetchData(pbmc1,vars = mtor_1_new)
for (i in 1:length(rownames(mtor1))) {
  mtor1$Average[i]<-mean(as.numeric(mtor1[i,]))
}

Total<-FetchData(pbmc1,vars = rownames(pbmc1))
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
control<-FetchData(pbmc1,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc1@meta.data$Mtor_1<-mtor1$Average-control$Average

mtor_2<-read.table("I:/MAC/Memory_B cells/PENG_RAPAMYCIN_RESPONSE_DN_new.txt",sep = "\t",header = TRUE)
mtor_2_new<-mtor_2[2:length(rownames(mtor_2)),]
mtor_2_new<-tolower(mtor_2_new)
library(Hmisc)
mtor_2_new<-capitalize(mtor_2_new)
mtor2<-FetchData(pbmc1,vars = mtor_2_new)
for (i in 1:length(rownames(mtor2))) {
  mtor2$Average[i]<-mean(as.numeric(mtor2[i,]))
}

Total<-FetchData(pbmc1,vars = rownames(pbmc1))
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
control<-FetchData(pbmc1,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc1@meta.data$Mtor_2<-mtor2$Average-control$Average

mtor_3<-read.table("I:/MAC/Memory_B cells/MTORC1_hallmarks.txt",sep = "\t",header = TRUE)
mtor_3_new<-mtor_3[2:length(rownames(mtor_3)),]
mtor_3_new<-tolower(mtor_3_new)
library(Hmisc)
mtor_3_new<-capitalize(mtor_3_new)
mtor3<-FetchData(pbmc1,vars = mtor_3_new)
for (i in 1:length(rownames(mtor3))) {
  mtor3$Average[i]<-mean(as.numeric(mtor3[i,]))
}

Total<-FetchData(pbmc1,vars = rownames(pbmc1))
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
control<-FetchData(pbmc1,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc1@meta.data$Mtor_3<-mtor3$Average-control$Average

Glutamine_1<-read.table("I:/MAC/Memory_B cells/GO_BP_Glutamine.txt",sep = "\t",header = TRUE)
Glutamine_1_new<-Glutamine_1[2:length(rownames(Glutamine_1)),]
Glutamine_1_new<-tolower(Glutamine_1_new)
library(Hmisc)
Glutamine_1_new<-capitalize(Glutamine_1_new)
Glutamine1<-FetchData(pbmc1,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc1,vars = rownames(pbmc1))
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
control<-FetchData(pbmc1,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc1@meta.data$Glutamine_1<-Glutamine1$Average-control$Average

remove(control,Data,Data_T,Glutamine_1,Glutamine1,Glutamine_1_new,mtor_1,mtor_1_new,mtor_2,mtor_2_new,mtor_3,mtor_3_new,mtor1,mtor2,mtor3)
remove(Total,Sortdata,Sortdatahi,Sortdatalo,rownamesDataT,all.genes,average_Glutamine1,average_mtor1,average_mtor2,average_mtor3,higene,logene,i,j,k,p)

library(monocle)
data<-as.matrix(pbmc1@assays$RNA@counts)
pd <-pbmc1@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,expressionFamily=negbinomial())

monocle_cds <- estimateSizeFactors(monocle_cds)


monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(monocle_cds),
                                    num_cells_expressed >= 10))
Myc_id <- row.names(subset(fData(monocle_cds), gene_short_name == "Myc"))
Irf4_id <- row.names(subset(fData(monocle_cds),
                            gene_short_name == "Irf4"))
Ly75_id <- row.names(subset(fData(monocle_cds), gene_short_name == "Ly75"))
cth <- newCellTypeHierarchy()
Myc<-FetchData(pbmc,vars = "Myc")
Irf4<-FetchData(pbmc,vars = "Irf4")
Ly75<-FetchData(pbmc,vars = "Ly75")
cth <- addCellType(cth, "Myc", classify_func =
                     function(x) { x[Myc_id,] >= mean(Myc$Myc) })
cth <- addCellType(cth, "Irf4", classify_func =
                     function(x) { x[Irf4_id,] >= mean(Irf4$Irf4) })
cth <- addCellType(cth, "Ly75", classify_func =
                     function(x) { x[Ly75_id,] >= mean(Ly75$Ly75) })



monocle_cds <- classifyCells(monocle_cds, cth)
marker_diff <- markerDiffTable(monocle_cds[expressed_genes,],
                               cth,
                               residualModelFormulaStr = "~ num_genes_expressed",
                               cores = 20)

candidate_clustering_genes <-
  row.names(subset(marker_diff, qval < 0.01))
marker_spec <-
  calculateMarkerSpecificity(monocle_cds[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
monocle_cds <- setOrderingFilter(monocle_cds, semisup_clustering_genes)
plot_pc_variance_explained(monocle_cds, return_all = F)
monocle_cds <- reduceDimension(monocle_cds, max_components = 3, num_dim = 6,
                               norm_method = 'log',
                               reduction_method = 'tSNE',
                               residualModelFormulaStr = "~ num_genes_expressed",
                               verbose = T)
monocle_cds <- clusterCells(monocle_cds,
                            num_clusters = 12,
                            frequency_thresh = 0.1,
                            cell_type_hierarchy = cth,
                            method = "densityPeak",
                            clustering_genes = "Myc")
plot_cell_clusters(monocle_cds, 2, 3, color = "Cluster",
                   markers = c("Myc", "Irf4","Ly75","Kdm6b"))
plot_cell_clusters(monocle_cds,2,3,color_by = "Cluster" )
Idents(pbmc1)<-monocle_cds@phenoData@data[["Cluster"]]
VlnPlot(pbmc1,features = c("Myc","Irf4","Kdm6b","Ly75"),sort = TRUE,pt.size = 0)
Myc1<-subset(pbmc1,idents="8")
pbmc2<-subset(pbmc1,idents = c("2"))
