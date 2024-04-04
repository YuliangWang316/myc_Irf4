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
# plan("multisession", workers = 8)
# options(future.globals.maxSize = 8000 * 1024^2)
GCB1.data<-Read10X("D:/GCB1_2.2.filtered_feature_bc_matrix/")
GCB2.data<-Read10X("D:/GCB2_2.2.filtered_feature_bc_matrix/")

GCB1.data<-as.data.frame(GCB1.data)
GCB2.data<-as.data.frame(GCB2.data)

for (i in 1:length(colnames(GCB1.data))) {
  colnames(GCB1.data)[i] <- paste(colnames(GCB1.data)[i],"GCB1",i,sep = "-")  
}

for (i in 1:length(colnames(GCB2.data))) {
  colnames(GCB2.data)[i] <- paste(colnames(GCB2.data)[i],"GCB2",i,sep = "-")  
}

GCB1.metadata<-data.frame(colnames(GCB1.data),rep("GCB1",length(colnames(GCB1.data))))
GCB2.metadata<-data.frame(colnames(GCB2.data),rep("GCB2",length(colnames(GCB2.data))))
colnames(GCB1.metadata)<-c("barcode","group")
colnames(GCB2.metadata)<-c("barcode","group")
rownames(GCB1.metadata)<-GCB1.metadata[,1]
rownames(GCB2.metadata)<-GCB2.metadata[,1]
pbmc.metadata<-rbind(GCB1.metadata,GCB2.metadata)

pbmc.data<-cbind(GCB1.data,GCB2.data)

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB",meta.data = pbmc.metadata,min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 50)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000, verbose=F)
VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F,model.use = "poisson")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F,npcs = 100)
ElbowPlot(pbmc,ndims = 100)

mtor_1<-read.table("D:/PENG_RAPAMYCIN_RESPONSE_DN.txt",sep = "\t",header = TRUE)
mtor_1_new<-mtor_1[2:length(rownames(mtor_1)),]
mtor_1_new<-tolower(mtor_1_new)
library(Hmisc)
mtor_1_new<-capitalize(mtor_1_new)
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

mtor_2<-read.table("D:/PENG_LEUCINE_DEPRIVATION_DN.txt",sep = "\t",header = TRUE)
mtor_2_new<-mtor_2[2:length(rownames(mtor_2)),]
mtor_2_new<-tolower(mtor_2_new)
library(Hmisc)
mtor_2_new<-capitalize(mtor_2_new)
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

Glutamine_1<-read.table("D:/GOBP_Glutamine.txt",sep = "\t",header = TRUE)
Glutamine_1_new<-Glutamine_1[2:length(rownames(Glutamine_1)),]
Glutamine_1_new<-tolower(Glutamine_1_new)
library(Hmisc)
Glutamine_1_new<-capitalize(Glutamine_1_new)
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

pbmc_Totalgene<-as.data.frame(t(as.data.frame(pbmc@assays$RNA@data)[c("Myc","Irf4","Ly75","Kdm6b"),]))
pbmc_Myc<-pbmc_Totalgene[which(pbmc_Totalgene$Myc > 0 ),]
pbmc_Irf4<-pbmc_Totalgene[which(pbmc_Totalgene$Irf4 > 0),]
pbmc_Ly75<-pbmc_Totalgene[which(pbmc_Totalgene$Ly75 > 0),]
pbmc_Kdm6b<-pbmc_Totalgene[which(pbmc_Totalgene$Kdm6b > 0),]

write.table(pbmc_Irf4,"pbmc_Irf4.txt",sep = "\t")
write.table(pbmc_Kdm6b,"pbmc_Kdm6b.txt",sep = "\t")
write.table(pbmc_Ly75,"pbmc_Ly75.txt",sep = "\t")
write.table(pbmc_Myc,"pbmc_Myc.txt",sep = "\t")

pbmc_Irf4<-rownames(pbmc_Irf4)
pbmc_Kdm6b<-rownames(pbmc_Kdm6b)
pbmc_Ly75<-rownames(pbmc_Ly75)
pbmc_Myc<-rownames(pbmc_Myc)

pbmc@meta.data$Myc <- rep("neg",length(rownames(pbmc@meta.data)))
pbmc@meta.data$Irf4 <- rep("neg",length(rownames(pbmc@meta.data)))
pbmc@meta.data$Ly75 <- rep("neg",length(rownames(pbmc@meta.data)))
pbmc@meta.data$Kdm6b <- rep("neg",length(rownames(pbmc@meta.data)))

for (i in pbmc_Irf4) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(i == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$Irf4[j] <- "Irf4"
    }
  }
}

for (i in pbmc_Kdm6b) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(i == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$Kdm6b[j] <- "Kdm6b"
    }
  }
}

for (i in pbmc_Ly75) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(i == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$Ly75[j] <- "Ly75"
    }
  }
}

for (i in pbmc_Myc) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(i == rownames(pbmc@meta.data)[j]){
      pbmc@meta.data$Myc[j] <- "Myc"
    }
  }
}

pbmc@meta.data$Mix<- paste(pbmc@meta.data$Myc,pbmc@meta.data$Irf4,pbmc@meta.data$Kdm6b,pbmc@meta.data$Ly75,sep = "_")
# Idents(pbmc)<-pbmc$Mix
# VlnPlot(pbmc,features = c("rna_Myc"),pt.size = 0,sort = TRUE) # + theme(axis.text.x = element_text(angle = 90))
# VlnPlot(pbmc,features = c("rna_Irf4"),pt.size = 0,sort = TRUE) # + theme(axis.text.x = element_text(angle = 90))
# VlnPlot(pbmc,features = c("rna_Kdm6b"),pt.size = 0,sort = TRUE) # + theme(axis.text.x = element_text(angle = 90))
# VlnPlot(pbmc,features = c("rna_Ly75"),pt.size = 0,sort = TRUE)  #+ theme(axis.text.x = element_text(angle = 90))
# VlnPlot(pbmc,features = c("Mtor_1"),pt.size = 0,sort = TRUE)  #+ theme(axis.text.x = element_text(angle = 90))
# VlnPlot(pbmc,features = c("Mtor_2"),pt.size = 0,sort = TRUE)  #+ theme(axis.text.x = element_text(angle = 90))
# VlnPlot(pbmc,features = c("Glutamine_1"),pt.size = 0,sort = TRUE)  #+ theme(axis.text.x = element_text(angle = 90))
# table(Idents(pbmc))

# elements <- unique(Idents(pbmc))[2:16]
# for (i in 1:15 ) {
#   all_combinations <- combn(elements, m = i, simplify = FALSE)
#   for (j in 1:length(all_combinations )) {
#     Myc1<-subset(pbmc,idents = all_combinations[[j]])
#     if (length(Myc1@meta.data[which(Myc1@meta.data$Ly75 == "pos"),])/length(rownames(Myc1@meta.data)) > 0.7){
#       if(length(Myc1@meta.data[which(Myc1@meta.data$Irf4 == "pos"),])/length(rownames(Myc1@meta.data)) > 0.7){
#         if(length(Myc1@meta.data[which(Myc1@meta.data$Ly75 == "pos"),])/length(rownames(Myc1@meta.data)) > 0.7){
#           if(length(Myc1@meta.data[which(Myc1@meta.data$Myc == "pos"),])/length(rownames(Myc1@meta.data)) > 0.7){
#             print(all_combinations[[j]])
#           }
#         }
#       }
#     }
#   }
# }
remove(all_combinations,control,Data,Data_T,Glutamine_1,Glutamine1,mtor_1,mtor_2,mtor1,mtor2,mtor_1_new,mtor_2_new,
       i,j,k,logene,higene,pbmc_Irf4,pbmc_Kdm6b,pbmc_Ly75,pbmc_Myc)
remove(p,elements,Glutamine_1_new,average_Glutamine1,average_mtor1,average_mtor2,all.genes,Total,Sortdata,Sortdatahi,Sortdatalo)
remove(rownamesDataT,Myc1,pbmc_Totalgene)

library(data.table)
library(dplyr)
library(ggplot2)
data<-pbmc@assays$RNA@data[c("Myc","Irf4","Kdm6b","Ly75"),]
data<-as.data.frame(t(as.data.frame(data)))
data<-cbind(data,pbmc@meta.data[rownames(data),][,c("Mtor_1","Mtor_2","Glutamine_1")])


num_genes <- 6
num_cells <- 15993
num_groups <- 11


population_size <- 100
max_generations <- 100
crossover_rate <- 1
mutation_rate <- 1


initialize_population <- function() {
  population <- matrix(0, nrow = num_cells, ncol = population_size)
  for (i in 1:population_size) {
    population[, i] <- sample(1:num_groups, num_cells,prob = c(0.15,0.14,0.1,0.1,0.1,0.09,0.08,0.07,0.06,0.06,0.04), replace = TRUE) 
  }
  return(population)
}
calculate_fitness <- function(population) {
  state<-cbind(data,population)
  gene_expression_means  <- aggregate(. ~ population, data = state, FUN = mean)
  max_expression_groups <- apply(gene_expression_means[,-1],2,function (x){return(gene_expression_means$population[which.max(x)])})
  max_expression_groups<-data.frame(max_expression_groups)
  colnames(max_expression_groups)<-"Group"
  fitness <- all(max_expression_groups$Group == max_expression_groups$Group[1])
  if(fitness){
    cells_above_zero <-  apply(state[which(state$population == max_expression_groups$Group[1]),c("Myc","Irf4","Kdm6b","Ly75")],2,function(gene) sum(gene>0))
    target_cells_percentage <- cells_above_zero / nrow(state[which(state$population == max_expression_groups$Group[1]),c("Myc","Irf4","Kdm6b","Ly75")])
    gene_expression_means_2<-gene_expression_means[order(gene_expression_means$Myc,decreasing = TRUE),]
    cells_above_zero2 <-  apply(data.frame(state[which(state$population ==gene_expression_means_2$population[2]),c("Myc")]),2,function(gene) sum(gene>0))
    target_cells_percentage2 <- cells_above_zero2 / nrow(data.frame(state[which(state$population ==gene_expression_means_2$population[2]),c("Myc")]))
    max2_expression_groups<-as.data.frame(apply(gene_expression_means[,2:8],2 , function(x){ return(gene_expression_means[order(x,decreasing = TRUE),][2,1])}))
    colnames(max2_expression_groups)<-"Group"
    fitness2 <- any(max2_expression_groups$Group != max2_expression_groups$Group[1])
    fitness <- fitness * fitness2
    fitness <- fitness * min(target_cells_percentage2,target_cells_percentage) +1
    
    
  }else{
    fitness <- 1
  }
  
  
  return(fitness)
}


selection <- function(population, fitness) {
  
  selected_indices <- sample(1:dim(population)[2], size = 10, prob = fitness / sum(fitness),replace = TRUE) #selection method need correction
  selected <- population[, selected_indices, drop = FALSE]
  return(selected)
}


crossover <- function(parent1, parent2) {
  
  crossover_point <- sample(1:15992, size = 1)
  child1 <- c(parent1[1:crossover_point], parent2[(crossover_point + 1):length(rownames(as.data.frame(parent1)))])
  child2 <- c(parent2[1:crossover_point], parent1[(crossover_point + 1):length(rownames(as.data.frame(parent1)))])
  return(list(child1, child2))
}


mutation <- function(individual) {
  
  mutation_point <- sample(1:length(rownames(as.data.frame(individual))), size = 1)
  individual[mutation_point] <- as.numeric(sample(1:num_groups, size = 1))
  return(individual)
}

t<-0
population <- initialize_population()
repeat{
  
  fitness <- apply(population, 2, calculate_fitness)
  if(any(fitness > 1.5)){
    break
  }
  new_population <- matrix(0, nrow = num_cells, ncol = 10)
  
  
  selected <- selection(population, fitness)
  
  
  for (i in seq(1, 10, by = 2)) {
    if (runif(1) < crossover_rate) {
      children <- crossover(selected[, i], selected[, i + 1])
      new_population[, i] <- children[[1]]
      new_population[, i + 1] <- children[[2]]
    } else {
      new_population[, i] <- selected[, i]
      new_population[, i + 1] <- selected[, i + 1]
    }
  }
  
  
  for (i in 1:seq(1, 10, by = 2)) {
    if (runif(1) < mutation_rate) {
      new_population[, i] <- mutation(new_population[, i])
    }
  }
  t=t+1
  population_new<-initialize_population()
  population <- cbind(population_new,new_population)
}


best_individual <- population[, which.max(apply(population, 2, calculate_fitness)), drop = FALSE]