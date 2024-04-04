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
pbmc <- NormalizeData(pbmc,normalization.method = "RC")
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

# 定义强化学习环境
rl_env <- function(data, target_genes) {
  state_space <- data[,target_genes]
  action_space <- c("A","B","C","D","E","F","G","H","I","J","K")
  
  reset <- function() {
    state <- state_space
    state$Group <- sample(action_space, nrow(state), replace = TRUE)
    return(state)
  }
  
  step <- function(state, action) {
    # 将一个细胞分到某个群体
    state$Group[state$Group == action] <- sample(action_space, sum(state$Group == action), replace = TRUE) ###这段需要重写
    # # 将剩下的细胞分到剩下的群体
    # remaining_groups <- setdiff(action_space, action)
    # state$Group[state$Group %in% remaining_groups] <- sample(remaining_groups, sum(state$Group %in% remaining_groups), replace = TRUE)
    gene_expression_means  <- aggregate(. ~ Group, data = state, FUN = mean)
    max_expression_groups <- apply(gene_expression_means[,-1],2,function (x){return(gene_expression_means$Group[which.max(x)])})
    max_expression_groups<-data.frame(max_expression_groups)
    colnames(max_expression_groups)<-"Group"
    reward <- all(max_expression_groups$Group == max_expression_groups$Group[1])
    
    cells_above_zero <- apply(state[, target_genes], 2, function(gene) sum(gene > 0))
    target_cells_percentage <- cells_above_zero / nrow(state)
    
    reward <- reward * min(target_cells_percentage)
    
    return(list(state = state, reward = reward))
  }
  
  return(list(reset = reset, step = step,action_space=action_space))
}


# 创建强化学习环境
env <- rl_env(data = data.frame(data), target_genes = c("Myc","Irf4","Kdm6b","Ly75","Mtor_1","Glutamine_1"))

q_learning <- function(env, num_episodes, alpha, gamma, epsilon) {
  Q <- matrix(0, nrow = length(env$action_space), ncol = length(env$reset()))
  
  for (episode in 1:num_episodes) {
    state <- env$reset()
    done <- FALSE
    
    while (!done) {
      # 选择动作
      if (length(env$action_space) > 0) {
        if (runif(1) < epsilon) {
          action <- sample(env$action_space, 1)
        } else {
          state_index <- as.numeric(paste(state, collapse = ""))
          action <- env$action_space[which.max(Q[, state_index])]
        }
      } else {
        action <- NULL  # No valid actions
      }
      
      # 执行动作并获取奖励
      result <- env$step(state, action)
      new_state <- result$state
      new_state_index <- as.numeric(paste(new_state, collapse = ""))
      reward <- result$reward
      
      # 更新Q值
      if (!is.null(action)) {
        Q[action, state_index] <- Q[action, state_index] + alpha * (reward + gamma * max(Q[, new_state_index]) - Q[action, state_index])
      }
      
      # 更新状态
      state <- new_state
      
      # 判断是否完成
      done <- sum(state) == length(env$reset())
    }
  }
  
  return(Q)
}


# 训练Q-learning代理
Q <- q_learning(env, num_episodes = 1000, alpha = 0.1, gamma = 0.9, epsilon = 0.1)

# 获取最优策略
optimal_policy <- apply(Q, 2, which.max)

