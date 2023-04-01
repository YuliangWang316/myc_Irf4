library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

counts <- Read10X_h5(filename = "D:/Bcells_singlecellATAC/aggr/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "D:/Bcells_singlecellATAC/aggr/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'D:/Bcells_singlecellATAC/aggr/outs/fragments.tsv.gz',
  min.cells = 1,
  
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(pbmc) <- annotations
pbmc <- NucleosomeSignal(object = pbmc)
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 1 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
brain <- pbmc
remove(pbmc)
brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- RunSVD(object = brain)
DepthCor(brain)

brain <- RunUMAP(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
S<-data.frame(0,0)
colnames(S)<-c("i","j")

for (i in 1:50) {
  for (j in seq(0.1,3,by=0.1)) {
brains <- FindNeighbors(
  object = brain,
  reduction = 'lsi',
  dims = 2:i
)
brains <- FindClusters(
  object = brains,
  algorithm = 3,
  resolution = j,
  verbose = FALSE
)

# DimPlot(object = brain, label = TRUE) + NoLegend()
gene.activities <- GeneActivity(brains)
pbmc.data<-as.data.frame(gene.activities)

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB1_filtered", min.cells = 3, min.features = 200)
rm(pbmc.data)
pbmcs <- NormalizeData(pbmc)




    
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
    
    
    
    
    
    if (levels(Kdm6b$ident)[1] == levels(Irf4$ident)[1]){
      if (levels(Kdm6b$ident)[1] == levels(Ly75$ident)[1]){
        if(levels(Kdm6b$ident)[1] == levels(Myc$ident)[1]){
          
          
              
              P<-data.frame(i,j)
              colnames(P)<-c("i","j")
              S<-rbind(S,P)
            
        }
      }
    }
  }
}

