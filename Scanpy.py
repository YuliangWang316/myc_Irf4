C:/Users/xjmik/anaconda3/python.exe
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce
import scrublet as scr
import numpy as np
import gseapy as gp
import sigc
import csv
import scvelo
np.random.seed(1)

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

adata = sc.read_10x_mtx(
    'D:/GCB1_2.2.filtered_feature_bc_matrix/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)    
sc.external.pp.scrublet(adata)
sce.pl.scrublet_score_distribution(adata)
adata = adata[adata.obs['predicted_doublet']==False, :]
sc.pl.highest_expr_genes(adata, n_top=20, )
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('mt-') 
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")
 # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=True, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo','pct_counts_hb'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
adata = adata[adata.obs['n_genes_by_counts'] < 6000]
adata = adata[adata.obs['pct_counts_mt']<30]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='auto')
sc.pl.pca(adata, color='Irf4')
sc.pl.pca_variance_ratio(adata, log=True,n_pcs=50)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40, metric= "cosine")
sc.tl.leiden(adata,resolution=1.2,n_iterations=2)
sc.tl.louvain(adata,flavor="igraph")
sc.tl.paga(adata)
sc.tl.dpt(adata)
sc.pl.paga(adata, plot=True)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')
sc.pl.umap(adata, color=['Myc', 'Irf4', 'Kdm6b','Ly75','Gls'])
sc.pl.umap(adata, color=['leiden', 'louvain'])
sc.pl.violin(adata, ['Myc', 'Irf4', 'Kdm6b','Ly75'], groupby='louvain')
sc.pl.dotplot(adata, [ 'Irf4', 'Kdm6b','Ly75','Prdm1'], groupby='louvain')
sc.pl.dotplot(adata, ['Irf4', 'Kdm6b','Ly75','Prdm1'], groupby='leiden')

sc.pl.stacked_violin(adata, ['Myc', 'Irf4', 'Kdm6b','Ly75'], groupby='louvain', rotation=90);
sc.pl.heatmap(adata, ['Myc', 'Irf4', 'Kdm6b','Ly75'], groupby='louvain', cmap='viridis', dendrogram=True)
sc.pl.dotplot(adata, ['Prdm1'], groupby='louvain')
with open("d:/LZ.txt") as f:
    LZ = [x.strip() for x in list(f)]

with open("d:/DZ.txt") as f:
    DZ = [x.strip() for x in list(f)]

with open("d:/GSE60927.txt") as f:
    Plasma = [x.strip() for x in list(f)]

with open("d:/Glutamine.txt") as f:
    Glutamine = [x.strip() for x in list(f)]

with open("d:/Mtor1.txt") as f:
    Mtor1 = [x.strip() for x in list(f)]

with open("d:/Mtor2.txt") as f:
    Mtor2 = [x.strip() for x in list(f)]




sc.tl.score_genes(adata,LZ,score_name="LZ")
sc.tl.score_genes(adata,DZ,score_name="DZ")
sc.tl.score_genes(adata,Plasma,score_name="Plasma")
sc.tl.score_genes(adata,Glutamine,score_name="Glutamine")
sc.tl.score_genes(adata,Mtor1,score_name="Mtor1")
sc.tl.score_genes(adata,Mtor2,score_name="Mtor2")

res = gp.gsva(data=adata.to_df().T, # row -> genes, column-> samples
        gene_sets="D:/Total.gmt",
        outdir="d:/GSVA",
        min_size=0,
        seed=1,
        threads=20)
a=adata.obs
a.to_csv("d:/GCB1metadata.txt", sep='\t')

sc.pl.umap(adata,color = "LZ")
sc.pl.umap(adata,color = "DZ")
sc.pl.umap(adata,color = "Plasma", vcenter='p55', vmin='p15',vmax='p95',sort_order = True,cmap="RdYlGn")
sc.pl.dotplot(adata, ['Myc', 'Irf4', 'Kdm6b','Ly75'], groupby='louvain')
sc.pl.dotplot(adata, ['Plasma'], groupby='leiden')
sc.pl.dotplot(adata, ['Plasma'], groupby='louvain')
sc.pl.dotplot(adata, ['DZ',"LZ"], groupby='louvain')
sc.pl.dotplot(adata, ["Glutamine"], groupby='louvain')
adata2 = sc.read_10x_mtx(
    'D:/GCB2_2.2.filtered_feature_bc_matrix/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)    
sc.external.pp.scrublet(adata2)
sce.pl.scrublet_score_distribution(adata2)
adata2 = adata2[adata2.obs['predicted_doublet']==False, :]
sc.pl.highest_expr_genes(adata2, n_top=20, )
sc.pp.filter_cells(adata2, min_genes=200)
sc.pp.filter_genes(adata2, min_cells=3)
adata2.var['mt'] = adata2.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata2, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata2, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata2, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata2, x='total_counts', y='n_genes_by_counts')
adata2 = adata2[adata2.obs['n_genes_by_counts'] < 6000]
adata2 = adata2[adata2.obs['pct_counts_mt']<35]
sc.pp.normalize_total(adata2, target_sum=1e4)
sc.pp.log1p(adata2)
sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata2.raw = adata2
adata2 = adata2[:, adata2.var.highly_variable]
sc.pp.regress_out(adata2, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata2, max_value=10)
sc.tl.pca(adata2, svd_solver='arpack')
sc.pl.pca(adata2, color='Myc')
sc.pl.pca_variance_ratio(adata2, log=True,n_pcs=50)
sc.pp.neighbors(adata2, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata2)
sc.tl.louvain(adata2)
sc.tl.paga(adata2)
sc.pl.paga(adata2, plot=True)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata2, init_pos='spectral')
sc.pl.umap(adata2, color=['Myc', 'Irf4', 'Kdm6b','Ly75'])
sc.pl.umap(adata2, color=['louvain','leiden'])
sc.pl.violin(adata2, ['Myc', 'Irf4', 'Kdm6b','Ly75'], groupby='leiden')
sc.pl.dotplot(adata2, ['Myc', 'Irf4', 'Kdm6b','Ly75','Prdm1'], groupby='louvain')
sc.pl.dotplot(adata2, ['Prdm1'], groupby='louvain')
sc.pl.stacked_violin(adata2, ['Myc', 'Irf4', 'Kdm6b','Ly75'], groupby='leiden', rotation=90);

sc.tl.score_genes(adata2,LZ,score_name="LZ")
sc.tl.score_genes(adata2,DZ,score_name="DZ")
sc.tl.score_genes(adata2,Plasma,score_name="Plasma")
sc.tl.score_genes(adata2,Glutamine,score_name="Glutamine")
sc.tl.score_genes(adata2,Mtor1,score_name="Mtor1")
sc.tl.score_genes(adata2,Mtor2,score_name="Mtor2")

sc.pl.umap(adata2,color = "LZ")
sc.pl.umap(adata2,color = "DZ")
sc.pl.umap(adata2,color = "Plasma",sort_order=True)
sc.pl.dotplot(adata2, ['Myc', 'Irf4', 'Kdm6b','Ly75'], groupby='louvain')
sc.pl.dotplot(adata2, ['Plasma'], groupby='louvain')
sc.pl.dotplot(adata2, ['DZ',"LZ"], groupby='louvain')
sc.pl.dotplot(adata2, ['Mtor1',"Mtor2","Glutamine"], groupby='louvain')

res2 = gp.gsva(data=adata2.to_df().T, # row -> genes, column-> samples
        gene_sets="D:/Total.gmt",
        outdir="d:/GSVA2",
        min_size=0,
        seed=1,
        threads=20)

b=adata2.obs
b.to_csv("d:/GCB2metadata.txt", sep='\t')


var_names = adata2.var_names.intersection(adata.var_names)
adata2 = adata2[:, var_names]
adata = adata[:, var_names]

sc.tl.ingest(adata, adata2, obs='louvain')
adata.uns['louvain_colors'] = adata_ref.uns['louvain_colors']