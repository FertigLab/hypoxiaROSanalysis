import pickle
import scanpy as sc
import pandas as pd

with open("./data/spatialhypoxia4patterns30k.pkl", 'rb') as fp:
    hypoxiaresult4pattern = pickle.load(fp)

# load CoGAPS result object
path = "/data/outs/filtered_feature_bc_matrix/matrix.mtx"
rawdata = sc.read_mtx(path)
rawdata.X=rawdata.X.todense()
features = pd.read_csv("/data/outs/filtered_feature_bc_matrix/features.tsv", sep="\t", header=None)
cell_labels = pd.read_csv("/data/outs/filtered_feature_bc_matrix/barcodes.tsv",
                       sep="\t", header=None)
rawdata.var_names=cell_labels[0]
rawdata.obs_names=features[1]

sc.pp.log1p(rawdata)

adata = hypoxiaresult4pattern.T
sc.pp.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]

# scale data and compute PCA
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

# find neighbor embeddings and run UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

patterns = list(adata.obs.columns)
# plot pattern amplitude on UMAP
sc.pl.umap(adata, color=patterns)

from PyCoGAPS.analysis_functions import *

pm = patternMarkers(hypoxiaresult4pattern, threshold="cut")

hypoxia_hallmarks = pd.read_table("hypoxiagenes", header=None)

genes = list(hypoxia_hallmarks[0])
genes = list(set(hypoxiaresult4pattern.obs_names).intersection(genes))