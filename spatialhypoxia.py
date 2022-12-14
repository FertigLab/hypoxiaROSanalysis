if __name__ == "__main__":
    from PyCoGAPS import *
    import pickle
    from PyCoGAPS.parameters import *
    from PyCoGAPS.pycogaps_main import CoGAPS
    import scanpy as sc
    import pandas as pd
    import anndata
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

    params = CoParams(path)

    setParams(params, {
        'nIterations': 50000,
        'seed': 42,
        'nPatterns': 4,
        'useSparseOptimization': True,
        'distributed': "genome-wide"
    })

    params.setDistributedParams(nSets=8)
    params.printParams()
    start = time.time()
    result = CoGAPS(rawdata, params)
    end = time.time()
    print("TIME:", end - start)

    print("Pickling...")
    pickle.dump(result, open("./data/spatialhypoxia4patterns30k.pkl", "wb"))
    print("Pickling complete!")
