import pandas as pd
import scanpy as sc
import scConnect as cn
import os
import anndata as ad
import multiprocessing as mp

def run_scConnect(mat_ori, label_ori, path_result_CCC, name_mat):
    path_write_curr = os.path.join(path_result_CCC,'scConnect')
    if not os.path.isdir(path_write_curr):
        os.mkdir(path_write_curr)
    path_write_curr_error = os.path.join(path_write_curr, 'ERROR')
    if not os.path.isdir(path_write_curr_error):
        os.mkdir(path_write_curr_error)
    if os.path.isfile(os.path.join(path_write_curr, name_mat)):
        return(0)
    if os.path.isfile(os.path.join(path_write_curr_error, name_mat)):
        return(0)
    try:
        adata = ad.AnnData(mat_ori).T
        # label_ori.index = label_ori['Cell']
        adata.obs = label_ori
        
        # gene calling
        adata_tissue = cn.genecall.meanExpression(adata, groupby="cell_type", normalization=False, use_raw=False,
                                                  transformation="log1p")
        adata_tissue = cn.connect.ligands(adata_tissue,organism='hsapiens')
        adata_tissue = cn.connect.receptors(adata_tissue,organism='hsapiens')

        # Calculating specificity of each ligand and recepor
        adata_tissue = cn.connect.specificity(adata_tissue, n=100, groupby="cell_type",organism='hsapiens')
        # Constructing the graph
        edges = cn.connect.interactions(emitter=adata_tissue, target=adata_tissue, self_reference=False,organism='hsapiens')
        table = pd.concat([pd.DataFrame({'Lcell': [edge[0]],
                                         'Rcell': [edge[1]],
                                         'interaction': [edge[2]['interaction']],
                                         'ligand': [edge[2]['ligand']],
                                         'ligand_zscore': [edge[2]['ligand_zscore']],
                                         'ligand_pval': [edge[2]['ligand_pval']],
                                         'receptor': [edge[2]['receptor']],
                                         'receptor_zscore': [edge[2]['receptor_zscore']],
                                         'receptor_pval': [edge[2]['receptor_pval']],
                                         'score': [edge[2]['score']],
                                         'log_score': [edge[2]['log_score']],
                                         'specificity': [edge[2]['specificity']],
                                         'importance': [edge[2]['importance']]
                                         }) for edge in edges], axis=0)
        ###
        table.index = list(range(len(table)))
        table.to_csv(os.path.join(path_write_curr, name_mat))
    
    except Exception as e:  # most generic exception you can catch
        logf = open(os.path.join(path_write_curr_error, name_mat), "w")
        logf.write(str(e))
        logf.close()
        return(0)

