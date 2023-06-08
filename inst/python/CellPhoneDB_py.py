import pandas as pd
import os
from cellphonedb.utils import db_utils, file_utils

from typing import Tuple
import numpy as np
import pickle
from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.exceptions.AllCountsFilteredException import AllCountsFilteredException
from cellphonedb.src.core.exceptions.NoInteractionsFound import NoInteractionsFound
from cellphonedb.src.core.methods import cpdb_statistical_analysis_helper, cpdb_statistical_analysis_complex_method
from cellphonedb.src.core.models.complex import complex_helper

import warnings

warnings.simplefilter("ignore", UserWarning)
from functools import partial
from multiprocessing.pool import Pool
import numpy_groupies as npg
from cellphonedb.src.core.core_logger import core_logger
from cellphonedb.src.core.models.complex import complex_helper
from tqdm.std import tqdm


def cpdb_statistical_analysis_helper_shuffled_analysis(iterations: int,
                                                       meta: pd.DataFrame,
                                                       counts: pd.DataFrame,
                                                       interactions: pd.DataFrame,
                                                       cluster_combinations: list,
                                                       complex_to_protein_ids: dict,
                                                       real_mean_analysis: pd.DataFrame,
                                                       threads: int,
                                                       separator: str) -> list:
    """
    Shuffles meta and calculates the means for each and saves it in a list.

    If threads > 1, runs it in a multiple threads to run it faster
    """
    results = []
    if threads == 1:
        # NB. At the time of writing, our parallelisation method does not work when Python is called
        # via Reticulate from R Studio on Windows. The only option then is to run with threads set to 1.
        # See: https://github.com/ventolab/CellphoneDB/issues/102
        progress_step = round(iterations / 100, 0)
        for i in range(iterations):
            result = cpdb_statistical_analysis_helper._statistical_analysis(cluster_combinations,
                                                                            counts,
                                                                            interactions,
                                                                            meta,
                                                                            complex_to_protein_ids,
                                                                            separator,
                                                                            real_mean_analysis,
                                                                            i)

            # if i % progress_step == 0:
            #     # Poor man's progress reporting
            #     print("{}%".format(round(i / progress_step), 0), end=' ')
            results.append(result)
    else:
        with Pool(processes=threads) as pool:
            statistical_analysis_thread = partial(cpdb_statistical_analysis_helper._statistical_analysis,
                                                  cluster_combinations,
                                                  counts,
                                                  interactions,
                                                  meta,
                                                  complex_to_protein_ids,
                                                  separator,
                                                  real_mean_analysis)
            for result in tqdm(pool.imap(statistical_analysis_thread, range(iterations)),
                               total=iterations):
                results.append(result)
    return results


def cpdb_statistical_analysis_complex_method_call(meta: pd.DataFrame,
                                                  counts: pd.DataFrame,
                                                  counts_data: str,
                                                  interactions: pd.DataFrame,
                                                  genes: pd.DataFrame,
                                                  complexes: pd.DataFrame,
                                                  complex_compositions: pd.DataFrame,
                                                  microenvs: pd.DataFrame,
                                                  pvalue: float,
                                                  separator: str = '|',
                                                  iterations: int = 1000,
                                                  threshold: float = 0.1,
                                                  threads: int = 4,
                                                  debug_seed: int = -1,
                                                  result_precision: int = 3,
                                                  debug: bool = False,
                                                  output_path: str = '',
                                                  ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # core_logger.info(
    #     '[Cluster Statistical Analysis] '
    #     'Threshold:{} Iterations:{} Debug-seed:{} Threads:{} Precision:{}'.format(threshold,
    #                                                                               iterations,
    #                                                                               debug_seed,
    #                                                                               threads,
    #                                                                               result_precision))
    if debug_seed >= 0:
        np.random.seed(debug_seed)
        # core_logger.warning('Debug random seed enabled. Set to {}'.format(debug_seed))

    # get reduced interactions (drop duplicates)
    interactions_reduced = interactions[['multidata_1_id', 'multidata_2_id']].drop_duplicates()

    # add multidata id and means to counts
    counts, counts_relations = cpdb_statistical_analysis_helper.add_multidata_and_means_to_counts(
        counts, genes, counts_data)

    if counts.empty:
        raise AllCountsFilteredException(hint='Are you using human data?')

    interactions_filtered, counts_filtered, complex_composition_filtered = \
        cpdb_statistical_analysis_helper.prefilters(interactions_reduced,
                                                    counts,
                                                    complexes,
                                                    complex_compositions)

    if interactions_filtered.empty:
        raise NoInteractionsFound()

    meta = meta.loc[counts.columns]

    complex_to_protein_row_ids = complex_helper.map_complex_to_protein_row_ids(complex_composition_filtered,
                                                                               counts_filtered)
    clusters = cpdb_statistical_analysis_helper.build_clusters(meta, counts_filtered, complex_to_protein_row_ids,
                                                               skip_percent=False)
    # core_logger.info('Running Real Analysis')
    cluster_combinations = cpdb_statistical_analysis_helper.get_cluster_combinations(clusters['names'], microenvs)
    base_result = cpdb_statistical_analysis_helper.build_result_matrix(interactions_filtered,
                                                                       cluster_combinations,
                                                                       separator)

    real_mean_analysis = cpdb_statistical_analysis_helper.mean_analysis(interactions_filtered,
                                                                        clusters,
                                                                        cluster_combinations,
                                                                        separator)

    real_percents_analysis = cpdb_statistical_analysis_helper.percent_analysis(clusters,
                                                                               threshold,
                                                                               interactions_filtered,
                                                                               cluster_combinations,
                                                                               separator)

    # core_logger.info('Running Statistical Analysis')
    statistical_mean_analysis = cpdb_statistical_analysis_helper_shuffled_analysis(iterations,
                                                                                   meta,
                                                                                   counts_filtered,
                                                                                   interactions_filtered,
                                                                                   cluster_combinations,
                                                                                   complex_to_protein_row_ids,
                                                                                   real_mean_analysis,
                                                                                   threads,
                                                                                   separator)

    result_percent = cpdb_statistical_analysis_helper.build_percent_result(real_mean_analysis,
                                                                           real_percents_analysis,
                                                                           statistical_mean_analysis,
                                                                           interactions_filtered,
                                                                           cluster_combinations,
                                                                           base_result,
                                                                           separator)

    if debug:
        with open(f"{output_path}/debug_intermediate.pkl", "wb") as fh:
            pickle.dump({
                "genes": genes,
                "interactions": interactions,
                "interactions_filtered": interactions_filtered,
                "interactions_reduced": interactions_reduced,
                "complex_compositions": complex_compositions,
                "counts": counts,
                "counts_relations": counts_relations,
                "clusters_means_percents": clusters,
                "cluster_combinations": cluster_combinations,
                "base_result": base_result,
                "real_mean_analysis": real_mean_analysis,
                "real_percent_analysis": real_percents_analysis,
                "statistical_mean_analysis": statistical_mean_analysis,
                "result_percent": result_percent}, fh)

    pvalues_result, means_result, significant_means, deconvoluted_result = cpdb_statistical_analysis_complex_method.build_results(
        interactions_filtered,
        interactions,
        counts_relations,
        real_mean_analysis,
        result_percent,
        clusters['means'],
        complex_composition_filtered,
        counts,
        genes,
        result_precision,
        pvalue,
        counts_data
    )
    return pvalues_result, means_result, significant_means, deconvoluted_result


def run_cellphoneDB(mat_ori, label_ori, path_result_CCC, path_cpdb, name_mat):
    path_write_temp = os.path.join(path_result_CCC, 'CellPhoneDB')
    if not os.path.isdir(path_write_temp): os.mkdir(path_write_temp)
    path_write_curr = os.path.join(path_write_temp, name_mat)
    if not os.path.isdir(path_write_curr): os.mkdir(path_write_curr)
    if os.path.isfile(path_write_curr + '//' + 'statistical_analysis_significant_means_.txt'):
        return (0)

    # label_ori.index = label_ori['Cell']

    cpdb_file_path = '//'.join([path_cpdb, 'cellphonedb.zip'])
    counts_data = 'gene_name'
    output_path = path_write_curr
    microenvs_file_path = None
    iterations = 1000
    threshold = 0.1
    threads = 1
    debug_seed = -1
    result_precision = 3
    pvalue = 0.05
    separator: str = '|'
    debug = False
    output_suffix = ''

    # Load into memory CellphoneDB data
    interactions, genes, complex_composition, complex_expanded, gene_synonym2gene_name = \
        db_utils.get_interactions_genes_complex(cpdb_file_path)

    # Load user files into memory
    counts = mat_ori
    meta = label_ori
    microenvs = pd.DataFrame()
    degs = pd.DataFrame()

    pvalues, means, significant_means, deconvoluted = \
        cpdb_statistical_analysis_complex_method_call(meta.copy(),
                                                      counts.copy(),
                                                      counts_data,
                                                      interactions,
                                                      genes,
                                                      complex_expanded,
                                                      complex_composition,
                                                      microenvs,
                                                      pvalue,
                                                      separator,
                                                      iterations,
                                                      threshold,
                                                      threads,
                                                      debug_seed,
                                                      result_precision,
                                                      debug,
                                                      output_path
                                                      )
    max_rank = significant_means['rank'].max()
    significant_means['rank'] = significant_means['rank'].apply(lambda rank: rank if rank != 0 else (1 + max_rank))
    significant_means.sort_values('rank', inplace=True)
    file_utils.save_dfs_as_tsv(output_path, output_suffix, "statistical_analysis",
                               {"deconvoluted": deconvoluted,
                                "means": means,
                                "pvalues": pvalues,
                                "significant_means": significant_means})
