# similarity of communication between replicates
import os
import pandas as pd
import numpy as np
from read_results_and_extract_LRlist_py import read_result, extract_LRlist
import rbo
import warnings
warnings.filterwarnings("ignore")

def jaccard(set1, set2):
    down = len(set(set1) | set(set2))
    if down ==0:
        return(0)
    else:
        return (float(len(set(set1) & set(set2))) / down)

def top_or_nonTop(method):
    if method in ['CellPhoneDB','CellChat','iCellNet','iTALK','NATMI','scConnect',
                  'SingleCellSignalR','CellCall','CytoTalk','NicheNet','Kumar']:
        return('top')
    if method in ['scMLnet','Skelly','Zhou']:
        return('nonTop')

def similarity_between_replicates(path_result_CCC,path_result_robust, name_mat_1,name_mat_2, cellType1, cellType2, list_method):
    dict_similairity = {}
    dict_i = 0
    for method in list_method:
        topType = top_or_nonTop(method)
        # read results and extract list of LR from replicate 1
        name_result_1 = name_mat_1
        table_result_1 = read_result(path_result_CCC, name_result_1, method)
        if len(table_result_1) == 0:
            list_LR_1 = np.array([])
        else:
            list_LR_1 = extract_LRlist(table_result_1, cellType1, cellType2, method)

        # read results and extract list of LR from replicate 2
        name_result_2 = name_mat_2
        table_result_2 = read_result(path_result_CCC, name_result_2, method)
        if len(table_result_2) == 0:
            list_LR_2 = np.array([])
        else:
            list_LR_2 = extract_LRlist(table_result_2, cellType1, cellType2, method)
        # similarity between two results
        if topType == 'top':
            for top in [20000,'rbo']:
                if top != 'rbo':
                    dict_similairity[dict_i] = {'method': method, 'cellType1': cellType1, 'cellType2': cellType2,
                                                'similar_type': 'Jaccard',
                                                'numLR_R1': min(top, len(list_LR_1)),'numLR_R2': min(top, len(list_LR_2)),
                                                'similarity': jaccard(list_LR_1[:min(top, len(list_LR_1))],
                                                                   list_LR_2[:min(top, len(list_LR_2))])}
                else:
                    dict_similairity[dict_i] = {'method': method, 'cellType1': cellType1, 'cellType2': cellType2,
                                                'similar_type': 'RBO',
                                                'numLR_R1': len(list_LR_1),'numLR_R2': len(list_LR_2),
                                                'similarity': rbo.RankingSimilarity(list_LR_1, list_LR_2).rbo()}
                dict_i = dict_i + 1
        else:  # nonTop
            dict_similairity[dict_i] = {'method': method, 'cellType1': cellType1,'cellType2': cellType2,'similar_type': 'Jaccard',
                                        'numLR_R1': len(list_LR_1),'numLR_R2': len(list_LR_2),
                                        'similarity': jaccard(list_LR_1, list_LR_2)}
            dict_i = dict_i + 1
    table_similairity = pd.DataFrame.from_dict(dict_similairity, 'index')
    table_similairity.to_csv(os.path.join(path_result_robust, '_'.join([name_mat_1, name_mat_2,'bioReplicate'])))

def similarity_between_simulated_data(path_result_CCC,path_result_robust, simulated_type, name_mat, cellType1, cellType2, list_method, list_rate, num_of_rep):
    dict_similairity = {}
    dict_i = 0
    num_of_rep = int(num_of_rep)
    list_rate = [int(i) for i in list_rate]
    for method in list_method:
        topType = top_or_nonTop(method)
        for rate in list_rate:
            for rand_i in range(num_of_rep):
                # read results and extract list of LR from the 1st result file
                name_result_1 = '_'.join([name_mat,simulated_type,str(rate),str(rand_i)])
                table_result_1 = read_result(path_result_CCC, name_result_1, method)
                if len(table_result_1) ==0:
                    list_LR_1 = np.array([])
                else:
                    list_LR_1 = extract_LRlist(table_result_1, cellType1, cellType2, method)
                for rand_j in range(num_of_rep):
                    # read results and extract list of LR from the 2nd result file
                    if rand_i >= rand_j: continue
                    name_result_2 = '_'.join([name_mat, simulated_type, str(rate), str(rand_j)])
                    table_result_2 = read_result(path_result_CCC, name_result_2, method)
                    if len(table_result_2) == 0:
                        list_LR_2 = np.array([])
                    else:
                        list_LR_2 = extract_LRlist(table_result_2, cellType1, cellType2, method)
                    # similarity between two results
                    if topType == 'top':
                        for top in [20000,'rbo']:
                            if top != 'rbo':
                                dict_similairity[dict_i] = {'method':method,'cellType1':cellType1,'cellType2':cellType2,
                                                            'rate':rate, 'rand_i':rand_i, 'rand_j':rand_j, 'similar_type': 'Jaccard',
                                                            'numLR_R1': min(top, len(list_LR_1)),'numLR_R2': min(top, len(list_LR_2)),
                                                            'similarity': jaccard(list_LR_1[:min(top, len(list_LR_1))],
                                                                               list_LR_2[:min(top, len(list_LR_2))])}
                            else:
                                dict_similairity[dict_i] = {'method': method, 'cellType1': cellType1,'cellType2': cellType2,
                                                            'rate': rate, 'rand_i': rand_i, 'rand_j': rand_j,
                                                            'similar_type': 'RBO',
                                                            'numLR_R1': len(list_LR_1),
                                                            'numLR_R2': len(list_LR_2),
                                                            'similarity': rbo.RankingSimilarity(list_LR_1, list_LR_2).rbo()}
                            dict_i = dict_i + 1
                    else:  # nonTop
                        dict_similairity[dict_i] = {'method': method, 'cellType1': cellType1,
                                                    'cellType2': cellType2,
                                                    'rate': rate, 'rand_i': rand_i, 'rand_j': rand_j,
                                                    'similar_type': 'Jaccard',
                                                    'numLR_R1': len(list_LR_1),
                                                    'numLR_R2': len(list_LR_2),
                                                    'similarity': jaccard(list_LR_1, list_LR_2)}
                        dict_i = dict_i + 1
    table_similairity = pd.DataFrame.from_dict(dict_similairity, 'index')
    table_similairity.to_csv(os.path.join(path_result_robust, '_'.join([name_mat, simulated_type])))

def similarity_between_simulated_and_original_data(path_result_CCC,path_result_robust, simulated_type, name_mat, cellType1, cellType2, list_method, list_rate, num_of_rep):
    dict_similairity = {}
    dict_i = 0
    num_of_rep = int(num_of_rep)
    list_rate = [int(i) for i in list_rate]
    for method in list_method:
        topType = top_or_nonTop(method)
        # read results and extract list of LR from the 1st result file
        name_result_ori = name_mat
        table_result_ori = read_result(path_result_CCC, name_result_ori, method)
        if len(table_result_ori) == 0:
            list_LR_ori = np.array([])
        else:
            list_LR_ori = extract_LRlist(table_result_ori, cellType1, cellType2, method)
        for rate in list_rate:
            for rand_i in range(num_of_rep):
                # read results and extract list of LR from the 2nd result file
                name_result_simulate = '_'.join([name_mat, simulated_type, str(rate), str(rand_i)])
                table_result_simulate = read_result(path_result_CCC, name_result_simulate, method)
                if len(table_result_simulate) == 0:
                    list_LR_simulate = np.array([])
                else:
                    list_LR_simulate = extract_LRlist(table_result_simulate, cellType1, cellType2, method)

                # similarity between two results
                if topType == 'top':
                    for top in [20000,'rbo']:
                        if top != 'rbo':
                            dict_similairity[dict_i] = {'method':method,'cellType1':cellType1,'cellType2':cellType2,
                                                        'rate':rate, 'rand_i':rand_i, 'similar_type': 'Jaccard',
                                                        'numLR_simulate': min(top, len(list_LR_simulate)),'numLR_ori': min(top, len(list_LR_ori)),
                                                        'similarity': jaccard(list_LR_simulate[:min(top, len(list_LR_simulate))],
                                                                           list_LR_ori[:min(top, len(list_LR_ori))])}
                        else:
                            dict_similairity[dict_i] = {'method': method, 'cellType1': cellType1,'cellType2': cellType2,
                                                        'rate': rate, 'rand_i': rand_i,
                                                        'similar_type': 'RBO',
                                                        'numLR_simulate': len(list_LR_simulate),
                                                        'numLR_ori': len(list_LR_ori),
                                                        'similarity': rbo.RankingSimilarity(list_LR_simulate, list_LR_ori).rbo()}
                        dict_i = dict_i + 1
                else:  # nonTop
                    dict_similairity[dict_i] = {'method': method, 'cellType1': cellType1,
                                                'cellType2': cellType2,
                                                'rate': rate, 'rand_i': rand_i,
                                                'similar_type': 'Jaccard',
                                                'numLR_simulate': len(list_LR_simulate),
                                                'numLR_ori': len(list_LR_ori),
                                                'similarity': jaccard(list_LR_simulate, list_LR_ori)}
                    dict_i = dict_i + 1
    table_similairity = pd.DataFrame.from_dict(dict_similairity, 'index')
    table_similairity.to_csv(os.path.join(path_result_robust, '_'.join([name_mat, simulated_type])))
