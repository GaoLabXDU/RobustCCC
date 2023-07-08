import os
import numpy as np
import pandas as pd
from read_results_and_extract_LRlist_py import read_result, extract_LRlist

def aggregate_result(path_result_CCC, name_mat, cellType1, cellType2, list_method, top):
    top = int(top)
    dict_aggregate = {}
    dict_i = 0
    for method in list_method:
        # read results and extract list of LR
        name_result = name_mat
        table_result = read_result(path_result_CCC, name_result, method)
        if len(table_result) == 0:
            list_LR = np.array([])
        else:
            list_LR = extract_LRlist(table_result, cellType1, cellType2, method)

        # extract top ligand-receptor pairs and merge
        list_LR_top = list_LR[:min(top, len(list_LR))]
        table_temp = pd.DataFrame({'method': [method] * len(list_LR_top),
                                   'cellType1': [cellType1] * len(list_LR_top),
                                   'cellType2': [cellType2] * len(list_LR_top),
                                   'list_LR_top': list_LR_top})
        table_temp.index = list(range(dict_i, dict_i + len(list_LR_top)))
        dict_temp = table_temp.to_dict(orient='index')
        dict_aggregate.update(dict_temp)
        dict_i = dict_i + len(table_temp)
    table_aggregate = pd.DataFrame.from_dict(dict_aggregate, 'index')
    table_aggregate.to_csv(os.path.join(path_result_CCC, 'aggregate_result_top'+str(top) +'.txt'))
