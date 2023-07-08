import os
import pandas as pd
import numpy as np
from scipy.stats import combine_pvalues
import warnings
warnings.filterwarnings("ignore")

def read_result(path_result, name_result, method):
    # file path is?
    if method == 'CellPhoneDB':
        path_file = os.path.join(path_result, method, name_result, 'statistical_analysis_pvalues_.txt')
    elif method == 'NATMI':
        path_file = os.path.join(path_result, method, name_result, 'Edges_lrc2p.csv')
    elif method == 'CytoTalk':
        path_file = os.path.join(path_result, method, name_result, 'Crosstalk_TypATypB.txt')
    else:
        path_file = os.path.join(path_result, method, name_result)

    # file exist?
    if not os.path.isfile(path_file):
        return(pd.DataFrame())

    # read result
    if method in ['CellPhoneDB', 'CytoTalk','CellCall']:
        table_result = pd.read_csv(path_file, sep='\t', index_col=False)  ## CellPhoneDB: pvalue<1;
    elif method == 'scConnect':
        table_result = pd.read_csv(path_file, index_col=0)
    else:
        table_result = pd.read_csv(path_file)

    if len(table_result) == 0:
        return (pd.DataFrame())
    return(table_result)

def extract_LRlist_cellPhoneDB(table_result, cellType1, cellType2):
    if len(table_result) == 0: return ([])
    col = cellType1 + '|' + cellType2
    cpdbTable_temp = table_result.loc[table_result.loc[:, col] < 0.05,]
    if len(cpdbTable_temp) ==0: return([])
    cpdbTable_sort = cpdbTable_temp.sort_values(by=col,ascending=True)
    LRlist = cpdbTable_sort['interacting_pair'].to_list()
    LRlist = np.array(LRlist)
    return (LRlist)

def extract_LRlist_CellChat(table_result):
    if len(table_result) ==0: return([])
    table_sort = table_result.sort_values(by='pval', ascending=True)
    LRlist = table_sort['LRname'].to_list()
    LRlist = np.array(LRlist)
    return (LRlist)

def extract_LRlist_iCellNet(table_result):
    if len(table_result) ==0: return([])
    table_sort = table_result.sort_values(by='score', ascending=False)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_iTALK(table_result, cellType1, cellType2):
    if len(table_result) ==0: return([])
    table_result = table_result.loc[table_result['cell_from'].str.contains(cellType1) &
                                    table_result['cell_to'].str.contains(cellType2),]
    if len(table_result) == 0: return ([])

    table_result['score'] = table_result['cell_from_mean_exprs'] * table_result['cell_to_mean_exprs']
    table_result['lr'] = table_result['ligand'] + '_' + table_result['receptor']
    table_result = table_result.loc[table_result['score'] > 0,]
    table_sort = table_result.sort_values(by='score', ascending=False)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_NATMI(table_result, cellType1, cellType2):
    if len(table_result) ==0: return([])
    table_result = table_result.loc[table_result['Sending cluster'].str.contains(cellType1) &
                                    table_result['Target cluster'].str.contains(cellType2),]
    if len(table_result) == 0: return ([])
    table_result['lr'] = table_result['Ligand symbol'] + '_' + table_result['Receptor symbol']
    table_sort = table_result.sort_values(by='Edge average expression weight', ascending=False)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_scConnect(table_result, cellType1, cellType2):
    if len(table_result) ==0: return([])
    table_result = table_result.loc[table_result['Lcell'].str.contains(cellType1) &
                                    table_result['Rcell'].str.contains(cellType2),]
    if len(table_result) == 0: return ([])
    table_sort = table_result.sort_values(by='importance', ascending=False)
    LRlist = table_sort['interaction'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_SingleCellSignalR(table_result):
    if len(table_result) ==0: return([])
    table_result['lr'] = table_result['cluster 1'] + '_' + table_result['cluster 2']
    table_sort = table_result.sort_values(by='LRscore',ascending=False)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_CellCall(table_result,cellType1, cellType2):
    if len(table_result) == 0: return ([])
    col = cellType1 + '-' + cellType2
    cpdbTable_temp = table_result.loc[table_result.loc[:, col] > 0,]
    if len(cpdbTable_temp) ==0: return([])
    cpdbTable_sort = cpdbTable_temp.sort_values(by=col,ascending=False)
    LRlist = cpdbTable_sort['LR'].to_list()
    LRlist = np.array(LRlist)
    return (LRlist)

def extract_LRlist_CytoTalk(table_result):
    if len(table_result) ==0: return([])
    table_result = table_result.loc[table_result['score'] > 0,]
    if len(table_result) == 0: return ([])
    table_result = table_result.loc[table_result['ligand'].str.contains('TypA')&
                                        table_result['receptor'].str.contains('TypB'),]
    if len(table_result) == 0: return ([])
    table_result['lr'] = table_result['ligand'] + '_' + table_result['receptor']
    table_sort = table_result.sort_values(by='score',ascending=False)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_NicheNet(table_result):
    if len(table_result) ==0: return([])
    table_result['lr'] = table_result['from'] + '_' + table_result['to']
    table_sort = table_result.sort_values(by='weight',ascending=False)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_scMLnet(table_result):
    if len(table_result) ==0: return([])
    list_Pva_Lig = table_result['Pva_Lig'].tolist()
    list_Pva_Rec = table_result['Pva_Rec'].tolist()
    table_result['commuPvalue'] = [combine_pvalues([list_Pva_Lig[pva_i], list_Pva_Rec[pva_i]])[1] for pva_i in
                                     range(len(list_Pva_Lig))]
    table_result['lr'] = table_result['Lig'] + '_' + table_result['Rec']
    table_sort = table_result.sort_values(by='commuPvalue', ascending=True)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_Kumar(table_result):
    if len(table_result) ==0: return([])
    table_result = table_result.loc[table_result['product_LR'] > 0,]
    if len(table_result) == 0: return ([])
    table_result['lr'] = table_result['Ligand'] + '_' + table_result['Receptor']
    table_sort = table_result.sort_values(by='product_LR', ascending=False)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_Skelly(table_result):
    if len(table_result) ==0: return([])
    table_result = table_result.loc[((table_result['rate_Lcells_Expressed']>0.2)&
                                         (table_result['rate_Rcells_Expressed']>0.2)),]
    if len(table_result) == 0: return ([])
    table_result['commuScore'] = table_result['rate_Lcells_Expressed'] * table_result['rate_Rcells_Expressed']
    table_result['lr'] = table_result['Ligand'] + '_' + table_result['Receptor']
    table_sort = table_result.sort_values(by='commuScore', ascending=False)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist_Zhou(table_result):
    if len(table_result) ==0: return([])
    table_result = table_result.loc[((table_result['LogFC_Lig']>1)&
                                         (table_result['LogFC_Rec']>1)&
                                         (table_result['Pva_Lig'] <0.05)&
                                         (table_result['Pva_Rec'] <0.05)),]
    if len(table_result) == 0: return ([])
    list_Pva_Lig = table_result['Pva_Lig'].tolist()
    list_Pva_Rec = table_result['Pva_Rec'].tolist()
    table_result['commuPvalue'] = [combine_pvalues([list_Pva_Lig[pva_i], list_Pva_Rec[pva_i]])[1] for pva_i in
                                   range(len(list_Pva_Lig))]
    table_result['lr'] = table_result['Ligand'] + '_' + table_result['Receptor']
    table_sort = table_result.sort_values(by='commuPvalue', ascending=True)
    LRlist = table_sort['lr'].to_list()
    LRlist = np.array(LRlist)
    return(LRlist)

def extract_LRlist(table_result,cellType1, cellType2, method):
    if method == 'CellPhoneDB':
        LRlist = extract_LRlist_cellPhoneDB(table_result, cellType1, cellType2)
    if method == 'CellChat':
        LRlist = extract_LRlist_CellChat(table_result)
    if method == 'ICELLNET':
        LRlist = extract_LRlist_iCellNet(table_result)
    if method == 'iTALK':
        LRlist = extract_LRlist_iTALK(table_result,cellType1, cellType2)
    if method == 'NATMI':
        LRlist = extract_LRlist_NATMI(table_result,cellType1, cellType2)
    if method == 'scConnect':
        LRlist = extract_LRlist_scConnect(table_result, cellType1, cellType2)
    if method == 'SingleCellSignalR':
        LRlist = extract_LRlist_SingleCellSignalR(table_result)
    if method == 'CellCall':
        LRlist = extract_LRlist_CellCall(table_result, cellType1, cellType2)
    if method == 'CytoTalk':
        LRlist = extract_LRlist_CytoTalk(table_result)
    if method == 'NicheNet':
        LRlist = extract_LRlist_NicheNet(table_result)
    if method == 'scMLnet':
        LRlist = extract_LRlist_scMLnet(table_result)
    if method == 'Kumar':
        LRlist = extract_LRlist_Kumar(table_result)
    if method == 'Skelly':
        LRlist = extract_LRlist_Skelly(table_result)
    if method == 'Zhou':
        LRlist = extract_LRlist_Zhou(table_result)
    return(np.unique(LRlist))
