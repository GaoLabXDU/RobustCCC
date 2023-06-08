# Data preprocessing
import os
import pandas as pd

def load_mo2hu(path_code_depend):
    table = pd.read_csv(os.path.join(path_code_depend,'HOM_MouseHumanSequence.txt'), sep='\t')
    geneList_hu = list(table['Human_symbol'])
    geneList_mo = list(table['Mouse_symbol'])
    dict_mo2hu = {}
    for i in range(len(geneList_hu)):
        gene_hu = geneList_hu[i]
        gene_mo = geneList_mo[i]
        geneL_hu =gene_hu.split(',')
        geneL_mo = gene_mo.split(',')
        if(len(geneL_mo)==1):
            dict_mo2hu[gene_mo] = gene_hu
        else:
            for mo_i in range(len(geneL_mo)):
                dict_mo2hu[geneL_mo[mo_i]] = gene_hu
    return(dict_mo2hu)

# convert gene symbol, mouse -> human
def convert_gene_symbol_mouse2human(path_data, name_mat, path_code_depend):
    mat_mouse = pd.read_csv(os.path.join(path_data, name_mat), index_col=0)
    dict_mo2hu = load_mo2hu(path_code_depend)
    list_gene_mouse = list(mat_mouse.index)
    list_gene_human = []
    list_gene_index = []
    for i in range(len(list_gene_mouse)):
        if (list_gene_mouse[i] in dict_mo2hu):
            gene_human = dict_mo2hu[list_gene_mouse[i]]
            list_temp_human = gene_human.split(',')
            if (list_temp_human == 1):
                list_gene_human.append(gene_human)
                list_gene_index.append(i)
            else:
                for temp_human in list_temp_human:
                    list_gene_human.append(temp_human)
                    list_gene_index.append(i)
    mat_human = mat_mouse.iloc[list_gene_index, :]
    mat_human.index = list_gene_human
    mat_human = mat_human[~mat_human.index.duplicated(keep='last')]
    mat_human.to_csv(os.path.join(path_data,'_'.join([name_mat,'human'])))

def split_mat_to_cell_type_pairs(path_data, name_mat, name_label, ct1, ct2):
    mat = pd.read_csv(os.path.join(path_data, name_mat), index_col=0)
    label = pd.read_csv(os.path.join(path_data, name_label), index_col=0)
    label_ct1 = label.loc[label['cell_type']==ct1,]
    mat_ct1 = mat.loc[:,list(label_ct1.index)]
    label_ct2 = label.loc[label['cell_type']==ct2,]
    mat_ct2 = mat.loc[:,list(label_ct2.index)]
    mat_ct12 = pd.concat([mat_ct1,mat_ct2],axis=1)
    label_ct12 = pd.concat([label_ct1,label_ct2],axis=0)
    mat_ct12.to_csv(os.path.join(path_data, '_'.join([name_mat, ct1, ct2])), index=True)
    label_ct12.to_csv(os.path.join(path_data, '_'.join([name_label, ct1, ct2])), index=True)

