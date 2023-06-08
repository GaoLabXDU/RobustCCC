# generating simulated data
import os
import pandas as pd
import numpy as np

def generate_simulated_data(simulated_type, path_data, list_rate, num_of_rep, name_mat, mat, name_label, label, ct1, ct2):
    if simulated_type == 'simuReplicate':
        generate_simulated_replicates(path_data, list_rate, num_of_rep, name_mat, mat, name_label, label, ct1, ct2)
    if simulated_type == 'GaussianNoise':
        generate_Gaussian_noise(path_data, list_rate, num_of_rep, name_mat, mat)
    if simulated_type == 'dropout':
        generate_dropout(path_data, list_rate, num_of_rep, name_mat, mat)
    if simulated_type == 'selfNoised':
        generate_self_noised_data(path_data, list_rate, num_of_rep, name_mat, mat)
    if simulated_type == 'cellTypePermu':
        generate_cell_type_permutaion(path_data, list_rate, num_of_rep, name_label, label, ct1, ct2)
    if simulated_type == 'ligRecPermu':
        generate_ligand_receptor_permutation(path_data, list_rate, num_of_rep, name_mat, mat)
    if simulated_type == 'expPermuted':
        generate_expression_permuted_data(path_data, list_rate, num_of_rep, name_mat, mat)


def generate_simulated_replicates(path_data, list_rate, num_of_rep, name_mat, mat, name_label, label, ct1, ct2):
    # generate matrix and label of cell type 1
    label_ct1 = label.loc[label['cell_type'] == ct1,]
    mat_ct1 = mat.loc[:, list(label_ct1.index)]
    list_cell_ct1 = list(mat_ct1.columns)

    # generate matrix and label of cell type 2
    label_ct2 = label.loc[label['cell_type'] == ct2,]
    mat_ct2 = mat.loc[:, list(label_ct2.index)]
    list_cell_ct2 = list(mat_ct2.columns)
    list_rate = [int(i) for i in list_rate]
    num_of_rep = int(num_of_rep)
    for rate in list_rate:
        for rand_i in range(num_of_rep):
            # sample cells in cell type 1
            np.random.seed(seed=rate*10+rand_i)
            list_cell_ct1_rand = np.random.choice(list_cell_ct1, replace=False,
                                                  size=int(len(list_cell_ct1) * (rate / 100)))
            list_cell_ct1_rand = sorted(list_cell_ct1_rand)
            mat_ct1_rand = mat_ct1.loc[:, list_cell_ct1_rand]

            # sample cells in cell type 2
            np.random.seed(seed=rate*10+rand_i)
            list_cell_ct2_rand = np.random.choice(list_cell_ct2, replace=False,
                                                  size=int(len(list_cell_ct2) * (rate / 100)))
            list_cell_ct2_rand = sorted(list_cell_ct2_rand)
            mat_ct2_rand = mat_ct2.loc[:, list_cell_ct2_rand]

            # merge sampled matrix of cell type 1 and 2, and save
            mat_rand = pd.concat([mat_ct1_rand, mat_ct2_rand],axis=1)
            mat_rand.to_csv(os.path.join(path_data, '_'.join([name_mat,'simuReplicate',str(rate),str(rand_i)])))

            # generate label of sampled cells of cell type 1 and 2, and save
            label_rand = pd.DataFrame({'Cell': list_cell_ct1_rand+list_cell_ct2_rand,
                                       'cell_type': [ct1] * len(list_cell_ct1_rand)+
                                                    [ct2] * len(list_cell_ct2_rand)})
            label_rand.to_csv(os.path.join(path_data, '_'.join([name_label,'simuReplicate',str(rate),str(rand_i)])),
                               header=True, index=False)
            ##


def generate_Gaussian_noise(path_data, list_rate, num_of_rep, name_mat, mat):
    list_cell = list(mat.columns)
    list_gene = list(mat.index)
    num_of_cells = len(list_cell)
    num_of_genes = len(list_gene)
    list_mean_genes = np.mean(mat, axis=1)
    list_std_genes = np.std(mat, axis=1)
    list_rate = [int(i) for i in list_rate]
    num_of_rep = int(num_of_rep)
    for rate in list_rate:
        for rand_i in range(num_of_rep):
            # generate noised mat and save
            np.random.seed(seed=rate * 10 + rand_i)
            noise = np.array([np.random.normal(list_mean_genes[row_i], list_std_genes[row_i], num_of_cells) for row_i in list(range(num_of_genes))])
            mat_noise = mat * (1 - rate / 100) + noise * (rate / 100)
            mat_noise[mat_noise < 0] = 0
            mat_noise.to_csv(os.path.join(path_data, '_'.join([name_mat,'GaussianNoise',str(rate),str(rand_i)])))

def generate_dropout(path_data, list_rate, num_of_rep, name_mat, mat):
    list_cell = list(mat.columns)
    list_gene = list(mat.index)
    list_rate = [int(i) for i in list_rate]
    num_of_rep = int(num_of_rep)
    for rate in list_rate:
        for rand_i in range(num_of_rep):
            # generate noised mat and save
            list_mat = np.array(mat).flatten()
            list_mat_idx_non0 = np.where(list_mat > 0)[0]
            np.random.seed(seed=rate * 10 + rand_i)
            list_select_idx_non0 = np.random.choice(list_mat_idx_non0, replace=False, size=int(len(list_mat_idx_non0) * (rate / 100)))
            list_mat[list_select_idx_non0] = [0] * len(list_select_idx_non0)
            mat_noise = list_mat.reshape(mat.shape)
            mat_noise = pd.DataFrame(mat_noise)
            mat_noise.index = list_gene
            mat_noise.columns = list_cell
            mat_noise.to_csv(os.path.join(path_data, '_'.join([name_mat, 'dropout', str(rate), str(rand_i)])))

def generate_self_noised_data(path_data, list_rate, num_of_rep, name_mat, mat):
    list_rate = [int(i) for i in list_rate]
    num_of_rep = int(num_of_rep)
    for rate in list_rate:
        for rand_i in range(num_of_rep):
            # generate noised mat and save
            np.random.seed(seed=rate * 10 + rand_i)
            noise = mat.copy()
            noise = np.array(noise).flatten()
            np.random.shuffle(noise)
            noise = noise.reshape(mat.shape)
            mat_noise = mat * (1 - rate / 100) + noise * (rate / 100)
            mat_noise[mat_noise < 0] = 0
            mat_noise.to_csv(os.path.join(path_data, '_'.join([name_mat,'selfNoised',str(rate),str(rand_i)])))


def generate_cell_type_permutaion(path_data, list_rate, num_of_rep, name_label, label, ct1, ct2):
    dict_permu = {ct1: ct2, ct2: ct1}
    list_cell = list(label.index)
    num_of_cells = len(label)
    label['Cell'] = list_cell
    label = label[['Cell','cell_type']]
    label.index = list(range(num_of_cells))
    list_rate = [int(i) for i in list_rate]
    num_of_rep = int(num_of_rep)
    for rate in list_rate:
        for rand_i in range(num_of_rep):
            np.random.seed(seed=rate * 10 + rand_i)
            list_cellIdx_permu = sorted(
                np.random.choice(list(range(num_of_cells)), replace=False, size=int(num_of_cells * (rate / 100))))
            list_ct_permu = [dict_permu[x] for x in label.loc[list_cellIdx_permu, 'cell_type'].tolist()]
            label_permu = label.copy()
            label_permu.loc[list_cellIdx_permu, 'cell_type'] = list_ct_permu
            label_permu.to_csv(
                os.path.join(path_data, '_'.join([name_label, 'cellTypePermu', str(rate), str(rand_i)])),
                header=True, index=False)

def generate_ligand_receptor_permutation(path_data, list_rate, num_of_rep, name_mat, mat):
    list_gene = np.array(list(mat.index))
    num_genes = len(list_gene)
    list_rate = [int(i) for i in list_rate]
    num_of_rep = int(num_of_rep)
    for rate in list_rate:
        for rand_i in range(num_of_rep):
            # generate noised mat and save
            np.random.seed(seed=rate * 10 + rand_i)
            list_geneIdx_select = sorted(np.random.choice(list(range(num_genes)), replace=False, size=int(num_genes * (rate / 100))))
            list_gene_select = list_gene[list_geneIdx_select]
            np.random.seed(seed=rate * 10 + rand_i)
            np.random.shuffle(list_gene_select)
            list_gene_permuted = list_gene.copy()
            list_gene_permuted[list_geneIdx_select] = list_gene_select
            mat_permuted = mat.copy()
            mat_permuted.index = list_gene_permuted
            mat_permuted.to_csv(os.path.join(path_data, '_'.join([name_mat, 'ligRecPermu', str(rate), str(rand_i)])))

def generate_expression_permuted_data(path_data, list_rate, num_of_rep, name_mat, mat):
    list_rate = [int(i) for i in list_rate]
    num_of_rep = int(num_of_rep)
    list_cell = list(mat.columns)
    list_gene = list(mat.index)
    for rate in list_rate:
        for rand_i in range(num_of_rep):
            list_permu = np.array(mat).flatten()
            list_permu_idx_0 = np.where(list_permu == 0)[0]
            list_permu_idx_non0 = np.where(list_permu > 0)[0]
            np.random.seed(seed=rate * 10 + rand_i)
            list_select_idx_0 = np.random.choice(list_permu_idx_0, replace=False,
                                                 size=int(len(list_permu_idx_0) * (rate / 100)))
            np.random.seed(seed=rate * 10 + rand_i)
            list_select_idx_non0 = np.random.choice(list_permu_idx_non0, replace=False,
                                                    size=int(len(list_permu_idx_non0) * (rate / 100)))
            list_select_idx = sorted(np.concatenate((list_select_idx_0, list_select_idx_non0)))
            list_permu_select = list_permu[list_select_idx]

            np.random.seed(seed=rate * 10 + rand_i)
            np.random.shuffle(list_permu_select)
            list_permu[list_select_idx] = list_permu_select
            mat_permu = list_permu.reshape(mat.shape)

            mat_permu = pd.DataFrame(mat_permu)
            mat_permu.index = list_gene
            mat_permu.columns = list_cell
            mat_permu.to_csv(os.path.join(path_data, '_'.join([name_mat, 'expPermuted', str(rate), str(rand_i)])))

