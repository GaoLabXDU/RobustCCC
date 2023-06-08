# RobustCCC
A user-friendly tool to evaluate the robustness of cell-cell communication methods
![image](https://github.com/chenxing-zhang/RobustCCC/blob/main/schematic.png)


## Descriptions
We develop a user-friendly tool, RobustCCC, which facilitates the robustness evaluation of cell-cell communication methods. RobustCCC offers the following capabilities: 1) generating simulated data, including replicated data, transcriptomic data noise and prior knowledge noise, with a single step, 2) installing and executing 14 cell-cell communication methods in a single step and 3) generating robustness evaluation reports in tabular form for easy interpretation.

## Installation

RobustCCC R package can be installed from Github using devtools:  

```
devtools::install_github('chenxing-zhang/RobustCCC')
```

RobustCCC also needs python environment using anaconda:

```
conda create -n RobustCCC_env python=3.8
```

## Tutorials
**0. activate environment and set path**
```
# 0.0 load RobustCCC and reticulate
library(RobustCCC)
library(reticulate) # use python

# 0.1 activate conda environment
path_conda_env = '//path//of//your//conda_env' 
reticulate::use_condaenv(path_conda_env) 

# 0.2 path of data, including scRNA-seq data and cell type anotation
path_data = '//path//of//your//data' 
# path_data = system.file('data',package='RobustCCC')

# 0.3 path of result. The folders named CCCmethods and similarity  will be created, which are used to save the method results and the similarity of the results between repelicate data respectively.
path_result <- '//path//of//result'
```

**1. check and install cell-cell communication methods**
```
# 1.0 select methods
list_methods = c("CellCall","CellChat","CytoTalk","iCellNet","iTALK","NicheNet","scMLnet","SingleCellSignalR","Zhou","Skelly","Kumar","NATMI","scConnect","CellPhoneDB")

# 1.1 check method avaliable
list_methods_not_installed = check_CCC_packages(list_methods) 

# 1.2 install method in list_methods
install_CCC_packages(list_methods_not_installed) 

# 1.3 check again
list_methods_not_installed = check_CCC_packages(list_methods)
```

**2. data pre-processing**
```
# 2.0 data information
name_mat = 'F001'  # 'F001' or 'F002'
name_label = 'F001_label'  # 'F001_label' or 'F002_label'
species = 'mouse'
Lcell='Endo'
Rcell='Micro'

# 2.1 data pre-processing: change mouse gene to human gene
name_mat_human = run_convert_gene_symbol_mouse2human(path_data, name_mat, species)

# 2.2 data pre-processing: split whole matrix to cell-type-pairs matrix based on label
run_split_mat_to_cell_type_pairs(path_data, name_mat_human, name_label, Lcell, Rcell)
```

**3. load mat and label**
```
# 3.0 data information, CTP: cell type pair
name_mat_CTP = paste(name_mat_human, Lcell, Rcell,sep='_')
name_label_CTP = paste(name_label, Lcell, Rcell,sep='_')

# 3.1 load mat and label
mat_ori = read.table(paste(path_data,name_mat_CTP,sep = '//'),sep=',',header = T,row.names = 1)
mat_ori = mat_ori[!duplicated(rownames(mat_ori)),]
label_ori = read.table(paste(path_data,name_label_CTP,sep='//'),sep=',',header = TRUE, row.names = 1)
```

**4. run cell-cell communication**
```
run_CCC_methods(name_mat_CTP, Lcell, Rcell, mat_ori, label_ori, list_methods, path_result)
```

**5. generate simulated data**
```
# 5.1 select simulated type
list_simulated_type = c('simuReplicate', 'GaussianNoise','dropout', 'cellTypePermu', 'ligRecPermu')

# 5.2 generate simulated data
run_generate_simulated_data(list_simulated_type, path_data, name_mat_CTP, mat_ori, name_label_CTP, label_ori, Lcell, Rcell)
```

**6. load simulated mat and label and run cell-cell communication**
```
run_CCC_methods_simulated(list_simulated_type, name_mat_CTP, Lcell, Rcell, list_methods, path_data, path_result)
```

**7. evaluate robustness**
```
# 7.1 select simulated type
list_simulated_type = c('bioReplicate', 'simuReplicate', 'GaussianNoise','dropout', 'cellTypePermu', 'ligRecPermu')

# 7.2 evaluate robustness by calculating similarity
# name_mat_1 and name_mat_2 are needed if data_type=='bioReplicate'
run_evaluate_robustness(path_result, list_simulated_type, list_methods, name_mat, name_mat_1=NULL, name_mat_2=NULL)
```
