#' Check CCC packages
#'
#' Check R and python packages of cell-cell communication methods
#'
#' @param list_methods List of cell-cell communication methods
#'
#' @return
#' @export
#'
#' @examples
check_CCC_packages <- function(list_methods){
  # all supported methods
  list_methods_all_R = c("cellcall","CellChat","CytoTalk","icellnet","iTALK","nichenetr","scMLnet","SingleCellSignalR")
  list_methods_all_python = c("scConnect","CellPhoneDB")

  # user selected methods
  list_methods_select_R = intersect(list_methods, list_methods_all_R)
  list_methods_select_python = intersect(list_methods, list_methods_all_python)

  # stats not-installed methods
  list_methods_not_installed_R = setdiff(list_methods_select_R, rownames(installed.packages()))
  list_methods_not_installed_python = setdiff(tolower(list_methods_select_python), reticulate::py_list_packages()$package)

  list_methods_not_installed = c(list_methods_not_installed_R, list_methods_not_installed_python)

  if(length(list_methods_not_installed)>1){
    cat(paste(paste(list_methods_not_installed[1:length(list_methods_not_installed)-1],collapse = ', '),
              'and',list_methods_not_installed[length(list_methods_not_installed)] ,'is/are not installed\n',sep=' '))
  }
  else if (length(list_methods_not_installed)==1) {
  cat(paste(paste(list_methods_not_installed[1:length(list_methods_not_installed)-1],collapse = ', '),
              list_methods_not_installed[length(list_methods_not_installed)] ,'is/are not installed\n',sep=' '))
  }
  else{cat('All selected methods are installed\n')}
  return(list_methods_not_installed)
}


#' Install CCC packages
#'
#' Install R and python packages of cell-cell communication methods
#'
#' @param List of cell-cell communication methods
#'
#' @return
#' @export
#'
#' @examples
install_CCC_packages <- function(list_methods){
  # all supported methods
  list_methods_all_R = c("cellcall","CellChat","CytoTalk","icellnet","iTALK","nichenetr","scMLnet","SingleCellSignalR")
  list_methods_all_python = c("scConnect","cellphonedb")

  # user selected but not-installed methods
  list_methods_select_R = intersect(list_methods, list_methods_all_R)
  list_methods_select_python = intersect(tolower(list_methods), tolower(list_methods_all_python))

  # install methods
  for(method in list_methods_select_R){install_r_packages(method)}
  for(method in list_methods_select_python){install_python_packages(method)}

}

#' @noRd
install_r_packages <- function(method){
  # install github R packages

  # cellcall
  # https://github.com/ShellyCoder/cellcall
  if(method == "cellcall"){
    cat('start install cellcall\n')
    BiocManager::install(c('clusterProfiler', 'ComplexHeatmap', 'enrichplot', 'DOSE'))
    devtools::install_github("ShellyCoder/cellcall",type = "source")
  }
  # CellChat
  # https://github.com/sqjin/CellChat
  if(method == "CellChat"){
    cat('start install CellChat\n')
    BiocManager::install(c('ComplexHeatmap', 'BiocNeighbors'))
    devtools::install_github("sqjin/CellChat")
  }

  # CytoTalk
  # https://github.com/tanlabcode/CytoTalk
  if(method == "CytoTalk"){
    cat('start install CytoTalk\n')
    devtools::install_github("tanlabcode/CytoTalk@c00132c")
  }

  # iCellNet
  # https://github.com/soumelis-lab/ICELLNET
  if(method == "icellnet"){
    cat('start install icellnet\n')
    BiocManager::install(c('AnnotationDbi', 'hgu133plus2.db'))
    devtools::install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")
  }

  # iTALK
  # https://github.com/Coolgenome/iTALK
  if(method == "iTALK"){
    cat('start install iTALK\n')
    devtools::install_github("Coolgenome/iTALK", build_vignettes = FALSE, force = TRUE)
  }

  # NicheNet
  # https://github.com/saeyslab/nichenetr
  if(method == "nichenetr"){
    cat('start install nichenetr\n')
    BiocManager::install(c('ComplexHeatmap'))
    devtools::install_github("saeyslab/nichenetr")
    if(!file.exists(paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'ligand_target_matrix.rds',sep='//'))){
      ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
      # ligand_target_matrix = readRDS(url("https://media.githubusercontent.com/media/chenxing-zhang/RobustCCC/main/inst/depend_cell_cell_communication_methods/NicheNet/ligand_target_matrix.rds"))
      saveRDS(ligand_target_matrix,paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'ligand_target_matrix.rds',sep='//'))
    }else if(file.size(paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'ligand_target_matrix.rds',sep='//'))<100000000){
      ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
      # ligand_target_matrix = readRDS(url("https://media.githubusercontent.com/media/chenxing-zhang/RobustCCC/main/inst/depend_cell_cell_communication_methods/NicheNet/ligand_target_matrix.rds"))
      saveRDS(ligand_target_matrix,paste(system.file('depend_cell_cell_communication_methods/NicheNet',package='RobustCCC'),'ligand_target_matrix.rds',sep='//'))
    }
  }

  # scMLnet
  # https://github.com/SunXQlab/scMLnet
  if(method == "scMLnet"){
    cat('start install scMLnet\n')
    devtools::install_github("YUZIXD/scMLnet")
  }

  # SingleCellSignalR
  # https://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html
  if(method == "SingleCellSignalR"){
    cat('start install SingleCellSignalR\n')
    BiocManager::install("SingleCellSignalR")
  }
}


#' @noRd
install_python_packages <- function(method){
  # NOTE: The scconnect method calls the DataFrame.iterrows() function, which is no longer supported in pandas (>1.5.3)
  reticulate::py_install('pandas==1.5.3', pip = TRUE)
  reticulate::py_install('fsspec', pip = TRUE)
  reticulate::py_install('rbo', pip = TRUE)

  # NATMI, no installation required, already included in 'NATMI.py'
  # https://github.com/asrhou/NATMI
  # git clone https://github.com/asrhou/NATMI.git
  # if(method == 'natmi'){
  #   reticulate::conda_install(conda_env,'fsspec')
  # }

  # scConnect
  # https://github.com/JonETJakobsson/scConnect
  if(method == 'scconnect'){
    cat('start install scConnect\n')
    # reticulate::conda_install(conda_env,'scConnect')
    reticulate::py_install('scConnect', pip = TRUE)
  }

  # CellPhoneDB
  # https://github.com/Teichlab/cellphonedb
  if(method == 'cellphonedb'){
    cat('start install cellphonedb\n')
    # reticulate::conda_install(conda_env,'cellphonedb')
    reticulate::py_install('cellphonedb', pip = TRUE)
  }
}



