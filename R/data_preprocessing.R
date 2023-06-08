#' Convert gene symbol
#'
#' Convert mouse gene symbol to mouse gene symbol
#'
#' @param path_data Path of data
#' @param name_mat Name of scRNA-seq data
#' @param species Species name of scRNA-seq data
#'
#' @return
#' @export
#'
#' @examples
run_convert_gene_symbol_mouse2human <- function(path_data, name_mat, species){
  # 2.1 data pre-processing: change mouse gene to human gene
  if(species == 'mouse'){
    if(file.exists(paste(path_data,paste(name_mat, '_human', sep=""),sep='//'))){
      cat(paste(name_mat, '_human already exists', sep=""))
      return(paste(name_mat, '_human', sep=""))
    }
    path_code_depend_dp = system.file('depend_data_preprocessing',package='RobustCCC')
    dp <- import_from_path('data_preprocessing_py', path=system.file('python', package='RobustCCC'), convert=TRUE)
    dp$convert_gene_symbol_mouse2human(path_data, name_mat, path_code_depend_dp)
    name_mat = paste(name_mat, '_human', sep="")
  }
  return(name_mat)
}

#' Split original data
#'
#' Split original scRNA-seq data to cell type pair scRNA-seq data
#'
#' @param path_data Path of data
#' @param name_mat Name of scRNA-seq data
#' @param name_label Name of cell label
#' @param Lcell Name of ligand(sender) cell
#' @param Rcell Name of receptor(receiver) cell
#'
#' @return
#' @export
#'
#' @examples
run_split_mat_to_cell_type_pairs <- function(path_data, name_mat, name_label, Lcell, Rcell){
  # 2.2 data pre-processing: split whole matrix to cell-type-pairs matrix based on label
  if(file.exists(paste(path_data,paste(name_mat,Lcell, Rcell, sep="_"),sep='//'))&
     file.exists(paste(path_data,paste(name_label, Lcell, Rcell, sep="_"),sep='//'))){
    cat(paste(paste(name_mat,Lcell, Rcell, sep="_"),'and', paste(name_label,Lcell, Rcell, sep="_"),'already exists', sep=" "))
  }else{
    dp <- import_from_path('data_preprocessing_py', path=system.file('python', package='RobustCCC'), convert=TRUE)
    dp$split_mat_to_cell_type_pairs(path_data, name_mat, name_label, Lcell, Rcell)
  }
}
