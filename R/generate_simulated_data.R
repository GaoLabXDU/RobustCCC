
#' Generate simulated data
#'
#' @param list_simulated_type List of simulated type
#' @param path_data Path of data
#' @param name_mat_CTP Name of scRNA-seq data for cell type pairs
#' @param mat_ori Original expression matrix
#' @param name_label_CTP Name of cell label for cell type pairs
#' @param label_ori Original cell label
#' @param Lcell Name of ligand(sender) cell
#' @param Rcell Name of receptor(receiver) cell
#'
#' @return
#' @export
#'
#' @examples
run_generate_simulated_data <- function(list_simulated_type, path_data, name_mat_CTP, mat_ori, name_label_CTP, label_ori, Lcell, Rcell){
  gsd <- import_from_path('generate_simulated_data_py', path=system.file('python', package='RobustCCC'), convert=TRUE)
  num_of_rep = 3
  for(simulated_type in list_simulated_type){
    if(simulated_type=="simuReplicate"){list_rate = c(95,90,85)}else{list_rate = c(5,10,15)}
    gsd$generate_simulated_data(simulated_type, path_data, list_rate, num_of_rep, name_mat_CTP, mat_ori, name_label_CTP, label_ori, Lcell, Rcell)
  }
}
