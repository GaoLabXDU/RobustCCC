#' Evaluate robustness
#'
#' Evaluate robustness of cell-cell communication by calculating the similarity
#' between inferred LR pairs in the original data and simulated data. The name_mat_1
#' and name_mat_2 are needed if data_type=='bioReplicate'
#'
#' @param path_result Path of result
#' @param list_simulated_type List of simulated type
#' @param list_methods List of cell-cell communication methods
#' @param name_mat Name of scRNA-seq data
#' @param Lcell Name of ligand(sender) cell
#' @param Rcell Name of receptor(receiver) cell
#' @param name_mat_1 Name of the first scRNA-seq replicate data
#' @param name_mat_2 Name of the second scRNA-seq replicate data
#'
#' @return
#' @export
#'
#' @examples
run_evaluate_robustness <- function(path_result, list_simulated_type, list_methods, name_mat, Lcell, Rcell, name_mat_1=NULL, name_mat_2=NULL){
  path_result_robust <- paste(path_result,'robustness',sep='//')
  if(!dir.exists(path_result_robust)){dir.create(path_result_robust)}
  path_result_CCC <- paste(path_result,'CCCmethods',sep='//')
  socbr <- import_from_path('evaluate_robustness_py', path=system.file('python', package='RobustCCC'), convert=TRUE)
  num_of_rep = 3
  for(simulated_type in list_simulated_type){
    if(simulated_type=="bioReplicate"){
      if(is.null(name_mat_1)|is.null(name_mat_2)){
        cat('name_mat_1 and name_mat_2 needed if data_type==bioReplicate')
        next
      }
      socbr$similarity_between_replicates(path_result_CCC,path_result_robust, name_mat_1, name_mat_2, Lcell, Rcell, list_methods)
    }else if(simulated_type=='simuReplicate'){
      list_rate = c(95,90,85)
      socbr$similarity_between_simulated_data(path_result_CCC,path_result_robust, simulated_type, name_mat, Lcell, Rcell , list_methods, list_rate, num_of_rep)
    }else{
      list_rate = c(5,10,15)
      socbr$similarity_between_simulated_and_original_data(path_result_CCC,path_result_robust, simulated_type, name_mat, Lcell, Rcell , list_methods, list_rate, num_of_rep)
    }
  }
}
