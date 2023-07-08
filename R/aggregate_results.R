#' Aggregate CCC results
#'
#' merge the top inferred results of each method sorted by communication score or P-value
#'
#' @param name_mat Name of scRNA-seq data
#' @param Lcell Name of ligand(sender) cell
#' @param Rcell Name of receptor(receiver) cell
#' @param list_methods List of cell-cell communication methods
#' @param path_result Path of result
#' @param top Number of top inferred results
#'
#' @return
#' @export
#'
#' @examples
run_aggregate_results <- function(name_mat, Lcell, Rcell, list_methods, path_result, top){
  path_result_CCC <- paste(path_result,'CCCmethods',sep='//')
  ar <- import_from_path('aggregate_results_py', path=system.file('python', package='RobustCCC'), convert=TRUE)
  ar$aggregate_result(path_result_CCC, name_mat, Lcell, Rcell, list_methods, top)
}
