## MODELS

#' #' description
#' #' @return GRB model
#' #' @seealso 
#' #' @export
#' #' @examples
#' #' 
#' GRB_ecoli_model <- function(){
#'   data(Ec_core);
#'   model=Ec_core;
#'   model <- changeBounds(model, 11, lb = 0)
#'   # model <- changeBounds(model, 13, lb = 0, ub = 0)
#'   model <- rmReact(model = model, react = 13)
#'   for (i in findExchReact(model)@react_pos){
#'     model <- changeBounds(model, i, lb = -1000, ub = 1000)
#'     # if (model@lowbnd[i] == 0){
#'     #   model <- changeBounds(model, i, lb = -1000)
#'     # }
#'   }
#' 
#'   ecoli <- as_GRBmodel(model)
#'   ecoli$show_output(FALSE)
#'   return(ecoli)
#' }
#' 
#' #' description
#' #' @return GRB falcon model
#' #' @seealso 
#' #' @export
#' #' @examples
#' #' 
#' GRB_ecoli_falcon_model <- function(){
#'   sybil_ecoli <- get_ecoli_model()
#'   ecoli_falcon_model <- GRB_generate_falcon_model(sybil_ecoli)
#'   return(ecoli_falcon_model)
#' }
#' 
#' #' Gurobi S. cerevisiae model
#' #' @param model filename of RData file for sybil model; defaults to 'yeast_model.RData'
#' #' @return GRB model
#' #' @seealso 
#' #' @export
#' #' @examples
#' #' 
#' GRB_yeast_model <- function(model = 'yeast_model.RData'){
#'   #maranas_model_exch_add_biom_rm/maranas_model_lipid_exch
#'   load(model)
#' 
#'   # setwd("~/GitHub/PathwayMining/")
#' 
#'   yeast <- as_GRBmodel(yeast_model)
#'   yeast$show_output(FALSE)
#' 
#'   return(yeast)
#' }
#' 
#' #' Gurobi S. cerevisiae falcon model
#' #' @param model filename of RData file for sybil model; defaults to 'yeast_model.RData'
#' #' @return GRB falcon model
#' #' @seealso 
#' #' @export
#' #' @examples
#' #' 
#' GRB_yeast_falcon_model <- function(model = 'yeast_model.RData'){
#' 
#'   load(model)
#' 
#'   sybil_yeast <- yeast_model
#'   yeast_falcon_model <- GRB_generate_falcon_model(sybil_yeast)
#'   return(yeast_falcon_model)
#' }
#' 
#' #' Gurobi S. mutans model
#' #' @param model filename of RData file for sybil model; defaults to 'mutans_model.RData'
#' #' @return GRB model
#' #' @seealso 
#' #' @export
#' #' @examples
#' #' 
#' GRB_mutans_model <- function(model = 'mutans_model.RData'){
#' 
#'   load(model)
#'  
#'   mutans <- as_GRBmodel(mutans)
#'   mutans$show_output(FALSE)
#' 
#'   return(mutans)
#' }
#' 
#' #' Gurobi S. mutans falcon model
#' #' @param model filename of RData file for sybil model; defaults to 'mutans_model.RData'
#' #' @return GRB falcon model
#' #' @seealso 
#' #' @export
#' #' @examples
#' #' 
#' GRB_mutans_falcon_model <- function(model = 'mutans_model.RData'){
#'   load(model)
#' 
#'   sybil_mutans <- mutans
#'   mutans_falcon_model <- GRB_generate_falcon_model(sybil_mutans)
#'   return(mutans_falcon_model)
#' }
#' 
#' #' Gurobi P. aeruginosa model
#' #' @param model filename of RData file for sybil model; defaults to 'pao_model.RData'
#' #' @return GRB model
#' #' @seealso 
#' #' @export
#' #' @examples
#' #' 
#' GRB_pao_model <- function(model = 'pao_model.RData'){
#'   load(model)
#' 
#'   pao <- as_GRBmodel(pao_model)
#'   pao$show_output(FALSE)
#' 
#'   return(pao)
#' }
#' 
#' #' Gurobi P. aeruginosa falcon model
#' #' @param model filename of RData file for sybil model; defaults to 'pao_model.RData'
#' #' @return GRB falcon model
#' #' @seealso 
#' #' @export
#' #' @examples
#' #' 
#' GRB_pao_falcon_model <- function(model = 'pao_model.RData'){
#'   load(model)
#' 
#'   pao_falcon <- GRB_generate_falcon_model(pao_model)
#'   pao_falcon$show_output(FALSE)
#' 
#'   return(pao_falcon)
#' }

## FUNCTIONS

#' Convert the output of GRB_generate_set_lists_cluster to a coupling matrix. G sets or R sets are fully coupled to each other in G* or R* sets
#' @param coupling_list list of list of integers which each indicate the couplings induced by each suppression in PACT
#' @param n_react number of reactions
#' @param vars names of the reactions; length should equal n_react
#' @param init_sets known sets to indicate coupling for in the matrix
#' @return boolean matrix indicating coupling between reactions
#' @seealso 
#' @export
#' @examples
#' 
full_coupling_matrix_from_coupling_vector_list <- function(coupling_list, n_react, vars, init_sets = NULL){

  coupling_matrix <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
  rownames(coupling_matrix) <- vars
  colnames(coupling_matrix) <- vars
  if (!is.null(init_sets)){
    coupling_matrix <- fill_coupling_matrix_from_sets(coupling_matrix, init_sets)
  }
  
  for (i in 1:length(coupling_list)){
    if (is.null(coupling_list[[i]])){next}
    coupling_matrix[coupling_list[[i]]] <- TRUE
  }
  
  coupling_matrix <- fill_coupling_matrix(coupling_matrix)
  
  return(coupling_matrix)
}

#' Convert the output of GRB_generate_set_lists_cluster to a coupling matrix
#' @param coupling_list list of list of integers which each indicate the couplings induced by each suppression in PACT
#' @param n_react number of reactions
#' @param vars names of the reactions; length should equal n_react
#' @param init_sets known sets to indicate coupling for in the matrix
#' @return boolean matrix indicating coupling between reactions
#' @seealso 
#' @export
#' @examples
#' 
coupling_matrix_from_coupling_vector_list <- function(coupling_list, n_react, vars, init_sets){
  
  coupling_matrix <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
  rownames(coupling_matrix) <- vars
  colnames(coupling_matrix) <- vars
  for (i in 1:length(coupling_list)){
    if (is.null(coupling_list[[i]])){next}
    print(i)
    intermediate_mtx <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
    rownames(intermediate_mtx) <- vars
    colnames(intermediate_mtx) <- vars
    intermediate_mtx[coupling_list[[i]]] <- TRUE
    intermediate_mtx <- fill_coupling_matrix_from_sets(intermediate_mtx, init_sets)
    coupling_matrix[which(intermediate_mtx)] <- TRUE
  }
  
  return(coupling_matrix)
}

#' identify which G sets and R sets within the same G* or R* set never couple to each other
#' @param full_coupling_mtx boolean matrix indicating fully coupled G* or R* sets
#' @param coupling_mtx  boolean matrix indicating only observed couplings
#' @param n_react number of reactions
#' @return boolean matrix with TRUE values indicating which reactions in the same G* or R* set never couple
#' @seealso 
#' @export
#' @examples
#' 
identify_intermediate_uncoupled <- function(full_coupling_mtx, coupling_mtx, n_react){
  uncoupled_matrix <- Matrix(data = FALSE, nrow = n_react, ncol = n_react)
  rownames(uncoupled_matrix) <- rownames(coupling_mtx)
  colnames(uncoupled_matrix) <- colnames(coupling_mtx)
  
  uncoupled <- which(full_coupling_mtx & !coupling_mtx)
  uncoupled_matrix[uncoupled] <- TRUE
  
  return(uncoupled_matrix)
}
