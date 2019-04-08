# library(dplyr)

#' generate GPR dataframe
#' @param model sybil model
#' @return dataframe including association between reactions, genes, and GPR
#' @seealso 
#' @export
#' @examples
#' 
generate_gpr <- function(model){
  GPR <- c()
  GPR$react <- model@react_id
  GPR$GENE <- model@genes
  GPR$GPR <- model@gpr

  return(GPR)
}

#' identify all genes associated with a reaction
#' @param GPR GPR dataframe
#' @param rxn_id string indicating the reaction of interest
#' @return list of strings of genes associated with the reaction of interest
#' @seealso 
#' @export
#' @examples
#' 
get_genes_from_react_id <- function(GPR, rxn_id){
  idx <- which(GPR$react == rxn_id)
  if (length(idx) == 0){
    return()
  }
  # print(paste(idx, GPR$GENE[idx]))
  return(GPR$GENE[idx][[1]])
}

#' extrapolate a gene set from a reaction set based on GPR
#' @param GPR GPR dataframe
#' @param rxn_set_list list of reaction sets
#' @return list of gene sets
#' @seealso 
#' @export
#' @examples
#' 
gene_set_from_rxn_set <- function(GPR, rxn_set_list){
  gene_set_list <- c()

  for (i in 1:length(rxn_set_list)){
    #print(i)
    gene_set <- c()
    for (rxn in rxn_set_list[[i]]){
      gene_set <- union(gene_set, get_genes_from_react_id(GPR,rxn))
    }
    if (length(gene_set) > 0){
      gene_set_list[i] <- list(gene_set)
    }
    else{
      print('empty gene set')
    }
  }

  return(gene_set_list)
}
