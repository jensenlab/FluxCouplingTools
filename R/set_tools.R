# functions related to identification and analysis of reaction/gene sets

#' get position of the set containing an element in a list of sets. may return more than one position
#' @param rxn element to search for
#' @param sets list of vectors to search through
#' @return list of sets containing rxn
#' @seealso 
#' @export
#' @examples
#' 
get_set_idx <- function(rxn, sets){
  
  idx <- grep(rxn, sets)
  sets_idxs <- c()
  for (j in idx){
    if (rxn %in% sets[[j]]){
      sets_idxs <- c(sets_idxs, j)
    }
  }
  return(sets_idxs)
}

#' get the position of the set containing all members of a group of elements
#' @param rxns list of elements to search for
#' @param sets list of vectors to search through
#' @return integer index of set containing rxns
#' @seealso 
#' @export
#' @examples
#' 
get_containing_set_idx <- function(rxns, sets){
  idxs <- get_set_idx(rxns[1], sets)
  
  for (i in idxs){
    if (all(rxns %in% sets[[i]])){
      return(i)
    }
  }
  return(0)
}

map_elements_to_set <- function(elem_list, set_list){
  elem_map <- matrix(data = NA, nrow = 1, ncol = length(elem_list))

  for (i in 1:length(elem_list)){
    idx <- get_set_idx(elem_list[i], set_list)
    # print(elem_list[i])
    # print(any(grepl(elem_list[i], mutans@allGenes)))
    # print(idx)
    if (length(idx) > 0){
      elem_map[i] <- idx
    }
  }

  return(elem_map)
}

#' check to see if set list contains elements of interest
#' @param rxns list if elements to check for
#' @param sets list of vectors to search through
#' @return boolean value indicating whether or not all elementst are present
#' @seealso 
#' @export
#' @examples
#' 
check_sets_for_containing <- function(rxns, sets){

  for (i in 1:length(sets)){
    if (all(rxns %in% sets[[i]])){
      return(TRUE)
    }
  }

  return(FALSE)
}

check_all_sets_for_containing <- function(comparison_sets, target_sets){
  set_error <- c()
  for (set_idx in 1:length(comparison_sets)){
    if (!check_sets_for_containing(comparison_sets[[set_idx]], target_sets)){
      # print(set)
      set_error <- c(set_error, set_idx)
      # return(FALSE)
    }
  }
  return(set_error)
  # return(TRUE)
}

check_sets_for_deletion <- function(rxn_list, sets){
  deleted <- c()
  for (i in rxn_list){
    if (!check_sets_for_containing(i, sets)){
      deleted <- c(deleted, i)
    }
  }

  return(deleted)
}

check_set_list_for_deletion <- function(rxn_list, set_list){
  deletions <- c()
  for (i in 1:length(set_list)){
    deletions[[i]] <- check_sets_for_deletion(rxn_list, set_list[[i]])
  }

  return(deletions)
}

find_all_sets_for_rxn <- function(rxn_id, set_lists){
  sets <- c()

  for (i in 1:length(set_lists)){
    set <- set_lists[[i]][get_set_idx(rxn_id, set_lists[[i]])]

    if (length(set) > 0){
      sets[i] <- set
    }
  }

  return(sets)
}

#' return T/F based on whether or not two set-lists are equivalent
#' @param set_1 first set 
#' @param set_2 second set
#' @return boolean indicating whether or not sets are identical
#' @seealso 
#' @export
#' @examples
#' 
compare_sets <- function(set_1, set_2){
  if (length(set_1) != length(set_2)){
    return(FALSE)
  }

  for (i in 1:length(set_1)){
    if (!all.equal(set_1[[i]], set_2[[i]])){
      return(FALSE)
    }
  }

  return(TRUE)
}

#' find the sets contained by 'rxns'
#' @param rxns list of elements in a set of interest
#' @param sets list of vectors indicating sets
#' @return list of integers indicating the indexes of the sets which are contained in the list of elements
#' @seealso 
#' @export
#' @examples
#' 
find_composing_sets <- function(rxns, sets){
  composition <- c()
  for (i in 1:length(sets)){
    if (length(sets[[i]]) == 0){next}
    if (any(sets[[i]] %in% rxns)){
      composition <- c(composition, i)
    }
  }

  return(composition)
}

optimize_suppression_idxs <- function(model, sets){
  idxs <- c()
  for (i in sets){
    idxs <- c(idxs, GRB_get_rxn_idx(model, i[[1]]))
  }

  return(unique(idxs))
}

#' generate sets from coupling mtx
#' @param coupling_mtx boolean matrix in which TRUE values indicate coupling between reaction or gene on rows and columns
#' @param init_mtx boolean matrix of same size and elements as coupling_mtx indicating additional couplings to add; defaults to NULL
#' @return list of sets
#' @seealso 
#' @export
#' @examples
#' 
get_list_of_sets_from_mtx <- function(coupling_mtx, init_mtx = NULL){ #2d columns

  if (!all.equal(rownames(coupling_mtx), colnames(coupling_mtx))){
    print("names mismatch in coupling matrix")
    break
  }

  if (!is.null(init_mtx)){
    if (!all.equal(rownames(coupling_mtx), rownames(init_mtx)) & !all.equal(colnames(coupling_mtx), colnames(init_mtx))){
      coupling_mtx[which(init_mtx)] <- TRUE
    }
  }
  coupling_mtx <- fill_coupling_matrix(coupling_mtx)
  sets <- c()

  rxns_list <- rownames(coupling_mtx)
  active <- matrix(data = TRUE, nrow = length(rxns_list), ncol = 1)

  for (i in 1:length(rxns_list)){
    if (!active[i]){next}
    if (!coupling_mtx[i,i]){active[i] <- FALSE; next}

    rxns_in_set <- which(coupling_mtx[i,])
    active[rxns_in_set] <- FALSE
    sets <- c(sets, list(rxns_list[rxns_in_set]))
  }

  return(sets)
}

core_rxn_id <- function(rxn_id){ # rxn id with parenthesis
  return(strsplit(rxn_id, split = "\\(")[[1]][1])
}

optimize_set_elements <- function(sets){
  elems <- unique(unlist(sets))
  n = length(elems)
  vec <- matrix(data = TRUE, nrow = 1, ncol = n)
  names(vec) <- elems
  for (i in sets){
    vec[i] <- FALSE
    vec[i[1]] <- TRUE
  }
  
  return(names(vec[which(vec)]))
}
