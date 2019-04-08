# library(pryr)
# library(data.tree)

#' Generate FALCON model from base metabolic Sybil model
#' @param model original Sybil model for organism metabolism
#' @param gene_sets known gene sets; defaults to empty
#' @param rxn_sets known rxn sets; defaults to empty
#' @return Sybil model expanded to include enzymatic activity
#' @seealso 
#' @export
#' @examples
generate_falcon_model <- function(model, gene_sets = c(), rxn_sets = c()){
  
  og_dim <- dim(model@S)
  og_react_id <- model@react_id # model$get_names()$VarName
  og_met_id <- model@met_id
  og_genes <- model@genes # may need to pass this in
  og_allGenes <- model@allGenes # may need to pass this in
  rxnGeneMat <- model@rxnGeneMat # may not need
  colnames(rxnGeneMat) <- model@allGenes # may not need
  rownames(rxnGeneMat) <- model@react_id # may not need
  
  marked_genes <- matrix(FALSE, nrow = length(og_allGenes), ncol = 1)
  marked_rxns <- matrix(FALSE, nrow = length(og_react_id), ncol = 1)
  rownames(marked_genes) <- model@allGenes
  rownames(marked_rxns) <- model@react_id
  
  # add all exchange reactions
  for (gene in og_allGenes){
    new_met <- paste('a', gene, sep = '_')
    model <- addExchReact(model, met = new_met, lb = -1000, ub = 1000)
  }
  
  genes_from_path <- function(path, rxn_idx){
    gene_idxs <- gprRule_to_idx(path)
    genes <- og_genes[rxn_idx]
    return(genes[[1]][gene_idxs])
  }
  
  # helper function for simple_add and or_add
  normal_add <- function(model, new_met_list, rxn_id, simple = FALSE, addExch = FALSE, identifier = NULL){
    
    rxn_idx <- which(model@react_id == rxn_id)
    exch <- findExchReact(model)
    
    # metabolites of existing reaction
    old_met_idxs <- which(model@S[(1:og_dim[1]), rxn_idx] != 0)
    old_met_list <- og_met_id[old_met_idxs]
    old_met_coeff <- model@S[old_met_idxs, rxn_idx]
    
    met_list <- c(unlist(old_met_list), unlist(new_met_list))
    met_list <- met_list[!is.na(met_list)]
    
    if (!simple){ # add reverse reaction if needed # change from 0 -> 1000 & -1000 -> 0
      lowbnd <- model@lowbnd[rxn_idx]
      uppbnd <- model@uppbnd[rxn_idx]
      model <- addReact(model, paste(rxn_id, identifier, 'fwd', sep = '_'), met = met_list,
                        Scoef = c(old_met_coeff, rep(-1, length(new_met_list))), lb = 0, ub = uppbnd, reversible = FALSE)
      model <- addReact(model, paste(rxn_id, identifier, 'rev', sep = '_'), met = met_list,
                        Scoef = c(old_met_coeff, rep(1, length(new_met_list))), lb = lowbnd, ub = 0, reversible = FALSE)
      
      ## NEED CONSTRAINTS TO PREVENT MODEL FROM PUSHING FLUX THROUGH BOTH DIRECTIONS AT ONCE
    }
    else { # changed from [-1000 to 1000] to [lowbnd to uppbnd]
      lowbnd <- model@lowbnd[rxn_idx]
      uppbnd <- model@uppbnd[rxn_idx]
      model <- addReact(model, rxn_id, met = met_list,
                        Scoef = c(old_met_coeff, rep(-1, length(new_met_list))), lb = lowbnd, ub = uppbnd, reversible = TRUE)
    }
    
    # rxn_removal_ids <- c(rxn_removal_ids, rxn_id)
    return(model)
  }
  
  # used for direct addition of gene/enzyme to reaction; not used when reaction_activity is needed (or case)
  simple_add <- function(model, new_met_list, rxn_id, simple = TRUE){
    # print(c(rxn_id, ":", new_met_list))
    model <- normal_add(model, new_met_list, rxn_id, simple = TRUE)
    return(model)
  }
  
  # used for multiple gene combinations; reaction_activity needed
  or_add <- function(model, path_list, rxn_id, simple = TRUE, split = TRUE){
    rxn_idx <- which(og_react_id == rxn_id)
    rxn_activity <- paste('a', rxn_id, sep = "_")
    
    # add conversion for path to react_activity
    identifier <- 1 # need to differentiate breakdown of rxn in name
    for (mets in path_list){ # CHECK THIS FUNCTION IF EVERYTHING BREAKS
      genes <- genes_from_path(mets, rxn_idx)
      new_mets <- paste('a_', genes, sep = '')
      met_list <- c(unlist(new_mets), rxn_activity) # genes first, then reaction_activity
      met_list <- met_list[!is.na(met_list)] # need to get rid of this
      coeff_list <- c(unlist(rep(-1, length(new_mets))), 1)
      if (split){ # genes are always consumed, reaction_activity is consumed or produced; reactions proceed forwards
        model <- addReact(model, paste(rxn_activity, 'fwd_conversion', identifier, sep = '_'), met = met_list,
                          Scoef = c(unlist(rep(-1, length(new_mets))), 1), lb = 0, ub = 1000, reversible = FALSE)
        model <- addReact(model, paste(rxn_activity, 'rev_conversion', identifier, sep = '_'), met = met_list,
                          Scoef = c(unlist(rep(-1, length(new_mets))), -1), lb = 0, ub = 1000, reversible = FALSE)
      }
      else {
        model <- addReact(model, paste(rxn_activity, 'conversion', identifier, sep = '_'), met = met_list,
                          Scoef = coeff_list, lb = -1000, ub = 1000, reversible = TRUE)
      }
      identifier <- identifier + 1
    }
    
    # add activity specific to react
    model <- normal_add(model, new_met_list = c(rxn_activity), rxn_id, simple = simple, addExch = FALSE)
    
    return(model)
  }
  
  # list of all rxns w which a gene participates
  gene_rxn_recurrence <- c()
  for (i in 1:length(model@allGenes)){
    gene_rxn_recurrence[i] <- list(which(rxnGeneMat[,i] == TRUE))
  }
  names(gene_rxn_recurrence) <- model@allGenes
  gene_rxn_promiscuity <- sapply(gene_rxn_recurrence, function(x) length(x))
  
  # loyal genes and reactions
  
  #print('RXN EXCLUSIVITY')
  
  rxn_exclusive_genes <- which(gene_rxn_promiscuity == 1)
  loyal_rxns <- c()
  for (i in gene_rxn_recurrence[rxn_exclusive_genes]){
    if (all(sapply(model@genes[i], function(x) x %in% names(rxn_exclusive_genes)))){
      loyal_rxns <- c(loyal_rxns, i)
    }
  }
  loyal_rxn_idxs <- unique(loyal_rxns)
  
  # all reactions which are comprised of non-promiscuous genes
  for (i in loyal_rxn_idxs){
    rxn_id <- og_react_id[i]
    gpr_rule <- model@gprRules[i]
    gpr_paths <- find_gpr_paths(gpr_rule)
    if (length(gpr_paths) == 1){
      genes <- genes_from_path(gpr_paths[[1]], i)
      model <- simple_add(model, new_met_list = paste('a_', genes, sep = ''), rxn_id)
      marked_genes[which(model@allGenes %in% genes)] <- TRUE
      marked_rxns[i] <- TRUE
    }
    else {
      genes <- c()
      for (path in gpr_paths){
        genes <- c(genes, genes_from_path(path, i))
      }
      genes <- unlist(genes)
      if (genes[1] == ''){next}
      model <- or_add(model, gpr_paths, og_react_id[i], split = TRUE)
      marked_genes[which(model@allGenes %in% genes)] <- TRUE
      marked_rxns[i] <- TRUE
    }
  }
  
  remaining_genes <- model@allGenes[-c(which(marked_genes == TRUE))]
  remaining_rxns <- og_react_id[-c(which(marked_rxns == TRUE))]
  
  # NORMAL ADD FOR ALL REMAINING REACTIONS
  for (rxn_idx in which(marked_rxns == FALSE)){
    rxn_id <- og_react_id[rxn_idx]
    gpr_rule <- model@gprRules[rxn_idx]
    gpr_paths <- find_gpr_paths(gpr_rule)
    if (nchar(gpr_paths[1]) == 0){next}
    model <- or_add(model, gpr_paths, rxn_id, split = TRUE)
    marked_rxns[rxn_idx] <- TRUE
  }
  
  return(model)
}

gprRule_to_idx <- function(list){
  for (i in 1:length(list)){
    list[i] <- gsub("[^0-9]", "", list[i])
  }

  return(as.numeric(list))
}

#' process names of coupled reactions to replace the names of enzymatic activity reactions to just indicate the genes
#' @param set_list list of lists containing sets of coupled reactions
#' @return processed list of lists
#' @seealso 
#' @export
#' @examples
clean_rxn_names_in_set <- function(set_list){

  clean_ex_a <- function(name){
    new_name <- strsplit(name, 'Ex_a_')[[1]]
    return(new_name[2])
  }
  
  new_sets <- c()

  for (i in 1:length(set_list)){
    k <- 1
    new_set <- c()
    for (j in 1:length(set_list[[i]])){
      new_name <- clean_ex_a(set_list[[i]][j])
      if (!is.na(new_name)){
        # set_list[[i]][j] <- clean_ex_a(set_list[[i]][j])
        new_set[k] <- clean_ex_a(set_list[[i]][j])
        k <- k+1
      }
      else {
        # set_list[[i]][j] <- ""
        next
      }
    }
    if (length(new_set) > 0){
      new_sets[i] <- list(new_set)
    }
  }

  return(new_sets)
}

#' reduce coupling matrix to only contain the reactions for enzymatic activity
#' @param coupling_matrix Matrix of TRUE/FALSE values
#' @return Matrix of TRUE/FALSE values
#' @seealso 
#' @export
#' @examples
isolate_gene_matrix <- function(coupling_matrix){
  row_genes <- which(grepl('Ex_a', rownames(coupling_matrix)))
  col_genes <- which(grepl('Ex_a', colnames(coupling_matrix)))

  gene_matrix <- coupling_matrix[row_genes, col_genes]

  return(gene_matrix)
}

#' reduce coupling matrix to only the metabolic reactions
#' @param coupling_matrix Matrix of TRUE/FALSE values
#' @return Matrix of TRUE/FALSE values
#' @seealso 
#' @export
#' @examples
isolate_rxn_matrix <- function(coupling_matrix){
  row_rxns <- which(!grepl('a_', rownames(coupling_matrix)))
  col_rxns <- which(!grepl('a_', colnames(coupling_matrix)))

  rxn_matrix <- coupling_matrix[row_rxns, col_rxns]

  return(rxn_matrix)
}

str_to_char <- function(string){
  chars <- strsplit(string, "")[[1]]
  remove <- which(chars == " ")
  if (length(remove) > 0){
    chars <- chars[-remove]
  }
  return(chars)
}

char_to_str <- function(chars){
  string <- paste(chars,sep = '', collapse = '')
  return(string)
}

extract_logical_from_base <- function(chars){
  and <- which(chars == "&")
  or <- which(chars == "|")
  
  if ((length(and) > 0) & (length(or) > 0)){
    print('not base')
    return()
  }
  
  if ((length(and) > 0)){
    return('&')
  }
  
  if ((length(or) > 0)){
    return('|')
  }
}

find_gpr_paths <- function(gprRule){
  
  gpr_char <- str_to_char(gprRule)
  len_logical <- which(gpr_char %in% c('&', '|'))
  
  if (length(len_logical) == 0){
    return(gsub("[()]", "", gprRule))
  }
  
  ast <- capture.output(call_tree(parse(text = gprRule)))
  
  logic_flag <- FALSE
  i <- 0 #starting idx
  while (!logic_flag){
    i <- i + 1
    if (('&' %in% str_to_char(ast[i])) | ('|' %in% str_to_char(ast[i]))){
      logic_flag <- TRUE
    }
  }
  
  gpr_tree <- tree_building(ast, i)
  paths <- get_paths_from_gpr(gpr_tree)
  
  return(paths)
}

tree_building <- function(list, index){
  
  chars <- strsplit(list[index], "")[[1]]
  og_depth <- (which(strsplit(list[index], "")[[1]] == "\\")-1)/2
  og_id <- chars[length(chars)]
  id1 <- "_" #og_id
  id2 <- "_"
  
  new_node <- Node$new(og_id)
  
  start_1 <- index + 1
  chars <- strsplit(list[start_1], "")[[1]]
  depth1 <- (which(chars == "\\")-1)/2
  id1 <- chars[length(chars)]
  
  start_2 <- start_1 + 1
  chars <- strsplit(list[start_2], "")[[1]]
  depth2 <- (which(chars == "\\")-1)/2
  id2 <- chars[length(chars)]
  while (depth2 != depth1 | id2 != id1){
    start_2 <- start_2 + 1
    chars <- strsplit(list[start_2], "")[[1]]
    depth2 <- (which(chars == "\\")-1)/2
    id2 <- chars[length(chars)]
  }
  index1 <- start_1

  while (!(id1 %in% c( '|', '&', 'x'))){
    index1 <- index1 + 1
    chars <- strsplit(list[index1], "")[[1]]
    depth1 <- (which(chars == "\\")-1)/2
    id1 <- chars[length(chars)]
  }
  
  index2 <- start_2
  while (!(id2 %in% c( '|', '&', 'x'))){
    index2 <- index2 + 1
    chars <- strsplit(list[index2], "")[[1]]
    depth2 <- (which(chars == "\\")-1)/2
    id2 <- chars[length(chars)]
  }
  
  left_node <- Node$new("_")
  right_node <- Node$new("_")
  if (id1 == 'x'){
    chars <- strsplit(list[index1+1], "")[[1]]
    
    len <- length(chars)
    id1 <- chars[len]
    while (id1 == " "){
      len <- len - 1
      id1 <- chars[len]
    }
    if (chars[len-1] != " "){
      id1 <- paste(chars[len-1], chars[len], sep = '')
    }
    
    left_node <- Node$new(id1)
  }
  else{
    left_node <- tree_building(list, index1)
  }
  if (id2 == 'x'){
    chars <- strsplit(list[index2+1], "")[[1]]
    
    len <- length(chars)
    id2 <- chars[len]
    while (id2 == " "){
      len <- len - 1
      id2 <- chars[len]
    }
    if (chars[len-1] != " "){
      id2 <- paste(chars[len-1], chars[len], sep = '')
    }
    
    right_node <- Node$new(id2)
  }
  else {
    right_node <- tree_building(list, index2)
  }
  
  if (right_node$name == left_node$name){
    right_node$name <- paste(right_node$name, 2, sep = '')
  }
  
  new_node$AddChildNode(left_node)
  new_node$AddChildNode(right_node)
  return(new_node)
}

get_paths_from_gpr <- function(node){
  paths <- c()
  name <- node$name
  # print(name)
  name <- strsplit(name, '')[[1]][1]
  # if (length(name) == 0){
  #   return(c())
  # }
  
  if (name != '&' & name != '|'){
    paths <- c(paste('x', '[', node$name, ']', sep = ''))
    return(paths)
  }
  
  left_paths <- get_paths_from_gpr(node$children[[1]])
  right_paths <- get_paths_from_gpr(node$children[[2]])
  
  if (name == '&'){
    idx <- 0
    for (i in 1:length(left_paths)){
      for (j in 1:length(right_paths)){
        idx <- idx + 1
        paths[idx] <- list(c(unlist(left_paths[i]), unlist(right_paths[j])))
      }
    }
  }
  
  if (name == '|'){
    for (i in 1:length(left_paths)){
      paths[i] <- list(unlist(left_paths[i]))
    }
    
    for (i in 1:length(right_paths)){
      paths[length(left_paths) + i] <- list(unlist(right_paths[i]))
    }
  }
  
  return(paths)
}

get_all_gpr_paths <- function(gprRules){
  paths <- c()
  
  for (i in gprRules){
    paths <- c(paths, list(find_gpr_paths(i)))
  }
  
  return(paths)
}
