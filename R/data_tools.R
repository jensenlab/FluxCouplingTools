## data tools

#' combine two tables
#' @param table_1 matrix
#' @param table_2 matrix
#' @param iter_col_1 which column of table_1 to merge
#' @param iter_col_2 which column of table_2 to merge
#' @return two-column matrix
#' @seealso 
#' @examples

combine_tables <- function(table_1, table_2, iter_col_1 = 1, iter_col_2 = 1){ # iterate down table 1, specify row
  aligned_mtx <- matrix(data = NA, nrow = nrow(table_1), ncol = ncol(table_2))

  for (i in 1:nrow(table_1)){
    it <- table_1[i, iter_col_1]
    pos <- which(table_2[, iter_col_2] == it)
    if (length(pos) > 0){
      aligned_mtx[i,] <- table_2[pos,]
    }
  }

  return(cbind(table_1, aligned_mtx[,-c(iter_col_2)]))
}

#' Read in coupling file
#' @param filename name of file which PACT output was written to
#' @param coupling_vector optional coupling vector, if warm-starting
#' @param completed_idxs optional copmpleted indexes, if warm starting
#' @return list of coupling vectors (list of lists containing positions of TRUE values in coupling matrix) and completed indexes (list of integers of suppressions already read)
#' @seealso 
#' @export
#' @examples

read_coupling_csv <- function(filename, coupling_vector = c(), completed_idxs = c()){
  lines <- readLines(filename)
  n_lines <- length(lines)

  for (line in lines){
    nums <- as.numeric(strsplit(line, split = ',')[[1]])
    idx <- nums[1]
    if (idx %in% completed_idxs){next}
    coupling_idxs <- nums[2:length(nums)]

    completed_idxs <- c(completed_idxs, idx)
    coupling_vector[[idx]] <- coupling_idxs
  }


  list(
    coupling_vector = coupling_vector,
    completed_idxs = completed_idxs
  )
}
