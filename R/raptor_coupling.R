# new coupling protocol adapted from jensenlab/raptor

# SUB-ROUTINES

#' change the bounds on a Gurobi model
#' @param model Gurobi model
#' @param rxn string indicating the reaction whose bounds should be altered
#' @param UB new upper bound for reaction
#' @param LB new lower bound for reaction
#' @return new model with altered bounds
#' @seealso 
#' @export
#' @examples
#' 
set_model_bounds <- function(model, rxn, UB, LB){
  model$setattr("UB", setNames(UB, rxn))
  model$setattr("LB", setNames(LB, rxn))

  return(model)
}

#' minimize or maximize flux through a reaction in a given Gurobi model
#' @param model Gurobi model
#' @param rxn string indicating reaction to optimize flux for
#' @param max boolean value indicating whether or not to maximize flux. FALSE indicates minimization 
#' @return fluxes for all reactions in models
#' @seealso 
#' @export
#' @examples
#' 
optimize_rxn <- function(model, rxn, max){
  model$setattr("Obj", setNames(1.0, rxn))
  if (max){
    model$set_model_sense(maximize=TRUE)
  }
  else {
    model$set_model_sense(minimize=TRUE)
  }
  model$optimize()
  sol <- model$get_solution()

  return(sol)
}

# MAIN FUNCTION

#' #' determine flux coupling between reactions in a Gurobi model
#' #' @param model Gurobi model
#' #' @param min_fva minimum correlation value to test for coupling
#' #' @param fix_frac const value used in fixing flux at non-zero value 
#' #' @param fix_tol_frac maximum allowed error when determining whether flux for a reaction is fixed
#' #' @param bnd_tol allowed error in comparing max & min flux
#' #' @param stored_obs maximum number of flux instancecs to store
#' #' @param cor_iter number of iterations after which correlation is considered in checking coupling
#' #' @param reaction_indexes indexes for reactions to check; check all if empty
#' #' @param compare_mtx boolean to indicate whether or not to couple known sets together if any individual reactions are coupled; FALSE by default
#' #' @param known_set_mtx boolean matrix of known coupling between reactions
#' #' @return boolean matrix indicating whether reactions are coupled to each other (reactions are assigned to rows and columns)
#' #' @seealso 
#' #' @export
#' #' @examples
#' #' 
#' flux_coupling_raptor <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01,
#'       bnd_tol = 0.1, stored_obs = 4000, cor_iter = 5, cor_check = TRUE,
#'       reaction_indexes = c(), compare_mtx = FALSE, known_set_mtx = Matrix(data = FALSE, nrow = 1, ncol = 1, sparse = TRUE)) {
#' 
#'   # min_fva_cor is minimum correlation between fluxes
#'   # bnd_tol is allowed error in comparing max & min flux
#'   # fix_frac is const value used in fixing flux at non-zero value
#'   # fix_tol_frac is error allowed in determining whether flux is fixed
#'   # stored_obs is # not flux values to be stored
#'   # cor_iter is number of iterations after which correlation is considered in checking coupling
#' 
#'   if (!cor_check){
#'     stored_obs = 1
#'   }
#' 
#'   n <- model$get_sizes()$NumVars
#' 
#'   if (is.null(known_set_mtx)){
#'     compare_mtx <- FALSE
#'   }
#'   else {
#'     if ((nrow(known_set_mtx) < n) | (ncol(known_set_mtx) < n)){compare_mtx <- FALSE}
#'   }
#' 
#'   # if empty set, then assume all reactions are to be inspected
#'   if (length(reaction_indexes) == 0){
#'     reaction_indexes <- c(1:n)
#'   }
#' 
#'   vars <- model$get_names()$VarName
#'   prev_obj <- model$getattr("Obj")
#'   model$setattr("Obj", setNames(numeric(n), vars)) # clear the objective
#'   prev_sense <- model$getattr("ModelSense")
#' 
#'   global_max <- rep(0, n)
#'   global_min <- rep(0, n)
#' 
#'   coupled <- Matrix(FALSE, nrow=n, ncol=n, dimnames=list(vars, vars), sparse = TRUE)
#' 
#'   blocked <- rep(FALSE, n)
#'   active <- rep(TRUE, n)
#' 
#'   not_fixed <- function(x,y) { # check for variability in flux, return TRUE or FALSE
#'     !is.infinite(x) && !is.infinite(y) && abs(x - y) > fix_tol_frac*max(abs(x), abs(y))
#'   }
#' 
#'   rxn_fix <- function(max_, min_){
#' 
#'     avg <- min_ + fix_frac*(max_ - min_)
#'     ct = 0
#'     while (near(avg, 0) & ct < 5){
#'       print(c(avg, max_, min_))
#'       avg <- avg  + fix_frac*(max_ - min_)
#'       ct = ct + 1
#'     }
#' 
#'     return(avg)
#'   }
#' 
#'   correlation_check <- function(flux, i, j){ # return true if correlation is high, or no correlation --> continue comparing flux
#'     # if false, skip in depth comparison
#' 
#'     n_entries <- length(which(!near(flux[,i], 0, tol = bnd_tol) | !near(flux[,j], 0, tol = bnd_tol)))
#'     if (n_entries > cor_iter){
#'       C <- cor(flux[,i], flux[,j])
#' 
#'       if ((is.na(C)) | (abs(C) < min_fva_cor)){
#'         return(FALSE)
#'       }
#'     }
#'     return(TRUE)
#'   }
#' 
#'   flux <- matrix(c(0), nrow = stored_obs, ncol = n)
#'   lp_calls <- 0
#' 
#'   update_flux <- function(flux_, idx, sol){
#'     if (stored_obs > 0){
#'       flux_[idx,] <- sol
#'     }
#'     return(flux_)
#'   }
#' 
#'   # if j is coupled to i, couple the rxns known to be coupled to j to i as well
#'   known_set_coupling <- function(i, j, coupled, active){
#'     set <- which(known_set_mtx[j,])
#' 
#'     coupled[i, set] <- TRUE
#'     coupled[set, set] <- TRUE
#'     active[set] <- FALSE
#' 
#'     list(coupled = coupled, active = active)
#'   }
#' 
#'   for (idx in 1:(length(reaction_indexes))) { # (i in 1:(n-1))
#'     # iterate over passed in idxs instead (idx in 1:length(reaction_indexes)); i <-  reaction_indexes[idx]
#'     i <-  reaction_indexes[idx]
#' 
#'     if (!active[i] | blocked[i]) next
#'     sub_max <- rep(-Inf, n)
#'     sub_min <- rep(Inf, n)
#'     prev_ub <- model$getattr("UB")[vars[i]]
#'     prev_lb <- model$getattr("LB")[vars[i]]
#' 
#'     if (!near(global_max[i], 0, tol = bnd_tol) | !near(global_min[i], 0, tol = bnd_tol)){
#'       fixed_val <- rxn_fix(global_max[i], global_min[i])
#'     }
#'     else {
#'       if (!near(model$getattr("UB")[vars[i]], 0)){ #model$getattr("UB")[vars[i]] > tol
#' 
#'         model$setattr("Obj", setNames(1.0, vars[i]))
#'         model$set_model_sense(maximize=TRUE)
#'         model$optimize()
#'         lp_calls <- lp_calls + 1
#'         sol <- model$get_solution()
#'         global_max <- pmax(global_max, sol$X)
#'         global_min <- pmin(global_min, sol$X)
#'         flux[lp_calls%%stored_obs,] <- sol$X
#' 
#'         if (!near(global_max[i], 0) | !near(global_min[i], 0)){
#'           fixed_val <- rxn_fix(global_max[i], global_min[i])
#'         }
#' 
#'         model$setattr("Obj", setNames(0.0, vars[i]))
#'       }
#'       if (!near(model$getattr("LB")[vars[i]], 0)){ #model$getattr("LB")[vars[i]] < (-1*tol_)
#' 
#'         model$setattr("Obj", setNames(1.0, vars[i]))
#'         model$set_model_sense(minimize=TRUE)
#'         model$optimize()
#'         lp_calls <- lp_calls + 1
#'         sol <- model$get_solution()
#'         global_max <- pmax(global_max, sol$X)
#'         global_min <- pmin(global_min, sol$X)
#'         flux[lp_calls%%stored_obs,] <- sol$X
#' 
#'         if (!near(global_min[i], 0) | !near(global_min[i], 0)){
#'           fixed_val <- rxn_fix(global_max[i], global_min[i])
#'         }
#' 
#'         model$setattr("Obj", setNames(0.0, vars[i]))
#'       }
#'       if (near(global_max[i], 0) & near(global_min[i], 0)){ #(abs(global_max[i]) < tol) & (abs(global_min[i]) < tol)
#'         blocked[i] <- TRUE
#' 	      active[i] <- FALSE
#'         next
#'       }
#'     }
#' 
#'     # set new bounds for selected rxn (temporarily)
#'     model$setattr("UB", setNames(fixed_val + 0.0*fix_tol_frac*abs(fixed_val), vars[i]))
#'     model$setattr("LB", setNames(fixed_val - 0.0*fix_tol_frac*abs(fixed_val), vars[i]))
#' 
#'     # couple reaction to itself if not blocked
#'     if (!blocked[i]){
#'       coupled[i,i] <- TRUE
#'       active[i] <- FALSE
#'     }
#' 
#'     if (idx == length(reaction_indexes)){break}
#' 
#'     for (idx2 in (idx+1):length(reaction_indexes)) { # (j in (i+1):n)
#'       # also keep this in passed in idxs (idx2 in (idx+1):length(reaction_indexes)); j <-  reaction_indexes[idx2]
#'       j <-  reaction_indexes[idx2]
#'       # check for fixed or blocked
#'       if (!active[j] | blocked[j]){next}
#'       if (not_fixed(sub_max[j], sub_min[j])){next}
#'       
#'       # check for uncoupled via correlation
#'       if (cor_check){
#'         if (!correlation_check(flux, i, j)){next}
#'       
#'       }
#' 
#'       skip <- FALSE
#' 
#'       max <- 0
#'       min <- 0
#' 
#'       model$setattr("Obj", setNames(1.0, vars[j]))
#'       if (!skip) {
#'         model$set_model_sense(maximize=TRUE)
#'         model$optimize()
#'         lp_calls <- lp_calls + 1
#'         sol <- model$get_solution()
#'         global_max <- pmax(global_max, sol$X)
#'         global_min <- pmin(global_min, sol$X)
#'         sub_max <- pmax(sub_max, sol$X)
#'         sub_min <- pmin(sub_min, sol$X)
#' 
#'         flux[lp_calls%%stored_obs,] <- sol$X
#' 
#'         max <- sol$X[j]
#' 
#'         skip <- not_fixed(sub_max[j], sub_min[j])
#'       }
#' 
#'       if (!skip) {
#'         model$set_model_sense(minimize=TRUE)
#'         model$optimize()
#'         lp_calls <- lp_calls + 1
#'         sol <- model$get_solution()
#'         global_max <- pmax(global_max, sol$X)
#'         global_min <- pmin(global_min, sol$X)
#'         sub_max <- pmax(sub_max, sol$X)
#'         sub_min <- pmin(sub_min, sol$X)
#' 
#'         flux[lp_calls%%stored_obs,] <- sol$X
#' 
#'         min <- sol$X[j]
#' 
#'         skip <- not_fixed(sub_max[j], sub_min[j])
#'       }
#' 
#'       if (near(max, 0) & near(min, 0)){skip = TRUE}
#' 
#'       if (!skip) { # finally label as coupled
#'         coupled[i,j] <- TRUE
#'         coupled[j,j] <- TRUE # make sure the reaction doesn't look blocked
#'         active[j] <- FALSE
#' 
#'         if (compare_mtx){
#'           output <- known_set_coupling(i, j, coupled, active)
#'           coupled <- output$coupled
#'           active <- output$active
#'         }
#'       }
#' 
#'       model$setattr("Obj", setNames(0.0, vars[j]))
#'     }
#' 
#' 
#'     # unfix i
#'     model$setattr("UB", prev_ub)
#'     model$setattr("LB", prev_lb)
#'   }
#' 
#'   model$setattr("Obj", prev_obj)
#'   model$setattr("ModelSense", prev_sense)
#' 
#'   print(lp_calls)
#' 
#'   list(
#'     coupled = coupled,
#'     lp_calls = lp_calls
#'   )
#' }

#' given a boolean coupling matrix, ensure that all reactions in each set are indicated to be coupled to each other
#' @param coupled boolean matrix
#' @return processed boolean matrix
#' @seealso 
#' @export
#' @examples
#' 
fill_coupling_matrix <- function(coupled){
  rows <- nrow(coupled)

  for (i in 1:nrow(coupled)){
    #identify set
    # if (!coupled[i,i]){next}
    set <- which(coupled[i,]) # true values in row
    if (length(set) < 1){next}
    set <- unique(c(i, set))
    coupled[set,set] <- TRUE

  }

  coupled[lower.tri(coupled)] <- FALSE
  return(coupled)
}

#' given a boolean coupling matrix, ensure that all reactions in each set are indicated to be coupled to each other based on known sets
#' @param mtx boolean coupling matrix to fill
#' @param sets sets to indicate in the matrix
#' @return processed boolean matrix
#' @seealso 
#' @export
#' @examples
#' 
fill_coupling_matrix_from_sets <- function(mtx, sets){
  for (set in sets){
    mtx[set,set] <- TRUE
  }
  
  mtx <- fill_coupling_matrix(mtx)
  
  return(mtx)
}

set_vector <- function(coupled){
  active <- matrix(data = TRUE, nrow = 1, ncol = ncol(coupled))
  set_num <- matrix(data = 0, nrow = 1, ncol = ncol(coupled))

  set_iter <- 1
  for (i in nrow(coupled)){
    if (!active[i]){next} # skip if already in a set
    set <- which(coupled[i,])
    if (length(set) < 1){next} # skip if blocked reaction (will not be coupled to itself)
    active[set] <- FALSE
    set_num[set] <- set_iter # enter same set # into all indexes in the same set
    set_iter <- set_iter + 1
  }

  return(set_num)
}
