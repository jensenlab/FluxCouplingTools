#library(readr)
# scripts for loading different models

# HOW TO HANDLE YEAST MODEL
# read txt to data frame
# use minval to convert to tsv for sybil (writeTSVmod)
# use sybil to read tsv (readTSVmod)


## NECESSARY PRESETS FOR USAGE:
# model
# S

## ECOLI MODEL

#' description
#' @param 
#' @return
#' @seealso 
#' @export
#' @examples
#' 
get_ecoli_model <- function(){
  data(Ec_core);
  model=Ec_core;
  model <- changeBounds(model, 11, lb = 0) # this and next lead to no PGI in GND coset
  # model <- changeBounds(model, 13, lb = 0, ub = 0)
  model <- rmReact(model = model, react = 13)
  for (i in findExchReact(model)@react_pos){
    model <- changeBounds(model, i, lb = -1000, ub = 1000)
    # if (model@lowbnd[i] == 0){
    #   model <- changeBounds(model, i, lb = -1000)
    # }
  }
  return(model)
}


#solver="cplexAPI"
# solver="glpkAPI"
# W=8000
# warmup=8000
# nPnts=3000
# steps=2

# opt <- fluxVar(model, percentage = 99)
# model_fva <- opt@lp_obj
# fva_min <- model_fva[1:95]
# fva_max <- model_fva[96:190]

## YEAST MODEL

# USE THIS FUNCTION FOR YEAST MODEL (MAIN)
#' generate Saccharomyces cerevisiae model
#' @param reactList filename for list of reactions
#' @param metList filename for list of metabolites
#' @return Sybil model
#' @seealso 
#' @export
#' @examples
#' 
get_yeast_maranas_model <- function(reactList = "S7 Model iSce926.tsv", metList = "S7 Model iSce926_met.tsv"){

  yeast_model <- readTSVmod(reactList = reactList, metList = metList)
                            # , remUnusedMetReact = FALSE, balanceReact = TRUE)

  # get rid of biomass and lipid pseudo reactions except for one
  remove <- which(grepl("pseudoreaction", yeast_model@react_name)) # should have 4 entries; removing first 3

  # yeast_model <- rmReact(model = yeast_model, react = remove[4])
  yeast_model <- rmReact(model = yeast_model, react = remove[3])
  yeast_model <- rmReact(model = yeast_model, react = remove[2])
  # yeast_model <- rmReact(model = yeast_model, react = remove[1]) # lipid

  growth_idx <- which(yeast_model@react_name == 'growth')
  biomass_idx <- which(yeast_model@react_name == "yeast 8 biomass pseudoreaction")
  lipid_idx <- which(yeast_model@react_name == "lipid pseudoreaction")

  yeast_model <- add_exch_rxns_to_model(yeast_model, lipid_idx)
  yeast_model <- add_exch_rxns_to_model(yeast_model, biomass_idx)

  yeast_model <- rmReact(model = yeast_model, react = biomass_idx, rm_met = TRUE) # last biomass reaction
  yeast_model <- rmReact(model = yeast_model, react = growth_idx, rm_met = TRUE) # growth reaction (objective function; biomass exchange)
  yeast_model <- rmReact(model = yeast_model, react = lipid_idx, rm_met = TRUE) # lipid reaction (objective function; biomass exchange)

  return(yeast_model)
}

## MUTANS MODEL
#' generate Streptococcus mutans model
#' @param reactList filename for list of reactions
#' @param metList filename for list of metabolites
#' @return Sybil model
#' @seealso 
#' @export
#' @examples
#' 
get_mutans_model <- function(reactList = "mutans_model.csv", metList = "mutans_model_met.csv"){
  mutans_model <- readTSVmod(reactList = reactList, metList = metList)

  mutans_model@met_name[401] <- "DAP-type peptidoglycan"
  mutans_model@met_name[422] <- "Lys-type peptidoglycan"

  # make sure biomass split is implemented
  exch_idxs <- which(mutans_model@S[,477] != 0) # biomass
  biom_consumed <- which(mutans_model@S[,477] < 0)
  biom_produced <- which(mutans_model@S[,477] > 0)

  exch <- findExchReact(mutans_model)
  add_exch <- c()
  for (i in exch_idxs){
  
    if (mutans_model@met_id[i] %in% c('C00002', 'C00044')){
      next
    }

    else {
      add_exch <- c(add_exch, i)
    }
  }

  biom_consumed <- intersect(add_exch, which(mutans_model@S[,477] < 0))
  biom_produced <- intersect(add_exch, which(mutans_model@S[,477] > 0))

  # add outlets for reactants in biomass rxn
  mutans_model <- addExchReact(mutans_model, met <- mutans_model@met_id[biom_consumed],
                               lb <- rep(0, length(biom_consumed)), ub <- rep(1000, length(biom_consumed)))

  # add new inlets for products with energy cost
  mutans_model <- addReact(mutans_model, 'ATP_decomp_new', met = c('C00002', 'C00008', 'C00009', 'C00080'),
                           Scoef = c(-1, 1, 1, 1), reversible = FALSE, lb = 0, ub = 1000)
  mutans_model <- addReact(mutans_model, 'GTP_decomp_new', met = c('C00044', 'C00035', 'C00013'),
                           Scoef = c(-1, 1, 1), reversible = FALSE, lb = 0, ub = 1000)
  mutans_model <- rmReact(model = mutans_model, react = 477)

  # remove duplicate reactions
  print(mutans_model@react_name[616])
  mutans_model <- rmReact(model = mutans_model, react = 616)
  print(mutans_model@react_name[366])
  mutans_model <- rmReact(model = mutans_model, react = 366)
  print(mutans_model@react_name[364])
  mutans_model <- rmReact(model = mutans_model, react = 364)

  return(mutans_model)
}

get_mutans_model_w_obj <- function(){
  mutans_model <- readTSVmod(reactList = "mutans_model.csv", metList = "mutans_model_met.csv")
  
  # print('biomass')
  # print(mutans_model@S[which(mutans_model@S[,477] != 0), 477])
  
  mutans_model@met_name[401] <- "DAP-type peptidoglycan"
  mutans_model@met_name[422] <- "Lys-type peptidoglycan"
  
  mutans_model@obj_coef[477] <- 1
  
  return(mutans_model)
}

#' generate Pseudomonas Aeruginosa model
#' @param reactList filename for list of reactions
#' @return Sybil model
#' @seealso 
#' @export
#' @examples
#' 
get_PAO_model <- function(reactList = "PAO1recon1_v23_reactions.csv"){
  
  # need to change file names
  pao_model <- readTSVmod(reactList = reactList, fielddelim = "\t") # metList = "PAO1recon1_v23_metabolites.csv", 
  
  reactions_to_replace <- c('PAO1_Biomass', 'PA_Biomass_v13ub', 'PA_Biomass_v13', 'PA_Biomass_v4')
  replacement_idxs <- which(pao_model@react_id %in% reactions_to_replace)
  
  exch_idxs <- c() # biomass
  consumed_met_idxs <- c()
  produced_met_idxs <- c()
  
  for (rxn in replacement_idxs){
    exch_idxs <- c(exch_idxs, which(pao_model@S[,rxn] != 0)) # biomass
    consumed_met_idxs <- c(consumed_met_idxs, which(pao_model@S[,rxn] < 0))
    produced_met_idxs <- c(produced_met_idxs, which(pao_model@S[,rxn] > 0))
  }
  consumed_mets <- pao_model@met_id[consumed_met_idxs]
  produced_mets <- pao_model@met_id[produced_met_idxs]
  
  exch <- findExchReact(pao_model)
  
  # add outlets for consumed metabolites
  for (met in consumed_mets){
    if (met %in% exch@met_id){
      exch_idx <- exch@react_pos[which(exch@met_id == met)]
      pao_model@uppbnd[exch_idx] <- 1000
      next
    }
    pao_model <- addExchReact(pao_model, met <- met, lb <- 0, ub <- 1000)
  }
  # add inlets or produced metabolites
  for (met in produced_mets){
    if (met %in% exch@met_id){
      exch_idx <- exch@react_pos[which(exch@met_id == met)]
      pao_model@lowbnd[exch_idx] <- -1000
      next
    }
    pao_model <- addExchReact(pao_model, met <- met, lb <- -1000, ub <- 0)
  }
  
  # remove biomass reactions
  pao_model <- rmReact(model = pao_model, react = replacement_idxs)
  
  return(pao_model)
}

get_blocked <- function(model){
  lb <- which(model@lowbnd > -1000)
  ub <- which(model@uppbnd < 1000)

  potential_blocked <- intersect(lb, ub)
  blocked <- c()
  for (i in potential_blocked){
    if (model@lowbnd[i] == 0 & model@uppbnd[i] == 0){
      blocked <- c(blocked, i)
    }
  }

  return(blocked)
}

generate_exch_rxn <- function(model, rxn_idx){
  rxn_id <- model@met_id[rxn_idx]
  rxn_name <- model@met_name[rxn_idx]
  paste("exc00000	", rxn_name," exchange	", " <=> ", rxn_id, "					1	0	0	0		Unknown		", sep = "")
}

#' add exchange reactions for metabolites in indicated reaction to sybil model
#' @param model sybil model to update
#' @param rxn_idx integer indicating reaction whose metabolites should have exchange reactions
#' @return sybil model
#' @seealso 
#' @export
#' @examples
#' 
add_exch_rxns_to_model <- function(model, rxn_idx){

  exch_idxs <- which(model@S[,rxn_idx] != 0)
  consumed <- which(model@S[,rxn_idx] < 0)
  produced <- which(model@S[,rxn_idx] > 0)

  exch <- findExchReact(model)
  add_exch <- c()
  for (i in exch_idxs){
    if ((model@met_id[i] %in% exch@met_id) | (paste(model@met_id[i], '[e]', sep = "") %in% exch@met_id)){
      print(paste('existing exch:', model@met_id[i]))
      next
    }

    else {
      add_exch <- c(add_exch, i)
    }
  }

  consumed <- intersect(add_exch, which(model@S[,rxn_idx] < 0))
  produced <- intersect(add_exch, which(model@S[,rxn_idx] > 0))

  # add outlets for reactants in biomass rxn
  model <- addExchReact(model, met <- model@met_id[consumed],
                               lb <- rep(0, length(consumed)), ub <- rep(1000, length(consumed)))
  # add inlets for products
  model <- addExchReact(model, met <- model@met_id[produced],
                               lb <- rep(-1000, length(produced)), ub <- rep(0, length(produced)))
  return(model)
}
