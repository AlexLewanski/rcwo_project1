###################################
### INITIAL PEDIGREE PROCESSING ###
###################################
ped_add_missing_indivs <- function(ped, id, fid, mid, sex, founder_par_val = '0') {
  
  ped_col_names <- colnames(ped)
  parent_vec <- unique(c(ped[, fid, drop = TRUE],
                         ped[, mid, drop = TRUE]))
  
  parent_vec <- parent_vec[!parent_vec %in% founder_par_val]
  
  absent_parent <- parent_vec[!parent_vec %in% ped[, id, drop = TRUE]]
  
  if (length(absent_parent) == 0) return(ped)
  
  sex_vec <- sapply(absent_parent, function(id, df) {c(1:2)[apply(df, 2, function(x) id %in% x)]},
                    df = ped[,c('fid', 'mid')])
  
  new_id_df <- stats::setNames(data.frame(absent_parent, founder_par_val, founder_par_val, sex_vec),
                               nm = c(id, fid, mid, sex))
  
  new_id_df[ped_col_names[!ped_col_names %in% colnames(new_id_df)]] <- NA
  
  new_id_df_order <- new_id_df[, match(ped_col_names, colnames(new_id_df))]
  
  return(rbind(ped, new_id_df_order))
}


col_process_internal <- function(ped,
                                 id_col,
                                 sire_col,
                                 dam_col,
                                 sex_col,
                                 keep_extra_cols) {
  
  colname_vec <- stats::setNames(c('id', 'fid', 'mid', 'sex'),
                                 nm = c(id_col, sire_col, dam_col, sex_col))
  
  if (any(!names(colname_vec) %in% colnames(ped))) {
    stop('The following column(s) are missing from the pedigree: ',
         paste(names(colname_vec)[!names(colname_vec) %in% colnames(ped)], collapse = ', '))
  }
  
  #rename columns
  colnames(ped)[match(names(colname_vec), colnames(ped))] <- colname_vec
  
  if (keep_extra_cols) {
    colname_vec <- c(colname_vec, colnames(ped)[!colnames(ped) %in% colname_vec])
  }
  
  return(ped[,colname_vec,drop=FALSE])
}

sex_eval_internal <- function(processed_ped) {
  
  ### dam/sire overlap ###
  sex_overlap <- intersect(processed_ped$fid, processed_ped$mid)
  
  if (length(sex_overlap[sex_overlap != '0']) > 0) {
    stop('The following individuals are recorded as both a sire and dam: ',
         paste(unique(sex_overlap[sex_overlap != '0']), collapse = ', '))
  }
  
  ### identifying conflicts between parent status and specified sex ###
  males_unknown <- processed_ped$id[processed_ped$sex %in% c(1, 0)]
  females_unknown <- processed_ped$id[processed_ped$sex %in% c(2, 0)]
  
  male_unk_dams <- males_unknown[males_unknown %in% processed_ped$mid]
  female_unk_sires <- females_unknown[females_unknown %in% processed_ped$fid]
  combined_sex_mismatch <- c(male_unk_dams, female_unk_sires)
  
  if (length(combined_sex_mismatch) > 0) {
    
    sex_correct_vec <- c(rep(2, length(male_unk_dams)),
                         rep(1, length(female_unk_sires)))
    
    processed_ped$sex[match(combined_sex_mismatch, processed_ped$id)] <- sex_correct_vec
  }
  
  return(processed_ped)
}


# Process a pedigree
#
#ped The unprocessed dataframe. It needs to contain the following columns: id_col (the ID of each individual), sire_col (the ID of the sire), dam_col (the ID of the dam), and sex_col (the sex of the individual).
#id_col The name of the id column.
#sire_col The name of the sire (father) column.
#dam_col The name of the dam (mother) column.
#sex_col The name of the sex column.
#founder_val One or more (stored as a vector) values that represent the value for founder parents.
#sex_vals A list containing the values representing each sex in sex_col. The list should be in the following format: list(male = c('male_val1', male_val2'), female = 'female_val', unknown = c('u1', 'u2', NA)).
#keep_extra_cols Logical value (TRUE or FALSE) indicating whether to keep extra columns in the pedigree beyond the id, sire, dam, and sex columns.
#disable_sex_check Logical value (TRUE or FALSE) indicating whether sex information should be evaluated.

process_ped <- function(ped,
                        id_col,
                        sire_col,
                        dam_col,
                        sex_col,
                        founder_val,
                        sex_vals,
                        keep_extra_cols = TRUE,
                        disable_sex_check = FALSE) {
  
  #===============
  #== OVERVIEW ===
  #===============
  
  #make the ped object a dataframe (if it isn't already)
  #rename and reorder the id, parent, and sex columns
  #recode the founder parent values to '0'
  #add individuals who are parents but aren't included in the id column
  #process and evaluate sex info:
  #   - recode sex categories to 1 (male), 2 (female), and 0 (unknown)
  #   - make sure single individuals aren't recorded as both sire and dams
  #   - if an individual has as unknown sex but is included as a sire or dam,
  #     update the sex info based on the sex of the parent type that it was
  #     recorded as
  
  
  #================================================
  #== CHECKING AND INITIAL PROCESSING OF INPUTS ===
  #================================================
  if (!identical(class(ped), 'data.frame')) ped <- as.data.frame(ped)
  stopifnot(identical(class(ped), 'data.frame'))
  
  for (i in list(id_col, sire_col, dam_col, sex_col)) {
    stopifnot(length(i) == 1 && class(i) == 'character')
  }
  
  stopifnot(is.list(sex_vals) && identical(sort(names(sex_vals)), c("female", "male", "unknown") ))
  
  #this should get all the argument inputs in a list but idk if this a good idea
  #c(as.list(environment()), list(...))
  logic_list <- list(keep_extra_cols, disable_sex_check)
  names(logic_list) <- c('keep_extra_cols', 'disable_sex_check')
  
  for (i in names(logic_list)) {
    if (! (isFALSE(logic_list[[i]]) | isTRUE(logic_list[[i]])) )
      stop('The following argument is not a single TRUE or FALSE value: ', i)
  }
  
  if (any(duplicated(ped[,id_col,drop = TRUE]))) {
    stop('The following individuals have >1 entry in the id column: ',
         paste(
           unique(ped[,id_col,drop = TRUE][duplicated(ped[,id_col,drop = TRUE])]),
           collapse = ', ')
    )
  }
  
  
  #========================
  #== COLUMN PROCESSING ===
  #========================
  
  #rename and reorder columns
  ped_processed <- col_process_internal(ped = ped,
                                        id_col = id_col,
                                        sire_col = sire_col,
                                        dam_col = dam_col,
                                        sex_col = sex_col,
                                        keep_extra_cols = keep_extra_cols)
  
  #if they aren't already, coerce id, fid, and mid to characters
  for (COL in c('id', 'fid', 'mid')) {
    if (!is.character(ped_processed[,COL,drop=TRUE])) {
      ped_processed[,COL] <- as.character(ped_processed[,COL,drop=TRUE])
    }
  }
  
  
  #===============================
  #== PEDIGREE INFO PROCESSING ===
  #===============================
  
  for (i in c('fid', 'mid')) {
    ped_processed[,i][ped_processed[,i,drop = TRUE] %in% founder_val] <- "0"
  }
  
  if (!all(ped_processed$sex %in% unlist(sex_vals)))
    stop('Sex is encoded with different(or additional) values than those in the sex_vals list')
  ped_processed$sex_old <- ped_processed$sex
  ped_processed$sex <- NA
  
  for (i in names(sex_vals) ) {
    ped_processed$sex[ped_processed$sex_old %in% sex_vals[[i]]] <- switch(i,
                                                                          male = 1,
                                                                          female = 2,
                                                                          unknown = 0)
  }
  
  ped_processed <- ped_add_missing_indivs(ped = ped_processed,
                                          id = 'id',
                                          fid = 'fid',
                                          mid = 'mid',
                                          sex = 'sex',
                                          founder_par_val = '0')
  
  if (!disable_sex_check) {
    ped_processed <- sex_eval_internal(processed_ped = ped_processed)
  }
  
  return(ped_processed)
}



#############################
### ALL-PURPOSE FUNCTIONS ###
#############################
harmonic_mean <- function(x) length(x)/sum(1/x)

#######################################
### PEDIGREE MANIPULATON AND CHECKS ###
#######################################
# Reorder pedigree so that all offspring follow their parents
#
#ped a pedigree
#id_col name or index of id column
#sire_col name or index of sire column
#dam_col name or index of dam column
#
# a pedigree ordered so that offspring follow their parents

reorder_ped <- function(ped,
                        id_col,
                        sire_col,
                        dam_col) {
  
  #Zhang et al. 2009
  #An algorithm to sort complex pedigrees chronologically without birthdates
  
  id_ind <- seq_along(ped[[id_col]])
  id_sire <- match(ped[[sire_col]], ped[[id_col]])
  id_dam <- match(ped[[dam_col]], ped[[id_col]])
  
  gen_vec <- rep(0, nrow(ped))
  
  while (TRUE) {
    empty_list <- list()
    
    for (i in id_ind) {
      
      if (!is.na(id_sire[i]) && !(id_sire[i] %in% unlist(empty_list))) {
        empty_list[[length(empty_list) + 1]] <- id_sire[i]
        gen_vec[id_sire[i]] <- gen_vec[id_sire[i]] + 1
      }
      
      if (!is.na(id_dam[i]) && !(id_dam[i] %in% unlist(empty_list))) {
        empty_list[[length(empty_list) + 1]] <- id_dam[i]
        gen_vec[id_dam[i]] <- gen_vec[id_dam[i]] + 1
      }
      
    }
    
    if (length(empty_list) == 0) break
    
    id_ind <- unlist(empty_list)
  }
  
  return(ped[order(gen_vec, decreasing = TRUE),])
}


# Evaluate whether the pedigree is ordered. NOTE: this has not been thoroughly vetted.
#
#ped a ped
#id_col name or index of id column
#sire_col name or index of sire column
#dam_col name or index of dam column
#
# A logical value indicating whether or not the pedigree is ordered with offspring after parents.
#
#
check_order_ped <- function(ped,
                            id_col,
                            sire_col,
                            dam_col) {
  
  fid_index <- match(ped[[sire_col]], ped[[id_col]])
  mid_index <- match(ped[[dam_col]], ped[[id_col]])
  
  for (i in seq_len(nrow(ped))) {
    if (isTRUE(i < fid_index[i]) | isTRUE(i < mid_index[i])) {
      return(FALSE)
    }
  }
  return(TRUE)
}



# Create indexed version of the pedigree
#
#ped pedigree (stored in a dataframe) with the following organization: column 1 --> id, column 2 --> sire id, column3 --> dam id. Founders should have parents coded as 0.
#
# A dataframe containing the pedigree with ids transformed to the index values of each individual (and founder parents coded as 0s).
#
#
index_pedigree <- function(ped) {
  
  return(
    cbind(seq_along(ped[,1,drop=TRUE]),
          match(ped[,2,drop=TRUE], ped[,1,drop=TRUE], nomatch = 0),
          match(ped[,3,drop=TRUE], ped[,1,drop=TRUE], nomatch = 0))
  )
}




########################################
### MISCELLANEOUS PEDIGREE FUNCTIONS ###
########################################

# Switch character IDs to numeric (for the ID, sire, and dam columns)
#
#ped pedigree
#id_col name or index of id column
#sire_col name or index of sire column
#dam_col name or index of dam column
#founder_val the value that represents a founder in the sire and dam cols
#
# a pedigree stored in a dataframe with updated individual ID values
#
#
character2numeric_id <- function(ped, id_col, sire_col, dam_col, founder_val = 0) {
  
  id_map <- data.frame(unique(ped[,id_col, drop = TRUE]),
                       seq_along(unique(ped[,id_col, drop = TRUE])))
  
  for (i in c(id_col, sire_col, dam_col)) {
    colnames(id_map) <- c(i, paste0(i, '_numeric'))
    ped <- dplyr::left_join(ped, id_map, by = i)
    ped[,paste0(i, '_numeric')][is.na(ped[,paste0(i, '_numeric')])] <- founder_val
  }
  
  return(ped)
}




#######################################
### EXTRACTION OF PEDIGREE ELEMENTS ###
#######################################
# Extract ancestors for one or more individuals
#
#ped a pedigree object
#indiv_vec the individuals that for which to extract ancestors
#include_repeats if individuals are duplicate ancestors, should duplicates be removed?
#
# a vector of ancestors
#
#
get_ancestors <- function(ped, indiv_vec, include_repeats = TRUE) {
  
  #if (any(c('tbl_df', "tbl") %in% class(ped))) ped <- as.data.frame(ped)
  
  #extract all the ancestors for each individual in indiv_vec
  anc_vec <- unlist(ped[match(indiv_vec, ped[,1, drop = TRUE]), c(2, 3), drop = FALSE],
                    use.names = FALSE)
  #anc_vec_removeNA <- anc_vec[!is.na(anc_vec)] #remove NAs (these are founders)
  anc_vec_removeNA <- anc_vec[anc_vec != "0"]
  #if len of anc vec is 0 after removing NAs, no ancestors exist for the indivs
  if (length(anc_vec_removeNA) == 0) return(NULL)
  
  #output --> include_repeats != TRUE, return vector of ancestors without reps
  if (!isTRUE(include_repeats)) return(unique(anc_vec_removeNA))
  return(anc_vec_removeNA)
}


get_offspring <- function(ped, indiv_vec, include_repeats = TRUE) {
  
  offspring_vec <- ped[,1,drop = TRUE][ped[,2,drop = TRUE] %in% indiv_vec | ped[,3,drop = TRUE] %in% indiv_vec]
  if (length(offspring_vec) == 0) return(NULL)
  return(offspring_vec)
}


# Extract the founders from a pedigree
#
#ped a pedigree object
#founder_vals the founder values
#
# a list with the founders and partial founders (the individuals with only one known parent)
#
#
get_founders <- function(ped, founder_vals = '0') {
  founder_list <- list()
  founder_list[['founders']] <- ped[apply(ped[,c(2,3)], 1, function(x) all(x %in% '0')),1,drop=TRUE]
  founder_list[['partial_founders']] <- ped[apply(ped[,c(2,3)], 1, function(x) sum(x %in% '0') == 1),1,drop=TRUE]
  founder_list[sapply(founder_list, function(x) length(x) == 0)] <- 'none'
  
  return(founder_list)
}



# Subset pedigree down to a focal set of individuals and their descendants
#
#ped pedigree (stored in dataframe) with the following organization: column 1 --> id, column 2 --> sire id, column 3 --> dam. Founder parents should be coded as 0s.
#indivs the individuals for whom you want to subset the pedigree
#
# a dataframe containing pedigree information for the focal individuals and their descendants based on the input pedigree
subset_ped_descendants <- function(ped, indivs) {
  
  indiv_list <- list(indivs)
  counter <- 1
  while (TRUE) {
    
    indiv_vec <- get_offspring(ped, indiv_list[[counter]])
    
    if (is.null(indiv_vec)) break
    
    indiv_list[[counter + 1]] <- indiv_vec[!indiv_vec %in% unlist(indiv_list)]
    
    counter <- counter + 1
  }
  
  return(unlist(indiv_list))
  
}


# Subset pedigree down to a focal set of individuals and their ancestors
#
#ped pedigree (stored in dataframe) with the following organization: column 1 --> id, column 2 --> sire id, column 3 --> dam. Founder parents should be coded as 0s.
#indivs the individuals for whom you want to subset the pedigree
#
# a dataframe containing pedigree information for the focal individuals and their ancestors based on the input pedigree
#
subset_ped_ancestors <- function(ped, indivs) {
  
  indiv_list <- list(indivs)
  counter <- 1
  while (TRUE) {
    
    indiv_vec <- get_ancestors(ped, indiv_list[[counter]],
                               include_repeats = FALSE)
    
    if (is.null(indiv_vec)) break
    
    indiv_list[[counter + 1]] <- indiv_vec[!indiv_vec %in% unlist(indiv_list)]
    
    counter <- counter + 1
  }
  
  return(unlist(indiv_list))
  
}

###########################
############################
### GENE DROP SIMULATION ###
############################
############################

#library(pedigree)
#library(tidyverse) #specifically tibble


#################
### FUNCTIONS ###
#################

#Single locus gene drop on a pedigree
#
#ped pedigree on which to perform gene drop
#id_col name or index of id column
#sire_col name or index of sire column
#dam_col name or index of dam column
#sex_col name or index of sex column
#sims number of gene drop simulations
#pop_freq the frequency of the allele in the starting population
#fixed_founder_genotypes should founder genotypes be fixed?
#calc_contribution should the ancestry contributions of each founder be recorded?
#report_progress should the simulation progress be shown?
simple_gene_drop <- function(ped,
                             id_col,
                             sire_col,
                             dam_col,
                             sex_col,
                             sims = 10,
                             pop_freq = 0.5,
                             fixed_founder_genotypes = FALSE,
                             calc_contribution = FALSE,
                             report_progress = TRUE) {
  
  #==========================
  #=== INITIAL PROCESSING ===
  #==========================
  
  ### CHECK INPUTS ###
  report_progress <- isTRUE(report_progress)
  fixed_founder_genotypes <- isTRUE(fixed_founder_genotypes)
  #OTHER CHECKS
  
  if (report_progress) message('Initial pedigree processing and simulation preparation.')
  
  ### INITIAL PROCESSING OF PEDIGREE ###
  #On my computer, feeding a tibble to pedigree's orderPed results in the R session
  #being aborted. This doesn't happen if ped is a dataframe
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)
  
  # (1) subset down to id, sire, and dam cols
  # (2) rename the columns
  # (3) reorder pedigree so that parents appear before offspring
  # (4) record size of pedigree
  ped_reorder <- ped[,c(id_col, sire_col, dam_col, sex_col)]
  colnames(ped_reorder) <- c('id', 'sire', 'dam', 'sex')
  
  #ordered_ped <- ped_reorder[order(pedigree::orderPed(ped_reorder)),]
  ordered_ped <- reorder_ped(ped_reorder,
                             id_col = 'id',
                             sire_col = 'sire',
                             dam_col = 'dam')
  ped_size <- nrow(ped_reorder)
  
  #record the index of sire and dams in new cols (this will make retrieval
  #of sire and dam genotype and genotype origin information easier and faster)
  ordered_ped$sire_id <- match(ordered_ped$sire, ordered_ped$id)
  ordered_ped$dam_id <- match(ordered_ped$dam, ordered_ped$id)
  
  
  ### CREATE GENOTYPE AND GENOTYPE ORIGIN MATRICES ###
  geno_mat <- matrix(data = NA, nrow = ped_size, ncol = 2,
                     dimnames = list(ordered_ped$id, c('sire_geno', 'dam_geno')))
  
  geno_orig_mat <- matrix(data = NA, nrow = ped_size, ncol = 2,
                          dimnames = list(ordered_ped$id, c('sire_geno_origin', 'dam_geno_origin')))
  
  if (isTRUE(calc_contribution)) {
    contr_mat <- matrix(data = 0, nrow = ped_size, ncol = 2,
                        dimnames = list(ordered_ped$id, c('sire_contr', 'dam_contr')))
    
    contr_list <- list()
  }
  
  
  ### CREATE GENOTYPE INFORMATION FOR FOUNDERS###
  #*** THIS IS AN OBVIOUS PLACE FOR IMPROVEMENT. HOW DO OTHERS DO THIS?
  
  #indices of founders
  founder_indices <- which(apply(ordered_ped[, c('sire', 'dam')], 1, function(x) all(x == '0')))
  
  #IF A SINGLE FIXED SET OF FOUNDER ALLELES ACROSS SIMULATIONS
  if (fixed_founder_genotypes) {
    #add genotype information for founders by two independent draws from the binomial
    #geno_mat[founder_indices,] <- t(replicate(length(founder_indices), rbinom(n = 2, size = 1, prob = 0.5), simplify = TRUE))
    geno_mat[founder_indices,] <- matrix(stats::rbinom(n = length(founder_indices)*2, size = 1, prob = pop_freq), ncol = 2)
  }
  
  #for both alleles in each founder, record the origin as the founder id (not sure if this is right)
  geno_orig_mat[founder_indices,] <- t(sapply(rownames(geno_orig_mat)[founder_indices], function(x) paste(x, 1:2, sep = "_"))) #separate allele identities for each founder
  #geno_orig_mat[founder_indices,] <- rownames(geno_orig_mat)[founder_indices] #if you want to have the same identity for both alleles of each founder
  
  #=============================
  #=== GENE DROP SIMULATIONS ===
  #=============================
  
  #initiate progress bar
  if (report_progress) {
    message('Starting gene drop simulation', ifelse(sims > 1, 's.', '.'))
    prog_bar <- utils::txtProgressBar(min = 0, max = sims, initial = 0, char = "*", style = 3)
  }
  
  gene_drop_list <- list() #list to store each gene drop simulation
  
  for (i in seq_len(sims)) {
    
    #IF FOUNDERS SHOULD BE RANDOMLY DRAWN FOR EACH SIMULATION (i.e. if fixed_founder_genotypes is FALSE)
    if (!fixed_founder_genotypes) {
      #add genotype information for founders by two independent draws from the binomial
      geno_mat[founder_indices,] <- matrix(stats::rbinom(n = length(founder_indices)*2,
                                                         size = 1,
                                                         prob = pop_freq),
                                           ncol = 2)
    }
    
    #transmission probabilities for each allele from each individual's parents
    #for example, the ith element in sire_draw_vec indicates whether the ith individual will inherit
    #allele 1 (sire allele) or allele 2 from its father.
    sire_draw_vec <- stats::rbinom(n = ped_size, size = 1, prob = 0.5) + 1L
    dam_draw_vec <- stats::rbinom(n = ped_size, size = 1, prob = 0.5) + 1L
    
    f_vec_mat <- matrix(data = 0, nrow = ped_size, ncol = 2,
                        dimnames = list(ordered_ped$id, c('sire_ibd_count', 'dam_ibd_count')))
    
    #perform gene drop
    gene_drop_list[[paste0('sim', i)]] <- single_gene_drop(ped_size = ped_size,
                                                           pedigree = ordered_ped,
                                                           sire_draw_vec = sire_draw_vec,
                                                           dam_draw_vec = dam_draw_vec,
                                                           geno_mat = geno_mat,
                                                           geno_orig_mat = geno_orig_mat,
                                                           f_vec_mat = f_vec_mat)
    
    if (isTRUE(calc_contribution)) {
      contr_list[[paste0('sim', i)]] <- single_genetic_contr(ped_size = ped_size,
                                                             pedigree = ordered_ped,
                                                             sire_draw_vec = sire_draw_vec,
                                                             dam_draw_vec = dam_draw_vec,
                                                             contr_mat = contr_mat)
    }
    
    if (report_progress) utils::setTxtProgressBar(prog_bar, i)
  }
  
  
  
  #=====================================
  #=== PROCESS AND RETURN SIM OUTPUT ===
  #=====================================
  #THIS CAN BE PRETTY SLOW FOR LARGE NUMBERS OF SIMULATIONS
  #THINK MORE ABOUT HOW TO SPEED THIS UP
  
  if (report_progress) message('\nProcessing simulation output.')
  output_list <- list()
  
  output_list[['geno']] <- lapply(gene_drop_list, function(x) {
    tibble::rownames_to_column(as.data.frame(x[['geno_mat']]), var = "id")
  } ) %>%
    dplyr::bind_rows(.id = 'sim')
  #dplyr::bind_rows(., .id = 'sim')
  
  output_list[['geno_origin']] <- lapply(gene_drop_list, function(x) {
    tibble::rownames_to_column(as.data.frame(x[['geno_orig_mat']]), var = "id")
  }) %>%
    dplyr::bind_rows(.id = 'sim')
  #dplyr::bind_rows(., .id = 'sim')
  
  output_list[['f_count']] <- lapply(gene_drop_list, function(x) {
    tibble::rownames_to_column(as.data.frame(x[['f_vec_mat']]), var = "id")
  }) %>%
    dplyr::bind_rows(.id = 'sim')
  #dplyr::bind_rows(., .id = 'sim')
  
  if (isTRUE(calc_contribution)) {
    output_list[['contr_list']] <- contr_list %>%
      dplyr::bind_rows(.id = 'sim')
    #dplyr::bind_rows(., .id = 'sim')
  }
  
  return(output_list)
  
}



single_gene_drop <- function(ped_size,
                             pedigree,
                             sire_draw_vec,
                             dam_draw_vec,
                             geno_mat,
                             geno_orig_mat,
                             f_vec_mat) {
  
  #iterate through every individual and perform the following tasks:
  # (1) randomly draw an allele from each parent
  # (2) copy the origin of the chosen allele
  
  #notes:
  # - for purposes of speed, the choice of which allele to transmit from
  #   each parent is pre-calculated and fed into the function (sire_draw_vec,
  #   dam_draw_vec)
  # -
  
  for (i in seq_len(ped_size)) {
    
    #if an individual doesn't have a sire or dam id, it is a founder
    #founders should already have genotypes so skip to next individual
    #THIS COULD BE IMPROVED (OR CONSIDERED MORE THOUGHTFULLY). THIS
    #WILL PRODUCE WEIRDNESS IF AN INDIVIDUAL HAS A SINGLE KNOWN PARENT
    if (any(is.na(pedigree[i,c('sire_id', 'dam_id')]))) next
    
    #retrieve sire and dam index for focal individual
    sire_index <- pedigree[i, 'sire_id']
    dam_index <- pedigree[i, 'dam_id']
    
    #get random index (1 or 2) for which sire and dam ID to transmit
    #sire_draw <- rbinom(n = 1, size = 1, prob = 0.5) + 1L #sum(runif(1, 0, 1) > 0.5) + 1L
    #dam_draw <- rbinom(n = 1, size = 1, prob = 0.5) + 1L #sum(runif(1, 0, 1) > 0.5) + 1L
    sire_draw <- sire_draw_vec[i]
    dam_draw <- dam_draw_vec[i]
    
    #add chosen sire and dam alleles for the individual's genotype
    geno_mat[i,1] <- geno_mat[sire_index, sire_draw]
    geno_mat[i,2] <- geno_mat[dam_index, dam_draw]
    
    #record the origin of each allele for the focal individual
    geno_orig_mat[i,1] <- geno_orig_mat[sire_index, sire_draw]
    geno_orig_mat[i,2] <- geno_orig_mat[dam_index, dam_draw]
    
    sire_ibd <- as.integer(length(unique(geno_orig_mat[sire_index,])) == 1)
    dam_ibd <- as.integer(length(unique(geno_orig_mat[dam_index,])) == 1)
    f_vec_mat[i,1] <- f_vec_mat[sire_index, sire_draw] + sire_ibd
    f_vec_mat[i,2] <- f_vec_mat[dam_index, dam_draw] + dam_ibd
    
  }
  
  return(list(geno_mat = geno_mat,
              geno_orig_mat = geno_orig_mat,
              f_vec_mat = f_vec_mat))
  
}




single_genetic_contr <- function(ped_size,
                                 pedigree,
                                 sire_draw_vec,
                                 dam_draw_vec,
                                 contr_mat) {
  
  contribute_list <- list()
  
  for (i in seq_len(ped_size)) {
    contribute_list[[pedigree$id[i]]] <- contr_mat
    contribute_list[[i]][i,] <- 1
    
    if (i == ped_size) {
      contribute_list[[i]] <- contribute_list[[i]][apply(contribute_list[[i]], 1, function(x) sum(x) != 0),,drop = FALSE] %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = 'id')
      break
    }
    
    for (x in (i + 1):ped_size) {
      
      if (any(is.na(pedigree[x,c('sire_id', 'dam_id')]))) next
      
      #retrieve sire and dam index for focal individual
      sire_index <- pedigree[x, 'sire_id']
      dam_index <- pedigree[x, 'dam_id']
      
      sire_draw <- sire_draw_vec[i]
      dam_draw <- dam_draw_vec[i]
      
      contribute_list[[i]][x,1] <- contribute_list[[i]][sire_index, sire_draw]
      contribute_list[[i]][x,2] <- contribute_list[[i]][dam_index, dam_draw]
    }
    
    contribute_list[[i]] <- contribute_list[[i]][apply(contribute_list[[i]], 1, function(x) !all(x == 0)),,drop = FALSE] %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = 'id')
  }
  
  return(
    contribute_list %>%
      #dplyr::bind_rows(., .id = 'focal_ind') %>%
      dplyr::bind_rows(.id = 'focal_ind') %>%
      dplyr::mutate(geno_sum = .data$sire_contr + .data$dam_contr) %>%
      dplyr::select(.data$focal_ind, .data$id, .data$geno_sum) %>%
      `rownames<-`( NULL )
  )
}





### ALTERNATIVE VERSION OF GENE DROP THAT USES MATRICES ###


# Single locus gene drop on a pedigree (using matrices to record the alleles)
#
#ped pedigree (stored in dataframe) with the following organization: column 1 --> id, column 2 --> sire id, column 3 --> dam. Founder parents should be coded as 0s.
#sims the number of gene drops to perform
#report_progress should the simulation progress be reported?

gene_drop_matrix <- function(ped, sims = 10, report_progress = TRUE) {
  ped_reorder <- reorder_ped(ped, 1, 2, 3)
  ped_reorder_index <- index_pedigree(ped_reorder)
  
  ped_size <- nrow(ped_reorder_index)
  
  init_allele_mat <- matrix(data = 0, nrow = ped_size, ncol = 2)
  
  #allele_mat_sire <- matrix(data = 0, nrow = ped_size, ncol = sims)
  #allele_mat_dam <- matrix(data = 0, nrow = ped_size, ncol = sims)
  
  allele_mat_list <- stats::setNames(replicate(2,
                                               matrix(data = 0,
                                                      nrow = ped_size,
                                                      ncol = sims,
                                                      dimnames = list(ped_reorder[,1,drop=TRUE], NULL)),
                                               simplify = FALSE),
                                     nm = c('sire', 'dam'))
  
  sire_founder_index <- which(ped_reorder_index[,2,drop=TRUE] == 0)
  dam_founder_index <- which(ped_reorder_index[,3,drop=TRUE] == 0)
  
  init_allele_mat[sire_founder_index,1] <- sire_founder_index
  init_allele_mat[dam_founder_index,2] <- -dam_founder_index
  
  #initiate progress bar
  if (report_progress) {
    message('Starting gene drop simulation', ifelse(sims > 1, 's.', '.'))
    prog_bar <- utils::txtProgressBar(min = 0, max = sims, initial = 0, char = "*", style = 3)
  }
  
  for (SIM in seq_len(sims)) {
    
    completed_mat <- single_drop_matrix(ped = ped_reorder_index,
                                        ped_size = ped_size,
                                        allele_mat = init_allele_mat)
    
    for (i in seq_len(2)) allele_mat_list[[i]][,SIM] <- completed_mat[,i]
    #allele_mat_sire[,SIM] <- completed_mat[,1]
    #allele_mat_dam[,SIM] <- completed_mat[,2]
    
    if (report_progress) utils::setTxtProgressBar(prog_bar, SIM)
    
  }
  
  allele_mat_list[['reordered_pedigree']] <- ped_reorder
  allele_mat_list[['indexed_pedigree']] <- ped_reorder_index
  allele_mat_list[['founder_alleles']] <- data.frame(id = c(ped_reorder[sire_founder_index,1,drop=TRUE],
                                                            ped_reorder[dam_founder_index,1,drop=TRUE]),
                                                     allele = c(sire_founder_index,
                                                                -dam_founder_index))
  allele_mat_list[['sim_count']] <- sims
  
  return(allele_mat_list)
}


single_drop_matrix <- function(ped,
                               ped_size,
                               allele_mat) {
  
  fid_draw <- stats::rbinom(n = ped_size, size = 1, prob = 0.5) + 1L
  mid_draw <- stats::rbinom(n = ped_size, size = 1, prob = 0.5) + 1L
  
  for (i in seq_len(ped_size)) {
    #sire
    if (ped[i,2,drop=TRUE] != 0) {
      allele_mat[i,1] <- allele_mat[ped[i,2,drop=TRUE], fid_draw[i]]
    }
    
    #dam
    if (ped[i,3,drop=TRUE] != 0) {
      allele_mat[i,2] <- allele_mat[ped[i,3,drop=TRUE], mid_draw[i]]
    }
  }
  
  return(allele_mat)
}


#########################
### COANCESTRY MATRIX ###
#########################
# Calculate kinship/coancestry matrix or additive relationship matrix from a pedigree
#
#ped the pedigree from which to calculate the kinship matrix
#id_col the name or index of the id column
#sire_col the name or index of the sire column
#dam_col the name or index of the dam column
#type Whether the kinship matrix or additive relationship matrix should be calculated
#
#the kinship or additive relationship matrix
calc_pedmat <- function(ped,
                        id_col = 1,
                        sire_col = 2,
                        dam_col = 3,
                        type = c('kinship', 'add_rel')) {
  #resources:
  #https://doi.org/10.3389/fgene.2021.655638
  #https://github.com/mayoverse/kinship2/blob/master/R/kinship.R
  
  #========================
  #== PREPARATION STEPS ===
  #========================
  
  #On my computer, feeding a tibble to pedigree's orderPed results in the R session
  #being aborted. This doesn't happen if ped is a dataframe
  #perhaps this isn't needed anymore because a custom function is now used...
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)
  
  #reorder pedigree so parents occur before offspring
  #ped_reorder <- ped[order(pedigree::orderPed(ped[,c(id_col, sire_col, dam_col)])),]
  ped_reorder <- reorder_ped(ped,
                             id_col = id_col,
                             sire_col = sire_col,
                             dam_col = dam_col)
  
  id_vec <- ped_reorder[,id_col, drop = TRUE] #id vector
  n_plus_one <- length(id_vec) + 1 #number of individuals in pedigree plus 1
  
  #vectors of parents: vectors hold the index of each parent in the id_vec to
  #facilitate easy retrieval of the parents' kinship values from the kinmat
  #if the parent isn't in id_vec, it is a founder so assign it a value of n_plus_one
  sire_vec <- match(ped_reorder[, sire_col, drop = TRUE], id_vec, nomatch = n_plus_one)
  dam_vec <- match(ped_reorder[, dam_col, drop = TRUE], id_vec, nomatch = n_plus_one)
  
  #empty matrix to store kinship values
  kinmat <- matrix(data = 0,
                   nrow = n_plus_one, ncol = n_plus_one,
                   dimnames = list(c(id_vec, n_plus_one), c(id_vec, n_plus_one))
  )
  
  
  #===========================
  #== KINSHIP CALCULATIONS ===
  #===========================
  for (x in seq_along(id_vec)) {
    kinmat[,x] <- kinmat[x,] <- (kinmat[dam_vec[x],] + kinmat[sire_vec[x],])/2
    kinmat[x,x] <- (1 + kinmat[dam_vec[x], sire_vec[x]])/2
  }
  
  #return the kinship matrix (without the extra founder column and row)
  return(
    switch(type,
           kinship = {kinmat[-(n_plus_one), -(n_plus_one)]},
           add_rel = {2*kinmat[-(n_plus_one), -(n_plus_one)]})
  )
}





# Calculate partial kinship matrix for a given pedigree and focal founder
#
#ped the pedigree from which to calculate the partial kinship matrix
#id_col the name or index of the id column
#sire_col the name or index of the sire column
#dam_col the name or index of the dam column
#founder_val the value in the sire and dam columns that represent founder
#focal_founder the founder for which to calculate the partial kinship matrix
#
#a partial kinship matrix for the focal founder
# @export
#
calc_partial_kinmat <- function(ped,
                                id_col = 1,
                                sire_col = 2,
                                dam_col = 3,
                                founder_val = 0,
                                focal_founder) {
  #resources:
  #https://doi.org/10.3389/fgene.2021.655638
  #https://github.com/mayoverse/kinship2/blob/master/R/kinship.R
  #https://doi.org/10.1016/j.livsci.2006.04.007
  
  #========================
  #== PREPARATION STEPS ===
  #========================
  
  #On my computer, feeding a tibble to pedigree's orderPed results in the R session
  #being aborted. This doesn't happen if ped is a dataframe
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)
  
  #reorder pedigree so parents occur before offspring
  #ped_reorder <- ped[order(pedigree::orderPed(ped[,c(id_col, sire_col, dam_col)])),] #12/13 --> check that this works
  ped_reorder <- reorder_ped(ped,
                             id_col = id_col,
                             sire_col = sire_col,
                             dam_col = dam_col)
  
  founders <- ped_reorder[ped_reorder[,sire_col] == founder_val & ped_reorder[,dam_col] == founder_val, id_col, drop = TRUE]
  non_focal_founders <- founders[founders != focal_founder]
  
  id_vec <- ped_reorder[,id_col] #id vector
  n_plus_one <- length(id_vec) + 1 #number of individuals in pedigree plus 1
  
  #vectors of parents: vectors hold the index of each parent in the id_vec to
  #facilitate easy retrieval of the parents' kinship values from the kinmat
  #if the parent isn't in id_vec, it is a founder so assign it a value of n_plus_one
  sire_vec <- match(ped_reorder[, sire_col], id_vec, nomatch = n_plus_one)
  dam_vec <- match(ped_reorder[, dam_col], id_vec, nomatch = n_plus_one)
  founder_index <- which(id_vec == focal_founder)
  
  #matrix of 0s to store kinship values
  partial_kinmat <- matrix(data = 0,
                           nrow = n_plus_one, ncol = n_plus_one,
                           dimnames = list(c(id_vec, n_plus_one), c(id_vec, n_plus_one))
  )
  
  
  #===================================
  #== PARTIAL KINSHIP CALCULATIONS ===
  #===================================
  for (x in seq_along(id_vec)) {
    if (id_vec[x] %in% non_focal_founders) {
      partial_kinmat[,x] <- partial_kinmat[x,] <- 0
    } else {
      partial_kinmat[,x] <- partial_kinmat[x,] <- (partial_kinmat[dam_vec[x],] + partial_kinmat[sire_vec[x],])/2
      
      if (x == founder_index) {
        partial_kinmat[x,x] <- 0.5
      } else {
        partial_kinmat[x,x] <- partial_kinmat[x,founder_index] + 0.5*partial_kinmat[dam_vec[x], sire_vec[x]]
      }
    }
  }
  
  return(partial_kinmat[-(n_plus_one), -(n_plus_one)])
}






#############################################
### INBREEDING BASED ON COANCESTRY MATRIX ###
#############################################

# Calculate inbreeding from a kinship/coancestry matrix
#
#ped the pedigree from which to calculate inbreeding
#kin_mat the kinship matrix from to which to calculate inbreeding. If this is not provided, the matrix is calculated internally.
#id_col the name or index of the id column
#sire_col the name or index of the sire column
#dam_col the name or index of the dam column
#
#a dataframe with inbreeding values for each individual
# @export
#
inbr_from_kinmat <- function(ped,
                             kin_mat = NULL,
                             id_col = 1,
                             sire_col = 2,
                             dam_col = 3) {
  
  #On my computer, feeding a tibble to orderPed results in the R session being
  #aborted. This doesn't happen if ped is a dataframe
  #MAYBE THIS ISN'T NECESSARY ANYMORE
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)
  
  #if a kinship matrix isn't provided, calculate one from the input pedigree
  if (is.null(kin_mat)) {
    kin_mat <- calc_pedmat(ped = ped,
                           id_col = id_col,
                           sire_col = sire_col,
                           dam_col = dam_col,
                           type = 'kinship')
  }
  
  if (!identical(colnames(kin_mat), rownames(kin_mat)))
    stop('The column and row names for kin_mat must be identical (including identical order)')
  
  ### extract the indices in kin_mat for sire and dam
  #(the dam and sire indices for founders will return NA)
  dam_vec <- match(ped[,dam_col, drop = TRUE], colnames(kin_mat))
  sire_vec <- match(ped[,sire_col, drop = TRUE], colnames(kin_mat))
  
  inbr_df <- data.frame(id = ped[,id_col, drop = TRUE],
                        f_ped = NA)
  
  for (i in seq_len(nrow(inbr_df))) {
    
    if ( any(is.na(c(dam_vec[i], sire_vec[i]))) ) {#check to make sure this works and is correct
      inbr_df[i,'f_ped'] <- 0
    } else {
      inbr_df[i,'f_ped'] <- kin_mat[dam_vec[i], sire_vec[i]]
    }
  }
  
  return(inbr_df)
}



# Calculate the partial inbreeding value for each individual derived from each founder
#
#ped the pedigree from which to calculate the partial kinship matrix
#id_col the name or index of the id column
#sire_col the name or index of the sire column
#dam_col the name or index of the dam column
#founder_val the value in the sire and dam columns that represent founder
#
#a list that contains a list of partial coancestry matrices and a dataframe of partial inbreeding values
# @export
#
partial_inbreeding <- function(ped,
                               id_col = 1,
                               sire_col = 2,
                               dam_col = 3,
                               founder_val = 0) {
  
  #On my computer, feeding a tibble to pedigree's orderPed results in the R session
  #being aborted. This doesn't happen if ped is a dataframe
  #THIS MIGHT NOT BE NECESSARY ANYMORE: 12/13/23
  if (!identical("data.frame", class(ped))) ped <- as.data.frame(ped)
  
  founders <- ped[ped[,sire_col] == founder_val & ped[,dam_col] == founder_val, id_col, drop = TRUE]
  
  partial_coanc_list <- lapply(stats::setNames(nm = founders), function(x) {
    calc_partial_kinmat(ped = ped,
                        id_col = id_col,
                        sire_col = sire_col,
                        dam_col = dam_col,
                        focal_founder = x)
  })
  
  partial_inbr_list <- lapply(partial_coanc_list, function(COANC_MAT, ped) {
    inbr_from_kinmat(ped = ped,
                     kin_mat = COANC_MAT,
                     id_col = id_col,
                     sire_col = sire_col,
                     dam_col = dam_col)
  }, ped = ped)
  
  #return list containing 2 elements:
  #(1) list of partial coancestry matrices
  #(2) dataframe of partial inbreeding values
  return(
    list(partial_coanc_list = partial_coanc_list,
         partial_inbr = dplyr::bind_rows(partial_inbr_list, .id = "focal_founder"))
  )
  
}


# Calculate the altered genetic relatedness matrix of Hunter et al. (2015)
#
#ped a pedigree
#cohort the cohort (containing one or more individuals) used to create the truncated pedigree
#
#a matrix containing the additive genetic relatedness value between all individuals calculated based on the truncated pedigree
# @export
#
calc_Astar <- function(ped, cohort) {
  ped[ped[,1,drop=TRUE] %in% cohort,c(2, 3)] <- '0'
  return(calc_pedmat(ped, type = "add_rel"))
}


# Calculate the expected genetic contribution of a set of individuals to another set based on the Hunter et al. (2019) method
#
#ped pedigree dataframe
#contributors the cohort for which to calculate genetic contribution
#recipients the "recipient" group for which to calculate the genetic contribution of the individual(s) included in contributors
#focal_contemporaries the individuals alive at the same time as the contributor cohort. This is only necessary if standardize is set to two.
#standardize the method used to standardize the genetic contribution values.
# This includes:
#   none: unstandardized value is returned
#   one: standardize the contribution value by the size of the recipients cohort
#   two: standardize the contribution value by the expected contribution value of all individuals alive at the same time as the contributors cohort. For this standardization, the ids of individuals in the population that exists at the same time as the contributors cohort must be input using the focal_contemporaries argument
#
#the expected genetic contribution of the contributors cohort to the recipients cohort
# @export
#
calc_gen_contr <- function(ped,
                           contributors,
                           recipients,
                           focal_contemporaries = NULL,
                           standardize = c('none', 'one', 'two')) {
  
  if (standardize == 'two') {
    if (is.null(focal_contemporaries))
      stop("if standardize is set to 'two', contemporaries of contributors need to be input via the focal_contemporaries argument")
    
    Astar_mat_contemps <- calc_Astar(ped, focal_contemporaries)
  }
  
  Astar_mat <- calc_Astar(ped = ped, cohort = contributors)
  
  unstd_contr <- sum(Astar_mat[rownames(Astar_mat) %in% recipients, colnames(Astar_mat) %in% contributors])
  
  return(
    switch(standardize,
           'none' = unstd_contr,
           'one' = unstd_contr/length(recipients),
           'two' = {
             unstd_contr/sum(Astar_mat_contemps[rownames(Astar_mat_contemps) %in% recipients, colnames(Astar_mat_contemps) %in% focal_contemporaries])
           }
    )
  )
}





################################################
### FUNCTIONS BASED ON GENE DROP SIMULATIONS ###
################################################
# Calculation of inbreeding from gene dropping
#
#gdrop_output output from the simple_gene_drop function (specifically, the geno_origin dataframe)
#
#a dataframe containing the inbreeding values
# @export
#
fped_gdrop <- function(gdrop_output) {
  return(
    gdrop_output %>%
      dplyr::group_by(.data$id) %>%
      dplyr::mutate(f_status = .data$sire_geno_origin == .data$dam_geno_origin) %>%
      dplyr::summarize(f_ped = sum(.data$f_status)/dplyr::n())
  )
}


# Calculation of partial inbreeding from gene dropping
#
#gdrop_output output from the simple_gene_drop function (specifically, the geno_origin dataframe)
#
#a dataframe containing each individual's inbreeding values contributed by each founder
# @export
#
partial_founder_fped_gdrop <- function(gdrop_output) {
  return(
    gdrop_output %>%
      dplyr::group_by(.data$id) %>%
      dplyr::mutate(f_status = .data$sire_geno_origin == .data$dam_geno_origin) %>%
      dplyr::mutate(fped_founder_id = ifelse(.data$f_status, .data$sire_geno_origin, NA)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(.data$id, .data$fped_founder_id) %>%
      dplyr::summarize(#id = first(id),
        #founder_id = first(fped_founder_id),
        founder_count = sum(.data$f_status),
        .groups = 'drop') %>%
      dplyr::group_by(.data$id) %>%
      dplyr::mutate(fped_partial_ancestor = .data$founder_count/sum(.data$founder_count)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(.data$fped_founder_id))
  )
}


# Calculate coancestry/kinship from gene dropping (the gene_drop_matrix function). NOTE: THIS FUNCTION HAS NOT BEEN TESTED
#
#gdrop_output output from the simple_gene_drop function (specifically, the geno_origin dataframe)
#samples the samples for which you want to calculate kinship values. If "all" is specified, all individuals in the pedigree will be included in the calculations.
#show_progress specifies whether the function should report the progress of the calculations
#
#a dataframe with all pairwise combos of individuals and their associated kinship value
# @export
#
kinship_gdrop <- function(gdrop_output, samples = 'all', show_progress = TRUE) {
  
  if (samples == 'all') samples <- unique(gdrop_output$id)
  
  sample_combo_df <- as.data.frame(t(utils::combn(samples, 2))) %>%
    dplyr::rename(samp1 = "V1",
                  samp2 = "V2") %>%
    dplyr::mutate(kinship = NA)
  
  show_progress <- isTRUE(show_progress)
  n_comp <- nrow(sample_combo_df)
  
  if (show_progress) p_bar <- utils::txtProgressBar(0, n_comp, style = 3, char = '*')
  
  for (i in seq_len(n_comp)) {
    sample_combo_df[i,3] <- gdrop_output %>%
      dplyr::filter(.data$id %in% unlist(sample_combo_df[i,-3])) %>%
      #pivot_longer(cols = ends_with('origin'), values_to = 'allele_id') %>%
      dplyr::group_by(.data$sim) %>%
      #arrange(id) %>%
      dplyr::summarize(kin_val = 0.25 * ( (.data$sire_geno_origin[1] == .data$sire_geno_origin[2]) +
                                            (.data$sire_geno_origin[1] == .data$dam_geno_origin[2]) +
                                            (.data$dam_geno_origin[1] == .data$sire_geno_origin[2]) +
                                            (.data$dam_geno_origin[1] == .data$dam_geno_origin[2]) )) %>%
      dplyr::pull(.data$kin_val) %>%
      mean()
    
    if (show_progress) utils::setTxtProgressBar(p_bar, i)
  }
  
  return(sample_combo_df)
}



# Calculate inbreeding from gene dropping (the gene_drop_matrix function)
#
#gdrop_mat_output gene dropping output from the gene_drop_matrix function
#
#a dataframe containing the estimated inbreeding for each individual
# @export
#
fped_gdrop_mat <- function(gdrop_mat_output) {
  
  return(
    data.frame(id = gdrop_mat_output$reordered_pedigree[,1,drop=TRUE],
               f_ped = apply(gdrop_mat_output$sire == gdrop_mat_output$dam, 1, sum)/ncol(gdrop_mat_output$sire))
  )
  
}


# Calculate the partial founder inbreeding values from gene dropping (the gene_drop_matrix function)
#
#gdrop_mat_output gene dropping output from the gene_drop_matrix function
#
#a dataframe containing the partial founder inbreeding values (partial_founder_fped) and the proportion of inbreeding that can be attributed to each founder (fped_prop)
# @export
#
partial_founder_fped_gdrop_mat <- function(gdrop_mat_output) {
  
  #get the indices of the individuals/sim combos that represent IBD
  inbreeding_indices <- which(gdrop_mat_output$sire == gdrop_mat_output$dam,
                              arr.ind = TRUE)
  
  ### NOTES: ###
  #STEP 1: dataframe where each row represents an IBD instance where ID is the indiv
  #where the IBD event occurred and allele_origin is the founder where the allele
  #originated. Because each allele copy from a founder is coded as x and -x, getting
  #the origin of the allele is as simple as taking the absolute value of allele copy
  #These IBD events are extracted via the inbreeding_indices matrix created above.
  #The sire matrix is arbitrarily used to obtain the allele origin for the IBD
  #alleles. But since we are working with instances of IBD, it would equivalent
  #to use the dam matrix for this purpose
  
  #STEP 2: sum up the number of times IBD occurs for each id and allele origin
  #combo (e.g., indiv x has an instance of IBD due to an allele provided by
  #founder y)
  
  #STEP 3: divide the counts by the total count to get the proportion of inbreeding
  #that can be attributed to each founder (fped_prop). Divide the counts by the
  #total number of simulations (extracted as column count from sire matrix) to
  #get the partial founder inbreeding (the probability that an individual will
  #be IBD at a locus due to a particular founder). This quantity is stored in
  #partial_founder_fped.
  return(
    #Step 1
    data.frame(id = rownames(inbreeding_indices),
               #sim = as.character(inbreeding_indices[,2,drop=TRUE]),
               allele_origin = as.character(abs(gdrop_mat_output$sire)[inbreeding_indices])) %>%
      dplyr::group_by(.data$id, .data$allele_origin) %>%
      dplyr::summarize(count = dplyr::n(), #Step 2
                       .groups = 'drop') %>%
      dplyr::group_by(.data$id) %>%
      dplyr::mutate(partial_founder_fped = .data$count/ncol(gdrop_mat_output$sire), #Step 3
                    fped_prop = .data$count/sum(.data$count)) %>% #Step 3 (ctd)
      dplyr::ungroup()
  )
}




# Calculate coancestry/kinship from gene dropping (the gene_drop_matrix function)
#
#gdrop_mat_output gene dropping output from the gene_drop_matrix function
#samples the samples for which you want to calculate kinship values. If "all" is specified, all individuals in the pedigree will be included in the calculations.
#show_progress specifies whether the function should report the progress of the calculations
#output should the kinship values be output as a dataframe or a matrix?
#
#a dataframe with all pairwise combos of individuals and their associated kinship value
# @export
#
kinship_gdrop_mat <- function(gdrop_mat_output,
                              samples = 'all',
                              show_progress = TRUE,
                              output = c('matrix', 'dataframe')) {
  
  output <- match.arg(output, several.ok = FALSE)
  
  if (identical(samples,'all')) {
    samples <- gdrop_mat_output$indexed_pedigree[,1,drop=TRUE]
  } else {
    samples <- which(gdrop_mat_output$reordered_pedigree[,1,drop=TRUE] %in% samples)
    if (length(samples) == 0) stop('None of input samples exist in the pedigree.')
  }
  
  #sample_combo_df <- as.data.frame(t(utils::combn(samples, 2))) %>%
  sample_combo_df <- rbind(t(utils::combn(samples, 2)),
                           cbind(samples, samples)) %>%
    as.data.frame() %>%
    dplyr::rename(samp1 = 1,
                  samp2 = 2) %>%
    #dplyr::rename(samp1 = "V1",
    #              samp2 = "V2") %>%
    dplyr::mutate(kinship = NA)
  
  show_progress <- isTRUE(show_progress)
  n_comp <- nrow(sample_combo_df)
  #sim_count <- ncol(gdrop_mat_output$sire)
  
  if (show_progress) p_bar <- utils::txtProgressBar(0, n_comp, style = 3, char = '*')
  
  for (i in seq_len(n_comp)) {
    
    #the IBD counts across all four inter-individual allele copy comparisons
    #this is processed outside of the for loop to calculate the kinship value
    sample_combo_df[i,'kinship'] <- (sum(gdrop_mat_output$sire[sample_combo_df$samp1[i],] == gdrop_mat_output$sire[sample_combo_df$samp2[i],]) +
                                       sum(gdrop_mat_output$sire[sample_combo_df$samp1[i],] == gdrop_mat_output$dam[sample_combo_df$samp2[i],]) +
                                       sum(gdrop_mat_output$dam[sample_combo_df$samp1[i],] == gdrop_mat_output$sire[sample_combo_df$samp2[i],]) +
                                       sum(gdrop_mat_output$dam[sample_combo_df$samp1[i],] == gdrop_mat_output$dam[sample_combo_df$samp2[i],]))
    
    if (show_progress) utils::setTxtProgressBar(p_bar, i)
  }
  
  sample_combo_df[,'samp1_name'] <- gdrop_mat_output$reordered_pedigree$id[sample_combo_df$samp1]
  sample_combo_df[,'samp2_name'] <- gdrop_mat_output$reordered_pedigree$id[sample_combo_df$samp2]
  
  #calculate kinship from the IBD counts by multiplying the summed IBD count by 0.25 (because 4 different
  #allele copy comparisons are being made) and then divide by the number of simulations used in gene dropping
  sample_combo_df[,'kinship'] <- sample_combo_df[,'kinship'] * 0.25 * (1/ncol(gdrop_mat_output$sire))
  
  if (output == 'dataframe') {
    return(sample_combo_df)
  } else {
    #output as a matrix
    #currently, self kinship isn't included (i.e., the diagonals are NA)
    kinship_mat <- tidyr::pivot_wider(sample_combo_df[,c('samp1_name', 'samp2_name', 'kinship')],
                                      names_from = 'samp2_name',
                                      values_from = 'kinship') %>%
      tibble::column_to_rownames(var = 'samp1_name') %>%
      as.matrix()
    
    kinship_mat_reorder <- kinship_mat[,match(rownames(kinship_mat), colnames(kinship_mat))] #reorder columns to match row order
    kinship_mat_reorder[lower.tri(kinship_mat_reorder)] <- t(kinship_mat_reorder)[lower.tri(kinship_mat_reorder)] #copy upper triangle to lower triangle
    
    return(kinship_mat_reorder)
  }
}



# Quantify founder ancestry from gene dropping (the gene_drop_matrix function)
#
#gdrop_mat_output gene dropping output from the gene_drop_matrix function
#founder_group_info the group membership of each founder. The dataframe should have one column named "id" that includes the founder ID and another column named "group" that indicates group membership
#
#a dataframe containing the ancestry proportion from each group
# @export
#
quantify_ancestry <- function(gdrop_mat_output,
                              founder_group_info) {
  
  if (!identical(rownames(gdrop_mat_output$sire), rownames(gdrop_mat_output$dam)))
    stop('The gene drop matrices must have the same individuals in the same order')
  
  
  ### Switched to inner_join from left_join
  #allele_group_map <- dplyr::left_join(gdrop_mat_output$founder_alleles,
  #                                     founder_group_info, by = 'id')
  
  allele_group_map <- dplyr::inner_join(gdrop_mat_output$founder_alleles,
                                        founder_group_info, by = 'id')
  
  #STEP 1:
  #-combine sire and dam gene drop matrices (cbind func)
  #-for each row (representing each individual), tally up each allele (table func)
  #-for the list of named vectors outputted from table, convert each vector to a dataframe
  # with a column of "values" (allele counts) and a column indicating the allele
  # labelled "ind" (stack func). Each dataframe corresponds to the alleles for
  #each individual
  
  #STEP 2: combine the dataframes into a single dataframe (bind_rows command)
  #        The newly added "ind" column indicates the individual that each
  #        allele count corresponds to
  
  #STEP 3: filter the considered alleles down to those originating from the
  #        groups included in the founder_group_info input
  
  #STEP 4: convert allele to integer (stack func outputs it as a factor)
  
  #STEP 5: join the allele count df with information about the group from which
  #        each allele is sourced. This information is provided by the user
  #        via the founder_group_info argument. The dataframes are joined along
  #        the allele column
  
  #STEP 6: count up the number of alleles belonging to each individual/group subset
  #        For example, the number of alleles that individual A inherited from
  #        group 1
  
  #STEP 7: standardize the allele counts by dividing by twice the number of simulations
  #        used in gene dropping. The doubling of the sim count is because we are
  #        working with diploid organisms in these simulations
  return(
    allele_group_membership_df <- lapply( #STEP 1
      apply(cbind(gdrop_mat_output$sire, gdrop_mat_output$dam), 1, table, simplify = FALSE), #STEP 1
      utils::stack) %>% #STEP 1
      dplyr::bind_rows(.id = 'id') %>% #STEP 1
      dplyr::rename("allele" = "ind") %>% #STEP 2
      dplyr::filter(.data$allele %in% unique(allele_group_map$allele)) %>% #STEP 3
      dplyr::mutate("allele" = as.integer(levels(.data$allele)[.data$allele])) %>% #STEP 4
      dplyr::left_join(allele_group_map[,c('group', 'allele')], #STEP 5
                       by = 'allele') %>% #STEP 5
      dplyr::group_by(.data$id, .data$group) %>% #STEP 6
      dplyr::summarize(allele_count = sum(.data$values), #STEP 6
                       .groups = 'drop') %>% #STEP 6
      dplyr::mutate(anc_prop = .data$allele_count/(gdrop_mat_output$sim_count*2) ) #STEP 7
  )
}
