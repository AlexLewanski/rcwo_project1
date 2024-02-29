############################
############################
### RCW CUSTOM FUNCTIONS ###
############################
############################

#######################
### DATA PROCESSING ###
#######################
add_intervening_years <- function(data, id_col, year_col) {
  
  #split the dataframe into list of dataframes by id_col, and create a dataframe
  #that includes all the years that are not found between the first and last
  #recorded years (i.e. years where the individual was alive but not found during
  #the census)
  missing_years_list <- lapply(split(data, data[, id_col]), function(x) {
    
    year_range <- range(x[, year_col])
    year_vec <- year_range[1]:year_range[2]
    
    #explicitly drop to vector because if this a tibble, then it will remain a
    #dataframe/tibble and all elements in missing_years_eval will be FALSE
    missing_years_eval <- year_vec %in% x[, year_col, drop = TRUE]
    if (all(missing_years_eval)) return(data.frame())
    new_df <- data.frame(id_col = x[1, id_col],
                         year_col = year_vec[!missing_years_eval])
    colnames(new_df) <- c(id_col, year_col)
    
    return(new_df)
  })
  
  #remove empty entries of the list, which results when an individual is recorded in all the years
  missing_years_list_empty_eval <- sapply(missing_years_list, function(x) nrow(x) != 0)
  if (any(missing_years_list_empty_eval)) {
    missing_years_df <- bind_rows(missing_years_list[missing_years_list_empty_eval])
  } else {
    #if no missing years exist across all indivs, return input dataframe
    return(data)
  }
  
  ### add empty columns for missing columns in the missing years dataframe ###
  column_names <- colnames(data) #columns in original dataframe
  missing_names <- column_names[!column_names %in% colnames(missing_years_df)] #missing columns
  for (i in missing_names) missing_years_df[, i] <- NA #add missing columns
  missing_years_df[,order(colnames(missing_years_df), column_names)] #reorder new df to match original
  
  data$row_origin <- 'observed'
  missing_years_df$row_origin <- 'deduced'
  
  #return dataframe with original data and the dataframes with the missing years
  return(bind_rows(data, missing_years_df)) 
}




ped_add_dummy_parents <- function(ped, id, fid, mid, sex, nest_id, founder_parent_vals) {
  
  colname_vec <- setNames(c('id', 'MaleID', 'FemaleID', 'sex', 'NatalNest'),
                          nm = c(id, fid, mid, sex, nest_id))
  
  for (i in names(colname_vec)) colnames(ped)[colnames(ped) == i] <- colname_vec[i]
  
  single_parent_nest_df <- ped %>%
    mutate(row_index = 1:n()) %>% 
    rowwise() %>% 
    filter(sum(c(MaleID, FemaleID) %in% founder_parent_vals) == 1)
  #filter(sum(is.na(c(MaleID, FemaleID))) == 1)
  
  parent_dummy_id_df <- single_parent_nest_df %>% 
    group_by(NatalNest) %>% 
    slice_head(n = 1) %>% 
    ungroup() %>% 
    mutate(missing_parent = if_else(MaleID %in% founder_parent_vals, 'MaleID', 'FemaleID')) %>% 
    group_split(missing_parent) %>% 
    map_df( ~ .x %>% 
              mutate(dummy_parent_id = paste(missing_parent, 1:n(), sep = "_") )) %>% 
    select(NatalNest, dummy_parent_id)
  
  parental_replacement_id_df <- left_join(single_parent_nest_df, 
                                          parent_dummy_id_df, 
                                          by = 'NatalNest') %>% 
    rowwise() %>% 
    mutate(column_index = if_else(FemaleID %in% founder_parent_vals, 
                                  which(colnames(single_parent_nest_df) == 'FemaleID'), 
                                  which(colnames(single_parent_nest_df) == 'MaleID'))) %>% 
    select(row_index, column_index, dummy_parent_id)
  
  ped[as.matrix(parental_replacement_id_df[,c('row_index', 'column_index')])] <- parental_replacement_id_df$dummy_parent_id
  
  for (i in colname_vec) colnames(ped)[colnames(ped) == i] <- names(colname_vec)[colname_vec == i]
  
  dummy_ids <- unique(parental_replacement_id_df$dummy_parent_id)
  dummy_info_df = data.frame(dummy_ids,
                             ifelse(dummy_ids %in% ped[,fid,drop = TRUE], 1, 2))
  colnames(dummy_info_df) <- c(id, sex)
  
  dummy_info_df[, fid] <- founder_parent_vals[1]
  dummy_info_df[, mid] <- founder_parent_vals[1]
  for (COLNAME in colnames(ped)[!colnames(ped) %in% colnames(dummy_info_df)]) dummy_info_df[, COLNAME] <- NA
  
  return(list(ped = bind_rows(ped, dummy_info_df[match(colnames(ped), colnames(dummy_info_df))]),
              dummy_ids = dummy_info_df[, c(id, sex)])
  )
}



#####################
### VISUALIZATION ###
#####################
circle_radius_area_prop <- function(max_r = 10, val_vec, max_val = 'observed') {
  
  if (max_val == 'observed') 
    max_val <- max(val_vec)
  #1
  max_area <- pi*(max_r^2)
  
  #2
  area_vec <- (max_area*val_vec)/max_val
  
  #3
  return(sqrt(area_vec/pi))
}


circle_origins <- function(radius_vec, stack_order = 'observed', circle_space = 1) {
  if (identical(stack_order, 'observed')) stack_order <- names(radius_vec)
  
  radius_vec_reorder <- radius_vec[match(stack_order, names(radius_vec))]
  
  return(cumsum((radius_vec_reorder + c(0, radius_vec_reorder[-length(radius_vec_reorder)]))) + cumsum(c(0, rep(circle_space, length(radius_vec_reorder) - 1))) - radius_vec_reorder[1])
}



generate_circle_coords <- function(radius_vec, origins_vec, coord_count = 100) {
  if (!identical(names(radius_vec), names(origins_vec))) radius_vec <- radius_vec[match(names(origins_vec), names(radius_vec))]
  
  coord_list <- list()
  
  for (i in names(origins_vec)) {
    coord_list[[i]] <- data.frame(group = i,
                                  x = 0 + radius_vec[i]*cos((0:vertices_count)*2*pi/vertices_count),
                                  y = origins_vec[i] + radius_vec[i]*sin((0:vertices_count)*2*pi/vertices_count))
  }
  
  return(do.call(rbind, coord_list))
}


