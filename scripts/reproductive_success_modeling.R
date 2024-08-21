##############################################################################################
### SCRIPT NAME: reproductive_success_modeling.R
### PURPOSE: identify predictors of RCW reproductive success and evaluate the relationships
###          between translocation-related variables and reproductive success
### PRODUCTS:
###  Translocation ancestry models:
###     lrs_mod_var_list1_named.rds: model formulas for the translocation ancestry mods
###     repro_success_info.rds: processed data for running the translocation ancestry mods
###     rcw_addrel.rds: additive genetic relationship matrix for the RCWs
###     transloc_anc_repro_mod_multipan.png (main text): lifetime repro success vs. transloc
###                                                      ancestry plot
###     mean_group_size_repro_mod_multipan.png (supplement): lifetime repro success vs. mean 
###                                                          group size plot
###     first_year_scaled_repro_mod_multipan.png (supplement): lifetime repro success vs. first 
###                                                            calendar breeding year plot
###     anc_count_repro_mod_multipan.png (supplement): lifetime repro success vs. indiv. ancestry 
###                                                    group count plot
###     fped_repro_mod_multipan.png (supplement): lifetime repro success vs. indiv. inbreeding plot
###     lrs_mod_fit_table.txt: table of fit info for the translocation ancestry mods 
###     repro_success_mods_params_table.txt: parameter estimates for the top translocation ancestry
###                                          mod and auxiliary models
###     repro_success_pp_check_multipan.png: posterior predictive check for the top translocation 
###                                          ancestry mod and auxiliary models
###
###  Translocation models:
###     transloc_mod_plot_info.RDS: list of results and info from the top translocation mod for 
###                                 plotting it in figure 1 of the main text
###     transloc_mod_pp_check.png: posterior predictive check for the top translocation mod
###     lrs_mod_transloceffects_fit_table.txt: table of fit info for the translocation mods 
###     repro_success_transloc_effects_mods_params_table.txt: parameter estimates for the top 
###                                                           translocation mod
##############################################################################################


#####################
### SCRIPT SET-UP ###
#####################

### Loading libraries ###
library(here)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(brms)
library(tidybayes)
library(kableExtra)
library(modelr)
library(see)
library(gridExtra)
library(grid)
#library(glmmTMB)
#library(GMMAT)

### OPTIONS ###
save_hpcc_info <- FALSE
#run_all_mods <- FALSE


### Custom functions ###
source(here('scripts', 'rcw_project_custom_functions.R'))
source(here('scripts', 'ped_functions.R'))
#devtools::load_all('/Users/alexlewanski/Documents/r_packages/pedutils')


### Loading data ###
file_names <- list.files(here('data', 'feb2024_databasemarch2023_processed'))
rcw_processed_list <- lapply(setNames(file_names[grepl('csv$', file_names)], nm = gsub(pattern = '\\.csv', '', file_names[grepl('csv$', file_names)])), function(x) {
  read.csv(here('data', 'feb2024_databasemarch2023_processed', x))
})


#repro_success_df <- read.csv(here('data', 'feb2024_databasemarch2023_processed', 'repro_success_dataset_processed.csv'))
#ancestry_df <- read.csv(here('data', 'feb2024_databasemarch2023_processed', 'ancestry_dataset_processed.csv'))

results_names <- setNames(nm = c("rcws_founder_info",
                                 "rcws_inbr", 
                                 "transloc_info_color", 
                                 "rcw_inbr_merge", 
                                 "rcws_ancestry_info",
                                 "rcws_ancestry_source_by_year",
                                 "rcw_partial_founder_summed_processed",
                                 "rcw_contr_info_df",
                                 "fped_prop_by_group",
                                 "inbr_founder_color")
)

results_list <- lapply(results_names, function(x) read.csv(here('results', paste0(x, '.csv'))))


#census_processed <- read.csv(here('data', 'feb2024_databasemarch2023_processed', 'census_processed.csv'))
pop_info_list <- lapply(setNames(nm = paste0('Scenario', 1:4)), function(X, DF) {
  DF[DF[,X,drop = TRUE],!grepl(pattern = 'Scenario', colnames(DF))]
}, DF = rcw_processed_list$census_processed)
pop_dat <- pop_info_list$Scenario3


### NOTES ###
#lifetime reproductive success based on fledgling

#censor variable: number of breeding opportunities/whether or not an individual is still alive
#inbreeding
#translocation ancestry
#year of first breeding
#sex
#average(group size)
#

#first breeding year
#number of opportunities that an individual had to breed (censored data)
#translocation ancestry proportion 
#inbreeding?
#Sex
#average(group size)?
#variability in group size?


### PLOTTING FUNCTIONS ###
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}


#######################
### DATA PROCESSING ###
#######################

### Kinship/additive genetic relatedness matrix ###
rcw_kinmat <- calc_pedmat(rcw_processed_list$ped_processed, type = 'kinship')
rcw_addrel <- calc_pedmat(rcw_processed_list$ped_processed, type = 'add_rel')

### Breeder dataframe for reproductive success modeling ###
#step 1: 
repro_success_df_sexpivot <- rcw_processed_list$repro_success_dataset_processed %>% 
  rename(fid_transloc_anc = transloc_anc_fid,
         mid_transloc_anc = transloc_anc_mid) %>% 
  #select(fid_dummy, mid_dummy, fid_compl_index, mid_compl_index) %>% 
  tidyr::pivot_longer(cols = starts_with('mid') | starts_with('fid'),
                      names_to = c("sex", ".value"),
                      names_pattern = "([A-Za-z]+)_(.+)")


#identify all the individuals with any missing informatioon on flegldin
repro_success_df_sexpivot_missinfo <- repro_success_df_sexpivot %>% 
  filter(is.na(dummy) |
           is.na(fldg_num) |
           (fldg_num > hatch_num)  |
           is.na(group_size)
         ) %>% 
  pull(dummy) %>% 
  unique()


#nests with 0 fledged offspring recorded in the nesting table but where there is
#at least one individual that is recorded coming from that nest
nests_with_offspring_in_pop <- rcw_processed_list$rcws %>% 
  filter(RCWid %in% unique(pop_dat$RCWid)) %>% 
  filter(!is.na(NatalNest)) %>% 
  pull(NatalNest) %>% 
  unique()

problem_nests <- rcw_processed_list$nests_processed %>% 
  filter(FldgNum == 0) %>% 
  filter(NatalNest_Temp %in% nests_with_offspring_in_pop) %>%
  pull(NatalNest_Temp) %>% 
  unique()

problem_ids <- rcw_processed_list$rcws %>% 
  filter(NatalNest %in% problem_nests) %>% 
  pull(RCWid)

processed_anc_info <- results_list$rcws_ancestry_info %>% 
  select(id, group, anc_prop) %>% 
  pivot_wider(names_from = group,
              values_from = anc_prop,
              values_fill = 0) %>% 
  rowwise() %>% 
  mutate(#compl_simpson = simpson(prop_vec = c(`non-transloc`, ANF, FTB, ONF, FTS, CBJTC, `WSF-CITRUS`), 
        #                         version = c('regular', 'complement')[2]),
         anc_count = sum(c(`non-transloc`, ANF, FTB, ONF, FTS, CBJTC, `WSF-CITRUS`) > 0),
  ) %>% 
  ungroup() %>% 
  rename(dummy = id)

repro_success_df_sexpivot1 <- repro_success_df_sexpivot %>% 
  # rowwise() %>% 
  # filter(#full_grandparent_info == 'yes' &
  #          !is.na(dummy) &
  #          !is.na(fldg_num) &
  #          #(fldg_num <= hatch_num) &
  #          identical(fldg_num <= hatch_num, TRUE) &
  #           !is.na(group_size)) %>% 
  rowwise() %>% 
  filter(identical(fldg_num <= hatch_num, TRUE)) %>% 
  filter(!dummy %in% repro_success_df_sexpivot_missinfo) %>% 
  ungroup() %>% 
  filter(!dummy %in% problem_ids) %>% 
  group_by(dummy) %>% 
  summarize(#transloc = first(transloc),
    full_grandparent_info = first(full_grandparent_info),
    sex = first(sex),
    mean_group_size = mean(group_size),
    transloc_anc = first(transloc_anc),
    fped = first(fped),
    first_year = min(year),
    last_year = max(year),
    life_fldg = sum(fldg_num),
    mean_fldg_per_nest = mean(fldg_num),
    mean_fldg_per_year = sum(fldg_num)/length(unique(year)),
    nesting_years = length(unique(year)),
    total_nests = n(),
    .groups = 'drop') %>% 
  mutate(alive = if_else(dummy %in% unique(pop_dat[pop_dat$year == 2022,]$'RCWid'), 
                         'yes', 'no'),
         birth_during_monitoring = if_else(dummy %in% rcw_processed_list$rcws[rcw_processed_list$rcws$MinAge >= 1994,]$RCWid,
                                           'yes', 'no'),
         transloc = if_else(dummy %in% results_list$transloc_info_color$RCWid,
                                           'yes', 'no'),
         first_year_scaled = first_year - min(first_year)) %>% 
  left_join(.,
            pop_dat %>% 
              group_by(RCWid) %>% 
              summarize(age = max(year) - min(year)) %>% 
              rename(dummy = RCWid),
            by = 'dummy') %>% 
  left_join(., processed_anc_info,
            by = 'dummy') %>% 
  mutate(mean_group_size_z = (mean_group_size - mean(mean_group_size))/sd(mean_group_size),
         transloc_anc_z = (transloc_anc - mean(transloc_anc))/sd(transloc_anc),
         fped_z = (fped - mean(fped))/sd(fped),
         first_year_scaled_z = (first_year_scaled - mean(first_year_scaled))/sd(first_year_scaled),
         #compl_simpson_z = (compl_simpson - mean(compl_simpson))/sd(compl_simpson),
         anc_count_z = (anc_count - mean(anc_count))/sd(anc_count),
         ) %>% 
  left_join(., rcw_processed_list$rcws %>% 
              select(RCWid, MinAge) %>% 
              rename(dummy = RCWid),
            by = 'dummy') %>% 
  mutate(birth_firstbreed_dif = first_year - MinAge,
         birth_lastbreed_dif = last_year - MinAge)



#############################################
### LIFETIME REPRO. SUCCESS MODEL FITTING ###
#############################################

### create lists of predictors for the models ###
lrs_var <- c('mean_group_size', 'first_year_scaled', 'sex', 'transloc_anc', 'anc_count', 'fped')
lrs_mod_var_list <- all_subset_combn(lrs_var)

lrs_mod_var_list1 <- lapply(lrs_mod_var_list, function(x) paste0('life_fldg | cens(censored) ~ ', paste(x, collapse = ' + '), ' + (1|gr(dummy,cov=A))'))


lrs_mod_var_list1_named <- setNames(lrs_mod_var_list1, 
                                    nm = sapply(lrs_mod_var_list, function(x) paste(x, collapse = '_')))


### SAVING OBJECTS FOR RUNNING THINGS ON HPCC ###
if (isTRUE(save_hpcc_info)) {
  saveRDS(lrs_mod_var_list1_named,
          here('data', 'data_hpc', 'lrs_mod_var_list1_named.rds'))
  
  saveRDS(repro_success_df_sexpivot1 %>%
            mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
            as.data.frame(),
          here('data', 'data_hpc', 'repro_success_info.rds'))
  
  saveRDS(rcw_addrel,
          here('data', 'data_hpc', 'rcw_addrel.rds'))
}



mod_fit_life_success <- read.csv(here('results', 'mod_fit_life_success.csv'))
mod_formulas_life_success <- read.delim(here('results', 'life_mod_formulas.txt'), 
                                        header = FALSE)


mod_fit_life_success_processed <- mod_fit_life_success %>% 
  mutate(mod_index = as.integer(sub('.*\\[\\[([0-9]+)\\]\\]', '\\1', X) )) %>% 
  mutate(rank = 1:n()) %>% 
  mutate(predictors = mod_formulas_life_success$V1[mod_index])

top_predictors_lifemod <- mod_fit_life_success_processed %>% 
  filter(elpd_diff >= -4) %>% 
  pull(predictors)

top_predictors_lifemod_vec <- unique(unlist(strsplit(top_predictors_lifemod, 
         split = ' \\+ ')))

top_predictors_lifemod_modstrig <- paste(top_predictors_lifemod_vec, collapse = ' + ')



#######################################################
### TOP TRANSLOCATION ANCESTRY MOD + AUXILIARY MODS ###
#######################################################

beta_prior_poisson <- c(prior_string("normal(0,5)", class = "b", coef = "anc_count"),
            prior_string("normal(0,5)", class = "b", coef = "first_year_scaled"),
            prior_string("normal(0,5)", class = "b", coef = "mean_group_size"),
            prior_string("normal(0,5)", class = "b", coef = "transloc_anc"),
            prior_string("normal(0,5)", class = "b", coef = "sexmid"),
            prior_string("normal(0,5)", class = "b", coef = "fped"),
            prior_string("normal(0,5)", class = "Intercept")
            #prior_string("normal(0,0.001)", class = "sd", coef = 'dummy'),
)


# focal_predictors <- strsplit(x = sub(".*~ ([a-z]+.*) \\+ \\(.*", "\\1", lrs_mod_var_list1_named[[50]]), split = ' \\+ ')[[1]]
# focal_predictors_update <- unname(sapply(focal_predictors, function(x) if (x == 'sex') 'sexmid' else x ))
# 
# c(prior_beta1[prior_beta1$coef %in% focal_predictors_update,],
#   prior_string("normal(0,5)", class = "Intercept"))
# 

beta_prior_hurdle_gamma <- c(prior_string("normal(0,5)", class = "b", coef = "anc_count"),
                        prior_string("normal(0,5)", class = "b", coef = "first_year_scaled"),
                        prior_string("normal(0,5)", class = "b", coef = "mean_group_size"),
                        prior_string("normal(0,5)", class = "b", coef = "transloc_anc"),
                        prior_string("normal(0,5)", class = "b", coef = "sexmid"),
                        prior_string("normal(0,5)", class = "b", coef = "fped"),
                        prior_string("normal(0,5)", class = "Intercept"),
                       prior_string("normal(0,5)", class = "b", dpar = 'hu', coef = "first_year_scaled"),
                       prior_string("normal(0,5)", class = "b", dpar = 'hu', coef = "transloc_anc"),
                       prior_string("normal(0,5)", class = "b", dpar = 'hu', coef = "mean_group_size"),
                       prior_string("normal(0,5)", class = "Intercept", dpar = 'hu')
                        #prior_string("normal(0,0.001)", class = "sd", coef = 'dummy'),
)


repro_mod_list <- list()

for (MOD_OUTCOME in c('life_fldg', 'nesting_years', 'birth_firstbreed_dif', 'birth_lastbreed_dif')) {
  #message('starting model: ', MOD_OUTCOME)
  
  #mod_string <- switch(MOD_OUTCOME,
  #                     birth_firstbreed_dif = {gsub('life_fldg \\| cens\\(censored\\)', MOD_OUTCOME, as.character(mod_fit_life_fldg1[[mod_fit_df$mod_index[mod_fit_df$rank == 1]]]$formula)[1])},
  #                     gsub('life_fldg', MOD_OUTCOME, as.character(mod_fit_life_fldg1[[mod_fit_df$mod_index[mod_fit_df$rank == 1]]]$formula)[1]))
  
  mod_string <- switch(MOD_OUTCOME,
                       birth_firstbreed_dif = {paste(MOD_OUTCOME, ' ~ ', top_predictors_lifemod_modstrig, '+ (1 | gr(dummy, cov = A))')},
                       paste(MOD_OUTCOME, '| cens(censored) ~ ', top_predictors_lifemod_modstrig, '+ (1 | gr(dummy, cov = A))'))
  message('starting model: ', mod_string)
  
  repro_mod_list[[MOD_OUTCOME]] <- brm(
    as.formula(mod_string),
    data = repro_success_df_sexpivot1 %>%
      filter(birth_during_monitoring == 'yes') %>%
      #filter(age != 0) %>%
      mutate(censored = if_else(alive == 'yes', 'right', 'none')),
    data2=list(A = rcw_addrel),
    family = poisson(link = "log"),
    #cov_ranef = list(A = rcw_addrel),
    chains = 4, cores = 2, iter = 10000,
    control = list(adapt_delta = 0.99),
    prior = beta_prior_poisson
  )
}


repro_mod_list[['mean_fldg_per_year']] <- brm(
    bf(
      as.formula(paste('mean_fldg_per_year | cens(censored) ~ ', top_predictors_lifemod_modstrig, '+ (1 | gr(dummy, cov = A))')),
      hu ~ first_year_scaled + transloc_anc + mean_group_size + (1 | gr(dummy, cov = A))
    ),
    data = repro_success_df_sexpivot1 %>%
      filter(birth_during_monitoring == 'yes') %>%
      #filter(age != 0) %>%
      mutate(censored = if_else(alive == 'yes', 'right', 'none')),
    data2=list(A = rcw_addrel),
    family = hurdle_gamma(),
    #cov_ranef = list(A = rcw_addrel),
    chains = 4, cores = 2, iter = 10000,
    control = list(adapt_delta = 0.999),
    prior = beta_prior_hurdle_gamma
  )


plot_mod_list <- setNames(repro_mod_list[names(repro_mod_list) != 'mean_fldg_per_year'],
         nm = c('life_fldg', 'nestyear_mod', 'firstbreed', 'lastbreed'))

# test_without_group <- brm(
#   as.formula("life_fldg | cens(censored) ~ first_year_scaled + sex + transloc_anc + anc_count + fped + (1 | gr(dummy, cov = A))"),
#   data = repro_success_df_sexpivot1 %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     #filter(age != 0) %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"),
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.99)
# )



##########################################
### VIZ OF TRANSLOCATION ANCESTRY MODS ###
##########################################

transloc_mod_color_vec <- setNames(c('#3f0b70', '#04ab80', '#d62828', '#3478e5'),
                                   nm = c("life_fldg", "nestyear_mod", "firstbreed", "lastbreed"))

transloc_mod_color_post <- setNames(c('#c7b6d8', '#9addcc', '#eea9a9', '#adc9f4'),
                                    nm = c("life_fldg", "nestyear_mod", "firstbreed", "lastbreed"))

repro_multipan_list <- list()

pred_vec <- c('transloc_anc', 'mean_group_size', 'first_year_scaled', 'anc_count', 'fped')
repro_mod_title_vec <- setNames(c('Translocation ancestry', 'Mean group size', 'First breeding year', 'Ancestry count', 'Inbreeding'),
                                nm = pred_vec)


response_vec <- setNames(c("life_fldg", "nesting_years", 'birth_firstbreed_dif', 'birth_lastbreed_dif'),
                         nm = c("life_fldg", "nestyear_mod", "firstbreed", "lastbreed"))

repro_plot_list <- list()

for (PRED in pred_vec) {
  message('starting ', PRED)
  mod_pred_list <- lapply(plot_mod_list, 
                          function(MOD_OBJ, dat, PRED) {
                            
                            dat %>% 
                              data_grid(transloc_anc = case_when(PRED == 'transloc_anc' ~ seq_range(transloc_anc, n = 100),
                                                                 TRUE ~ quantile(transloc_anc, prob = 0.5)),
                                        mean_group_size = case_when(PRED == 'mean_group_size' ~ seq_range(mean_group_size, n = 100),
                                                                    TRUE ~ quantile(mean_group_size, prob = 0.5)),
                                        sex = 'mid',
                                        fped =  case_when(PRED == 'fped'~ seq_range(fped, n = 100),
                                                            TRUE ~ quantile(fped, prob = 0.5)),
                                        first_year_scaled = case_when(PRED == 'first_year_scaled'~ seq_range(first_year_scaled, n = 100),
                                                                      TRUE ~ quantile(first_year_scaled, prob = 0.5)),
                                        anc_count = case_when(PRED == 'anc_count' ~ seq_range(anc_count, n = 100),
                                                                TRUE ~ quantile(anc_count, prob = 0.5))
                              ) %>% 
                              add_epred_draws(MOD_OBJ,
                                              ndraws = 500,
                                              re_formula = NA)
                          },
                          PRED = PRED,
                          dat = repro_success_df_sexpivot1 %>% 
                            #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
                            filter(birth_during_monitoring == 'yes')  %>%
                            mutate(censored = if_else(alive == 'yes', 'right', 'none'),
                                   birth_year = MinAge - min(MinAge))
  )
  
  mod_plot_list <- lapply(setNames(nm = names(mod_pred_list)), 
                          function(PRED_OBJ, mod, color_vec, color_post, FOCAL_VAR, dat, response_vec) {
                            
                            var_name <- switch(FOCAL_VAR,
                                               life_fldg = {''})
                            
                            mod[[PRED_OBJ]] %>% 
                              group_by(.data[[FOCAL_VAR]]) %>% 
                              summarize(lower25 = quantile(.epred, prob = 0.025),
                                        lower10 = quantile(.epred, prob = 0.10),
                                        lower20 = quantile(.epred, prob = 0.20),
                                        mean = mean(.epred),
                                        median = quantile(.epred, prob = 0.5),
                                        upper80 = quantile(.epred, prob = 0.80),
                                        upper90 = quantile(.epred, prob = 0.90),
                                        upper975 = quantile(.epred, prob = 0.975),
                                        .groups = 'drop') %>% 
                              ggplot(aes(x = .data[[FOCAL_VAR]], y = mean)) + 
                              geom_point(data = dat,
                                         aes(x = .data[[PRED]], y = .data[[response_vec[PRED_OBJ]]]),
                                         color = '#818181', alpha = 0.5, size = 1.2) +
                              geom_line(data = mod[[PRED_OBJ]],
                                        #aes(x = .data[[FOCAL_VAR]], y = `.epred`, group = .draw),
                                        aes_string(x = FOCAL_VAR, y = '.epred', group = '.draw'),
                                        linewidth = 0.25, 
                                        color = color_post[PRED_OBJ], 
                                        alpha = 0.4) +
                              geom_ribbon(aes(.data[[FOCAL_VAR]],
                                              ymin = lower25, ymax = upper975),
                                          color = color_vec[PRED_OBJ], linetype = 'dashed', linewidth = 1.25,
                                          lineend = 'round', fill = NA,
                                          #fill = '#3f0b70', #'#c7b6d8',
                                          alpha = 0.2) +
                              geom_line(linewidth = 2.5, color = color_vec[PRED_OBJ],
                                        lineend = 'round') +
                              theme_bw() +
                              theme(panel.grid.minor = element_blank(),
                                    panel.grid.major = element_line(linetype = 'dashed', 
                                                                    linewidth = 0.4,
                                                                    color = '#d8d8d8'))
                          }, mod = mod_pred_list, 
                          color_vec = transloc_mod_color_vec, 
                          color_post = transloc_mod_color_post,
                          FOCAL_VAR = PRED,
                          dat = repro_success_df_sexpivot1 %>% 
                            #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
                            filter(birth_during_monitoring == 'yes')  %>%
                            mutate(censored = if_else(alive == 'yes', 'right', 'none'),
                                   birth_year = MinAge - min(MinAge)),
                          response_vec = response_vec
                          )
  
  repro_combined_plot <- cowplot::plot_grid(mod_plot_list$life_fldg +
                                              theme(axis.title.x = element_blank(),
                                                    panel.border = element_blank(),
                                                    axis.line = element_line(linewidth = 0.5)) +
                                              ylab('Lifetime reproductive success'),
                                            mod_plot_list$nestyear_mod +
                                              theme(axis.title.x = element_blank(),
                                                    panel.border = element_blank(),
                                                    axis.line = element_line(linewidth = 0.5)) +
                                              ylab('Total nesting years'),
                                            mod_plot_list$firstbreed +
                                              theme(axis.title.x = element_blank(),
                                                    panel.border = element_blank(),
                                                    axis.line = element_line(linewidth = 0.5)) +
                                              ylab('First nesting age'),
                                            mod_plot_list$lastbreed +
                                              theme(axis.title.x = element_blank(),
                                                    panel.border = element_blank(),
                                                    axis.line = element_line(linewidth = 0.5)) +
                                              ylab('Last nesting age'),
                                            ncol = 2,
                                            align = 'hv',
                                            axis = 'lr',
                                            labels = c('A', 'B', 'C', 'D'),
                                            label_fontface = 'bold')
  
  repro_draws_list <- lapply(setNames(nm = names(plot_mod_list)), 
                             function(MOD, mod_list) {
                               output_list <- list()
                               #message('starting ')
                               beta_draws <- as_draws_df(mod_list[[MOD]])
                               
                               ci_95 <- quantile(beta_draws[[paste0('b_', PRED)]],
                                                 prob = c(0.025, 0.975))
                               output_list[['ci_95']] <- ci_95
                               output_list[['draws']] <- beta_draws %>% 
                                 select(all_of(paste0('b_', PRED))) %>% 
                                 mutate(in_ci = if_else(.data[[paste0('b_', PRED)]] < ci_95[1] | .data[[paste0('b_', PRED)]] > ci_95[2], 'out', 'in'),
                                        group = 'one') %>% 
                                 mutate(mod = case_when(MOD == 'life_fldg' ~ 'Lifetime\nreproductive\nsuccess',
                                                        MOD == 'nestyear_mod' ~ 'Total\nnesting\nyears',
                                                        MOD == 'firstbreed' ~ 'First\nnesting\nage',
                                                        MOD == 'lastbreed' ~ 'Last\nnesting\nage'))
                               
                               return(output_list)
                             }, mod_list = plot_mod_list) 
  
  repro_draws_only <- lapply(repro_draws_list, function(x) x$draws) %>% 
    bind_rows()
  
  
  
  repro_base_violin_plot <- repro_draws_only %>%
    mutate(mod = factor(mod, levels = rev(c("Lifetime\nreproductive\nsuccess", "Total\nnesting\nyears", "First\nnesting\nage", "Last\nnesting\nage")))) %>% 
    ggplot() +
    geom_violin(aes(x = factor(mod),
                    y = .data[[paste0('b_', PRED)]]))
  
  repro_base_violin_plot_build <- ggplot2::ggplot_build(repro_base_violin_plot)$data[[1]]
  repro_base_violin_plot_build <- transform(repro_base_violin_plot_build,
                                            xminv = x - violinwidth * (x - xmin),
                                            xmaxv = x + violinwidth * (xmax - x))
  
  repro_base_violin_plot_build <- rbind(plyr::arrange(transform(repro_base_violin_plot_build, x = xminv), y),
                                        plyr::arrange(transform(repro_base_violin_plot_build, x = xmaxv), -y))
  
  
  
  #1: last nesting age
  #2: first nesting age
  #3: total nesting years
  #4 lifetime reproductive success
  
  repro_base_violin_plot_build <- repro_base_violin_plot_build %>% 
    mutate(fill_group = case_when(group == 1 & y > repro_draws_list$lastbreed$ci_95['2.5%'] & y < repro_draws_list$lastbreed$ci_95['97.5%'] ~ 'inside',
                                  group == 1 & (y < repro_draws_list$lastbreed$ci_95['2.5%'] | y > repro_draws_list$lastbreed$ci_95['97.5%']) ~ 'outside',
                                  group == 2 & y > repro_draws_list$firstbreed$ci_95['2.5%'] & y < repro_draws_list$firstbreed$ci_95['97.5%'] ~ 'inside',
                                  group == 2 & (y < repro_draws_list$firstbreed$ci_95['2.5%'] | y > repro_draws_list$firstbreed$ci_95['97.5%']) ~ 'outside',
                                  group == 3 & y > repro_draws_list$nestyear_mod$ci_95['2.5%'] & y < repro_draws_list$nestyear_mod$ci_95['97.5%'] ~ 'inside',
                                  group == 3 & (y < repro_draws_list$nestyear_mod$ci_95['2.5%'] | y > repro_draws_list$nestyear_mod$ci_95['97.5%']) ~ 'outside',
                                  group == 4 & y > repro_draws_list$life_fldg$ci_95['2.5%'] & y < repro_draws_list$life_fldg$ci_95['97.5%'] ~ 'inside',
                                  group == 4 & (y < repro_draws_list$life_fldg$ci_95['2.5%'] | y > repro_draws_list$life_fldg$ci_95['97.5%']) ~ 'outside'))
  
  
  #  $fill_group <- ifelse(p_build$y >= 0,'Above','Below')
  #This is necessary to ensure that instead of trying to draw
  # 3 polygons, we're telling ggplot to draw six polygons
  repro_base_violin_plot_build$group1 <- with(repro_base_violin_plot_build,
                                              interaction(factor(group),factor(fill_group)))
  
  
  
  violin_post_plot <- repro_base_violin_plot_build %>% 
    filter(fill_group == 'inside') %>% 
    ggplot() + 
    geom_violin(data = repro_draws_only %>%
                  mutate(mod = factor(mod, levels = rev(c("Lifetime\nreproductive\nsuccess", "Total\nnesting\nyears", "First\nnesting\nage", "Last\nnesting\nage")))), 
                aes(x = mod, y = .data[[paste0('b_', PRED)]], fill = mod),
                color = NA, alpha = 0.8) +
    scale_fill_manual(values = unname(rev(transloc_mod_color_post))) +
    ggnewscale::new_scale_fill() +
    geom_polygon(aes(x = x,y = y, group = group1, fill = group1)) +
    scale_fill_manual(values = unname(rev(transloc_mod_color_vec))) +
    theme(legend.position = 'none',
          axis.line.x = element_line(color = 'black', linewidth = 0.5),
          panel.background = element_blank(),
          #panel.grid.major.y = element_line(color = '#d8d8d8', linetype = 'dashed', linewidth = 0.4),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(color = rev(transloc_mod_color_vec), face = 'bold')) +
    coord_flip() +
    ylab('ylab') +
    #scale_x_discrete(name = NULL, sec.axis = sec_axis(~., name = "Y-axis Label on right side")) +
    guides(y = "none") +
    guides(y.sec = guide_axis_label_trans(~.x)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black')
  
  
  repro_multipan_list[[PRED]] <- cowplot::plot_grid(
    gridExtra::grid.arrange(arrangeGrob(repro_combined_plot, 
                             bottom = textGrob(repro_mod_title_vec[PRED], 
                                               gp=gpar(fontface="plain", col="black", fontsize = 12))
    )),
    gridExtra::grid.arrange(arrangeGrob(violin_post_plot + #mod_coefficient_sina +
                               theme(axis.title.x = element_blank()), 
                             bottom = textGrob(paste0(repro_mod_title_vec[PRED], " coefficent"), 
                                               gp=gpar(fontface="plain", col="black", fontsize = 12))
    )),
    ncol = 2, rel_widths = c(0.6, 0.2),
    labels = c('', 'E'), label_fontface = 'bold')
  
}


for (i in names(repro_multipan_list)) {
  if (i == "transloc_anc") {
    output_path <- here('figures', 'main_paper', paste0(i, '_repro_mod_multipan.png') )
  } else {
    output_path <- here('figures', 'supplement', paste0(i, '_repro_mod_multipan.png') )
  }
  
  cowplot::ggsave2(filename = output_path,
                   plot = repro_multipan_list[[i]],
                   width = 9*1.25, height = 5*1.25, bg = 'white')
}



###################################################
### MOD FIT INFO OF TRANSLOCATION ANCESTRY MODS ###
###################################################

lrs_mod_fit_table <- mod_fit_life_success %>% 
  mutate(mod_index = as.integer(sub('.*\\[\\[([0-9]+)\\]\\]', '\\1', X) )) %>% 
  mutate(rank = 1:n()) %>% 
  mutate(predictors = mod_formulas_life_success$V1[mod_index],
         across(where(is.numeric), function(x) round(x, 3))) %>% 
  mutate(predictors = gsub('mean_group_size', 'Mean group size', predictors),
         predictors = gsub('first_year_scaled', 'First breeding year', predictors),
         predictors = gsub('sex', 'Sex', predictors),
         predictors = gsub('anc_count', 'Ancestry count', predictors),
         predictors = gsub('transloc_anc', 'Transloc. ancestry', predictors),
         predictors = gsub('fped', 'Fped', predictors)) %>% 
  #mutate(across(where(is.numeric), function(x) round(x, 3))) %>% 
  select(rank, predictors, elpd_diff, se_diff) %>%
  rename(Rank = rank, 
         `Fixed effects` = predictors, 
         `ELPD diff` = elpd_diff, 
         `SE ELPD diff` = se_diff) %>% 
  kbl('latex', booktabs = TRUE, 
      align = "c",
      linesep = "") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position"),
                position = "center",
                font_size = 7)


cat(lrs_mod_fit_table, 
    file = here('tables', 'lrs_mod_fit_table.txt'), append = FALSE)


table_list_repro_mods <- list()
for (MOD in names(repro_mod_list)) {
  
  table_list_repro_mods[[MOD]] <- summary(repro_mod_list[[MOD]])$fixed %>%
    rbind(., summary(repro_mod_list[[MOD]])$random$dummy) %>% 
    mutate(Response = MOD) %>% 
    rownames_to_column(var = 'Parameter_orig') %>% 
    mutate(Parameter = case_when(
      Parameter_orig == "mean_group_size" ~ "Mean group size",
      Parameter_orig == "first_year_scaled" ~ "First breeding year",
      Parameter_orig == "sexmid" ~ "Sex (female)",
      Parameter_orig == "anc_count" ~ "Ancestry count",
      Parameter_orig == "transloc_anc" ~ "Transloc. ancestry",
      Parameter_orig == "fped" ~ "Fped",
      Parameter_orig == "hu_transloc_anc" ~ "Transloc. ancestry (hurdle)",
      Parameter_orig == "hu_first_year_scaled" ~ "First breeding year (hurdle)",
      Parameter_orig == "hu_mean_group_size" ~ "Mean group size (hurdle)",
      Parameter_orig == "hu_Intercept" ~ "Intercept (hurdle)",
      Parameter_orig == "sd(Intercept)" ~ "sd(Intercept) [rand. effect]",
      Parameter_orig == "sd(hu_Intercept)" ~ "sd(hurdle Intercept) [rand. effect]",
      TRUE ~ Parameter_orig
    )) %>% 
    select(Response, Parameter, Estimate, `l-95% CI`, `u-95% CI`, Rhat, Bulk_ESS, Tail_ESS) %>% 
    rename(`lower 95% CI` = `l-95% CI`,
           `upper 95% CI` = `u-95% CI`,
           `Bulk ESS` = Bulk_ESS, 
           `Tail ESS` = Tail_ESS) 
}

repro_success_mods_params_table <- table_list_repro_mods %>% 
  bind_rows() %>% 
  mutate(across(where(is.numeric), function(x) round(x, 3))) %>% 
  kbl('latex', booktabs = TRUE,
      align = "c", linesep = "",
      escape = TRUE) %>% 
  collapse_rows(columns = 1:2,
                valign = "middle",
                latex_hline = "major") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position"),
                position = "center")


cat(repro_success_mods_params_table, 
    file = here('tables', 'repro_success_mods_params_table.txt'), append = FALSE)


repro_mod_list_ppcheck <- setNames(repro_mod_list,
                                   nm = c('Life repro. success', 'Total nesting years', 'First nesting age', 'Last nesting age', 'Mean fledglings per year'))


repro_success_pp_check_list <- lapply(setNames(nm = names(repro_mod_list_ppcheck)), function(x, MOD, col_vec) {
  
  pp_plot <- pp_check(MOD[[x]], ndraws = 1000) +
    ggtitle(x) +
    theme_classic() +
    theme(legend.position = 'none',
          plot.title = element_text(size = 15))
  
  pp_plot$layers[[1]]$aes_params$colour <- col_vec[[x]]
  return(pp_plot)
  
}, MOD = repro_mod_list_ppcheck,
col_vec = setNames(c('#c7b6d8', '#9addcc', '#eea9a9', '#adc9f4', 'gray'),
                   nm = c("Life repro. success", "Total nesting years", "First nesting age", "Last nesting age", "Mean fledglings per year")))


repro_success_pp_check_multipan <- cowplot::plot_grid(repro_success_pp_check_list$`Life repro. success`,
                                                      repro_success_pp_check_list$`Total nesting years`,
                                                      repro_success_pp_check_list$`Mean fledglings per year`,
                                                      repro_success_pp_check_list$`First nesting age`,
                                                      repro_success_pp_check_list$`Last nesting age`, nrow = 1)

cowplot::ggsave2(filename = here('figures', 'supplement', 'repro_success_pp_check_multipan.png' ),
                 plot = repro_success_pp_check_multipan,
                 width = 9*1.25, height = 2*1.25, bg = 'white')



### Model fit table ###
# lrs_mod_fit_table <- mod_fit_df %>% 
#   select(rank, predictors, elpd_diff, se_diff) %>% 
#   mutate(across(where(is.numeric), function(x) round(x, 3))) %>% 
#   #mutate(elpd_diff = round(elpd_diff, 2),
#   #       se_diff = round(se_diff, 2)) %>% 
#   rename(Rank = rank, 
#          `Fixed effects` = predictors, 
#          `ELPD diff` = elpd_diff, 
#          `SE ELPD diff` = se_diff) %>% 
#   kbl('latex', booktabs = TRUE, align = "c") %>% 
#   kable_styling(latex_options = c("scale_down", "hold_position"),
#                 position = "center")

# mod_fit_life_success %>% 
#   mutate(mod_index = as.integer(sub('.*\\[\\[([0-9]+)\\]\\]', '\\1', X) )) %>% 
#   mutate(rank = 1:n()) %>% 
#   mutate(predictors = mod_formulas_life_success$V1[mod_index],
#          across(where(is.numeric), function(x) round(x, 3))) %>% 
#   mutate(predictors = gsub('mean_group_size', 'Mean group size', predictors),
#          predictors = gsub('first_year_scaled', 'First breeding year', predictors),
#          predictors = gsub('sex', 'Sex', predictors),
#          predictors = gsub('anc_count', 'Ancestry count', predictors),
#          predictors = gsub('transloc_anc', 'Transloc. ancestry', predictors),
#          predictors = gsub('fped', 'Fped', predictors)) %>% 
#   #mutate(across(where(is.numeric), function(x) round(x, 3))) %>% 
#   select(rank, predictors, elpd_diff, se_diff) %>%
#   rename(Rank = rank, 
#          `Fixed effects` = predictors, 
#          `ELPD diff` = elpd_diff, 
#          `SE ELPD diff` = se_diff)

### Model parameter table ###
# table_mod_list <- setNames(c(mod_fit_life_fldg1[mod_fit_df$mod_index[mod_fit_df$rank == 1]], 
#                              auxiliary_mod_list),
#                           nm = c('Life repro. success', 'Total nesting years', 'First nesting age', 'Last nesting age', 'Mean fledglings per year'))



############################
### TRANSLOCATION MODELS ###
############################

transloc_var <- c('mean_group_size', 'first_year_scaled', 'sex', 'transloc')
transloc_mod_var_list <- all_subset_combn(transloc_var)

transloc_mod_var_list1 <- lapply(transloc_mod_var_list, function(x) paste0('life_fldg | cens(censored) ~ ', paste(x, collapse = ' + '), ' + (1|gr(dummy,cov=A))'))


# #saving objects for the hpcc
# 
# transloc_mod_var_list1_named <- setNames(transloc_mod_var_list1, 
#                                          nm = sapply(transloc_mod_var_list, function(x) paste(x, collapse = '_')))


bprior <- c(prior_string("normal(0,0.42)", class = "b", coef = "anc_count"),
            prior_string("normal(0,0.08)", class = "b", coef = "first_year_scaled"),
            prior_string("normal(0,0.42)", class = "b", coef = "mean_group_size"),
            prior_string("normal(0,2)", class = "b", coef = "transloc_anc"),
            prior_string("normal(0, 1.5)", class = "Intercept"),
            #prior_string("normal(0,0.001)", class = "sd", coef = 'dummy'),
            prior_string("normal(0,0.02)", class = "sd", coef = 'Intercept', group = 'dummy')
)



mod_fit_life_fldg_transloc <- lapply(transloc_mod_var_list1, function(x, dat, add_rel) {
  message('Starting', x)
  
  return(
    brm(
      as.formula(x), 
      data = dat, 
      data2=list(A = add_rel),
      family = poisson(link = "log"), 
      #cov_ranef = list(A = rcw_addrel),
      chains = 4, cores = 3, iter = 8000,
      prior = prior_string("normal(0,5)", class = "b"),
      save_pars = save_pars(all = TRUE)
    )
  )
  
}, dat = repro_success_df_sexpivot1 %>% 
  filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
  mutate(censored = if_else(alive == 'yes', 'right', 'none')),
add_rel = rcw_addrel
)



mod_fit_life_fldg_transloc <- lapply(mod_fit_life_fldg_transloc, function(x) {
  message('finished')
  return(
    add_criterion(x, c("waic", "loo"))
  )
})


life_fldg_transloc_loocompare <- loo_compare(mod_fit_life_fldg_transloc[[1]],
                                             mod_fit_life_fldg_transloc[[2]],
                                             mod_fit_life_fldg_transloc[[3]],
                                             mod_fit_life_fldg_transloc[[4]],
                                             mod_fit_life_fldg_transloc[[5]],
                                             mod_fit_life_fldg_transloc[[6]],
                                             mod_fit_life_fldg_transloc[[7]],
                                             mod_fit_life_fldg_transloc[[8]],
                                             mod_fit_life_fldg_transloc[[9]],
                                             mod_fit_life_fldg_transloc[[10]],
                                             mod_fit_life_fldg_transloc[[11]],
                                             mod_fit_life_fldg_transloc[[12]],
                                             mod_fit_life_fldg_transloc[[13]],
                                             mod_fit_life_fldg_transloc[[14]],
                                             mod_fit_life_fldg_transloc[[15]])



repro_mod_draws <- repro_success_df_sexpivot1 %>% 
  filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
  mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
  data_grid(transloc = c('yes', 'no'),
            mean_group_size = quantile(mean_group_size, prob = 0.5),
            sex = 'mid',
            first_year_scaled = quantile(first_year_scaled, prob = 0.5)
  ) %>% 
  add_epred_draws(mod_fit_life_fldg_transloc[[15]],
                  ndraws = 500,
                  re_formula = NA)



color_info <- repro_success_df_sexpivot1 %>% 
  filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
  mutate(censored = if_else(alive == 'yes', 'right', 'none'),
         transloc = factor(transloc, levels = c('no', 'yes'))) %>% 
  left_join(., 
            results_list$transloc_info_color %>% 
              select(RCWid, source_color) %>% 
              rename(dummy = RCWid),
            by = 'dummy') %>% 
  mutate(source_color = if_else(is.na(source_color), "#d8d8d8", source_color)) %>% 
  select(dummy, source_color)

transloc_repro_mod_output <- list(mod_plotting = repro_mod_draws,
                                  color_info = color_info,
                                  raw_dat = repro_success_df_sexpivot1 %>% 
                                    filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
                                    mutate(transloc_final = factor(case_when(transloc == 'no' ~ 'non-translocated\n(0% transloc. ancestry)',
                                                                             transloc == 'yes' ~ 'translocated'),
                                                                   levels = c('translocated', 'non-translocated\n(0% transloc. ancestry)')))
)

saveRDS(transloc_repro_mod_output,
        file = here('results', 'transloc_mod_plot_info.RDS'))


pp_plot_transloc_mod <- pp_check(mod_fit_life_fldg_transloc[[15]], ndraws = 1000) +
  theme_classic() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 15))
pp_plot_transloc_mod$layers[[1]]$aes_params$colour <- "#ad5aad"


cowplot::ggsave2(filename = here('figures', 'supplement', 'transloc_mod_pp_check.png' ),
                 plot = pp_plot_transloc_mod,
                 width = 5*1.25, height = 4*1.25, bg = 'white')


lrs_mod_fit_table_transloceffects <- life_fldg_transloc_loocompare %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'mod_name') %>% 
  mutate(mod_index = as.integer(sub('.*\\[\\[([0-9]+)\\]\\]', '\\1', mod_name) )) %>% 
  mutate(rank = 1:n()) %>% 
  mutate(predictors = sapply(transloc_mod_var_list, function(x) paste(x, collapse = ' + '))[mod_index],
         across(where(is.numeric), function(x) round(x, 3))) %>% 
  mutate(predictors = gsub('mean_group_size', 'Mean group size', predictors),
         predictors = gsub('first_year_scaled', 'First breeding year', predictors),
         predictors = gsub('sex', 'Sex', predictors),
         predictors = gsub('anc_count', 'Ancestry count', predictors),
         predictors = gsub('transloc_anc', 'Transloc. ancestry', predictors),
         predictors = gsub('fped', 'Fped', predictors)) %>% 
  #mutate(across(where(is.numeric), function(x) round(x, 3))) %>% 
  select(rank, predictors, elpd_diff, se_diff) %>%
  rename(Rank = rank, 
         `Fixed effects` = predictors, 
         `ELPD diff` = elpd_diff, 
         `SE ELPD diff` = se_diff) %>% 
  kbl('latex', booktabs = TRUE, 
      align = "c",
      linesep = "") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position"),
                position = "center",
                font_size = 7)


cat(lrs_mod_fit_table_transloceffects, 
    file = here('tables', 'lrs_mod_transloceffects_fit_table.txt'), append = FALSE)


top_mod_transloc_effects <- as.integer(sub('.*\\[\\[([0-9]+)\\]\\]', '\\1', rownames(life_fldg_transloc_loocompare)[1]))

life_fldg_transloc_effects <- summary(mod_fit_life_fldg_transloc[[top_mod_transloc_effects]])$fixed %>%
  rbind(., summary(mod_fit_life_fldg_transloc[[top_mod_transloc_effects]])$random$dummy) %>% 
  mutate(Response = 'lifetime repro. success') %>% 
  rownames_to_column(var = 'Parameter_orig') %>% 
  mutate(Parameter = case_when(
    Parameter_orig == "mean_group_size" ~ "Mean group size",
    Parameter_orig == "first_year_scaled" ~ "First breeding year",
    Parameter_orig == "sexmid" ~ "Sex (female)",
    Parameter_orig == "anc_count" ~ "Ancestry count",
    Parameter_orig == "transloc_anc" ~ "Transloc. ancestry",
    Parameter_orig == "translocyes" ~ "Translocated",
    Parameter_orig == "fped" ~ "Fped",
    Parameter_orig == "hu_transloc_anc" ~ "Transloc. ancestry (hurdle)",
    Parameter_orig == "hu_first_year_scaled" ~ "First breeding year (hurdle)",
    Parameter_orig == "hu_mean_group_size" ~ "Mean group size (hurdle)",
    Parameter_orig == "hu_Intercept" ~ "Intercept (hurdle)",
    Parameter_orig == "sd(Intercept)" ~ "sd(Intercept) [rand. effect]",
    Parameter_orig == "sd(hu_Intercept)" ~ "sd(hurdle Intercept) [rand. effect]",
    TRUE ~ Parameter_orig
  )) %>% 
  select(Response, Parameter, Estimate, `l-95% CI`, `u-95% CI`, Rhat, Bulk_ESS, Tail_ESS) %>% 
  rename(`lower 95% CI` = `l-95% CI`,
         `upper 95% CI` = `u-95% CI`,
         `Bulk ESS` = Bulk_ESS, 
         `Tail ESS` = Tail_ESS) 

repro_success_transloc_effects_mods_params_table <- life_fldg_transloc_effects %>% 
  #bind_rows() %>% 
  mutate(across(where(is.numeric), function(x) round(x, 3))) %>% 
  kbl('latex', booktabs = TRUE,
      align = "c", linesep = "",
      escape = TRUE) %>% 
  collapse_rows(columns = 1:2,
                valign = "middle",
                latex_hline = "major") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position"),
                position = "center")


cat(repro_success_transloc_effects_mods_params_table, 
    file = here('tables', 'repro_success_transloc_effects_mods_params_table.txt'), append = FALSE)


# repro_mod_compare <- repro_mod_draws %>% 
#   #mutate(transloc = factor(transloc, levels = c('no', 'yes'))) %>% 
#   mutate(transloc_final = factor(case_when(transloc == 'no' ~ 'non-translocated\n(0% transloc. ancestry)',
#                                            transloc == 'yes' ~ 'translocated'),
#                                  levels = c('translocated', 'non-translocated\n(0% transloc. ancestry)'))) %>% 
#   ungroup() %>% 
#   ggplot()  +
#   #geom_hline(yintercept = 0, 
#   #           color = "#4c4c4c", 
#   #           linewidth = 1.5) +
#   #geom_violin(aes(x = transloc_final, y = .epred),
#   #            linewidth = 1.1, fill = "#595959", color = "#595959") +
#   geom_point(data = repro_success_df_sexpivot1 %>% 
#                filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
#                mutate(transloc_final = factor(case_when(transloc == 'no' ~ 'non-translocated\n(0% transloc. ancestry)',
#                                                         transloc == 'yes' ~ 'translocated'),
#                                               levels = c('translocated', 'non-translocated\n(0% transloc. ancestry)'))),
#              aes(x = transloc_final, y = life_fldg, color = dummy),
#              position = position_jitter(seed = 42,
#                                         width = 0.15,
#                                         height = 0),
#              #color = "#b9b9b9",
#              alpha = 0.6,
#              size = 1.4) +
#   geom_violin(aes(x = transloc_final, y = .epred),
#               linewidth = 1.1, fill = NA, color = "#595959") +
#   theme_bw() +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(linetype = 'dashed', color = '#d8d8d8'),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line.y = element_line(color = "#4c4c4c", linewidth = 1.5),
#         legend.position = 'none',
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 10)) +
#   xlab("Translocated") +
#   ylab("Life repro. success") +
#   scale_y_continuous(expand = c(0, 0)) +
#   #scale_color_manual(values = c("#b9b9b9", "#595959")) +
#   coord_cartesian(clip = 'off') +
#   scale_color_manual(values = setNames(color_info$source_color,
#                                        nm = color_info$dummy))
# 
# 
# repro_mod_compare <- repro_mod_draws %>% 
#   #mutate(transloc = factor(transloc, levels = c('no', 'yes'))) %>% 
#   mutate(transloc_final = factor(case_when(transloc == 'no' ~ 'Non-translocated\n(0% transloc. ancestry)',
#                                            transloc == 'yes' ~ 'Translocated'),
#                                  levels = c('Translocated', 'Non-translocated\n(0% transloc. ancestry)'))) %>% 
#   ungroup() %>% 
#   ggplot()  +
#   #geom_hline(yintercept = 0, 
#   #           color = "#4c4c4c", 
#   #           linewidth = 1.5) +
#   #geom_violin(aes(x = transloc_final, y = .epred),
#   #            linewidth = 1.1, fill = "#595959", color = "#595959") +
#   # geom_jitter(data = repro_success_df_sexpivot1 %>% 
#   #              filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
#   #              mutate(transloc_final = factor(case_when(transloc == 'no' ~ 'non-translocated',
#   #                                                       transloc == 'yes' ~ 'translocated'),
#   #                                             levels = c('non-translocated', 'translocated'))),
#   #            aes(x = transloc_final, y = life_fldg, color = dummy),
# #            position = position_nudge(x = -0.5, y = 0),
# #            #color = "#b9b9b9",
# #            alpha = 0.85,
# #            size = 1.2,
# #            width = 0.15,
# #            height = 0) +
# geom_point(data = repro_success_df_sexpivot1 %>%
#              filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>%
#              mutate(transloc_final = factor(case_when(transloc == 'no' ~ 'Non-translocated\n(0% transloc. ancestry)',
#                                                       transloc == 'yes' ~ 'Translocated'),
#                                             levels = c('Translocated', 'Non-translocated\n(0% transloc. ancestry)'))),
#            aes(x = transloc_final, y = life_fldg, color = dummy, fill = dummy, shape = sex),
#            position = ggpp::position_jitternudge(seed = 42,
#                                                  width = 0.19,
#                                                  height = 0,
#                                                  x = -0.24,
#                                                  nudge.from = 'jittered'),
#            #color = "#b9b9b9",
#            alpha = 0.6,
#            size = 2.5) +
#   geom_violinhalf(aes(x = transloc_final, y = .epred),
#                   linewidth = 1.1, fill = "#595959", color = "#595959") +
#   theme_bw() +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(linetype = 'dashed', color = '#d8d8d8'),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line.y = element_line(color = "#4c4c4c", linewidth = 1.5),
#         legend.position = 'none',
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 10.5),
#         axis.title.y = element_text(size = 13)) +
#   #xlab("Translocated") +
#   ylab("Life repro. success") +
#   scale_y_continuous(expand = c(0, 0)) +
#   #scale_color_manual(values = c("#b9b9b9", "#595959")) +
#   coord_cartesian(clip = 'off') +
#   scale_color_manual(values = setNames(color_info$source_color,
#                                        nm = color_info$dummy)) +
#   scale_fill_manual(values = setNames(color_info$source_color,
#                                       nm = color_info$dummy)) +
#   scale_shape_manual(values = setNames(c(22, 21),
#                                        nm = c('fid', 'mid')))
#
#
# mod_fit_life_success %>% 
#   mutate(mod_index = as.integer(sub('.*\\[\\[([0-9]+)\\]\\]', '\\1', X) )) %>% 
#   mutate(rank = 1:n()) %>% 
#   mutate(predictors = mod_formulas_life_success$V1[mod_index],
#          across(where(is.numeric), function(x) round(x, 3))) %>% 
#   mutate(predictors = gsub('mean_group_size', 'Mean group size', predictors),
#          predictors = gsub('first_year_scaled', 'First breeding year', predictors),
#          predictors = gsub('sex', 'Sex', predictors),
#          predictors = gsub('anc_count', 'Ancestry count', predictors),
#          predictors = gsub('transloc_anc', 'Transloc. ancestry', predictors),
#          predictors = gsub('fped', 'Fped', predictors)) %>% 
#   #mutate(across(where(is.numeric), function(x) round(x, 3))) %>% 
#   select(rank, predictors, elpd_diff, se_diff) %>%
#   rename(Rank = rank, 
#          `Fixed effects` = predictors, 
#          `ELPD diff` = elpd_diff, 
#          `SE ELPD diff` = se_diff)


### Model parameter table ###
# table_mod_list <- setNames(c(mod_fit_life_fldg1[mod_fit_df$mod_index[mod_fit_df$rank == 1]], 
#                              auxiliary_mod_list),
#                           nm = c('Life repro. success', 'Total nesting years', 'First nesting age', 'Last nesting age', 'Mean fledglings per year'))



###########################################
### CODE FOR RUNNING THE MODELS LOCALLY ###
###########################################
# ### fit all combos of models ###
# if (isTRUE(run_all_mods)) {
#   mod_fit_life_fldg <- lapply(lrs_mod_var_list1, function(x, dat, add_rel) {
#     message('Starting', x)
#     
#     return(
#       brm(
#         as.formula(x), 
#         data = dat, 
#         data2=list(A = add_rel),
#         family = poisson(link = "log"), 
#         #cov_ranef = list(A = rcw_addrel),
#         chains = 3, cores = 3, iter = 10000,
#         save_pars = save_pars(all = TRUE)
#       )
#     )
#     
#   }, dat = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#   add_rel = rcw_addrel
#   )
#   
#   
#   test_mod_with_fp <- brm(
#     life_fldg | cens(censored) ~ mean_group_size + first_year_scaled + sex + transloc_anc + anc_count + fped + (1|gr(dummy,cov=A)), 
#     data = repro_success_df_sexpivot1 %>%
#       #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#       filter(birth_during_monitoring == 'yes') %>%
#       mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#     data2=list(A = rcw_addrel),
#     family = poisson(link = "log"), 
#     #cov_ranef = list(A = rcw_addrel),
#     chains = 3, cores = 3, iter = 10000,
#     save_pars = save_pars(all = TRUE)
#   )
#   
#   # lrs_mod_var_list1_named[[57]]
#   # lrs_mod_var_list1_named[[42]]
#   # lrs_mod_var_list1_named[[58]]
#   # lrs_mod_var_list1_named[[45]]
#   
#   mod_fit <- loo_compare(mod_list_addcriterion[[1]],
#                          mod_list_addcriterion[[2]],
#                          mod_list_addcriterion[[3]],
#                          mod_list_addcriterion[[4]],
#                          mod_list_addcriterion[[5]],
#                          mod_list_addcriterion[[6]],
#                          mod_list_addcriterion[[7]],
#                          mod_list_addcriterion[[8]],
#                          mod_list_addcriterion[[9]],
#                          mod_list_addcriterion[[10]],
#                          mod_list_addcriterion[[11]],
#                          mod_list_addcriterion[[12]],
#                          mod_list_addcriterion[[13]],
#                          mod_list_addcriterion[[14]],
#                          mod_list_addcriterion[[15]],
#                          mod_list_addcriterion[[16]],
#                          mod_list_addcriterion[[17]],
#                          mod_list_addcriterion[[18]],
#                          mod_list_addcriterion[[19]],
#                          mod_list_addcriterion[[20]],
#                          mod_list_addcriterion[[21]],
#                          mod_list_addcriterion[[22]],
#                          mod_list_addcriterion[[23]],
#                          mod_list_addcriterion[[24]],
#                          mod_list_addcriterion[[25]],
#                          mod_list_addcriterion[[26]],
#                          mod_list_addcriterion[[27]],
#                          mod_list_addcriterion[[28]],
#                          mod_list_addcriterion[[29]],
#                          mod_list_addcriterion[[30]],
#                          mod_list_addcriterion[[31]],
#                          mod_list_addcriterion[[32]],
#                          mod_list_addcriterion[[33]],
#                          mod_list_addcriterion[[34]],
#                          mod_list_addcriterion[[35]],
#                          mod_list_addcriterion[[36]],
#                          mod_list_addcriterion[[37]],
#                          mod_list_addcriterion[[38]],
#                          mod_list_addcriterion[[39]],
#                          mod_list_addcriterion[[40]],
#                          mod_list_addcriterion[[41]],
#                          mod_list_addcriterion[[42]],
#                          mod_list_addcriterion[[43]],
#                          mod_list_addcriterion[[44]],
#                          mod_list_addcriterion[[45]],
#                          mod_list_addcriterion[[46]],
#                          mod_list_addcriterion[[47]],
#                          mod_list_addcriterion[[48]],
#                          mod_list_addcriterion[[49]],
#                          mod_list_addcriterion[[50]],
#                          mod_list_addcriterion[[51]],
#                          mod_list_addcriterion[[52]],
#                          mod_list_addcriterion[[53]],
#                          mod_list_addcriterion[[54]],
#                          mod_list_addcriterion[[55]],
#                          mod_list_addcriterion[[56]],
#                          mod_list_addcriterion[[57]],
#                          mod_list_addcriterion[[58]],
#                          mod_list_addcriterion[[59]],
#                          mod_list_addcriterion[[60]],
#                          mod_list_addcriterion[[61]],
#                          mod_list_addcriterion[[62]],
#                          mod_list_addcriterion[[63]],
#                          criterion = "loo")
#   
#   # write.csv(as.data.frame(mod_fit),
#   #           row.names = TRUE)
#   
#   mod_fit_life_fldg1 <- lapply(mod_fit_life_fldg, function(x) {
#     message('finished')
#     return(
#       add_criterion(x, c("waic", "loo"))
#     )
#   })
#   
#   
#   mod_fit <- loo_compare(mod_fit_life_fldg1[[1]],
#                          mod_fit_life_fldg1[[2]],
#                          mod_fit_life_fldg1[[3]],
#                          mod_fit_life_fldg1[[4]],
#                          mod_fit_life_fldg1[[5]],
#                          mod_fit_life_fldg1[[6]],
#                          mod_fit_life_fldg1[[7]],
#                          mod_fit_life_fldg1[[8]],
#                          mod_fit_life_fldg1[[9]],
#                          mod_fit_life_fldg1[[10]],
#                          mod_fit_life_fldg1[[11]],
#                          mod_fit_life_fldg1[[12]],
#                          mod_fit_life_fldg1[[13]],
#                          mod_fit_life_fldg1[[14]],
#                          mod_fit_life_fldg1[[15]],
#                          mod_fit_life_fldg1[[16]],
#                          mod_fit_life_fldg1[[17]],
#                          mod_fit_life_fldg1[[18]],
#                          mod_fit_life_fldg1[[19]],
#                          mod_fit_life_fldg1[[20]],
#                          mod_fit_life_fldg1[[21]],
#                          mod_fit_life_fldg1[[22]],
#                          mod_fit_life_fldg1[[23]],
#                          mod_fit_life_fldg1[[24]],
#                          mod_fit_life_fldg1[[25]],
#                          mod_fit_life_fldg1[[26]],
#                          mod_fit_life_fldg1[[27]],
#                          mod_fit_life_fldg1[[28]],
#                          mod_fit_life_fldg1[[29]],
#                          mod_fit_life_fldg1[[30]],
#                          mod_fit_life_fldg1[[31]],
#                          criterion = "loo")
#   
#   mod_fit_df <- as.data.frame(mod_fit) %>% 
#     rownames_to_column(var = 'mod_index') %>% 
#     mutate(mod_index = as.integer(sub('.*\\[\\[([0-9]+)\\]\\]', '\\1', mod_index) )) %>% 
#     mutate(rank = 1:n())
#   
#   mod_fit_df$full_formula <- sapply(mod_fit_df$mod_index, function(index, mod) {as.character(mod[[index]]$formula)[1]}, mod = mod_fit_life_fldg1)
#   mod_fit_df$predictors <- sub(".*~ ([a-z]+.*) \\+ \\(.*", "\\1", mod_fit_df$full_formula)
#   
#   # GGally::ggpairs(repro_success_df_sexpivot1 %>%
#   #                   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   #                   filter(birth_during_monitoring == 'yes') %>%
#   #                   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#   #                          birth_year = MinAge - min(MinAge)) %>% 
#   #                   select(mean_group_size, first_year_scaled, sex, transloc_anc, anc_count_z))
#   # 
# }
 


#################################
### CODE NOT CURRENTLY IN USE ###
#################################
# 
# 
# transloc_first_year <- repro_success_df_sexpivot1 %>% 
#   filter(transloc == 'yes') %>% 
#   pull(first_year_scaled) %>% 
#   unique()
# 
# test_prior_basic <- brm(
#   #life_fldg | cens(censored) ~ mean_group_size + transloc + sex + first_year_scaled + (1 | gr(dummy, cov = A)),
#   life_fldg | cens(censored) ~ transloc + (1 | gr(dummy, cov = A)),
#   data = repro_success_df_sexpivot1 %>% 
#     filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
#     filter(first_year_scaled >= 8),
#     #filter(first_year_scaled %in% transloc_first_year),
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"),
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, iter = 3000#,
#   #prior = bprior,
#  # sample_prior = 'only'
# )
# 
# test_prior1 <- brm(
#   #life_fldg | cens(censored) ~ mean_group_size + transloc + sex + first_year_scaled + (1 | gr(dummy, cov = A)),
#   life_fldg | cens(censored) ~ transloc + first_year_scaled + (1 | gr(dummy, cov = A)),
#   data = repro_success_df_sexpivot1 %>% 
#     filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
#     filter(first_year_scaled >= 8),
#   #filter(first_year_scaled %in% transloc_first_year),
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"),
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, iter = 3000#,
#   #prior = bprior,
#   # sample_prior = 'only'
# )
# 
# loo_compare(add_criterion(test_prior_basic, criterion = 'loo'),
#             add_criterion(test_prior1, criterion = 'loo'))
# 
# repro_success_df_sexpivot1 %>% 
#   filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
#   mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
#   filter(first_year_scaled %in% transloc_first_year) %>% 
#   ggplot() +
#   geom_violin(aes(x = transloc, y = first_year_scaled))
# 
# 
# 
# test_draws <- repro_success_df_sexpivot1 %>% 
#   filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
#   mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
#   data_grid(transloc = c('yes', 'no'),
#             mean_group_size = quantile(mean_group_size, prob = 0.5),
#             sex = 'mid',
#             first_year_scaled = quantile(first_year_scaled, prob = 0.5)
#   ) %>% 
#   add_epred_draws(test_prior,
#                   ndraws = 500,
#                   re_formula = NA)
# 
# 
# test_draws %>% 
#   ggplot() +
#   geom_violin(aes(x = transloc, y = .epred)) +
#   geom_point(data = repro_success_df_sexpivot1 %>% 
#                filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
#                mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#              aes(x = transloc, y = life_fldg))
# 
# 
# test_draws %>% 
#   ggplot()  +
#   geom_violin(data = repro_success_df_sexpivot1 %>% 
#                 filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
#                 mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#               aes(x = transloc, y = life_fldg)) +
#   geom_violin(aes(x = transloc, y = .epred))
# 
# repro_success_df_sexpivot1 %>% 
#   filter( (transloc_anc < 1e-13 & birth_during_monitoring == 'yes') | transloc == 'yes') %>% 
#   mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
#   ggplot() +
#   geom_boxplot(aes(x = transloc, y = life_fldg))
# 
# 
# 
# 
# test_prior <- brm(
#   life_fldg | cens(censored) ~ mean_group_size + anc_count + first_year_scaled + transloc_anc + (1 | gr(dummy, cov = A)),
#   data = repro_success_df_sexpivot1 %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     #filter(age != 0) %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"),
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 1, iter = 2000,
#   prior = bprior,
#   sample_prior = 'only'
# )
# 
# pp_check(test_prior, ndraws = 100) +
#   xlim(0, 15)
# 
# 
# default_prior(life_fldg | cens(censored) ~ mean_group_size + first_year_scaled + transloc_anc + anc_count_z + (1 | gr(dummy, cov = A)),
#               data = repro_success_df_sexpivot1 %>%
#                 filter(birth_during_monitoring == 'yes') %>%
#                 #filter(age != 0) %>%
#                 mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#               data2=list(A = rcw_addrel),
#               family = poisson(link = 'log'))
# 
# 
# default_prior(mean_fldg_per_year | cens(censored) ~ mean_group_size + first_year_scaled + transloc_anc + anc_count_z + (1 | gr(dummy, cov = A)),
#               data = repro_success_df_sexpivot1 %>%
#                 filter(birth_during_monitoring == 'yes') %>%
#                 #filter(age != 0) %>%
#                 mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#               data2=list(A = rcw_addrel),
#               family = hurdle_gamma())
# 
# 
# stancode(mean_fldg_per_year | cens(censored) ~ mean_group_size + first_year_scaled + transloc_anc + anc_count_z + (1 | gr(dummy, cov = A)),
#          data = repro_success_df_sexpivot1 %>%
#            filter(birth_during_monitoring == 'yes') %>%
#            #filter(age != 0) %>%
#            mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#          data2=list(A = rcw_addrel),
#          family = hurdle_gamma())
# 
# default_prior(mean_fldg_per_year | cens(censored) ~ mean_group_size + first_year_scaled + transloc_anc + anc_count_z + (1 | gr(dummy, cov = A)),
#               data = repro_success_df_sexpivot1 %>%
#                 filter(birth_during_monitoring == 'yes') %>%
#                 #filter(age != 0) %>%
#                 mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#               data2=list(A = rcw_addrel),
#               family = hurdle_gamma())
# 
# 
# 
# 
# get_prior(mean_fldg_per_year | cens(censored) ~ mean_group_size + first_year_scaled + transloc_anc + anc_count_z + (1 | gr(dummy, cov = A)),
#           data = repro_success_df_sexpivot1 %>%
#             filter(birth_during_monitoring == 'yes') %>%
#             #filter(age != 0) %>%
#             mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#           data2=list(A = rcw_addrel),
#           family = hurdle_gamma())
# 
# pp_check(test_prior)
# 
# 
# stancode(birth_firstbreed_dif | cens(censored) ~ mean_group_size + first_year_scaled + transloc_anc + anc_count_z + (1 | gr(dummy, cov = A)),
#          data = repro_success_df_sexpivot1 %>%
#            filter(birth_during_monitoring == 'yes') %>%
#            #filter(age != 0) %>%
#            mutate(censored = if_else(alive == 'yes', 'right', 'none')),
#          family = poisson(link = "log"),
#          data2=list(A = rcw_addrel),
#          prior = bprior)
#
#
# ### TOTAL NESTING YEARS MODELS ###
# auxiliary_mod_list[['nestyear_mod']] <- brm(
#   as.formula(gsub('life_fldg', 'nesting_years', as.character(mod_fit_life_fldg1[[mod_fit_df$mod_index[mod_fit_df$rank == 1]]]$formula)[1])), 
#   data = repro_success_df_sexpivot1 %>% 
#     filter(birth_during_monitoring == 'yes') %>% 
#     #filter(age != 0) %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.9)
# )
# 
# 
# ### TOTAL NESTING YEARS MODELS ###
# auxiliary_mod_list[['mean_fldg_mod']] <- brm(
#   as.formula(gsub('life_fldg', 'mean_fldg_per_year', as.character(mod_fit_life_fldg1[[31]]$formula)[1])), 
#   data = repro_success_df_sexpivot1 %>% 
#     filter(birth_during_monitoring == 'yes') %>% 
#     #filter(age != 0) %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = hurdle_gamma(), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.95)
# )
# 
# mean_fldg_mod_alt <- brm(
#   bf(
#     as.formula(gsub('life_fldg', 'mean_fldg_per_year', as.character(mod_fit_life_fldg1[[31]]$formula)[1])),
#     hu ~ first_year_scaled + transloc_anc + (1 | gr(dummy, cov = A))
#   ),
#   data = repro_success_df_sexpivot1 %>% 
#     filter(birth_during_monitoring == 'yes') %>% 
#     #filter(age != 0) %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = hurdle_gamma(), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.95)
# )
#
#
# mean_fldg_mod_alt1 <- brm(
#   bf(
#     as.formula(gsub('life_fldg', 'mean_fldg_per_year', as.character(mod_fit_life_fldg1[[31]]$formula)[1])),
#     hu ~ first_year_scaled + transloc_anc + (1 | gr(dummy, cov = A))
#   ),
#   data = repro_success_df_sexpivot1 %>% 
#     filter(birth_during_monitoring == 'yes') %>% 
#     #filter(age != 0) %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = hurdle_gamma(), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.95)
# )


# ### TOTAL NESTING YEARS MODELS ###
# auxiliary_mod_list[['firstbreed_mod']] <- brm(
#   as.formula(gsub('life_fldg', 'birth_firstbreed_dif', as.character(mod_fit_life_fldg1[[31]]$formula)[1])),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
#
# test_pp$layers[[1]]$computed_mapping<- 'black'
# 
# 
# test_pp$layers[[1]]$aes_params$alpha <- 0.7
# test_pp$layers[[1]]$aes_params$colour <- 'gray'
# 
# 
# 
# 
# test_pp <- brms::pp_check(mean_fldg_mod_alt1,
#          ndraws = 1000)
# 
# pp_check(lastbreed_mod,
#          ndraws = 1000)


# #MODEL 3: MIN NESTING YEAR
# firstbreed_mod <- brm(
#   as.formula(gsub('life_fldg', 'birth_firstbreed_dif', as.character(mod_fit_life_fldg1[[31]]$formula)[1])),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
# 
# 
# firstbreed_mod1 <- brm(
#   birth_firstbreed_dif | cens(censored) ~ mean_group_size + first_year_scaled + transloc_anc + anc_count_z + (1 | gr(dummy, cov = A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
# 
# 
# firstbreed_mod2 <- brm(
#   birth_firstbreed_dif | cens(censored) ~ first_year_scaled + transloc_anc + anc_count_z + (1 | gr(dummy, cov = A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
# 
# #MODEL 4: MAX NESTING YEAR
# lastbreed_mod <- brm(
#   as.formula(gsub('life_fldg', 'birth_lastbreed_dif', as.character(mod_fit_life_fldg1[[31]]$formula)[1])),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
# 
# lapply(c(mod_fit_life_fldg1[[31]], nestyear_mod, firstbreed_mod, lastbreed_mod), function(MOD) message(class(MOD)))
# 
# repro_draws_list <- lapply(setNames(nm = names(repro_mod_list)), 
#                            function(MOD, mod_list) {
#                              output_list <- list()
#                              #message('starting ')
#                              beta_draws <- as_draws_df(mod_list[[MOD]])
#                              
#                              ci_95 <- quantile(beta_draws$b_transloc_anc,
#                                                prob = c(0.025, 0.975))
#                              output_list[['ci_95']] <- ci_95
#                              output_list[['draws']] <- beta_draws %>% 
#                                select(b_transloc_anc) %>% 
#                                mutate(in_ci = if_else(b_transloc_anc < ci_95[1] | b_transloc_anc > ci_95[2], 'out', 'in'),
#                                       group = 'one') %>% 
#                                mutate(mod = case_when(MOD == 'life_fldg' ~ 'Lifetime\nrepro.\nsuccess',
#                                                       MOD == 'nestyear_mod' ~ 'Total\nnesting\nyears',
#                                                       MOD == 'firstbreed' ~ 'First\nnesting\nage',
#                                                       MOD == 'lastbreed' ~ 'Last\nnesting\nage'))
#                              
#                              return(output_list)
#                            }, mod_list = repro_mod_list) 
# 
# repro_draws_only <- lapply(repro_draws_list, function(x) x$draws) %>% 
#   bind_rows()
# 
# 
# mod_pred_list_transloc <- lapply(setNames(list(mod_fit_life_fldg1[[31]], nestyear_mod, firstbreed_mod, lastbreed_mod),
#                                    nm = c('life_fldg', 'nestyear_mod', 'firstbreed', 'lastbreed')), 
#                         function(MOD_OBJ, dat) {
#                           
#                           dat %>% 
#                             data_grid(transloc_anc = seq_range(transloc_anc, n = 100),
#                                       mean_group_size = quantile(mean_group_size, prob = 0.5),
#                                       sex = 'mid',
#                                       first_year_scaled = quantile(first_year_scaled, prob = 0.5),
#                                       anc_count_z = quantile(anc_count_z, prob = 0.5)
#                                       ) %>% 
#                             add_epred_draws(MOD_OBJ,
#                                             ndraws = 500,
#                                             re_formula = NA)
#                         },
#                         dat = repro_success_df_sexpivot1 %>% 
#                           #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#                           filter(birth_during_monitoring == 'yes')  %>%
#                           mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#                                  birth_year = MinAge - min(MinAge))
#                         )
# 
# # test <- repro_success_df_sexpivot1 %>% 
# #   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
# #   filter(birth_during_monitoring == 'yes')  %>%
# #   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
# #          birth_year = MinAge - min(MinAge)) %>% 
# #   data_grid(transloc_anc = seq_range(transloc_anc, n = 100),
# #             mean_group_size = quantile(mean_group_size, prob = 0.5),
# #             sex = 'mid',
# #             first_year_scaled = quantile(first_year_scaled, prob = 0.5),
# #             anc_count_z = quantile(anc_count_z, prob = 0.5)
# #   ) %>% 
# #   add_epred_draws(nestyear_mod,
# #                   ndraws = 100,
# #                   re_formula = NA,
# #                   seed = 5093) %>% 
# #   mutate(draw_val = `.draw`)
# # 
# # test2 <- spread_draws(nestyear_mod, b_transloc_anc, 
# #                       ndraws = 100, seed = 5093) %>% 
# #   rownames_to_column(var = 'draw_val') %>% 
# #   mutate(draw_val = as.numeric(draw_val))
# # 
# # test3 <- test %>% 
# #   left_join(., test2,
# #             by = 'draw_val')
# # 
# # test3 %>% 
# #   ggplot() +
# #   geom_line(aes(x = transloc_anc, y = `.epred`, group = draw_val,
# #                 color = b_transloc_anc),
# #             linewidth = 1, 
# #             #color = '#c7b6d8', 
# #             alpha = 0.6) +
# #   scale_color_gradient2(midpoint = 0.56,
# #                         low = '#c5625b', 
# #                         mid = '#f2f2f2', 
# #                         high = '#326934'#,
# #                         #limits = c(lower_lim_sina, upper_lim_sina)
# #   )
# # 
# 
# transloc_base_violin_plot <- repro_draws_only %>%
#   mutate(mod = factor(mod, levels = rev(c("Lifetime\nrepro.\nsuccess", "Total\nnesting\nyears", "First\nnesting\nage", "Last\nnesting\nage")))) %>% 
#   ggplot() +
#   geom_violin(aes(x = factor(mod),
#                   y = b_transloc_anc))
# 
# transloc_base_violin_plot_build <- ggplot2::ggplot_build(transloc_base_violin_plot)$data[[1]]
# transloc_base_violin_plot_build <- transform(transloc_base_violin_plot_build,
#                                              xminv = x - violinwidth * (x - xmin),
#                                              xmaxv = x + violinwidth * (xmax - x))
# 
# transloc_base_violin_plot_build <- rbind(plyr::arrange(transform(transloc_base_violin_plot_build, x = xminv), y),
#                                          plyr::arrange(transform(transloc_base_violin_plot_build, x = xmaxv), -y))
# 
# #1: last nesting age
# #2: first nesting age
# #3: total nesting years
# #4 lifetime reproductive success
# 
# transloc_base_violin_plot_build <- transloc_base_violin_plot_build %>% 
#   mutate(fill_group = case_when(group == 1 & y > repro_draws_list$lastbreed$ci_95['2.5%'] & y < repro_draws_list$lastbreed$ci_95['97.5%'] ~ 'inside',
#                                 group == 1 & (y < repro_draws_list$lastbreed$ci_95['2.5%'] | y > repro_draws_list$lastbreed$ci_95['97.5%']) ~ 'outside',
#                                 group == 2 & y > repro_draws_list$firstbreed$ci_95['2.5%'] & y < repro_draws_list$firstbreed$ci_95['97.5%'] ~ 'inside',
#                                 group == 2 & (y < repro_draws_list$firstbreed$ci_95['2.5%'] | y > repro_draws_list$firstbreed$ci_95['97.5%']) ~ 'outside',
#                                 group == 3 & y > repro_draws_list$nestyear_mod$ci_95['2.5%'] & y < repro_draws_list$nestyear_mod$ci_95['97.5%'] ~ 'inside',
#                                 group == 3 & (y < repro_draws_list$nestyear_mod$ci_95['2.5%'] | y > repro_draws_list$nestyear_mod$ci_95['97.5%']) ~ 'outside',
#                                 group == 4 & y > repro_draws_list$life_fldg$ci_95['2.5%'] & y < repro_draws_list$life_fldg$ci_95['97.5%'] ~ 'inside',
#                                 group == 4 & (y < repro_draws_list$life_fldg$ci_95['2.5%'] | y > repro_draws_list$life_fldg$ci_95['97.5%']) ~ 'outside'))
# 
# 
# #  $fill_group <- ifelse(p_build$y >= 0,'Above','Below')
# #This is necessary to ensure that instead of trying to draw
# # 3 polygons, we're telling ggplot to draw six polygons
# transloc_base_violin_plot_build$group1 <- with(transloc_base_violin_plot_build,
#                                                interaction(factor(group),factor(fill_group)))
# 
# 
# transloc_mod_color_vec <- setNames(c('#3f0b70', '#04ab80', '#d62828', '#3478e5'),
#                                    nm = names(mod_pred_list_transloc))
# 
# transloc_mod_color_post <- setNames(c('#c7b6d8', '#9addcc', '#eea9a9', '#adc9f4'),
#                                     nm = names(mod_pred_list_transloc))
# 
# 
# transloc_violin_post_plot <- transloc_base_violin_plot_build %>% 
#   filter(fill_group == 'inside') %>% 
#   ggplot() + 
#   geom_violin(data = repro_draws_only %>%
#                 mutate(mod = factor(mod, levels = rev(c("Lifetime\nrepro.\nsuccess", "Total\nnesting\nyears", "First\nnesting\nage", "Last\nnesting\nage")))), 
#               aes(x = mod, y = b_transloc_anc, fill = mod),
#               color = NA, alpha = 0.8) +
#   scale_fill_manual(values = unname(rev(transloc_mod_color_post))) +
#   ggnewscale::new_scale_fill() +
#   geom_polygon(aes(x = x,y = y, group = group1, fill = group1)) +
#   scale_fill_manual(values = unname(rev(transloc_mod_color_vec))) +
#   theme(legend.position = 'none',
#         axis.line.x = element_line(color = 'black', linewidth = 0.5),
#         panel.background = element_blank(),
#         #panel.grid.major.y = element_line(color = '#d8d8d8', linetype = 'dashed', linewidth = 0.4),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(color = rev(transloc_mod_color_vec), face = 'bold')) +
#   coord_flip() +
#   ylab('Translocation coefficient') +
#   #scale_x_discrete(name = NULL, sec.axis = sec_axis(~., name = "Y-axis Label on right side")) +
#   guides(y = "none") +
#   guides(y.sec = guide_axis_label_trans(~.x)) +
#   geom_hline(yintercept = 0, linetype = 'dashed', color = 'black')
# 
# 
# transloc_mod_color_vec <- setNames(c('#3f0b70', '#04ab80', '#d62828', '#3478e5'),
#                           nm = names(mod_pred_list_transloc))
# 
# transloc_mod_color_post <- setNames(c('#c7b6d8', '#9addcc', '#eea9a9', '#adc9f4'),
#                                    nm = names(mod_pred_list_transloc))
# 
# 
# mod_plot_list_transloc <- lapply(setNames(nm = names(mod_pred_list_transloc)), 
#                                  function(PRED_OBJ, mod, color_vec, color_post) {
#   mod[[PRED_OBJ]] %>% 
#     group_by(transloc_anc) %>% 
#     summarize(lower25 = quantile(.epred, prob = 0.025),
#               lower10 = quantile(.epred, prob = 0.10),
#               lower20 = quantile(.epred, prob = 0.20),
#               median = quantile(.epred, prob = 0.5),
#               upper80 = quantile(.epred, prob = 0.80),
#               upper90 = quantile(.epred, prob = 0.90),
#               upper975 = quantile(.epred, prob = 0.975),
#               .groups = 'drop') %>% 
#     ggplot(aes(x = transloc_anc, y = median)) + 
#     geom_line(data = mod[[PRED_OBJ]],
#               aes(x = transloc_anc, y = `.epred`, group = .draw),
#               linewidth = 0.25, 
#               color = color_post[PRED_OBJ], 
#               alpha = 0.4) +
#     geom_ribbon(aes(transloc_anc,
#                     ymin = lower25, ymax = upper975),
#                 color = color_vec[PRED_OBJ], linetype = 'dashed', linewidth = 1.5,
#                 lineend = 'round', fill = NA,
#                 #fill = '#3f0b70', #'#c7b6d8',
#                 alpha = 0.2) +
#     geom_line(linewidth = 3, color = color_vec[PRED_OBJ],
#               lineend = 'round') +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank(),
#           panel.grid.major = element_line(linetype = 'dashed', 
#                                           linewidth = 0.4,
#                                           color = '#d8d8d8'))
# }, mod = mod_pred_list_transloc, color_vec = transloc_mod_color_vec, color_post = transloc_mod_color_post)
# 
# 
# 
# repro_combined_plot <- cowplot::plot_grid(mod_plot_list_transloc$life_fldg +
#                      theme(axis.title.x = element_blank(),
#                            panel.border = element_blank(),
#                            axis.line = element_line(linewidth = 0.5)) +
#                      ylab('Life reproductive success'),
#                    mod_plot_list_transloc$nestyear_mod +
#                      theme(axis.title.x = element_blank(),
#                            panel.border = element_blank(),
#                            axis.line = element_line(linewidth = 0.5)) +
#                      ylab('Total nesting years'),
#                    mod_plot_list_transloc$firstbreed +
#                      theme(axis.title.x = element_blank(),
#                            panel.border = element_blank(),
#                            axis.line = element_line(linewidth = 0.5)) +
#                      ylab('First nesting age'),
#                    mod_plot_list_transloc$lastbreed +
#                      theme(axis.title.x = element_blank(),
#                            panel.border = element_blank(),
#                            axis.line = element_line(linewidth = 0.5)) +
#                      ylab('Last nesting age'),
#                    ncol = 2,
#                    align = 'hv',
#                    axis = 'lr',
#                    labels = c('A', 'B', 'C', 'D'),
#                    label_fontface = 'bold')
# 
# 
# cowplot::plot_grid(repro_combined_plot,
#                    mod_coefficient_sina,
#                    ncol = 2, rel_widths = c(0.8, 0.2))
# 
# 
# library(grid)
# library(gridExtra)
# 
# x_title <- textGrob("Expected translocation ancestry", 
#                    gp=gpar(fontface="plain", col="black", fontsize = 12))
# 
# #add to plot
# 
# transloc_repro_mod_multipan <- cowplot::plot_grid(
#   grid.arrange(arrangeGrob(repro_combined_plot, 
#                            bottom = textGrob("Expected translocation ancestry", 
#                                             gp=gpar(fontface="plain", col="black", fontsize = 12))
#                            )),
#   grid.arrange(arrangeGrob(transloc_violin_post_plot + #mod_coefficient_sina +
#                              theme(axis.title.x = element_blank()), 
#                            bottom = textGrob("Translocation coefficent", 
#                                              gp=gpar(fontface="plain", col="black", fontsize = 12))
#                            )),
#   ncol = 2, rel_widths = c(0.6, 0.2),
#   labels = c('', 'E'), label_fontface = 'bold')
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'transloc_repro_mod_multipan.png' ),
#                  plot = transloc_repro_mod_multipan,
#                  width = 9*1.25, height = 5*1.25, bg = 'white')
# 
# 
# posterior_samples(repro_mod_list[[1]], add_chain = TRUE)
# 
# ggplot() +
#   geom_line(data = as_draws_array(repro_mod_list[[1]], variable = "b_transloc_anc")[,1,] %>% 
#               as.data.frame() %>% 
#               mutate(step = 1:n()),
#             aes(x = step, y = `1.b_transloc_anc`))
# 
# 
# 
# 
# performance::check_collinearity(repro_mod_list[[1]])
# performance::check_collinearity(repro_mod_list[[2]])
# performance::check_collinearity(repro_mod_list[[3]])
# performance::check_collinearity(repro_mod_list[[4]])
# 
# names(repro_mod_list)
# 
# 
# 
# 
# repro_draws_list_firstyear <- lapply(setNames(nm = names(repro_mod_list)), 
#                            function(MOD, mod_list) {
#                              output_list <- list()
#                              #message('starting ')
#                              beta_draws <- as_draws_df(mod_list[[MOD]])
#                              
#                              ci_95 <- quantile(beta_draws$b_first_year_scaled,
#                                                prob = c(0.025, 0.975))
#                              output_list[['ci_95']] <- ci_95
#                              output_list[['draws']] <- beta_draws %>% 
#                                select(b_first_year_scaled) %>% 
#                                mutate(in_ci = if_else(b_first_year_scaled < ci_95[1] | b_first_year_scaled > ci_95[2], 'out', 'in'),
#                                       group = 'one') %>% 
#                                mutate(mod = case_when(MOD == 'life_fldg' ~ 'Lifetime\nrepro.\nsuccess',
#                                                       MOD == 'nestyear_mod' ~ 'Total\nnesting\nyears',
#                                                       MOD == 'firstbreed' ~ 'First\nnesting\nage',
#                                                       MOD == 'lastbreed' ~ 'Last\nnesting\nage'))
#                              
#                              return(output_list)
#                            }, mod_list = repro_mod_list) 
# 
# repro_draws_only_firstyear <- lapply(repro_draws_list_firstyear, function(x) x$draws) %>% 
#   bind_rows()
# 
# 
# mod_pred_list_firstyear <- lapply(setNames(list(mod_fit_life_fldg1[[31]], nestyear_mod, firstbreed_mod, lastbreed_mod),
#                                           nm = c('life_fldg', 'nestyear_mod', 'firstbreed', 'lastbreed')), 
#                                  function(MOD_OBJ, dat) {
#                                    
#                                    dat %>% 
#                                      data_grid(transloc_anc = quantile(transloc_anc, prob = 0.5),
#                                                mean_group_size = quantile(mean_group_size, prob = 0.5),
#                                                sex = 'mid',
#                                                first_year_scaled = seq_range(first_year_scaled, n = 100),
#                                                anc_count_z = quantile(anc_count_z, prob = 0.5)
#                                      ) %>% 
#                                      add_epred_draws(MOD_OBJ,
#                                                      ndraws = 500,
#                                                      re_formula = NA)
#                                  },
#                                  dat = repro_success_df_sexpivot1 %>% 
#                                    #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#                                    filter(birth_during_monitoring == 'yes')  %>%
#                                    mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#                                           birth_year = MinAge - min(MinAge))
# )
# 
# 
# 
# firstyear_base_violin_plot <- repro_draws_only_firstyear %>%
#   mutate(mod = factor(mod, levels = rev(c("Lifetime\nrepro.\nsuccess", "Total\nnesting\nyears", "First\nnesting\nage", "Last\nnesting\nage")))) %>% 
#   ggplot() +
#   geom_violin(aes(x = factor(mod),
#                   y = b_first_year_scaled))
# 
# firstyear_base_violin_plot_build <- ggplot2::ggplot_build(firstyear_base_violin_plot)$data[[1]]
# firstyear_base_violin_plot_build <- transform(firstyear_base_violin_plot_build,
#                                              xminv = x - violinwidth * (x - xmin),
#                                              xmaxv = x + violinwidth * (xmax - x))
# 
# firstyear_base_violin_plot_build <- rbind(plyr::arrange(transform(firstyear_base_violin_plot_build, x = xminv), y),
#                                          plyr::arrange(transform(firstyear_base_violin_plot_build, x = xmaxv), -y))
# 
# 
# #1: last nesting age
# #2: first nesting age
# #3: total nesting years
# #4 lifetime reproductive success
# 
# firstyear_base_violin_plot_build <- firstyear_base_violin_plot_build %>% 
#   mutate(fill_group = case_when(group == 1 & y > repro_draws_list_firstyear$lastbreed$ci_95['2.5%'] & y < repro_draws_list_firstyear$lastbreed$ci_95['97.5%'] ~ 'inside',
#                                 group == 1 & (y < repro_draws_list_firstyear$lastbreed$ci_95['2.5%'] | y > repro_draws_list_firstyear$lastbreed$ci_95['97.5%']) ~ 'outside',
#                                 group == 2 & y > repro_draws_list_firstyear$firstbreed$ci_95['2.5%'] & y < repro_draws_list_firstyear$firstbreed$ci_95['97.5%'] ~ 'inside',
#                                 group == 2 & (y < repro_draws_list_firstyear$firstbreed$ci_95['2.5%'] | y > repro_draws_list_firstyear$firstbreed$ci_95['97.5%']) ~ 'outside',
#                                 group == 3 & y > repro_draws_list_firstyear$nestyear_mod$ci_95['2.5%'] & y < repro_draws_list_firstyear$nestyear_mod$ci_95['97.5%'] ~ 'inside',
#                                 group == 3 & (y < repro_draws_list_firstyear$nestyear_mod$ci_95['2.5%'] | y > repro_draws_list_firstyear$nestyear_mod$ci_95['97.5%']) ~ 'outside',
#                                 group == 4 & y > repro_draws_list_firstyear$life_fldg$ci_95['2.5%'] & y < repro_draws_list_firstyear$life_fldg$ci_95['97.5%'] ~ 'inside',
#                                 group == 4 & (y < repro_draws_list_firstyear$life_fldg$ci_95['2.5%'] | y > repro_draws_list_firstyear$life_fldg$ci_95['97.5%']) ~ 'outside'))
# 
# 
# #  $fill_group <- ifelse(p_build$y >= 0,'Above','Below')
# #This is necessary to ensure that instead of trying to draw
# # 3 polygons, we're telling ggplot to draw six polygons
# firstyear_base_violin_plot_build$group1 <- with(firstyear_base_violin_plot_build,
#                                                interaction(factor(group),factor(fill_group)))
# 
# 
# firstyear_mod_color_vec <- setNames(c('#3f0b70', '#04ab80', '#d62828', '#3478e5'),
#                                    nm = names(mod_pred_list_firstyear))
# 
# firstyear_mod_color_post <- setNames(c('#c7b6d8', '#9addcc', '#eea9a9', '#adc9f4'),
#                                     nm = names(mod_pred_list_firstyear))
# 
# 
# 
# 
# firstyear_violin_post_plot <- firstyear_base_violin_plot_build %>% 
#   filter(fill_group == 'inside') %>% 
#   ggplot() + 
#   geom_violin(data = repro_draws_only_firstyear %>%
#                 mutate(mod = factor(mod, levels = rev(c("Lifetime\nrepro.\nsuccess", "Total\nnesting\nyears", "First\nnesting\nage", "Last\nnesting\nage")))), 
#               aes(x = mod, y = b_first_year_scaled, fill = mod),
#               color = NA, alpha = 0.8) +
#   scale_fill_manual(values = unname(rev(firstyear_mod_color_post))) +
#   ggnewscale::new_scale_fill() +
#   geom_polygon(aes(x = x,y = y, group = group1, fill = group1)) +
#   scale_fill_manual(values = unname(rev(firstyear_mod_color_vec))) +
#   theme(legend.position = 'none',
#         axis.line.x = element_line(color = 'black', linewidth = 0.5),
#         panel.background = element_blank(),
#         #panel.grid.major.y = element_line(color = '#d8d8d8', linetype = 'dashed', linewidth = 0.4),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(color = rev(firstyear_mod_color_vec), face = 'bold')) +
#   coord_flip() +
#   ylab('First breeding year coefficient') +
#   #scale_x_discrete(name = NULL, sec.axis = sec_axis(~., name = "Y-axis Label on right side")) +
#   guides(y = "none") +
#   guides(y.sec = guide_axis_label_trans(~.x)) +
#   geom_hline(yintercept = 0, linetype = 'dashed', color = 'black')
# 
# 
# 
# mod_pred_list_firstyear$nestyear_mod  %>% 
#   group_by(first_year_scaled) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop')
# 
# mod_plot_list_firstyear <- lapply(setNames(nm = names(mod_pred_list_firstyear)), 
#                                  function(PRED_OBJ, mod, color_vec, color_post) {
#                                    mod[[PRED_OBJ]] %>% 
#                                      group_by(first_year_scaled) %>% 
#                                      summarize(lower25 = quantile(.epred, prob = 0.025),
#                                                lower10 = quantile(.epred, prob = 0.10),
#                                                lower20 = quantile(.epred, prob = 0.20),
#                                                median = quantile(.epred, prob = 0.5),
#                                                upper80 = quantile(.epred, prob = 0.80),
#                                                upper90 = quantile(.epred, prob = 0.90),
#                                                upper975 = quantile(.epred, prob = 0.975),
#                                                .groups = 'drop') %>% 
#                                      ggplot(aes(x = first_year_scaled, y = median)) + 
#                                      geom_line(data = mod[[PRED_OBJ]],
#                                                aes(x = first_year_scaled, y = `.epred`, group = .draw),
#                                                linewidth = 0.25, 
#                                                color = color_post[PRED_OBJ], 
#                                                alpha = 0.4) +
#                                      geom_ribbon(aes(first_year_scaled,
#                                                      ymin = lower25, ymax = upper975),
#                                                  color = color_vec[PRED_OBJ], linetype = 'dashed', linewidth = 1.5,
#                                                  lineend = 'round', fill = NA,
#                                                  #fill = '#3f0b70', #'#c7b6d8',
#                                                  alpha = 0.2) +
#                                      geom_line(linewidth = 3, color = color_vec[PRED_OBJ],
#                                                lineend = 'round') +
#                                      theme_bw() +
#                                      theme(panel.grid.minor = element_blank(),
#                                            panel.grid.major = element_line(linetype = 'dashed', 
#                                                                            linewidth = 0.4,
#                                                                            color = '#d8d8d8'))
#                                  }, mod = mod_pred_list_firstyear, color_vec = firstyear_mod_color_vec, color_post = firstyear_mod_color_post)
# 
# 
# repro_combined_plot_firstyear <- cowplot::plot_grid(mod_plot_list_firstyear$life_fldg +
#                                             theme(axis.title.x = element_blank(),
#                                                   panel.border = element_blank(),
#                                                   axis.line = element_line(linewidth = 0.5)) +
#                                             ylab('Life reproductive success'),
#                                           mod_plot_list_firstyear$nestyear_mod +
#                                             theme(axis.title.x = element_blank(),
#                                                   panel.border = element_blank(),
#                                                   axis.line = element_line(linewidth = 0.5)) +
#                                             ylab('Total nesting years'),
#                                           mod_plot_list_firstyear$firstbreed +
#                                             theme(axis.title.x = element_blank(),
#                                                   panel.border = element_blank(),
#                                                   axis.line = element_line(linewidth = 0.5)) +
#                                             ylab('First nesting age'),
#                                           mod_plot_list_firstyear$lastbreed +
#                                             theme(axis.title.x = element_blank(),
#                                                   panel.border = element_blank(),
#                                                   axis.line = element_line(linewidth = 0.5)) +
#                                             ylab('Last nesting age'),
#                                           ncol = 2,
#                                           align = 'hv',
#                                           axis = 'lr',
#                                           labels = c('A', 'B', 'C', 'D'),
#                                           label_fontface = 'bold')
# 
# 
# #cowplot::plot_grid(repro_combined_plot,
# #                   mod_coefficient_sina,
# #                   ncol = 2, rel_widths = c(0.8, 0.2))
# 
# #library(grid)
# #library(gridExtra)
# 
# x_title_firstyear <- textGrob("First breeding year", 
#                     gp=gpar(fontface="plain", col="black", fontsize = 12))
# 
# #add to plot
# 
# firstyear_repro_mod_multipan <- cowplot::plot_grid(
#   grid.arrange(arrangeGrob(repro_combined_plot_firstyear, 
#                            bottom = textGrob("First breeding year", 
#                                              gp=gpar(fontface="plain", col="black", fontsize = 12))
#   )),
#   grid.arrange(arrangeGrob(firstyear_violin_post_plot + #mod_coefficient_sina +
#                              theme(axis.title.x = element_blank()), 
#                            bottom = textGrob("First breeding year coefficent", 
#                                              gp=gpar(fontface="plain", col="black", fontsize = 12))
#   )),
#   ncol = 2, rel_widths = c(0.6, 0.2),
#   labels = c('', 'E'), label_fontface = 'bold')
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'firstyear_repro_mod_multipan.png' ),
#                  plot = firstyear_repro_mod_multipan,
#                  width = 9*1.25, height = 5*1.25, bg = 'white')
# 


# ci_95 <- quantile(as_draws_df(nestyear_mod)$b_transloc_anc, 
#                   prob = c(0.025, 0.975))
# 
# 
# test_ridge <- as_draws_df(nestyear_mod) %>% 
#   select(b_transloc_anc) %>% 
#   mutate(in_ci = if_else(b_transloc_anc < ci_95[1] | b_transloc_anc > ci_95[2], 'out', 'in'),
#          group = 'one') %>% 
#   ggplot(aes(x = b_transloc_anc, y = group, group = group,
#              color = b_transloc_anc)) + 
#   ggridges::geom_density_ridges2(fill = NA,
#                                  #color = NA,
#                                  jittered_points = TRUE,
#                                  aes(point_color = b_transloc_anc)
#                                  ) +
#   ggridges::scale_point_color_gradient(#midpoint = 0,
#                                        low = '#c5625b', 
#                                        #mid = '#f2f2f2', 
#                                        high = '#326934') +
#   scale_color_gradient2(midpoint = 0,
#                         low = '#c5625b', 
#                         mid = '#f2f2f2', 
#                         high = '#326934'#,
#                         #limits = c(lower_lim_sina, upper_lim_sina)
#   ) +
#   theme_void() +
#   theme(legend.position = 'none')
# 
# 
# as_draws_df(mod_fit_life_fldg1[[31]])
# 
# beta_draws <- as_draws_df(mod_fit_life_fldg1[[31]])
# 
# ci_95 <- quantile(beta_draws$b_transloc_anc,
#                   prob = c(0.025, 0.975))
# 
# str_replace(as.character(mod_fit_life_fldg1[[31]]$formula$formula)[2], '| cens(censored)', '')
# 
# repro_mod_list <- setNames(list(mod_fit_life_fldg1[[31]], nestyear_mod, firstbreed_mod, lastbreed_mod),
#                            nm = c('life_fldg', 'nestyear_mod', 'firstbreed', 'lastbreed'))
# repro_draws_list <- lapply(setNames(nm = names(repro_mod_list)), 
#                           function(MOD, mod_list) {
#                             output_list <- list()
#                             #message('starting ')
#                             beta_draws <- as_draws_df(mod_list[[MOD]])
#                             
#                             ci_95 <- quantile(beta_draws$b_transloc_anc,
#                                               prob = c(0.025, 0.975))
#                             output_list[['ci_95']] <- ci_95
#                             output_list[['draws']] <- beta_draws %>% 
#                               select(b_transloc_anc) %>% 
#                               mutate(in_ci = if_else(b_transloc_anc < ci_95[1] | b_transloc_anc > ci_95[2], 'out', 'in'),
#                                      group = 'one') %>% 
#                               mutate(mod = case_when(MOD == 'life_fldg' ~ 'Lifetime\nrepro.\nsuccess',
#                                                      MOD == 'nestyear_mod' ~ 'Total\nnesting\nyears',
#                                                      MOD == 'firstbreed' ~ 'First\nnesting\nage',
#                                                      MOD == 'lastbreed' ~ 'Last\nnesting\nage'))
#                             
#                             return(output_list)
#                           }, mod_list = repro_mod_list) 
# 
# repro_draws_only <- lapply(repro_draws_list, function(x) x$draws) %>% 
#   bind_rows()
# 
# 
# 
# 
# 
# mod_coefficient_sina <- repro_draws_list %>%
#   mutate(mod = factor(mod, levels = rev(c("Lifetime\nrepro.\nsuccess", "Total\nnesting\nyears", "First\nnesting\nage", "Last\nnesting\nage")))) %>% 
#   ggplot() +
#   geom_hline(yintercept = 0, linewidth = 0.5) +
#   ggforce::geom_sina(aes(x = mod, 
#                          y = b_transloc_anc,
#                          color = b_transloc_anc,
#                          group = mod),
#                      fill = 'white',
#                      shape = 21,
#                      stroke = 0.4,
#                      method = "density",
#   ) +
#   ggforce::geom_sina(aes(x = mod, y = b_transloc_anc,
#                          fill = b_transloc_anc,
#                          color = b_transloc_anc,
#                          alpha = in_ci,
#                          group = mod),
#                      shape = 21,
#                      stroke = 0.4,
#                      #position = 'jitter',
#                      method = "density",
#                      #color = NA
#   ) +
#   scale_alpha_manual(values = c(1, 0)) +
#   scale_fill_gradient2(midpoint = 0,
#                         low = '#c5625b', 
#                         mid = '#dadddb', 
#                         high = '#326934'#,
#                         #limits = c(lower_lim_sina, upper_lim_sina)
#   ) +
#   scale_color_gradient2(midpoint = 0,
#                        low = '#c5625b',
#                        mid = '#dadddb',
#                        high = '#326934'#,
#                        #limits = c(lower_lim_sina, upper_lim_sina)
#   ) +
#   #theme_void() +
#   theme(legend.position = 'none',
#         axis.line.x = element_line(color = 'black', linewidth = 0.5),
#         panel.background = element_blank(),
#         panel.grid.major.y = element_line(color = '#d8d8d8', linetype = 'dashed', linewidth = 0.4),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(color = rev(transloc_mod_color_vec), face = 'bold')) +
#   coord_flip() +
#   ylab('Translocation coefficient') +
#   #scale_x_discrete(name = NULL, sec.axis = sec_axis(~., name = "Y-axis Label on right side")) +
#   guides(y = "none") +
#   guides(y.sec = guide_axis_label_trans(~.x))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# p <- ggplot() + 
#   geom_violin(data = dat,aes(x = factor(x),y = y))
# p_build <- ggplot2::ggplot_build(p)$data[[1]]
# 
# p_build <- transform(p_build,
#                      xminv = x - violinwidth * (x - xmin),
#                      xmaxv = x + violinwidth * (xmax - x))
# 
# p_build <- rbind(plyr::arrange(transform(p_build, x = xminv), y),
#                  plyr::arrange(transform(p_build, x = xmaxv), -y))
# 
# 
# p_build$fill_group <- ifelse(p_build$y >= 0,'Above','Below')
# #This is necessary to ensure that instead of trying to draw
# # 3 polygons, we're telling ggplot to draw six polygons
# p_build$group1 <- with(p_build,interaction(factor(group),factor(fill_group)))
# 
# repro_draws_only %>%
#   mutate(mod = factor(mod, levels = rev(c("Lifetime\nrepro.\nsuccess", "Total\nnesting\nyears", "First\nnesting\nage", "Last\nnesting\nage")))) %>% 
# 
# ggplot() +
#   geom_violin(data = dat, aes(x = x, y = y, group = x, fill = as.factor(x))) +
#   scale_fill_manual(values = c('red', 'purple', 'white')) +
#   ggnewscale::new_scale_fill() +
#   geom_polygon(data = p_build %>% 
#                  filter(grepl('Below', group1)) %>% 
#                  mutate(fill_val = case_when(group == 1 ~ 'gray',
#                                          group == 2 ~ 'black',
#                                          group == 3 ~ 'blue')),
#                aes(x = x,y = y,group = group1, fill = fill_val)) +
#   #scale_fill_identity() #+
#   scale_fill_manual(values = c('gray', 'black', 'blue') )
# 
# 
# repro_sina_list <- lapply(setNames(list(mod_fit_life_fldg1[[31]], nestyear_mod, firstbreed_mod, lastbreed_mod),
#                                    nm = c('life_fldg', 'nestyear_mod', 'firstbreed', 'lastbreed')), 
#                           function(MOD) {
#   ci_95 <- quantile(as_draws_df(MOD)$b_transloc_anc,
#                     prob = c(0.025, 0.975))
#   as_draws_df(MOD) %>% 
#     select(b_transloc_anc) %>% 
#     mutate(in_ci = if_else(b_transloc_anc < ci_95[1] | b_transloc_anc > ci_95[2], 'out', 'in'),
#            group = 'one') %>% 
#     ggplot() +
#     ggforce::geom_sina(aes(x = 'a', y = b_transloc_anc,
#                            color = b_transloc_anc,
#                            alpha = in_ci,
#                            size = in_ci,
#                            group = 'single'),
#                        #position = 'jitter',
#                        method = "density"
#     ) +
#     scale_size_manual(values = c(0.5, 0.3)) +
#     scale_alpha_manual(values = c(1, 0.1)) +
#     scale_color_gradient2(midpoint = 0,
#                           low = '#c5625b', 
#                           mid = '#f2f2f2', 
#                           high = '#326934'#,
#                           #limits = c(lower_lim_sina, upper_lim_sina)
#     ) +
#     # geom_hline(data = data.frame(val = quantile(as_draws_df(nestyear_mod)$b_transloc_anc, 
#     #                                             prob = c(0.025, 0.975))),
#     #            aes(yintercept = val,
#     #                color = ),
#     #            linewidth = 1.5) +
#     theme_void() +
#     theme(legend.position = 'none')
# })
#   
# test_sina <- as_draws_df(nestyear_mod) %>% 
#   select(b_transloc_anc) %>% 
#   mutate(in_ci = if_else(b_transloc_anc < ci_95[1] | b_transloc_anc > ci_95[2], 'out', 'in'),
#          group = 'one') %>% 
#   ggplot() +
#   ggforce::geom_sina(aes(x = 'a', y = b_transloc_anc,
#                          color = b_transloc_anc,
#                          alpha = in_ci,
#                          size = in_ci,
#                          group = 'single'),
#                      #position = 'jitter',
#                      method = "density"
#                      ) +
#   scale_size_manual(values = c(0.5, 0.3)) +
#   scale_alpha_manual(values = c(1, 0.1)) +
#   scale_color_gradient2(midpoint = 0,
#                         low = '#c5625b', 
#                         mid = '#f2f2f2', 
#                         high = '#326934'#,
#                         #limits = c(lower_lim_sina, upper_lim_sina)
#                         ) +
#   # geom_hline(data = data.frame(val = quantile(as_draws_df(nestyear_mod)$b_transloc_anc, 
#   #                                             prob = c(0.025, 0.975))),
#   #            aes(yintercept = val,
#   #                color = ),
#   #            linewidth = 1.5) +
#   theme_void() +
#   theme(legend.position = 'none')
# 
# 
# ylab_vec <- setNames(c('Lifetime reproductive success', 'Nesting years', 'First breeding age', 'Last breeding age'),
#                      nm = names(mod_plot_list_transloc))
# repro_success_plot_inset_list <- lapply(names(mod_plot_list_transloc), function(MOD, main_plot, sina, ylab_title) {
#   
#   cowplot::ggdraw() +
#     cowplot::draw_plot(
#       main_plot[[MOD]] +
#         theme(axis.title.x = element_blank(),
#               panel.border = element_blank(),
#               axis.line = element_line(linewidth = 1),
#               panel.grid.major = element_blank()) +
#         ylab(ylab_title[MOD])
#     ) +
#     cowplot::draw_plot(sina[[MOD]] + coord_flip(), 
#                        x = 0.37, #ifelse(MOD == 'firstbreed', 0.85, 0.2), 
#                        y = .8, 
#                        hjust = 0.5,
#                        width = .35, height = .2)
#   
# }, main_plot = mod_plot_list_transloc, sina = repro_sina_list, ylab_title = ylab_vec)
# 
# 
# x_title <- textGrob("Expected translocation ancestry", 
#                     gp=gpar(fontface="plain", col="black", fontsize = 12))
# 
# #add to plot
# repro_combined_plot <- cowplot::plot_grid(repro_success_plot_inset_list[[1]],
#                                           repro_success_plot_inset_list[[2]],
#                                           repro_success_plot_inset_list[[3]],
#                                           repro_success_plot_inset_list[[4]],
#                                           nrow = 2,
#                                           align = 'hv',
#                                           axis = 'lr')
# grid.arrange(arrangeGrob(repro_combined_plot, bottom = x_title))
# 


# mod_pred_list_transloc$life_fldg %>% 
#   group_by(transloc_anc) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = transloc_anc, y = median)) + 
#   geom_line(data = mod_pred_list_transloc$life_fldg,
#             aes(x = transloc_anc, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(transloc_anc,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed'))


# repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(transloc_anc = seq_range(first_year_scaled, n = 100),
#             mean_group_size = quantile(mean_group_size, prob = 0.5),
#             sex = 'mid',
#             first_year_scaled = quantile(first_year_scaled, prob = 0.5),
#             anc_count_z = quantile(anc_count_z, prob = 0.5)
#   ) %>% 
#   add_epred_draws(nestyear_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# 
# nestyear_mod_pred_first_year_scaled <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(transloc_anc = quantile(transloc_anc, prob = 0.5),
#             first_year_scaled = seq_range(first_year_scaled, n = 100)) %>% 
#   add_epred_draws(nestyear_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# 
# plot_nestyear_mod_pred_transloc_anc <- nestyear_mod_pred_transloc_anc %>% 
#   group_by(transloc_anc) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = transloc_anc, y = median)) + 
#   geom_line(data = nestyear_mod_pred_transloc_anc,
#             aes(x = transloc_anc, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(transloc_anc,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Expected translocated ancestry') + 
#   ylab('Number of nesting years')



# repro_success_df_sexpivot1 %>% 
#   filter(birth_during_monitoring == 'yes') %>% 
#   #filter(age != 0) %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
#   ggplot() +
#   geom_histogram(aes(mean_fldg_per_year), bins = 15)
# 


# ### LIFETIME REPRODUCTIVE SUCCESS MODEL ###
# 
# life_fldg_mod <- brm(
#   life_fldg | cens(censored) ~ mean_group_size + transloc_anc + first_year_scaled + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000
# )
# 
# life_fldg_formula_list <- list(
#   c('mean_group_size', 'first_year_scaled', 'sex'),
#   c('mean_group_size', 'first_year_scaled', 'sex', 'transloc_anc'),
#   c('mean_group_size', 'first_year_scaled', 'sex', 'anc_count_z'),
#   c('mean_group_size', 'first_year_scaled', 'sex', 'transloc_anc', 'anc_count_z')
# )
# 
# 
# 
# life_fldg_formula_list1 <- lapply(life_fldg_formula_list, function(x) paste0('life_fldg | cens(censored) ~ ', paste(x, collapse = ' + '), ' + (1|gr(dummy,cov=A))'))
# 
# mod_fit_life_fldg <- lapply(life_fldg_formula_list1, function(x, dat, add_rel) {
#   message('finish')
#   
#   return(
#     brm(
#       as.formula(x), 
#       data = dat, 
#       data2=list(A = add_rel),
#       family = poisson(link = "log"), 
#       #cov_ranef = list(A = rcw_addrel),
#       chains = 3, cores = 2, iter = 10000,
#       save_pars = save_pars(all = TRUE)
#     )
#   )
#   
# }, dat = repro_success_df_sexpivot1 %>%
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes') %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none')),
# add_rel = rcw_addrel
# )
# 
# 
# 
# mod_fit_life_fldg1 <- lapply(mod_fit_life_fldg, function(x) {
#   message('finished')
#   return(
#     add_criterion(x, c("waic"))
#   )
# })
# 
# mod_fit_life_fldg1[[1]]
# 
# 
# model_simple2.1_censored <- add_criterion(model_simple2.1_censored, 
#                                           c("loo"), moment_match = TRUE)
# 
# 
# 
# 
# 


#!/usr/bin/Rscript

#####################
### SCRIPT SET-UP ###
#####################

# array_index <- as.numeric(commandArgs(trailingOnly=TRUE)[[1]])
# 
# ### paths ###
# output_path <- "/mnt/home/lewanski/rcw_project1/output_results/life_fldg_mods/"
# input_path <- "/mnt/home/lewanski/rcw_project1/input_data/"
# 
# ### packages ##
# library(brms)
# 
# 
# ### load data and info ###
# predictor_list <- readRDS()
# 
# #data
# dat <- readRDS()
# 
# #relationship matrix
# rcw_addrel <- readRDS()
# 
# 
# #predictor list
# 
# focal_mod <- predictor_list[array_index]
# focal_mod_name <- names(predictor_list)[array_index]
# 
# mod_object <- add_criterion(
#   brm(
#     as.formula(x), 
#     data = repro_success_df_sexpivot1 %>%
#       #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#       filter(birth_during_monitoring == 'yes') %>%
#       mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#     data2 = list(A = rcw_addrel),
#     family = poisson(link = "log"), 
#     chains = 4, iter = 8000
#     ),
#   c("loo"), moment_match = TRUE)
# 
# 
# saveRDS(mod_object, file = paste0(output_path, focal_mod_name, '.rds'))
# 
# 
# 
# 
# 
# 
# mod_fit <- print(
#   loo_compare(mod_fit_life_fldg1[[1]],
#               mod_fit_life_fldg1[[2]],
#               mod_fit_life_fldg1[[3]],
#               mod_fit_life_fldg1[[4]],
#               criterion = "waic"),
#   simplify = FALSE
# )
# 
# 
# 
# 
# 
# loo_compare(mod_fit_life_fldg1,
#             criterion = "waic")
# 
# 
# length(mod_fit_life_fldg1)
# 
# 
# 
# 
# 
# life_fldg_mod2 <- brm(
#   life_fldg | cens(censored) ~ mean_group_size + transloc_anc + first_year_scaled + compl_simpson_z + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000
# )
# 
# 
# 
# life_fldg_mod1 <- brm(
#   life_fldg | cens(censored) ~ mean_group_size + transloc_anc + first_year_scaled + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     filter(first_year < 2012) %>% 
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000,
#   save_pars = save_pars(all = TRUE)
# )
# 
# stancode(life_fldg_mod1)
# waic(life_fldg_mod)
# 
# test_loo <- loo_compare(add_criterion(life_fldg_mod, criterion = "loo"),
#             add_criterion(life_fldg_mod2, criterion = "loo"),
#             criterion = "loo")
# 
# 
# life_fldg_mod2
# 
# 
# # loo_compare(loo(life_fldg_mod, moment_match = TRUE),
# #             loo(life_fldg_mod1, moment_match = TRUE)
# #             )
# 
# 
# 
# 
# 
# ### MODELS EXPLORING WHY LRS IS RELATED TO TIME ###
# 
# #MODEL 1: NESTING YEARS
# nestyear_mod <- brm(
#   nesting_years | cens(censored) ~ first_year_scaled + transloc_anc + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>% 
#     filter(birth_during_monitoring == 'yes') %>% 
#     #filter(age != 0) %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.9)
# )
# #performance::check_collinearity(nestyear_mod)
# 
# 
# 
# #MODEL 2:MEAN FLEDGE/NEST
# mean_fldg_mod <- brm(
#   mean_fldg_per_year | cens(censored) ~ transloc_anc + first_year_scaled + (1|gr(dummy,cov=A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = hurdle_gamma(), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.97)
# )
# 
# 
# 
# #MODEL 3: MIN NESTING YEAR
# firstbreed_mod <- brm(
#   birth_firstbreed_dif ~ transloc_anc + birth_year  + sex + (1|gr(dummy,cov=A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
# 
# #MODEL 4: MAX NESTING YEAR
# lastbreed_mod <- brm(
#   birth_lastbreed_dif | cens(censored) ~ transloc_anc + birth_year + sex + (1|gr(dummy,cov=A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
# 
# lastbreed_mod1 <- brm(
#   birth_lastbreed_dif | cens(censored) ~ transloc_anc + birth_year + sex + (1|gr(dummy,cov=A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     filter(first_year <= 2012) %>% 
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
# performance::check_collinearity(lastbreed_mod)
# 
# 
# ###################################
# ### POSTERIOR PREDICTIVE CHECKS ###
# ###################################
# 
# 
# 
# 
# ###########################
# ### MODEL VISUALIZATION ###
# ###########################
# 
# ### MODEL 1:MEAN FLEDGE/NEST ###
# lrs_mod_pred_transloc_anc <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(transloc_anc = seq_range(transloc_anc, n = 100),
#             mean_group_size = quantile(mean_group_size, prob = 0.5),
#             first_year_scaled = quantile(first_year_scaled, prob = 0.5)) %>% 
#   add_epred_draws(life_fldg_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# lrs_mod_pred_mean_group_size <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(transloc_anc = quantile(transloc_anc, prob = 0.5),
#             mean_group_size = seq_range(mean_group_size, n = 100),
#             first_year_scaled = quantile(first_year_scaled, prob = 0.5)) %>% 
#   add_epred_draws(life_fldg_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# lrs_mod_pred_first_year_scaled <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(transloc_anc = quantile(transloc_anc, prob = 0.5),
#             mean_group_size = quantile(mean_group_size, prob = 0.5),
#             first_year_scaled = seq_range(first_year_scaled, n = 100)) %>% 
#   add_epred_draws(life_fldg_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# 
# plot_lrs_mod_pred_transloc_anc <- lrs_mod_pred_transloc_anc %>% 
#   group_by(transloc_anc) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = transloc_anc, y = median)) + 
#   geom_line(data = lrs_mod_pred_transloc_anc,
#             aes(x = transloc_anc, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(transloc_anc,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Expected translocated ancestry') + 
#   ylab('Lifetime reproductive success')
# 
# 
# plot_lrs_mod_pred_mean_group_size <- lrs_mod_pred_mean_group_size %>% 
#   group_by(mean_group_size) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = mean_group_size, y = median)) + 
#   geom_line(data = lrs_mod_pred_mean_group_size,
#             aes(x = mean_group_size, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(mean_group_size,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Mean group size') + 
#   ylab('Lifetime reproductive success')
# 
# plot_lrs_mod_pred_first_year_scaled <- lrs_mod_pred_first_year_scaled %>% 
#   group_by(first_year_scaled) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = first_year_scaled, y = median)) + 
#   geom_line(data = lrs_mod_pred_first_year_scaled,
#             aes(x = first_year_scaled, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(first_year_scaled,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('First breeding year (shifted)') + 
#   ylab('Lifetime reproductive success')
# 
# 
# lrs_multipan <- cowplot::plot_grid(plot_lrs_mod_pred_transloc_anc,
#                    plot_lrs_mod_pred_mean_group_size,
#                    plot_lrs_mod_pred_first_year_scaled,
#                    nrow = 1)
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'lrs_multipan.png' ),
#                  plot = lrs_multipan,
#                  width = 10.5*1, height = 3*1, bg = 'white')
# 
# 
# 
# ### MODEL 2:MEAN FLEDGE/NEST ###
# nestyear_mod_pred_transloc_anc <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(transloc_anc = seq_range(transloc_anc, n = 100),
#             first_year_scaled = quantile(first_year_scaled, prob = 0.5)) %>% 
#   add_epred_draws(nestyear_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# nestyear_mod_pred_first_year_scaled <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(transloc_anc = quantile(transloc_anc, prob = 0.5),
#             first_year_scaled = seq_range(first_year_scaled, n = 100)) %>% 
#   add_epred_draws(nestyear_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# 
# plot_nestyear_mod_pred_transloc_anc <- nestyear_mod_pred_transloc_anc %>% 
#   group_by(transloc_anc) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = transloc_anc, y = median)) + 
#   geom_line(data = nestyear_mod_pred_transloc_anc,
#             aes(x = transloc_anc, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(transloc_anc,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Expected translocated ancestry') + 
#   ylab('Number of nesting years')
# 
# 
# plot_nestyear_mod_pred_first_year_scaled <- nestyear_mod_pred_first_year_scaled %>% 
#   group_by(first_year_scaled) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = first_year_scaled, y = median)) + 
#   geom_line(data = nestyear_mod_pred_first_year_scaled,
#             aes(x = first_year_scaled, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(first_year_scaled,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('First breeder year (shifted)') + 
#   ylab('Number of nesting years')
# 
# 
# nestyear_multipan <- cowplot::plot_grid(plot_nestyear_mod_pred_transloc_anc,
#                                         plot_nestyear_mod_pred_first_year_scaled,
#                                         nrow = 1)
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'nestyear_multipan.png' ),
#                  plot = nestyear_multipan,
#                  width = 7*1.15, height = 3*1.15, bg = 'white')
# 
# 
# 
# 
# 
# 
# 
# ### MODEL 3:MEAN FLEDGE/NEST ###
# mean_fldg_mod_pred_transloc_anc <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(transloc_anc = seq_range(transloc_anc, n = 100),
#             first_year_scaled = quantile(first_year_scaled, prob = 0.5)) %>% 
#   add_epred_draws(mean_fldg_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# mean_fldg_mod_pred_first_year_scaled <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(transloc_anc = quantile(transloc_anc, prob = 0.5),
#             first_year_scaled = seq_range(first_year_scaled, n = 100)) %>% 
#   add_epred_draws(mean_fldg_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# plot_mean_fldg_mod_pred_transloc_anc <- mean_fldg_mod_pred_transloc_anc %>% 
#   group_by(transloc_anc) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = transloc_anc, y = median)) + 
#   geom_line(data = mean_fldg_mod_pred_transloc_anc,
#             aes(x = transloc_anc, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(transloc_anc,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Expected translocated ancestry') + 
#   ylab('Mean(fledgling/nesting year)')
# 
# 
# plot_mean_fldg_mod_pred_first_year_scaled <- mean_fldg_mod_pred_first_year_scaled %>% 
#   group_by(first_year_scaled) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = first_year_scaled, y = median)) + 
#   geom_line(data = mean_fldg_mod_pred_first_year_scaled,
#             aes(x = first_year_scaled, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(first_year_scaled,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('First breeder year (shifted)') + 
#   ylab('Mean(fledgling/nesting year)')
# 
# mean_fldg_multipan <- cowplot::plot_grid(plot_mean_fldg_mod_pred_transloc_anc,
#                                          plot_mean_fldg_mod_pred_first_year_scaled,
#                                    nrow = 1)
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'mean_fldg_multipan.png' ),
#                  plot = mean_fldg_multipan,
#                  width = 7*1.15, height = 3*1.15, bg = 'white')
# 
# 
# 
# 
# 
# ### MODEL 3: MIN NESTING YEAR ###
# model_firstbreed_pred_birthyear <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(birth_year = seq_range(birth_year, n = 100),
#             transloc_anc = quantile(transloc_anc, prob = 0.5),
#             sex = 'mid') %>% 
#   add_epred_draws(firstbreed_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# model_firstbreed_pred_transloc <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(birth_year = quantile(birth_year, prob = 0.5),
#             transloc_anc = seq_range(transloc_anc, n = 100),
#             sex = 'mid') %>% 
#   add_epred_draws(firstbreed_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# 
# plot_model_firstbreed_pred_birthyear <- model_firstbreed_pred_birthyear %>% 
#   group_by(birth_year) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = birth_year, y = median)) + 
#   geom_line(data = model_firstbreed_pred_birthyear,
#             aes(x = birth_year, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(birth_year,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Birth year') + 
#   ylab('Age at first nesting year')
# 
# 
# plot_model_firstbreed_pred_transloc <- model_firstbreed_pred_transloc %>% 
#   group_by(transloc_anc) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = transloc_anc, y = median)) + 
#   geom_line(data = model_firstbreed_pred_transloc,
#             aes(x = transloc_anc, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(transloc_anc,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Expected translocation ancestry') +
#   ylab('Age at first nesting year')
# 
# firstbreed_multipan <- cowplot::plot_grid(plot_model_firstbreed_pred_birthyear,
#                                           plot_model_firstbreed_pred_transloc,
#                                          nrow = 1)
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'firstbreed_multipan.png' ),
#                  plot = firstbreed_multipan,
#                  width = 7*1.15, height = 3*1.15, bg = 'white')
# 
# 
# 
# ### MODEL 4: MAX NESTING YEAR ###
# model_lastbreed_pred_birthyear <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(birth_year = seq_range(birth_year, n = 100),
#             transloc_anc = quantile(transloc_anc, prob = 0.5),
#             sex = 'mid') %>% 
#   add_epred_draws(lastbreed_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# model_lastbreed_pred_transloc <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(birth_year = quantile(birth_year, prob = 0.5),
#             transloc_anc = seq_range(transloc_anc, n = 100),
#             sex = 'mid') %>% 
#   add_epred_draws(lastbreed_mod,
#                   ndraws = 1000,
#                   re_formula = NA)
# 
# 
# plot_model_lastbreed_pred_birthyear <- model_lastbreed_pred_birthyear %>% 
#   group_by(birth_year) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = birth_year, y = median)) + 
#   geom_line(data = model_lastbreed_pred_birthyear,
#             aes(x = birth_year, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(birth_year,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Birth year') + 
#   ylab('Age at last nesting year')
# 
# 
# plot_model_lastbreed_pred_transloc <- model_lastbreed_pred_transloc %>% 
#   group_by(transloc_anc) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = transloc_anc, y = median)) + 
#   geom_line(data = model_lastbreed_pred_transloc,
#             aes(x = transloc_anc, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(transloc_anc,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round', fill = NA,
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Expected translocation ancestry') + 
#   ylab('Age at last nesting year')
# 
# lastbreed_multipan <- cowplot::plot_grid(plot_model_lastbreed_pred_birthyear,
#                                           plot_model_lastbreed_pred_transloc,
#                                           nrow = 1)
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'lastbreed_multipan.png' ),
#                  plot = lastbreed_multipan,
#                  width = 7*1.15, height = 3*1.15, bg = 'white')

# library(brms)
# library(loo)
# 
# model_simple2 <- brm(
#   life_fldg ~ mean_group_size, 
#   data = repro_success_df_sexpivot1 %>%
#     filter(alive == 'no' & birth_during_monitoring == 'yes'), 
#   #data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   chains = 2, cores = 2, iter = 1000
# )
# 
# stancode(model_simple2)
# 
# model_simple2.1 <- brm(
#   life_fldg ~ mean_group_size + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     filter(alive == 'no' & birth_during_monitoring == 'yes'), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 1000
# )
#
# model_simple2.1 <- brm(
#   breeder ~ transloc_anc + min_year_scaled + Sex + (1|gr(RCWid,cov=A)), 
#   data = processed_breed_info1 %>%
#     filter(birth_during_monitoring == 'yes' & alive == 'no') %>%
#     mutate(min_year_scaled = min_year - min(min_year)) %>% 
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = bernoulli(link = "logit"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.9)
# )
#
# life_mod <- brm(
#   age | cens(censored) ~ first_year_scaled + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>% 
#     filter(birth_during_monitoring == 'yes') %>% 
#     filter(age != 0) %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.9)
# )
#
# nestcount_mod <- brm(
#   nesting_years | cens(censored) ~ first_year_scaled + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>% 
#     filter(birth_during_monitoring == 'yes') %>% 
#     filter(age != 0) %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.9)
# )
#
# model_nesting_count <- brm(
#   nesting_years | cens(censored) ~ transloc_anc + first_year_scaled + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 10000,
#   control = list(adapt_delta = 0.9)
# )
#
# model_mean_fldg <- brm(
#     mean_fldg_per_year | cens(censored) ~ transloc_anc + first_year_scaled + (1|gr(dummy,cov=A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = hurdle_gamma(), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
#
# model_firstbreed <- brm(
#   birth_firstbreed_dif ~ transloc_anc + birth_year + (1|gr(dummy,cov=A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
#
# model_lastbreed <- brm(
#   birth_lastbreed_dif | cens(censored) ~ transloc_anc + birth_year + sex + (1|gr(dummy,cov=A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
#
# model_lastbreed_2012 <- brm(
#   birth_lastbreed_dif ~ transloc_anc + birth_year + sex + (1|gr(dummy,cov=A)),
#   data = repro_success_df_sexpivot1 %>%
#     filter(first_year <= 2012) %>% 
#     filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     #filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )
#
# repro_success_df_sexpivot1 %>%
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes') %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
#   ggplot() +
#   geom_histogram(aes(nesting_years))

# model_simple2.1_censored <- brm(
#   life_fldg | cens(censored) ~ mean_group_size + transloc_anc + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 2, cores = 2, iter = 5000
# )
# 
# repro_success_df_sexpivot1 %>%
#   filter(alive == 'no' & birth_during_monitoring == 'yes') %>% 
#   nrow()
# 
# repro_success_df_sexpivot1 %>%
#   #filter(birth_during_monitoring == 'yes') %>% 
#   nrow()
# 
# model_simple3.1_censored <- brm(
#   life_fldg | cens(censored) ~ mean_group_size + transloc_anc + first_year_scaled + anc_count + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000
# )
# 
# model_simple3.3_censored <- brm(
#   life_fldg | cens(censored) ~ mean_group_size + transloc_anc + first_year_scaled + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000
# )
# 
# model_simple3.2_censored <- brm(
#   life_fldg | cens(censored) ~ mean_group_size + first_year_scaled + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000
# )
# model_binom <- brm(
#   life_fldg | cens(censored) ~ mean_group_size + transloc_anc + first_year_scaled + (1|gr(dummy,cov=A)), 
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none')), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = "log"), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000
# )

#library(tidybayes)
# test_pred <- repro_success_df_sexpivot1 %>%
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes') %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none')) %>% 
#   data_grid(transloc_anc = seq_range(transloc_anc, n = 100),
#             mean_group_size = mean(mean_group_size),
#             first_year_scaled = quantile(first_year_scaled , prob = 0.5)) %>% 
#   add_epred_draws(model_simple3.3_censored,
#                   ndraws = 1000,
#                   re_formula = NA)

# test_pred %>% 
#   group_by(transloc_anc) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = transloc_anc, y = median)) + 
#   geom_line(data = test_pred,
#             aes(x = transloc_anc, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   # geom_ribbon(aes(transloc_anc,
#   #                 ymin = lower25, ymax = upper975),
#   #             fill = '#ece6f2', alpha = 0.2) +
#   #geom_ribbon(aes(transloc_anc,
#   #                ymin = lower10, ymax = upper90),
#   #            fill = '#c7b6d8', color = '#c7b6d8', alpha = 0.2) +
#   #geom_ribbon(aes(transloc_anc,
#   #                ymin = lower25, ymax = upper975),
#   #            #color = '#c7b6d8', linetype = 'dashed', linewidth = 2.5
#   #            fill = '#3f0b70', #'#c7b6d8', 
#   #            alpha = 0.2) +
#   # geom_ribbon(aes(transloc_anc,
#   #                 ymin = lower10, ymax = upper90),
#   #             color = '#c7b6d8', fill = NA, 
#   #             linetype = 'solid', linewidth = 2) +
#   # #geom_ribbon(aes(transloc_anc,
#   #                ymin = lower20, ymax = upper80),
#   #            fill = '#a386be', alpha = 0.2) +
#   geom_ribbon(aes(transloc_anc,
#                 ymin = lower25, ymax = upper975),
#             color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#             lineend = 'round',
#             #fill = '#3f0b70', #'#c7b6d8',
#             alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Translocated ancestry') + 
#   ylab('Lifetime reproductive success')
#   #geom_line(aes(group = `.draw`)) +
#   #stat_lineribbon()

# model_firstbreed_pred <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(birth_year = seq_range(birth_year, n = 100),
#             transloc_anc = quantile(transloc_anc, prob = 0.5),
#             sex = 'mid') %>% 
#   add_epred_draws(model_firstbreed,
#                   ndraws = 1000,
#                   re_formula = NA)

# model_lastbreed_2012_pred <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(first_year <= 2012) %>% 
#   filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(birth_year = seq_range(birth_year, n = 100),
#             transloc_anc = quantile(transloc_anc, prob = 0.5),
#             sex = 'mid') %>% 
#   add_epred_draws(model_lastbreed_2012,
#                   ndraws = 1000,
#                   re_formula = NA)

# model_lastbreed_pred <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge)) %>% 
#   data_grid(birth_year = seq_range(birth_year, n = 100),
#             transloc_anc = quantile(transloc_anc, prob = 0.5),
#             sex = 'mid') %>% 
#   add_epred_draws(model_lastbreed,
#                   ndraws = 1000,
#                   re_formula = NA)

# model_firstbreed_pred %>% 
#   group_by(birth_year) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = birth_year, y = median)) + 
#   geom_line(data = model_firstbreed_pred,
#             aes(x = birth_year, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(birth_year,
#                  ymin = lower25, ymax = upper975),
#              color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#              lineend = 'round',
#              #fill = '#3f0b70', #'#c7b6d8',
#              alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Birth year') + 
#   ylab('Age at first nesting year') +
#   ylim(0, 3.5)
#
# model_lastbreed_2012_pred %>% 
#   group_by(birth_year) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = birth_year, y = median)) + 
#   geom_line(data = model_lastbreed_2012_pred,
#             aes(x = birth_year, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(birth_year,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round',
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
#   geom_line(linewidth = 3, color = '#3f0b70',
#             lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Birth year') + 
#   ylab('Age at final nesting year') +
#   ylim(0, 6.5)
#
# model_lastbreed_pred %>% 
#   group_by(birth_year) %>% 
#   summarize(lower25 = quantile(.epred, prob = 0.025),
#             lower10 = quantile(.epred, prob = 0.10),
#             lower20 = quantile(.epred, prob = 0.20),
#             median = quantile(.epred, prob = 0.5),
#             upper80 = quantile(.epred, prob = 0.80),
#             upper90 = quantile(.epred, prob = 0.90),
#             upper975 = quantile(.epred, prob = 0.975),
#             .groups = 'drop') %>% 
#   ggplot(aes(x = birth_year, y = median)) + 
#   geom_line(data = model_lastbreed_pred,
#             aes(x = birth_year, y = `.epred`, group = .draw),
#             linewidth = 0.24, color = '#c7b6d8', alpha = 0.6) +
#   geom_ribbon(aes(birth_year,
#                   ymin = lower25, ymax = upper975),
#               color = '#3f0b70', linetype = 'dashed', linewidth = 1.5,
#               lineend = 'round',
#               #fill = '#3f0b70', #'#c7b6d8',
#               alpha = 0.2) +
# geom_line(linewidth = 3, color = '#3f0b70',
#           lineend = 'round') +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linetype = 'dashed')) +
#   xlab('Birth year') + 
#   ylab('Age at final nesting year') +
#   ylim(0, 6.5)
#
# dataset_explore_multipan <- GGally::ggpairs(repro_success_df_sexpivot1 %>% 
#                   #select(!dummy & !first_year) %>% 
#                   #select(!matches("_z$")) %>% 
#                   filter(birth_during_monitoring == 'yes') %>% 
#                   select(sex, mean_group_size, transloc_anc, fped, first_year_scaled, anc_count)
#                 )
# 
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'dataset_explore_multipan.png' ),
#                  plot = dataset_explore_multipan,
#                  width = 14*0.8, height = 10*0.8, bg = 'white')


# adults_notalive <- pop_dat %>% 
#   group_by(RCWid) %>% 
#   summarize(year = max(year),
#             #Origin_update = first(Origin_update)
#             ) %>% 
#   ungroup() %>% 
#   left_join(.,
#             rcw_processed_list$rcws %>% 
#               select(RCWid, MinAge, Sex),
#             by = 'RCWid') %>% 
#   filter(year > MinAge #& Origin_update == 'Local'
#          ) %>% 
#   filter(year < 2022) %>% 
#   pull(RCWid) %>% 
#   unique()

# library(pedigreemm)
# 
# df_pedigreemm <- repro_success_df_sexpivot1 %>% 
#   filter(alive == 'no' & full_grandparent_info == 'no') %>% 
#   as.data.frame()

# test_pedigreemm <- pedigreemm::pedigreemm(life_fldg ~ transloc_anc + 
#                                             sex + mean_group_size + 
#                                             first_year_scaled + 
#                                             #fped +
#                                             (1|dummy), 
#                        data = repro_success_df_sexpivot1 %>% 
#                          filter(alive == 'no' & full_grandparent_info == 'yes'), 
#                        family = poisson(link = "log"), 
#                        REML = TRUE,
#                        pedigree = list(dummy = rcw_ped_object)
#                        )
# 
# summary(pedigreemm::pedigreemm(life_fldg ~ transloc_anc + sex +
#                                  fped +
#                                  (1|dummy), 
#                                data = repro_success_df_sexpivot1 %>% 
#                                  filter(alive == 'no' & full_grandparent_info == 'yes'), 
#                                family = poisson(link = "log"), 
#                                REML = TRUE,
#                                pedigree = list(dummy = rcw_ped_object)
# ))
# 
# MuMIn::r.squaredGLMM(test_lifetime_repro)
# library(DHARMa)
# library(lme4)
# #library(ggeffects)
# 
# variable_vec <- c("fped_z",
#                   "transloc_anc_z",
#                   "mean_group_size_z",
#                   "sex",
#                   "first_year_scaled_z",
#                   "compl_simpson",
#                   "anc_count_z")
# 
# 
# variable_combo_list <- unlist(lapply(1:length(variable_vec), combn, 
#               x = variable_vec, simplify = FALSE), 
#        recursive = FALSE)
# 
# #variable_combo_list[[length(variable_combo_list) + 1]] <- "1"
# 
# variable_combo_list1 <- lapply(variable_combo_list, function(x) paste0('life_fldg ~ ', paste(x, collapse = ' + '), ' + (1|dummy)'))
# 
# 
# mod_list <- lapply(variable_combo_list1, function(x, dat, ped) {
#   message(x)
#   return(pedigreemm::pedigreemm(as.formula(x), 
#                                 data = dat, 
#                                 family = poisson(link = "log"), 
#                                 REML = TRUE,
#                                 pedigree = list(dummy = ped),
#                                 control = glmerControl(optimizer= "Nelder_Mead", #Nelder_Mead",
#                                                        boundary.tol = 1e-15,
#                                                        optCtrl = list(maxfun = 1000000))
#   ))
#   
# }, dat = repro_success_df_sexpivot1 %>% 
#   filter(alive == 'no' & birth_during_monitoring == 'yes'),
# ped = rcw_ped_object)
# 
# 
# mod_list_r2 <- lapply(variable_combo_list1, function(x, dat, ped) {
#   
#   mod <- pedigreemm::pedigreemm(as.formula(x), 
#                                 data = dat, 
#                                 family = poisson(link = "log"), 
#                                 REML = TRUE,
#                                 pedigree = list(dummy = ped),
#                                 control = glmerControl(optimizer= "Nelder_Mead", #Nelder_Mead",
#                                              boundary.tol = 1e-15,
#                                              optCtrl = list(maxfun = 1000000))
#                                 )
#   MuMIn::r.squaredGLMM(mod)
#   
#   
# }, dat = repro_success_df_sexpivot1 %>% 
#   filter(alive == 'no' & birth_during_monitoring == 'yes'),
# ped = rcw_ped_object)
# 
# 
# mod_summary_df <- data.frame(mod_id = 1:length(mod_list),
#            mod_formula = unlist(variable_combo_list1),
#            AICc = sapply(mod_list, MuMIn::AICc),
#            R2m = sapply(mod_list_r2, function(x) x['trigamma', 'R2m']),
#            R2c = sapply(mod_list_r2, function(x) x['trigamma', 'R2c']))
# 
# 
# summary_test_topmod <- summary(mod_list[[30]])
# coefs_vec <- rownames(summary_test_topmod$coefficients)
# 
# 
# predict_list <- list()
# 
# for (i in coefs_vec[!coefs_vec %in% c("(Intercept)", 'sexmid')]) {
#   predictor_range <- range(repro_success_df_sexpivot1[repro_success_df_sexpivot1$alive == 'no' & repro_success_df_sexpivot1$birth_during_monitoring == 'yes',i])
#   
#   predict_list[[i]][['range']] <- predictor_range
#   predict_list[[i]][['emmean']] <- as.data.frame(emmeans::emmeans(mod_list[[30]], 
#                                                                   i,
#                                                                   at = setNames(list(seq(predictor_range[1], 
#                                                                                          predictor_range[2], 
#                                                                                          length.out = 100)), 
#                                                                                 nm = i),
#                                                                   type = "response")
#   )
# }
# 
# transloc_anc_plot <- predict_list$transloc_anc[['emmean']] %>% 
#   ggplot() +
#   geom_point(data = repro_success_df_sexpivot1 %>%
#                filter(alive == 'no' & repro_success_df_sexpivot1$birth_during_monitoring == 'yes'),
#              aes(x = transloc_anc_z,
#                  y = life_fldg,
#                  color = first_year_scaled),
#              size = 1,
#              alpha = 0.6) +
#   geom_ribbon(aes(transloc_anc_z,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = transloc_anc_z, y = rate),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw()
# 
# mean_group_size_plot <- predict_list$mean_group_size_z[['emmean']] %>% 
#   ggplot() +
#   geom_point(data = repro_success_df_sexpivot1 %>%
#                filter(alive == 'no' & repro_success_df_sexpivot1$birth_during_monitoring == 'yes'),
#              aes(x = mean_group_size_z,
#                  y = life_fldg,
#                  color = first_year_scaled),
#              size = 1,
#              alpha = 0.6) +
# geom_ribbon(aes(mean_group_size_z,
#                 ymin = asymp.LCL, ymax = asymp.UCL),
#             fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = mean_group_size_z, y = rate),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw()
# 
# first_year_scaled_plot <- predict_list$first_year_scaled_z[['emmean']] %>% 
#   ggplot() +
#   geom_point(data = repro_success_df_sexpivot1 %>%
#                filter(alive == 'no' & repro_success_df_sexpivot1$birth_during_monitoring == 'yes'),
#              aes(x = first_year_scaled_z,
#                  y = life_fldg,
#                  color = first_year_scaled),
#              size = 1,
#              alpha = 0.6) +
#   geom_ribbon(aes(first_year_scaled_z,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = first_year_scaled_z, y = rate),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw()
#
# estimated_relationship_bestfit <- cowplot::plot_grid(transloc_anc_plot,
#                    mean_group_size_plot,
#                    first_year_scaled_plot,
#                    ncol = 1)
# 
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'estimated_relationship_bestfit.png' ),
#                  plot = estimated_relationship_bestfit,
#                  width = 7*0.8, height = 10*0.8, bg = 'white')


# life_fldg_pedmod <- pedigreemm::pedigreemm(life_fldg ~ 
#                          transloc_anc +
#                          mean_group_size +
#                          fped + 
#                          sex + 
#                          first_year_scaled +
#                          (1|dummy), 
#                        data = repro_success_df_sexpivot1 %>% 
#                          filter(alive == 'no'), 
#                        family = poisson(link = "log"), 
#                        REML = TRUE,
#                        pedigree = list(dummy = rcw_ped_object),
#                        control = glmerControl(optimizer="Nelder_Mead",
#                                               boundary.tol = 1e-10))
# 
# fped_mod <- pedigreemm::pedigreemm(life_fldg ~ 
#                          transloc_anc +
#                          mean_group_size +
#                          #fped + 
#                          sex +
#                          first_year_scaled +
#                          (1|dummy), 
#                        data = repro_success_df_sexpivot1 %>% 
#                          filter(alive == 'no' & birth_during_monitoring == 'yes'), #& full_grandparent_info == 'yes'), 
#                        family = poisson(link = "log"), 
#                        REML = TRUE,
#                        pedigree = list(dummy = rcw_ped_object),
#                        control = glmerControl(optimizer="Nelder_Mead",
#                                              boundary.tol = 1e-10))
# 
# 
# dharma_test_lifetime_repro <- DHARMa::simulateResiduals(fittedModel = fped_mod,
#                                                         plot = FALSE)
# 
# plot(dharma_test_lifetime_repro)
# 
# #without year: 1779.91
# #with year: 1701.984
# #without sex: 1701.984
# 
# test_lifetime_repro <- pedigreemm::pedigreemm(life_fldg ~ 
#                                  transloc_anc +
#                                  mean_group_size +
#                                  sex + 
#                                  first_year_scaled +
#                                  (1|dummy), 
#                                data = repro_success_df_sexpivot1 %>% 
#                                  filter(alive == 'no' & birth_during_monitoring == 'yes'), 
#                                family = poisson(link = "log"), 
#                                REML = TRUE,
#                                pedigree = list(dummy = rcw_ped_object))
# 
# summary_test_lifetime_repro <- summary(test_lifetime_repro)
# coefs_vec <- rownames(summary_test_lifetime_repro$coefficients)
# 
# 
# predict_list <- list()
# 
# for (i in coefs_vec[!coefs_vec %in% c("(Intercept)", 'sexmid')]) {
#   predictor_range <- range(repro_success_df_sexpivot1[repro_success_df_sexpivot1$alive == 'no',i])
#   
#   predict_list[[i]][['range']] <- predictor_range
#   predict_list[[i]][['emmean']] <- as.data.frame(emmeans::emmeans(test_lifetime_repro, 
#                                                     i,
#                                                     at = setNames(list(seq(predictor_range[1], 
#                                                                            predictor_range[2], 
#                                                                            length.out = 100)), 
#                                                                   nm = i),
#                                                     type = "response")
#                                                  )
# }

# predict_list$transloc_anc[['emmean']] %>% 
#   ggplot() +
#   geom_point(data = repro_success_df_sexpivot1 %>%
#                filter(alive == 'no'),
#              aes(x = transloc_anc,
#                  y = life_fldg,
#                  color = first_year_scaled),
#              size = 0.5,
#              alpha = 0.6) +
#   geom_point(data = repro_success_df_sexpivot1 %>%
#                filter(alive == 'no') %>%
#                mutate(transloc_interval = cut(transloc_anc, seq(0, 1, by = 0.1)),
#                       interval_bin = match(transloc_interval, levels(transloc_interval)),
#                       interval_bin = if_else(transloc_anc == 0, 1, interval_bin),
#                       interval_mdpt = c((seq(0, 0.9, 0.1) + seq(0.1, 1, 0.1))/2)[interval_bin]) %>%
#                group_by(interval_mdpt) %>%
#                summarize(life_fldg_25 = quantile(life_fldg, prob = 0.25),
#                          life_fldg_75 = quantile(life_fldg, prob = 0.25),
#                          life_fldg_mean = mean(life_fldg),
#                          count = n()),
#              aes(x = interval_mdpt, y = life_fldg_mean, size = count)) +
#   geom_segment(data = repro_success_df_sexpivot1 %>%
#                  filter(alive == 'no') %>% 
#                  mutate(transloc_interval = cut(transloc_anc, seq(0, 1, by = 0.1)),
#                         interval_bin = match(transloc_interval, levels(transloc_interval)),
#                         interval_bin = if_else(transloc_anc == 0, 1, interval_bin),
#                         interval_mdpt = c((seq(0, 0.9, 0.1) + seq(0.1, 1, 0.1))/2)[interval_bin]) %>% 
#                  group_by(interval_mdpt) %>% 
#                  summarize(life_fldg_25 = quantile(life_fldg, prob = 0.25),
#                            life_fldg_75 = quantile(life_fldg, prob = 0.75),
#                            life_fldg_mean = mean(life_fldg),
#                            count = n()),
#                aes(x = interval_mdpt, 
#                    y = life_fldg_25, 
#                    xend = interval_mdpt, 
#                    yend = life_fldg_75
#   )) +
#   geom_ribbon(aes(transloc_anc,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = transloc_anc, y = rate),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw()
# 
# 
# repro_success_df_sexpivot1 %>%
#   filter(alive == 'no') %>% 
#   mutate(rounded_anc = round(transloc_anc, 2)) %>% 
#   pull(rounded_anc)
# 
# 
# test_vec <- runif(100, 0, 1)
# transloc_interval <- cut(test_vec, seq(0, 1, by = 0.1))
# 
# test_df <- data.frame(transloc_anc = test_vec,
#            interval = transloc_interval,
#            interval_bin = match(transloc_interval, levels(transloc_interval)),
#            interval_mdpt = c((seq(0, 0.9, 0.1) + seq(0.1, 1, 0.1))/2)[test_df$interval_bin])
# predict_list$first_year_scaled[['emmean']] %>% 
#   ggplot() +
#   geom_point(data = repro_success_df_sexpivot1 %>%
#                filter(alive == 'no'),
#              aes(x = first_year_scaled,
#                  y = life_fldg)) +
#   geom_ribbon(aes(first_year_scaled,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = first_year_scaled, y = rate),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw()
# predict_list$mean_group_size[['emmean']] %>% 
#   ggplot() +
#   geom_point(data = repro_success_df_sexpivot1 %>%
#                filter(alive == 'no'),
#              aes(x = mean_group_size,
#                  y = life_fldg)) +
#   geom_ribbon(aes(mean_group_size,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = mean_group_size, y = rate),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw()
# testemmeans %>% 
#   as.data.frame() %>% 
#   ggplot() +
#   geom_point(data = repro_success_df_sexpivot1 %>% 
#                filter(alive == 'no'),
#              aes(x = mean_group_size,
#                  y = life_fldg)) +
#   geom_ribbon(aes(mean_group_size,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = mean_group_size, y = rate ),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw()
# testemmeans <- emmeans::emmeans(test_lifetime_repro, 
#                                 'transloc_anc',
#                                 at = list(transloc_anc = seq(from = 0, to = 1, length.out = 100)),
#                                 type = "response")
# 
# testemmeans %>% 
#   as.data.frame() %>% 
#   ggplot() +
#   geom_point(data = repro_success_df_sexpivot1 %>%
#               filter(alive == 'no'),
#             aes(x = transloc_anc,
#                 y = life_fldg),
#             size = 1, alpha = 0.5) +
#   geom_ribbon(aes(transloc_anc,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = transloc_anc, y = rate ),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw()
# testemmeans <- emmeans::emmeans(test_lifetime_repro, 
#                                 'first_year_scaled',
#                                 at = list(first_year_scaled = seq(from = 0, to = 27, length.out = 100)),
#                                 type = "response")
# testemmeans %>% 
#   as.data.frame() %>% 
#   ggplot() +
#   geom_point(data = repro_success_df_sexpivot1 %>%
#                filter(alive == 'no'),
#              aes(x = first_year_scaled,
#                  y = life_fldg),
#              size = 1, alpha = 0.5) +
#   geom_ribbon(aes(first_year_scaled,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = first_year_scaled, y = rate ),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw()

#pedigreemm::pedigree(sire = rcw_processed_list$ped_processed$fid, 
#                     dam = rcw_processed_list$ped_processed$mid, 
#                     label = rcw_processed_list$ped_processed$id)
#example("pedigree-class")
#
#test_ped <- pedigree(sire = c(NA,NA,1, 1,4,5),
#                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
#
#ped_mod_process <- rcw_processed_list$ped_processed %>% 
#  mutate(fid_na = if_else(fid == '0', NA, fid),
#         mid_na = if_else(fid == '0', NA, mid))
#
#new("pedigree",
#    sire = ped_mod_process$fid_na,
#    dam  = ped_mod_process$mid_na,
#    label = ped_mod_process$id)

#ped_reorder <- reorder_ped(rcw_processed_list$ped_processed,
#            1, 2, 3)
# age_firstyear_lifefledge_mulitpan <- cowplot::plot_grid(
#   repro_success_df_sexpivot1 %>% 
#     filter(alive == 'no' & birth_during_monitoring == 'yes') %>% 
#     ggplot() +
#     geom_point(aes(x = first_year, y = age)) +
#     theme_bw() +
#     ggtitle('cor = -0.13; p = 0.03675'),
#   repro_success_df_sexpivot1 %>% 
#     filter(alive == 'no' & birth_during_monitoring == 'yes') %>% 
#     ggplot() +
#     geom_point(aes(x = age, y = life_fldg)) +
#     theme_bw() +
#     ggtitle('cor = 0.71; p < 2.2e-16'),
#   ncol = 2
# )
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'age_firstyear_lifefledge_mulitpan.png' ),
#                  plot = age_firstyear_lifefledge_mulitpan,
#                  width = 14*0.7, height = 6*0.7, bg = 'white')

# processed_breed_info <- pop_dat %>% 
#   left_join(.,
#             rcw_processed_list$rcws %>% 
#               select(RCWid, MinAge, Sex),
#             by = 'RCWid') %>% 
#   group_by(RCWid) %>%
#   summarize(max_year = max(year),
#             min_year = min(year),
#             MinAge = min(MinAge),
#             Sex = first(Sex),
#             .groups = 'drop') %>% 
#   mutate(max_age = max_year - min_year) %>% 
#   mutate(alive = if_else(RCWid %in% unique(pop_dat[pop_dat$year == 2022,]$RCWid), 
#                          'yes', 'no'),
#          birth_during_monitoring = if_else(RCWid %in% rcw_processed_list$rcws[rcw_processed_list$rcws$MinAge >= 1994,]$RCWid,
#                                            'yes', 'no'),
#          breeder = if_else(RCWid %in% c(rcw_processed_list$nests_processed$fid_dummy, rcw_processed_list$nests_processed$mid_dummy),
#                            1, 0)
#   ) %>% 
#   left_join(.,
#             results_list$rcws_ancestry_info %>% 
#               select(id, group, anc_prop)  %>%
#               pivot_wider(names_from = 'group',
#                           values_from = 'anc_prop',
#                           values_fill = 0) %>% 
#               mutate(transloc_anc = 1 - `non-transloc`) %>%
#               rename(RCWid = id) %>%
#               select(RCWid, transloc_anc),
#             by = 'RCWid')  %>% 
#   left_join(.,
#             results_list$rcws_inbr %>% 
#               rename(RCWid = id) %>% 
#               select(RCWid, f_ped),
#             by = 'RCWid') %>% 
#   mutate(min_year_z = (min_year - mean(min_year))/sd(min_year),
#          max_age_z = (max_age - mean(max_age))/sd(max_age),
#          transloc_anc_z = (transloc_anc - mean(transloc_anc))/sd(transloc_anc),
#          f_ped_z = (f_ped - mean(f_ped))/sd(f_ped))
# 
# processed_breed_info1 <- processed_breed_info %>% 
#   #group_by(RCWid) %>% 
#   filter(Sex %in% c('M', 'F')) %>% 
#   filter(max_year > MinAge) #%>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes')
# pop_dat %>% 
#   left_join(.,
#             rcw_processed_list$rcws %>% 
#               select(RCWid, MinAge, Sex),
#             by = 'RCWid') %>% 
#   group_by(RCWid) %>%
#   summarize(max_year = max(year),
#             min_year = min(year),
#             MinAge = min(MinAge),
#             Sex = first(Sex),
#             .groups = 'drop') %>% 
#   mutate(max_age = max_year - min_year) %>% 
#   mutate(alive = if_else(RCWid %in% unique(pop_dat[pop_dat$year == 2022,]$RCWid), 
#                          'yes', 'no'),
#          birth_during_monitoring = if_else(RCWid %in% rcw_processed_list$rcws[rcw_processed_list$rcws$MinAge >= 1994,]$RCWid,
#                                            'yes', 'no'),
#          breeder = if_else(RCWid %in% c(rcw_processed_list$nests_processed$fid_dummy, rcw_processed_list$nests_processed$mid_dummy),
#                            1, 0)
#   ) %>% 
#   group_by(RCWid) %>%
#   summarize(count = n()) %>% 
#   pull(count)

# variable_vec <- c("max_age_z",
#                   "min_year_z",
#                   "Sex",
#                   "transloc_anc_z",
#                   "f_ped_z")
# 
# variable_combo_list <- unlist(lapply(1:length(variable_vec), combn, 
#                                      x = variable_vec, simplify = FALSE), 
#                               recursive = FALSE)
# 
# #variable_combo_list[[length(variable_combo_list) + 1]] <- "1"
# 
# variable_combo_list1 <- lapply(variable_combo_list, function(x) paste0('breeder ~ ', paste(x, collapse = ' + '), ' + (1|RCWid)'))

# mod_list_breederbinom <- lapply(variable_combo_list1, function(x, dat, ped) {
#   message(x)
#   return(pedigreemm::pedigreemm(as.formula(x), 
#                                 data = dat, 
#                                 family = binomial(link = "logit"), 
#                                 REML = TRUE,
#                                 pedigree = list(RCWid = ped),
#                                 control = glmerControl(optimizer= "Nelder_Mead", #Nelder_Mead",
#                                                        boundary.tol = 1e-15,
#                                                        optCtrl = list(maxfun = 1000000))
#   ))
#   
# }, dat = processed_breed_info1,
# ped = rcw_ped_object)
# mod_list_r2 <- lapply(variable_combo_list1, function(x, dat, ped) {
#   
#   mod <- pedigreemm::pedigreemm(as.formula(x), 
#                                 data = dat, 
#                                 family = binomial(link = "logit"),
#                                 REML = TRUE,
#                                 pedigree = list(RCWid = ped),
#                                 control = glmerControl(optimizer= "Nelder_Mead", #Nelder_Mead",
#                                                        boundary.tol = 1e-15,
#                                                        optCtrl = list(maxfun = 1000000))
#   )
#   MuMIn::r.squaredGLMM(mod)
#   
#   
# }, dat = processed_breed_info1, ped = rcw_ped_object)
# mod_summary_df <- data.frame(mod_id = 1:length(mod_list_breederbinom),
#                              mod_formula = unlist(variable_combo_list1),
#                              AICc = sapply(mod_list_breederbinom, MuMIn::AICc),
#                              R2m = sapply(mod_list_r2, function(x) x['delta', 'R2m']),
#                              R2c = sapply(mod_list_r2, function(x) x['delta', 'R2c'])
#                              )
# # summary_test_topmod_binom <- summary(mod_list_breederbinom[[22]])
# coefs_vec_binom <- rownames(summary_test_topmod_binom$coefficients)
# 
# predict_list_binom <- list()
# 
# for (i in coefs_vec_binom[!coefs_vec_binom %in% c("(Intercept)", 'SexM')]) {
#   predictor_range <- range(processed_breed_info1[,i])
#   
#   predict_list_binom[[i]][['range']] <- predictor_range
#   predict_list_binom[[i]][['emmean']] <- as.data.frame(emmeans::emmeans(mod_list_breederbinom[[22]], 
#                                                                   i,
#                                                                   at = setNames(list(seq(predictor_range[1], 
#                                                                                          predictor_range[2], 
#                                                                                          length.out = 100)), 
#                                                                                 nm = i),
#                                                                   type = "response")
#   )
# }


# binom_breeder_multipan <- cowplot::plot_grid(
#   predict_list_binom$transloc_anc_z$emmean %>% 
#     ggplot() +
#     # geom_point(data = processed_breed_info1,
#     #            aes(x = transloc_anc_z,
#     #                y = breeder),
#     #            size = 1,
#     #            alpha = 0.6) +
#     geom_jitter(data = processed_breed_info1,
#                 aes(x = transloc_anc_z,
#                     y = breeder),
#                 position = position_jitter(height = 0.02),
#                 size = 1,
#                 alpha = 0.6) +
#     geom_ribbon(aes(transloc_anc_z,
#                     ymin = asymp.LCL, ymax = asymp.UCL),
#                 fill = 'gray', alpha = 0.4) +
#     geom_line(aes(x = transloc_anc_z, y = prob),
#               linewidth = 2,
#               color = 'gray') +
#     theme_bw() +
#     ylab('P(breeder)'),
#   predict_list_binom$min_year_z$emmean %>% 
#     ggplot() +
#     #geom_point(data = processed_breed_info1,
#     #           aes(x = min_year_z,
#     #               y = breeder),
#     #           size = 1,
#     #           alpha = 0.6) +
#     geom_jitter(data = processed_breed_info1,
#                 aes(x = min_year_z,
#                     y = breeder),
#                 position = position_jitter(height = 0.02),
#                 size = 1,
#                 alpha = 0.6) +
#     geom_ribbon(aes(min_year_z,
#                     ymin = asymp.LCL, ymax = asymp.UCL),
#                 fill = 'gray', alpha = 0.4) +
#     geom_line(aes(x = min_year_z, y = prob),
#               linewidth = 2,
#               color = 'gray') +
#     theme_bw() +
#     ylab('P(breeder)'),
#   ncol = 1
# )
# 
# cowplot::ggsave2(filename = here('figures', 'test_figs', 'figures', 'binom_breeder_multipan.png' ),
#                  plot = binom_breeder_multipan,
#                  width = 8*0.8, height = 12*0.8, bg = 'white')
# 
#library("lme4qtl")
# testemmeans <- emmeans::emmeans(mod_list[[22]], 
#                                 'transloc_anc_z',
#                                 at = list(transloc_anc_z = seq(from = range(processed_breed_info1$transloc_anc_z)[1], 
#                                                                to = range(processed_breed_info1$transloc_anc_z)[2], 
#                                                                length.out = 100)),
#                                 type = "response")
# 
# # as.data.frame(testemmeans) %>% 
#   ggplot() +
#   geom_point(data = processed_breed_info1,
#                            aes(x = transloc_anc_z,
#                                y = breeder),
#                            size = 1,
#                            alpha = 0.6) +
#   geom_ribbon(aes(transloc_anc_z,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = transloc_anc_z, y = prob),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw() +
#   ylab('P(breeder)')
# 
# testemmeans <- emmeans::emmeans(mod_list_breederbinom[[22]], 
#                                 'min_year_z',
#                                 at = list(min_year_z = seq(from = range(processed_breed_info1$min_year_z)[1], 
#                                                                to = range(processed_breed_info1$min_year_z)[2], 
#                                                                length.out = 100)),
#                                 type = "response")
# 
# as.data.frame(testemmeans) %>% 
#   ggplot() +
#   geom_point(data = processed_breed_info1,
#              aes(x = min_year_z,
#                  y = breeder),
#              size = 1,
#              alpha = 0.6) +
#   geom_ribbon(aes(min_year_z,
#                   ymin = asymp.LCL, ymax = asymp.UCL),
#               fill = 'gray', alpha = 0.4) +
#   geom_line(aes(x = min_year_z, y = prob),
#             linewidth = 2,
#             color = 'gray') +
#   theme_bw() +
#   ylab('P(breeder)')

# 
# binom_pedigreemm <- pedigreemm::pedigreemm(breeder ~ max_age_z + min_year_z + Sex + anc_prop_z + f_ped_z + (1|RCWid), 
#                                           data = processed_breed_info1,
#                                           family = binomial(link = "logit"),
#                                           REML = TRUE,
#                                           pedigree = list(RCWid = rcw_ped_object),
#                                           control = glmerControl(optimizer= "Nelder_Mead", #Nelder_Mead",
#                                                                  boundary.tol = 1e-15,
#                                                                  optCtrl = list(maxfun = 1000000)))

# data.frame(id = c(rcw_processed_list$nests_processed$fid_dummy, rcw_processed_list$nests_processed$mid_dummy),
#            sex1 = c(rep('M', length(rcw_processed_list$nests_processed$fid_dummy)), rep('F', length(rcw_processed_list$nests_processed$mid_dummy)))) %>% 
#   group_by(id) %>% 
#   slice_head(n = 1) %>% 
#   ungroup()


# nests_with_offspring_in_pop <- rcw_processed_list$rcws %>% 
#   filter(RCWid %in% unique(pop_dat$RCWid)) %>% 
#   filter(!is.na(NatalNest)) %>% 
#   pull(NatalNest) %>% 
#   unique()
# 
# 
# 
# indiv_by_birth_year <- split(rcw_processed_list$rcws$RCWid, rcw_processed_list$rcws$MinAge)
# 
# 
# lapply(indiv_by_birth_year, function(x, last_year_info) {
#   sum(x %in% last_year_info)
# }, last_year_info = pop_dat$RCWid[pop_dat$year == 2022])
# 
# 
# library(brms)
# 

# cor.test(
#   repro_success_df_sexpivot1[repro_success_df_sexpivot1$full_grandparent_info == 'yes',]$fped,
#   repro_success_df_sexpivot1[repro_success_df_sexpivot1$full_grandparent_info == 'yes',]$transloc_anc,
#   method = 'spearman'
#   )
#          
# 
# repro_success_df_sexpivot1 %>% 
#   filter(full_grandparent_info == 'yes') %>% 
#   filter(fped > 0.12)
# 
# 
# repro_success_df_sexpivot1 %>% 
#   filter(full_grandparent_info == 'yes') %>% 
#   ggplot() +
#   geom_point(aes(y = fped,
#                  x = transloc_anc))
# 


# repro_success_df_sexpivot_process <- repro_success_df_sexpivot %>% 
#   # rowwise() %>% 
#   # filter(#full_grandparent_info == 'yes' &
#   #          !is.na(dummy) &
#   #          !is.na(fldg_num) &
#   #          #(fldg_num <= hatch_num) &
#   #          identical(fldg_num <= hatch_num, TRUE) &
#   #           !is.na(group_size)) %>% 
#   rowwise() %>% 
#   filter(identical(fldg_num <= hatch_num, TRUE)) %>% 
#   filter(!dummy %in% repro_success_df_sexpivot_missinfo) %>% 
#   ungroup() %>% 
#   filter(!dummy %in% problem_ids) %>% 
#   # group_by(dummy) %>% 
#   # summarize(#transloc = first(transloc),
#   #   full_grandparent_info = first(full_grandparent_info),
#   #   sex = first(sex),
#   #   mean_group_size = mean(group_size),
#   #   transloc_anc = first(transloc_anc),
#   #   fped = first(fped),
#   #   first_year = min(year),
#   #   life_fldg = sum(fldg_num),
#   #   .groups = 'drop') %>% 
#   mutate(alive = if_else(dummy %in% unique(pop_dat[pop_dat$year == 2022,]$'RCWid'), 
#                          'yes', 'no'),
#          birth_during_monitoring = if_else(dummy %in% rcw_processed_list$rcws[rcw_processed_list$rcws$MinAge >= 1994,]$RCWid,
#                                            'yes', 'no')#,
#          #first_year_scaled = first_year - min(first_year)
#          ) %>% 
#   left_join(.,
#             pop_dat %>% 
#               group_by(RCWid) %>% 
#               summarize(age = max(year) - min(year)) %>% 
#               rename(dummy = RCWid),
#             by = 'dummy') %>% 
#   left_join(., processed_anc_info,
#             by = 'dummy') %>%
#   group_by(dummy) %>% 
#   mutate(first_year = min(year),
#          scaled_year = year - first_year) %>% 
#   ungroup()
#   
# %>% 
#  mutate(mean_group_size_z = (mean_group_size - mean(mean_group_size))/sd(mean_group_size),
#         transloc_anc_z = (transloc_anc - mean(transloc_anc))/sd(transloc_anc),
#         fped_z = (fped - mean(fped))/sd(fped),
#         first_year_scaled_z = (first_year_scaled - mean(first_year_scaled))/sd(first_year_scaled),
#         compl_simpson_z = (compl_simpson - mean(compl_simpson))/sd(compl_simpson),
#         anc_count_z = (anc_count - mean(anc_count))/sd(anc_count),
#  )


# summary(
#   pedigreemm::pedigreemm(fldg_num ~ aggregation +
#                            first_year +
#                            scaled_year +
#                            (1|dummy), 
#                          data = repro_success_df_sexpivot_process %>% 
#                            filter(birth_during_monitoring == 'yes'), 
#                          family = poisson(link = "log"), 
#                          REML = TRUE,
#                          pedigree = list(dummy = rcw_ped_object)
#   )
#   
# )



# df <- data.frame("data" = runif(1000))
# 
# dens <- density(df$data)
# new_df1 <- data.frame(y = c(dens$x[dens$x < 0.5], rev(dens$x[dens$x < 0.5])),
#                       x = c(-dens$y[dens$x < 0.5], rev(dens$y[dens$x < 0.5])),
#                       z = 'red2')
# new_df2 <- data.frame(y = c(dens$x[dens$x >= 0.5], rev(dens$x[dens$x >= 0.5])),
#                       x = c(-dens$y[dens$x >= 0.5], rev(dens$y[dens$x >= 0.5])),
#                       z = 'green3')
# 
# df %>% 
#   ggplot() +
#   geom_violin(aes(x, y))
# 
# ggplot(rbind(new_df1, new_df2), aes(x, y, fill = z)) + 
#   geom_polygon() +
#   scale_fill_identity() +
#   scale_x_continuous(breaks = 0, expand = c(1, 1), labels = 'DATA', name = '')


# #Data setup
# set.seed(123)
# dat <- data.frame(x = rep(1:3,each = 100),
#                   y = c(rnorm(100,-1),rnorm(100,0),rnorm(100,1)))


# p <- ggplot() + 
#   geom_violin(data = dat,aes(x = factor(x),y = y))
# p_build <- ggplot2::ggplot_build(p)$data[[1]]
# 
# p_build <- transform(p_build,
#                      xminv = x - violinwidth * (x - xmin),
#                      xmaxv = x + violinwidth * (xmax - x))
# 
# p_build <- rbind(plyr::arrange(transform(p_build, x = xminv), y),
#                  plyr::arrange(transform(p_build, x = xmaxv), -y))
# 
# 
# p_build$fill_group <- ifelse(p_build$y >= 0,'Above','Below')
# #This is necessary to ensure that instead of trying to draw
# # 3 polygons, we're telling ggplot to draw six polygons
# p_build$group1 <- with(p_build,interaction(factor(group),factor(fill_group)))
# 
# ggplot() + 
#   geom_violin(data = dat, aes(x = x, y = y, group = x)) +
#   geom_polygon(data = p_build,
#                aes(x = x,y = y,group = group1,fill = fill_group))

##################
##################
##################
# mod_dat <- repro_success_df_sexpivot1 %>% 
#   #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#   filter(birth_during_monitoring == 'yes')  %>%
#   mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#          birth_year = MinAge - min(MinAge))
# 
# var_vec <- c('intercept', 'mean_group_size', 'transloc_anc', 'first_year_scaled', 'anc_count', 'fped', 'sex')
# response_vec <- c('life_fldg', 'nesting_years', 'mean_fldg_per_year', 'birth_firstbreed_dif', 'birth_lastbreed_dif')
# 
# mod_dat_seq <- mod_dat %>% 
#   reframe(mean_group_size = seq(min(mean_group_size), max(mean_group_size), length.out = 50),
#           transloc_anc = seq(min(transloc_anc), max(transloc_anc), length.out = 50),
#           first_year_scaled = seq(min(first_year_scaled), max(first_year_scaled), length.out = 50),
#           anc_count = seq(min(anc_count), max(anc_count), length.out = 50),
#           fped = seq(min(fped), max(fped), length.out = 50),
#           sex = rep(c(0, 1), each = 25)) %>% 
#   mutate(intercept = 1)
#
# prior_func <- function(response, variable, draws = 100) {
#   if (response == 'life_fldg') {
#     
#     prior_vals <- switch(variable,
#                          intercept = {rnorm(draws, 0, 1.5)},
#                          mean_group_size = {rnorm(draws, 0, 0.42)},
#                          transloc_anc = {rnorm(draws, 0, 2)},
#                          first_year_scaled = {rnorm(draws, 0, 0.08)},
#                          anc_count = {rnorm(draws, 0, 0.42)},
#                          fped = {rnorm(draws, 0, 15)},
#                          sex = {rnorm(draws, 0, 1.9)},
#                          stop('not a valid variable'))
#     
#   } else if (response == 'nesting_years') {
#     
#     prior_vals <- switch(variable,
#                          intercept = {rnorm(draws, 0, 1.3)},
#                          mean_group_size = {rnorm(draws, 0, 0.41)},
#                          transloc_anc = {rnorm(draws, 0, 1.9)},
#                          first_year_scaled = {rnorm(draws, 0, 0.07)},
#                          anc_count = {rnorm(draws, 0, 0.41)},
#                          fped = {rnorm(draws, 0, 14)},
#                          sex = {rnorm(draws, 0, 1.7)},
#                          stop('not a valid variable'))
#     
#   } else if (response == 'mean_fldg_per_year') {
#     
#     prior_vals <- switch(variable,
#                          intercept = {rnorm(draws, 0, 1.4)},
#                          mean_group_size = {rnorm(draws, 0, 0.15)},
#                          transloc_anc = {rnorm(draws, 0, 1.25)},
#                          first_year_scaled = {rnorm(draws, 0, 0.1)},
#                          anc_count = {rnorm(draws, 0, 0.5)},
#                          fped = {rnorm(draws, 0, 1.25)},
#                          sex = {rnorm(draws, 0, 1.25)},
#                          stop('not a valid variable'))
#     
#   } else if (response == 'birth_firstbreed_dif') {
#     
#     prior_vals <- switch(variable,
#                          intercept = {rnorm(draws, 0, 1)},
#                          mean_group_size = {rnorm(draws, 0, 0.32)},
#                          transloc_anc = {rnorm(draws, 0, 1.5)},
#                          first_year_scaled = {rnorm(draws, 0, 0.05)},
#                          anc_count = {rnorm(draws, 0, 0.38)},
#                          fped = {rnorm(draws, 0, 13)},
#                          sex = {rnorm(draws, 0, 1.5)},
#                          stop('not a valid variable'))
#     
#   } else if (response == 'birth_lastbreed_dif') {
#     
#     prior_vals <- switch(variable,
#                          intercept = {rnorm(draws, 0, 1.3)},
#                          mean_group_size = {rnorm(draws, 0, 0.41)},
#                          transloc_anc = {rnorm(draws, 0, 1.9)},
#                          first_year_scaled = {rnorm(draws, 0, 0.07)},
#                          anc_count = {rnorm(draws, 0, 0.41)},
#                          fped = {rnorm(draws, 0, 14)},
#                          sex = {rnorm(draws, 0, 1.7)},
#                          stop('not a valid variable'))
#     
#   } else {
#     stop('Not a valid response')
#   }
#   
#   return(prior_vals)
# }
#
# param_draws <- 500
# sim_per_param_draw <- 50
# 
# prior_draws_list <- list()
# 
# for (RESPONSE in response_vec) {
#   var_list <- list()
#   for (VAR in var_vec) {
#     
#     prior_vec <- prior_func(response = RESPONSE, 
#                             variable = VAR, 
#                             draws = param_draws)
#     
#     var_list[[VAR]] <- lapply(prior_vec, function(param_val, dat, var, draw_count) {
#       data.frame(var = var,
#                  param_val = param_val,
#                  y = rpois(50, lambda = exp(param_val*dat[[var]])))
#     }, dat = mod_dat_seq, var = VAR, draw_count = sim_per_param_draw) %>% 
#       bind_rows() %>% 
#       mutate(type = 'simulated')
#     
#     var_list[[VAR]] <- var_list[[VAR]] %>% 
#       filter(y <= quantile(y, prob = 0.995) ) %>% 
#       rbind(., mod_dat %>% 
#               select(all_of(RESPONSE)) %>% 
#               rename(y = all_of(RESPONSE)) %>% 
#               mutate(type = 'obs',
#                      var = VAR,
#                      param_val = NA) %>% 
#               select(var, param_val, y, type)
#             )
#   }
#   
#   prior_draws_list[[RESPONSE]] <- bind_rows(var_list)
#   
# }
#
# prior_draws_list$life_fldg %>% 
#   #bind_rows(.id = 'response') %>% 
#   ggplot() +
#   geom_jitter(aes(x = 'a', y = y, color = type)) +
#   theme_bw() +
#   facet_wrap(~var, scales = 'free_y', nrow = 1)
# 
# prior_draws_list$nesting_years %>% 
#   #bind_rows(.id = 'response') %>% 
#   ggplot() +
#   geom_jitter(aes(x = 'a', y = y, color = type)) +
#   theme_bw() +
#   facet_wrap(~var, scales = 'free_y', nrow = 1)
# 
# prior_draws_list$birth_firstbreed_dif %>% 
#   #bind_rows(.id = 'response') %>% 
#   ggplot() +
#   geom_jitter(aes(x = 'a', y = y, color = type), height = 0.1) +
#   theme_bw() +
#   facet_wrap(~var, scales = 'free_y', nrow = 1)
# 
# prior_draws_list$birth_lastbreed_dif %>% 
#   #bind_rows(.id = 'response') %>% 
#   ggplot() +
#   geom_jitter(aes(x = 'a', y = y, color = type), height = 0.1) +
#   theme_bw() +
#   facet_wrap(~var, scales = 'free_y', nrow = 1)
#
# prior_draws_list$life_fldg %>% 
#   #bind_rows(.id = 'response') %>% 
#   ggplot() +
#   geom_jitter(aes(x = 'a', y = y, color = type)) +
#   theme_bw() +
#   facet_grid(var ~ response, scales = 'free_y')
# 
# prior_sim_df %>% 
#   select(y) %>% 
#   mutate(type = 'sim') %>% 
#   rbind(mod_dat %>% 
#           select(life_fldg) %>% 
#           rename(y = life_fldg) %>% 
#           mutate(type = 'obs')
#   ) %>% 
#   #filter(param == '-0.2048') %>% 
#   ggplot() +
#   geom_jitter(aes(x = 'a', y = y, color = type)) +
#   theme_bw()
#
# prior_intercept <- rnorm(500, 0, 1.25)
# prior_mean_group_size <- rnorm(500, 0, 0.35)
# prior_sim_df <- lapply(prior_mean_group_size, function(param_val, dat) {
#   data.frame(param = as.character(round(param_val, 4)),
#              x = mod_dat_seq$mean_group_size,
#              y = rpois(50, lambda = exp(param_val*mod_dat_seq$mean_group_size)))
# }, dat = mod_dat_seq) %>% 
#   bind_rows()
# 
# prior_sim_df %>% 
#   select(y) %>% 
#   mutate(type = 'sim') %>% 
#   rbind(mod_dat %>% 
#           select(life_fldg) %>% 
#           rename(y = life_fldg) %>% 
#         mutate(type = 'obs')
#         ) %>% 
#   #filter(param == '-0.2048') %>% 
#   ggplot() +
#   geom_jitter(aes(x = 'a', y = y, color = type)) +
#   theme_bw()
#
# prior_mean_group_size <- rnorm(500, 0, 1.4)
# 
# prior_sim_df <- lapply(prior_mean_group_size, function(param_val, dat) {
#   data.frame(param = as.character(round(param_val, 4)),
#              x = mod_dat_seq$mean_group_size,
#              y = rpois(50, lambda = exp(param_val*mod_dat_seq$transloc_anc)))
# }, dat = mod_dat_seq) %>% 
#   bind_rows()
#
# prior_sim_df %>% 
#   select(y) %>% 
#   mutate(type = 'sim') %>% 
#   rbind(mod_dat %>% 
#           select(life_fldg) %>% 
#           rename(y = life_fldg) %>% 
#           mutate(type = 'obs')
#   ) %>% 
#   #filter(param == '-0.2048') %>% 
#   ggplot() +
#   geom_jitter(aes(x = 'a', y = y, color = type)) +
#   theme_bw()
# prior_mean_group_size <- rnorm(500, 0, .08)
# 
# prior_sim_df <- lapply(prior_mean_group_size, function(param_val, dat) {
#   data.frame(param = as.character(round(param_val, 4)),
#              x = mod_dat_seq$mean_group_size,
#              y = rpois(50, lambda = exp(param_val*mod_dat_seq$first_year_scaled)))
# }, dat = mod_dat_seq) %>% 
#   bind_rows()
# 
# prior_sim_df %>% 
#   select(y) %>% 
#   mutate(type = 'sim') %>% 
#   rbind(mod_dat %>% 
#           select(life_fldg) %>% 
#           rename(y = life_fldg) %>% 
#           mutate(type = 'obs')
#   ) %>% 
#   #filter(param == '-0.2048') %>% 
#   ggplot() +
#   geom_jitter(aes(x = 'a', y = y, color = type)) +
#   theme_bw()


# firstbreed_mod1 <- brm(
#   birth_firstbreed_dif | cens(censored) ~ mean_group_size + first_year_scaled + transloc_anc + anc_count_z + (1 | gr(dummy, cov = A)),
#   data = repro_success_df_sexpivot1 %>%
#     #filter(alive == 'no' & birth_during_monitoring == 'yes') %>%
#     filter(birth_during_monitoring == 'yes') %>%
#     mutate(censored = if_else(alive == 'yes', 'right', 'none'),
#            birth_year = MinAge - min(MinAge)), 
#   data2=list(A = rcw_addrel),
#   family = poisson(link = 'log'), 
#   #cov_ranef = list(A = rcw_addrel),
#   chains = 3, cores = 2, iter = 10000, 
#   control = list(adapt_delta = 0.9)
# )

# (flat)         b                                                    default
# (flat)         b       anc_count_z                             (vectorized)
# (flat)         b first_year_scaled                             (vectorized)
# (flat)         b   mean_group_size                             (vectorized)
# (flat)         b      transloc_anc 

# intercept = {rnorm(draws, 0, 1.5)},
# mean_group_size = {rnorm(draws, 0, 0.42)},
# transloc_anc = {rnorm(draws, 0, 2)},
# first_year_scaled = {rnorm(draws, 0, 0.08)},
# anc_count = {rnorm(draws, 0, 0.42)},
# fped = {rnorm(draws, 0, 15)},
# sex = {rnorm(draws, 0, 1.9)}
# 
# bprior <- c(prior_string("normal(0,0.42)", class = "b", coef = "anc_count"),
#             prior_string("normal(0,0.08)", class = "b", coef = "first_year_scaled"),
#             prior_string("normal(0,0.42)", class = "b", coef = "mean_group_size"),
#             prior_string("normal(0,2)", class = "b", coef = "transloc_anc"),
#             prior_string("normal(0, 1.5)", class = "Intercept"),
#             #prior_string("normal(0,0.001)", class = "sd", coef = 'dummy'),
#             prior_string("normal(0,0.02)", class = "sd", coef = 'Intercept', group = 'dummy')
#             )
