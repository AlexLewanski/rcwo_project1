########################################################
########################################################
### RCW DATA EXPLORATION, EVALUATION, AND PROCESSING ###
########################################################
########################################################

#####################
### SCRIPT SET-UP ###
#####################

### PACKAGES ###
library(tidyverse)
library(readxl)
library(clock)
library(here)

here::here()


### LOADING CUSTOM FUNCTIONS AND DATA ###

#custom functions
devtools::load_all('/Users/alexlewanski/Documents/r_packages/pedutils')
source(here('scripts', 'data_processing_custom_functions.R'))


#data
file_names <- c('Nests.xlsx', 'RCWs.xlsx', 'RCWs_wout_format.xlsx', 'LastSeenDetail.xlsx', 'Translocation.xlsx')
raw_dat_list <- lapply(setNames(file_names, tolower(gsub('\\.xlsx', '', file_names))), function(x) {
  read_excel(here('data', 'feb2024_databasemarch2023', x))
})

#rcws_nest_query1 <- read_excel(here('data', 'feb2024_databasemarch2023', "rcws_nest_query1.xlsx"))




##############################
### CENSUS DATA PROCESSING ###
##############################

#filter census data down to the post-breeding census where birds were actually observed
summer_census_detected <- raw_dat_list$lastseendetail %>% 
  filter(Detected == 'Yes' & SurveyType == "Census")  %>% 
  filter(get_month(RecordDate) >= 6 & get_month(RecordDate) <= 8) %>% 
  mutate(year = get_year(RecordDate))  %>% 
  mutate(RCWid_old = RCWid) %>% 
  mutate(RCWid = case_when(year < 2005 & RCWid == 'GB=Z' ~ 'GB-Z',
                           TRUE ~ RCWid)) %>% 
  select(RCWid, Detected, SurveyType, year) %>% 
  group_by(RCWid, year) %>% 
  slice_head(n = 1) #if there are more than one observations/indiv/census, keep only one observation


#filter rcws info down to birth information
#add empty Detected column so that it can be combined with the census data
birth_year_info <- raw_dat_list$rcws %>% 
  mutate(Detected = NA,
         SurveyType = 'birth_info') %>% 
  rename(year = MinAge) %>% 
  select(RCWid, Detected, SurveyType, year) %>% 
  filter(!is.na(year))

#combine birth and census information
obs_info_combine <- bind_rows(summer_census_detected,
                              birth_year_info)

#add observations for any year that a bird was not observed when the observation
#was bookended by observation and/or birth
obs_info_combine_addint <- add_intervening_years(data = obs_info_combine,
                                                 id_col = 'RCWid',
                                                 year_col = 'year')

#combine the census/birth info with the information on each woodpecker
obs_info_combine_addint1 <- obs_info_combine_addint %>%   
  left_join(., raw_dat_list$rcws %>% 
              select(RCWid, Origin),
            by = 'RCWid') %>% 
  mutate(Origin_update = case_when(is.na(Origin) ~ 'unknown',
                                   TRUE ~ Origin))

#unknown_indivs <- unique(obs_info_combine_addint1[obs_info_combine_addint1$Origin_update == 'unknown',]$RCWid)


pop_filters_df <- lapply(split(obs_info_combine_addint1, obs_info_combine_addint1$RCWid), function(DF) {
  
  #1: CENSUS_ONLY
  #Only census observations
  DF$Scenario1 <- DF$SurveyType %in% 'Census'
  
  #2: CENSUS + DEDUCED
  #Census observations and filling in the the gaps between non-consecutive 
  #census observations
  if ('Census' %in% DF$SurveyType) {
    first_year <- min(DF$year[DF$SurveyType %in% 'Census'])
    DF$Scenario2 <- DF$year >= first_year
  } else {
    DF$Scenario2 <- FALSE
  }
  
  #3: CENSUS/BIRTH + DEDUCED (INCLUDE BIRTH YEAR IF FLEDGED, BORN IN THE POP, AND FOUND IN A SUBSEQUENT CENSUS)
  #Similar to 2. Include all census observations but for birds born in the population, record 
  #the start year as the year where the 
  if (DF$Origin_update[1] %in% c('Immigrant', 'unknown')) {
    
    if ('Census' %in% DF$SurveyType) {
      first_year <- min(DF$year[DF$SurveyType %in% 'Census'])
      DF$Scenario3 <- DF$year >= first_year  & DF$SurveyType %in% c('Census', NA)
    } else {
      DF$Scenario3 <- FALSE
    }
    
  } else {
    
    if (all(DF$SurveyType %in% 'birth_info') ) {
      
      DF$Scenario3 <- FALSE
    } else {
      first_year <- min(DF$year)
      DF$Scenario3 <- switch(DF$Origin_update[1],
                                      'Local' = DF$year >= first_year,
                                      'Translocated' = DF$year > first_year)
      
      if (sum(DF$year == min(DF$year)) > 1) {
        DF$Scenario3[which(DF$year %in% min(DF$year))[2]] <- FALSE
      }
    }
  }
  
  #4: CENSUS/BIRTH + DEDUCED (EXCLUDE BIRTH YEAR EVEN IF FLEDGED AND FOUND IN A SUBSEQUENT CENSUS)
  #Identical to 3 but exclude birth year. This would represent a focus on the breeder population
  if (DF$Origin_update[1] %in% c('Immigrant', 'unknown')) {
    
    if ('Census' %in% DF$SurveyType) {
      first_year <- min(DF$year[DF$SurveyType %in% 'Census'])
      DF$Scenario4 <- DF$year >= first_year & DF$SurveyType %in% c('Census', NA)
    } else {
      DF$Scenario4 <- FALSE
    }
    
  } else {
    
    if (all(DF$SurveyType %in% 'birth_info') ) {
      
      DF$Scenario4 <- FALSE
    } else {
      first_year <- min(DF$year)
      DF$Scenario4 <- DF$year > first_year
    }
  } 
  
  return(DF)
  
}) %>% 
  bind_rows() %>% 
  group_by(year, RCWid) %>% 
  mutate(Scenario2 = case_when(sum(Scenario2) > 1 ~ c(TRUE, rep(FALSE, length(Scenario2) - 1)),
                         TRUE ~ Scenario2))

scenarios_long_list <- lapply(setNames(nm = paste0('Scenario', 1:4)), function(x, pop_info) {
  pop_info[pop_info %>% pull(x),]
}, pop_info = pop_filters_df)


count_ind_year_scenario <- scenarios_long_list %>% 
  bind_rows(.id = 'scenario') %>% 
  group_by(scenario, year, RCWid) %>% 
  summarize(count = n(), .groups = 'drop') %>% 
  pull(count)


message("Are all individuals found at most once a year in the processed census info?: ", all(count_ind_year_scenario == 1))


year_seq_ind_year_scenario <- scenarios_long_list %>% 
  bind_rows(.id = 'scenario') %>% 
  group_by(scenario, RCWid) %>% 
  summarize(year_check = identical(as.double(sort(year)), as.double(min(year):max(year))),
            .groups = 'drop') %>% 
  filter(!year_check & !(scenario %in% 'Scenario1')) %>% 
  pull(year_check) %>% 
  all()

message("Are all individuals found over a continuous sequence of years in the processed census info?: ", year_seq_ind_year_scenario)



###########################
### PEDIGREE PROCESSING ###
###########################

rcw_pedigree_info <- raw_dat_list$nests %>% 
  select(ID, MaleID, FemaleID) %>% 
  #left_join(., raw_dat_list$rcws_wout_format[,c('ID', 'RCWid')] %>% 
  #            rename(MaleID = ID), by = 'MaleID') %>% 
  #select(-MaleID) %>% 
  #rename(MaleID = RCWid)  %>% 
  #left_join(., raw_dat_list$rcws_wout_format[,c('ID', 'RCWid')] %>% 
  #            rename(FemaleID = ID), by = 'FemaleID') %>% 
  #select(-FemaleID) %>% 
  rename(#FemaleID = RCWid,
         NatalNest = ID) %>% 
  left_join(raw_dat_list$rcws_wout_format[,c('RCWid', 'NatalNest', 'Sex')], 
            .,
            by = 'NatalNest')

rcws_ped <- process_ped(ped = rcw_pedigree_info,
                        id_col = 'RCWid',
                        sire_col = "MaleID",
                        dam_col = "FemaleID",
                        sex_col= 'Sex',
                        founder_val = NA,
                        sex_vals = list(male = "M", 
                                        female = 'F', 
                                        unknown = c("U", NA)),
                        keep_extra_cols = TRUE,
                        disable_sex_check = FALSE)

rcws_ped_dummy_pars <- ped_add_dummy_parents(ped = rcws_ped, 
                                             id = 'id',
                                             fid = 'fid', 
                                             mid = 'mid',
                                             sex = 'sex',
                                             nest_id = 'NatalNest', 
                                             founder_parent_vals = '0')


#pedigree checks
sire_check <- all(rcws_ped_dummy_pars$ped$fid %in% rcws_ped_dummy_pars$ped$id | rcws_ped_dummy_pars$ped$fid == '0')
dam_check <- all(rcws_ped_dummy_pars$ped$mid %in% rcws_ped_dummy_pars$ped$id | rcws_ped_dummy_pars$ped$mid == '0')

message('Do all sires have an entry (i.e., are found in the column)? ', sire_check)
message('Do all dams have an entry (i.e., are found in the column)? ', dam_check)
message('Are there any repeated entries for individuals in the pedigree file? ', any(duplicated(rcws_ped_dummy_pars$ped$id)))



###########################
### PEDIGREE PROCESSING ###
###########################

### ADD DUMMY PARENTS TO THE NEST INFORMATION
raw_dat_list$nests$fid_dummy <- raw_dat_list$nests$MaleID
raw_dat_list$nests$mid_dummy <- raw_dat_list$nests$FemaleID

dummy_male_info <- rcws_ped_dummy_pars$ped[rcws_ped_dummy_pars$ped$fid %in% rcws_ped_dummy_pars$dummy_ids$id[rcws_ped_dummy_pars$dummy_ids$sex == 1],]
dummy_female_info <- rcws_ped_dummy_pars$ped[rcws_ped_dummy_pars$ped$mid %in% rcws_ped_dummy_pars$dummy_ids$id[rcws_ped_dummy_pars$dummy_ids$sex == 2],]

raw_dat_list$nests$fid_dummy[match(dummy_male_info$NatalNest, raw_dat_list$nests$ID)] <- dummy_male_info$fid
raw_dat_list$nests$mid_dummy[match(dummy_female_info$NatalNest, raw_dat_list$nests$ID)] <- dummy_female_info$mid



#################################
### OUTPUTTING PROCESSED DATA ###
#################################

write.csv(raw_dat_list$translocation,
          here('data', 'feb2024_databasemarch2023_processed', 'translocation.csv'),
          row.names = FALSE)

write.csv(pop_filters_df,
          here('data', 'feb2024_databasemarch2023_processed', 'census_processed.csv'),
          row.names = FALSE)

write.csv(raw_dat_list$nests,
          here('data', 'feb2024_databasemarch2023_processed', 'nests_processed.csv'),
          row.names = FALSE)

write.csv(raw_dat_list$rcws,
          here('data', 'feb2024_databasemarch2023_processed', 'rcws.csv'),
          row.names = FALSE)

write.csv(rcws_ped_dummy_pars$ped,
          here('data', 'feb2024_databasemarch2023_processed', 'ped_processed.csv'),
          row.names = FALSE)

write.csv(rcws_ped_dummy_pars$dummy_ids,
          here('data', 'feb2024_databasemarch2023_processed', 'dummy_parents_info.csv'),
          row.names = FALSE)





as.data.frame(raw_dat_list$translocation[ raw_dat_list$translocation$Type == 'Inter-population',])








scenarios_long_list %>% 
  bind_rows(.id = 'scenario') %>% 
  group_by(scenario, year) %>% 
  summarize(`Census count` = n(),
            .groups = 'drop') %>% 
  ggplot(aes(x = year, y = `Census count`, color = scenario)) +
  geom_line(alpha = 1.75, linewidth = 0.5) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c('#2a9d8f', '#e76f51', '#a64ca6', 'gray')) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  xlab('Year')



































































# observation_type_val <- (table(summer_census_detected_addint$row_origin)/nrow(summer_census_detected_addint))*100
# 
# 
# missing_census_info <- lapply(split(summer_census_detected, summer_census_detected$RCWid), function(x) {
#   dif_vec <- sort(x$year, decreasing = TRUE)[-length(x$year)] - sort(x$year, decreasing = TRUE)[-1]
#   dif_vec1 <- dif_vec[dif_vec > 1] - 1
#   if (length(dif_vec1) == 0) dif_vec1 <- 0
#   
#   return(data.frame(gap_count = dif_vec1))
# }) %>% 
#   bind_rows(.id = 'RCWid')
# 
# 
# census_summarized <- summer_census_detected_addint %>% 
#   group_by(year) %>% 
#   summarize('observed + deduced' = n(),
#             'observed' = sum(row_origin == 'observed'), 
#             .groups = 'drop') %>% 
#   #full_join(., census_minage_addint, by = )
#   pivot_longer(cols = c("observed + deduced", "observed"), values_to = 'Census count') %>% 
#   bind_rows(.,  
#             census_minage_addint_summarize)
# 
# 
# 
# census_summarized
# 
# 
# 
# cor_mat <- pivot_wider(census_summarized, names_from = 'name', values_from = "Census count") %>% 
#   select(observed, `observed + deduced`, `observed + deduced (min age)`) %>% 
#   na.omit() %>% 
#   cor(method = 'pearson')
# 
# 
# 
# 
# census_plot_scenarios <- census_summarized %>% 
#   ggplot(aes(x = year, y = `Census count`, color = name)) +
#   geom_line(alpha = 1.75, size = 0.5) +
#   geom_point(size = 2.5) +
#   scale_color_manual(values = c('#2a9d8f', '#e76f51', '#a64ca6')) +
#   theme_classic() +
#   theme(legend.title = element_blank()) +
#   xlab('Year') +
#   annotate('text', x = 1997, y = 150, 
#            label = paste0('obs vs. obs + deduced: cor = ', round(cor_mat[1,2], 3)),
#            hjust = 0,
#            size = 4) +
#   annotate('text', x = 1997, y = 160, 
#            label = paste0('obs vs. obs + deduced (minage): cor = ', round(cor_mat[1,3], 3)),
#            hjust = 0,
#            size = 4)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# deduced_vs_observed_plot <- summer_census_detected_addint %>% 
#   group_by(year, row_origin) %>% 
#   summarize(count = n(), .groups = 'drop') %>% 
#   group_by(year) %>%
#   mutate(count_perc = (count/sum(count))*100 ) %>% 
#   ungroup() %>% 
#   mutate(row_origin = factor(row_origin, levels = c('observed', 'deduced'))) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = row_origin),
#            #color = 'gray', size = 0.2,
#            position = "stack", stat="identity", width = 0.98, alpha = 0.7) +
#   scale_fill_manual(name = 'Observation type', 
#                     values = c('#2a9d8f', '#e76f51')) + #c('#e5e5e5', 'purple')) +
#   geom_text(aes(x = year, y = count, label = round(count_perc, 1) ), 
#             #position = position_dodge(width = 0.9), 
#             vjust = -0.25) +
#   annotate('text', x = 1998, y = 155, 
#            label = paste0('Overall deduced: ', round(observation_type_val['deduced'], 2), '%\nOverall observed: ', round(observation_type_val['observed'], 2), '%'),
#            hjust = 0,
#            size = 5) +
#   theme_classic() +
#   xlab('Year') + ylab('Census count')
# 
# ggsave(plot = deduced_vs_observed_plot,
#        here('data_exploration', 'deduced_vs_observed_plot.png'),
#        width = 8*1.5, height = 4*1.5)
# 
# 
# 
# 
# 
# 
# census_gaps_barplot <- missing_census_info %>% 
#   group_by(gap_count) %>% 
#   #summarize(Count = n(), .groups = 'drop') %>% 
#   mutate(gap_count = factor(gap_count)) %>% 
#   ggplot(aes(x = gap_count)) +
#   geom_bar() +
#   geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
#   xlab('Census separation') +
#   ylab('Year') +
#   theme_classic() +
#   theme(plot.margin = margin(15, 0, 0, 0)) +
#   coord_cartesian(clip = 'off') +
#   annotate('text', x = '1', y = 600, 
#            label = ifelse(any(table(missing_census_info$RCWid) > 1),
#                           'Some individuals have\nmore than one gap',
#                           'Each individual has\nat most one gap'),
#            hjust = 0,
#            size = 5)
# 
# ggsave(plot = census_gaps_barplot,
#        here('data_exploration', 'census_gaps_barplot.png'),
#        width = 6, height = 4)
# 
# 
# 
# 
# 
# 
# 
# #when are individuals first detected?
# #difference between census 1st observation date and rcw data sheet 1st observation sheet
# 
# first_detected_summary <- summer_census_detected_addint %>% 
#   group_by(RCWid) %>% 
#   arrange(year) %>% 
#   slice_head(n = 1) %>% 
#   ungroup() %>% 
#   group_by(year) %>% 
#   summarize(count_firstdetect = n())
# 
# birth_year_summary <- rcws %>% 
#   filter(RCWid %in% summer_census_detected_addint$RCWid) %>% 
#   select(RCWid, MinAge) %>% 
#   group_by(RCWid) %>% 
#   slice_head(n = 1) %>% 
#   ungroup() %>% 
#   group_by(MinAge) %>% 
#   summarize(count_birthyear = n()) %>% 
#   rename(year = MinAge)
# 
# 
# 
# full_join(first_detected_summary, birth_year_summary, by = 'year') %>% 
#   pivot_longer(cols = c('count_firstdetect', 'count_birthyear'), 
#                names_to = 'count_type', values_to = 'count')  %>% 
#   filter(!is.na(count) & !is.na(year)) %>% 
#   ggplot() +
#   geom_point(aes(x = year, y = count, color = count_type), size = 3) +
#   theme_classic() +
#   theme(panel.grid.major.x = element_line(),
#         panel.grid.minor.x = element_line()) +
#   scale_x_continuous(minor_breaks = seq(1993, 2022, 1)) +
#   xlab('Year') + ylab('Count')
# 
# 
# 
# identical(sort(unique(first_detected_summary$year)),
#           min(first_detected_summary$year):max(first_detected_summary$year))
# 
# # first_detected_summary %>% 
# #   ggplot() +
# #   geom_point(aes(x = year, y = count), size = 3) +
# #   theme_classic() +
# #   geom_hline(yintercept = 0, color = 'red') +
# #   xlab('Year first detected') +
# #   ylab('Count') +
# #   ggtitle('Census data')
# 
#  
#   
# 
# left_join(., rcws[,c('RCWid', 'MinAge')], by = 'RCWid') %>% 
# 
# #   select(RCWid, first_census_year, MinAge) %>% 
# #   #pivot_longer(cols = c('first_census_year', 'MinAge'), names_to = 'year_type', values_to = 'year') %>% 
# 
# # first_detected_summary %>% 
# #   left_join(., rcws[,c('RCWid', 'MinAge')], by = 'RCWid')
# #   rename(first_census_year = year) %>%
# #   select(RCWid, first_census_year, MinAge) %>% 
# #   pivot_longer(cols = c('first_census_year', 'MinAge'), names_to = 'year_type', values_to = 'year') %>% 
# #   ggplot() +
# #   geom_point(aes(x = year, y = count), size = 3) +
# #   theme_classic() +
# #   geom_hline(yintercept = 0, color = 'red') +
# #   xlab('Year first detected') +
# #   ylab('Count') +
# #   ggtitle('Census data')
# # 
# 
# compare_minage_vs_censussize <- first_observation_compare %>% 
#   #filter(!is.na(dif_firstcens_minage)) %>% 
#   mutate(dif_firstcens_minage = factor(dif_firstcens_minage)) %>% 
#   #group_by(dif_firstcens_minage) %>%
#   #summarize(dif_firstcens_minage_count = n(),
#   #          .groups = 'drop') %>% 
#   ggplot(aes(x = dif_firstcens_minage)) +
#   geom_bar() +
#   geom_text(stat='count', aes(label=..count..), vjust=-1) +
#   theme_classic() +
#   theme(plot.margin = margin(15, 0, 0, 0)) +
#   coord_cartesian(clip = 'off') +
#   xlab('1st census year - minimum recorded age') +
#   ylab('Count')
# 
# ggsave(plot = compare_minage_vs_censussize,
#        here('data_exploration', 'compare_minage_vs_censussize.png'),
#        width = 6, height = 4)
# 
# 
# # first_observation_compare %>% 
# #   arrange(dif_firstcens_minage) %>% 
# #   mutate(row_ind = 1:n(),
# #          RCWid = factor(RCWid, levels = RCWid)) %>% 
# #   rename(first_census_year = year) %>% 
# #   select(RCWid, first_census_year, MinAge) %>% 
# #   #pivot_longer(cols = c('first_census_year', 'MinAge'), names_to = 'year_type', values_to = 'year') %>% 
# #   ggplot() +
# #   geom_segment(aes(x = first_census_year, xend = MinAge,
# #                    y = RCWid, yend=RCWid))
# 
# 
# 
# 
# 
#   
#   
# #characterize gaps between individuals in the census vs. individuals in the rcws data sheet
# 
# 
# #can we categorize individuals into non-migrant founders, natural migrants, residents, and translocated individuals
# 
# table(rcws$Translocation)
# table(rcws$Origin)
# 
# 
# trans_indivs <- rcws %>% 
#   filter(Origin == 'Translocated') %>% 
#   pull(RCWid)
# 
# 
# compare_minage_vs_censussize_translocated <- summer_census_detected_addint %>% 
#   filter(RCWid %in% trans_indivs) %>% 
#   group_by(RCWid) %>% 
#   arrange(year) %>% 
#   slice_head(n = 1) %>% 
#   ungroup() %>% 
#   full_join(., 
#             rcws %>% 
#               filter(Origin == 'Translocated') %>% 
#               select(RCWid, MinAge),
#             by = 'RCWid') %>% 
#   mutate(year_dif = as.factor(year - MinAge)) %>%
#   ggplot(aes(x = year_dif)) +
#   geom_bar() +
#   geom_text(stat='count', aes(label=..count..), vjust=-1) +
#   theme_classic() +
#   theme(plot.margin = margin(15, 0, 0, 0)) +
#   coord_cartesian(clip = 'off') +
#   xlab('1st census year - minimum recorded age') +
#   ylab('Count')
# 
# ggsave(plot = compare_minage_vs_censussize_translocated,
#        here('data_exploration', 'compare_minage_vs_censussize_translocated.png'),
#        width = 6, height = 4)
# 
# 
# 
# 
# rcws %>% 
#   filter(Translocation == 'Inter-population') %>% 
#   print(n = 70)
# 
# 
# rcws %>% 
#   filter(Origin == 'Translocated') %>% 
#   select(RCWid, Translocation, Sex, MinAge) %>% 
#   group_by(MinAge) %>% 
#   summarize(count = n())
# 
# 
# 
# #can we confidently say when each translocated indiv. was introduced to the population? I.e., are their gaps between
# #between when an individual is know to be introduced into the population (perhaps recorded in the rcws sheet) and
# #when it is first detected in the census?
# 
# 
# #difference between inter- and intra-population translocation? Frequency of each translocation type?
# 
# #Potential breeders?
# 
# #a cluster that contains a male and female breeder, regardless of 
# #the number of helpers, is designated as a potential breeding group
# #(PBG). Individuals may become breeders after their hatch year or 
# #remain as helpers for many years, and breeders can hold their status
# #within a PBG throughout their entire lifespan, which can be up to 16
# #years of age for males and 17 years of age for females (USFWS 2003)
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
# 
# 
# rcws_nest_query1$NatalNest
# colnames(rcws)
# left_join(rcws_wout_format %>% rename(id_nest = NatalNest), 
#           nests %>% rename(id_nest = ID),
#           by = 'id_nest')
# 
# 
# 
# 
# test_injoin <- inner_join(nests %>% rename(id_nest = ID),
#            rcws_wout_format %>% rename(id_nest = NatalNest), 
#           by = 'id_nest')
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
# rcws_nest_query1_reord <- as.data.frame(rcws_nest_query1[,c('RCWid', "MaleID", "FemaleID", "NatalNest")])
# 
# 
# rcws_ped <- process_ped(ped = test_injoin[,c('RCWid', "MaleID", "FemaleID", 'Sex')],
#                         id_col = 'RCWid',
#                         sire_col = "MaleID",
#                         dam_col = "FemaleID",
#                         sex_col= 'Sex',
#                         founder_val = NA,
#                         sex_vals = list(male = "M", female = 'F', unknown = c("U", NA)),
#                         keep_extra_cols = TRUE,
#                         disable_sex_check = FALSE)
# 
# 
# rcws_ped_query <- process_ped(ped = rcws_nest_query1 %>% 
#                                 left_join(., rcws[,c('RCWid', 'Sex')], ),
#                         id_col = 'RCWid',
#                         sire_col = "MaleID",
#                         dam_col = "FemaleID",
#                         sex_col= 'Sex',
#                         founder_val = NA,
#                         sex_vals = list(male = "M", female = 'F', unknown = c("U", NA)),
#                         keep_extra_cols = FALSE,
#                         disable_sex_check = FALSE)
# 
# rcws_ped_query %>% 
#   filter(grepl('^GB[-=]Z$', RCWid))
# 
# 
# rcws_ped_query_reord <- reorder_ped(ped = rcws_ped_query, 
#             id_col = 1, 
#             sire_col = 2, 
#             dam_col = 3)
# 
# check_order_ped(rcws_ped_query_reord, 
#                 id_col = 1, 
#                 sire_col = 2, 
#                 dam_col = 3)
# 
# rcws_nest_query1 %>% 
#   left_join(., rcws[,c('RCWid', 'Sex')], ) %>% 
#   filter(RCWid %in% get_founders(rcws_ped_query)$partial_founders)
#   
# 
# 
# summer_census_detected <- lastseendetail %>% 
#   filter(Detected == 'Yes' & SurveyType == "Census")  %>% 
#   filter(get_month(RecordDate) >= 6 & get_month(RecordDate) <= 8) %>% 
#   mutate(year = get_year(RecordDate))  %>% 
#   mutate(RCWid_old = RCWid) %>% 
#   mutate(RCWid = case_when(year < 2005 & RCWid == 'GB=Z' ~ 'GB-Z',
#                            TRUE ~ RCWid))
# 
# 
# summer_census_detected_addint <- add_intervening_years(data = summer_census_detected, 
#                                                        id_col = 'RCWid', 
#                                                        year_col = 'year')
# 
# observation_type_val <- (table(summer_census_detected_addint$row_origin)/nrow(summer_census_detected_addint))*100
# 
# summer_census_detected_addint %>% 
#   group_by(year, row_origin) %>% 
#   summarize(count = n(), .groups = 'drop') %>% 
#   group_by(year) %>%
#   mutate(count_perc = (count/sum(count))*100 ) %>% 
#   ungroup() %>% 
#   mutate(row_origin = factor(row_origin, levels = c('observed', 'deduced'))) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = row_origin),
#            #color = 'gray', size = 0.2,
#            position = "stack", stat="identity", width = 0.98) +
#   scale_fill_manual(name = 'Observation type', 
#                     values = c('#e5e5e5', 'purple')) +
#   geom_text(aes(x = year, y = count, label = round(count_perc, 1) ), 
#             #position = position_dodge(width = 0.9), 
#             vjust = -0.25) +
#   annotate('text', x = 1998, y = 155, 
#            label = paste0('Overall deduced: ', round(observation_type_val['deduced'], 2), '%\nOverall observed: ', round(observation_type_val['observed'], 2), '%'),
#            hjust = 0,
#            size = 5) +
#   theme_classic() +
#   xlab('Year') + ylab('Count')
# 
# missing_census_info %>% 
#   group_by(gap_count) %>% 
#   #summarize(Count = n(), .groups = 'drop') %>% 
#   mutate(gap_count = factor(gap_count)) %>% 
#   ggplot(aes(x = gap_count)) +
#   geom_bar() +
#   geom_text(stat='count', aes(label=..count..), vjust=-1) +
#   xlab('Census separation') +
#   ylab('Year') +
#   theme_classic() +
#   theme(plot.margin = margin(15, 0, 0, 0)) +
#   coord_cartesian(clip = 'off')
# 
# any(table(missing_census_info$RCWid) > 1)
# 
# 
# 
# 
# # summer_census_detected %>% 
# #   mutate(year = get_year(RecordDate)) %>% 
# #   group_by(year) %>% 
# #   summarize(count = n()) %>% 
# #   ggplot() +
# #   geom_point(aes(x = year, y = count))
# #   
# 
# 
# # summer_census_detected %>% 
# #   mutate(year = get_year(RecordDate)) %>% 
# #   group_by(RCWid) %>% 
# #   filter( (max(year) - min(year) + 1) != n()) %>% 
# #   filter(RCWid == 'Z-SKY') %>% 
# #   arrange(year)
# 
# 
# lastseen_rcws <- unique(lastseendetail$RCWid)
# rcws_id <- unique(rcws$RCWid)
# 
# lastseen_rcws[!lastseen_rcws %in% rcws_id]
# 
# rcws %>% 
#   filter(RCWid %in% rcws_id[!rcws_id %in% summer_census_detected$RCWid]) %>% 
#   print(n = 400)
# 
# rcws %>% 
#   filter(RCWid %in% rcws_id[!rcws_id %in% lastseen_rcws]) %>% 
#   print(n = 120)
# 
# 
# 
# 
#   
# 
# unique(lastseendetail$RCWid)[!unique(lastseendetail$RCWid) %in% rcws$RCWid]
# 
# 
# table(rcws$Fledged)
# 
# rcws_fledged <- rcws %>% 
#   filter(Fledged != 'N')
# 
# 
# unique(rcws$RCWid)[!unique(rcws$RCWid) %in% unique(lastseendetail$RCWid)]
# 
# unknowns <- unique(rcws_fledged$RCWid)[!unique(rcws_fledged$RCWid) %in% unique(lastseendetail$RCWid)]
# 
# rcws_fledged[rcws_fledged$RCWid %in% unknowns,] %>% 
#   print(n = 20)
# 
# 
# 
# 
# unique(rcws_ped$id)[!unique(rcws_ped$id) %in% unique(summer_census_detected$RCWid)]
# 
# unique(summer_census_detected$RCWid)[!unique(summer_census_detected$RCWid) %in% unique(rcws_ped$id)]
# 
# 
# 
# 
# missing_census_info <- lapply(split(summer_census_detected, summer_census_detected$RCWid), function(x) {
#   dif_vec <- sort(x$year, decreasing = TRUE)[-length(x$year)] - sort(x$year, decreasing = TRUE)[-1]
#   dif_vec1 <- dif_vec[dif_vec > 1] - 1
#   if (length(dif_vec1) == 0) dif_vec1 <- 0
#   
#   return(data.frame(gap_count = dif_vec1))
# }) %>% 
#   bind_rows(.id = 'RCWid')
# 
# 
# missing_census_info
# 
# 
# 
# (table(missing_census_info$gap_count)/length(missing_census_info$gap_count))*100
# 
# 
# missing_census_info %>% 
#   group_by(RCWid) %>% 
#   summarize(gap_count = n(), .groups = 'drop') %>% 
#   pull(gap_count) %>% 
#   max()
# 
# 
# 
# missing_census_info %>% 
#   filter(gap_count > 1)
# 
# summer_census_detected %>% 
#   filter(RCWid == 'GB=Z')
# 
# rcws %>% 
#   filter(RCWid == 'GB=Z')
# 
# rcws %>% 
#   filter(RCWid == 'GB-Z')
# 
# 
# 
# missing_census_info
# 
# 
# 
# 
# 
# test_df <- data.frame(year = c(sample(c(1, 2, 3, 4, 6, 7, 8, 9, 10)), sample(c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 15, 16)), c(3)),
#            group = c(rep('a', 9), rep('b', 12), c('c') ))
# 
# test_df_rand <- test_df[sample(nrow(test_df)),]
# rownames(test_df_rand) <- NULL
# 
# lapply(split(test_df_rand, test_df_rand$group), function(x) {
#   dif_vec <- sort(x$year, decreasing = TRUE)[-length(x$year)] - sort(x$year, decreasing = TRUE)[-1]
#   dif_vec1 <- dif_vec[dif_vec > 1] - 1
#   if (length(dif_vec1) == 0) dif_vec1 <- 0
#   
#   return(data.frame(gap = dif_vec1))
# }) %>% 
#   bind_rows(.id = 'group')
# 
# dif_vec <- sort(test_df$year, decreasing = TRUE)[-length(test_df$year)] - sort(test_df$year, decreasing = TRUE)[-1]
# dif_vec1 <- dif_vec[dif_vec > 1] - 1
# 
# 
# 
# #how many birds fledge (and thus should be out in the pop) but aren't ever in the census?
# #compare first year of sighting recorded in rcws vs. the census
# #are birds found in the year that they fledge?
# 
# 
# 
# 
# 
# 
# 
# 
# test_df <- data.frame(year = c(200, 200, 200, 201, 201, 201, 201),
#            type = c('a', 'b', 'a', 'b', 'b', 'b', 'a'))
# 
# test_df_rug <- test_df %>% 
#   filter(type == 'a') %>% 
#   mutate(color = c('red', 'blue', 'red'))
# 
# test_df %>%
#   group_by(year, type) %>% 
#   summarize(count = n(), .groups = 'drop') %>% 
#   ggplot(aes(fill=type, y=count, x=year)) + 
#   geom_bar(position="stack", stat="identity") +
#   geom_rug(data = test_df_rug,
#            aes(x = year, y = 0), sides = 'b', size= 1.5,
#            color = test_df_rug$color,
#            position = position_jitterdodge(dodge.width = 0,
#                                          jitter.width = 0.1,
#                                          jitter.height = 0))
# 
# 
# 
# 
#   geom_jitter(data = test_df %>% filter(type == 'a'),
#               aes(y = - 0.1),
#               size = 1.5, alpha = 0.7,
#               position = position_jitterdodge(dodge.width = 0,
#                                               jitter.width = 0,
#                                               jitter.height = 0.5))
# 
# raw_dat_list$nests %>% 
#   group_by(Year) %>% 
#   summarize(n())
# 
# table(raw_dat_list$nests$Year)
# 
# length(raw_dat_list$nests$ID)
# length(unique(raw_dat_list$nests$ID))
# 
# raw_dat_list$nests$Cluster
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
# summer_census_detected <- lastseendetail %>% 
#   filter(Detected == 'Yes' & SurveyType == "Census")  %>% 
#   filter(get_month(RecordDate) >= 6 & get_month(RecordDate) <= 8) %>% 
#   mutate(year = get_year(RecordDate))  %>% 
#   mutate(RCWid_old = RCWid) %>% 
#   mutate(RCWid = case_when(year < 2005 & RCWid == 'GB=Z' ~ 'GB-Z',
#                            TRUE ~ RCWid))
# 
# summer_census_detected_addint <- add_intervening_years(data = summer_census_detected, 
#                                                        id_col = 'RCWid', 
#                                                        year_col = 'year') %>% 
#   left_join(., 
#             rcws %>% 
#               select(RCWid, Origin),
#             by = 'RCWid')
# 
# summer_census_detected_addint_join <- left_join(summer_census_detected_addint, 
#                                                 rcws %>% 
#                                                   select(RCWid, Origin))
# 
# 
# 
# census_minage1 <- summer_census_detected_addint %>% 
#   select(RCWid, year) %>% 
#   mutate(type = 'census') %>% 
#   bind_rows(.,
#             rcws %>% 
#               filter(RCWid %in% summer_census_detected_addint$RCWid) %>% 
#               select(RCWid, MinAge) %>% 
#               rename(year = MinAge) %>% 
#               filter(!is.na(year)) %>% 
#               mutate(type = 'minage')
#             )
# 
# 
# 
# 
# 
# 
# 
# 
# #if (sum(test_subset$SurveyType == 'birth_info') > 1)
# #  warning('There are more than one birth_info entries for ', IND)
# 
# 
# #Four potential filtering 
# 
# #1: CENSUS_ONLY
# #Only census observations
# 
# #YEAR1    YEAR2    YEAR3    YEAR4    YEAR5
# # B                  X                 X
# #include Y3, Y5
# 
# 
# #2: CENSUS + DEDUCED
# #Census observations and filling in the the gaps between non-consecutive census observations
# 
# #YEAR1    YEAR2    YEAR3    YEAR4    YEAR5
# # B                  X                 X
# #include Y3, Y4, Y5
# 
# 
# #3: CENSUS/BIRTH + DEDUCED (INCLUDE BIRTH YEAR IF FLEDGED, BORN IN THE POP, AND FOUND IN A SUBSEQUENT CENSUS)
# #Similar to 2. Include all census observations but for birds born in the population, record 
# #the start year as the year where the 
# 
# #YEAR1    YEAR2    YEAR3    YEAR4    YEAR5
# # B                                   
# #include none if found in year1, don't include if never refound
# 
# #YEAR1    YEAR2    YEAR3    YEAR4    YEAR5
# # B                  X                 X
# #include Y1, Y2, Y3, Y4, Y5
# 
# 
# #4: CENSUS/BIRTH + DEDUCED (EXCLUDE BIRTH YEAR EVEN IF FLEDGED AND FOUND IN A SUBSEQUENT CENSUS)
# #Identical to 3 but exclude birth year. This would represent a focus on the breeder population
# 
# #YEAR1    YEAR2    YEAR3    YEAR4    YEAR5
# # B                  X                 X
# #include Y3, Y5
# 
# julian(raw_dat_list$nests$StartDate)
# 
# test_nest <- raw_dat_list$nests %>% 
#   mutate(julian = lubridate::yday(raw_dat_list$nests$StartDate) ) %>% 
#   filter(!is.na(julian) & !is.na(HatchNum))
# 
# 
# raw_dat_list$nests$
# 
# summary(lm(test_nest$julian ~ test_nest$Year))
# summary(lm(test_nest$HatchNum ~ test_nest$julian))
# summary(glm(test_nest$HatchNum ~ test_nest$julian, family = 'poisson'))
# 
# raw_dat_list$nests %>% 
#   mutate(julian = lubridate::yday(raw_dat_list$nests$StartDate) ) %>% 
#   filter(!is.na(julian) & !is.na(HatchNum)) %>% 
#   ggplot(aes(y = HatchNum, x = julian)) +
#   geom_point() +
#   geom_smooth(method='lm', formula = y~x)

