##############################################################################################
### SCRIPT NAME: translocation_analysis_viz.R
### PURPOSE: Summarize and visualize the performances of the translocated birds
### PRODUCTS:
###     transloc_explore_multipanel1.png: Fig. 1 in paper main text
###     translocation_info_table.txt: table summarizing translocations included
###                                   in supplement
###     translocation_outcomes_info.txt: file containing info. on translocations
###                                      reported in the main text
##############################################################################################

#####################
### SCRIPT SET-UP ###
#####################

### PACKAGES ###
library(tidyverse)
library(readxl)
library(here)
library(ggnetwork)
library(network)
library(clock)
library(see)
library(kableExtra)
#library(sna)
#library(egg)

source(here('scripts', 'ped_functions.R'))


here::here()


### LOADING DATA AND RESULTS ###
#population composition data
census_processed <- read.csv(here('data', 'feb2024_databasemarch2023_processed', 'census_processed.csv'))

translocations <- read_excel(here('data', 'feb2024_databasemarch2023', 'Translocation.xlsx'))
nests <- read_excel(here('data', 'feb2024_databasemarch2023', 'Nests.xlsx'))
rcws <- read_excel(here('data', 'feb2024_databasemarch2023', 'RCWs.xlsx'))

results_names <- setNames(nm = c("rcws_founder_info",
                                 "rcws_inbr", 
                                 "transloc_info_color", 
                                 "ne_f_df", "rcw_inbr_merge", 
                                 "rcws_ancestry_info",
                                 "rcws_ancestry_source_by_year",
                                 "rcw_partial_founder_summed_processed",
                                 "rcw_contr_info_df",
                                 "fped_prop_by_group",
                                 "inbr_founder_color")
)
results_list <- lapply(results_names, function(x) read.csv(here('results', paste0(x, '.csv'))))

#results for mark--recapture model of translocated vs. non-translocated survival
cjs_translocation_top_mod <- read.csv(here('results', 'cjs_translocation_top_mod.csv'))

transloc_mod_plot_info <- readRDS(file = here('results', 'transloc_mod_plot_info.RDS'))


### aesthetic info for plotting
group_color_vec <- setNames(c("#f0597d", "#ffdb19", "#595959", "#06d6a0", "#118ab2", "#800080", "#ededed"),
                            nm = c("ANF", "CBJTC", "FTB", "FTS", "ONF", "WSF-CITRUS", "non-transloc"))



#######################
### DATA PROCESSING ###
#######################

#POPULATION COMPOSITION
pop_info_list <- lapply(setNames(nm = paste0('Scenario', 1:4)), function(X, DF) {
  DF[DF[,X,drop = TRUE],!grepl(pattern = 'Scenario', colnames(DF))]
}, DF = census_processed)
pop_dat <- pop_info_list$Scenario3


file_names <- list.files(here('data', 'feb2024_databasemarch2023_processed'))
rcw_processed_list <- lapply(setNames(file_names[grepl('csv$', file_names)], nm = gsub(pattern = '\\.csv', '', file_names[grepl('csv$', file_names)])), function(x) {
  read.csv(here('data', 'feb2024_databasemarch2023_processed', x))
})


#info on translocations
translocations_inter <- translocations[translocations$Type == "Inter-population",]


nest_count_summary <- nests %>% 
  pivot_longer(cols = c('MaleID', 'FemaleID'), names_to = 'parent_type', values_to = 'RCWid') %>% 
  filter(RCWid %in% translocations_inter$RCWid) %>% 
  group_by(RCWid, Year) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  group_by(RCWid) %>% 
  summarize(nest_years = n())


#how long were they in the population

translocations_inter1 <- pop_dat %>% 
  filter(RCWid %in% translocations_inter$RCWid) %>% 
  group_by(RCWid) %>% 
  summarize(year_count = n(),
            #first_year = min(Year),
            #last_year = max(Year),
            first_year = min(year),
            last_year = max(year),
            .groups = 'drop') %>% 
  right_join(., translocations_inter,
             by = 'RCWid') %>% 
  mutate(year_count = if_else(is.na(year_count), 0, year_count)) %>% 
  left_join(., rcws, by = 'RCWid') %>% 
  left_join(., results_list$transloc_info_color[,c('RCWid', 'sex_color', 'source_color')], by = 'RCWid') %>% 
  left_join(., nest_count_summary, by = 'RCWid') %>% 
  mutate(nest_years = if_else(is.na(nest_years), 0, nest_years))

sex_col_list <- list()
source_col_list <- list()

for (i in unique(translocations_inter1$Sex)) {
  rgb <- col2rgb(translocations_inter1$sex_color[translocations_inter1$Sex == i])
  lab <- convertColor(t(rgb), 'sRGB', 'Lab')
  
  sex_col_list[[i]] <- data.frame(sex_color = translocations_inter1$sex_color[translocations_inter1$Sex == i][order(lab[, 'L'])],
                              sex_order = seq_along(translocations_inter1$sex_color[translocations_inter1$Sex == i])
                              )
}

for (i in unique(translocations_inter1$Source)) {
  rgb <- col2rgb(translocations_inter1$source_color[translocations_inter1$Source == i])
  lab <- convertColor(t(rgb), 'sRGB', 'Lab')
  
  source_col_list[[i]] <- data.frame(source_color = translocations_inter1$source_color[translocations_inter1$Source == i][order(lab[, 'L'])],
                             source_order = seq_along(translocations_inter1$source_color[translocations_inter1$Source == i])
  )
}


translocations_inter2 <- translocations_inter1 %>% 
  left_join(., bind_rows(sex_col_list), by = 'sex_color') %>% 
  left_join(., bind_rows(source_col_list), by = 'source_color') %>% 
  arrange(Source, source_order) %>% 
  mutate(id_factor_source = factor(RCWid, levels = RCWid)) %>% 
  arrange(Sex, sex_order) %>% 
  mutate(id_factor_sex = factor(RCWid, levels = RCWid))




# nest_info_processed_barplot <- transloc_nests %>%
#   mutate(female_status = if_else(FemaleID %in% translocations_inter$RCWid, 'transloc', 'nontransloc'),
#          male_status = if_else(MaleID %in% translocations_inter$RCWid, 'transloc', 'nontransloc')) %>% 
#   rowwise() %>% 
#   mutate(mate_combo = paste(sort(c(female_status, male_status)), collapse = '\u2013')) %>% 
#   ungroup() %>% 
#   group_by(mate_combo) %>% 
#   summarize(Nests = n(),
#             Fledglings = sum(fldg_count, na.rm = TRUE)
#   ) %>% 
#   pivot_longer(cols = c('Nests', 'Fledglings'), names_to = 'var', values_to = 'val') %>% 
#   group_by(var) %>% 
#   mutate(val_perc = paste0(round((val/sum(val)) * 100, 1), '%')) %>% 
#   ungroup() %>% 
#   mutate(val_plot = if_else(mate_combo == 'nontransloc\u2013transloc', -1*val, val))



transloc_hist_multipanel <- translocations_inter2 %>% 
  pivot_longer(cols = c(year_count, nest_years),
               names_to = 'variable',
               values_to = 'value') %>% 
  mutate(variable_updated = factor(case_when(variable == 'nest_years' ~ 'Nesting years',
                                      TRUE ~ 'Years in population'),
                                   levels = c('Years in population', 'Nesting years')) ) %>% 
  ggplot() +
  geom_bar(aes(x = value, fill = id_factor_source)) +
  scale_fill_manual(values = setNames(results_list$transloc_info_color[[paste0('source', '_color')]],
                                      nm = results_list$transloc_info_color$RCWid)) +
  theme_bw() +
  theme(#strip.background = element_blank(), 
        #strip.text = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white'), 
    strip.text = element_text(size = 13, face = 'bold'),
    axis.text = element_text(size = 10),
        legend.position = 'none',
        panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = '#4c4c4c', linewidth = 1.5),
        panel.border = element_blank(),
        panel.spacing.y = unit(1, "lines"),
    axis.title.x = element_text(size = 13),
    panel.background = element_rect(fill = 'transparent', color = 'transparent')
        #axis.line.x.top = element_blank(),
        #axis.ticks.x.top = element_blank(),
        #axis.title.x.top = element_blank()
    ) +
  facet_wrap(~variable_updated,
             scales = 'free_y',
             nrow = 1) +
  xlab('Count') +
  ylab('Number of individuals') +
  coord_cartesian(clip = 'off') +
  #geom_text(aes(label = variable_updated), x = 1.5, y = Inf, 
  #          hjust = 0, vjust = 0.5, check_overlap = TRUE,
  #          position = position_nudge(x = 100), size = 5,
  #          fontface = 'bold') +
  scale_y_continuous(expand = c(0, 0))#+
  #scale_x_continuous(sec.axis = dup_axis())



# transloc_nests1 <- nests %>% 
#   #pivot_longer(cols = c('MaleID', 'FemaleID'), names_to = 'parent_type', values_to = 'RCWid') %>% 
#   filter(MaleID %in% translocations_inter$RCWid | FemaleID %in% translocations_inter$RCWid) %>% 
#   select(FemaleID, MaleID) %>% 
#   group_by(FemaleID, MaleID) %>% 
#   slice_head(n = 1) %>% 
#   mutate(val = 1L) %>% 
#   ungroup()
# 
# test_mat <- matrix(rep(0, length(c(transloc_nests1$FemaleID, transloc_nests1$MaleID))^2),
#        nrow = length(c(transloc_nests1$FemaleID, transloc_nests1$MaleID)))
# colnames(test_mat) <- c(transloc_nests1$FemaleID, transloc_nests1$MaleID)
# rownames(test_mat) <- c(transloc_nests1$FemaleID, transloc_nests1$MaleID)


transloc_nests <- nests %>% 
  #pivot_longer(cols = c('MaleID', 'FemaleID'), names_to = 'parent_type', values_to = 'RCWid') %>% 
  filter(MaleID %in% translocations_inter$RCWid | FemaleID %in% translocations_inter$RCWid) %>% 
  select(FemaleID, MaleID, FldgNum) %>% 
  group_by(FemaleID, MaleID) %>% 
  summarize(weight = n(),
            fldg_count = sum(FldgNum, na.rm = TRUE),
            .groups = 'drop') %>% 
  mutate(fldg_yes = if_else(fldg_count > 0, 'yes', 'no'))

# n <- network(transloc_nests, directed = FALSE)

# simple_edge_df <- data.frame(
#   from = c("b", "c", "c", "d", "a"),
#   to = c("a", "b", "a", "a", "b"),
#   weight = c(1, 1, 2, 2, 3),
#   stringsAsFactors = FALSE
# )

# simple_vertex_df <- data.frame(
#   name = letters[1:5],
#   residence = c("urban", "rural", "suburban", "suburban", "rural"),
#   stringsAsFactors = FALSE
# )


test_vertex <- transloc_nests %>% 
  pivot_longer(cols = c("MaleID", "FemaleID"), names_to = 'ID') %>% 
  mutate(sex = if_else(ID == 'MaleID', 'M', 'F'),
         transloc = if_else(value %in% translocations_inter$RCWid, 'yes', 'no')) %>% 
  rename(name = value) %>% 
  select(name, sex, transloc) %>% 
  group_by(name) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  left_join(., results_list$transloc_info_color %>% 
              rename(name = RCWid) %>% 
              select(name, source_color),
            by = 'name') %>% 
  mutate(source_color = if_else(is.na(source_color), 'gray', source_color),
         name1 = name)


#add color values for non-translocated birds
updated_color_df <- results_list$transloc_info_color %>% 
  select(RCWid, source_color) %>% 
  bind_rows(., 
            data.frame(RCWid = test_vertex$name[test_vertex$transloc == 'no'],
                       source_color = '#b9b9b9')
  )


test_net <- as.network(transloc_nests, vertices = test_vertex,
                       directed = FALSE)

#set.seed(91091)
set.seed(91096)
test_net_plot <- ggplot(test_net, 
       aes(x = x, y = y, xend = xend, yend = yend)) +
  #geom_edges(aes(linewidth = fldg_count,
  #               linetype = fldg_yes), color = "#dcdcdc") +
  geom_edges(aes(linewidth = fldg_count,
                 linetype = fldg_yes), color = "#595959") +
  geom_nodes(aes(color = name1, shape = sex), 
             size = 5.8) +
  theme_blank() +
  scale_linewidth_continuous(name = 'Fledgling\ncount',
                             breaks = c(1, 6, 12), 
                             range = c(0.75, 4.25)) +
  scale_linetype_manual(values = c('dotted', 'solid'),
                        guide = "none") +
  #scale_size_manual(values = c(2.5, 5.5)) +
  scale_color_manual(values = setNames(updated_color_df$source_color,
                                      nm = updated_color_df$RCWid),
                     guide = "none") +
  scale_shape_manual(name = 'Sex',
                     values = c(19, 15)) + #female: circle; male: square
  theme(
    #legend.position = 'none',
        #plot.margin = unit(c(-0.08, 8.1, 0.48, 0.45), "cm"), #original
        plot.margin = unit(c(0.25, 9, 1, 0.9), "cm"),
        #legend.position = c(0.94, .1),
        legend.position = c(0.09, 0.06),
        legend.background = element_rect(fill = NA, color = NA),
        legend.key = element_rect(fill = NA, color = NA),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.direction = "vertical", 
        legend.box = "horizontal",
        legend.spacing.x = unit(0.02, "cm")) +
  coord_cartesian(clip = 'off') +
  guides(shape = guide_legend(order = 1),
         linewidth = guide_legend(order = 2))

# test_net_plot1 <- test_net_plot + 
#   annotation_custom(
#     ggplotGrob(establishment_barplot),
#     xmin = 1.01, xmax = 1.43, ymin = 0.52, ymax = 1.02
#   ) +
#   annotation_custom(
#     ggplotGrob(examp_bar),
#     xmin = 1.03, xmax = 1.46, ymin = -0.04, ymax = 0.48
#   )
  
# test_net_plot1 <- test_net_plot + 
#   annotation_custom(
#     ggplotGrob(establishment_barplot),
#     xmin = 1.01, xmax = 1.43, ymin = 0.52, ymax = 1.02
#   ) +
#   annotation_custom(
#     ggplotGrob(test_compare),
#     xmin = 1.15, xmax = 1.46, ymin = -0.04, ymax = 0.48
#   )


# transloc_hist_multipan <- cowplot::plot_grid(year_hist,
#                    nest_hist,
#                    nrow = 2)
# 
# 
# transloc_explore_multipanel <- cowplot::plot_grid(test_net_plot1,
#                    transloc_hist_multipanel +
#                      theme(plot.margin = unit(c(0, 0.4, 0, 0.4), "cm")),
#                    ncol = 1,
#                    rel_heights = c(0.67, 0.33))
# 
# 
# transloc_explore_multipanel_label <- transloc_explore_multipanel +
#   annotate("text", 
#            x = c(0.025, 0.75, 0.75, 0.025, 0.503), 
#            y = c(0.955, 0.955, 0.58, 0.32, 0.32), 
#            label = LETTERS[1:5], 
#            fontface = 2, size = 7.25,
#            colour= "black",
#            angle = c(0), hjust = 0.5)
# 
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'transloc_explore_multipanel.png'),
#                  plot = transloc_explore_multipanel_label,
#                  width = 10*1, height = 9.5*1, bg = 'white')




#Pairings: transloc--transloc vs. nontransloc--transloc

#nests
#nestlings


# simple_bar_df <- data.frame(prop = c(0.4, -0.6, 0.8, -0.2),
#            group = c('Nests', 'Nests', 'Fledglings', 'Fledglings'),
#            group1 = c('1', '2', '1', '2'),
#            label = c('40%', '60%', '80%', '20%'))
# 
# examp_bar <- nest_info_processed_barplot %>% 
#   mutate(var = factor(var, levels = c('Nests', "Fledglings"))) %>% 
#   ggplot(aes(x = var, y = val_plot, fill = mate_combo)) +
#   geom_bar(stat = "identity", position = "identity") +
#   geom_hline(yintercept = 0, color = 'black', linewidth = 1) +
#   geom_text(aes(label = val_perc),
#             nudge_y = c(-20, -20, 20, 20),
#             size = 4,
#             color = 'black') +
#   theme_bw() +
#   theme(axis.ticks.x = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 12, face = 'bold'),
#         axis.line.y = element_line(color = 'black', linewidth = 1),
#         legend.position = 'none',
#         panel.background = element_rect(fill = 'transparent', color = 'transparent'),
#         plot.background = element_rect(fill = 'transparent', color = 'transparent')) + 
#   annotate("text", y = c(100, -100), x = 2.65, 
#            label = c("transloc\u2013transloc","transloc\u2013nontransloc"), 
#            fontface=1, size = 4,
#            colour= "black",
#            angle = c(90 + 180), hjust = 0.5) +
#   coord_cartesian(clip = 'off') +
#   ylab('Count') +
#   scale_x_discrete(position = "top") +
#   #ylim(-190, 190) +
#   scale_y_continuous(limits = c(-190, 185),
#                      breaks = c(-160, -120, -80, -40, 0, 40, 80, 120, 160),
#                      labels = as.character(abs(c(-160, -120, -80, -40, 0, 40, 80, 120, 160)))) +
#   scale_fill_manual(values = c("#bfbfbf", "#595959"))



# examp_bar <- nest_info_processed_barplot %>% 
#   mutate(var = factor(var, levels = c('Nests', "Fledglings"))) %>% 
#   ggplot(aes(x = var, y = val_plot, fill = mate_combo)) +
#   geom_bar(stat = "identity", position = "identity") +
#   geom_hline(yintercept = 0, color = 'black', linewidth = 1) +
#   geom_text(aes(label = val_perc),
#             nudge_y = c(-15, -15, 15, 15),
#             size = 4,
#             color = 'black') +
#   theme_bw() +
#   theme(axis.ticks.x = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 13.5, face = 'bold'),
#         axis.line.y = element_line(color = 'black', linewidth = 1),
#         legend.position = 'none',
#         panel.background = element_rect(fill = 'transparent', color = 'transparent'),
#         plot.background = element_rect(fill = 'transparent', color = 'transparent')) + 
#   annotate("text", y = c(50, -210), x = 1.5, 
#            label = c("transloc\u2013transloc","transloc\u2013nontransloc"), 
#            fontface = 1, size = 4.5,
#            colour= "black",
#            angle = c(0), hjust = 0.5) +
#   coord_cartesian(clip = 'off') +
#   ylab('Count') +
#   scale_x_discrete(position = "top") +
#   #ylim(-190, 190) +
#   scale_y_continuous(limits = c(-210, 60),
#                      breaks = c(-180, 160, -120, -80, -40, 0, 40),
#                      labels = as.character(abs(c(-180 , -160, -120, -80, -40, 0, 40)))) +
#   scale_fill_manual(values = c("#b9b9b9", "#595959"))
# 
# #+
#   #coord_flip()



# source_test_df <- data.frame(source = c('S1', 'S1', 'S2', 'S2' , 'S3', 'S3'),
#                              established = c('Yes', 'No', 'Yes', 'No', "Yes", "No"),
#                              count = c(3, -1, 4, 0, 2, -6))

#source_test_df %>% 
  
establishment_barplot <- translocations_inter2 %>% 
  mutate(establishment = if_else(year_count > 0, 'Yes', 'No'),
         Source = factor(Source, levels = rev(c('ONF', 'FTS', 'WSF-CITRUS', 'FTB', 'ANF', 'CBJTC')))) %>% 
  group_by(Source, establishment) %>% 
  summarize(Count = n(),
            .groups = 'drop') %>% 
  mutate(count_update = if_else(establishment == 'Yes', -1*Count, Count)) %>% 
  ggplot(aes(x = Source, y = count_update 
             )) +
  geom_bar(fill = 'white',
           stat = "identity", position = "identity", linewidth = 1) +
  geom_bar(aes(fill = Source, 
               #color = source,
               alpha = establishment), 
           stat = "identity", position = "identity", linewidth = 1) +
  geom_hline(yintercept = 0, color = "#4c4c4c", linewidth = 1.5) +
  coord_flip(clip = 'off') +
  theme_bw() +
  scale_alpha_manual(values = c(0.24, 1), guide = 'none') +
  theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = 'transparent'), 
        axis.text.x = element_text(size = 9),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "#4c4c4c", linewidth = 1.5),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold',
                                  margin=unit(c(0, 0, 0.65, 0), "cm")
                                  ) 
        )+
  ggtitle('Establishment by source')+
  annotate("text", y = c(-1.5, 1.5), x = c(6.5, 6.5),
           hjust = c(1, 0),
           label = c('Yes', 'No'), 
           fontface = 1, size = 4.5,
           colour= "black",
           angle = c(0), 
           #hjust = 0.5, 
           vjust = - 0.5) +
  ylab('Count') +
  scale_y_continuous(breaks = c(-15, -10, -5, 0, 5),
                     labels = as.character(abs(c(-15, -10, -5, 0, 5)))) +
  scale_fill_manual(values = group_color_vec)


# year_cohort_df <- data.frame(year = c(1994, 1994, 2005, 2005, 2008, 2009, 2009, 2012, 2012, 2016),
#            count = c(5,   6,    3,    1,    7,    3,    4,    2,    1,    8),
#            establish = c('Y', 'N', 'Y', 'N', 'Y', 'Y', 'N', 'Y', 'N', 'Y'),
#            source = c('a', 'a', 'b', 'b', 'c', 'a', 'a', 'a', 'a', 'b'))
# 
# year_cohort_df_update <- rbind(year_cohort_df,
#                                data.frame(year = c(1994:2022)[!c(1994:2022) %in% year_cohort_df$year],
#                                           count = 0,
#                                           establish = 'Y',
#                                           source = 'a')
#                                )

# translocations_inter2 %>% 
#   mutate(year = clock::get_year(TranslocationDate),
#          establishment = if_else(year_count > 0, 'Yes', 'No')) %>% 
#   group_by(year, establishment) %>% 
#   summarize(indiv_count = n(),
#             Source = first(Source),
#             .groups = 'drop') %>% 
#   mutate(
#     Source = factor(Source, levels = rev(c('ONF', 'FTS', 'WSF-CITRUS', 'FTB', 'ANF', 'CBJTC'))),
#     indiv_count_update = if_else(establishment == 'Yes', indiv_count, indiv_count*-1)
#   ) %>% 
#   rbind(.,
#         data.frame(year = c(1994:2022)[!c(1994:2022) %in% clock::get_year(translocations_inter2$TranslocationDate)],
#                    establishment = NA,
#                    indiv_count = 0,
#                    Source = NA,
#                    indiv_count_update = 0)
#   ) %>% 
#   print(n = 50)


year_establishment_barplot <- translocations_inter2 %>% 
  mutate(year = clock::get_year(TranslocationDate),
         establishment = if_else(year_count > 0, 'Yes', 'No')) %>% 
  group_by(year, establishment) %>% 
  summarize(indiv_count = n(),
            Source = first(Source),
            .groups = 'drop') %>% 
  mutate(
    Source = factor(Source, levels = rev(c('ONF', 'FTS', 'WSF-CITRUS', 'FTB', 'ANF', 'CBJTC'))),
    indiv_count_update = if_else(establishment == 'Yes', indiv_count, indiv_count*-1)
    ) %>% 
  # rbind(.,
  #       data.frame(year = c(1994:2022)[!c(1994:2022) %in% clock::get_year(translocations_inter2$TranslocationDate)],
  #                  establishment = 'Yes',
  #                  indiv_count = 0,
  #                  Source = 'ANF',
  #                  indiv_count_update = 0)
  #       ) %>% 
  ggplot(aes(x = year, y = indiv_count_update 
  )) +
  geom_bar(fill = 'white',
           stat = "identity", position = "identity", linewidth = 1) +
  geom_bar(aes(fill = Source, 
               #color = source,
               alpha = establishment),
           stat = "identity", position = "identity", linewidth = 1) +
  geom_hline(yintercept = 0, 
             color = "#4c4c4c", 
             linewidth = 1.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.grid.minor = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = 'none',
        axis.text.y = element_text(size = 10),
        axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
        plot.title = element_text(hjust = 0.5, size = 13.5, face = 'bold',
                                  margin=unit(c(0, 0, 0.65, 0), "cm")
        ) ,
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.title.x.top = element_blank(),
        axis.ticks.x = element_blank()
  ) +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(sec.axis = dup_axis(),
                     limits = c(1994, 2022),
                     breaks = c(2000, 2005, 2010, 2015, 2020),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6),
                     labels = abs(c(-4, 2, 0, 2, 4, 6)),
                     limits = c(-4.1, 6.1)
                     ) +
  ylab('Count') +
  scale_alpha_manual(values = c(0.24, 1), guide = 'none') +
  scale_fill_manual(values = group_color_vec) +
  annotate("text", 
           x = 1994.5, 
           y = 6, 
           label = 'Establishment\nby year',
           lineheight = .85,
           fontface = 2, size = 4,
           colour= "black",
           angle = c(0), hjust = 0) +
  annotate("text", 
           x = 1994.5, 
           y = 1.75, 
           label = 'Yes', 
           fontface = 1, size = 4.35,
           colour= "black",
           angle = c(0), hjust = 0) +
  annotate("text", 
           x = 1994.5, 
           y = -1.75, 
           label = 'No', 
           fontface = 1, size = 4.35,
           colour= "black",
           angle = c(0), hjust = 0) +
  annotate("text", 
           x = 2021.5, 
           y = 1.2, 
           label = 'Year', 
           fontface = 1, size = 3.9,
           colour= "black",
           angle = c(0), hjust = 0.5)#+
  #scale_x_continuous(limits = c(1994, 2022),
  #                   expand = c(0, 0)) #+
  #scale_y_continuous(expand = c(0, 0))





cjs_transloc_plot <- cjs_translocation_top_mod %>%
  mutate(sex_update = if_else(sex == 0, 'female', 'male'),
         transloc_update = factor(if_else(transloc == 0, 'Non-translocated\n(0% transloc. ancestry)', 'Translocated'),
                                  levels = c('Translocated', 'Non-translocated\n(0% transloc. ancestry)'))) %>%
  filter(param_type == 'Phi') %>%
  ggplot(aes(x = sex_update, y = estimate)) +
  # geom_segment(aes(x = level2, y = low95, yend = high95, color = level1), 
  #              linewidth = 6, lineend = 'round') +
  # geom_point(aes(shape = level2), 
  #            size = 3.5, color = 'white') +
  geom_segment(aes(x = sex_update, y = lcl, yend = ucl, color = transloc_update), 
               linewidth = 3.7, lineend = 'round', alpha = 1) +
  geom_point(aes(shape = sex_update, fill = transloc_update), 
             size = 6, stroke = 0) + 
  scale_shape_manual(name = 'Sex',
                     values = c(21, 22)) + #female: circle; male: square
  facet_wrap(~transloc_update, strip.position = "top", nrow = 1) +
  scale_color_manual(values = c('#818181', "#d8d8d8")) +
  scale_fill_manual(values = c('#4c4c4c', "#acacac")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(size = 10, vjust = 1,
                                  margin = margin(t = 0, b = -8)),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = 'transparent', color = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = 'transparent'), 
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "#4c4c4c", linewidth = 1.5),
        #axis.ticks.y = element_blank(),
        legend.position = 'none',
        axis.text.y = element_text(size = 10),
        strip.background = element_blank()
  ) +
  ylab('Apparent survival') +
  coord_cartesian(clip = 'off')


repro_mod_compare <- transloc_mod_plot_info$mod_plotting %>% 
  #mutate(transloc = factor(transloc, levels = c('no', 'yes'))) %>% 
  mutate(transloc_final = factor(case_when(transloc == 'no' ~ 'Non-translocated\n(0% transloc. ancestry)',
                                           transloc == 'yes' ~ 'Translocated'),
                                 levels = c('Translocated', 'Non-translocated\n(0% transloc. ancestry)'))) %>% 
  ungroup() %>% 
  ggplot()  +
  geom_point(data = transloc_mod_plot_info$raw_dat %>% 
               #mutate(transloc = factor(transloc, levels = c('no', 'yes'))) %>% 
               mutate(transloc_final = factor(case_when(transloc == 'no' ~ 'Non-translocated\n(0% transloc. ancestry)',
                                                        transloc == 'yes' ~ 'Translocated'),
                                              levels = c('Translocated', 'Non-translocated\n(0% transloc. ancestry)'))),
           aes(x = transloc_final, y = life_fldg, color = dummy, fill = dummy, shape = sex),
           position = ggpp::position_jitternudge(seed = 42,
                                                 width = 0.19,
                                                 height = 0,
                                                 x = -0.24,
                                                 nudge.from = 'jittered'),
           #color = "#b9b9b9",
           alpha = 0.6,
           size = 2.5) +
  geom_violinhalf(aes(x = transloc_final, y = .epred),
                  linewidth = 1.1, fill = "#595959", color = "#595959") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(color = "#4c4c4c", linewidth = 1.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  #xlab("Translocated") +
  ylab("Life repro. success") +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_color_manual(values = c("#b9b9b9", "#595959")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual(values = setNames(transloc_mod_plot_info$color_info$source_color,
                                       nm = transloc_mod_plot_info$color_info$dummy)) +
  scale_fill_manual(values = setNames(transloc_mod_plot_info$color_info$source_color,
                                      nm = transloc_mod_plot_info$color_info$dummy)) +
  scale_shape_manual(values = setNames(c(22, 21),
                                       nm = c('fid', 'mid')))



test_net_plot1 <- test_net_plot + 
  annotation_custom(
    ggplotGrob(establishment_barplot),
    #xmin = 0.99, xmax = 1.518, ymin = 0.485, ymax = 1.015
    xmin = 0.99 + 0.05, xmax = 1.518 + 0.06, ymin = 0.67, ymax = 1.07
  ) +
  annotation_custom(
    ggplotGrob(cjs_transloc_plot),
    #xmin = 1.045, xmax = 1.525, ymin = -0.1, ymax = 0.43
    xmin = 1.045 + 0.05, xmax = 1.525 + 0.06, ymin = 0.26, ymax = 0.66
  ) +
  annotation_custom(
    ggplotGrob(repro_mod_compare),
    #xmin = 1.045, xmax = 1.525, ymin = -0.1, ymax = 0.43
    xmin = 1.045 + 0.05, xmax = 1.525 + 0.06, ymin = -0.17, ymax = 0.23
  )


# test_net_plot1 <- test_net_plot + 
#   annotation_custom(
#     ggplotGrob(establishment_barplot),
#     #xmin = 0.99, xmax = 1.518, ymin = 0.485, ymax = 1.015
#     xmin = 0.99 + 0.05, xmax = 1.518 + 0.05, ymin = 0.485, ymax = 1.015
#   ) +
#   annotation_custom(
#     ggplotGrob(repro_mod_compare),
#     #xmin = 1.045, xmax = 1.525, ymin = -0.1, ymax = 0.43
#     xmin = 1.045 + 0.05, xmax = 1.525 + 0.05, ymin = -0.1, ymax = 0.43
#   )



transloc_explore_multipanel <- cowplot::plot_grid(year_establishment_barplot +
                                                    theme(plot.margin = unit(c(0.5, 0.4, 0.15, 0.4), "cm")),
                                                  test_net_plot1,
                                                  transloc_hist_multipanel +
                                                    theme(plot.margin = unit(c(0.22, 0.4, -0.09, 0.4), "cm"),
                                                          plot.background = element_rect(fill = 'transparent', color = 'transparent')),
                                                  nrow = 3,
                                                  rel_heights = c(0.2, 0.75, 0.23))


transloc_explore_multipanel_label <- transloc_explore_multipanel +
  annotate("text", 
           x = c(0.025, 0.025, 0.71, 0.71, 0.71, 0.025, 0.506), 
           y = c(0.985, 0.80, 0.80, 0.61, 0.405, 0.20, 0.20), 
           label = LETTERS[1:7], 
           fontface = 2, size = 7.45,
           colour= "black",
           angle = c(0), hjust = 0.5)

cowplot::ggsave2(filename = here('figures', 'main_paper', 'transloc_explore_multipanel1.png'),
                 plot = transloc_explore_multipanel_label +
                   theme(plot.margin = unit(c(0.1, 0.2, 0.15, 0), "cm")),
                 width = 10*1.1, height = 8.5*1.1, bg = 'white')

# cowplot::plot_grid(year_establishment_barplot,
#                    transloc_explore_multipanel_label,
#                    nrow = 2,
#                    rel_heights = c(0.2, 0.8))



#############################################
### SUMMARIES FOR PAPER'S RESULTS SECTION ###
#############################################

transloc_info_doc_path <- here('results', 'text_results', 'translocation_outcomes_info.txt')

cat('##########################\n### SUMMARY OF TRANSLOCATIONS ###\n##########################',
    file = transloc_info_doc_path, sep="\n", append = FALSE)

cat('Translocation info:',
    file = transloc_info_doc_path, sep="\n", append = FALSE)

translocations_inter2_process <- translocations_inter2 %>% 
  mutate(establish = if_else(year_count > 0, 1, 0),
         release_year = clock::get_year(TranslocationDate)) %>% 
  select(RCWid, Sex, establish, Source, nest_years, release_year)

cohort_info_df <- translocations_inter2_process %>% 
  group_by(release_year) %>% 
  summarize(Source = first(Source),
            cohort_size = n()) %>% 
  arrange(release_year)

write.table(cohort_info_df, 
            transloc_info_doc_path,
            append = TRUE,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

cat('\n\n',
    file = transloc_info_doc_path, append = TRUE)





fisher_test_list <- lapply(setNames(nm = c('Sex', 'Source', 'release_year')), function(x, dat) {
  
  dat %>% 
    select(all_of(c(x, 'establish'))) %>% 
    table() %>% 
    fisher.test()
  
}, dat = translocations_inter2_process)

cat("### Fisher's Exact Tests comparing establishment proportions ###",
    file = transloc_info_doc_path, sep="\n", append = TRUE)

for (i in names(fisher_test_list)) {
  cat(i,
      file = transloc_info_doc_path, sep="\n", append = TRUE)
  cat(paste0('Estimate = ', round(fisher_test_list$Sex$estimate, 3) ),
      file = transloc_info_doc_path, sep="\n", append = TRUE)
  cat(paste0('P = ', round(fisher_test_list$Sex$p.value, 3)),
      file = transloc_info_doc_path, sep="\n", append = TRUE)
  cat('\n',
      file = transloc_info_doc_path, sep="", append = TRUE)
}



#establishment percentage by source population
establish_perc_range_sourcepop <- translocations_inter2_process %>% 
  group_by(Source) %>% 
  summarize(establish_percent = round(100*(sum(establish)/n()), 2) ) %>% 
  pull(establish_percent) %>% 
  range()

cat(paste0('Establishment percentage range by source population: ', paste(establish_perc_range_sourcepop, collapse = '-'), '%'),
    file = transloc_info_doc_path, sep="\n", append = TRUE)


#establishment percentage by translocation cohort (birds released in a particular year)
establish_perc_range_releaseyear <- translocations_inter2_process %>% 
  group_by(release_year) %>% 
  summarize(establish_percent = round(100*(sum(establish)/n()), 2) ) %>% 
  pull(establish_percent) %>% 
  range()

cat(paste0('Establishment percentage range by release year: ', paste(establish_perc_range_releaseyear, collapse = '-'), '%'),
    file = transloc_info_doc_path, sep="\n", append = TRUE)

cat('\n',
    file = transloc_info_doc_path, append = TRUE)


#Basic breeding summary:
# -number of translocated birds that nested
# -percentage of translocated birds that nested
breeding_summary_table <- translocations_inter2_process %>% 
  filter(establish == 1) %>% 
  summarize(breed_count = sum(nest_years > 0),
            breed_prop = round(100*(sum(nest_years > 0)/n()), 2)) %>% 
  as.data.frame()

cat('Number of transloc. birds that nested: ', breeding_summary_table$breed_count,
    file = transloc_info_doc_path, sep="\n", append = TRUE)

cat('Percentage of transloc. birds that nested: ', paste0(breeding_summary_table$breed_prop, '%'),
    file = transloc_info_doc_path, sep="\n", append = TRUE)

cat(file = transloc_info_doc_path, sep="\n", append = TRUE)

#for translocated birds that nested, was is the range median values for nesting years
nesting_year_range_median <- translocations_inter2_process %>% 
  filter(nest_years > 0) %>% 
  pull(nest_years) %>% 
  quantile(., prob = c(0, 0.5, 1))

cat('Median and range of breeding years for translocated birds that nested:',
    file = transloc_info_doc_path, sep="\n", append = TRUE)
cat('range: ', paste0(nesting_year_range_median['0%'], '--', nesting_year_range_median['100%']),
    file = transloc_info_doc_path, sep="\n", append = TRUE)
cat('median: ', nesting_year_range_median['50%'],
    file = transloc_info_doc_path, sep="\n", append = TRUE)


parent_offspring_in_pop <- rcw_processed_list$ped_processed %>% 
  pivot_longer(cols = c(fid, mid),
               values_to = 'parent_id',
               names_to = 'parent_type') %>% 
  filter(parent_id != '0') %>% 
  rename(RCWid = id) %>% 
  inner_join(., pop_dat %>% 
               group_by(RCWid) %>% 
               summarize(years_in_pop = n(),
                        .groups = 'drop'), 
             by = 'RCWid',
             relationship = "many-to-many") %>% 
  group_by(parent_id) %>% 
  summarize(offspring_in_pop = n())

translocations_inter2_process1 <- translocations_inter2_process %>%
  left_join(., parent_offspring_in_pop %>% rename(RCWid = parent_id),
            by = 'RCWid', val) %>%
  mutate(offspring_in_pop = if_else(is.na(offspring_in_pop), 0, offspring_in_pop)) %>%
  left_join(., pop_dat %>%
              group_by(RCWid) %>%
              summarize(years_in_pop = n(), .groups = 'drop'),
            by = 'RCWid')

# translocations_inter2_process1 %>% 
#   filter(establish == 1) %>% 
#   summarize(total_years_in_pop = sum(years_in_pop),
#             lower_years_in_pop = min(years_in_pop),
#             upper_years_in_pop = max(years_in_pop),
#             median_years_in_pop = quantile(years_in_pop, prob = 0.5),
#             #total_offspring_in_pop = sum(offspring_in_pop), #this double counts some birds
#             lower_range_offspring_in_pop = min(offspring_in_pop),
#             upper_range_offspring_in_pop = max(offspring_in_pop),
#             median_offspring_in_pop = quantile(offspring_in_pop, prob = 0.5))


#the number of offspring from translocated birds that made it into a census
transloc_offspring_in_census <- pop_dat %>% 
  filter(RCWid %in% get_offspring(rcw_processed_list$ped_processed, translocations_inter2_process1$RCWid)) %>% 
  pull(RCWid) %>% 
  unique() %>% 
  length()

cat('Number of offspring from translocated birds in a census: ', transloc_offspring_in_census,
    file = transloc_info_doc_path, sep="\n", append = TRUE)

cat(file = transloc_info_doc_path, sep="\n", append = TRUE)


#the number of translocated--translocated and translocated--nontranslocated breeding pairs 
breeding_pair_types_table <- transloc_nests %>% 
  mutate(pairing_type = case_when(FemaleID %in% translocations_inter$RCWid & MaleID %in% translocations_inter$RCWid ~ 'transloc_transloc',
                                  (FemaleID %in% translocations_inter$RCWid & !MaleID %in% translocations_inter$RCWid) | (!FemaleID %in% translocations_inter$RCWid & MaleID %in% translocations_inter$RCWid) ~ 'transloc_local',
                                  TRUE ~ 'alt_pairing')) %>% 
  group_by(pairing_type) %>% 
  summarize(pairing_count = n(),
            fldg_count_pairing = sum(fldg_count))

cat("Breakdown of breeding pairs that involve at least one translocated bird:",
    file = transloc_info_doc_path, sep="\n", append = TRUE)
write.table(breeding_pair_types_table, 
            transloc_info_doc_path,
            append = TRUE,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)



###################################################
### TRANSLOCATION SUMMARY TABLE (IN SUPPLEMENT) ###
###################################################
translocation_info_table <- translocations_inter1 %>% 
  rename(`Donor abbr.` = Source) %>% 
  mutate(`Translocation year` = get_year(TranslocationDate),
         transloc_month_num = get_month(TranslocationDate),
         `Donor population` = case_when(`Donor abbr.` == 'ONF' ~ 'Osceola National Forest',
                                        `Donor abbr.` == 'FTS' ~ 'Fort Stewart',
                                        `Donor abbr.` == 'WSF-CITRUS' ~ 'Withlacoochee State Forest---Citrus',
                                        `Donor abbr.` == 'FTB' ~ 'Fort Benning',
                                        `Donor abbr.` == 'ANF' ~ 'Apalachicola National Forest',
                                        `Donor abbr.` == 'CBJTC' ~ 'Camp Blanding Joint Training Center'),
         State = if_else(`Donor abbr.` %in% c('FTB', 'FTS'), 'GA',  'FL')) %>% 
  select(`Donor population`,
         `Donor abbr.`,
         State,
         `Translocation year`, 
         transloc_month_num, 
         Sex) %>% 
  group_by(`Donor population`, `Translocation year`) %>% 
  summarize(`Donor abbr.` = first(`Donor abbr.`),
            State = first(State),
            transloc_month_num = first(transloc_month_num),
            `Total count` = n(),
            Male = sum(Sex == 'M'),
            Female = sum(Sex == 'F'),
            .groups = 'drop') %>% 
  mutate(`Translocation month` = month.abb[transloc_month_num]) %>% 
  select(`Donor population`, 
         `Donor abbr.`,
         State,
         `Translocation year`, 
         `Translocation month`, 
         Male, 
         Female,
         `Total count`) %>% 
  arrange(`Donor population`, `Translocation year`) %>% 
   kbl('latex', booktabs = TRUE,
       align = "c", linesep = "",
       escape = TRUE) %>% 
   collapse_rows(columns = 1:2,
                 valign = "middle",
                 latex_hline = "major") %>% 
   kable_styling(latex_options = c("scale_down", "hold_position"),
                 position = "center")
 
#output table
cat(translocation_info_table, 
    file = here('tables', 'translocation_info_table.txt'), append = FALSE)
 

### INFO ABOUT DONOR POPULATION ABBREVIATIONS ###
#APAFR: Avon Park
#ONF: Osceola National Forest
#FTS: Fort Stewart
#WSF-CITRUS: Withlacoochee State Forest - Citrus
#FTB: Fort Benning
#ANF: Apalachicola National Forest
#CBJTC: Camp Blanding Joint Training Center


#translocations_inter1 %>% 
#  mutate(`Translocation year` = get_year(TranslocationDate),
#         `Translocation month` = get_month(TranslocationDate)) %>% 
#  filter(Source == 'FTS')

#translocations_inter1 %>% 
#  mutate(Year = clock::get_year(TranslocationDate)) %>% 
#  group_by(Year) %>% 
#  summarize(Source = first(Source),
#  )



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# year_hist <- translocations_inter2 %>% 
#   ggplot() +
#   geom_bar(aes(x = year_count, fill = id_factor_source)) +
#   scale_fill_manual(values = setNames(results_list$transloc_info_color[[paste0('source', '_color')]],
#                                       nm = results_list$transloc_info_color$RCWid)) +
#   theme_bw() +
#   theme(legend.position = 'none')


# nest_hist <- translocations_inter2 %>% 
#   ggplot() +
#   geom_bar(aes(x = nest_years, fill = id_factor_source)) +
#   scale_fill_manual(values = setNames(results_list$transloc_info_color[[paste0('source', '_color')]],
#                                       nm = results_list$transloc_info_color$RCWid)) +
#   theme_bw() +
#   theme(legend.position = 'none')
# 

#male_counts <- apply(table(nests$MaleID, nests$Cluster) > 0, 1, sum)
#female_counts <- apply(table(nests$FemaleID, nests$Cluster) > 0, 1, sum)

#male_counts[male_counts > 1]

#female_counts[female_counts > 1]


#https://stackoverflow.com/questions/61193516/how-to-sort-colours-in-r
