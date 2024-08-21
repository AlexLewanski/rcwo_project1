##############################################################################################
### SCRIPT NAME: pop_composition_results_output.R
### PURPOSE: 
### PRODUCTS:
###     translocation.csv: info on translocated individuals
###     census_processed.csv: processed census data
###     nests_processed.csv: info. on all nesting events
###     rcws.csv: info on each individual RCW that has been documented in the population
###     ped_processed.csv: processed version of the population pedigree
###     dummy_parents_info.csv: info on the dummy parents that were added to the pedigree
###     pop_count_scenario_plot_w_cor.png: plot of pop size based on the different approaches
###                                        to processing the census data
###     obs_vs_deduced_barplot.png: plot showing the percentage of observed and deduced
###                                 observations found each year in the scenario 3
###                                 census data (this is the processed version that is
###                                 used in the paper).
##############################################################################################


#####################
### SCRIPT SET-UP ###
#####################

#source(here('scripts', 'ped_functions.R'))
#source(here('scripts', 'rcw_project_custom_functions.R'))

### PACKAGES ###
library(here)
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(egg)
#library(ggridges)
library(magick)
library(patchwork)
library(clock)

here::here()


output_figs <- TRUE


#color palettes
green_colfunc <- colorRampPalette(c("#add6d6", "#329999", "#287a7a", "#1e5b5b"))
red_colfunc <- colorRampPalette(c("#ea9999", "#cc0000", "#a30000", "#7a0000"))

pink_colfunc <- colorRampPalette(c("#f1698a", "#d85070", "#a83e57", "#90354b"))
yellow_colfunc <- colorRampPalette(c("#ffdb19", "#e5c100"))  
mint_colfunc <- colorRampPalette(c("#69e6c6", "#06d6a0", "#038060", "#1e5b5b"))   
purple_colfunc <- colorRampPalette(c("#d8b2d8", "#b266b2", "#800080", "#590059"))
gray_colfunc <- colorRampPalette(c("#b2b2b2", "#808080", "#595959", "#333333"))
blue_colfunc <- colorRampPalette(c("#9fd0e0", "#58adc9", "#118ab2", "#0b607c"))

group_color_vec <- setNames(c("#f0597d", "#ffdb19", "#595959", "#06d6a0", "#118ab2", "#800080", "#ededed"),
                            nm = c("ANF", "CBJTC", "FTB", "FTS", "ONF", "WSF-CITRUS", "non-transloc"))

#RCW illustrations
flying_rcw <- image_flop(image_read(here('figures', 'main_paper', 'bird_illustrations', 'rcw_flying1.png')))
perched_rcw <- image_rotate(image_read(here('figures', 'main_paper', 'bird_illustrations', 'rcw_perched1.png')),
                            degrees = 0)

results_names <- setNames(nm = c("rcws_founder_info",
                                 "rcws_inbr", 
                                "transloc_info_color", 
                                #"ne_f_df", 
                                "rcw_inbr_merge", 
                                "rcws_ancestry_info",
                                "rcws_ancestry_source_by_year",
                                "rcw_partial_founder_summed_processed",
                                "rcw_contr_info_df",
                                "fped_prop_by_group",
                                "inbr_founder_color")
                         )

census_processed <- read.csv(here('data', 'feb2024_databasemarch2023_processed', 'census_processed.csv'))
pop_info_list <- lapply(setNames(nm = paste0('Scenario', 1:4)), function(X, DF) {
  DF[DF[,X,drop = TRUE],!grepl(pattern = 'Scenario', colnames(DF))]
}, DF = census_processed)
pop_dat <- pop_info_list$Scenario3

results_list <- lapply(results_names, function(x) read.csv(here('results', paste0(x, '.csv'))))

lastseendetail <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'LastSeenDetail.xlsx'))
clusterdetail <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'ClusterDetail.xlsx'))
socialclassdetail <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'SocialClassDetail.xlsx'))
nests <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'Nests.xlsx'))
rcws <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'RCWs.xlsx')) 


###################################################################
### FIG 1: POPULATION OVERVIEW (POP SIZE, INBREEDING, ANCESTRY) ###
###################################################################


# pop_dat_with_anc <- pop_dat %>% 
#   rename(id = RCWid) %>% 
#   left_join(., results_list$rcws_ancestry_info,
#             by = 'id',
#             relationship = "many-to-many") %>% 
#   select(id, year, group, anc_prop) %>% 
#   pivot_wider(values_from = anc_prop,
#               names_from = group,
#               values_fill = 0) %>% 
#   group_by(year) %>% 
#   #arrange(!!!rlang::syms(c(anc_vec))) %>% 
#   arrange(`non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`, 
#           .by_group = TRUE) %>% 
#   #group_by(year) %>%
#   mutate(rank_ind = 1:n()) %>% 
#   pivot_longer(cols = c(`non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`), # c(`1`, `2`, `3`),
#                names_to = 'anc_group',
#                values_to = 'anc_val') %>%
#   filter(anc_val > 0) %>% 
#   ungroup()

pop_dat_with_anc <- pop_dat %>% 
  rename(id = RCWid) %>% 
  left_join(., results_list$rcws_ancestry_info,
            by = 'id',
            relationship = "many-to-many")

all_anc_vec <- c("CBJTC","ANF", "FTB", "WSF-CITRUS", "FTS", "ONF", "non-transloc")

pop_dat_with_anc1 <- lapply(split(pop_dat_with_anc, pop_dat_with_anc$year), function(x, all_anc_vec) { 
  anc_vec <- all_anc_vec[all_anc_vec %in% unique(x$group)]
  
  x %>% 
    select(id, year, group, anc_prop) %>% 
    pivot_wider(values_from = anc_prop,
                names_from = group,
                values_fill = 0) %>% 
    #group_by(year) %>% 
    arrange(!!!rlang::syms(c(anc_vec))) %>% 
    #arrange(`non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`, 
    #        .by_group = TRUE) %>% 
    #group_by(year) %>%
    mutate(rank_ind = 1:n()) %>% 
    pivot_longer(cols = all_of(anc_vec), # c(`1`, `2`, `3`),
                 names_to = 'anc_group',
                 values_to = 'anc_val') %>%
    filter(anc_val > 0) %>% 
    ungroup()
}, all_anc_vec = all_anc_vec) %>% 
  bind_rows()

test <- pop_dat_with_anc1 %>% 
  #filter(year %in% c(2015, 2016) ) %>% 
  select(year, id, year, anc_group, anc_val, rank_ind) %>% 
  mutate(year_fac = year) %>% 
  rename(Individuals = rank_ind) %>% 
  ggplot() +
  geom_hline(yintercept = 0.5, 
             color = c(rep('white', 6), '#d8d8d8', rep('white', 9), '#d8d8d8', rep('white', 9), '#d8d8d8', 'white', 'white'),
             linetype = 'dashed') +
  geom_bar(aes(y = anc_val, x = Individuals, fill = anc_group),
           position = "stack", stat = "identity",
           linewidth = 0, width = 1) +
  scale_fill_manual(values = group_color_vec) +
  scale_color_manual(values = group_color_vec) +
  coord_flip(clip = 'off') +
  facet_grid(~year_fac, switch = "x") +
  theme_bw() +
  theme(
    #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
    #strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    
    #strip.background = element_rect(fill = NA, color = "white"),
    strip.text = element_text(size = 7.5),
    panel.spacing.x = unit(0,"cm"),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
    axis.line.x = element_blank(),
    axis.ticks = element_line(color = '#808080', linewidth = 0.7),
    legend.title = element_blank(),
    legend.position = 'none',
    plot.margin = unit(c(-0.01, 0.1, -0.01, 1), "cm"),
    axis.title.y = element_text(size = 10)
    #axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
  ) +
  scale_x_continuous(expand = c(0.01, 0)#,
                     #breaks = c(50, 100, 150)
  ) +
  scale_y_continuous(lim = c(-0.1, 1.1), expand = c(0, 0)) +
  ylab('Individuals') +
  geom_text(data = results_list$transloc_info_color %>% 
              group_by(year) %>%
              summarize(Source = first(Source),
                        count = n(),
                        .groups = 'drop') %>% 
              mutate(year_fac = year) %>% 
              left_join(.,
                        pop_dat %>% 
                          group_by(year) %>% 
                          summarize(pop_count = n()),
                        by = 'year'),
            aes(x = pop_count + 13, y = 0.5, label = count, color = Source),
            fontface = 'bold',
            size = 4.75)


lastseendetail_subset <- lastseendetail %>% 
  filter(Detected == 'Yes' & SurveyType == 'Census' & clock::get_month(lastseendetail$RecordDate) %in% c(6, 7, 8))

lastseendetail_subset$most_recent_date <- lubridate::mdy_hms(NA)
lastseendetail_subset$ClusterID <- NA
for (i in seq_len(nrow(lastseendetail_subset))) {
  
  focal_date <- lastseendetail_subset$RecordDate[i]
  if (is.na(focal_date)) next
  if (!lastseendetail_subset$RCWid[i] %in% clusterdetail$RCWid) stop('STOP!')
  dates_vec <- clusterdetail$RecordDate[clusterdetail$RCWid == lastseendetail_subset$RCWid[i]]
  
  if (any(dates_vec <= focal_date)) {
    lastseendetail_subset$most_recent_date[i] <- max(dates_vec[dates_vec <= focal_date])[1]
  } else {
    lastseendetail_subset$most_recent_date[i] <- min(dates_vec)[1]
  }
  
  lastseendetail_subset$ClusterID[i] <- clusterdetail$ClusterID[clusterdetail$RCWid == lastseendetail_subset$RCWid[i] & clusterdetail$RecordDate == lastseendetail_subset$most_recent_date[i]][1]
}



lastseendetail_subset2 <- lastseendetail_subset %>% 
  #filter(RecordDate %in% censusdata$RecordDate[censusdata$SurveySeason == 'Post-season']) %>% 
  left_join(rcws[,c('RCWid', 'Sex', 'MinAge', 'AgeKnown')], by = 'RCWid') %>% 
  mutate(census_year = clock::get_year(RecordDate)) %>% 
  filter(MinAge < census_year | is.na(MinAge)) %>% 
  mutate(ClusterID_Year = paste0(ClusterID, "_", census_year)) %>% 
  group_by(ClusterID_Year) %>% 
  mutate(both_sex = if_else(all(c('M', 'F') %in% Sex), 1, 0)) %>% 
  ungroup() %>% 
  mutate(transloc = if_else(RCWid %in% results_list$transloc_info_color$RCWid, 'yes', 'no') )


pbg_summary <- lastseendetail_subset2 %>% 
  group_by(ClusterID_Year) %>%
  summarize(ClusterID = first(ClusterID),
            year = first(census_year),
            both_sex = if_else(all(c('M', 'F') %in% Sex), 1, 0),
            `Group\ncomposition` = case_when(all(transloc == 'yes') ~ 'only transloc.',
                                             all(transloc == 'no') ~ 'no transloc.',
                                             TRUE ~ 'combined'),
            .groups = 'drop') %>%
  filter(both_sex == 1) #%>%
##filter(!grepl('^RR', ClusterID)) %>% 
#group_by(year) %>% 
#summarize(pbg_count = n(),
#          .groups = 'drop')



bottom_bar_pbg <- pbg_summary %>% 
  mutate(year_fac = year,
         `Group\ncomposition` = factor(`Group\ncomposition`, 
                                       levels = rev(c('no transloc.', 'combined', 'only transloc.')))
  ) %>% 
  ggplot() +
  geom_vline(xintercept = 'a',
             color = c(rep('white', 6), '#d8d8d8', rep('white', 9), '#d8d8d8', rep('white', 9), '#d8d8d8', 'white', 'white'),
             linetype = 'dashed') +
  #geom_bar(aes(x = 'a', y = count, fill = anc_prop),
  #         position = "stack", stat = "identity", 
  #         linewidth = 0, width = 1) +
  geom_bar(aes(x = 'a', fill = `Group\ncomposition`),
           position = "stack", stat = 'count', width = 0.98) +
  scale_fill_manual(values = c("#6a6a6a", '#bdbdbd', "#ededed")) + #c('#333333', "#999999", "#ededed")) +
  #facet_grid(~year_fac, switch = "x") +
  facet_grid(~year_fac) +
  theme_bw() +
  ylab('Potential breeding groups') +
  theme(
    legend.title = element_text(size = 9.75),
    #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
    #strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    #strip.background = element_rect(fill = NA, color = "white"),
    strip.text = element_text(size = 7.5),
    panel.spacing.x = unit(0,"cm"),
    #axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
    #axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.5),
    axis.line.x = element_blank(),
    axis.ticks = element_line(color = '#808080', linewidth = 0.7),
    #legend.title = element_blank(),
    #legend.position = 'none',
    plot.margin = unit(c(-0.05, 0.1, -0.3, 1), "cm"),
    axis.title.y = element_text(size = 10, margin = margin(0, 2.5, 0, 0, unit = "mm"))
    #axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
  ) +
  #scale_y_reverse(expand = expansion(mult = c(0, 0), 
  #                                     add = c(1, 0))) +
  theme(plot.margin=unit(c(-0.01, 1, -0.01,1), "cm"),
        #axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()#,
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_x_continuous(expand = c(0, 0)#, 
  #                   #limits = c(1993.5, 2022.5)
  #                   ) +
  coord_cartesian(clip = 'off')


rcw_inbr_merge_processed <- results_list$rcw_inbr_merge %>% 
  group_by(year) %>% 
  summarize(perc_10 = quantile(f_ped, prob = 0.1),
            mean = mean(f_ped),
            perc_90 = quantile(f_ped, prob = 0.9),
            .groups = 'drop') %>% 
  arrange(year) %>% 
  mutate(change_color = case_when(mean > lag(mean) ~ "#f0597d",
                                  mean < lag(mean) ~ "#118ab2",
                                  TRUE ~ '#4c4c4c'))

fped_plot_facet <- results_list$rcw_inbr_merge %>% 
  ggplot() +
  geom_vline(xintercept = 'a',
             color = c(rep('white', 6), '#d8d8d8', rep('white', 9), '#d8d8d8', rep('white', 9), '#d8d8d8', 'white', 'white'),
             linetype = 'dashed') +
  geom_point(aes(x = 'a', y = f_ped),
             shape = 21, colour = '#b7b7b7', fill = '#e3e3e3', 
             size = 1, alpha = 0.8,
             position = position_jitter(w = 0.16, h = 0)) +
  geom_segment(data = rcw_inbr_merge_processed,
               aes(x = 'a', xend = 'a', 
                   y = perc_10, yend = perc_90),
               colour = '#4c4c4c', linewidth = 1.15) +
  geom_point(data = rcw_inbr_merge_processed,
             aes(x = 'a', y = mean),
             shape = 21, color = rcw_inbr_merge_processed$change_color, 
             fill = 'white', 
             size = 2.25, stroke = 1.5) +
  theme_bw() +
  facet_wrap(~year,
             nrow = 1) +
  theme(
    #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
    #strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    #strip.background = element_rect(fill = NA, color = "white"),
    strip.text = element_text(size = 7.5),
    panel.spacing.x = unit(0,"cm"),
    #axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
    #axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.5),
    axis.line.x = element_blank(),
    axis.ticks = element_line(color = '#808080', linewidth = 0.7),
    #legend.title = element_blank(),
    #legend.position = 'none',
    plot.margin = unit(c(-0.05, 0.1, -0.3, 1), "cm"),
    axis.title.y = element_text(size = 10)
    #axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
  ) +
  #scale_y_reverse(expand = expansion(mult = c(0, 0), 
  #                                   add = c(1, 0))) +
  theme(plot.margin=unit(c(-0.01, 1, -0.01,1), "cm"),
        #axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()#,
  ) +
  #scale_y_continuous(lim = c(-0.1, 1.1), expand = c(0, 0)) +
  #scale_x_continuous(expand = c(0, 0)#, 
  #                   #limits = c(1993.5, 2022.5)
  #                   ) +
  ylab(expression(italic(F)["P"])) +
  coord_cartesian(clip = 'off')

transloc_ancestry_facet <- results_list$rcws_ancestry_source_by_year %>%
  mutate(group = factor(group, 
                        levels = rev(c('non-transloc', 'ANF', 'CBJTC', 'FTB', 'FTS', 'ONF', 'WSF-CITRUS')))) %>% 
  ggplot() +
  geom_hline(yintercept = 0.5, 
             color = c(rep('white', 6), '#d8d8d8', rep('white', 9), '#d8d8d8', rep('white', 9), '#d8d8d8', 'white', 'white'),
             linetype = 'dashed') +
  geom_bar(aes(x = 'a', y = ancestry_prop, fill = group), 
           color = 'white', linewidth = 0.15, width = 0.72,
           position = "stack", stat = "identity") + 
  #scale_fill_manual(values = c('#ededed', '#4c4c4c')) +
  scale_fill_manual(values = group_color_vec) +
  facet_grid(~year) +
  theme_bw() +
  theme(
    #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
    #strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    #strip.background = element_rect(fill = NA, color = "white"),
    strip.text = element_text(size = 7.5),
    panel.spacing.x = unit(0,"cm"),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
    axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.5),
    #axis.line.x = element_blank(),
    axis.ticks = element_line(color = '#808080', linewidth = 0.7),
    #legend.title = element_blank(),
    legend.position = 'none',
    plot.margin = unit(c(-0.05, 0.1, -0.1, 1), "cm"),
    axis.title.y = element_text(size = 10, margin = margin(0, 1.8, 0, 0, unit = "mm"))
    #axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0.0", "0.5", "1.0"),
                     expand = expansion(c(0, 0))) +
  xlab('Year') +
  ylab("Ancestry") +
  theme(legend.key.size = unit(0.4, 'cm')) +
  coord_cartesian(clip = 'off')

### custom legends ###
fped_plot_legend <- ggdraw(cowplot::get_legend(
  data.frame(x = 1,
             y = 3:1,
             label = factor(c('increase', 'no change', 'decrease'), 
                            levels = c('increase', 'no change', 'decrease')) ) %>% 
    ggplot() +
    geom_point(aes(x = x,
                   y = y,
                   color = label),
               shape = 21,
               fill = 'white',
               size = 2.1, stroke = 1.5) +
    theme_bw() +
    theme(legend.key = element_blank(),
          legend.title = element_blank(),
          legend.spacing.y = unit(-1.5, 'mm'),
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent", color = "transparent")) +
    scale_color_manual(substitute(paste('Change in mean(', italic(F)["p"], ")")),
                       values = c("#f0597d", '#4c4c4c', "#118ab2")) +
    guides(color = guide_legend(byrow = TRUE)) #required to adjust vertical spacing of legend
  
)) +
  theme(plot.margin=unit(c(0,-5, 0, -5), "cm"),
        panel.background = element_rect(fill = "transparent", color = "transparent"))


founder_prop_legend <- results_list$rcws_founder_info %>%
  group_by(group) %>% 
  summarize(count = n(),
            .groups = 'drop') %>%
  mutate(group = factor(group, 
                        levels = rev(c('CBJTC', 'ANF', 'FTB', 'WSF-CITRUS', 'FTS', 'ONF', 'non-transloc')))) %>% 
  mutate(x = 0.95) %>%
  arrange(desc(group)) %>% 
  mutate(cumsum_count = cumsum(count) - (0.5*count)) %>% 
  mutate(cumsum_count_adjust = cumsum_count + c(-13, -5, 0, 10, 17, 24, 0)) %>% 
  mutate(update_label = paste0(group, ' (', count, ')')) %>% 
  ggplot(aes(x = x, y = count, fill = group, label = count)) +
  geom_segment(aes(x = 1.47,
                   xend = 2,
                   y = cumsum_count, 
                   yend = cumsum_count_adjust,
                   color = group),
               linewidth = 1.1) +
  geom_bar(position = 'stack', stat = "identity", color = 'white', linewidth = 0.15) +
  scale_color_manual(values = group_color_vec) +
  scale_fill_manual(values = group_color_vec) +
  geom_text(aes(x = 2.1, y = cumsum_count_adjust, label = update_label),
            size = 2.6,
            hjust = 0) +
  theme_void() +
  theme(legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  xlim(-0.32, 4) +
  geom_text(aes(x = 0.5, y = 220, label = "Founder\nproportions"),
            size = 3.5,
            hjust = 0,
            lineheight = .9)


rcw_pop_overview_facet <- egg::ggarrange(fped_plot_facet +
                                           theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm")),
                                         bottom_bar_pbg +
                                           theme(plot.margin = unit(c(0.2, 0, 0, 0.1), "cm")),
                                         test +
                                           theme(plot.margin = unit(c(0.2, 0, 0, 0.1), "cm")),
                                         transloc_ancestry_facet +
                                           theme(plot.margin = unit(c(0.2, 0, 0, 0.1), "cm")),
                                         nrow = 4,
                                         labels = c('A', 'B', 'C', 'D'),
                                         label.args = list(gp = grid::gpar(font = 2, cex = 1.2, hjust = -0.1)),
                                         heights = c(30, 50, 50, 22))


test_multipan <- ggdraw(rcw_pop_overview_facet) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.1), "cm")) +
  annotate("text", x = c(0.23, 0.521, 0.814), y = -0.011, label = c("2000", "2010", "2020"),
           size = 3.25) +
  annotate("text", x = c(0.23, 0.521, 0.814), y = 1.015, label = c("2000", "2010", "2020"),
           size = 3.25)

rcw_pop_overview_rcw <- test_multipan +
  cowplot::draw_image(flying_rcw,
                      scale = 0.12,
                      x = 0.4, y = 0.46) +
  patchwork::inset_element(founder_prop_legend,
                           left = 0.918 - 0.033, bottom = 0.18, right = 1.015 - 0.033, top = 0.42
                           #left = 0.849 - 0.0358, bottom = 0.5, right = 0.98 - 0.0358, top = 0.70
                           #left = 0.849, bottom = -0.01, right = 0.972, top = 0.20
  ) +
  patchwork::inset_element(fped_plot_legend,
                           #left = 0, bottom = 0.9, right = 0.38, top = 0.95,
                           #left = 0.891 - 0.0358, bottom = 0.81, right = 0.94 - 0.0358, top = 0.88
                           #left = 0.891 - 0.032, bottom = 0.85, right = 0.94 - 0.032, top = 0.86
                           left = 0.965 - 0.032, bottom = 0.83, right = 0.99 - 0.032, top = 0.86
  )


#library(grid)
#library(gridExtra)

shared_axis_title <- textGrob("Year", 
                              gp = gpar(fontface = "plain", col = "black", fontsize=15),
                              x = unit(0.445, "npc"),
                              y = unit(0.8, "npc"))

#add to plot

rcw_pop_overview_rcw_title <- grid.arrange(arrangeGrob(ggdraw(rcw_pop_overview_rcw), 
                                                       bottom = shared_axis_title))

cowplot::ggsave2(filename = here('figures', 'main_paper', 'rcw_pop_overview_fig_TEST.png'),
                 #plot = ggdraw(rcw_pop_overview_rcw) + theme(plot.margin = unit(c(-0.3, 0.3, 0.25, -0.6), "cm")),
                 plot = ggdraw(rcw_pop_overview_rcw_title) +
                   theme(plot.margin = unit(c(-0.22, -0.2, -0.2, 0), "cm")), #ggdraw(rcw_pop_overview_rcw) +
                 #theme(plot.margin = unit(c(0, -0.9, 0, 0), "cm")),
                 #width = 10*1, height = 8.1*1, 
                 width = 10*1.2, height = 5.5*1.2, 
                 bg = 'white')



###################################################################
### FIG 2 ###
###################################################################

fped_contr_transloc_ids <- results_list$rcw_partial_founder_summed_processed %>% 
  filter(year == 2022) %>% 
  filter(group != 'non-transloc') %>% 
  pull(founder_id)


fig2_list <- list()

for (i in c('source', 'sex')[1]) { #currently only including coloration by source
  
  fig2_list[[i]][['contr_plot']] <- results_list$rcw_contr_info_df %>%
    group_by(id) %>%
    filter(any(contr != 0)) %>% #filter out the ids with all 0s
    filter(contr > 0 | (contr == 0 & year > max(year[contr > 0])  )) %>%
    left_join(., results_list$rcws_founder_info,
              by = 'id') %>%
    ggplot() +
    geom_line(data = . %>% 
                filter(group == 'non-transloc'),
              aes(x = year, y = contr, group = id, linewidth = `group`),
              color = '#e3e3e3', linewidth = 0.7) +
    geom_point(data = . %>%
                 group_by(id) %>%
                 filter(year == min(year)) %>%
                 slice_head(n = 1)%>% 
                 filter(group == 'non-transloc'),
               aes(x = year, y = contr),
               color = '#e3e3e3', size = 0.7) +
    geom_line(data = . %>% 
                filter(group != 'non-transloc'),
              aes(x = year, y = contr, group = id, 
                  linewidth = group, color = id),
              linewidth = 1.1) +
    geom_point(data = . %>%
                 group_by(id) %>%
                 filter(year == min(year)) %>%
                 slice_head(n = 1)%>% 
                 filter(group != 'non-transloc'),
               aes(x = year, y = contr, color = id),
               size = 1.5) +
    scale_color_manual(values = setNames(results_list$transloc_info_color[[paste0(i, '_color')]],
                                         nm = results_list$transloc_info_color$RCWid)) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
          panel.border = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
          axis.ticks = element_line(color = '#808080', linewidth = 0.7),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.32,1, 0, 0.5), "cm"),
          legend.position = 'none',
          axis.title.y = element_text(size = 12, margin = margin(0, 12, 0, 0))
    ) +
    scale_x_continuous(expand = c(0, 0),
                       #expand = c(0.01, 0), 
                       #limits = c(1993.5, 2022.5),
                       limits = c(1993.5, 2023.5),
                       position = "top") +
    xlab('Year') +
    ylab("Expected genetic contribution") +
    geom_text(data = results_list$rcw_contr_info_df %>%
                filter(year == 2022 & id %in% fped_contr_transloc_ids) %>%
                #filter(year == 2022 & id %in% "ZG-GWO") %>%
                mutate(y_pos = case_when(id == "OWO-AZ" ~ contr + 0.0005,
                                         id == "ZG-GWO" ~ contr - 0.003,
                                         id == 'WZ-YKA' ~ contr + 0.001,
                                         TRUE ~ contr)),
              aes(x = 2022.6, y = y_pos, label = id),
              size = 2.55,
              hjust = 0) +
    geom_segment(data = results_list$rcw_contr_info_df %>%
                   filter(year == 2022 & id %in% fped_contr_transloc_ids) %>%
                   #filter(id %in% c("OWO-AZ", "ZG-GWO", 'WZ-YKA')) %>%
                   mutate(yend = case_when(id == "OWO-AZ" ~ contr + 0.0005,
                                           id == "ZG-GWO" ~ contr - 0.003,
                                           id == 'WZ-YKA' ~ contr + 0.001,
                                           TRUE ~ contr)),
                 
                 aes(x = 2022.05, xend = 2022.5,
                     y = contr, yend = yend)) +
    coord_cartesian(clip = 'off')
  # ggrepel::geom_text_repel(data = results_list$rcw_contr_info_df %>%
  #             filter(year == 2022 & id %in% fped_contr_transloc_ids),
  #           aes(x = 2022.1, y = contr, label = id),
  #           size = 2.5,
  #           hjust = 0,
  #           min.segment.length = 0.3,
  #           box.padding = 0.4,
  #           seed = 4534,
  #           direction = 'y',
  #           nudge_y = setNames(c(      -0.001,   0.001,    -0.002,   -0.003,   0.0001), 
  #                              nm = c('ZG-GWO', 'OHA-ZK', 'OWO-AZ', 'WZ-YKA', 'OAG-HZ'))
  # )
  
  
  fig2_list[[i]][['fped_prop_by_indiv_founder_plot']] <- results_list$rcw_partial_founder_summed_processed %>% 
    ggplot() +
    geom_bar(aes(x = year, y = pfound_fsum_scale, fill = founder_id,
                 group = -pfound_fsum_scale),
             color = 'white', linewidth = 0.2,
             position = "stack", stat = "identity", width = 0.98, alpha = 0.85) +
    #geom_point(data = . %>% filter(group == 'transloc'),
    geom_point(data = . %>% filter(group != 'non-transloc'),
               aes(x = year, y = val_cumsum, color = founder_id),
               size = 1.5) +
    geom_line(#data = . %>% filter(group == 'transloc'),
      data = . %>% filter(group != 'non-transloc'),
      aes(x = year, y = val_cumsum, 
          group = founder_id, color = founder_id),
      linewidth = 0.7, linetype = 'solid') +
    scale_fill_manual(values = setNames(results_list$inbr_founder_color[[paste0(i, '_color')]],
                                        nm = results_list$inbr_founder_color$founder_id)) +
    scale_color_manual(values = setNames(results_list$inbr_founder_color[[paste0(i, '_color')]],
                                         nm = results_list$inbr_founder_color$founder_id)) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
      panel.border = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
      axis.ticks = element_line(color = '#808080', linewidth = 0.7),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = unit(c(0.25,1,-0.05,0.5), "cm"),
      legend.position = 'none',
      axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 4.5, b = 0, l = 0))) +
    ylab(substitute(paste("Total ", italic(F)["p"]))) +
    #ylab(substitute(italic(F)["P"])) +
    scale_x_continuous(expand = c(0, 0),
                       #expand = c(0.01, 0), 
                       limits = c(1993.5, 2023.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(clip = 'off') +
    geom_text(data = results_list$rcw_partial_founder_summed_processed %>% 
                filter(founder_id %in% fped_contr_transloc_ids & year == 2022),
              aes(x = 2022.6, y = val_cumsum, label = founder_id),
              size = 2.5,
              hjust = 0)
  
  fig2_list[[i]][['fped_prop_by_group_plot']] <- results_list$fped_prop_by_group %>%
    mutate(group = case_when(group == "non-transloc" ~ "non-transloc",
                             TRUE ~ 'transloc')) %>% 
    group_by(year, group) %>% 
    summarize(fped_group_prop = sum(fped_group_prop),
              .groups = 'drop') %>% 
    ggplot() +
    geom_bar(aes(x = year, y = fped_group_prop, fill = group), 
             color = 'white', linewidth = 0.35, width = 0.70,
             position = "stack", stat = "identity") + 
    scale_fill_manual(values = c('#ededed', '#4c4c4c')) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
          panel.border = element_blank(),
          #axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
          axis.line.x = element_blank(),
          axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
          plot.margin = unit(c(0.25,1,0,0.5), "cm"),
          axis.title.x = element_text(size = 18),
          legend.position = 'none',
          axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 14, b = 0, l = 0))
    ) +
    # scale_y_continuous(breaks = c(0, 0.5, 1), 
    #                    labels = c("0.0", "0.5", "1.0")) +
    # scale_x_continuous(expand = c(0.01, 0),
    #                    limits = c(1993.5, 2023.5)
    #                    #limits = c(1993.5, 2022.5)
    #                    ) +
    # xlab('Year') +
    # ylab(substitute(paste(italic(F)["p"], ' prop.')))
    scale_y_continuous(expand = c(0, 0),
                       breaks = c(0, 0.5, 1),
                       limits = c(-0.05, 1.02),
                       labels = c("0.0", "0.5", "1.0")) +
    scale_x_continuous(expand = c(0, 0),
                       #expand = c(0.01, 0),
                       limits = c(1993.5, 2023.5)
                       #limits = c(1993.5, 2022.5)
    ) +
    xlab('Year') +
    ylab(substitute(paste(italic(F)["P"], ' prop.'))) +
    geom_segment(aes(x = 1993.5, xend = 2022.5,
                     y = -0.05, yend = -0.05),
                 color = '#4c4c4c', linewidth = 1.1) +
    coord_cartesian(clip = 'off')
  
  fig2_list[[i]][['panel1']] <- egg::ggarrange(fig2_list[[i]][['contr_plot']],
                                               fig2_list[[i]][['fped_prop_by_indiv_founder_plot']],
                                               fig2_list[[i]][['fped_prop_by_group_plot']],
                                               nrow = 3,
                                               labels = c('A', 'B', 'C'),
                                               label.args = list(gp = grid::gpar(font = 2, cex = 1.3)),
                                               heights = c(65, 65, 12))
  
  fig2_list[[i]][['panel2']] <- cowplot::ggdraw(fig2_list[[i]][['panel1']]) +
    cowplot::draw_image(perched_rcw,
                        scale = 0.2,
                        #x = -0.3695, 
                        x = -0.405, 
                        y = -.13) +
    theme(plot.margin = unit(c(0.1, -0.7, 0, 0.1), "cm"))
  
  if (isTRUE(output_figs)) {
    cowplot::ggsave2(filename =  here('figures', 'main_paper', paste0('indiv_fped_gencontr_panel_', i, '.png')),
                     plot = fig2_list[[i]][['panel2']],
                     width = 11*1.1, height = 5.9*1.1, bg = 'white')
  }
}



###################################################################
### FIG 3 ###
###################################################################

#data components:
# population size each year

# contributions of each year's cohort

# file_names <- list.files(here('data', 'feb2024_databasemarch2023_processed'))
# rcw_processed_list <- lapply(setNames(file_names, nm = gsub(pattern = '\\.csv', '', file_names)), function(x) {
#   read.csv(here('data', 'feb2024_databasemarch2023_processed', x))
# })
  
# time between cohort and current year
# transloc_details <- results_list$transloc_info_color %>% 
#   rename(id = RCWid) %>% 
#   mutate(first_contr_year = year + 1,
#          time_lag = 2022 - first_contr_year) %>% 
#   left_join(., 
#             results_list$rcw_contr_info_df %>% 
#               filter(year == 2022) %>% 
#               rename(contr_year = year), 
#             by = 'id') %>% 
#   #group_by(first_contr_year) %>% 
#   # summarize(time_lag = first(time_lag),
#   #           cohort_contr = sum(contr),
#   #           transloc_count = n()) %>% 
#   left_join(., pop_info_list$Scenario4 %>% 
#               group_by(year) %>% 
#               summarize(pop_size = n()) %>% 
#               rename(first_contr_year = year),
#             by = 'first_contr_year')


contr_info_processed_year <- results_list$rcw_contr_info_df %>% 
  left_join(., results_list$transloc_info_color %>% 
              mutate(first_contr_year = year + 1) %>% 
              rename(id = RCWid) %>% 
              select(id, first_contr_year, Source),
            by = 'id') %>% 
  filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% 
  filter(year >= first_contr_year) %>% 
  #mutate(years_after_contr = first_contr_year - year) %>% 
  group_by(year, first_contr_year) %>% 
  summarize(Source = first(Source),
            cohort_contr = sum(contr),
            `Translocation\ncount` = length(unique(id)),
            .groups = 'drop')

contr_info_processed_source <- results_list$rcw_contr_info_df %>% 
  left_join(., results_list$transloc_info_color %>% 
              mutate(first_contr_year = year + 1) %>% 
              rename(id = RCWid) %>% 
              select(id, first_contr_year, Source),
            by = 'id') %>% 
  filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% 
  filter(year >= first_contr_year) %>% 
  #mutate(years_after_contr = first_contr_year - year) %>% 
  group_by(Source) %>% 
  mutate(year_cohort_count = length(unique(first_contr_year))) %>% 
  ungroup() %>% 
  group_by(year, Source) %>% 
  summarize(cohort_contr = sum(contr),
            year_cohort_count = first(year_cohort_count),
            .groups = 'drop')


# results_list$rcw_contr_info_df %>% 
#   left_join(., results_list$transloc_info_color %>% 
#               mutate(first_contr_year = year + 1) %>% 
#               rename(id = RCWid) %>% 
#               select(id, first_contr_year, Source),
#             by = 'id') %>% 
#   filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% 
#   filter(year >= first_contr_year) %>% 
#   filter(year == max(year)) %>% 
#   arrange(contr)
# 

# contr_info_processed_source %>% 
#   filter(Source == 'ANF') %>% 
#   ggplot() +
#   geom_line(aes(x = year, y = cohort_contr)) +
#   geom_line(data = contr_info_processed_year %>% 
#               filter(Source == 'ANF'),
#             aes(x = year,
#                 y = cohort_contr,
#                 color = as.character(first_contr_year)
#                 )
#             ) +
#   theme_bw()

# contr_info_processed_year %>% 
#   filter(Source == 'ANF') %>% 
#   ggplot() + 
#   geom_area(aes(x = year, 
#                 y = cohort_contr, 
#                 fill = as.character(first_contr_year))) +
#   geom_line(data = contr_info_processed_source %>% 
#               filter(Source == 'ANF'),
#             mapping = aes(x = year, y = cohort_contr)
#             )

processed_data <- contr_info_processed_year %>% 
  mutate(`Translocation\nyear` = first_contr_year - 1) %>% 
  group_by(first_contr_year) %>% 
  mutate(zeroed_year = year - min(year))

background_contr_df <- lapply(setNames(nm = unique(processed_data$Source)), function(x, dat) {
  return(
    dat %>% 
      mutate(Source1 = Source,
             Source = x)
  )
}, dat = processed_data) %>% 
  bind_rows()



cohort_contr_multipanel <- ggplot()  +
  geom_hline(yintercept = 0, linewidth = 1.2) +
  geom_line(data = background_contr_df,
            aes(x = year, y = cohort_contr,
                group = `Translocation\nyear`),
            linewidth = 0.4,
            alpha = 0.3,
            color = '#8c8c8c'
  ) +
  geom_point(data = background_contr_df,
             aes(x = year, y = cohort_contr,
                 group = `Translocation\nyear`),
             size = 0.9, 
             alpha = 0.3,
             color = '#8c8c8c'
  ) +
  geom_line(data = processed_data,
            aes(x = year, y = cohort_contr,
                group = `Translocation\nyear`,
                color = `Translocation\ncount`),
            linewidth = 1.85) +
  geom_point(data = processed_data,
             aes(x = year, y = cohort_contr,
                 group = `Translocation\nyear`,
                 color = `Translocation\ncount`),
             size = 2.25, alpha = 1) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),
    strip.background.x = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = 'dashed', 
                                      color = '#d8d8d8',
                                      linewidth = 0.3),
    panel.border = element_blank(),
    axis.line = element_line(color = '#4c4c4c', 
                             linewidth = 1.5),
    axis.line.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank()
  ) +
  scale_colour_gradientn(colors = c("#feb2b2", "#fe7f7f", "#fe0000", "#b10000", "#7f0000", "#650000"),
                         breaks = c(1, 4, 7, 10), 
                         labels = c("1", "4", "7", "10")) +
  xlab("Year") +
  ylab('Expected genetic contribution') +
  facet_wrap(~Source, 
             #scale = 'free_x', 
             ncol = 1) +
  geom_text(data = processed_data,
            aes(label = Source), x = 1998.5, y = Inf, 
            hjust = 0, vjust = 1.5, check_overlap = TRUE,
            position = position_nudge(x = 100), size = 5,
            fontface = 'bold') +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(sec.axis = dup_axis())

cowplot::ggsave2(filename = here('figures', 'main_paper', 'cohort_contr_multipanel.png'),
                 plot = cohort_contr_multipanel,
                 width = 5*1, height = 7.5*1, bg = 'white')



###################################################################
### FIG 4 ###
###################################################################

if (!exists('pop_dat_with_anc1')) {
  pop_dat_with_anc <- pop_dat %>% 
    rename(id = RCWid) %>% 
    left_join(., results_list$rcws_ancestry_info,
              by = 'id',
              relationship = "many-to-many")
  
  all_anc_vec <- c("CBJTC","ANF", "FTB", "WSF-CITRUS", "FTS", "ONF", "non-transloc")
  
  pop_dat_with_anc1 <- lapply(split(pop_dat_with_anc, pop_dat_with_anc$year), function(x, all_anc_vec) { 
    anc_vec <- all_anc_vec[all_anc_vec %in% unique(x$group)]
    
    x %>% 
      select(id, year, group, anc_prop) %>% 
      pivot_wider(values_from = anc_prop,
                  names_from = group,
                  values_fill = 0) %>% 
      #group_by(year) %>% 
      arrange(!!!rlang::syms(c(anc_vec))) %>% 
      #arrange(`non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`, 
      #        .by_group = TRUE) %>% 
      #group_by(year) %>%
      mutate(rank_ind = 1:n()) %>% 
      pivot_longer(cols = all_of(anc_vec), # c(`1`, `2`, `3`),
                   names_to = 'anc_group',
                   values_to = 'anc_val') %>%
      filter(anc_val > 0) %>% 
      ungroup()
  }, all_anc_vec = all_anc_vec) %>% 
    bind_rows()
}


anc_mixing_df <- data.frame(year = sort(unique(pop_dat$year)),
                            md = NA)

for (YEAR in anc_mixing_df$year) {
  #for (YEAR in 1999) {  
  focal_year_indivs <- pop_dat %>% 
    filter(year == YEAR) %>% 
    pull(RCWid)
  
  anc_info_focal_year <- results_list$rcws_ancestry_info %>% 
    filter(id %in% focal_year_indivs)
  
  total_anc_info <- anc_info_focal_year %>% 
    group_by(group) %>% 
    summarize(anc_sum = sum(anc_prop), .groups = 'drop') %>% 
    mutate(anc_prop = anc_sum/sum(anc_sum))
  
  if (nrow(total_anc_info) == 1) {
    anc_mixing_df$md[anc_mixing_df$year == YEAR] <- NaN
    next
  }
  
  var_vec <- setNames(rep(NA, length(unique(total_anc_info$group))), nm = unique(total_anc_info$group))
  for (i in names(var_vec)) {
    
    anc_subset <- anc_info_focal_year %>% 
      filter(group == i)
    
    var_vec[i] <- var(c(anc_subset$anc_prop, rep(0, length(unique(anc_info_focal_year$id[!anc_info_focal_year$id %in% anc_subset$id])))))
  }
  
  var_vec_remove <- var_vec[names(var_vec) != 'non-transloc']
  total_anc_info_remove <- total_anc_info %>% 
    filter(group != 'non-transloc')
  
  anc_mixing_df$md[anc_mixing_df$year == YEAR] <- 1 - (sum(var_vec_remove)/sum(sapply(total_anc_info_remove$anc_prop, function(x) x*(1 - x))))
  
  #anc_mixing_df$md[anc_mixing_df$year == YEAR] <- 1 - (sum(var_vec)/sum(sapply(total_anc_info$anc_prop, function(x) x*(1 - x))))
  
}

anc_by_year_quantile_summary <- pop_dat %>% 
  rename(id = RCWid) %>% 
  left_join(., results_list$rcws_ancestry_info,
            by = 'id',
            relationship = "many-to-many") %>% 
  group_by(year, id) %>% 
  summarize(ind_ancestry_count = n(), .groups = 'drop') %>% 
  group_by(year) %>% 
  summarize(perc_5 = quantile(ind_ancestry_count, prob = 0.05),
            perc_10 = quantile(ind_ancestry_count, prob = 0.10),
            perc_25 = quantile(ind_ancestry_count, prob = 0.25),
            perc_50 = quantile(ind_ancestry_count, prob = 0.50),
            perc_75 = quantile(ind_ancestry_count, prob = 0.75),
            perc_90 = quantile(ind_ancestry_count, prob = 0.90),
            perc_95 = quantile(ind_ancestry_count, prob = 0.95),
            .groups = 'drop')


anc_count_perc_plot <- anc_by_year_quantile_summary %>% 
  mutate(year_fac = factor(year)) %>% 
  ggplot() +
  geom_point(data = pop_dat %>%
               rename(id = RCWid) %>%
               left_join(., results_list$rcws_ancestry_info,
                         by = 'id',
                         relationship = "many-to-many") %>%
               group_by(year, id) %>%
               summarize(ind_ancestry_count = n(), .groups = 'drop') %>%
               mutate(year_fac = factor(year)),
             aes(x = factor(1), y = ind_ancestry_count),
             shape = 21, colour = '#b7b7b7', fill = '#e3e3e3',
             size = 0.5, alpha = 0.8,
             position = position_jitter(w = 0.4, h = 0.001)) +
  geom_segment(aes(x = factor(1), xend = factor(1), 
                   y = perc_10, yend = perc_90),
               colour = '#4c4c4c', linewidth = 1.15) +
  geom_point(aes(x = factor(1), y = perc_50),
             shape = 21, color = 'black', 
             fill = 'white', 
             size = 2, stroke = 1.4) +
  facet_grid(~year_fac) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        #strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(0,"cm"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.grid.major.x = element_line(linetype = 'dashed', 
        #                                  color = '#d8d8d8'),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
        axis.ticks = element_line(color = '#808080', linewidth = 0.7),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0, 1), "cm"),
        axis.title.y = element_text(margin = margin(0, 18, 0, 0))
  ) +
  #scale_x_continuous(limits = c(1993.5, 2022.5), position = "top",
  #                  expand = c(0.01, 0)) +
  ylab("Individual ancestry count") +
  coord_cartesian(clip = 'off')


anc_mixing_plot <- anc_mixing_df %>% 
  mutate(year_fac = factor(year)) %>% 
  filter(!is.nan(md)) %>% 
  ggplot() +
  #geom_hline(yintercept = c(0,1), linewidth = 0.4, color = '#b7b7b7') +
  geom_point(aes(x = 1, y = md),
             shape = 21, color = 'black', 
             fill = 'white', 
             size = 2, stroke = 1.4) +
  ylim(-0.1, 1) +
  facet_grid(~year_fac, drop = FALSE) +
  theme_bw() +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = NA, color = "white"),
    strip.text = element_text(size = 7.5),
    panel.spacing = unit(0,"cm"),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major.x = element_line(linetype = 'dashed', 
    #                                  color = '#d8d8d8'),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
    axis.ticks = element_line(color = '#808080', linewidth = 0.7),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.2, 0.1, 0, 1), "cm"),
    axis.title.y = element_text(margin = margin(0, 5.5, 0, 0))
  ) +
  #scale_x_continuous(limits = c(1993.5, 2022.5), position = "bottom",
  #                   expand = c(0.01, 0)) +
  ylab('Translocation ancestry mixing')

indiv_anc_barplot <- pop_dat_with_anc1 %>% 
  #filter(year %in% c(2015, 2016) ) %>% 
  select(year, id, year, anc_group, anc_val, rank_ind) %>% 
  mutate(year_fac = year) %>% 
  rename(Individuals = rank_ind) %>% 
  ggplot() +
  geom_bar(aes(y = anc_val, x = Individuals , fill = anc_group),
           position = "stack", stat = "identity",
           linewidth = 0, width = 1) +
  scale_fill_manual(values = group_color_vec) +
  coord_flip() +
  facet_grid(~year_fac, switch = "x") +
  theme_bw() +
  theme(
    #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
    strip.placement = "outside",
    strip.background = element_rect(fill = NA, color = "white"),
    strip.text = element_text(size = 7.5),
    panel.spacing.x = unit(0,"cm"),
    axis.title.x = element_text(size = 19),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = '#4c4c4c', linewidth = 1.5),
    axis.ticks = element_line(color = '#808080', linewidth = 0.7),
    legend.title = element_blank(),
    plot.margin = unit(c(0.3, 0.1, 0, 1), "cm"),
    axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(lim = c(-0.1, 1.1), expand = c(0, 0)) +
  ylab('Year')

# ancestry_mixing_multipanel <- cowplot::plot_grid(anc_mixing_plot,
#                    anc_count_perc_plot,
#                    indiv_anc_barplot,
#                                       nrow = 3,
#                                       align = "vh",
#                                       axis = 'rltb', #'tblr',
#                                       labels = c('A', 'B', 'C'),
#                                       label_x = 0.09,
#                                       label_y = c(0.89,1),
#                                       label_size = 23,
#                    rel_heights = c(0.3, 0.3, 0.55))

ancestry_mixing_multipanel <- egg::ggarrange(anc_mixing_plot,
                                             anc_count_perc_plot,
                                             indiv_anc_barplot,
                                             nrow = 3,
                                             labels = c('A', 'B', 'C'),
                                             label.args = list(gp = gpar(font = 2, cex = 1.4),
                                                               hjust = -0.5,
                                                               vjust = 1.5),
                                             heights = c(0.3, 0.3, 0.6))


if (isTRUE(output_figs)) {
  cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'ancestry_mixing_multipanel.png'),
                   plot = ancestry_mixing_multipanel,
                   width = 9.75*1.2, height = 8*1.2, bg = 'white')
}



###################################################################
### Suppmat. FIG 4 ###
###################################################################

per_year_growth_rate_plot <- stack(table(results_list$rcw_inbr_merge$year)) %>% 
  mutate(year = as.numeric(as.character(ind)),
         count = as.numeric(values)) %>% 
  select(year, count) %>% 
  arrange(year) %>% 
  mutate(year_boundaries = case_when(year <= min(year) + 1 ~ 'yes',
                                     year == max(year) ~ 'yes',
                                     TRUE ~ 'no'),
         perc_growth_rate = ((count - lag(count))/lag(count))*100) %>% 
  filter(year > 1994) %>% 
  ggplot() +
  geom_bar(aes(x = year, y = perc_growth_rate, alpha = year_boundaries), 
           stat = "identity", fill = '#808080') +
  geom_point(shape = "l", data = results_list$transloc_info_color,
             aes(x = year, y = 0, color = RCWid), size = 5.5,
             alpha = 1,
             position = position_jitter(width = 0.42, height = 0, seed = 1915),
             show.legend = FALSE) +
  scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
                                       nm = results_list$transloc_info_color$RCWid)) +
  geom_hline(yintercept = 0, color = '#4c4c4c', linewidth = 1.25) +
  #geom_line(aes(x = year, y = perc_growth_rate), color = '#808080', linewidth = 0.75) +
  #geom_point(aes(x = year, y = perc_growth_rate),
  #           size = 3) +
  scale_alpha_manual(values = c(1, 0.4)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.border = element_blank(),
        #axis.line.x = element_blank(),
        axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
        axis.ticks = element_line(color = '#808080', linewidth = 0.7),
        legend.position = 'none') +
  xlab('Year') +
  ylab('Per-year growth rate') +
  scale_y_continuous(breaks = seq(-20, 24, by = 5))

cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'per_year_growth_rate_plot.png'),
                 plot = per_year_growth_rate_plot,
                 width = 10.5*1, height = 7*1, bg = 'white')


### Summaries reported in paper ###
summary_doc_path <- here('results', 'text_results', 'ancestry_inbreeding_summary.txt')

cat('##########################\n### ANCESTRY SUMMARIES ###\n##########################',
    file = summary_doc_path, sep="\n", append = FALSE)

pop_dat_with_anc_summ <- pop_dat %>% 
  rename(id = RCWid) %>% 
  left_join(., results_list$rcws_ancestry_info,
            by = 'id',
            relationship = "many-to-many") %>% 
  select(id, year, group, anc_prop) %>% 
  pivot_wider(values_from = anc_prop,
              names_from = group,
              values_fill = 0) %>% 
  group_by(year) %>% 
  #arrange(!!!rlang::syms(c(anc_vec))) %>% 
  arrange(`non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`, 
          .by_group = TRUE) %>% 
  #group_by(year) %>%
  mutate(rank_ind = 1:n()) %>% 
  pivot_longer(cols = c(`non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`), # c(`1`, `2`, `3`),
               names_to = 'anc_group',
               values_to = 'anc_val') %>%
  filter(anc_val > 0) %>% 
  ungroup()


indiv_with_anc_2022 <- pop_dat_with_anc_summ %>% 
  filter(year == 2022) %>% 
  filter(anc_group != 'non-transloc') %>% 
  pull(id) %>% 
  unique() %>% 
  length()

total_2022_indiv <- pop_dat_with_anc_summ %>% 
  filter(year == 2022) %>% 
  pull(id) %>% 
  unique() %>% 
  length()

cat(paste0('2022 indivs with translocation ancestry: ', indiv_with_anc_2022), 
    file = summary_doc_path, 
    sep="\n", append = TRUE)
cat(paste0('total 2022 indiv count: ', total_2022_indiv), 
    file = summary_doc_path, 
    sep="\n", append = TRUE)

cat('\n', 
    file = summary_doc_path, append = TRUE)

indiv_transloc_anc_info <- results_list$rcw_inbr_merge %>% 
  #filter(year == 2022 & anc_prop > 0) %>%
  filter(year == 2022) %>%
  summarize(mean_transloc_anc = round(mean(anc_prop*100), 1),
            sd_transloc_anc = round(sd(anc_prop*100), 1) )

cat('2022 individual-level translocation ancestry:', file = summary_doc_path, sep="\n", append = TRUE)
cat(paste0('mean = ', indiv_transloc_anc_info$mean_transloc_anc), 
    file = summary_doc_path, sep = "\n", append = TRUE)
cat(paste0('sd = ', indiv_transloc_anc_info$sd_transloc_anc), 
    file = summary_doc_path, sep = "\n", append = TRUE)

cat('\n', 
    file = summary_doc_path, append = TRUE)

#pop-level ancestry 
pop_level_total_transloc <- results_list$rcws_ancestry_source_by_year %>% 
  filter(year == 2022) %>% 
  filter(group != "non-transloc") %>%
  pull(ancestry_prop) %>% 
  sum()
cat(paste0('population-level translocation ancestry = ', round(pop_level_total_transloc, 4)), 
    file = summary_doc_path, sep = "\n", append = TRUE)

cat('\n', 
    file = summary_doc_path, append = TRUE)

donor_transloc_info <- results_list$rcws_ancestry_source_by_year %>% 
  filter(year == 2022) %>% 
  filter(group != "non-transloc") %>%
  mutate(group_anc = round(group_anc, 2), 
         ancestry_prop = round(ancestry_prop, 2), 
         prop_transloc = round(100*(ancestry_prop/sum(ancestry_prop)), 2))


write.table(donor_transloc_info[,-1], 
            summary_doc_path,
            append = TRUE,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

cat('\n\n', 
    file = summary_doc_path, append = TRUE)


cat('############################\n### INBREEDING SUMMARIES ###\n############################',
    file = summary_doc_path, sep="\n", append = TRUE)

first_inbr_year <- results_list$rcw_inbr_merge %>% 
  filter(f_ped > 0) %>% 
  pull(year) %>% 
  min()

cat(paste0('first inbreeding year = ', first_inbr_year), 
    file = summary_doc_path, sep = "\n", append = TRUE)

cat(paste0('\nAVERAGE INCREASE IN MEAN Fped STARTING AT FIRST INBREEDING YEAR:'), 
    file = summary_doc_path, sep = "\n", append = TRUE)

fped_pace <- results_list$rcw_inbr_merge %>% 
  filter(year >= first_inbr_year) %>% 
  group_by(year) %>% 
  summarize(mean_fped = mean(f_ped)) %>% 
  arrange(year) %>% 
  mutate(average_fped_dif = mean_fped - lag(mean_fped)) %>% 
  filter(!is.na(average_fped_dif)) %>% 
  summarize(mean_average_fped_change = round(mean(average_fped_dif), 6),
            sd_average_fped_change = round(sd(average_fped_dif), 6))

write.table(fped_pace, 
            summary_doc_path,
            append = TRUE,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)


cat(paste0('\nMAX Fped INFO:'), 
    file = summary_doc_path, sep = "\n", append = TRUE)


max_fped_info <- results_list$rcw_inbr_merge %>% 
  filter(year >= first_inbr_year) %>% 
  group_by(year) %>% 
  summarize(mean_fped = mean(f_ped)) %>% 
  filter(mean_fped == max(mean_fped)) %>% 
  mutate(mean_fped = round(mean_fped, 4))

write.table(max_fped_info, 
            summary_doc_path,
            append = TRUE,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

# results_list$rcw_inbr_merge %>% 
#   filter(year >= first_inbr_year) %>% 
#   group_by(year) %>% 
#   summarize(non_zero_fped_perc = 100*(sum(f_ped > 0)/n())) %>% 
#   arrange(year) %>% 
#   mutate(non_zero_fped_perc_change = non_zero_fped_perc - lag(non_zero_fped_perc)) %>% 
#   filter(!is.na(non_zero_fped_perc_change)) %>% 
#   summarize(mean_non_zero_fped_perc_change = mean(non_zero_fped_perc_change),
#             sd_non_zero_fped_perc_change = sd(non_zero_fped_perc_change))





###################################################################
### FIG 1: POPULATION OVERVIEW (POP SIZE, INBREEDING, ANCESTRY) ###
###################################################################

# ne_plot <- ne_f_df %>%
#   filter(ne < Inf) %>%
#   ggplot() +
#   geom_segment(aes(x = year,
#                    y = ne - ne_std,
#                    xend = year,
#                    yend = ne + ne_std),
#                color = '#b7b7b7', linewidth = 1.5) +
#   geom_point(aes(x = year, y = ne),
#              shape = 21, colour = '#4c4c4c', fill = NA, 
#              size = 1.75, stroke = 1.2) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', 
#                                           color = '#d8d8d8'),
#         panel.border = element_blank(),
#         axis.line.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.margin = unit(c(0.1, 1, -0.05, 1), "cm"),
#         axis.title.y = element_text(margin = margin(0, 4.75, 0, 0))
#         ) +
#   scale_x_continuous(limits = c(1993.5, 2022.5), position = "top",
#                      expand = c(0.01, 0)) +
#   ylab("Ne") +
#   coord_cartesian(clip = 'off')

# rcw_inbr_merge_processed <- results_list$rcw_inbr_merge %>% 
#   group_by(year) %>% 
#   summarize(perc_10 = quantile(f_ped, prob = 0.1),
#             mean = mean(f_ped),
#             perc_90 = quantile(f_ped, prob = 0.9),
#             .groups = 'drop') %>% 
#   arrange(year) %>% 
#   mutate(change_color = case_when(mean > lag(mean) ~ '#ff4c4c',
#                                   mean < lag(mean) ~ '#329932',
#                                   TRUE ~ '#4c4c4c'))

# fped_plot <- results_list$rcw_inbr_merge %>% 
#   ggplot() +
#   geom_point(aes(x = year, y = f_ped),
#              shape = 21, colour = '#b7b7b7', fill = '#e3e3e3', 
#              size = 1, alpha = 0.8,
#              position = position_jitter(w = 0.16, h = 0)) +
#   geom_segment(data = rcw_inbr_merge_processed,
#                aes(x = year, xend = year, 
#                    y = perc_10, yend = perc_90),
#                colour = '#4c4c4c', linewidth = 1.15) +
#   geom_point(data = rcw_inbr_merge_processed,
#              aes(x = year, y = mean),
#              shape = 21, color = rcw_inbr_merge_processed$change_color, 
#              fill = 'white', 
#              size = 2, stroke = 1.4) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', 
#                                           color = '#d8d8d8'),
#         panel.border = element_blank(),
#         axis.line.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.margin = unit(c(0.1, 1, -0.05, 1), "cm"),
#         axis.title.y = element_text(margin = margin(0, 5.5, 0, 0))
#   ) +
#   scale_x_continuous(limits = c(1993.5, 2022.5), position = "top",
#                      expand = c(0.01, 0)) +
#   ylab(expression(italic(F)["p"])) +
#   coord_cartesian(clip = 'off')


# top_bar <- results_list$rcw_inbr_merge %>%
#   mutate(count = 1) %>%
#   group_by(year) %>%
#   arrange(f_ped) %>%
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = f_ped),
#            position = "stack", stat = "identity", width = 0.98) +
#   scale_fill_gradient(name = expression(italic(F)["p"]), 
#                       low = '#f4eaf4', high = '#730073') +
#   geom_point(shape = "l", 
#              #shape = 108, 
#              data = results_list$transloc_info_color,
#              aes(x = year, y = 0, color = RCWid), size = 5.5,
#              alpha = 1,
#              position = position_jitter(width = 0.44, height = 0, seed = 195),
#              show.legend = FALSE) +
#   scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
#                                        nm = results_list$transloc_info_color$RCWid)) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 1))) +
#   theme_bw() +
#   theme(legend.title = element_text(size = 10),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#         panel.border = element_blank(),
#         axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
#         axis.text.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.margin=unit(c(0.32,1,-0.05,1), "cm"),
#         axis.title.y = element_blank()
#   ) +
#   scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
#   coord_cartesian(clip = 'off')

# bottom_bar <- results_list$rcw_inbr_merge %>%
#   mutate(count = 1) %>%
#   group_by(year) %>%
#   arrange(anc_prop) %>%
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = anc_prop),
#            position = "stack", stat = "identity", width = 0.98) +
#   scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = '#ededed', high = '#4c4c4c') +
#   geom_point(shape = "l", data = results_list$transloc_info_color,
#              aes(x = year, y = 0, color = RCWid), size = 5.5,
#              alpha = 1,
#              position = position_jitter(width = 0.44, height = 0, seed = 195),
#              show.legend = FALSE) +
#   scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
#                                        nm = results_list$transloc_info_color$RCWid)) +
#   theme_bw() +
#   theme(legend.title = element_text(size = 10),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#         panel.border = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5)
#   ) +
#   scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(1, 0))) +
#   theme(plot.margin=unit(c(-0.05,1,-0.05,1), "cm"),
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank()#,
#   ) +
#   scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
#   coord_cartesian(clip = 'off')
# 
# 
# transloc_ancestry <- results_list$rcws_ancestry_source_by_year %>%
#   mutate(group = factor(group, 
#                 levels = rev(c('ANF', 'CBJTC', 'FTB', 'FTS', 'ONF', 'WSF-CITRUS', 'non-transloc')))) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = ancestry_prop, fill = group), 
#            color = 'white', linewidth = 0.15, width = 0.72,
#            position = "stack", stat = "identity") + 
#   #scale_fill_manual(values = c('#ededed', '#4c4c4c')) +
#   scale_fill_manual(values = group_color_vec) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#         panel.border = element_blank(),
#         axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#         plot.margin = unit(c(0.25,1,0,1), "cm"),
#         axis.title.x = element_text(size = 18),
#         legend.position = 'none',
#         legend.title = element_blank(),
#         axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0))
#   ) +
#   scale_y_continuous(breaks = c(0, 0.5, 1), 
#                      labels = c("0.0", "0.5", "1.0")) +
#   scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
#   xlab('Year') +
#   ylab("Ancestry") +
#   theme(legend.key.size = unit(0.4, 'cm')) +
#   coord_cartesian(clip = 'off')
# 
# ### custom legends ###
# fped_plot_legend <- ggdraw(cowplot::get_legend(
#   data.frame(x = 1,
#              y = 3:1,
#              label = factor(c('increase', 'no change', 'decrease'), 
#                             levels = c('increase', 'no change', 'decrease')) ) %>% 
#     ggplot() +
#     geom_point(aes(x = x,
#                    y = y,
#                    color = label),
#                shape = 21,
#                fill = 'white',
#                size = 2, stroke = 1.4) +
#     theme_bw() +
#     theme(legend.key = element_blank(),
#           legend.title = element_blank(),
#           legend.background = element_rect(fill = "transparent"),
#           panel.background = element_rect(fill = "transparent", color = "transparent")) +
#     scale_color_manual(substitute(paste('Change in mean(', italic(F)["p"], ")")),
#                        values = c('#ff4c4c', '#4c4c4c', '#329932'))
#   
# )) +
#   theme(plot.margin=unit(c(0,-5, 0, -5), "cm"),
#         panel.background = element_rect(fill = "transparent", color = "transparent"))
# 
# founder_prop_legend <- results_list$rcws_founder_info %>%
#   group_by(group) %>% 
#   summarize(count = n(),
#             .groups = 'drop') %>%
#   mutate(group = factor(group, 
#                         levels = rev(c('CBJTC', 'ANF', 'FTB', 'WSF-CITRUS', 'FTS', 'ONF', 'non-transloc')))) %>% 
#   mutate(x = 0.95) %>%
#   arrange(desc(group)) %>% 
#   mutate(cumsum_count = cumsum(count) - (0.5*count)) %>% 
#   mutate(cumsum_count_adjust = cumsum_count + c(-13, -5, 0, 10, 17, 24, 0)) %>% 
#   mutate(update_label = paste0(group, ' (', count, ')')) %>% 
#   ggplot(aes(x = x, y = count, fill = group, label = count)) +
#   geom_segment(aes(x = 1.47,
#                    xend = 2,
#                    y = cumsum_count, 
#                    yend = cumsum_count_adjust,
#                    color = group),
#                linewidth = 1.1) +
#   geom_bar(position = 'stack', stat = "identity", color = 'white', linewidth = 0.15) +
#   scale_color_manual(values = group_color_vec) +
#   scale_fill_manual(values = group_color_vec) +
#   geom_text(aes(x = 2.1, y = cumsum_count_adjust, label = update_label),
#             size = 2.8,
#             hjust = 0) +
#   theme_void() +
#   theme(legend.position = 'none') +
#   coord_cartesian(clip = 'off') +
#   xlim(-0.32, 4) +
#   geom_text(aes(x = 2.76, y = 207, label = "Founder proportions"),
#             size = 3.5)
# 
# 
# ### add in legends and RCW image ###
# rcw_pop_overview <- ggpubr::annotate_figure(egg::ggarrange(fped_plot,
#                                                            top_bar,
#                                                            bottom_bar,
#                                                            transloc_ancestry,
#                                                            nrow = 4,
#                                                            labels = c('A', 'B', '', 'C'),
#                                                            label.args = list(gp = grid::gpar(font = 2, cex = 1.2)),
#                                                            heights = c(30, 50, 50, 22)),
#                                             left = grid::textGrob("Individuals", rot = 90, hjust = 0.58, vjust = 5.2, gp = grid::gpar(cex = 1))) +
#   theme(plot.margin = unit(c(0.3, -0.8, 0, -3), "cm"))
# 
# rcw_pop_overview_rcw <- rcw_pop_overview +
#   cowplot::draw_image(flying_rcw,
#                       scale = 0.12,
#                       x = 0.325, y = 0.432) +
#   patchwork::inset_element(founder_prop_legend,
#                            left = 0.849, bottom = -0.01, right = 0.98, top = 0.20
#                            #left = 0.849, bottom = -0.01, right = 0.972, top = 0.20
#                            ) +
#   patchwork::inset_element(fped_plot_legend,
#                            #left = 0, bottom = 0.9, right = 0.38, top = 0.95,
#                            left = 0.891, bottom = 0.81, right = 0.94, top = 0.88)
# 
# 
# if (isTRUE(output_figs)) {
#   cowplot::ggsave2(filename = here('figures', 'main_paper', 'rcw_pop_overview_fig.png'),
#                    plot = ggdraw(rcw_pop_overview_rcw) + theme(plot.margin = unit(c(-0.3, 0.3, 0.25, -0.6), "cm")),
#                    width = 10*1, height = 8.1*1, bg = 'white')
#   
#   # cowplot::ggsave2(filename = here('figures', 'main_paper', 'rcw_pop_overview_fig.png'),
#   #                  plot = rcw_pop_overview_rcw,
#   #                  width = 10*1, height = 8.1*1, bg = 'white')
# }
# 
# 
# # pop_dat %>% 
# #   group_by(year) %>% 
# #   summarize(count = n()) %>% 
# #   print(n = 50)
# #################################################################
# ### FIG 2: GENETIC CONTRIBUTIONS AND MORE INFO. ON INBREEDING ###
# #################################################################
# fped_contr_transloc_ids <- results_list$rcw_partial_founder_summed_processed %>% 
#   filter(year == 2022) %>% 
#   filter(group != 'non-transloc') %>% 
#   pull(founder_id)
# 
# 
# 
# fig2_list <- list()
# 
# for (i in c('sex', 'source')) {
#   
#   fig2_list[[i]][['contr_plot']] <- results_list$rcw_contr_info_df %>%
#     group_by(id) %>%
#     filter(any(contr != 0)) %>% #filter out the ids with all 0s
#     filter(contr > 0 | (contr == 0 & year > max(year[contr > 0])  )) %>%
#     left_join(., results_list$rcws_founder_info,
#               by = 'id') %>%
#     ggplot() +
#     geom_line(data = . %>% 
#                 filter(group == 'non-transloc'),
#               aes(x = year, y = contr, group = id, linewidth = `group`),
#               color = '#e3e3e3', linewidth = 0.7) +
#     geom_point(data = . %>%
#                  group_by(id) %>%
#                  filter(year == min(year)) %>%
#                  slice_head(n = 1)%>% 
#                  filter(group == 'non-transloc'),
#                aes(x = year, y = contr),
#                color = '#e3e3e3', size = 0.7) +
#     geom_line(data = . %>% 
#                 filter(group != 'non-transloc'),
#               aes(x = year, y = contr, group = id, 
#                   linewidth = group, color = id),
#               linewidth = 1.1) +
#     geom_point(data = . %>%
#                  group_by(id) %>%
#                  filter(year == min(year)) %>%
#                  slice_head(n = 1)%>% 
#                  filter(group != 'non-transloc'),
#                aes(x = year, y = contr, color = id),
#                size = 1.5) +
#     scale_color_manual(values = setNames(results_list$transloc_info_color[[paste0(i, '_color')]],
#                                          nm = results_list$transloc_info_color$RCWid)) +
#     theme_bw() +
#     theme(panel.grid.major.y = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#           panel.border = element_blank(),
#           axis.line.x = element_blank(),
#           axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#           axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#           axis.title.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           plot.margin = unit(c(0.32,1, 0, 0.5), "cm"),
#           legend.position = 'none',
#           axis.title.y = element_text(size = 12, margin = margin(0, 12, 0, 0))
#     ) +
#     scale_x_continuous(expand = c(0, 0),
#                        #expand = c(0.01, 0), 
#                        #limits = c(1993.5, 2022.5),
#                        limits = c(1993.5, 2023.5),
#                        position = "top") +
#     xlab('Year') +
#     ylab("Expected genetic contribution") +
#     geom_text(data = results_list$rcw_contr_info_df %>%
#                                filter(year == 2022 & id %in% fped_contr_transloc_ids) %>%
#                                #filter(year == 2022 & id %in% "ZG-GWO") %>%
#                                mutate(y_pos = case_when(id == "OWO-AZ" ~ contr + 0.0005,
#                                                         id == "ZG-GWO" ~ contr - 0.003,
#                                                         id == 'WZ-YKA' ~ contr + 0.001,
#                                                         TRUE ~ contr)),
#                              aes(x = 2022.6, y = y_pos, label = id),
#                              size = 2.55,
#                              hjust = 0) +
#     geom_segment(data = results_list$rcw_contr_info_df %>%
#                    filter(year == 2022 & id %in% fped_contr_transloc_ids) %>%
#                    #filter(id %in% c("OWO-AZ", "ZG-GWO", 'WZ-YKA')) %>%
#                    mutate(yend = case_when(id == "OWO-AZ" ~ contr + 0.0005,
#                                             id == "ZG-GWO" ~ contr - 0.003,
#                                             id == 'WZ-YKA' ~ contr + 0.001,
#                                            TRUE ~ contr)),
#                  
#                  aes(x = 2022.05, xend = 2022.5,
#                      y = contr, yend = yend)) +
#     coord_cartesian(clip = 'off')
#     # ggrepel::geom_text_repel(data = results_list$rcw_contr_info_df %>%
#     #             filter(year == 2022 & id %in% fped_contr_transloc_ids),
#     #           aes(x = 2022.1, y = contr, label = id),
#     #           size = 2.5,
#     #           hjust = 0,
#     #           min.segment.length = 0.3,
#     #           box.padding = 0.4,
#     #           seed = 4534,
#     #           direction = 'y',
#     #           nudge_y = setNames(c(      -0.001,   0.001,    -0.002,   -0.003,   0.0001), 
#     #                              nm = c('ZG-GWO', 'OHA-ZK', 'OWO-AZ', 'WZ-YKA', 'OAG-HZ'))
#     # )
# 
#   
#   fig2_list[[i]][['fped_prop_by_indiv_founder_plot']] <- results_list$rcw_partial_founder_summed_processed %>% 
#     ggplot() +
#     geom_bar(aes(x = year, y = pfound_fsum_scale, fill = founder_id,
#                  group = -pfound_fsum_scale),
#              color = 'white', linewidth = 0.2,
#              position = "stack", stat = "identity", width = 0.98, alpha = 0.85) +
#     #geom_point(data = . %>% filter(group == 'transloc'),
#     geom_point(data = . %>% filter(group != 'non-transloc'),
#                aes(x = year, y = val_cumsum, color = founder_id),
#                size = 1.5) +
#     geom_line(#data = . %>% filter(group == 'transloc'),
#       data = . %>% filter(group != 'non-transloc'),
#       aes(x = year, y = val_cumsum, 
#           group = founder_id, color = founder_id),
#       linewidth = 0.7, linetype = 'solid') +
#     scale_fill_manual(values = setNames(results_list$inbr_founder_color[[paste0(i, '_color')]],
#                                         nm = results_list$inbr_founder_color$founder_id)) +
#     scale_color_manual(values = setNames(results_list$inbr_founder_color[[paste0(i, '_color')]],
#                                          nm = results_list$inbr_founder_color$founder_id)) +
#     theme_bw() +
#     theme(
#       panel.grid.major.y = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#       panel.border = element_blank(),
#       axis.line.x = element_blank(),
#       axis.text.x = element_blank(),
#       axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#       axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#       axis.title.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       plot.margin = unit(c(0.25,1,-0.05,0.5), "cm"),
#       legend.position = 'none',
#       axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 4.5, b = 0, l = 0))) +
#     ylab(substitute(paste("Total ", italic(F)["p"]))) +
#     #ylab(substitute(italic(F)["P"])) +
#     scale_x_continuous(expand = c(0, 0),
#                        #expand = c(0.01, 0), 
#                        limits = c(1993.5, 2023.5)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     coord_cartesian(clip = 'off') +
#     geom_text(data = results_list$rcw_partial_founder_summed_processed %>% 
#                 filter(founder_id %in% fped_contr_transloc_ids & year == 2022),
#               aes(x = 2022.6, y = val_cumsum, label = founder_id),
#               size = 2.5,
#               hjust = 0)
#   
#   fig2_list[[i]][['fped_prop_by_group_plot']] <- results_list$fped_prop_by_group %>%
#     mutate(group = case_when(group == "non-transloc" ~ "non-transloc",
#                              TRUE ~ 'transloc')) %>% 
#     group_by(year, group) %>% 
#     summarize(fped_group_prop = sum(fped_group_prop),
#               .groups = 'drop') %>% 
#     ggplot() +
#     geom_bar(aes(x = year, y = fped_group_prop, fill = group), 
#              color = 'white', linewidth = 0.35, width = 0.70,
#              position = "stack", stat = "identity") + 
#     scale_fill_manual(values = c('#ededed', '#4c4c4c')) +
#     theme_bw() +
#     theme(panel.grid.major.y = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#           panel.border = element_blank(),
#           #axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
#           axis.line.x = element_blank(),
#           axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#           plot.margin = unit(c(0.25,1,0,0.5), "cm"),
#           axis.title.x = element_text(size = 18),
#           legend.position = 'none',
#           axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 14, b = 0, l = 0))
#     ) +
#     # scale_y_continuous(breaks = c(0, 0.5, 1), 
#     #                    labels = c("0.0", "0.5", "1.0")) +
#     # scale_x_continuous(expand = c(0.01, 0),
#     #                    limits = c(1993.5, 2023.5)
#     #                    #limits = c(1993.5, 2022.5)
#     #                    ) +
#     # xlab('Year') +
#     # ylab(substitute(paste(italic(F)["p"], ' prop.')))
#     scale_y_continuous(expand = c(0, 0),
#                        breaks = c(0, 0.5, 1),
#                        limits = c(-0.05, 1.02),
#                        labels = c("0.0", "0.5", "1.0")) +
#     scale_x_continuous(expand = c(0, 0),
#                        #expand = c(0.01, 0),
#                        limits = c(1993.5, 2023.5)
#                        #limits = c(1993.5, 2022.5)
#     ) +
#     xlab('Year') +
#     ylab(substitute(paste(italic(F)["P"], ' prop.'))) +
#     geom_segment(aes(x = 1993.5, xend = 2022.5,
#                      y = -0.05, yend = -0.05),
#                  color = '#4c4c4c', linewidth = 1.1) +
#     coord_cartesian(clip = 'off')
#   
#   fig2_list[[i]][['panel1']] <- egg::ggarrange(fig2_list[[i]][['contr_plot']],
#                                                fig2_list[[i]][['fped_prop_by_indiv_founder_plot']],
#                                                fig2_list[[i]][['fped_prop_by_group_plot']],
#                                                nrow = 3,
#                                                labels = c('A', 'B', 'C'),
#                                                label.args = list(gp = grid::gpar(font = 2, cex = 1.3)),
#                                                heights = c(65, 65, 12))
#   
#   fig2_list[[i]][['panel2']] <- cowplot::ggdraw(fig2_list[[i]][['panel1']]) +
#     cowplot::draw_image(perched_rcw,
#                         scale = 0.2,
#                         #x = -0.3695, 
#                         x = -0.406, 
#                         y = -.13) +
#     theme(plot.margin = unit(c(0.1, -0.7, 0, 0.1), "cm"))
#   
#   if (isTRUE(output_figs)) {
#     cowplot::ggsave2(filename = here('figures', 'main_paper', paste0('indiv_fped_gencontr_panel_', i, '.png') ),
#                      plot = fig2_list[[i]][['panel2']],
#                      width = 11*1.1, height = 5.9*1.1, bg = 'white')
#   }
# }
# 
# 
# 
# 
# 
# census_processed <- read.csv(here('data', 'feb2024_databasemarch2023_processed', 'census_processed.csv'))
# pop_info_list <- lapply(setNames(nm = paste0('Scenario', 1:4)), function(X, DF) {
#   DF[DF[,X,drop = TRUE],!grepl(pattern = 'Scenario', colnames(DF))]
# }, DF = census_processed)
# pop_dat <- pop_info_list$Scenario3
# 
# 
# 
# 
# 
# results_list$fped_prop_by_group %>%
#   mutate(group = case_when(group == "non-transloc" ~ "non-transloc",
#                            TRUE ~ 'transloc')) %>% 
#   group_by(year, group) %>% 
#   summarize(fped_group_prop = sum(fped_group_prop),
#             .groups = 'drop') %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = fped_group_prop, fill = group), 
#            color = 'white', linewidth = 0.35, width = 0.70,
#            position = "stack", stat = "identity") + 
#   scale_fill_manual(values = c('#ededed', '#4c4c4c')) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#         panel.border = element_blank(),
#         #axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
#         axis.line.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#         plot.margin = unit(c(0.25,1,0,1), "cm"),
#         axis.title.x = element_text(size = 18),
#         legend.position = 'none',
#         axis.title.y = element_text(margin = margin(t = 0, r = 14, b = 0, l = 0))
#   ) +
#   scale_y_continuous(expand = c(0, 0.01),
#                      breaks = c(0, 0.5, 1),
#                      limits = c(-0.05, 1.02),
#                      labels = c("0.0", "0.5", "1.0")) +
#   scale_x_continuous(expand = c(0, 0),
#                      limits = c(1993.5, 2023.5)
#                      #limits = c(1993.5, 2022.5)
#   ) +
#   xlab('Year') +
#   ylab(substitute(paste(italic(F)["p"], ' prop.'))) +
#   geom_segment(aes(x = 1993.5, xend = 2023,
#                    y = -0.05, yend = -0.05),
#              color = '#4c4c4c', linewidth = 1.5) +
#   coord_cartesian(clip = 'off')
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
# fig2_list[[i]][['fped_prop_by_indiv_founder_plot']] +
#   geom_text(data = results_list$rcw_partial_founder_summed_processed %>% 
#               filter(founder_id %in% fped_contr_transloc_ids & year == 2022),
#             aes(x = 2025, y = val_cumsum, label = founder_id),
#             size = 3) +
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#   coord_cartesian(clip = 'off')


# examp_indivs <- pop_dat %>% 
#   filter(year == 2010) %>% 
#   pull(RCWid)
# 
# anc_info_test <- results_list$rcws_ancestry_info %>% 
#   filter(id %in% examp_indivs)
# 
# 
# total_anc_info <- anc_info_test %>% 
#   group_by(group) %>% 
#   summarize(anc_sum = sum(anc_prop), .groups = 'drop') %>% 
#   mutate(anc_prop = anc_sum/sum(anc_sum))
# 
# 
# 
# anc_mixing_df <- data.frame(year = sort(unique(pop_dat$year)),
#                             md = NA)
# 
# for (YEAR in anc_mixing_df$year) {
# #for (YEAR in 1999) {  
#   focal_year_indivs <- pop_dat %>% 
#     filter(year == YEAR) %>% 
#     pull(RCWid)
#   
#   anc_info_focal_year <- results_list$rcws_ancestry_info %>% 
#     filter(id %in% focal_year_indivs)
#   
#   total_anc_info <- anc_info_focal_year %>% 
#     group_by(group) %>% 
#     summarize(anc_sum = sum(anc_prop), .groups = 'drop') %>% 
#     mutate(anc_prop = anc_sum/sum(anc_sum))
#   
#   if (nrow(total_anc_info) == 1) {
#     anc_mixing_df$md[anc_mixing_df$year == YEAR] <- NaN
#     next
#   }
#   
#   var_vec <- setNames(rep(NA, length(unique(total_anc_info$group))), nm = unique(total_anc_info$group))
#   for (i in names(var_vec)) {
#     
#     anc_subset <- anc_info_focal_year %>% 
#       filter(group == i)
#     
#     var_vec[i] <- var(c(anc_subset$anc_prop, rep(0, length(unique(anc_info_focal_year$id[!anc_info_focal_year$id %in% anc_subset$id])))))
#   }
#   
#   var_vec_remove <- var_vec[names(var_vec) != 'non-transloc']
#   total_anc_info_remove <- total_anc_info %>% 
#     filter(group != 'non-transloc')
#   
#   anc_mixing_df$md[anc_mixing_df$year == YEAR] <- 1 - (sum(var_vec_remove)/sum(sapply(total_anc_info_remove$anc_prop, function(x) x*(1 - x))))
#   
#   #anc_mixing_df$md[anc_mixing_df$year == YEAR] <- 1 - (sum(var_vec)/sum(sapply(total_anc_info$anc_prop, function(x) x*(1 - x))))
#   
# }
# 
# anc_by_year_quantile_summary <- pop_dat %>% 
#   rename(id = RCWid) %>% 
#   left_join(., results_list$rcws_ancestry_info,
#           by = 'id',
#           relationship = "many-to-many") %>% 
#   group_by(year, id) %>% 
#   summarize(ind_ancestry_count = n(), .groups = 'drop') %>% 
#   group_by(year) %>% 
#   summarize(perc_5 = quantile(ind_ancestry_count, prob = 0.05),
#             perc_10 = quantile(ind_ancestry_count, prob = 0.10),
#             perc_25 = quantile(ind_ancestry_count, prob = 0.25),
#             perc_50 = quantile(ind_ancestry_count, prob = 0.50),
#             perc_75 = quantile(ind_ancestry_count, prob = 0.75),
#             perc_90 = quantile(ind_ancestry_count, prob = 0.90),
#             perc_95 = quantile(ind_ancestry_count, prob = 0.95),
#             .groups = 'drop')
# 
# 
# anc_count_perc_plot <- anc_by_year_quantile_summary %>% 
#   mutate(year_fac = factor(year)) %>% 
#   ggplot() +
#   geom_point(data = pop_dat %>%
#                rename(id = RCWid) %>%
#                left_join(., results_list$rcws_ancestry_info,
#                          by = 'id',
#                          relationship = "many-to-many") %>%
#                group_by(year, id) %>%
#                summarize(ind_ancestry_count = n(), .groups = 'drop') %>%
#                mutate(year_fac = factor(year)),
#              aes(x = factor(1), y = ind_ancestry_count),
#              shape = 21, colour = '#b7b7b7', fill = '#e3e3e3',
#              size = 0.5, alpha = 0.8,
#              position = position_jitter(w = 0.4, h = 0.001)) +
#   geom_segment(aes(x = factor(1), xend = factor(1), 
#                  y = perc_10, yend = perc_90),
#              colour = '#4c4c4c', linewidth = 1.15) +
#   geom_point(aes(x = factor(1), y = perc_50),
#              shape = 21, color = 'black', 
#              fill = 'white', 
#              size = 2, stroke = 1.4) +
#   facet_grid(~year_fac) +
#   theme_bw() +
#   theme(strip.background = element_blank(), 
#         strip.text = element_blank(),
#         axis.text.x = element_blank(),
#         #strip.background = element_rect(fill = NA, color = "white"),
#         panel.spacing = unit(0,"cm"),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.grid.major.x = element_line(linetype = 'dashed', 
#         #                                  color = '#d8d8d8'),
#         panel.grid.major.x = element_blank(),
#         panel.border = element_blank(),
#         axis.line.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.margin = unit(c(0.3, 0.1, 0, 1), "cm"),
#         axis.title.y = element_text(margin = margin(0, 18, 0, 0))
#   ) +
#   #scale_x_continuous(limits = c(1993.5, 2022.5), position = "top",
#    #                  expand = c(0.01, 0)) +
#   ylab("Individual ancestry count") +
#   coord_cartesian(clip = 'off')
# 
# 
# anc_mixing_plot <- anc_mixing_df %>% 
#   mutate(year_fac = factor(year)) %>% 
#   filter(!is.nan(md)) %>% 
#   ggplot() +
#   #geom_hline(yintercept = c(0,1), linewidth = 0.4, color = '#b7b7b7') +
#   geom_point(aes(x = 1, y = md),
#              shape = 21, color = 'black', 
#              fill = 'white', 
#              size = 2, stroke = 1.4) +
#   ylim(-0.1, 1) +
#   facet_grid(~year_fac, drop = FALSE) +
#   theme_bw() +
#   theme(
#     strip.placement = "outside",
#     strip.background = element_rect(fill = NA, color = "white"),
#     strip.text = element_text(size = 7.5),
#     panel.spacing = unit(0,"cm"),
#     panel.grid.major.y = element_blank(),
#     axis.text.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.grid.major.x = element_line(linetype = 'dashed', 
#         #                                  color = '#d8d8d8'),
#     panel.grid.major.x = element_blank(),
#         panel.border = element_blank(),
#         axis.line.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.margin = unit(c(0.2, 0.1, 0, 1), "cm"),
#         axis.title.y = element_text(margin = margin(0, 5.5, 0, 0))
#   ) +
#   #scale_x_continuous(limits = c(1993.5, 2022.5), position = "bottom",
#   #                   expand = c(0.01, 0)) +
#   ylab('Translocation ancestry mixing')
#   
# indiv_anc_barplot <- pop_dat_with_anc1 %>% 
#   #filter(year %in% c(2015, 2016) ) %>% 
#   select(year, id, year, anc_group, anc_val, rank_ind) %>% 
#   mutate(year_fac = year) %>% 
#   rename(Individuals = rank_ind) %>% 
#   ggplot() +
#   geom_bar(aes(y = anc_val, x = Individuals , fill = anc_group),
#            position = "stack", stat = "identity",
#            linewidth = 0, width = 1) +
#   scale_fill_manual(values = group_color_vec) +
#   coord_flip() +
#   facet_grid(~year_fac, switch = "x") +
#   theme_bw() +
#   theme(
#     #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
#         strip.placement = "outside",
#         strip.background = element_rect(fill = NA, color = "white"),
#         strip.text = element_text(size = 7.5),
#        panel.spacing.x = unit(0,"cm"),
#         axis.title.x = element_text(size = 19),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(color = '#4c4c4c', linewidth = 1.5),
#        axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#        legend.title = element_blank(),
#        plot.margin = unit(c(0.3, 0.1, 0, 1), "cm"),
#        axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
#        ) +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(lim = c(-0.1, 1.1), expand = c(0, 0)) +
#   ylab('Year')
# 
# # ancestry_mixing_multipanel <- cowplot::plot_grid(anc_mixing_plot,
# #                    anc_count_perc_plot,
# #                    indiv_anc_barplot,
# #                                       nrow = 3,
# #                                       align = "vh",
# #                                       axis = 'rltb', #tblr',
# #                                       labels = c('A', 'B', 'C'),
# #                                       label_x = 0.09,
# #                                       label_y = c(0.89,1),
# #                                       label_size = 23,
# #                    rel_heights = c(0.3, 0.3, 0.55))
# 
# ancestry_mixing_multipanel <- egg::ggarrange(anc_mixing_plot,
#                                              anc_count_perc_plot,
#                                              indiv_anc_barplot,
#                                              nrow = 3,
#                                              labels = c('A', 'B', 'C'),
#                                              label.args = list(gp = gpar(font = 2, cex = 1.4),
#                                                                hjust = -0.5,
#                                                                vjust = 1.5),
#                                               heights = c(0.3, 0.3, 0.6))
# 
#                    
# if (isTRUE(output_figs)) {
#   cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'ancestry_mixing_multipanel.png'),
#                    plot = ancestry_mixing_multipanel,
#                    width = 9.75*1.2, height = 8*1.2, bg = 'white')
# }
# 
# 
# test_df_1 <- data.frame(year = c(1, 1, 1, 1, 2, 2, 2),
#                         val = c(10, 30, 40, 1, 10, 11, 2),
#                         id = c(1, 2, 1, 2, 4, 1, 2))
#                         
# test_df_2 <- data.frame(year = 2,
#                         id = c('a', 'a', 'a', 'b', 'b', 'c'),
#                         anc = c('1', '2', '3', '1', '3', '1'),
#                         val = c(0.8, 0.7, 0.1, 0.5, 0.5, 1))
# 
# 
# 
# pop_indiv_anc_processed <- pop_dat %>% 
#   rename(id = RCWid) %>% 
#   left_join(., results_list$rcws_ancestry_info,
#             by = 'id',
#             relationship = "many-to-many") %>%
#   pivot_wider(values_from = anc_prop,
#               names_from = group,
#               values_fill = 0) %>% 
#   group_by(year) %>% 
#   arrange(year, `non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`, 
#           .by_group = TRUE) %>% 
#   group_by(year) %>%
#   mutate(rank_ind = 1:n()) %>% 
#   pivot_longer(cols = unique(results_list$rcws_ancestry_info$group), # c(`1`, `2`, `3`),
#                names_to = 'anc_group',
#                values_to = 'anc_val') %>%
#   filter(anc_val > 0) %>% 
#   ungroup()
# 
# 
# pop_dat_with_anc <- pop_dat %>% 
#   rename(id = RCWid) %>% 
#   left_join(., results_list$rcws_ancestry_info,
#             by = 'id',
#             relationship = "many-to-many")
# 
# data.frame(col1 = c('b', 'a'),
#            col2 = 2:1) %>% 
#   #arrange(.dots=c("col1", "col2"))
#   #arrange(col1, col2)
#   arrange(!!!rlang::syms(c("col1", "col2")))
# 
# 
# all_anc_vec <- c("CBJTC","ANF", "FTB", "WSF-CITRUS", "FTS", "ONF", "non-transloc")
# 
# pop_dat_with_anc1 <- lapply(split(pop_dat_with_anc, pop_dat_with_anc$year), function(x, all_anc_vec) { 
#   anc_vec <- all_anc_vec[all_anc_vec %in% unique(x$group)]
#   
#   x %>% 
#     select(id, year, group, anc_prop) %>% 
#   pivot_wider(values_from = anc_prop,
#               names_from = group,
#               values_fill = 0) %>% 
#   #group_by(year) %>% 
#   arrange(!!!rlang::syms(c(anc_vec))) %>% 
#   #arrange(`non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`, 
#   #        .by_group = TRUE) %>% 
#   #group_by(year) %>%
#   mutate(rank_ind = 1:n()) %>% 
#   pivot_longer(cols = all_of(anc_vec), # c(`1`, `2`, `3`),
#                names_to = 'anc_group',
#                values_to = 'anc_val') %>%
#   filter(anc_val > 0) %>% 
#   ungroup()
# }, all_anc_vec = all_anc_vec) %>% 
#   bind_rows()
# 
# pop_dat_with_anc <- pop_dat %>% 
#   rename(id = RCWid) %>% 
#   left_join(., results_list$rcws_ancestry_info,
#             by = 'id',
#             relationship = "many-to-many") %>% 
#   select(id, year, group, anc_prop) %>% 
#   pivot_wider(values_from = anc_prop,
#               names_from = group,
#               values_fill = 0) %>% 
#   group_by(year) %>% 
#   #arrange(!!!rlang::syms(c(anc_vec))) %>% 
#   arrange(`non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`, 
#           .by_group = TRUE) %>% 
#   #group_by(year) %>%
#   mutate(rank_ind = 1:n()) %>% 
#   pivot_longer(cols = c(`non-transloc`, `ANF`, `FTB`, `ONF`, `FTS`, `CBJTC`, `WSF-CITRUS`), # c(`1`, `2`, `3`),
#                names_to = 'anc_group',
#                values_to = 'anc_val') %>%
#   filter(anc_val > 0) %>% 
#   ungroup()
# 
# 
# 
#  
# # test_plot <- test_df %>%   
# #   ggplot() +
# #   geom_bar(aes(y = anc_prop, x = pos, fill = group, group = year),
# #            position="stack", stat="identity",
# #            linewidth = 1, width = 1) +
# #   coord_flip() +
# #   facet_grid(~year, switch = "x") +
# #   theme_bw() +
# #   theme(strip.placement = "outside",
# #         strip.background = element_rect(fill = NA, color = "white"),
# #         panel.spacing = unit(-.8,"cm"),
# #         panel.grid.major = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         #panel.grid.major.x = element_line(linetype = 'dashed', 
# #         #                                  color = '#d8d8d8'),
# #         panel.border = element_blank(),
# #         #axis.line.x = element_blank(),
# #         axis.line = element_line(color = '#4c4c4c', linewidth = 1.5),
# #         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
# #         axis.title.x = element_blank(),
# #         axis.ticks.x = element_blank(),
# #         axis.text.x = element_blank()
# #         ) +
# #   scale_x_continuous(expand = c(0.01, 0))
# #   #scale_x_continuous(limits = c(1993.5, 2022.5), position = "top",
# #   #                   expand = c(0.01, 0))
# # 
# # 
# # plot_grid(test_plot,
# #           ggplot(data.frame(year = 1:3,
# #                             val = 1:3,
# #                             group = 'x')) +
# #             geom_point(aes(x = 1, y = val)) +
# #             facet_grid(~year) +
# #             theme(panel.spacing = unit(-.8,"cm")),
# #           nrow = 2,
# #           align = "v",
# #           axis = 'tblr')
# # 
# 
# 
# per_year_growth_rate_plot <- stack(table(results_list$rcw_inbr_merge$year)) %>% 
#   mutate(year = as.numeric(as.character(ind)),
#          count = as.numeric(values)) %>% 
#   select(year, count) %>% 
#   arrange(year) %>% 
#   mutate(year_boundaries = case_when(year <= min(year) + 1 ~ 'yes',
#                                      year == max(year) ~ 'yes',
#                                      TRUE ~ 'no'),
#          perc_growth_rate = ((count - lag(count))/lag(count))*100) %>% 
#   filter(year > 1994) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = perc_growth_rate, alpha = year_boundaries), 
#            stat = "identity", fill = '#808080') +
#   geom_point(shape = "l", data = results_list$transloc_info_color,
#              aes(x = year, y = 0, color = RCWid), size = 5.5,
#              alpha = 1,
#              position = position_jitter(width = 0.42, height = 0, seed = 1915),
#              show.legend = FALSE) +
#   scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
#                                        nm = results_list$transloc_info_color$RCWid)) +
#   geom_hline(yintercept = 0, color = '#4c4c4c', linewidth = 1.25) +
#   #geom_line(aes(x = year, y = perc_growth_rate), color = '#808080', linewidth = 0.75) +
#   #geom_point(aes(x = year, y = perc_growth_rate),
#   #           size = 3) +
#   scale_alpha_manual(values = c(1, 0.4)) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#         panel.border = element_blank(),
#         #axis.line.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#         legend.position = 'none') +
#   xlab('Year') +
#   ylab('Per-year growth rate') +
#   scale_y_continuous(breaks = seq(-20, 24, by = 5))
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'per_year_growth_rate_plot.png'),
#                  plot = per_year_growth_rate_plot,
#                  width = 10.5*1, height = 7*1, bg = 'white')
# 
# 
# 
# 
# 
#  
# 
# transloc_abund_vs_ancestry_plot <- results_list$rcws_founder_info %>%
#   group_by(group) %>% 
#   summarize(count = n(),
#             .groups = 'drop') %>%
#   mutate(group = factor(group, 
#                         levels = rev(c('CBJTC', 'ANF', 'FTB', 'WSF-CITRUS', 'FTS', 'ONF', 'non-transloc')))) %>% 
#   left_join(., results_list$rcws_ancestry_source_by_year %>% filter(year == 2022) %>% 
#               select(!X),
#             by = c('group')) %>% 
#   arrange(desc(group)) %>% 
#   filter(group != 'non-transloc') %>% 
#   rename(Translocations = count,
#          `2022 estimated ancestry` = ancestry_prop) %>% 
#   mutate(`2022 estimated ancestry` = 100*`2022 estimated ancestry`) %>% 
#   pivot_longer(cols = c(Translocations, `2022 estimated ancestry`),
#                values_to = 'Percentage',
#                names_to = 'Type') %>% 
#   mutate(Type = factor(Type, levels = c("Translocations", "2022 estimated ancestry"))) %>% 
#   ggplot() +
#   geom_line(aes(x = Type, y = Percentage, group = group, color = group),
#             linewidth = 2) +
#   geom_point(aes(x = Type, y = Percentage, color = group),
#              size = 4) +
#   theme_bw() +
#   theme(panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#         panel.border = element_blank(),
#         #axis.line.x = element_blank(),
#         axis.line = element_line(color = '#4c4c4c', linewidth = 1.5),
#         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#         #legend.position = 'none',
#         legend.title = element_blank(),
#         axis.title.x = element_blank()) +
#   scale_color_manual(values = group_color_vec)
#   
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'transloc_abund_vs_ancestry_plot.png'),
#                  plot = transloc_abund_vs_ancestry_plot,
#                  width = 10*0.75, height = 7.5*0.75, bg = 'white')
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
# ########
# ########
# ########
# 
# #data components:
# # population size each year
# 
# # contributions of each year's cohort
# 
# file_names <- list.files(here('data', 'feb2024_databasemarch2023_processed'))
# rcw_processed_list <- lapply(setNames(file_names, nm = gsub(pattern = '\\.csv', '', file_names)), function(x) {
#   read.csv(here('data', 'feb2024_databasemarch2023_processed', x))
# })
# 
# # time between cohort and current year
# transloc_details <- results_list$transloc_info_color %>% 
#   rename(id = RCWid) %>% 
#   mutate(first_contr_year = year + 1,
#          time_lag = 2022 - first_contr_year) %>% 
#   left_join(., 
#             results_list$rcw_contr_info_df %>% 
#               filter(year == 2022) %>% 
#               rename(contr_year = year), 
#             by = 'id') %>% 
#   #group_by(first_contr_year) %>% 
#   # summarize(time_lag = first(time_lag),
#   #           cohort_contr = sum(contr),
#   #           transloc_count = n()) %>% 
#   left_join(., pop_info_list$Scenario4 %>% 
#               group_by(year) %>% 
#               summarize(pop_size = n()) %>% 
#               rename(first_contr_year = year),
#             by = 'first_contr_year')
# 
# 
# 
# contr_info_processed_year <- results_list$rcw_contr_info_df %>% 
#   left_join(., results_list$transloc_info_color %>% 
#               mutate(first_contr_year = year + 1) %>% 
#               rename(id = RCWid) %>% 
#               select(id, first_contr_year, Source),
#             by = 'id') %>% 
#   filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% 
#   filter(year >= first_contr_year) %>% 
#   #mutate(years_after_contr = first_contr_year - year) %>% 
#   group_by(year, first_contr_year) %>% 
#   summarize(Source = first(Source),
#             cohort_contr = sum(contr),
#             `Translocation\ncount` = length(unique(id)),
#             .groups = 'drop')
# 
# contr_info_processed_source <- results_list$rcw_contr_info_df %>% 
#   left_join(., results_list$transloc_info_color %>% 
#               mutate(first_contr_year = year + 1) %>% 
#               rename(id = RCWid) %>% 
#               select(id, first_contr_year, Source),
#             by = 'id') %>% 
#   filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% 
#   filter(year >= first_contr_year) %>% 
#   #mutate(years_after_contr = first_contr_year - year) %>% 
#   group_by(Source) %>% 
#   mutate(year_cohort_count = length(unique(first_contr_year))) %>% 
#   ungroup() %>% 
#   group_by(year, Source) %>% 
#   summarize(cohort_contr = sum(contr),
#             year_cohort_count = first(year_cohort_count),
#             .groups = 'drop')
# 
# 
# # results_list$rcw_contr_info_df %>% 
# #   left_join(., results_list$transloc_info_color %>% 
# #               mutate(first_contr_year = year + 1) %>% 
# #               rename(id = RCWid) %>% 
# #               select(id, first_contr_year, Source),
# #             by = 'id') %>% 
# #   filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% 
# #   filter(year >= first_contr_year) %>% 
# #   filter(year == max(year)) %>% 
# #   arrange(contr)
# # 
# 
# # contr_info_processed_source %>% 
# #   filter(Source == 'ANF') %>% 
# #   ggplot() +
# #   geom_line(aes(x = year, y = cohort_contr)) +
# #   geom_line(data = contr_info_processed_year %>% 
# #               filter(Source == 'ANF'),
# #             aes(x = year,
# #                 y = cohort_contr,
# #                 color = as.character(first_contr_year)
# #                 )
# #             ) +
# #   theme_bw()
# 
# # contr_info_processed_year %>% 
# #   filter(Source == 'ANF') %>% 
# #   ggplot() + 
# #   geom_area(aes(x = year, 
# #                 y = cohort_contr, 
# #                 fill = as.character(first_contr_year))) +
# #   geom_line(data = contr_info_processed_source %>% 
# #               filter(Source == 'ANF'),
# #             mapping = aes(x = year, y = cohort_contr)
# #             )
# 
# processed_data <- contr_info_processed_year %>% 
#   mutate(`Translocation\nyear` = first_contr_year - 1) %>% 
#   group_by(first_contr_year) %>% 
#   mutate(zeroed_year = year - min(year))
# 
# background_contr_df <- lapply(setNames(nm = unique(processed_data$Source)), function(x, dat) {
#   return(
#     dat %>% 
#       mutate(Source1 = Source,
#              Source = x)
#   )
# }, dat = processed_data) %>% 
#   bind_rows()
# 
# 
# 
# 
# 
# cohort_contr_multipanel <- ggplot()  +
#   geom_hline(yintercept = 0, linewidth = 1.2) +
#   geom_line(data = background_contr_df,
#             aes(x = year, y = cohort_contr,
#                 group = `Translocation\nyear`),
#             linewidth = 0.4,
#             alpha = 0.3,
#             color = '#8c8c8c'
#             ) +
#   geom_point(data = background_contr_df,
#              aes(x = year, y = cohort_contr,
#                  group = `Translocation\nyear`),
#              size = 0.9, 
#              alpha = 0.3,
#              color = '#8c8c8c'
#              ) +
#   geom_line(data = processed_data,
#             aes(x = year, y = cohort_contr,
#                 group = `Translocation\nyear`,
#                 color = `Translocation\ncount`),
#             linewidth = 1.85) +
#   geom_point(data = processed_data,
#              aes(x = year, y = cohort_contr,
#                  group = `Translocation\nyear`,
#                  color = `Translocation\ncount`),
#              size = 2.25, alpha = 1) +
#   theme_bw() +
#   theme(
#     axis.title = element_text(size = 15),
#       strip.background.x = element_blank(),
#       strip.text.x = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.grid.major.x = element_line(linetype = 'dashed', 
#                                         color = '#d8d8d8',
#                                         linewidth = 0.3),
#       panel.border = element_blank(),
#       axis.line = element_line(color = '#4c4c4c', 
#                                linewidth = 1.5),
#     axis.line.x.top = element_blank(),
#     axis.ticks.x.top = element_blank(),
#     axis.title.x.top = element_blank()
#         ) +
#   scale_colour_gradientn(colors = c("#feb2b2", "#fe7f7f", "#fe0000", "#b10000", "#7f0000", "#650000"),
#                         breaks = c(1, 4, 7, 10), 
#                         labels = c("1", "4", "7", "10")) +
#   xlab("Year") +
#   ylab('Expected genetic contribution') +
#   facet_wrap(~Source, 
#              #scale = 'free_x', 
#              ncol = 1) +
#   geom_text(data = processed_data,
#             aes(label = Source), x = 1998.5, y = Inf, 
#             hjust = 0, vjust = 1.5, check_overlap = TRUE,
#             position = position_nudge(x = 100), size = 5,
#             fontface = 'bold') +
#   scale_y_continuous(expand = c(0, 0)) +
#   coord_cartesian(clip = 'off') +
#   scale_x_continuous(sec.axis = dup_axis())
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'cohort_contr_multipanel.png'),
#                  plot = cohort_contr_multipanel,
#                  width = 5*1, height = 7.5*1, bg = 'white')
# 
# 
# #cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'cohort_contr_multipanel.png'),
# #                 plot = cohort_contr_multipanel,
# #                 width = 8.2*1.5, height = 3.5*1.5, bg = 'white')
# 
# #cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'cohort_contr_multipanel.png'),
# #                 plot = cohort_contr_multipanel,
# #                 width = 8.2*1, height = 10.5*1, bg = 'white')
# 
# # ggplot()  +
# #   geom_line(data = processed_data,
# #             aes(x = year, y = cohort_contr,
# #                 group = `Translocation\nyear`,
# #                 color = Source),
# #             linewidth = 1.95) +
# #   geom_point(data = processed_data,
# #              aes(x = year, y = cohort_contr,
# #                  group = `Translocation\nyear`,
# #                  color = Source),
# #              size = 2.5, alpha = 1) +
# #   theme_bw() +
# #   theme(
# #     axis.title = element_text(size = 15),
# #     strip.background.x = element_blank(),
# #     strip.text.x = element_blank(),
# #     panel.grid.major = element_blank(),
# #     panel.grid.minor = element_blank(),
# #     panel.grid.major.x = element_line(linetype = 'dashed', 
# #                                       color = '#d8d8d8',
# #                                       linewidth = 0.3),
# #     panel.border = element_blank(),
# #     axis.line = element_line(color = '#4c4c4c', 
# #                              linewidth = 1.5),
# #     axis.line.x.top = element_blank(),
# #     axis.ticks.x.top = element_blank(),
# #     axis.title.x.top = element_blank()
# #   )
# 
# ########
# ########
# ########
# nests <- read_excel(here('data', 'feb2024_databasemarch2023', 'Nests.xlsx'))
# 
# tol <- 1e-12
# 
# 
# results_list$rcw_contr_info_df %>% 
#   group_by(id) %>% 
#   filter(!all(contr <= tol)) %>% 
#   ungroup() %>% 
#   filter(year == 2022) %>% 
#   filter(contr <= tol) %>% 
#   pull(id)
# 
# 
# transloc_noncontributors <- results_list$rcw_contr_info_df %>% 
#   group_by(id) %>% 
#   filter(!all(contr <= tol)) %>% 
#   ungroup() %>% 
#   filter(year == 2022) %>% 
#   filter(contr <= tol) %>% 
#   filter(id %in% results_list$transloc_info_color$RCWid) %>%
#   pull(id)
# 
# 
# nests %>%
#   filter(MaleID %in% transloc_noncontributors | FemaleID %in% transloc_noncontributors) %>% 
#   pivot_longer(cols = c(FemaleID, MaleID), names_to = 'parent_type', values_to = 'RCWid') %>%
#   select(parent_type, RCWid, NatalNest_Temp)
# 
# file_names1 <- list.files(here('data', 'feb2024_databasemarch2023_processed'))
# file_names <- file_names1[!grepl('RDS$' ,file_names1)]
# rcw_processed_list <- lapply(setNames(file_names, nm = gsub(pattern = '\\.csv', '', file_names)), function(x) {
#   read.csv(here('data', 'feb2024_databasemarch2023_processed', x))
# })
# 
# 
# survived_offspring_info <- sapply(transloc_noncontributors, function(x, PED, POP_DAT) {
#   
#   offspring_vec <- get_offspring(PED, x)
#   
#   if (is.null(offspring_vec)) return(0)
#   
#   return(sum(any(offspring_vec %in% POP_DAT$RCWid)))
#   
#   }, PED = rcw_processed_list$ped_processed, POP_DAT = pop_dat)
# 
# message(sum(survived_offspring_info), ' translocated birds without 2022 gen. contributions had offspring that were found in a census')
# 
# 
# 
# rcw_processed_list$ped_processed
# 
# 
# sapply(setNames(nm = as.character(1994:2022)), function(x, contr_info, transloc_ids) {
#   sum(
#     contr_info %>% 
#       filter(year == as.numeric(x)) %>% 
#       arrange(desc(contr)) %>% 
#       slice_head(n = 5) %>% 
#       pull(id) %in% transloc_ids
#   )
# }, contr_info = results_list$rcw_contr_info_df, transloc_ids = results_list$transloc_info_color$RCWid)
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
# contr_info_processed %>% 
#   group_by(first_contr_year) %>% 
#   mutate(zeroed_year = year - min(year)) %>% 
#   group_by()
# 
# lapply(contr_info_processed)
# 
# contr_info_processed %>% 
#   filter(id == 'ZA-AHA') %>% 
#   arrange(year) %>% 
#   mutate(percent_change = abs((contr - lag(contr))/lag(contr))*100 )
# 
# 
# unique(results_list$rcw_contr_info_df$id)
# 
# results_list$transloc_info_color %>% 
#   rename(id = RCWid) %>% 
#   mutate(first_contr_year = year + 1,
#          time_lag = 2022 - first_contr_year) %>% 
#   left_join(., 
#              results_list$rcw_contr_info_df %>% 
#                filter(year == 2022) %>% 
#                rename(contr_year = year), 
#              by = 'id') %>% 
#   group_by(first_contr_year) %>% 
#   summarize(time_lag = first(time_lag),
#             cohort_contr = sum(contr),
#             transloc_count = n()) %>% 
#   left_join(., pop_info_list$Scenario4 %>% 
#               group_by(year) %>% 
#               summarize(pop_size = n()) %>% 
#               rename(first_contr_year = year),
#             by = 'first_contr_year')
# 
# 
# 
# 
# rcw_processed_list$nests_processed
# 
# 
# 
# transloc_details %>% 
#   ggplot() +
#   geom_point(aes(x = time_lag,
#                  y = contr))
# 
# transloc_details %>% 
#   ggplot() +
#   geom_point(aes(x = pop_size,
#                  y = contr)) 
# 
# cor.test(transloc_details$contr,
#          transloc_details$pop_size)
# 
# 
# 
# summary(lm(transloc_details$contr ~ transloc_details$time_lag))
# summary(lm(transloc_details$contr ~ transloc_details$pop_size))
# 
# 
# rcw_processed_list$census_processed %>% 
#   group_by(year) %>% 
#   summarize(pop_size = n())
# 
# pop_info_list <- lapply(setNames(nm = paste0('Scenario', 1:4)), function(X, DF) {
#   DF[DF[,X,drop = TRUE],!grepl(pattern = 'Scenario', colnames(DF))]
# }, DF = rcw_processed_list$census_processed)
# 
# 
# 
# 
# 
# 
# results_list$rcw_contr_info_df %>% 
#   filter(id %in% unique(results_list$transloc_info_color$RCWid)[1:5]) %>% 
#   group_by(id) %>% 
#   filter(any(contr > 0)) %>% 
#   mutate(first_contr_year = min(year[contr > 0]),
#          biggest_contr = max(contr)) %>% 
#   filter(year >= first_contr_year) %>% 
#   arrange(year, .by_group = TRUE) %>% 
#   mutate(dif = abs(contr - lag(contr)) ) %>% 
#   filter(!is.na(dif)) %>% 
#   mutate(dif_std = (dif - mean(dif))/sd(dif)) %>% 
#   ggplot() +
#   geom_line(aes(x = year, y = dif_std, color = id)) +
#   theme(legend.position = 'none')
#   
# 
# 
# 
# results_list$rcw_contr_info_df %>% 
#   filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% 
#   group_by(id) %>% 
#   filter(any(contr > 0)) %>% 
#   mutate(first_contr_year = min(year[contr > 0]),
#          biggest_contr = max(contr)) %>% 
#   filter(year >= first_contr_year) %>% 
#   arrange(year, .by_group = TRUE) %>% 
#   mutate(dif = abs(contr - lag(contr)) ) %>% 
#   filter(!is.na(dif)) %>% 
#   filter(lag(contr) != 0) %>% 
#   mutate(dif_std = (dif - mean(dif))/sd(dif)) %>% 
#   filter( n() > 1) %>% 
#   summarise(r = cor(year, dif_std, method = 'spearman')) %>% 
#   ggplot() +
#   geom_histogram(aes(x = r), bins = 15) +
#   geom_vline(xintercept = 0, color = 'red')
# 
# 
# 
# 
# ########
# ########
# ########
# library(broom)
# library(tidyr)
# 
# #KKY-WZ' --> look into this indiv
# processed_contr_by_year <- results_list$rcw_contr_info_df %>% 
#   filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% 
#   group_by(id) %>% 
#   #filter(id == 'YYH-GZ') %>% 
#   filter(!all(contr - 0 < 1e-15)) %>% 
#   mutate(first_contr_year = min(year[contr > 0]),
#          last_contr_year = max(year[contr > 0]),
#          biggest_contr = max(contr))%>% 
#   filter(year >= first_contr_year & (year <= last_contr_year + 1)) %>% 
#   #filter(year >= first_contr_year) %>% 
#   arrange(year, .by_group = TRUE) %>% 
#   mutate(dif = abs(contr - lag(contr)) ) %>% 
#   #filter(!is.na(dif)) %>% 
#   ##filter(id == 'KKY-WZ')
#   #filter(lag(contr) != 0) %>% 
#   mutate(dif_std = (dif - mean(dif, na.rm = TRUE))/sd(dif, na.rm = TRUE)) %>% 
#   filter(n() > 2) %>% 
#   ungroup()
# 
# model_year_contr_list <- lapply(split(processed_contr_by_year, processed_contr_by_year$id), function(x) {
#   lm(contr ~ year, data = x)$coefficients[['year']]
# })
# 
# 
# processed_contr_by_year2 <- processed_contr_by_year %>% 
#   left_join(., 
#             data.frame(id = names(model_year_contr_list),
#                        slope = unlist(model_year_contr_list)),
#             by = 'id')
# 
# 
# processed_contr_by_year3 <- processed_contr_by_year2 %>% 
#   group_by(id) %>% 
#   mutate(contemp_contr = if_else((2022 %in% year) & (all(contr > 0)), 'yes', 'no')) %>% 
#   summarise(#cor = cor(year, dif, method = 'spearman', na.rm = TRUE),
#             cor = cor.test(year, dif, method = 'spearman', 
#                            na.action=na.omit, exact = FALSE)$estimate,
#             slope = first(slope),
#             year_count = n(),
#             contemp_contr = first(contemp_contr))
# 
# 
# example_indivs_data_processed %>% 
#   filter(id == 'YYH-GZ')
# 
# results_list$rcw_contr_info_df %>% 
#   filter(id == 'YYH-GZ')
# 
# example_indivs <- c('ZA-AHA',
#                     'YYH-GZ',
#                     'OHA-ZK',
#                     'KOW-BZ')
# 
# year_vs_contr_dif_cor_plot <- processed_contr_by_year3 %>% 
#   ggplot() +
#   geom_hline(yintercept = 0,
#              color = '#4c4c4c',
#              linetype = 'solid',
#              alpha = 0.7,
#              linewidth = 1.2) +
#   geom_point(aes(x = year_count, y = cor,
#                  fill = slope,
#                  size = contemp_contr),
#              color = '#808080',
#              shape = 21) +
#   theme_bw() + 
#   theme(
#     axis.title = element_text(size = 15),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#     panel.border = element_blank(),
#     axis.line = element_line(color = '#4c4c4c', linewidth = 1.5),
#     axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#     #legend.title = element_blank(),
#     legend.position = 'bottom',
#     legend.margin = margin(0, 0.75, 0, 0, unit="cm"),
#     legend.key.size = unit(1, 'lines'),
#     legend.text = element_text(margin = margin(r = 0.25, unit = 'cm'), 
#                                size = 7),
#     legend.key.width = unit(0.75, 'cm')
#     ) +
#   xlab('Number of years') +
#   ylab('Correlation of contribution change vs. year') +
#   scale_fill_gradient2(low = "#3288bd",
#                         mid = "#ffffbf",
#                         high = "#d53e4f",
#                         midpoint = 0,
#                        name = 'Contr. vs. year\nslope') +
#   scale_size_manual(values = c(2.75, 6.75), name = "Final year\ncontributor") +
#   geom_text(data = . %>% 
#               filter(id %in% example_indivs) %>% 
#               arrange(year_count) %>% 
#               mutate(label = LETTERS[1:length(example_indivs)]),
#             aes(x = year_count, y = cor, label = label),
#             fontface = 'bold',
#             size = 5,
#             hjust = -0.9, 
#             vjust = 0.5)
# 
# 
# extract_layer_info <- layer_data(year_vs_contr_dif_cor_plot, 2)
# 
# focal_indivs_info <- processed_contr_by_year3 %>% 
#   filter(id %in% example_indivs) %>% 
#   left_join(., extract_layer_info %>% 
#               rename(year_count = x,
#                      cor = y) %>% 
#               select(cor, year_count, fill), 
#             by = c('cor', 'year_count')) %>% 
#   arrange(year_count) %>% 
#   mutate(letter_label = LETTERS[1:length(example_indivs)])
# 
# example_indivs_data_processed <- results_list$rcw_contr_info_df %>% 
#   filter(id %in% example_indivs) %>% 
#   #mutate(val = rank())
#   #arrange()
#   group_by(id) %>% 
#   #filter(any(contr > 0)) %>% 
#   mutate(first_contr_year = min(year[contr > 0]),
#          last_contr_year = max(year[contr > 0]),
#          biggest_contr = max(contr)) %>% 
#   filter(year >= first_contr_year & (year <= last_contr_year + 1)) %>% 
#   left_join(., focal_indivs_info, by = 'id')
# 
# 
# example_trajectory_multipanel <- example_indivs_data_processed %>% 
#   ggplot() +
#   geom_hline(yintercept = 0, color = '#4c4c4c', linewidth = 1.5) +
#   geom_line(aes(x = year, y = contr),
#             color = '#808080',
#             linewidth = 0.75) +
#   geom_point(aes(x = year, y = contr),
#              fill = 'white',
#              color = '#808080',
#              shape = 21,
#              size = 2.5
#   ) +
#   scale_color_manual(values = setNames(focal_indivs_info$fill,
#                                        nm = focal_indivs_info$letter_label)) +
#   #scale_fill_manual(values = setNames(focal_indivs_info$fill,
#   #                                     nm = focal_indivs_info$letter_label)) +
#   geom_smooth(aes(x = year, y = contr, 
#                   color = letter_label),
#               method = "lm", se = FALSE,
#               #color = '#808080',
#               linewidth = 2.5,
#               lineend = 'round') +
#   facet_wrap(~letter_label, ncol = 1) +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         legend.position = 'none') +
#   geom_text(data = example_indivs_data_processed,
#           aes(label = letter_label), x = 1999, y = Inf, 
#           hjust = 0, vjust = 1.5, check_overlap = TRUE,
#           position = position_nudge(x = 100), size = 5,
#           fontface = 'bold') +
#   theme(
#     axis.title = element_text(size = 15),
#     strip.background.x = element_blank(),
#     strip.text.x = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_line(linetype = 'dashed', 
#                                       color = '#d8d8d8',
#                                       size = 0.3),
#     panel.border = element_blank(),
#     axis.line = element_line(color = '#4c4c4c', linewidth = 1.5),
#     axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#     axis.line.x.top = element_blank(),
#     axis.ticks.x.top = element_blank(),
#     axis.title.x.top = element_blank()
#   ) +
#   scale_y_continuous(expand = c(0, 0)) +
#   xlab('Year') +
#   ylab('Expected genetic contribution') +
#   coord_cartesian(clip = 'off') +
#   scale_x_continuous(sec.axis = dup_axis())
# 
# 
# 
# year_vs_contr_dif_cor_multipanel <- plot_grid(year_vs_contr_dif_cor_plot,
#           example_trajectory_multipanel,
#           ncol = 2,
#           align = 'h',
#           axis = 'tblr',
#           rel_widths= c(0.7, 0.25))
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'year_vs_contr_dif_cor_plot.png'),
#                  plot = year_vs_contr_dif_cor_multipanel,
#                  width = 10.5*1, height = 7*1, bg = 'white')
# 
# 
# 
# 
# 
# processed_contr_by_year3 %>% 
#   arrange(year_count, cor) %>% 
#   print(n = 100)
# 
# extract_layer_info
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
# year_vs_contr_dif_cor_plot
# 
# 
# 
# extract_layer_info
# 
# test_var <- layer_data(year_vs_contr_dif_cor_plot)
# 
# 
# 
# 
# year_vs_contr_dif_cor_plot$layers[[2]]
# 
# 
# 
# cowplot::ggsave2(filename = here('figures', 'supplement', 'figures', 'year_vs_contr_dif_cor_plot.png'),
#                  plot = year_vs_contr_dif_cor_plot,
#                  width = 10.5*1, height = 7*1, bg = 'white')
# ########
# ########
# ########
# 
# 
# 
# 
# test_scale <- scale_fill_gradient2(low = "#3288bd",
#                      mid = "#ffffbf",
#                      high = "#d53e4f",
#                      midpoint = 0,
#                      name = 'Contr. vs. year\nslope')
# 
# test_scale(0.1)
# 
# test_ramp <- colorRampPalette(colors = c("red", "white", "blue"))
# 
# img <- function(obj, nam) {
#   image(1:length(obj), 1, as.matrix(1:length(obj)), col=obj, 
#         main = nam, ylab = "", xaxt = "n", yaxt = "n",  bty = "n")
# }
# 
# img(test_ramp(100), "red-white-blue")
# 
# test_ramp(c(0, 1, -1))
# 
# processed_contr_by_year2 %>% 
#   group_by(id) %>% 
#   summarise(cor = cor(year, dif, method = 'spearman'),
#             slope = first(slope),
#             year_count = n()
#   ) %>% 
#   filter(id == 'YGB-WZ')
# 
# processed_contr_by_year2 %>% 
#   filter(id == 'YGB-WZ') %>% 
#   ggplot() +
#   geom_point(aes(year, contr)) +
#   geom_line(aes(year, contr))
# 
# 
# 
# 
# 
# 
# 
# processed_contr_by_year2 %>% 
#   summarise(cor = cor(year, dif, method = 'spearman'),
#             slope = first(slope),
#             year_count = n()
#   ) %>% 
#   print(n = 40)
# #Correlation between year and absolute difference between adjacent 
# #year contributions
# 
# #interpretation: with an increasing number of years, the correlation
# #between year and magnitude of inter-year change decreases
# 
# #magnitude of change generally decreases over time and this 
# #relationship is stronger when more years are considered
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
# ### TROUBLESHOOTING ###
# 
# # for a couple individuals, ("KKY-WZ" "Z-YYY"  "ZS-GL"), they have at least one
# #year of no genetic contribution
# 
# 
# 
# 
# rcw_contr_df_filter <- lapply(split(results_list$rcw_contr_info_df, results_list$rcw_contr_info_df$id), function(x) {
#   
#   min_contr_year <- x %>% 
#     filter(contr > 0) %>% 
#     pull(year) %>% 
#     min()
#   
#   x$min_contr_year <- min_contr_year
#   
#   return(x[x$year >= min_contr_year,])
#   
#   } ) %>% 
#   bind_rows()
# 
# rcw_contr_df_filter1 <- lapply(split(rcw_contr_df_filter, rcw_contr_df_filter$id), function(x) {
#   
#   pos_years <- x %>% 
#     filter(contr > 0) %>% 
#     pull(year)
#   
#   zero_years <- x %>% 
#     filter(abs(contr - 0) < 1e-15 ) %>% 
#     pull(year)
#   
#   #positive processing
#   if (length(pos_years) > 0) {
#     x$smallest_pos_year <- min(pos_years)
#     x$biggest_pos_year <- max(pos_years)
#   } else {
#     x$smallest_pos_year <-  NA
#     x$biggest_pos_year <- NA
#   }
#   
#   #zero processing
#   if (length(zero_years) > 0) {
#     x$smallest_zero_year <- min(zero_years)
#     x$biggest_zero_year <- max(zero_years)
#   } else {
#     x$smallest_zero_year <-  NA
#     x$biggest_zero_year <- NA
#   }
#   
#   return(x)
# }) %>% 
#   bind_rows()
# 
# 
# 
# 
# 
# 
# 
# rcw_contr_df_filter2 <- rcw_contr_df_filter1 %>% 
#   group_by(id) %>% 
#   arrange(year, .by_group = TRUE) %>% 
#   mutate(lag_val = if_else(lag(contr) > 0, 'pos', 'zero'),
#          lead_val = if_else(lead(contr) > 0, 'pos', 'zero'),
#          zero = if_else(contr > 0, 'no', 'yes')) %>% 
#   filter(!is.na(lag_val) & !is.na(lead_val)) %>% 
#   ungroup()
#   
# problem_ids <- rcw_contr_df_filter2 %>% 
#   filter(zero == 'yes' & year > smallest_pos_year & year < biggest_pos_year) %>% 
#   pull(id) %>% 
#   unique()
# 
# 
# results_list$rcw_contr_info_df %>% 
#   filter(id %in% problem_ids) %>% 
#   group_by(id) %>% 
#   arrange(year, .by_group = TRUE) %>% 
#   ungroup() %>% 
#   print(n = 100)
# 
# 
# 
# pop_info_list$Scenario4 %>% 
#   filter(RCWid == 'Z-YYY')
# 
# 
# nests %>% 
#   filter(MaleID == 'Z-YYY' | FemaleID == 'Z-YYY')
# 
# 
# results_list$rcw_contr_info_df %>% 
#   filter(id == "ZP-AG")
# 
# results_list$rcw_contr_info_df %>% 
#   filter(id == "ZP-AG")
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
# # stack(table(results_list$transloc_info_color$year)) %>% 
# #   mutate(year = as.numeric(as.character(ind)),
# #          count = as.numeric(values)) %>% 
# #   select(year, count) %>% 
# #   right_join(., data.frame(year = 1994:2022), by = 'year') %>% 
# #   mutate(count = if_else(is.na(count), 0, count)) %>% 
# #   arrange(year) %>% 
# #   mutate(cumsum_transloc = cumsum(count))
# 
# 
# #################################
# ### CODE NOT CURRENTLY IN USE ###
# #################################
# 
# # ancestry_plot_lewontin <- results_list$rcw_inbr_merge %>%
# #   mutate(count = 1) %>%
# #   group_by(year) %>%
# #   arrange(anc_prop) %>%
# #   ggplot() +
# #   geom_bar(aes(x = year, y = count, fill = anc_prop),
# #            position = "stack", stat = "identity", width = 0.98) +
# #   scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = '#ededed', high = '#4c4c4c') +
# #   geom_point(shape = "l", data = results_list$transloc_info_color,
# #              aes(x = year, y = 0), size = 5.5, color = 'black',
# #              alpha = 1,
# #              position = position_jitter(width = 0.37, height = 0, seed = 5345),
# #              show.legend = FALSE) +
# #   scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
# #                                        nm = results_list$transloc_info_color$RCWid)) +
# #   theme_bw() +
# #   theme(
# #     panel.grid.major.y = element_blank(),
# #     panel.grid.minor = element_blank(),
# #     panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
# #     panel.border = element_blank(),
# #     #axis.line.x = element_blank(),
# #     #axis.ticks.x = element_blank(),
# #     axis.line = element_line(color = '#4c4c4c', linewidth = 1.5)
# #   ) +
# #   #scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(1, 0))) +
# #   theme(
# #     plot.margin=unit(c(-0.05,1, 0,1), "cm"),
# #     axis.title = element_text(size = 18),
# #     #axis.title.x = element_blank(),
# #     axis.text.x = element_text(size = 15)#,
# #   ) +
# #   ylab('Individuals') +
# #   xlab('Year') + 
# #   coord_cartesian(clip = 'off')
# # 
# # 
# # inbreeding_plot_lewontin <- results_list$rcw_inbr_merge %>% 
# #   ggplot() +
# #   geom_point(aes(x = year, y = f_ped),
# #              shape = 21, colour = '#b7b7b7', fill = '#e3e3e3', 
# #              size = 1, alpha = 0.8,
# #              position = position_jitter(w = 0.16, h = 0)) +
# #   geom_segment(data = rcw_inbr_merge_processed,
# #                aes(x = year, xend = year, 
# #                    y = perc_10, yend = perc_90),
# #                colour = '#4c4c4c', linewidth = 1.15) +
# #   geom_point(data = rcw_inbr_merge_processed,
# #              aes(x = year, y = mean),
# #              shape = 21, colour = rcw_inbr_merge_processed$change_color, 
# #              fill = 'white', 
# #              size = 2, stroke = 1.4) +
# #   theme_bw() +
# #   theme(axis.title.y = element_text(size = 18),
# #         axis.text.x = element_text(size = 15),
# #         panel.grid.major.y = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         panel.grid.major.x = element_line(linetype = 'dashed', 
# #                                           color = '#d8d8d8'),
# #         panel.border = element_blank(),
# #         axis.line.x = element_blank(),
# #         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
# #         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
# #         axis.title.x = element_blank(),
# #         axis.ticks.x = element_blank(),
# #         plot.margin = unit(c(0.1, 1, 0.2, 1), "cm"),
# #         
# #   ) +
# #   scale_x_continuous(limits = c(1993.5, 2022.5), position = "top",
# #                      expand = c(0.01, 0)) +
# #   ylab('Inbreeding') +
# #   coord_cartesian(clip = 'off')
# # 
# # 
# # lewontin_fig <- cowplot::plot_grid(inbreeding_plot_lewontin,
# #                    ancestry_plot_lewontin,
# #                    nrow = 2,
# #                    align = "v",
# #                    axis = 'tblr',
# #                    labels = c('A', 'B'),
# #                    label_x = 0.09,
# #                    label_y = c(0.89,1),
# #                    label_size = 23) +
# #   cowplot::draw_image(flying_rcw,
# #                       scale = 0.25,
# #                       x = 0.32, y = 0.35) +
# #   theme(plot.margin = unit(c(0, -1, 0, -0.7), "cm"))
# # 
# # lewontin_fig <- ggdraw(
# #   egg::ggarrange(inbreeding_plot_lewontin,
# #                  ancestry_plot_lewontin,
# #                  labels = c('A', 'B'),
# #                  label.args = list(gp = grid::gpar(font = 1, cex = 1.4)),
# #           #hjust = c(0.2, 0.9),
# #           nrow = 2,
# #           heights = c(40, 60))
# #           ) +
# #   cowplot::draw_image(flying_rcw,
# #                       scale = 0.25,
# #                       x = 0.32, y = 0.35) +
# #   theme(plot.margin = unit(c(0, -1, 0, -0.7), "cm"))
# # 
# # 
# # cowplot::ggsave2(filename = '/Users/alexlewanski/Documents/michigan_state/funding/lewontin_2024/lewontin_fig.png',
# #                  plot = lewontin_fig,
# #                  width = 9.6*1, height = 6.5*1, bg = 'white')
# 
# 
# 
# ### PLOTS FOR TEA PARTY ###
# 
# # tea_party_plot <- results_list$rcw_contr_info_df %>%
# #   group_by(id) %>%
# #   filter(any(contr != 0)) %>% #filter out the ids with all 0s
# #   filter(contr > 0 | (contr == 0 & year > max(year[contr > 0])  )) %>%
# #   left_join(., results_list$rcws_founder_info,
# #             by = 'id') %>%
# #   ggplot() +
# #   geom_line(data = . %>% 
# #               filter(group == 'non-transloc'),
# #             aes(x = year, y = contr, group = id, linewidth = `group`),
# #             color = '#e3e3e3', linewidth = 0.7) +
# #   geom_point(data = . %>%
# #                group_by(id) %>%
# #                filter(year == min(year)) %>%
# #                slice_head(n = 1)%>% 
# #                filter(group == 'non-transloc'),
# #              aes(x = year, y = contr),
# #              color = '#e3e3e3', size = 0.7) +
# #   geom_line(data = . %>% 
# #               filter(group != 'non-transloc'),
# #             aes(x = year, y = contr, group = id, 
# #                 linewidth = group, color = id),
# #             linewidth = 1.1) +
# #   geom_point(data = . %>%
# #                group_by(id) %>%
# #                filter(year == min(year)) %>%
# #                slice_head(n = 1)%>% 
# #                filter(group != 'non-transloc'),
# #              aes(x = year, y = contr, color = id),
# #              size = 1.5) +
# #   scale_color_manual(values = setNames(results_list$transloc_info_color[[paste0(i, '_color')]],
# #                                        nm = results_list$transloc_info_color$RCWid)) +
# #   theme_bw() +
# #   theme(panel.grid.major.y = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
# #         panel.border = element_blank(),
# #         #axis.line.x = element_blank(),
# #         axis.line = element_line(color = '#4c4c4c', linewidth = 1.5),
# #         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
# #         #axis.title.x = element_blank(),
# #         axis.ticks.x = element_blank(),
# #         axis.title = element_text(size = 15),
# #         plot.margin = unit(c(0.32,1, 0, 1), "cm"),
# #         legend.position = 'none',
# #         axis.title.y = element_text(margin = margin(0, 12.5, 0, 0))
# #   ) +
# #   scale_x_continuous(expand = c(0.01, 0), 
# #                      limits = c(1993.5, 2022.5),
# #                      position = "bottom") +
# #   xlab('Year') +
# #   ylab("Exp. genetic contribution")
# 
# # tea_party_plot2 <- results_list$rcw_inbr_merge %>%
# #   mutate(count = 1) %>%
# #   group_by(year) %>%
# #   arrange(anc_prop) %>%
# #   ggplot() +
# #   geom_bar(aes(x = year, y = count, fill = anc_prop),
# #            position = "stack", stat = "identity", width = 0.98) +
# #   scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = '#ededed', high = '#4c4c4c') +
# #   geom_point(shape = "l", data = results_list$transloc_info_color,
# #              aes(x = year, y = 0, color = RCWid), size = 5.5,
# #              alpha = 1,
# #              position = position_jitter(width = 0.44, height = 0, seed = 195),
# #              show.legend = FALSE) +
# #   scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
# #                                        nm = results_list$transloc_info_color$RCWid)) +
# #   theme_bw() +
# #   theme(legend.title = element_text(size = 10),
# #         panel.grid.major.y = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
# #         panel.border = element_blank(),
# #         #axis.line.x = element_blank(),
# #         axis.ticks.x = element_blank(),
# #         axis.line = element_line(color = '#4c4c4c', linewidth = 1.5)
# #   ) +
# #   #scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(1, 0))) +
# #   #theme(
# #     #plot.margin=unit(c(-0.05,1,-0.05,1), "cm"),
# #         #axis.title.y = element_blank(),
# #         #axis.title.x = element_blank(),
# #    #     axis.text.x = element_blank()#,
# #   #) +
# #   scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
# #   coord_cartesian(clip = 'off') +
# #   xlab('Year') +
# #   ylab('Individual count')
# 
# 
# # tea_party_migration_popsizeonly <- results_list$rcw_inbr_merge %>%
# #   mutate(count = 1) %>%
# #   group_by(year) %>%
# #   arrange(anc_prop) %>%
# #   ggplot() +
# #   geom_bar(aes(x = year, y = count, fill = anc_prop),
# #            position = "stack", stat = "identity", width = 0.98) +
# #   scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = '#ededed', high = '#ededed') +
# #   geom_point(shape = "l", data = results_list$transloc_info_color,
# #              aes(x = year, y = 0, color = RCWid), size = 5.5,
# #              alpha = 1,
# #              position = position_jitter(width = 0.44, height = 0, seed = 195),
# #              show.legend = FALSE) +
# #   scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
# #                                        nm = results_list$transloc_info_color$RCWid)) +
# #   theme_bw() +
# #   theme(legend.title = element_text(size = 10),
# #         panel.grid.major = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         #panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
# #         panel.border = element_blank(),
# #         #axis.line.x = element_blank(),
# #         axis.ticks.x = element_blank(),
# #         axis.line = element_line(color = '#4c4c4c', linewidth = 1.5)
# #   ) +
# #   #scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(1, 0))) +
# #   #theme(
# #   #plot.margin=unit(c(-0.05,1,-0.05,1), "cm"),
# #   #axis.title.y = element_blank(),
# #   #axis.title.x = element_blank(),
# #   #     axis.text.x = element_blank()#,
# #   #) +
# #   scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
# #   coord_cartesian(clip = 'off') +
# #   xlab('Year') +
# #   ylab('Individual count')
# 
# # tea_party_migration_timeline <- results_list$rcw_inbr_merge %>%
# #   mutate(count = 1) %>%
# #   group_by(year) %>%
# #   arrange(anc_prop) %>%
# #   ggplot() +
# #   geom_bar(aes(x = year, y = count, fill = anc_prop),
# #            position = "stack", stat = "identity", width = 0.98) +
# #   scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = 'white', high = 'white') +
# #   geom_point(shape = "l", data = results_list$transloc_info_color,
# #              aes(x = year, y = 0, color = RCWid), size = 5.5,
# #              alpha = 1,
# #              position = position_jitter(width = 0.44, height = 0, seed = 195),
# #              show.legend = FALSE) +
# #   scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
# #                                        nm = results_list$transloc_info_color$RCWid)) +
# #   theme_bw() +
# #   theme(legend.title = element_text(size = 10),
# #         panel.grid.major = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         #panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
# #         panel.border = element_blank(),
# #         #axis.line.x = element_blank(),
# #         axis.ticks.x = element_blank(),
# #         axis.line = element_line(color = '#4c4c4c', linewidth = 1.5)
# #   ) +
# #   #scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(1, 0))) +
# #   #theme(
# #   #plot.margin=unit(c(-0.05,1,-0.05,1), "cm"),
# #   #axis.title.y = element_blank(),
# #   #axis.title.x = element_blank(),
# #   #     axis.text.x = element_blank()#,
# #   #) +
# #   scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
# #   coord_cartesian(clip = 'off') +
# #   xlab('Year') +
# #   ylab('Individual count')
# 
# # 
# # ggsave(filename = '/Users/alexlewanski/Documents/michigan_state/presentations/PNAS_special_issue_tea_party_spring2024/genetic_contr.png',
# #        plot = tea_party_plot,
# #        width = 18*0.65, height = 6*0.65, bg = 'white')
# # 
# # ggsave(filename = '/Users/alexlewanski/Documents/michigan_state/presentations/PNAS_special_issue_tea_party_spring2024/population_ancestry_full.png',
# #        plot = tea_party_plot2,
# #        width = 14*0.65, height = 6*0.65, bg = 'white')
# # 
# # 
# # ggsave(filename = '/Users/alexlewanski/Documents/michigan_state/presentations/PNAS_special_issue_tea_party_spring2024/population_size.png',
# #        plot = tea_party_migration_popsizeonly,
# #        width = 14*0.65, height = 6*0.65, bg = 'white')
# # 
# # 
# # ggsave(filename = '/Users/alexlewanski/Documents/michigan_state/presentations/PNAS_special_issue_tea_party_spring2024/population_ancestry_timeline.png',
# #        plot = tea_party_migration_timeline,
# #        width = 14*0.65, height = 6*0.65, bg = 'white')
# 
# 
# # table(pop_dat_with_anc1$year)
# # 
# # pop_dat_with_anc1 %>% 
# #   group_by(year, id) %>%
# #   slice_head() %>% 
# #   group_by(year) %>% 
# #   summarize(count = n()) %>% 
# #   pull(count) %>% 
# #   max()
# 
# test <- pop_dat_with_anc1 %>% 
#   #filter(year %in% c(2015, 2016) ) %>% 
#   select(year, id, year, anc_group, anc_val, rank_ind) %>% 
#   mutate(year_fac = year) %>% 
#   rename(Individuals = rank_ind) %>% 
#   ggplot() +
#   geom_hline(yintercept = 0.5, 
#              color = c(rep('white', 6), '#d8d8d8', rep('white', 9), '#d8d8d8', rep('white', 9), '#d8d8d8', 'white', 'white'),
#              linetype = 'dashed') +
#   geom_bar(aes(y = anc_val, x = Individuals, fill = anc_group),
#            position = "stack", stat = "identity",
#            linewidth = 0, width = 1) +
#   scale_fill_manual(values = group_color_vec) +
#   scale_color_manual(values = group_color_vec) +
#   coord_flip(clip = 'off') +
#   facet_grid(~year_fac, switch = "x") +
#   theme_bw() +
#   theme(
#     #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
#     #strip.placement = "outside",
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     
#     #strip.background = element_rect(fill = NA, color = "white"),
#     strip.text = element_text(size = 7.5),
#     panel.spacing.x = unit(0,"cm"),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_blank(),
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#     axis.line.x = element_blank(),
#     axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#     legend.title = element_blank(),
#     legend.position = 'none',
#     plot.margin = unit(c(-0.01, 0.1, -0.01, 1), "cm"),
#     axis.title.y = element_text(size = 10)
#     #axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
#   ) +
#   scale_x_continuous(expand = c(0.01, 0)#,
#                      #breaks = c(50, 100, 150)
#                      ) +
#   scale_y_continuous(lim = c(-0.1, 1.1), expand = c(0, 0)) +
#   ylab('Individuals') +
#   geom_text(data = results_list$transloc_info_color %>% 
#               group_by(year) %>%
#               summarize(Source = first(Source),
#                         count = n(),
#                         .groups = 'drop') %>% 
#               mutate(year_fac = year) %>% 
#               left_join(.,
#                         pop_dat %>% 
#                           group_by(year) %>% 
#                           summarize(pop_count = n()),
#                         by = 'year'),
#             aes(x = pop_count + 13, y = 0.5, label = count, color = Source),
#             fontface = 'bold',
#             size = 4.75)
# 
# 
# 
# # bottom_bar_test <- results_list$rcw_inbr_merge %>%
# #   mutate(year_fac = year) %>%
# #   mutate(count = 1) %>%
# #   group_by(year) %>%
# #   arrange(anc_prop) %>%
# #   ggplot() +
# #   geom_vline(xintercept = 'a',
# #              color = c(rep('white', 6), '#d8d8d8', rep('white', 9), '#d8d8d8', rep('white', 9), '#d8d8d8', 'white', 'white'),
# #              linetype = 'dashed') +
# #   geom_bar(aes(x = 'a', y = count, fill = anc_prop),
# #            position = "stack", stat = "identity", 
# #            linewidth = 0, width = 1) +
# #   #facet_grid(~year_fac, switch = "x") +
# #   facet_grid(~year_fac) +
# #   scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = '#ededed', high = '#4c4c4c') +
# #   # geom_point(shape = "l", data = results_list$transloc_info_color,
# #   #            aes(x = year, y = 0, color = RCWid), size = 5.5,
# #   #            alpha = 1,
# #   #            position = position_jitter(width = 0.44, height = 0, seed = 195),
# #   #            show.legend = FALSE) +
# #   # scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
# #   #                                      nm = results_list$transloc_info_color$RCWid)) +
# #   theme_bw() +
# #   theme(
# #     #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
# #     #strip.placement = "outside",
# #     strip.background = element_blank(),
# #     strip.text.x = element_blank(),
# #     #strip.background = element_rect(fill = NA, color = "white"),
# #     strip.text = element_text(size = 7.5),
# #     panel.spacing.x = unit(0,"cm"),
# #     axis.title = element_blank(),
# #     axis.ticks.x = element_blank(),
# #     axis.text.x = element_blank(),
# #     panel.border = element_blank(),
# #     panel.grid.major = element_blank(),
# #     panel.grid.minor = element_blank(),
# #     axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
# #     #axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.5),
# #     axis.line.x = element_blank(),
# #     axis.ticks = element_line(color = '#808080', linewidth = 0.7),
# #     #legend.title = element_blank(),
# #     #legend.position = 'none',
# #     plot.margin = unit(c(-0.05, 0.1, -0.3, 1), "cm"),
# #     #axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
# #   ) +
# #   # theme(legend.title = element_text(size = 10),
# #   #       panel.grid.major.y = element_blank(),
# #   #       panel.grid.minor = element_blank(),
# #   #       panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
# #   #       panel.border = element_blank(),
# #   #       axis.line.x = element_blank(),
# #   #       axis.ticks.x = element_blank(),
# #   #       axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5)
# #   # ) +
# #   scale_y_reverse(expand = expansion(mult = c(0, 0), 
# #                                      add = c(1, 0))) +
# #   theme(plot.margin=unit(c(-0.01, 1, -0.01,1), "cm"),
# #         axis.title.y = element_blank(),
# #         axis.title.x = element_blank(),
# #         axis.text.x = element_blank()#,
# #   ) +
# #   #scale_y_continuous(lim = c(-0.1, 1.1), expand = c(0, 0)) +
# #   #scale_x_continuous(expand = c(0, 0)#, 
# #   #                   #limits = c(1993.5, 2022.5)
# #   #                   ) +
# #   coord_cartesian(clip = 'off')
# # 
# 
# # line_plot <- data.frame(x = 1994:2022,
# #                         y = rnorm(29, 0, 1)) %>% 
# #   ggplot() +
# #   geom_line(aes(x = x, y = y)) +
# #   geom_point(aes(x = x, y = y)) +
# #     theme_bw()
# 
# 
# 
# 
# set.seed(32)
# 
# test_1 <- test +
#   geom_point(data = results_list$transloc_info_color %>%
#                mutate(rand_val = runif(54, 0.01, 0.99),
#                year_fac = year,
#                       x = 'a'), 
#              shape = "l", 
#              aes(x = 2, y = rand_val, color = RCWid), size = 5.5,
#              alpha = 1,
#              #position = position_jitter(width = 0, height = 0, seed = 1915),
#              show.legend = FALSE) +
#   scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
#                                        nm = results_list$transloc_info_color$RCWid)
#                      ) +
#   theme(axis.line.x = element_blank())
# 
# 
#   geom_jitter(data = data.frame(year_fac = c(2022, 2022),
#                         x = c('a', 'a)')),
#              aes(x = x, y = 0))
# 
# 
# bottom_bar_test +
#   geom_point(data = data.frame(year_fac = factor(c(2022)),
#                                x = 'a',
#                                y = 1),
#              aes(x = x, y = 0))
# 
# 
# ggarrange(fped_plot_facet,
#           test, #+
#                       # geom_text(data = results_list$transloc_info_color %>%
#                       #             group_by(year) %>%
#                       #             summarize(count = n(),
#                       #                       .groups = 'drop') %>%
#                       #             mutate(year_fac = year),
#                       #           aes(x = 4, y = 0.5, label = count),
#                       #           fontface = 'bold'
#                       #),
#                     bottom_bar_pbg,
#           #heights = c(1 - (47/172), (47/172)),
#           heights = c(0.1, 0.5, 0.5),
#                     nrow = 3)
# 
# test_pan <- plot_grid(test +
#             geom_text(data = results_list$transloc_info_color %>% 
#                         group_by(year) %>% 
#                         summarize(count = n(),
#                                   .groups = 'drop') %>% 
#                         mutate(year_fac = year),
#                       aes(x = 4, y = 0.5, label = count),
#                       fontface = 'bold'
#             ), 
#             bottom_bar_pbg,
#           #bottom_bar_pbg,
#           
#           nrow = 2,
#           align = 'hv',
#           axis = 'tblr')
# 
# 
# test_pan + 
#   theme(plot.margin=unit(c(0, 0, 0.5, 0), "cm")) +
#   annotate("text", x = c(0.255, 0.47, 0.7), y = -0.02, label = c("2000", "2010", "2020"),
#            size = 3.5)
# 
# 
# 
#   geom_point(data = results_list$transloc_info_color %>%
#                mutate(rand_val = runif(54, 0.01, 0.99),
#                       year_fac = year,
#                       x = 'a'), 
#              shape = "l", 
#              aes(x = 2, y = rand_val, color = RCWid), size = 5.5,
#              alpha = 1,
#              #position = position_jitter(width = 0, height = 0, seed = 1915),
#              show.legend = FALSE)
# 
# 
# 
# 
# 
# 
# egg::ggarrange(fped_plot,
#                top_bar,
#                bottom_bar,
#                transloc_ancestry,
#                nrow = 4,
#                labels = c('A', 'B', '', 'C'),
#                label.args = list(gp = grid::gpar(font = 2, cex = 1.2)),
#                heights = c(30, 50, 50, 22))




# lastseendetail <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'LastSeenDetail.xlsx'))
# clusterdetail <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'ClusterDetail.xlsx'))
# socialclassdetail <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'SocialClassDetail.xlsx'))
# nests <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'Nests.xlsx'))
# rcws <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'RCWs.xlsx')) 

# clusterdetail %>% 
#   mutate(year = get_year(RecordDate),
#          month = get_month(RecordDate)) %>% 
#   group_by(RCWid, year) %>% 
#   filter(month == max(month)) %>% 
#   rename(maxrecorddate = RecordDate) %>% 
#   select(RCWid, maxrecorddate) %>% 
#   left_join(., lastseendetail %>% 
#               filter(SurveyType == 'Census'), 
#             by = c('RCWid'))
# 
# 
# testvec <- socialclassdetail %>% 
#   group_by(RCWid, get_year(RecordDate)) %>% 
#   filter(get_month(RecordDate) >= 6 & get_month(RecordDate) <= 8) %>% 
#   filter(get_month(RecordDate) == max(get_month(RecordDate))) %>% 
#   summarize(count_per_year = n()) %>% 
#   pull(RCWid)
# 
# socialclassdetail$id_year <- paste0(socialclassdetail$RCWid, get_year(socialclassdetail$RecordDate) )
# socialclassdetail$year <- get_year(socialclassdetail$RecordDate)
# socialclassdetail$month <- socialclassdetail$month <- get_month(socialclassdetail$RecordDate)

# socialclassdetail_filter <- lapply(split(socialclassdetail, socialclassdetail$id_year), function(x) {
#   if (any(x$month >= 6 & x$month <= 8)) {
#     subset_df <- x[x$month >= 6 & x$month <= 8,]
#   } else {
#     subset_df <- x
#   }
#   return(subset_df[subset_df$month == max(subset_df$month),][1,])
# }) %>% 
#   bind_rows()


#Z-AOY, HKP-Z, BWG-Z, GZ-OOK

# lastseen_socialclass_join <- lastseendetail %>% 
#   mutate(year = get_year(RecordDate)) %>% 
#   filter(SurveyType == 'Census' & Detected == 'Yes' & (get_month(RecordDate) >= 6 & get_month(RecordDate) <= 8 ) ) %>% 
#   left_join(., socialclassdetail_filter, by = c('RCWid', 'year')) %>% 
#   select(RCWid, RecordDate.x, RecordDate.y, year, SocialClass) %>%
#   filter(RCWid == 'GZ-OOK')

#library(clock)




# aggregationyear <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'AggregationYear.xlsx')) 
# censusdata <- readxl::read_xlsx(here('data', 'feb2024_databasemarch2023', 'CensusData.xlsx'))
# 


#census_cluster_subset <- clusterdetail[clusterdetail$MoveType %in% 'No move_Census' & get_month(clusterdetail$RecordDate) %in% c(6, 7, 8),]



#& identical(TRUE, AgeKnown == 'Y') )

# clusterdetail %>% 
#   filter(RecordDate %in% censusdata$RecordDate[censusdata$SurveySeason == 'Post-season']) %>% 
#   left_join(rcws[,c('RCWid', 'Sex', 'MinAge', 'AgeKnown')], by = 'RCWid') %>% 
#   mutate(census_year = clock::get_year(RecordDate)) %>% 
#   filter( MinAge < census_year  | is.na(MinAge) | AgeKnown != 'Y') %>% 
#   group_by(ClusterID, census_year) %>% 
#   summarize(both_sex = if_else(all(c('M', 'F') %in% Sex), 1, 0),
#             .groups = 'drop') %>% 
#   #filter(!grepl("^RR", ClusterID)) %>% 
#   group_by(census_year) %>% 
#   summarize(pbg = sum(both_sex)) %>% 
#   print(n = 30)

# 
# 
# clusterdetail %>% 
#   #filter(RecordDate %in% censusdata$RecordDate[censusdata$SurveySeason == 'Post-season']) %>% 
#   left_join(rcws[,c('RCWid', 'Sex', 'MinAge')], by = 'RCWid') %>% 
#   mutate(census_year = get_year(RecordDate)) %>% 
#   filter(census_year == 1997) %>% 
#   filter(grepl("^AP", ClusterID))
# #group_by(ClusterID) %>% 
# #summarize(indiv_count = n()) %>% 
# #print(n = 100)
# 


# clusterdetail %>% 
#   filter(get_year(RecordDate) == 1997) %>% 
#   left_join(., 
#             rcws %>% select(RCWid, MinAge, Sex), 
#             by = 'RCWid') %>% 
#   group_by(ClusterID) %>% 
#   summarize()
# 
# nests %>% 
#   filter(Year == 1999) %>% 
#   print(n = 100)
# 
# nests %>% 
#   group_by(Year, Cluster) %>% 
#   slice_head(n = 1) %>% 
#   group_by(Year) %>% 
#   summarize(realized_breeding_groups = n()) %>% 
#   print(n = 40)


# lastseendetail_subset <- lastseendetail %>% 
#   filter(Detected == 'Yes' & SurveyType == 'Census' & clock::get_month(lastseendetail$RecordDate) %in% c(6, 7, 8))
# 
# lastseendetail_subset$most_recent_date <- lubridate::mdy_hms(NA)
# lastseendetail_subset$ClusterID <- NA
# for (i in seq_len(nrow(lastseendetail_subset))) {
#   
#   focal_date <- lastseendetail_subset$RecordDate[i]
#   if (is.na(focal_date)) next
#   if (!lastseendetail_subset$RCWid[i] %in% clusterdetail$RCWid) stop('STOP!')
#   dates_vec <- clusterdetail$RecordDate[clusterdetail$RCWid == lastseendetail_subset$RCWid[i]]
#   
#   if (any(dates_vec <= focal_date)) {
#     lastseendetail_subset$most_recent_date[i] <- max(dates_vec[dates_vec <= focal_date])[1]
#   } else {
#     lastseendetail_subset$most_recent_date[i] <- min(dates_vec)[1]
#   }
#   
#   lastseendetail_subset$ClusterID[i] <- clusterdetail$ClusterID[clusterdetail$RCWid == lastseendetail_subset$RCWid[i] & clusterdetail$RecordDate == lastseendetail_subset$most_recent_date[i]][1]
# }
# 
# 
# 
# 
# 
# lastseendetail_subset2 <- lastseendetail_subset %>% 
#   #filter(RecordDate %in% censusdata$RecordDate[censusdata$SurveySeason == 'Post-season']) %>% 
#   left_join(rcws[,c('RCWid', 'Sex', 'MinAge', 'AgeKnown')], by = 'RCWid') %>% 
#   mutate(census_year = clock::get_year(RecordDate)) %>% 
#   filter(MinAge < census_year | is.na(MinAge)) %>% 
#   mutate(ClusterID_Year = paste0(ClusterID, "_", census_year)) %>% 
#   group_by(ClusterID_Year) %>% 
#   mutate(both_sex = if_else(all(c('M', 'F') %in% Sex), 1, 0)) %>% 
#   ungroup() %>% 
#   mutate(transloc = if_else(RCWid %in% results_list$transloc_info_color$RCWid, 'yes', 'no') )
# 
# 
# pbg_summary <- lastseendetail_subset2 %>% 
#   group_by(ClusterID_Year) %>%
#   summarize(ClusterID = first(ClusterID),
#             year = first(census_year),
#             both_sex = if_else(all(c('M', 'F') %in% Sex), 1, 0),
#             `Group\ncomposition` = case_when(all(transloc == 'yes') ~ 'only transloc.',
#                                    all(transloc == 'no') ~ 'no transloc.',
#                                    TRUE ~ 'combined'),
#             .groups = 'drop') %>%
#   filter(both_sex == 1) #%>%
#   ##filter(!grepl('^RR', ClusterID)) %>% 
#   #group_by(year) %>% 
#   #summarize(pbg_count = n(),
#   #          .groups = 'drop')
# 
# 
# 
# bottom_bar_pbg <- pbg_summary %>% 
#   mutate(year_fac = year,
#          `Group\ncomposition` = factor(`Group\ncomposition`, 
#                                        levels = rev(c('no transloc.', 'combined', 'only transloc.')))
#          ) %>% 
#   ggplot() +
#   geom_vline(xintercept = 'a',
#              color = c(rep('white', 6), '#d8d8d8', rep('white', 9), '#d8d8d8', rep('white', 9), '#d8d8d8', 'white', 'white'),
#              linetype = 'dashed') +
#   #geom_bar(aes(x = 'a', y = count, fill = anc_prop),
#   #         position = "stack", stat = "identity", 
#   #         linewidth = 0, width = 1) +
#   geom_bar(aes(x = 'a', fill = `Group\ncomposition`),
#            position = "stack", stat = 'count', width = 0.98) +
#   scale_fill_manual(values = c("#6a6a6a", '#bdbdbd', "#ededed")) + #c('#333333', "#999999", "#ededed")) +
#   #facet_grid(~year_fac, switch = "x") +
#   facet_grid(~year_fac) +
#   theme_bw() +
#   ylab('Potential breeding groups') +
#   theme(
#     legend.title = element_text(size = 9.75),
#     #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
#     #strip.placement = "outside",
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     #strip.background = element_rect(fill = NA, color = "white"),
#     strip.text = element_text(size = 7.5),
#     panel.spacing.x = unit(0,"cm"),
#     #axis.title = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_blank(),
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#     #axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.5),
#     axis.line.x = element_blank(),
#     axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#     #legend.title = element_blank(),
#     #legend.position = 'none',
#     plot.margin = unit(c(-0.05, 0.1, -0.3, 1), "cm"),
#     axis.title.y = element_text(size = 10, margin = margin(0, 2.5, 0, 0, unit = "mm"))
#     #axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
#   ) +
#   #scale_y_reverse(expand = expansion(mult = c(0, 0), 
#   #                                     add = c(1, 0))) +
#   theme(plot.margin=unit(c(-0.01, 1, -0.01,1), "cm"),
#         #axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank()#,
#   ) +
#   scale_y_continuous(expand = c(0, 0)) +
#   #scale_x_continuous(expand = c(0, 0)#, 
#   #                   #limits = c(1993.5, 2022.5)
#   #                   ) +
#   coord_cartesian(clip = 'off')
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
# rcw_inbr_merge_processed <- results_list$rcw_inbr_merge %>% 
#   group_by(year) %>% 
#   summarize(perc_10 = quantile(f_ped, prob = 0.1),
#             mean = mean(f_ped),
#             perc_90 = quantile(f_ped, prob = 0.9),
#             .groups = 'drop') %>% 
#   arrange(year) %>% 
#   mutate(change_color = case_when(mean > lag(mean) ~ "#f0597d",
#                                   mean < lag(mean) ~ "#118ab2",
#                                   TRUE ~ '#4c4c4c'))
# 
# fped_plot_facet <- results_list$rcw_inbr_merge %>% 
#   ggplot() +
#   geom_vline(xintercept = 'a',
#              color = c(rep('white', 6), '#d8d8d8', rep('white', 9), '#d8d8d8', rep('white', 9), '#d8d8d8', 'white', 'white'),
#              linetype = 'dashed') +
#   geom_point(aes(x = 'a', y = f_ped),
#              shape = 21, colour = '#b7b7b7', fill = '#e3e3e3', 
#              size = 1, alpha = 0.8,
#              position = position_jitter(w = 0.16, h = 0)) +
#   geom_segment(data = rcw_inbr_merge_processed,
#                aes(x = 'a', xend = 'a', 
#                    y = perc_10, yend = perc_90),
#                colour = '#4c4c4c', linewidth = 1.15) +
#   geom_point(data = rcw_inbr_merge_processed,
#              aes(x = 'a', y = mean),
#              shape = 21, color = rcw_inbr_merge_processed$change_color, 
#              fill = 'white', 
#              size = 2.25, stroke = 1.5) +
#   theme_bw() +
#   facet_wrap(~year,
#              nrow = 1) +
#   theme(
#     #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
#     #strip.placement = "outside",
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     #strip.background = element_rect(fill = NA, color = "white"),
#     strip.text = element_text(size = 7.5),
#     panel.spacing.x = unit(0,"cm"),
#     #axis.title = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_blank(),
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#     #axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.5),
#     axis.line.x = element_blank(),
#     axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#     #legend.title = element_blank(),
#     #legend.position = 'none',
#     plot.margin = unit(c(-0.05, 0.1, -0.3, 1), "cm"),
#     axis.title.y = element_text(size = 10)
#     #axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
#   ) +
#   #scale_y_reverse(expand = expansion(mult = c(0, 0), 
#   #                                   add = c(1, 0))) +
#   theme(plot.margin=unit(c(-0.01, 1, -0.01,1), "cm"),
#         #axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank()#,
#   ) +
#   #scale_y_continuous(lim = c(-0.1, 1.1), expand = c(0, 0)) +
#   #scale_x_continuous(expand = c(0, 0)#, 
#   #                   #limits = c(1993.5, 2022.5)
#   #                   ) +
#   ylab(expression(italic(F)["P"])) +
#   coord_cartesian(clip = 'off')
#   
# 
#  
# transloc_ancestry_facet <- results_list$rcws_ancestry_source_by_year %>%
#   mutate(group = factor(group, 
#                         levels = rev(c('non-transloc', 'ANF', 'CBJTC', 'FTB', 'FTS', 'ONF', 'WSF-CITRUS')))) %>% 
#   ggplot() +
#   geom_hline(yintercept = 0.5, 
#              color = c(rep('white', 6), '#d8d8d8', rep('white', 9), '#d8d8d8', rep('white', 9), '#d8d8d8', 'white', 'white'),
#              linetype = 'dashed') +
#   geom_bar(aes(x = 'a', y = ancestry_prop, fill = group), 
#            color = 'white', linewidth = 0.15, width = 0.72,
#            position = "stack", stat = "identity") + 
#   #scale_fill_manual(values = c('#ededed', '#4c4c4c')) +
#   scale_fill_manual(values = group_color_vec) +
#   facet_grid(~year) +
#   theme_bw() +
#   theme(
#     #panel.spacing = unit(c(20, rep(0, 27)), "lines"),
#     #strip.placement = "outside",
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     #strip.background = element_rect(fill = NA, color = "white"),
#     strip.text = element_text(size = 7.5),
#     panel.spacing.x = unit(0,"cm"),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_blank(),
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#     axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.5),
#     #axis.line.x = element_blank(),
#     axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#     #legend.title = element_blank(),
#     legend.position = 'none',
#     plot.margin = unit(c(-0.05, 0.1, -0.1, 1), "cm"),
#     axis.title.y = element_text(size = 10, margin = margin(0, 1.8, 0, 0, unit = "mm"))
#     #axis.title.y = element_text(margin = margin(0, 8.25, 0, 0))
#   ) +
#   scale_y_continuous(breaks = c(0, 0.5, 1), 
#                      labels = c("0.0", "0.5", "1.0"),
#                      expand = expansion(c(0, 0))) +
#   xlab('Year') +
#   ylab("Ancestry") +
#   theme(legend.key.size = unit(0.4, 'cm')) +
#   coord_cartesian(clip = 'off')
# 
# ### custom legends ###
# fped_plot_legend <- ggdraw(cowplot::get_legend(
#   data.frame(x = 1,
#              y = 3:1,
#              label = factor(c('increase', 'no change', 'decrease'), 
#                             levels = c('increase', 'no change', 'decrease')) ) %>% 
#     ggplot() +
#     geom_point(aes(x = x,
#                    y = y,
#                    color = label),
#                shape = 21,
#                fill = 'white',
#                size = 2.1, stroke = 1.5) +
#     theme_bw() +
#     theme(legend.key = element_blank(),
#           legend.title = element_blank(),
#           legend.spacing.y = unit(-1.5, 'mm'),
#           legend.background = element_rect(fill = "transparent"),
#           panel.background = element_rect(fill = "transparent", color = "transparent")) +
#     scale_color_manual(substitute(paste('Change in mean(', italic(F)["p"], ")")),
#                        values = c("#f0597d", '#4c4c4c', "#118ab2")) +
#     guides(color = guide_legend(byrow = TRUE)) #required to adjust vertical spacing of legend
# 
# )) +
#   theme(plot.margin=unit(c(0,-5, 0, -5), "cm"),
#         panel.background = element_rect(fill = "transparent", color = "transparent"))



# founder_prop_legend <- results_list$rcws_founder_info %>%
#   group_by(group) %>% 
#   summarize(count = n(),
#             .groups = 'drop') %>%
#   mutate(group = factor(group, 
#                         levels = rev(c('CBJTC', 'ANF', 'FTB', 'WSF-CITRUS', 'FTS', 'ONF', 'non-transloc')))) %>% 
#   mutate(x = 0.95) %>%
#   arrange(desc(group)) %>% 
#   mutate(cumsum_count = cumsum(count) - (0.5*count)) %>% 
#   mutate(cumsum_count_adjust = cumsum_count + c(-13, -5, 0, 10, 17, 24, 0)) %>% 
#   mutate(update_label = paste0(group, ' (', count, ')')) %>% 
#   ggplot(aes(x = x, y = count, fill = group, label = count)) +
#   geom_segment(aes(x = 1.47,
#                    xend = 2,
#                    y = cumsum_count, 
#                    yend = cumsum_count_adjust,
#                    color = group),
#                linewidth = 1.1) +
#   geom_bar(position = 'stack', stat = "identity", color = 'white', linewidth = 0.15) +
#   scale_color_manual(values = group_color_vec) +
#   scale_fill_manual(values = group_color_vec) +
#   geom_text(aes(x = 2.1, y = cumsum_count_adjust, label = update_label),
#             size = 2.6,
#             hjust = 0) +
#   theme_void() +
#   theme(legend.position = 'none') +
#   coord_cartesian(clip = 'off') +
#   xlim(-0.32, 4) +
#   geom_text(aes(x = 0.5, y = 220, label = "Founder\nproportions"),
#             size = 3.5,
#             hjust = 0,
#             lineheight = .9)

### add in legends and RCW image ###
# rcw_pop_overview <- ggpubr::annotate_figure(egg::ggarrange(fped_plot,
#                                                            top_bar,
#                                                            bottom_bar,
#                                                            transloc_ancestry,
#                                                            nrow = 4,
#                                                            labels = c('A', 'B', '', 'C'),
#                                                            label.args = list(gp = grid::gpar(font = 2, cex = 1.2)),
#                                                            heights = c(30, 50, 50, 22)),
#                                             left = grid::textGrob("Individuals", rot = 90, hjust = 0.58, vjust = 5.2, gp = grid::gpar(cex = 1))) +
#   theme(plot.margin = unit(c(0.3, -0.8, 0, -3), "cm"))


# rcw_pop_overview_facet <- egg::ggarrange(fped_plot_facet +
#                                            theme(plot.margin = unit(c(0, 0, 0, 0.1), "cm")),
#                                          bottom_bar_pbg +
#                                            theme(plot.margin = unit(c(0.2, 0, 0, 0.1), "cm")),
#                                          test +
#                                            theme(plot.margin = unit(c(0.2, 0, 0, 0.1), "cm")),
#                                          transloc_ancestry_facet +
#                                            theme(plot.margin = unit(c(0.2, 0, 0, 0.1), "cm")),
#                                          nrow = 4,
#                                          labels = c('A', 'B', 'C', 'D'),
#                                          label.args = list(gp = grid::gpar(font = 2, cex = 1.2, hjust = -0.1)),
#                                          heights = c(30, 50, 50, 22))
# 
# 
# test_multipan <- ggdraw(rcw_pop_overview_facet) +
#   theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.1), "cm")) +
#   annotate("text", x = c(0.23, 0.521, 0.814), y = -0.011, label = c("2000", "2010", "2020"),
#            size = 3.25) +
#   annotate("text", x = c(0.23, 0.521, 0.814), y = 1.015, label = c("2000", "2010", "2020"),
#            size = 3.25)
# 
# rcw_pop_overview_rcw <- test_multipan +
#   cowplot::draw_image(flying_rcw,
#                       scale = 0.12,
#                       x = 0.4, y = 0.46) +
#   patchwork::inset_element(founder_prop_legend,
#                            left = 0.918 - 0.033, bottom = 0.18, right = 1.015 - 0.033, top = 0.42
#                            #left = 0.849 - 0.0358, bottom = 0.5, right = 0.98 - 0.0358, top = 0.70
#                            #left = 0.849, bottom = -0.01, right = 0.972, top = 0.20
#   ) +
#   patchwork::inset_element(fped_plot_legend,
#                            #left = 0, bottom = 0.9, right = 0.38, top = 0.95,
#                            #left = 0.891 - 0.0358, bottom = 0.81, right = 0.94 - 0.0358, top = 0.88
#                            #left = 0.891 - 0.032, bottom = 0.85, right = 0.94 - 0.032, top = 0.86
#                            left = 0.965 - 0.032, bottom = 0.83, right = 0.99 - 0.032, top = 0.86
#                            )
#   
# 
# shared_axis_title <- textGrob("Year", 
#                               gp = gpar(fontface = "plain", col = "black", fontsize=15),
#                               x = unit(0.445, "npc"),
#                               y = unit(0.8, "npc"))
# 
# #add to plot
# 
# rcw_pop_overview_rcw_title <- grid.arrange(arrangeGrob(ggdraw(rcw_pop_overview_rcw), 
#                                                        bottom = shared_axis_title))
# 
# cowplot::ggsave2(filename = here('figures', 'main_paper', 'rcw_pop_overview_fig.png'),
#                  #plot = ggdraw(rcw_pop_overview_rcw) + theme(plot.margin = unit(c(-0.3, 0.3, 0.25, -0.6), "cm")),
#                  plot = ggdraw(rcw_pop_overview_rcw_title) +
#                    theme(plot.margin = unit(c(-0.22, -0.2, -0.2, 0), "cm")), #ggdraw(rcw_pop_overview_rcw) +
#                  #theme(plot.margin = unit(c(0, -0.9, 0, 0), "cm")),
#                  #width = 10*1, height = 8.1*1, 
#                  width = 10*1.2, height = 5.5*1.2, 
#                  bg = 'white')
# 
# 
# # cowplot::ggsave2(filename = here('figures', 'main_paper', 'rcw_pop_overview_fig_TEST.png'),
# #                  #plot = ggdraw(rcw_pop_overview_rcw) + theme(plot.margin = unit(c(-0.3, 0.3, 0.25, -0.6), "cm")),
# #                  plot = rcw_pop_overview_rcw_final, #ggdraw(rcw_pop_overview_rcw) +
# #                    #theme(plot.margin = unit(c(0, -0.9, 0, 0), "cm")),
# #                  width = 10*1, height = 8.1*1, bg = 'white')
# 
# #  annotate("text", x = c(0.255, 0.47, 0.7), y = -0.02, label = c("2000", "2010", "2020"),
# #           size = 3.5)
# 
# 
# 
# 
# 
# anc_mixing_df %>% 
#   filter(!is.nan(md)) %>% 
#   mutate(md_update = if_else(md < 0, 0, md)) %>% 
#   arrange(year) %>% 
#   mutate(md_update_dif = md_update - lag(md_update)) %>% 
#   filter(!is.na(md_update_dif)) %>%
#   summarize(mean_md_dif = round(mean(md_update_dif), 3),
#             sd_md_dif = round(sd(md_update_dif), 3)) %>% 
#   nrow()
# 
# round(anc_mixing_df$md[anc_mixing_df$year %in% 2022], 3)
# 
# table(results_list$rcw_inbr_merge$group)
# 
# 
# 
# 
# pop_dat %>% 
#   filter(year == 2022) %>% 
#   nrow()
# 
# ### Summaries reported in paper ###
# 
# pop_dat_with_anc <- pop_dat %>% 
#   rename(id = RCWid) %>% 
#   left_join(., results_list$rcws_ancestry_info,
#             by = 'id',
#             relationship = "many-to-many")
# 
# 
# pop_dat_with_anc %>% 
#   filter(year == 2022) %>% 
#   filter(anc_group != 'non-transloc') %>% 
#   pull(id) %>% 
#   unique() %>% 
#   length()
# 
# results_list$rcw_inbr_merge %>% 
#   filter(year == 2022 & anc_prop > 0) %>%
#   summarize(mean_transloc_anc = round(mean(anc_prop*100), 1),
#             sd_transloc_anc = round(sd(anc_prop*100), 1) )
# 
# results_list$rcw_inbr_merge %>% 
#   filter(year == 2022 & anc_prop > 0) %>%
#   summarize(mean_transloc_anc = round(mean(anc_prop), 2),
#             sd_transloc_anc = round(sd(anc_prop), 2) )
# 
# #pop-level ancestry 
# results_list$rcws_ancestry_source_by_year %>% 
#   filter(year == 2022) %>% 
#   filter(group != "non-transloc") %>%
#   pull(ancestry_prop) %>% 
#   sum()
# 
# results_list$rcws_ancestry_source_by_year %>% 
#   filter(year == 2022) %>% 
#   filter(group != "non-transloc") %>%
#   mutate(prop_transloc = 100*(ancestry_prop/sum(ancestry_prop)))
# 
# #inbreeding
# 
# first_inbr_year <- results_list$rcw_inbr_merge %>% 
#   filter(f_ped > 0) %>% 
#   pull(year) %>% 
#   min()
# 
# 
# results_list$rcw_inbr_merge %>% 
#   filter(year >= first_inbr_year) %>% 
#   group_by(year) %>% 
#   summarize(mean_fped = mean(f_ped)) %>% 
#   arrange(year) %>% 
#   mutate(average_fped_dif = mean_fped - lag(mean_fped)) %>% 
#   filter(!is.na(average_fped_dif)) %>% 
#   summarize(mean_average_fped_change = mean(average_fped_dif),
#             sd_average_fped_change = sd(average_fped_dif))
# 
# 
# results_list$rcw_inbr_merge %>% 
#   filter(year >= first_inbr_year) %>% 
#   group_by(year) %>% 
#   summarize(mean_fped = mean(f_ped)) %>% 
#   filter(mean_fped == max(mean_fped))
# 
# 
# results_list$rcw_inbr_merge %>% 
#   filter(year >= first_inbr_year) %>% 
#   group_by(year) %>% 
#   summarize(non_zero_fped_perc = 100*(sum(f_ped > 0)/n())) %>% 
#   arrange(year) %>% 
#   mutate(non_zero_fped_perc_change = non_zero_fped_perc - lag(non_zero_fped_perc)) %>% 
#   filter(!is.na(non_zero_fped_perc_change)) %>% 
#   summarize(mean_non_zero_fped_perc_change = mean(non_zero_fped_perc_change),
#             sd_non_zero_fped_perc_change = sd(non_zero_fped_perc_change))
#   
# 
# rcws %>% 
#   filter(RCWid == 'ZG-GWO')
# 
# first_transloc_contr_year <- min(results_list$rcw_partial_founder_summed_processed$year[results_list$rcw_partial_founder_summed_processed$group != 'non-transloc'])
#   
# results_list$rcw_partial_founder_summed_processed %>%
#   filter(group != 'non-transloc') %>% 
#   pull(founder_id) %>% 
#   unique() %>% 
#   length()
#   
# 
# results_list$rcw_partial_founder_summed_processed %>%
#   #filter(group != 'non-transloc') %>%
#   group_by(year) %>% 
#   mutate(rank = rank(-pfound_fsum_scale),
#          percentile = 100*(ecdf(pfound_fsum_scale)(pfound_fsum_scale)),
#          percentage_contr = 100*(pfound_fsum_scale/sum(pfound_fsum_scale))) %>% 
#   ungroup() %>% 
#   filter(group != 'non-transloc') %>% 
#   select(year, founder_id, percentage_contr, rank, percentile) %>% 
#   group_by(founder_id) %>% 
#   filter(rank == min(rank))
#   filter(percentage_contr == max(percentage_contr))
# 
# rank(-c(10, 9, 1))
# 
# ecdf(1:10)(1:10)
# 
# results_list$rcw_partial_founder_summed_processed %>%
#   #filter(group != 'non-transloc') %>%
#   group_by(year) %>% 
#   summarize(total_contr = 100*(sum(pfound_fsum_scale[group != 'non-transloc'])/sum(pfound_fsum_scale))) %>% 
#   filter(year >= first_transloc_contr_year) %>% 
#   filter(total_contr == max(total_contr) | total_contr == min(total_contr))
# 
# 
# results_list$rcw_partial_founder_summed_processed %>%
#   filter(group != 'non-transloc') %>%
#   group_by(year) %>%
#   filter(pfound_fsum_scale == max(pfound_fsum_scale)) %>% 
#   ungroup() %>% 
#   arrange(year)
#   summarize(fped = sum(pfound_fsum_scale)) %>% 
#   print(n = 50)
# 
# 
# results_list$rcw_partial_founder_summed_processed %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = pfound_fsum_scale, fill = founder_id,
#                group = -pfound_fsum_scale),
#            color = 'white', linewidth = 0.2,
#            position = "stack", stat = "identity", width = 0.98, alpha = 0.85)


# unique(indiv_transloc_info_weird$id) %in% pop_dat$RCWid
# 
# results_list$rcw_contr_info_df %>% 
#   left_join(., results_list$transloc_info_color %>% 
#               mutate(first_contr_year = year + 1) %>% 
#               rename(id = RCWid) %>% 
#               select(id, first_contr_year, Source),
#             by = 'id') %>% 
#   filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% 
#   filter(id %in% unique(indiv_transloc_info_weird$id)) %>% 
#   pull(contr) %>% 
#   range()


# indiv_transloc_info_establish <- results_list$rcw_contr_info_df %>% 
#   left_join(., results_list$transloc_info_color %>% 
#               mutate(first_contr_year = year + 1) %>% 
#               rename(id = RCWid) %>% 
#               select(id, first_contr_year, Source),
#             by = 'id') %>% 
#   filter(id %in% unique(results_list$transloc_info_color$RCWid)) %>% #filter down to translocated birds
#   filter(year >= first_contr_year) %>% #remove contr values from years before a bird is established
#   filter(id %in% translocations_inter2_process$RCWid[translocations_inter2_process$establish == 1]) #filter down to translocated birds that established
#   
# transloc_finalyear_noncontr <- indiv_transloc_info_establish %>% 
#   filter(year == 2022 & contr <= 1e-15) %>% 
#   pull(id)
# 
# length(transloc_finalyear_noncontr)
# 
# translocations_inter2_process %>% 
#   filter(RCWid %in% transloc_finalyear_noncontr) %>% 
#   mutate(nest = if_else(nest_years > 0, 'yes', 'no')) %>% 
#   group_by(nest) %>% 
#   summarize(count = n())
# 
# pop_dat %>% 
#   filter(RCWid %in% transloc_finalyear_noncontr)
# 
# translocations_inter2_process %>% 
#   filter(RCWid %in% transloc_finalyear_noncontr) %>% 
#   group_by(nest_years) %>% 
#   summarize(count = n())
# 
# pop_dat %>% 
#   filter(RCWid %in% transloc_finalyear_noncontr) %>% 
#   group_by(RCWid) %>% 
#   summarize(last_reported_year = max(year))
# 
# 
# contr_info_rank <- results_list$rcw_contr_info_df %>% 
#   select(!X) %>% 
#   mutate(transloc = if_else(id %in% results_list$transloc_info_color$RCWid, 'yes', 'no')) %>% 
#   #group_by(id) %>% 
#   filter(contr > tol) %>% 
#   #ungroup() %>% 
#   group_by(year) %>% 
#   mutate(rank = rank(-contr, ties.method = 'average'))
# 
# contr_info_rank %>% 
#   filter(year == 2011) %>% 
#   filter(rank <= 7)
# 
# contr_info_rank %>% 
#   filter(year %in% c(2010, 2011, 2012)) %>% 
#   filter(id %in% c('ZB-CA', 'ZB-YO'))
# 
# contr_info_rank %>% 
#   group_by(year) %>% 
#   filter(rank <= 5) %>% 
#   summarize(transloc_count = sum(transloc == 'yes'),
#             .groups = 'drop') %>% 
#   print(n = 100)
# 
# contr_info_rank %>% 
#   filter(id == 'OHA-ZK')
# 
# 
# contr_info_rank %>% 
#   filter(transloc == 'yes') %>% 
#   filter(rank == 1)
# 
# 
# quantify_admix <- function(x) {
#   return(length(x)/(sum(c(1,var(x)), na.rm = TRUE)))
# }
# 
