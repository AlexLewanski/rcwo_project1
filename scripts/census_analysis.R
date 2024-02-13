##############################################################################
##############################################################################
### ANALYSIS AND VISUALIZATION OF RCW POPULATION-LEVEL TRENDS AND PATTERNS ###
##############################################################################
##############################################################################

#####################
### SCRIPT SET-UP ###
#####################

### PACKAGES ###
library(here)
library(tidyverse)
#library(readxl)
library(clock)
library(cowplot)
library(grid)
library(gridExtra)
library(egg)
#library(ggridges)
library(magick)

here::here()


### LOADING CUSTOM FUNCTIONS AND DATA ###

#custom functions
devtools::load_all('/Users/alexlewanski/Documents/r_packages/pedutils')
source(here('scripts', 'data_processing_custom_functions.R'))

file_names <- list.files(here('data', 'feb2024_databasemarch2023_processed'))
rcw_processed_list <- lapply(setNames(file_names, nm = gsub(pattern = '\\.csv', '', file_names)), function(x) {
  read.csv(here('data', 'feb2024_databasemarch2023_processed', x))
})


### ANALYSIS OPTIONS ###
gd_count <- 20000L

### VISUALIZATION CONTENT ###
output_figs <- FALSE

#color palettes
green_colfunc <- colorRampPalette(c("#add6d6", "#329999", "#287a7a", "#1e5b5b"))
red_colfunc <- colorRampPalette(c("#ea9999", "#cc0000", "#a30000", "#7a0000"))

#RCW illustrations
flying_rcw <- image_flop(image_read(here('figures', 'main_paper', 'bird_illustrations', 'rcw_flying1.png')))
perched_rcw <- image_rotate(image_read(here('figures', 'main_paper', 'bird_illustrations', 'rcw_perched1.png')),
                            degrees = 0)



##################################
### INITIAL PROCESSING OF DATA ###
##################################

### extracting different scenarios for annual population composition ###
pop_info_list <- lapply(setNames(nm = paste0('Scenario', 1:4)), function(X, DF) {
  DF[DF[,X,drop = TRUE],!grepl(pattern = 'Scenario', colnames(DF))]
}, DF = rcw_processed_list$census_processed)


transloc_ids <- rcw_processed_list$census_processed %>% 
  filter(Origin_update == 'Translocated') %>% 
  pull(RCWid) %>% 
  unique()

rcws_founder_info <- data.frame(id = get_founders(rcw_processed_list$ped_processed)$founders) %>%
  mutate(group = case_when(id %in% transloc_ids ~ 'transloc',
                           TRUE ~ 'non-transloc'))

transloc_info <- rcw_processed_list$rcws %>% 
  filter(Origin == 'Translocated') %>% 
  select(RCWid, MinAge, Sex) %>%
  rename(year = MinAge)

transloc_info_color <- lapply(split(transloc_info, transloc_info$Sex), function(x) {
  x$sex_color <- switch(x$Sex[1],
                        'M' = {red_colfunc(nrow(x))},
                        'F' = {green_colfunc(nrow(x))})
  return(x)
}) %>% 
  bind_rows()


################
### ANALYSES ###
################

#inbreeding
rcws_inbr <- inbr_from_kinmat(rcw_processed_list$ped_processed[,1:3])

rcw_pop_info_inbr <- pop_info_list[['Scenario3']] %>% 
  left_join(rcws_inbr %>% rename(RCWid = id), by = 'RCWid')


### GENE DROPPING (AND DOWNSTREAM QUANTIFICATIONS) ###
set.seed(345993)
rcws_gdrop <- gene_drop_matrix(ped = rcw_processed_list$ped_processed[,1:3],
                               sims = gd_count,
                               report_progress = TRUE)
set.seed(NULL)

rcws_ancestry_info <- quantify_ancestry(gdrop_mat_output = rcws_gdrop,
                                        founder_group_info = rcws_founder_info)

rcws_ancestry_info_update <- rcws_ancestry_info %>%
  bind_rows(., rcws_ancestry_info %>%
              group_by(id) %>%
              filter(n() < 2) %>%
              filter(!'transloc' %in% group) %>%
              ungroup() %>%
              select(id) %>%
              mutate(group = 'transloc',
                     allele_count = gd_count*2,
                     anc_prop = 0)) %>%
  filter(group == 'transloc')

rcw_inbr_merge <- pop_info_list[['Scenario3']] %>%
  rename(id = RCWid) %>% 
  left_join(rcws_inbr, by = 'id')  %>%
  left_join(rcws_ancestry_info_update, by = 'id')

rcw_partial_founder_fped <- partial_founder_fped_gdrop_mat(gdrop_mat_output = rcws_gdrop)

rcw_processed_partial_founder <- left_join(rcw_partial_founder_fped %>%
                                             rename(allele = allele_origin),
                                           rcws_gdrop$founder_alleles %>%
                                             mutate(allele = as.character(allele)) %>%
                                             rename(founder_id = id),
                                           by = 'allele') %>%
  left_join(., rcws_founder_info %>% rename(founder_id = id),
            by = 'founder_id') %>%
  select(id, group, partial_founder_fped, fped_prop) %>%
  left_join(pop_info_list[['Scenario3']] %>% rename(id = RCWid), ., by = 'id', relationship = "many-to-many") %>%
  filter(!is.na(partial_founder_fped)) %>%
  group_by(year, group) %>%
  summarize(summed_partial_inbr = sum(partial_founder_fped, na.rm = TRUE),
            number_indivs = n(),
            .groups = 'drop') %>%
  left_join(pop_info_list[['Scenario3']] %>%
              group_by(year) %>%
              summarize(pop_size = n(), .groups = 'drop'),
            by = 'year') %>%
  mutate(summed_partial_inbr_std = summed_partial_inbr/pop_size)

rcw_partial_founder_summed <- left_join(rcw_partial_founder_fped %>%
                                          rename(allele = allele_origin),
                                        rcws_gdrop$founder_alleles %>%
                                          mutate(allele = as.character(allele)) %>%
                                          rename(founder_id = id),
                                        by = 'allele') %>%
  left_join(., rcws_founder_info %>% rename(founder_id = id),
            by = 'founder_id') %>%
  select(id, group, partial_founder_fped, fped_prop, founder_id) %>%
  left_join(pop_info_list[['Scenario3']] %>% rename(id = RCWid), ., by = 'id', relationship = "many-to-many") %>% 
  group_by(year) %>%
  mutate(pop_size = n()) %>% 
  ungroup() %>% 
  filter(!is.na(partial_founder_fped))

fped_prop_by_group <- rcw_partial_founder_summed %>% 
  group_by(group, year) %>% 
  summarize(fped_group = sum(partial_founder_fped),
            .groups = 'drop') %>% 
  group_by(year) %>% 
  mutate(fped_group_prop = fped_group/sum(fped_group)) %>% 
  ungroup()

inbr_founder_color <- rcw_partial_founder_summed %>% 
  group_by(founder_id) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  left_join(., transloc_info_color %>% rename(founder_id = RCWid),
            by = 'founder_id') %>% 
  select(founder_id, sex_color) %>% 
  mutate(sex_color = if_else(is.na(sex_color), '#e3e3e3', sex_color)) %>% 
  as.data.frame()


### EFFECTIVE POPULATION SIZE ###
rcws_eqn <- ped_summary_stats(ped = rcw_processed_list$ped_processed, 
                              summary = "equiv_gen")

ne_f_df <- data.frame(year = sort(unique(pop_info_list[['Scenario3']]$year)),
                      ne = NA,
                      ne_std = NA)

for (X in seq_len(nrow(ne_f_df))) {
  message(X)
  ne_val <- calc_ped_Ne(ped = rcw_processed_list$ped_processed,
                        Ne = c('Ne_f', 'Ne_c')[1],
                        id = pop_info_list[['Scenario3']][pop_info_list[['Scenario3']]$year == ne_f_df$year[X],'RCWid'],
                        inbr = rcws_inbr,
                        eqg = rcws_eqn)
  
  ne_f_df[X,'ne'] <- ne_val$Ne_f$Ne
  ne_f_df[X,'ne_std'] <- ne_val$Ne_f$stderr
}





rcw_founder_list_gencontr <- list()
rcw_contr_df <- data.frame(year = sort(unique(pop_info_list[['Scenario3']]$year)),
                           contr = NA)

for (ID in rcws_founder_info$id) {
  rcw_founder_list_gencontr [[ID]] <- rcw_contr_df
  for (year_index in seq_along(rcw_contr_df$year)) {
    
    rcw_founder_list_gencontr [[ID]]$contr[year_index] <- calc_gen_contr(
      ped = rcw_processed_list$ped_processed,
      contributors = ID,
      recipients = pop_info_list[['Scenario3']][pop_info_list[['Scenario3']]$year == rcw_contr_df$year[year_index],]$RCWid,
      focal_contemporaries = NULL,
      standardize = c('none', 'one', 'two')[2]
    )
  }
  message('finished ', ID)
}

rcw_contr_info_df <- bind_rows(rcw_founder_list_gencontr, .id = 'id')



###################################################################
### FIG 1: POPULATION OVERVIEW (POP SIZE, INBREEDING, ANCESTRY) ###
###################################################################

ne_plot <- ne_f_df %>%
  filter(ne < Inf) %>%
  ggplot() +
  geom_segment(aes(x = year,
                   y = ne - ne_std,
                   xend = year,
                   yend = ne + ne_std),
               color = '#b7b7b7', linewidth = 1.5) +
  geom_point(aes(x = year, y = ne),
             shape = 21, colour = '#4c4c4c', fill = NA, 
             size = 1.75, stroke = 1.2) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed', 
                                          color = '#d8d8d8'),
        panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
        axis.ticks = element_line(color = '#808080', linewidth = 0.7),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.1, 1, -0.05, 1), "cm"),
        axis.title.y = element_text(margin = margin(0, 4.75, 0, 0))
        ) +
  scale_x_continuous(limits = c(1993.5, 2022.5), position = "top",
                     expand = c(0.01, 0)) +
  ylab("Ne") +
  coord_cartesian(clip = 'off')


top_bar <- rcw_inbr_merge %>%
  mutate(count = 1) %>%
  group_by(year) %>%
  arrange(f_ped) %>%
  ggplot() +
  geom_bar(aes(x = year, y = count, fill = f_ped),
           position = "stack", stat = "identity", width = 0.98) +
  scale_fill_gradient(name = expression(italic(F)["p"]), 
                      low = '#e5cce5', high = '#730073') +
  geom_point(shape = 108, data = transloc_info_color,
             aes(x = year, y = 0, color = RCWid), size = 4.5,
             alpha = 1,
             position = position_jitter(width = 0.35, height = 0, seed = 5345),
             show.legend = FALSE) +
  scale_color_manual(values = setNames(transloc_info_color$sex_color,
                                       nm = transloc_info_color$RCWid)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 1))) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.border = element_blank(),
        axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
        axis.text.x = element_blank(),
        axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
        axis.ticks = element_line(color = '#808080', linewidth = 0.7),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0.32,1,-0.05,1), "cm"),
        axis.title.y = element_blank()
        ) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
  coord_cartesian(clip = 'off')

bottom_bar <- rcw_inbr_merge %>%
  mutate(count = 1) %>%
  group_by(year) %>%
  arrange(anc_prop) %>%
  ggplot() +
  geom_bar(aes(x = year, y = count, fill = anc_prop),
           position = "stack", stat = "identity", width = 0.98) +
  scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = '#ededed', high = '#4c4c4c') +
  geom_point(shape = 108, data = transloc_info_color,
             aes(x = year, y = 0, color = RCWid), size = 4.5,
             alpha = 1,
             position = position_jitter(width = 0.35, height = 0, seed = 5345),
             show.legend = FALSE) +
  scale_color_manual(values = setNames(transloc_info_color$sex_color,
                                       nm = transloc_info_color$RCWid)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5)
  ) +
  scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(1, 0))) +
  theme(plot.margin=unit(c(-0.05,1,-0.05,1), "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()#,
        ) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
  coord_cartesian(clip = 'off')

transloc_ancestry <- rcw_inbr_merge %>% 
  group_by(year) %>% 
  summarize(transloc = sum(anc_prop)/n(),
            .groups = 'drop') %>% 
  mutate(nontransloc = 1 - transloc) %>% 
  pivot_longer(cols = c(transloc, nontransloc),
               names_to = 'ancestry_type',
               values_to = 'ancestry_prop')  %>% 
  ggplot() +
  geom_bar(aes(x = year, y = ancestry_prop, fill = ancestry_type), 
           color = 'white', linewidth = 0.5, width = 0.72,
           position = "stack", stat = "identity") + 
  scale_fill_manual(values = c('#ededed', '#4c4c4c')) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.border = element_blank(),
        axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
        axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
        plot.margin = unit(c(0.25,1,0,1), "cm"),
        axis.title.x = element_text(size = 18),
        legend.position = 'none',
        axis.title.y = element_text(margin = margin(t = 0, r = 6.6, b = 0, l = 0))
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0.0", "0.5", "1.0")) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
  xlab('Year') +
  ylab("Ancestry")
  
rcw_pop_overview <- ggpubr::annotate_figure(egg::ggarrange(ne_plot,
                                       top_bar,
                                       bottom_bar,
                                       transloc_ancestry,
                                       nrow = 4,
                                       labels = c('A', 'B', '', 'C'),
                                       label.args = list(gp = grid::gpar(font = 2, cex = 1.2)),
                                       heights = c(35, 50, 50, 10)),
                        left = grid::textGrob("Individuals", rot = 90, hjust = 1.04, vjust = 5.2, gp = grid::gpar(cex = 1))) +
  theme(plot.margin = unit(c(0.3, -0.8, 0, -0.45), "cm"))

rcw_pop_overview_rcw <- rcw_pop_overview +
  cowplot::draw_image(flying_rcw,
                      scale = 0.18,
                      x = 0.315, y = 0.35)

if (isTRUE(output_figs)) {
  cowplot::ggsave2(filename = here('figures', 'main_paper', 'rcw_pop_overview_fig.png'),
                   plot = rcw_pop_overview_rcw,
                   width = 10*1, height = 7.25*1, bg = 'white')
}



#################################################################
### FIG 2: GENETIC CONTRIBUTIONS AND MORE INFO. ON INBREEDING ###
#################################################################

genetic_contr_plot <- rcw_contr_info_df %>%
  group_by(id) %>%
  filter(any(contr != 0)) %>% #filter out the ids with all 0s
  filter(contr > 0 | (contr == 0 & year > max(year[contr > 0])  )) %>%
  left_join(., rcws_founder_info,
            by = 'id') %>%
  ggplot() +
  geom_line(data = . %>% 
              filter(group == 'non-transloc'),
            aes(x = year, y = contr, group = id, linewidth = `Founder type`),
            color = '#e3e3e3', linewidth = 0.7) +
  geom_point(data = . %>%
               group_by(id) %>%
               filter(year == min(year)) %>%
               slice_head(n = 1)%>% 
               filter(group == 'non-transloc'),
             aes(x = year, y = contr),
             color = '#e3e3e3', size = 0.7) +
  geom_line(data = . %>% 
              filter(group == 'transloc'),
            aes(x = year, y = contr, group = id, 
                linewidth = `Founder type`, color = id),
            linewidth = 1.1) +
  geom_point(data = . %>%
               group_by(id) %>%
               filter(year == min(year)) %>%
               slice_head(n = 1)%>% 
               filter(group == 'transloc'),
             aes(x = year, y = contr, color = id),
             size = 1.5) +
  scale_color_manual(values = setNames(transloc_info_color$sex_color,
                                       nm = transloc_info_color$RCWid)) +
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
        plot.margin = unit(c(0.32,1, 0, 1), "cm"),
        legend.position = 'none',
        axis.title.y = element_text(margin = margin(0, 12.5, 0, 0))
  ) +
  scale_x_continuous(expand = c(0.01, 0), 
                     limits = c(1993.5, 2022.5),
                     position = "top") +
  xlab('Year') +
  ylab("Exp. genetic contribution")


fped_prop_by_indiv_founder_plot <- rcw_partial_founder_summed %>% 
  group_by(founder_id, year) %>%
  summarize(partial_founder_fped_sum = sum(partial_founder_fped),
            group = first(group),
            pop_size = first(pop_size),
            .groups = 'drop') %>%
  group_by(year) %>%
  mutate(partial_founder_fped_sum_prop = partial_founder_fped_sum/sum(partial_founder_fped_sum)) %>% 
  arrange(partial_founder_fped_sum_prop) %>% 
  mutate(year_rank = factor(n():1, levels = n():1) ) %>% 
  ungroup() %>% 
  mutate(partial_founder_fped_sum_scaled = partial_founder_fped_sum/pop_size) %>% 
  ggplot() +
  geom_bar(aes(x = year, y = partial_founder_fped_sum_scaled, group = year_rank, fill = founder_id),
           color = 'white', linewidth = 0.2,
           position = "stack", stat = "identity", width = 0.98, alpha = 0.9) +
  scale_fill_manual(values = setNames(inbr_founder_color$sex_color,
                                       nm = inbr_founder_color$id)) +
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
    plot.margin = unit(c(0.25,1,-0.05,1), "cm"),
    legend.position = 'none') +
  ylab(substitute(paste("Scaled ", italic(F)["p"]))) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = 'off')


fped_prop_by_group_plot <- fped_prop_by_group %>% 
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
        axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
        axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
        plot.margin = unit(c(0.25,1,0,1), "cm"),
        axis.title.x = element_text(size = 18),
        legend.position = 'none',
        axis.title.y = element_text(margin = margin(t = 0, r = 14, b = 0, l = 0))
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0.0", "0.5", "1.0")) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
  xlab('Year') +
  ylab(substitute(paste(italic(F)["p"], ' prop.')))


indiv_fped_gencontr_panel <- egg::ggarrange(genetic_contr_plot,
               fped_prop_by_indiv_founder_plot,
               fped_prop_by_group_plot,
               nrow = 3,
               labels = c('A', 'B', 'C'),
               label.args = list(gp = grid::gpar(font = 2, cex = 1.2)),
               heights = c(65, 65, 10))

indiv_fped_gencontr_panel_rcw <- cowplot::ggdraw(indiv_fped_gencontr_panel) +
  cowplot::draw_image(perched_rcw,
                      scale = 0.195,
                      x = -0.3695, y = -.15) +
  theme(plot.margin = unit(c(0.2, -0.7, 0, 0.1), "cm"))

if (isTRUE(output_figs)) {
  cowplot::ggsave2(filename = here('figures', 'main_paper', 'indiv_fped_gencontr_panel.png'),
                   plot = indiv_fped_gencontr_panel_rcw,
                   width = 10*1, height = 6.5*1, bg = 'white')
}




#################################
### CODE NOT CURRENTLY IN USE ###
#################################
#rcw_pop_info_s3 <- rcw_processed_list$census_processed[rcw_processed_list$census_processed$Scenario3,]

# rcw_pop_overview <- ggpubr::annotate_figure(egg::ggarrange(ne_plot,
#                                        top_bar,
#                                        bottom_bar,
#                                        nrow = 3,
#                                        heights = c(35, 50, 50)),
#                         left = grid::textGrob("Individuals", rot = 90, hjust = 1, vjust = 5.2, gp = grid::gpar(cex = 1))) +
#   theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm"))
# 
# cowplot::ggsave2(filename = here('figures', 'main_paper', 'rcw_pop_overview_fig.png'),
#                  plot = rcw_pop_overview,
#                  width = 10*1, height = 7*1, bg = 'white')


# bottom_bar <- rcw_inbr_merge %>%
#   mutate(count = 1) %>%
#   group_by(year) %>%
#   arrange(anc_prop) %>%
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = anc_prop),
#            position = "stack", stat = "identity", width = 0.98) +
#   scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = '#ededed', high = '#4c4c4c') +
#   geom_point(shape = 108, data = transloc_info_color,
#              aes(x = year, y = 0, color = RCWid), size = 4.5,
#              alpha = 1,
#              position = position_jitter(width = 0.35, height = 0, seed = 5345)) +
#   scale_color_manual(values = setNames(transloc_info_color$sex_color,
#                                        nm = transloc_info_color$RCWid)) +
#   theme_bw() +
#   theme(panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#         panel.border = element_blank(),
#         axis.line.x = element_line(color = '#4c4c4c', size = 1.1),
#         axis.line.y = element_line(color = '#4c4c4c', size = 1.5),
#   ) +
#   scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(1, 0))) +
#   theme(plot.margin=unit(c(-0.05,1,1,1), "cm"),
#         axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 18),
#         legend.position = 'none') +
#   xlab('Year') + xlim(1993.5, 2022.5) +
#   coord_cartesian(clip = 'off')

# rcw_pop_overview <- ggpubr::annotate_figure(egg::ggarrange(ne_plot,
#                                                            top_bar,
#                                                            bottom_bar_alt,
#                                                            transloc_contr,
#                                                            nrow = 4,
#                                                            labels = c('A', 'B', '', 'C'),
#                                                            label.args = list(gp = grid::gpar(font = 2, cex = 1.2)),
#                                                            heights = c(35, 50, 50, 10)),
#                                             left = grid::textGrob("Individuals", rot = 90, hjust = 1.18, vjust = 5.2, gp = grid::gpar(cex = 1))) +
#   theme(plot.margin = unit(c(0.3, 0, 0, 0), "cm"))
# 
# rcw_pop_overview_rcw <- rcw_pop_overview +
#   cowplot::draw_image(perched_rcw,
#                       scale = 0.17,
#                       x = 0.315, y = 0.35)


# fped_decomp_transloc_plot <- rcw_processed_partial_founder %>%
#   ggplot() +
#   geom_bar(aes(x = year, y = summed_partial_inbr_std, fill = group),
#            #color = 'gray', size = 0.2,
#            position = "stack", stat = "identity", width = 0.98, alpha = 0.8) +
#   theme_bw() +
#   theme(#panel.grid.major = element_blank(),
#     #panel.grid.minor = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#     panel.border = element_blank(),
#     axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
#     #axis.line.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.1),
#     axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     plot.margin = unit(c(1,1,-0.05,1), "cm")) +
#   scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
#   ylab("Expected inbreeding") +
#   #scale_fill_manual(values = c("#00cdcd", '#e166e1'))
#   scale_fill_manual(name = "Founder type",
#                     values = c('#ededed', '#4c4c4c'))
#
# perind_partial_fped <- left_join(rcw_partial_founder_fped %>%
#                                    rename(allele = allele_origin),
#                                  rcws_gdrop$founder_alleles %>%
#                                    mutate(allele = as.character(allele)) %>%
#                                    rename(founder_id = id),
#                                  by = 'allele') %>%
#   left_join(., rcws_founder_info %>% rename(founder_id = id),
#             by = 'founder_id') %>%
#   select(id, group, partial_founder_fped, fped_prop, founder_id) %>%
#   left_join(pop_info_list[['Scenario3']] %>% rename(id = RCWid), ., by = 'id', relationship = "many-to-many") %>% 
#   #mutate(partial_founder_fped = case_when(is.na(partial_founder_fped) ~ 0,
#   #                                        TRUE ~ partial_founder_fped)) %>% 
#   #filter(!is.na(partial_founder_fped)) %>% 
#   group_by(id, year) %>% 
#   summarize(transloc_fped_perind = sum(partial_founder_fped[group == "transloc"], na.rm = TRUE),
#             total_fped = sum(partial_founder_fped, na.rm = TRUE),
#             .groups = 'drop')
# 
# perind_partial_fped %>% 
#   mutate(transloc_prop = transloc_fped_perind/total_fped,
#          transloc_prop = if_else(is.na(transloc_prop), 0, transloc_prop)) %>% 
#   arrange(transloc_prop) %>% 
#   ggplot(aes(x = year, y = total_fped, fill = transloc_prop)) +
#   geom_point(color = '#808080', pch = 21, size = 3.25, stroke = 0.25, alpha = 0.8,
#              position = position_jitter(width = 0.15)) +
#   theme_bw() +
#   theme(#panel.grid.major = element_blank(),
#     #panel.grid.minor = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#     panel.border = element_blank(),
#     axis.line.x = element_line(color = '#4c4c4c', linewidth = 1.1),
#     #axis.line.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.1),
#     axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank()) +
#   theme(plot.margin = unit(c(1,1,-0.05,1), "cm")
#         #axis.title.y = element_blank()
#   ) +
#   xlim(1993.5, 2022.5) +
#   ylab("Inbreeding") +
#   scale_fill_gradient2(low = "#00cdcd", mid="#f2f2f2", high= '#e166e1', midpoint = 0.5)
# 
# 
# 
# 
# 
# perind_partial_fped %>% 
#   mutate(transloc_prop = transloc_fped_perind/total_fped,
#          transloc_prop = if_else(is.na(transloc_prop), 0, transloc_prop)) %>% 
#   mutate(dif_trans_ped = (0.5*total_fped)  - transloc_prop) %>% 
#   group_by(id) %>% 
#   slice_head(n = 1) %>% 
#   ggplot() +
#   geom_point(aes(x = dif_trans_ped, y = total_fped)) +
#   xlim(-1, 1)



# 
# rcw_contr_info_processed <- rcw_contr_info_df %>%
#   group_by(id) %>%
#   filter(any(contr != 0)) %>% #filter out the ids with all 0s
#   filter(contr > 0 | (contr == 0 & year > max(year[contr > 0])  )) %>%
#   left_join(., rcws_founder_info,
#             by = 'id')
# 
# rcw_contr_info_processed_transloc <- rcw_contr_info_processed %>% 
#   filter(group == 'transloc') %>% 
#   left_join(., color_info, by = 'id')
# 
# 
# rcw_contr_info_processed_transloc %>% 
#   ggplot() +
#   geom_line(aes(x = year, y = contr, color = id, linewidth = `Founder type`), linewidth = 1.1) +
#   scale_color_manual(values = unique(rcw_contr_info_processed_transloc$green_col)) +
#   theme(legend.position = 'none')
# 

# rcw_contr_info_df %>%
#   group_by(id) %>%
#   filter(any(contr != 0)) %>% #filter out the ids with all 0s
#   filter(contr > 0 | (contr == 0 & year > max(year[contr > 0])  )) %>%
#   left_join(., rcws_founder_info,
#             by = 'id') %>%
#   #mutate(`Founder type` = factor(group, levels = c('transloc', 'non-transloc'))) %>%
#   ggplot() +
#   geom_line(data = . %>%
#               filter(group == 'non-transloc'),
#             aes(x = year, y = contr, group = id, linewidth = `Founder type`),
#             color = '#d8d8d8', linewidth = 0.75) +
#   geom_point(data = . %>%
#                group_by(id) %>%
#                filter(year == min(year)) %>%
#                slice_head(n = 1)%>%
#                filter(group == 'non-transloc'),
#              aes(x = year, y = contr),
#              color = '#d8d8d8', size = 0.75) +
#   geom_line(data = . %>%
#               filter(group == 'transloc'),
#             aes(x = year, y = contr, group = id, linewidth = `Founder type`),
#             color = '#6495ED', linewidth = 1.1) +
#   geom_point(data = . %>%
#                group_by(id) %>%
#                filter(year == min(year)) %>%
#                slice_head(n = 1)%>%
#                filter(group == 'transloc'),
#              aes(x = year, y = contr),
#              color = '#6495ED', size = 1.5) +
#   #scale_color_manual(values = c('#6495ED', '#d8d8d8')) +
#   #scale_size_manual(values = c(1.5, 0.75)) +
#   #scale_linewidth_manual(values = c(1.1, 0.75)) +
#   theme_classic() +
#   theme(panel.grid.major.x = element_line(color = '#d8d8d8',
#                                           linetype = 'dashed',
#                                           linewidth = 0.4)) +
#   xlab('Year') +
#   ylab("Expected genetic contribution")
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
# last_seen_detail <- read_excel(here('data',
#                                     'march_2023_database',
#                                     'LastSeenDetail.xlsx'))
# 
# rcws <- read_excel(here('data',
#                         'march_2023_database',
#                         'RCWs.xlsx'))
# 
# nests <- read_excel(here('data',
#                          'march_2023_database',
#                          'Nests.xlsx'))
#   
# 
# table(rcws$Translocation)
# 
# #the vast majority of summer censuses occur in July (7), but a few occur at the
# #end of June (6) or the beginning of August (8)
# last_seen_summer_census <- last_seen_detail %>% 
#   filter(Detected == 'Yes' & 
#          SurveyType == 'Census' & 
#            get_month(RecordDate) %in% c(6, 7, 8)) %>% 
#   mutate(year = get_year(RecordDate))
# 
# 
# 
# ######################
# ### QUALITY CHECKS ###
# ######################
# 
# #brief quality check: make sure that each indiv is only found once each year
# id_count_per_year <- last_seen_summer_census %>% 
#   group_by(RCWid, year) %>% 
#   summarize(count = n(), .groups = 'drop') %>%
#   pull(count)
# 
# #equivalent to:
# #all(table(last_seen_summer_census$RCWid,
# #      last_seen_summer_census$year) < 2)
# 
# if (all(id_count_per_year < 2)) {
#   message('Check passed: all IDs are only found once per year')
# } else {
#   stop('Check NOT passed: at least one ID is duplicated in a year')
# }
# 
# 
# survey_per_year_count <- last_seen_detail %>% 
#   filter(Detected == 'Yes' & 
#            SurveyType == 'Census' & 
#            get_month(RecordDate) %in% c(6, 7, 8)) %>% 
#   group_by(RecordDate) %>% 
#   summarize(survey_dates = first(RecordDate), .groups = 'drop') %>% 
#   mutate(year = get_year(RecordDate)) %>% 
#   group_by(year) %>% 
#   summarize(summer_survey_count = n()) %>%
#   pull(summer_survey_count)
# 
# if (all(survey_per_year_count == 1)) {
#   message('Check passed: each year only has one summer survey')
# } else {
#   stop('Check NOT passed: at least one year has more than one summer survey')
# }
# 
# 
# summer_census_interp <- add_intervening_years(data = last_seen_summer_census,
#                                                id_col = 'RCWid',
#                                                year_col = 'year')
# 
# id_count_per_year_post_interp <- summer_census_interp %>% 
#   group_by(RCWid, year) %>% 
#   summarize(count = n(), .groups = 'drop') %>%
#   pull(count)
# 
# if (all(id_count_per_year_post_interp == 1)) {
#   message('Check passed: each indiv is only found once a year year post-interpolation')
# } else {
#   stop('Check NOT passed: at least one individual is found >1 times per year')
# }
# 
# 
# 
# 
# 
# 
# check_census_year_df <- summer_census_interp %>% 
#   group_by(RCWid) %>% 
#   filter(year == min(year)) %>% 
#   rename(first_year = year) %>% 
#   slice_head(n = 1) %>% 
#   ungroup() %>% 
#   select(RCWid, first_year) %>% 
#   left_join(rcws,.,
#             by = 'RCWid') %>%
#   rowwise() %>% 
#   mutate(year_dif = case_when(is.na(MinAge) | is.na(first_year) ~ 100,
#                               TRUE ~ first_year - MinAge)) 
# 
# nrow(check_census_year_df[check_census_year_df$year_dif == 100,])/nrow(check_census_year_df)
# nrow(check_census_year_df[check_census_year_df$year_dif > 0 & check_census_year_df$year_dif != 100,])/nrow(check_census_year_df)
# 
# 
# 
# check_census_year_df %>% 
#   select(RCWid, year_dif, MinAge, DateNew, first_year) %>% 
#   filter(year_dif > 1 & year_dif != 100)
# 
# check_census_year_df %>% 
#   select(RCWid, year_dif, MinAge, DateNew, first_year) %>% 
#   filter(year_dif < 0)
# 
# summer_census_interp %>% 
#   filter(grepl('^GB[-=]Z$', RCWid) )
# 
# rcws %>% 
#   filter(grepl('^GB[-=]Z$', RCWid)) %>% 
#   select(RCWid, DateNew, MinAge, Origin)
# 
# nests %>% 
#   filter(grepl('^GB[-=]Z$', MaleID) | grepl('^GB[-=]Z$', FemaleID)) %>% 
#   select(NatalNest_Temp, StartDate, MaleID, FemaleID)
# 
# rcws %>% 
#   filter(RCWid == 'GB=Z')
# 
# check_census_year_df %>%
#   ggplot() +
#   geom_histogram(aes(x = year_dif), binwidth = 1) +
#   scale_x_continuous(limits=c(-12, 8), breaks = -12:8) +
#   theme_bw()
# 
# rcws %>% 
#   #select(RCWid, MinAge, DateNew, AgeKnown) %>% 
#   filter(is.na(MinAge))
# 
# 
# table(rcws$MinAge)
# 
# 
# 
# last_seen_summer_census_update <- last_seen_summer_census %>% 
#   mutate(RCWid_old = RCWid) %>% 
#   mutate(RCWid = case_when(get_year(RecordDate) < 2005 & RCWid == 'GB=Z' ~ 'GB-Z',
#                            TRUE ~ RCWid))
# 
# 
# summer_census_interp_idupdate <- add_intervening_years(data = last_seen_summer_census_update,
#                                                        id_col = 'RCWid',
#                                                        year_col = 'year')
# 
# 
# 
# 
# ### comparing the pre- and post-interpolation population sizes ###
# pop_size_dataframe <- left_join(
#   summer_census_interp_idupdate %>% 
#     group_by(year) %>% 
#     summarize(Interpolated = n(), .groups = 'drop'),
#   last_seen_summer_census %>% 
#     group_by(year) %>% 
#     summarize(Observed = n(), .groups = 'drop'),
#   by = 'year'
# )
# 
# pop_corval <- cor(pop_size_dataframe$Observed, pop_size_dataframe$Interpolated, 
#                   method = 'spearman')
# 
# 
# pop_size_joint_plot <- pop_size_dataframe %>% 
#   pivot_longer(cols = c(Observed, Interpolated),
#                names_to = 'count_type', values_to = 'Census count') %>% 
#   ggplot() +
#   geom_line(aes(x = year, y = `Census count`, color = count_type),
#              alpha = 0.4, linewidth = 0.9) +
#   geom_point(aes(x = year, y = `Census count`, color = count_type),
#              alpha = 0.6, size = 4) +
#   scale_color_manual(values = c('#2a9d8f', '#e76f51')) +
#   theme_bw() +
#   theme(panel.grid.major = element_line(linewidth = 0.35, color = '#f2f2f2'),
#         panel.grid.minor = element_line(linewidth = 0.35, linetype = 'dashed', color = '#f2f2f2'),
#         panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
#         legend.title = element_blank(),
#         axis.title = element_text(size = 15),
#         legend.position = "bottom",
#         axis.line = element_line(linewidth = 1, lineend = 'round')) +
#   xlab('Year') +
#   ylab('Census count')
# 
# 
# pre_post_pop_plot <- pop_size_dataframe %>% 
#   ggplot() +
#   geom_abline(slope = 1, intercept = 0, linetype = 'dashed', linewidth = 2, color = '#808080') +
#   geom_point(aes(x = Observed, y = Interpolated),
#              alpha = 0.6, size = 5, color = "#03045e") +
#   theme_bw() +
#   theme(panel.grid.major = element_line(linewidth = 0.35, color = '#f2f2f2'),
#         panel.grid.minor = element_line(linewidth = 0.35, linetype = 'dashed', color = '#f2f2f2'),
#         panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
#         axis.title = element_text(size = 15),
#         plot.title = element_text(size = 20),
#         axis.line = element_line(linewidth = 1, lineend = 'round')) +
#   xlab('Observed census count') +
#   ylab('Interpolation census count') +
#   ggtitle(paste0('Correlation: ', round(pop_corval, 4)) )
# 
# 
# (plot_compare_multipanel <- plot_grid(pre_post_pop_plot, pop_size_joint_plot, 
#                                      ncol = 2, align = 'hv', axis = "tblr"))
# 
# if (isTRUE(output_figs)) {
#   ggsave(filename = here("figures", "observed_vs_interpolated_popsize.png"),
#          plot = plot_compare_multipanel,
#          width = 10, height = 5)
# }
# 
# 
# 
# 
# 
# 
# 
# inbr_calcs <- lapply(seq(0, 0.11, by = 0.01), function(inbr_val) {
#   inbr_df <- stack_v2(
#     inbreeding_wrapper(ped = rcw_ped,
#                        clear_founder_inbreeding = TRUE,
#                        founder_inbreeding_val = inbr_val,
#                        founder_ids = NULL,
#                        inbreeding_ids = NULL, 
#                        Xchrom = FALSE)
#   ) %>% 
#     mutate(founder_fped_val = inbr_val)
#   rownames(inbr_df) <- NULL
#   return(inbr_df)
# }) %>% 
#   bind_rows() %>% 
#   rename(RCWid = ind,
#          fped = values) %>% 
#   pivot_wider(names_from = founder_fped_val, 
#               values_from = fped, 
#               names_prefix = 'founder_')
# 
# 
# combined_census_indiv_info <- left_join(summer_census_interp_idupdate,
#           rcws[,c('RCWid', "Translocation")],
#           by = 'RCWid') %>% 
#   left_join(., inbr_calcs,
#             by = 'RCWid')
# 
# 
# 
# 
# founders(rcw_ped)
# inbr_calcs %>% 
#   left_join(., nests_long %>% 
#               group_by(RCWid) %>% 
#               summarize(EggNum = sum(EggNum, na.rm = TRUE),
#                         HatchNum = sum(HatchNum, na.rm = TRUE),
#                         FldgNum = sum(FldgNum, na.rm = TRUE)),
#               by = 'RCWid') %>% 
#   mutate(EggNum = if_else(is.na(EggNum), 0, EggNum),
#          HatchNum = if_else(is.na(HatchNum), 0, HatchNum),
#          FldgNum = if_else(is.na(FldgNum), 0, FldgNum)) %>% 
#   filter(!RCWid %in% founders(rcw_ped)) %>% 
#   ggplot(aes(x = founder_0, y = FldgNum)) +
#   geom_point() +
#   geom_smooth(method='lm')
# 
# 
# 
# prac_inbr_depression <- inbr_calcs %>% 
#   left_join(., nests_long %>% 
#               group_by(RCWid) %>% 
#               summarize(EggNum = sum(EggNum, na.rm = TRUE),
#                         HatchNum = sum(HatchNum, na.rm = TRUE),
#                         FldgNum = sum(FldgNum, na.rm = TRUE)),
#             by = 'RCWid') %>% 
#   mutate(EggNum = if_else(is.na(EggNum), 0, EggNum),
#          HatchNum = if_else(is.na(HatchNum), 0, HatchNum),
#          FldgNum = if_else(is.na(FldgNum), 0, FldgNum)) %>% 
#   left_join(., rcws[,c('RCWid', 'Sex')], by = 'RCWid')
# 
# prac_inbr_depression %>% 
#   #filter(!RCWid %in% founders(rcw_ped)) %>% 
#   ggplot(aes(x = founder_0, y = HatchNum)) +
#   geom_point() +
#   geom_smooth(method = 'lm', method.args = list(family = "poisson"))
# 
# library(rsq)
# 
# female_mod <- glm(HatchNum ~ founder_0.01, 
#         family = 'poisson', 
#         data = prac_inbr_depression %>% 
#           filter(Sex == 'F')
#             )
# 
# male_mod <- glm(HatchNum ~ founder_0.01, 
#                   family = 'poisson', 
#                   data = prac_inbr_depression %>% 
#                     filter(Sex == 'M')
# )
# 
# mod_sex <- glm(FldgNum ~ founder_0.01 + Sex + founder_0.01*Sex, 
#     family = 'poisson', 
#     data = prac_inbr_depression %>% filter(Sex != 'U')
# )
# 
# rsq(mod_sex)
# 
# summary(mod_sex)
# summary(female_mod)
# 
# summary(glm(HatchNum ~ founder_0, family = 'poisson', data = prac_inbr_depression #%>% filter(!RCWid %in% founders(rcw_ped))
# ))
# 
# rsq(glm(HatchNum ~ founder_0, family = 'poisson', data = prac_inbr_depression %>% 
#           filter(!RCWid %in% founders(rcw_ped)) %>% 
#           filter(founder_0 < 0.2)
# ))
# 


# 
# 
# 
# inbreeding_facet_differentfounder <- combined_census_indiv_info %>% 
#   pivot_longer(cols = starts_with('founder'), names_to = 'Founder_inbreeding_val', values_to = 'f_ped') %>% 
#   group_by(year, Founder_inbreeding_val) %>% 
#   summarize(mean_fped = mean(f_ped),
#             lower_quantile = quantile(f_ped, probs = 0.025),
#             upper_quantile = quantile(f_ped, probs = 0.975)) %>%
#   ungroup() %>% 
#   ggplot(aes(x = year, y = mean_fped)) +
#   geom_hline(yintercept = 0, color = '#999999', size = 0.5) +
#   geom_line(alpha = 0.5, size = 1, color = "#03045e") +
#   geom_point(alpha = 0.6, size = 3, color = "#03045e") +
#   theme_bw() +
#   theme(#panel.grid.major = element_line(size = 0.35, color = '#f2f2f2'),
#         #panel.grid.minor = element_line(size = 0.35, linetype = 'dashed', color = '#f2f2f2'),
#         panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#         axis.title = element_text(size = 15),
#         plot.title = element_text(size = 20),
#         axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#         strip.background = element_rect(fill = "#737373"),
#         strip.text = element_text(color = "white", size = 12)) +
#   facet_wrap(~Founder_inbreeding_val, scales = "free_y") +
#   xlab('Year') +
#   ylab('Mean inbreeding')
# 
# ggsave(filename = here("figures", "inbreeding_facet_differentfounder_inbr.png"),
#        plot = inbreeding_facet_differentfounder,
#        width = 11, height = 7)
# 
# inbreeding_facet_differentfounder_observeddat <- combined_census_indiv_info %>% 
#   pivot_longer(cols = starts_with('founder'), 
#                names_to = 'Founder_inbreeding_val', 
#                values_to = 'f_ped') %>% 
#   ggplot(aes(x = year, y = f_ped)) +
#   geom_point(alpha = 0.5, size = 0.6, color = "gray") +
#   geom_hline(yintercept = 0, color = '#999999', size = 0.5) +
#   geom_line(data = . %>% group_by(year, Founder_inbreeding_val) %>% 
#               summarize(mean_fped = mean(f_ped),
#                         lower_quantile = quantile(f_ped, probs = 0.025),
#                         upper_quantile = quantile(f_ped, probs = 0.975)) %>%
#               ungroup(),
#             aes(x = year, y = mean_fped),
#             alpha = 0.5, size = 1, color = "#03045e") +
#   geom_point(data = . %>% group_by(year, Founder_inbreeding_val) %>% 
#                summarize(mean_fped = mean(f_ped),
#                          lower_quantile = quantile(f_ped, probs = 0.025),
#                          upper_quantile = quantile(f_ped, probs = 0.975)) %>%
#                ungroup(),
#              aes(x = year, y = mean_fped),
#              alpha = 0.6, size = 3, color = "#03045e") +
#   theme_bw() +
#   theme(#panel.grid.major = element_line(size = 0.35, color = '#f2f2f2'),
#     #panel.grid.minor = element_line(size = 0.35, linetype = 'dashed', color = '#f2f2f2'),
#     panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#     panel.grid.minor = element_blank(),
#     panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#     axis.title = element_text(size = 15),
#     plot.title = element_text(size = 20),
#     axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#     strip.background = element_rect(fill = "#737373"),
#     strip.text = element_text(color = "white", size = 12)) +
#   facet_wrap(~Founder_inbreeding_val, scales = "free_y") +
#   xlab('Year') +
#   ylab('Mean inbreeding')
# 
# ggsave(filename = here("figures", "inbreeding_facet_differentfounder_inbr_withobserved.png"),
#        plot = inbreeding_facet_differentfounder_observeddat,
#        width = 11, height = 7)
# 
# 
# combined_census_indiv_mig_status <- combined_census_indiv_info %>%
#   group_by(RCWid) %>% 
#   mutate(status = case_when(RCWid %in% translocated_ids & year == min(year) ~ 'migrant (first gen)', #init migrants
#                             RCWid %in% translocated_ids & year != min(year) ~ 'migrant (later gen)', #later year migrants
#                             RCWid %in% translocated_descendants ~ 'migrant descendant', #migrant descendants
#                             TRUE ~ 'non-migrant')) %>% 
#   ungroup()
# 
# 
# 
# 
# 
# top_bar <- combined_census_indiv_mig_status %>% 
#   mutate(count = 1) %>% 
#   group_by(year) %>% 
#   arrange(founder_0) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = founder_0), 
#            #color = 'gray', size = 0.2,
#            position = "stack", stat="identity", width = 0.98) +
#   scale_fill_gradient(low = '#e5cce5', high = '#730073',
#                       name = bquote(F[ped])) +
#   scale_x_continuous(expand = c(0.01, 0)) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 3))) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line.x = element_line(color = '#4c4c4c', size = 1.1, lineend = 'butt'),
#         axis.text.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', size = 1.5, lineend = 'butt'),
#         axis.ticks = element_line(color = '#808080', size = 0.7),
#         axis.title.x = element_blank(), 
#         axis.ticks.x = element_blank()) +
#   theme(plot.margin = unit(c(1, 1, -0.05,1), "cm"),
#         axis.title.y = element_blank())
# 
# 
# gray_vec <- c('#e4e4e4', '#c9c9c9', '#a6a6a6', '#737373')
# status_vec <- c('non-migrant', 'migrant descendant', 'migrant (later gen)', 'migrant (first gen)')
# bottom_bar <- combined_census_indiv_mig_status %>% 
#   mutate(count = 1) %>% 
#   mutate(migrant = factor(status, 
#                           levels = status_vec)) %>% 
#   #mutate(migrant = factor(migrant, levels = c('non-migrant',  'migrant descendant', 'migrant'))) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = migrant), 
#            #color = 'gray', size = 0.2,
#            position = "stack", stat = "identity", width = 0.98) +
#   scale_fill_manual(values = gray_vec) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line.x = element_line(color = '#4c4c4c', size = 1.1, lineend = 'butt'),
#         axis.line.y = element_line(color = '#4c4c4c', size = 1.5, lineend = 'butt'),
#         legend.title = element_blank()) +
#   scale_x_continuous(expand = c(0.01, 0)) +
#   scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(3, 0))) +
#   theme(plot.margin = unit(c(-0.05, 1, 1, 1), "cm"),
#         axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 18)) +
#   xlab('Year')
# 
# pop_descripton_multipanel <- ggpubr::annotate_figure(ggarrange(top_bar, bottom_bar, nrow = 2), 
#                         left = textGrob("Count", rot = 90, hjust = 0.1, vjust = 2.3, gp = gpar(cex = 1.3)))
# 
# 
# ggsave(filename = here("figures", "pop_descripton_multipanel.png"),
#        plot = pop_descripton_multipanel,
#        width = 11*0.8, height = 7*0.8, bg = "white")
# 
# # combined_census_indiv_mig_status %>% 
# #   mutate(year = factor(year),
# #          rand_val = rnorm(n(), 0.5*(as.numeric(as.character(year)) - min(as.numeric(as.character(year)))), 5)) %>% 
# #   ggplot(aes(x = founder_0, y = year, color = rand_val)) +
# #   #geom_density_ridges(fill = '#d8d8d8', color = '#b97fb9', size = 0.65, scale = 1, rel_min_height = 0.005) + #fill = '#bfbfbf'
# #   geom_density_ridges(aes(point_fill = rand_val), point_shape = 21,
# #                       alpha = 0, color = NA, jittered_points = TRUE, point_alpha = 1, scale = 0.85) +
# #   #scale_point_color_gradient(low = '#e5cce5', high = '#730073') +
# #   scale_point_fill_gradient(low = 'white', high = '#730073') +
# #   xlab('Fped') + ylab('Year') +
# #   theme_bw(base_size = 13) +
# #   theme(panel.grid.major.x = element_line(linetype = 'dashed', size = 0.65, color = '#f2f2f2'),
# #         #panel.grid.major.x = element_blank(),
# #         #panel.grid.minor.x = element_line(linetype = 'dashed', size = 0.5, color = '#f2f2f2'),
# #         panel.grid.minor.x = element_blank(),
# #         panel.grid.major.y = element_line(linetype = 'solid', size = 0.7, color = '#f2f2f2'),
# #         panel.border = element_blank(),
# #         axis.line = element_line(linetype = 'solid', size = 0.55, color = '#4c4c4c'),
# #         axis.ticks = element_line(color = '#808080', size = 0.7),
# #         axis.text = element_text(size = 9.5))
# 
# 
# str(prac_df)
# str(combined_census_indiv_mig_status)
# 
# 
# 
# 
# 
# 
# year_vec <- 1995:2020
# 
# prac_df <- data.frame(year = factor(rep(year_vec, times = sample(30:70, length(year_vec), replace = TRUE)), levels = as.character(year_vec)))
# 
# prac_df$f_ped <-  as.numeric(as.character(prac_df$year)) - min(as.numeric(as.character(prac_df$year))) + rnorm(nrow(prac_df), mean = 0, sd = 8)  #runif(nrow(prac_df), 0, 1)
# prac_df$count <- rep(1, nrow(prac_df))
# #prac_df$migrant <- if_else(rbinom(n = nrow(prac_df), size = 1, p = 0.1) == 1, 'migrant', 'nonmigrant')
# prac_df$migrant <- sample(c('non-migrant', 'migrant (first gen)', 'migrant (later gen)'), nrow(prac_df), c(0.8, 0.1, 0.1), replace = TRUE)
# 
# top_bar <- prac_df %>% 
#   group_by(year) %>% 
#   arrange(f_ped) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = f_ped), 
#            #color = 'gray', size = 0.2,
#            position = "stack", stat = "identity", width = 0.98) +
#   scale_fill_gradient(low = '#e5cce5', high = '#730073') +
#   scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 1))) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line.x = element_line(color = '#4c4c4c', size = 1.1),
#         axis.text.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', size = 1.5),
#         axis.ticks = element_line(color = '#808080', size = 0.7),
#         axis.title.x = element_blank(), 
#         axis.ticks.x = element_blank()) +
#   theme(plot.margin=unit(c(1,1,-0.05,1), "cm"),
#         axis.title.y = element_blank())
# 
# bottom_bar <- prac_df %>% 
#   mutate(migrant = factor(migrant, levels = c('non-migrant',  'migrant (later gen)', 'migrant (first gen)'))) %>% 
#   #mutate(migrant = factor(migrant, levels = c('non-migrant',  'migrant descendant', 'migrant'))) %>% 
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = migrant), 
#            #color = 'gray', size = 0.2,
#            position = "stack", stat = "identity", width = 0.98) +
#   scale_fill_manual(values = c('#d8d8d8', '#a6a6a6', '#808080')) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line.x = element_line(color = '#4c4c4c', size = 1.1),
#         axis.line.y = element_line(color = '#4c4c4c', size = 1.5),
#         legend.title = element_blank()) +
#   scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(1, 0))) +
#   theme(plot.margin=unit(c(-0.05,1,1,1), "cm"),
#         axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 18)) +
#   xlab('Year')
# #geom_hline(yintercept = c(5, 10, 15, 20), color = 'white', size = 0.25, linetype = 'dashed')
# 
# 
# 
# nests %>% 
#   rowwise() %>% 
#   filter(any(is.na(c(MaleID, FemaleID)))) %>% 
#   select(NatalNest_Temp, MaleID, FemaleID, Year)
# nests$MaleID[!nests$MaleID %in% ped3_sex_update$ped$RCWid]
# nests$FemaleID[!nests$FemaleID %in% ped3_sex_update$ped$RCWid]
# 
# nests$MaleID
# male_reproductive_info <- lapply(unique(nests$MaleID), function(x) {
#   summarize_nest_offspring_info(nest_df = nests, 
#                                 rcw_df = rcws, 
#                                 bird_id = x)
# }) %>% 
#   bind_rows()
# 
# reproductive_summary_list <- lapply(list(Male = unique(nests$MaleID),
#                                          Female = unique(nests$FemaleID)), 
#                                     function(x) {
#                                       lapply(x, function(x) {
#                                         summarize_nest_offspring_info(nest_df = nests, 
#                                                                       rcw_df = rcws, 
#                                                                       bird_id = x)
#                                       }) %>% 
#                                         bind_rows()
#                                       
#                                     })
# 
# translocated_ids <- rcws[rcws$Translocation == 'Inter-population',]$RCWid
# translocated_descendants <- document_descendants(ped = rcw_ped, 
#                                                  focal_ancestors = translocated_ids) %>%
#   rowwise() %>% 
#   filter(isTRUE(focal_ancestor_descendant)) %>% 
#   ungroup() %>% 
#   pull(id)
# 
# 
# reproductive_summary_df <- bind_rows(reproductive_summary_list, .id = 'sex') %>% 
#   rename(RCWid = id) %>% 
#   left_join(., rcws[,c('RCWid', "Translocation")],
#             by = 'RCWid') %>% 
#   mutate(status = case_when(RCWid %in% translocated_ids ~ 'migrant', #init migrants
#                             RCWid %in% translocated_descendants ~ 'migrant descendant', #migrant descendants
#                             TRUE ~ 'non-migrant')) %>% 
#   filter(!is.na(RCWid))
# 
# fledgling_hist_by_sex <- reproductive_summary_df %>% 
#   mutate(status = factor(status,
#                          levels = c('non-migrant', 'migrant descendant', 'migrant'))) %>% 
#   ggplot() +
#   #geom_point(data = . %>% 
#   #             group_by(sex, status) %>% 
#   #             summarize(mean_count = mean(fledglings)),
#   #           aes(x = mean_count, y = 0 - 1, color = status), shape = 21, stroke = 2, size = 2) +
#   geom_segment(data = . %>% 
#                  group_by(sex, status) %>% 
#                  summarize(mean_count = mean(fledglings)),
#                aes(x = mean_count, y = 0, xend = mean_count, yend = 35, color = status),
#                linetype = 'dashed', size = 1.5) +
#   geom_point(data = . %>% 
#                group_by(sex, status) %>% 
#                summarize(mean_count = mean(fledglings)),
#              aes(x = mean_count, y = 35, color = status), shape = 21, stroke = 2, size = 3, fill = 'white') +
#   geom_histogram(aes(x = fledglings, fill = status), position = "stack", binwidth = 1) +
#   scale_fill_manual(values = c('#e4e4e4', '#b4b4b4', '#737373')) +
#   scale_color_manual(values = c('#e4e4e4', '#b4b4b4', '#737373')) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#         axis.title = element_text(size = 15),
#         plot.title = element_text(size = 20),
#         axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#         strip.background = element_rect(fill = "#737373"),
#         strip.text = element_text(color = "white", size = 12)) +
#   facet_wrap(~sex) +
#   xlab('Fledgling count') +
#   ylab('Count')
# 
# # ggsave(filename = here("figures", ".png"),
# #        plot = example_inbreeding_sim_parent_reassign,
# #        width = 8, height = 5)
# 
# 
# 
# reproductive_summary_df %>% 
#   mutate(status = factor(status,
#                          levels = c('non-migrant', 'migrant descendant', 'migrant'))) %>% 
#   group_by(sex, status) %>% 
#   summarize(mean_count = mean(fledglings)) 
# 
# surviving_indivs <- combined_census_indiv_mig_status %>% 
#   filter(year == max(year)) %>% 
#   pull(RCWid)
# 
# point_line_height <- rev(rep(c(26, 29, 32), times = 2))
# fledgling_hist_by_sex_remove_alive <- reproductive_summary_df %>% 
#   filter(!RCWid %in% surviving_indivs) %>% 
#   mutate(status = factor(status,
#                          levels = c('non-migrant', 'migrant descendant', 'migrant'))) %>% 
#   ggplot() +
#   geom_segment(data = . %>% 
#                  group_by(sex, status) %>% 
#                  summarize(mean_count = mean(fledglings)),
#                aes(x = mean_count, 
#                    y = 0, 
#                    xend = mean_count, 
#                    yend = point_line_height, 
#                    color = status),
#                linetype = 'dashed', size = 1.5) +
#   geom_point(data = . %>% 
#                group_by(sex, status) %>% 
#                summarize(mean_count = mean(fledglings)),
#              aes(x = mean_count, y = point_line_height, color = status), size = 5) +
#   geom_histogram(aes(x = fledglings, fill = status), position = "stack", binwidth = 1) +
#   scale_color_manual(values = c('#e4e4e4', '#b4b4b4', '#737373')) +
#   scale_fill_manual(values = c('#e4e4e4', '#b4b4b4', '#737373')) + #c('#e4e4e4', '#c9c9c9', '#808080')) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#         axis.title = element_text(size = 15),
#         plot.title = element_text(size = 20),
#         axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#         strip.background = element_rect(fill = "#737373"),
#         strip.text = element_text(color = "white", size = 12)) +
#   facet_wrap(~sex) +
#   xlab('Fledgling count') +
#   ylab('Count')
# 
# ggsave(filename = here("figures", "fledgling_hist_by_sex_remove_alive.png"),
#        plot = fledgling_hist_by_sex_remove_alive,
#        width = 8, height = 4)
# 
# 
# 
# 
# point_line_height <- rev(rep(c(65, 70, 75), times = 2))
# offspring_breeders_hist_by_sex_remove_alive <- reproductive_summary_df %>% 
#   filter(!RCWid %in% surviving_indivs) %>% 
#   mutate(status = factor(status,
#                          levels = c('non-migrant', 'migrant descendant', 'migrant'))) %>% 
#   ggplot() +
#   geom_segment(data = . %>% 
#                  group_by(sex, status) %>% 
#                  summarize(mean_count = mean(total_breeders)),
#                aes(x = mean_count, y = 0, xend = mean_count, yend = point_line_height, color = status),
#                linetype = 'dashed', size = 1.5) +
#   geom_point(data = . %>% 
#                group_by(sex, status) %>% 
#                summarize(mean_count = mean(total_breeders)),
#              aes(x = mean_count, y = point_line_height, color = status), size = 5) +
#   geom_histogram(aes(x = total_breeders, fill = status), position = "stack", binwidth = 1) +
#   scale_color_manual(values = c('#e4e4e4', '#b4b4b4', '#737373')) +
#   scale_fill_manual(values = c('#e4e4e4', '#b4b4b4', '#737373')) + #c('#e4e4e4', '#c9c9c9', '#808080')) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#         axis.title = element_text(size = 15),
#         plot.title = element_text(size = 20),
#         axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#         strip.background = element_rect(fill = "#737373"),
#         strip.text = element_text(color = "white", size = 12)) +
#   facet_wrap(~sex) +
#   xlab('Number of offspring breeders') +
#   ylab('Count')
# 
# ggsave(filename = here("figures", "offspring_breeders_hist_by_sex_remove_alive.png"),
#        plot = offspring_breeders_hist_by_sex_remove_alive,
#        width = 8, height = 4)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# reproductive_summary_df %>% 
#   filter(!RCWid %in% surviving_indivs) %>% 
#   arrange(total_breeding_years) %>% 
#   slice_head(n = 10)
# 
# point_line_height <- rev(rep(c(37, 40, 43), times = 2))
# reproductive_summary_df %>% 
#   filter(!RCWid %in% surviving_indivs) %>% 
#   mutate(status = factor(status,
#                          levels = c('non-migrant', 'migrant descendant', 'migrant'))) %>% 
#   ggplot() +
#   geom_segment(data = . %>% 
#                  group_by(sex, status) %>% 
#                  summarize(mean_count = mean(total_breeding_years)),
#                aes(x = mean_count, y = 0, xend = mean_count, yend = point_line_height, color = status),
#                linetype = 'dashed', size = 1.5) +
#   geom_point(data = . %>% 
#                group_by(sex, status) %>% 
#                summarize(mean_count = mean(total_breeding_years)),
#              aes(x = mean_count, y = point_line_height, color = status), size = 5) +
#   geom_histogram(aes(x = total_breeding_years, fill = status), position = "stack", binwidth = 1) +
#   scale_color_manual(values = c('#e4e4e4', '#b4b4b4', '#737373')) +
#   scale_fill_manual(values = c('#e4e4e4', '#b4b4b4', '#737373')) + #c('#e4e4e4', '#c9c9c9', '#808080')) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#         axis.title = element_text(size = 15),
#         plot.title = element_text(size = 20),
#         axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#         strip.background = element_rect(fill = "#737373"),
#         strip.text = element_text(color = "white", size = 12)) +
#   facet_wrap(~sex) +
#   xlab('Number of breeding years') +
#   ylab('Count')
# 
# 
# 
# 
# 
# 
# nests_processed <- nests %>%
#   left_join(., rcws[,c('RCWid', "Translocation")] %>% rename(MaleID = RCWid),
#             by = 'MaleID') %>% 
#   rename(male_translocation_info = Translocation) %>% 
#   left_join(., rcws[,c('RCWid', "Translocation")] %>% rename(FemaleID = RCWid),
#             by = 'FemaleID') %>% 
#   rename(female_translocation_info = Translocation) %>% 
#   mutate(male_status = case_when(MaleID %in% translocated_ids ~ 'migrant', #init migrants
#                             MaleID %in% translocated_descendants ~ 'migrant descendant', #migrant descendants
#                             TRUE ~ 'non-migrant'),
#          female_status = case_when(FemaleID %in% translocated_ids ~ 'migrant', #init migrants
#                                    FemaleID %in% translocated_descendants ~ 'migrant descendant', #migrant descendants
#                                    TRUE ~ 'non-migrant')) %>% 
#   filter(!is.na(MaleID) & !is.na(FemaleID)) %>% 
#   rowwise() %>% 
#   mutate(parental_combo1 = paste(sort(c(male_status, female_status)),collapse = '/'),
#          parental_combo2 = case_when(all(c(male_status, female_status) == 'migrant') ~ 'mig./mig.',
#                                      all(c(male_status, female_status) == 'non-migrant') ~ 'non-mig./non-mig.',
#                                      TRUE ~ 'mig. desc. mating'))
# 
# 
# nests_processed %>% 
#   pivot_longer(cols = c(EggNum, HatchNum, Helpers), names_to = "variable", values_to = "value") %>% 
#   filter(!is.na(value)) %>% 
#   ggplot() +
#   geom_jitter(aes(x = parental_combo2, y = value), 
#               position = position_jitter(height = 0.1),
#               size = 0.7, alpha = 0.8, color = 'gray') +
#   geom_boxplot(aes(x = parental_combo2, y = value), 
#                alpha = 0.5, color = '#B67380', fill = '#B67380',
#                outlier.shape = NA, width = 0.5) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#         axis.title = element_text(size = 15),
#         plot.title = element_text(size = 20),
#         axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#         strip.background = element_rect(fill = "#737373"),
#         strip.text = element_text(color = "white", size = 12)) +
#   facet_wrap(~variable, scales = 'free_y', nrow = 1) +
#   xlab('Parental origin') + ylab('Value')
# 
# 
# 
# 
# 
# 
# nests_processed %>% 
#   ggplot() +
#   geom_boxplot(aes(x = parental_combo2, y = HatchNum)) +
#   theme_bw()
# 
# nests_processed %>% 
#   filter(!is.na(EggNum)) %>%
#   ggplot() +
#   geom_jitter(aes(x = parental_combo2, y = EggNum), 
#               position = position_jitter(height = 0.1),
#               size = 1, alpha = 0.4) +
#   geom_boxplot(aes(x = parental_combo2, y = EggNum), alpha = 0.5) +
#   theme_bw() +
#   theme(panel.grid.minor.y = element_blank()) +
#   xlab('Parent migration status') + 
#   ylab('Egg Count')
# 
# 
# 
# plot_list_sex <- lapply(setNames(nm = c('male', 'female')), function(sex) {
#   nests_processed %>% 
#     # group_by(!!! rlang::syms(capitalize_first(paste0(sex, 'ID')))) %>% 
#     #  summarize(nest_count = n(),
#     #            total_egg_count = sum(EggNum, na.rm = TRUE),
#     #            total_hatch_count = sum(HatchNum, na.rm = TRUE),
#     #            mean_helper = mean(Helpers, na.rm = TRUE),
#     #            origin = first(!!! rlang::syms(paste0(sex, '_status')))) %>%
#     pivot_longer(cols = c(EggNum, HatchNum, Helpers), 
#                  names_to = "variable", 
#                  values_to = "value") %>%
#     filter(!is.na(value)) %>%
#     ggplot() +
#     geom_jitter(aes_string(x = paste0(sex, '_status'), y = 'value'), 
#                 position = position_jitter(height = 0.1),
#                 size = 0.7, alpha = 0.8) +
#     geom_boxplot(aes_string(x = paste0(sex, '_status'), y = 'value'), 
#                  alpha = 0.5, color = '#B67380', fill = '#B67380',
#                  outlier.shape = NA, width = 0.5) +
#     facet_wrap(~variable, scales = 'free_y', nrow = 1) +
#     ggtitle(capitalize_first(sex)) +
#     xlab('Origin') + ylab('Value') +
#     theme_bw() +
#     theme(legend.title = element_blank(),
#           panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#           panel.grid.minor = element_blank(),
#           panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#           axis.title = element_text(size = 15),
#           plot.title = element_text(size = 20),
#           axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#           strip.background = element_rect(fill = "#737373"),
#           strip.text = element_text(color = "white", size = 12))
# })
# 
# 
# 
# nests_processed %>% 
#   group_by(!!! rlang::syms(capitalize_first(paste0('male', 'ID')))) %>% 
#   summarize(nest_count = n(),
#             total_egg_count = sum(EggNum, na.rm = TRUE),
#             total_hatch_count = sum(HatchNum, na.rm = TRUE),
#             mean_helper = mean(Helpers, na.rm = TRUE),
#             origin = first(!!! rlang::syms(paste0('male', '_status')))) %>%
#   pivot_longer(cols = c(total_egg_count, total_hatch_count, mean_helper), 
#                names_to = "variable", 
#                values_to = "value") %>%
#   filter(!is.na(value))
# 
# 
# 
# nests_processed %>% 
#   filter((!!! rlang::syms(capitalize_first(paste0('male', 'ID'))) %in% surviving_indivs))
# 
# nests_processed %>% 
#   filter(!FemaleID %in% surviving_indivs)
# 
# 
# plot_list_sex <- lapply(setNames(nm = c('male', 'female')), function(sex) {
#   
#   nest_processed_remove_surviving <- nests_processed[!nests_processed[,capitalize_first(paste0(sex, 'ID')),drop = TRUE] %in% surviving_indivs,]
#   
#   nest_processed_remove_surviving %>% 
#     filter(!(!!! rlang::syms(capitalize_first(paste0(sex, 'ID'))) %in% surviving_indivs)) %>% 
#     group_by(!!! rlang::syms(capitalize_first(paste0(sex, 'ID')))) %>%
#      summarize(nest_count = n(),
#                total_egg_count = sum(EggNum, na.rm = TRUE),
#                total_hatch_count = sum(HatchNum, na.rm = TRUE),
#                mean_helper = mean(Helpers, na.rm = TRUE),
#                origin = first(!!! rlang::syms(paste0(sex, '_status')))) %>%
#     pivot_longer(cols = c(total_egg_count, total_hatch_count, mean_helper), 
#                  names_to = "variable", 
#                  values_to = "value") %>%
#     filter(!is.na(value)) %>%
#     ggplot() +
#     geom_jitter(aes_string(x = 'origin', y = 'value'), 
#                 position = position_jitter(height = 0.1),
#                 size = 0.7, alpha = 0.8) +
#     geom_boxplot(aes_string(x = 'origin', y = 'value'), 
#                  alpha = 0.5, color = '#B67380', fill = '#B67380',
#                  outlier.shape = NA, width = 0.5) +
#     facet_wrap(~variable, scales = 'free_y', nrow = 1) +
#     ggtitle(capitalize_first(sex)) +
#     xlab('Origin') + ylab('Value') +
#     theme_bw() +
#     theme(legend.title = element_blank(),
#           panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#           panel.grid.minor = element_blank(),
#           panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#           axis.title = element_text(size = 15),
#           plot.title = element_text(size = 20),
#           axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#           strip.background = element_rect(fill = "#737373"),
#           strip.text = element_text(color = "white", size = 12))
# })
# 
# 
# 
# cowplot::plot_grid(plot_list_sex$male,
#                    plot_list_sex$female,
#                    nrow = 2)

# prac_census_year <- data.frame(
#   year = rep(1:10, times = 3),
#   id = paste0('id_', rep(1:3, each = 10)),
#   extra_col = sample(x = 1:400, 30, replace = TRUE)
# )
# 
# 
# prac_census_year_subset <- prac_census_year[-sample(1:nrow(prac_census_year), 10),]
# 
# add_intervening_years(prac_census_year_subset, 'id', 'year') %>%
#   arrange(id, year)
# 
# 
# 
# cor(summer_census_add_int %>% 
#       group_by(year) %>%
#       summarize(count = n(), .groups = 'drop') %>% pull(count),
#     last_seen_summer_census %>% 
#       group_by(year) %>%
#       summarize(count = n(), .groups = 'drop') %>% pull(count))
# 
# data.frame(interp_census_size = summer_census_add_int %>% 
#         group_by(year) %>%
#         summarize(count = n(), .groups = 'drop') %>% pull(count),
#         no_interp_census_size = last_seen_summer_census %>% 
#         group_by(year) %>%
#         summarize(count = n(), .groups = 'drop') %>% pull(count)) %>% 
#   mutate(year = 1:29) %>% 
#   ggplot() +
#   geom_point(aes(x = year, y = no_interp_census_size), color = 'blue') +
#   geom_point(aes(x = year, y = interp_census_size), color = 'gray') +
#   ylim(0, 190) +
#   theme_bw()
# 
# 
# summer_census_add_int %>% 
#   group_by(year) %>%
#   summarize(count = n(), .groups = 'drop') %>% 
#   ggplot() +
#   geom_point(aes(x = year, y = count)) +
#   theme_bw()
# 
# last_seen_summer_census %>% 
#   group_by(year) %>%
#   summarize(count = n(), .groups = 'drop') %>% 
#   ggplot() +
#   geom_point(aes(x = year, y = count)) +
#   theme_bw()
# 
# 
# #figure out 
# 
# 
# nests %>% 
#   select(Year, MaleID, FemaleID) %>% 
#   pivot_longer(cols = !Year, names_to = 'sex', values_to = 'RCWid') %>% 
#   left_join(., rcws[,c('RCWid', "Translocation")],
#             by = 'RCWid') %>% 
#   group_by(Translocation) %>% 
#   slice_head(n = 1)
# 
# summer_census_interp_idupdate
# 
# nests_sim_info <- nests %>% 
#   select(Year, MaleID, FemaleID) %>% 
#   pivot_longer(cols = !Year, names_to = 'sex', values_to = 'RCWid') %>% 
#   left_join(., rcws[,c('RCWid', "Translocation")],
#             by = 'RCWid') %>% 
#   mutate(status = case_when(Translocation == 'Inter-population' ~ 'migrant',
#                             TRUE ~ 'non-migrant'),
#          sex = if_else(sex == 'MaleID', 'M', 'F'))
# 
# 
# nests_sim_info
# 
#   
# census1 <- summer_census_interp_idupdate %>% 
#   left_join(., rcws[,c('RCWid', "Translocation", 'MinAge')],
#             by = 'RCWid') %>% 
#   mutate(status = case_when(Translocation == 'Inter-population' ~ 'migrant',
#                             TRUE ~ 'non-migrant')) %>% 
#   filter(status != 'migrant' & (year != MinAge) )
#   
# 
# potential_breeders_list <- list()
# loop_index <- 1
# for (i in sort(unique(census1$year))) {
#   
#   existing_breeders_yeari <- nests_sim_info %>% 
#     filter(Year == i) %>% 
#     pull(RCWid)
#   
#   potential_breeders_list[[loop_index]] <- census1 %>% 
#     filter(year == i) %>% 
#     filter(!RCWid %in% existing_breeders_yeari)
#   
#   loop_index <- loop_index + 1
# }
# 
# potential_breeders_df1 <- bind_rows(potential_breeders_list)
# potential_breeders_df <- left_join(potential_breeders_df1, 
#                                    rcws[,c('RCWid', 'Sex')],
#                                    by = 'RCWid')
# 
# 
# ped_sim_annotated <- left_join(ped3_sex_update$ped, 
#           rcws[,c('RCWid', "Translocation")],
#           by = 'RCWid') %>% 
#   rename(focal_translocation_info = Translocation) %>%
#   left_join(., rcws[,c('RCWid', "Translocation")] %>% rename(MaleID = RCWid),
#                                                              by = 'MaleID') %>% 
#   rename(male_translocation_info = Translocation) %>% 
#   left_join(., rcws[,c('RCWid', "Translocation")] %>% rename(FemaleID = RCWid),
#             by = 'FemaleID') %>% 
#   rename(female_translocation_info = Translocation) %>% 
#   mutate(focal_status = case_when(focal_translocation_info == 'Inter-population' ~ 'migrant',
#                                   TRUE ~ 'non-migrant'),
#          male_status = case_when(male_translocation_info == 'Inter-population' ~ 'migrant',
#                                   TRUE ~ 'non-migrant'),
#          female_status = case_when(female_translocation_info == 'Inter-population' ~ 'migrant',
#                                   TRUE ~ 'non-migrant'))
# 
# 
# 
# replace_migrant_parents <- function(ped, 
#                                     potential_breeders_df,
#                                     migrant_ids,
#                                     iters = 10) {
#   
#   ped_list <- list()
#   
#   male_parent_migrants_index_df <- 
#     data.frame(index = which(ped$male_status == 'migrant'),
#                year = ped$Year[ped$male_status == 'migrant'])
#   
#   female_parent_migrants_index_df <- 
#     data.frame(index = which(ped$female_status == 'migrant'),
#                year = ped$Year[ped$female_status == 'migrant'])
#   
#   for (i in seq_len(iters)) {
#     
#     ped[,paste0('sim', i,'_MaleID')] <- ped[,'MaleID',drop = TRUE]
#     ped[,paste0('sim', i,'_FemaleID')] <- ped[,'FemaleID',drop = TRUE]
#     
#     ### MALE REPLACEMENT ###
#     for (YEAR in unique(male_parent_migrants_index_df$year)) {
#       year_male_subset <- male_parent_migrants_index_df %>% 
#         filter(year == YEAR)#filter(year == YEAR)
#       
#       potential_male_breeders <- potential_breeders_df %>% 
#         filter(year == YEAR & Sex == 'M') %>% 
#         pull(RCWid)
#       
#       if (nrow(year_male_subset) > length(potential_male_breeders))
#         message('not enough potential male breeders to sample without replacement')
#       
#       year_male_subset$parent_replacement <- sample(potential_male_breeders, 
#                                                     size = nrow(year_male_subset), 
#                                                     replace = TRUE)
#       
#       ped[year_male_subset$index,paste0('sim', i,'_MaleID')] <- year_male_subset$parent_replacement
#       
#     }
#     
#     ### FEMALE REPLACEMENT ###
#     for (YEAR in unique(female_parent_migrants_index_df$year)) {
#       year_female_subset <- female_parent_migrants_index_df %>% 
#         filter(year == YEAR)#filter(year == YEAR)
#       
#       
#       potential_female_breeders <- potential_breeders_df %>% 
#         filter(year == YEAR & Sex == 'F') %>% 
#         pull(RCWid)
#       
#       #message('number of female migrants: ', nrow(year_female_subset))
#       #message('number of potential replacements: ', length(potential_female_breeders))
#       
#       #if (nrow(year_female_subset) > length(potential_female_breeders))
#       #  message('not enough potential female breeders to sample without replacement')
#       
#       year_female_subset$parent_replacement <- sample(potential_female_breeders, 
#                                                       size = nrow(year_female_subset), 
#                                                       replace = TRUE)
#       
#       ped[year_female_subset$index,paste0('sim', i,'_FemaleID')] <- year_female_subset$parent_replacement
#       
#     }
#     
#     ped_remove_migrants <- ped %>% 
#       filter(!RCWid %in% migrant_ids)
#     ped_list[[i]] <- ped(id = ped_remove_migrants$RCWid, 
#                          fid = ped_remove_migrants[,paste0('sim', i,'_MaleID'),drop = TRUE], 
#                          mid = ped_remove_migrants[,paste0('sim', i,'_FemaleID'),drop = TRUE], 
#                          sex = ped_remove_migrants$sex_corrected)
#     
#     message('finished ', i)
#   }
#   
#   return(list(ped_df = ped,
#          ped_list = ped_list)
#          )
# }
# 
# 
# replace_migrant_parents_output <- replace_migrant_parents(ped = ped_sim_annotated, 
#                                                           potential_breeders_df = potential_breeders_df, 
#                                                           migrant_ids = translocated_ids,
#                                                           iters = 500)
# 
# inbr_calc_list <- list()
# for (SIM in 1:length(replace_migrant_parents_output$ped_list)) {
#   inbr_calc_list[[SIM]] <- lapply(seq(0, 0.03, by = 0.01), function(inbr_val) {
#     inbr_df <- stack_v2(
#       inbreeding_wrapper(ped = replace_migrant_parents_output$ped_list[[SIM]],
#                          clear_founder_inbreeding = TRUE,
#                          founder_inbreeding_val = inbr_val,
#                          founder_ids = NULL,
#                          inbreeding_ids = NULL, 
#                          Xchrom = FALSE)
#     ) %>% 
#       mutate(founder_fped_val = inbr_val)
#     rownames(inbr_df) <- NULL
#     return(inbr_df)
#   }) %>% 
#     bind_rows() %>% 
#     rename(RCWid = ind,
#            fped = values) %>% 
#     pivot_wider(names_from = founder_fped_val, 
#                 values_from = fped, 
#                 names_prefix = 'founder_') %>% 
#     mutate(sim = paste0('sim_', SIM))
#   message(SIM)
# }
# 
# 
# inbr_combined_info_list <- lapply(inbr_calc_list, function(INBR_DF) {
#   left_join(summer_census_interp_idupdate,
#             rcws[,c('RCWid', "Translocation")],
#             by = 'RCWid') %>% 
#     left_join(., INBR_DF,
#               by = 'RCWid')
# })
# 
# inbr_combined_info_df <- bind_rows(inbr_combined_info_list)
# 
# 
# inbr_combined_info_df_remove_migrants <- inbr_combined_info_df %>%
#   filter(!RCWid %in% translocated_ids)
# 
# 
# inbr_combined_info_df_remove_migrants
# 
# 
# combined_census_indiv_info <- left_join(summer_census_interp_idupdate,
#                                         rcws[,c('RCWid', "Translocation")],
#                                         by = 'RCWid') %>% 
#   left_join(., inbr_calcs,
#             by = 'RCWid') %>% 
#   select(ID, RCWid, RecordDate, Detected, RoostConfirmed, Comments,
#          SurveyType, year, RCWid_old, Translocation, founder_0, founder_0.01, 
#          founder_0.02, founder_0.03) %>% 
#   mutate(sim = 'observed')
# 
# inbr_combined_info_df_remove_migrants_incl_observed <- rbind(
#   inbr_combined_info_df_remove_migrants,
#   combined_census_indiv_info %>% filter(!RCWid %in% translocated_ids)
# )
# 
# 
# sim_count <- 500
# 
# example_inbreeding_sim_parent_reassign <- inbr_combined_info_df_remove_migrants_incl_observed %>% 
#   #filter(sim %in% c('observed', paste0('sim_', 1:sim_count))) %>% 
#   mutate(sim = factor(sim, levels = c(paste0('sim_', 1:sim_count), 'observed'))) %>% 
#   pivot_longer(cols = starts_with('founder'), names_to = 'Founder_inbreeding_val', values_to = 'f_ped') %>% 
#   group_by(sim, year, Founder_inbreeding_val) %>% 
#   summarize(sim = first(sim),
#             mean_fped = mean(f_ped),
#             lower_quantile = quantile(f_ped, probs = 0.025),
#             upper_quantile = quantile(f_ped, probs = 0.975)) %>%
#   ungroup() %>% 
#   ggplot(aes(x = year, y = mean_fped, color = sim, size = sim, alpha = sim)) +
#   geom_hline(yintercept = 0, color = '#999999', size = 0.5) +
#   geom_line(alpha = 0.5, lineend = 'round') + #color = "#03045e"
#   geom_point(alpha = 0.6) +
#   scale_color_manual(values = c(rep('gray', sim_count), '#3a0ca3') ) + 
#   scale_size_manual(values = c(rep(0.4, sim_count), 3.5) ) + 
#   scale_alpha_manual(values = c(rep(0.2, sim_count), 0.8) ) + 
#     theme_bw() +
#     theme(#panel.grid.major = element_line(size = 0.35, color = '#f2f2f2'),
#       #panel.grid.minor = element_line(size = 0.35, linetype = 'dashed', color = '#f2f2f2'),
#       panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#       panel.grid.minor = element_blank(),
#       panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#       axis.title = element_text(size = 15),
#       plot.title = element_text(size = 20),
#       axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#       strip.background = element_rect(fill = "#737373"),
#       strip.text = element_text(color = "white", size = 12),
#       legend.position = 'none') +
#     facet_wrap(~Founder_inbreeding_val, scales = "free_y") +
#     xlab('Year') +
#     ylab('Mean inbreeding')
#   
# ggsave(filename = here("figures", "example_inbreeding_sim_parent_reassign.png"),
#        plot = example_inbreeding_sim_parent_reassign,
#        width = 8, height = 5)
#   
# 
# inbreeding_facet_differentfounder <- combined_census_indiv_info %>% 
#   pivot_longer(cols = starts_with('founder'), names_to = 'Founder_inbreeding_val', values_to = 'f_ped') %>% 
#   group_by(year, Founder_inbreeding_val) %>% 
#   summarize(mean_fped = mean(f_ped),
#             lower_quantile = quantile(f_ped, probs = 0.025),
#             upper_quantile = quantile(f_ped, probs = 0.975)) %>%
#   ungroup() %>% 
#   ggplot(aes(x = year, y = mean_fped)) +
#   geom_hline(yintercept = 0, color = '#999999', size = 0.5) +
#   geom_line(alpha = 0.5, size = 1, color = "#03045e") +
#   geom_point(alpha = 0.6, size = 3, color = "#03045e") +
#   theme_bw() +
#   theme(#panel.grid.major = element_line(size = 0.35, color = '#f2f2f2'),
#     #panel.grid.minor = element_line(size = 0.35, linetype = 'dashed', color = '#f2f2f2'),
#     panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = '#f2f2f2'),
#     panel.grid.minor = element_blank(),
#     panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#     axis.title = element_text(size = 15),
#     plot.title = element_text(size = 20),
#     axis.line = element_line(size = 1, lineend = 'round', color = '#737373'),
#     strip.background = element_rect(fill = "#737373"),
#     strip.text = element_text(color = "white", size = 12)) +
#   facet_wrap(~Founder_inbreeding_val, scales = "free_y") +
#   xlab('Year') +
#   ylab('Mean inbreeding')
# 


# 
# male_parent_migrants_index_df <- 
#   data.frame(index = which(ped_sim_annotated$male_status == 'migrant'),
#              year = ped_sim_annotated$Year[ped_sim_annotated$male_status == 'migrant'])
# 
# female_parent_migrants_index_df <- 
#   data.frame(index = which(ped_sim_annotated$female_status == 'migrant'),
#              year = ped_sim_annotated$Year[ped_sim_annotated$female_status == 'migrant'])
# 
# 
# ped_sim_annotated[,'sim1_MaleID'] <- ped_sim_annotated[,'MaleID',drop = TRUE]
# ped_sim_annotated[,'sim1_FemaleID'] <- ped_sim_annotated[,'FemaleID',drop = TRUE]
# 
# ### MALE REPLACEMENT ###
# for (YEAR in unique(male_parent_migrants_index_df$year)) {
#   year_male_subset <- male_parent_migrants_index_df %>% 
#     filter(year == YEAR)#filter(year == YEAR)
#   
#   
#   potential_male_breeders <- potential_breeders_df %>% 
#     filter(year == YEAR & Sex == 'M') %>% 
#     pull(RCWid)
#   
#   if (nrow(year_male_subset) > length(potential_male_breeders))
#     message('not enough potential male breeders to sample without replacement')
#   
#   year_male_subset$parent_replacement <- sample(potential_male_breeders, 
#                                                 size = nrow(year_male_subset), 
#                                                 replace = TRUE)
#   
#   ped_sim_annotated[year_male_subset$index,'sim1_MaleID'] <- year_male_subset$parent_replacement
#   
# }
# 
# ### FEMALE REPLACEMENT ###
# for (YEAR in unique(female_parent_migrants_index_df$year)) {
#   year_female_subset <- female_parent_migrants_index_df %>% 
#     filter(year == YEAR)#filter(year == YEAR)
#   
#   
#   potential_female_breeders <- potential_breeders_df %>% 
#     filter(year == YEAR & Sex == 'F') %>% 
#     pull(RCWid)
#   
#   #message('number of female migrants: ', nrow(year_female_subset))
#   #message('number of potential replacements: ', length(potential_female_breeders))
#   
#   #if (nrow(year_female_subset) > length(potential_female_breeders))
#   #  message('not enough potential female breeders to sample without replacement')
#   
#   year_female_subset$parent_replacement <- sample(potential_female_breeders, 
#                                                 size = nrow(year_female_subset), 
#                                                 replace = TRUE)
#   
#   ped_sim_annotated[year_female_subset$index,'sim1_FemaleID'] <- year_female_subset$parent_replacement
#   
# }
# 
# 
# new_ped <- ped(id = ped_sim_annotated$RCWid, 
#                fid = ped_sim_annotated[,'sim1_MaleID',drop = TRUE], 
#                mid = ped_sim_annotated[,'sim1_FemaleID',drop = TRUE], 
#                sex = ped_sim_annotated$sex_corrected)
# 
# 
# 
# prac_df <- data.frame(a = 1:10,
#                       b = 1:10,
#                       c = letters[1:10])
# 
# prac_df[c(2, 5, 6, 8),'b'] <- c(200, 500, 600, 800)
# 
# 
# ped_sim_annotated
# 
# which(ped_sim_annotated$male_status == 'migrant')
# which(ped_sim_annotated$female_status == 'migrant')
# 
# 
# 
# 
# table(nests_processed$Year)
# 
# fledging_count_multipanel <- nests_processed %>% 
#   pivot_longer(cols = c(MaleID, FemaleID), names_to = 'parent', values_to = 'parent_id') %>% 
#   mutate(class = if_else(parent == 'MaleID', male_status, female_status)) %>% 
#   mutate(parent_sex = if_else(parent == 'MaleID', 'Male', 'Female') ) %>% 
#   group_by(Year, class, parent_sex) %>% 
#   summarize(total_fledgling_count = sum(FldgNum, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   ggplot() +
#   geom_bar(aes(x = Year, y = total_fledgling_count, fill = class), 
#            #color = 'gray', size = 0.2,
#            position = "stack", stat = "identity", width = 0.98) +
#   scale_fill_manual(values = rev(c('#e4e4e4', '#b4b4b4', '#737373'))) +
#   theme_bw() +
#   scale_x_continuous(expand = expansion(mult = c(0, 0), add = c(0, 1))) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 1))) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "gray", fill = NA, size = 0.5),
#         axis.line.x = element_line(color = '#4c4c4c', size = 1.1),
#         axis.line.y = element_line(color = '#4c4c4c', size = 1.5),
#         axis.ticks = element_line(color = '#808080', size = 0.7),
#         axis.title = element_text(size = 14),
#         strip.background = element_rect(fill = "#737373"),
#         strip.text = element_text(color = "white", size = 12)) +
#   theme(plot.margin=unit(c(1,1,1,1), "cm")) +
#   xlab('Year') +
#   ylab('Number of fledglings') +
#   facet_wrap(~parent_sex)
# 
# ggsave(filename = here("figures", "fledging_count_multipanel.png"),
#        plot = fledging_count_multipanel,
#        width = 10, height = 5)
#   
#   
# nests_processed %>% 
#   pivot_longer(cols = c(MaleID, FemaleID), names_to = 'parent', values_to = 'parent_id') %>% 
#   mutate(class = if_else(parent == 'MaleID', male_status, female_status)) %>% 
#   group_by(Year, class, parent) %>% 
#   summarize(total_fledgling_count = sum(FldgNum, na.rm = TRUE))  %>% 
#   print(n = 150)
# 
# 
# 
# 
# 
# 
# 
# 
# candidate_breeders <- summer_census_interp_idupdate %>% 
#   left_join(., rcws[,c('RCWid', 'Sex', 'MinAge', 'Translocation')],
#             by = 'RCWid') %>% 
#   filter(year != MinAge)
# 
# nests_long <- pivot_longer(nests, 
#                            cols = c('MaleID', 'FemaleID'),
#                            names_to = 'Parent',
#                            values_to = 'RCWid') %>% 
#   select(RCWid, Parent, Year, EggNum, HatchNum, FldgNum)
# 
# breeders_df_list <- list()
# index <- 1
# for (i in unique(nests_long$Year)) {
#   year_subset <- nests_long %>% 
#     filter(Year == i)
#   
#   candidate_breeders_subset <- candidate_breeders %>% 
#     filter(year == i) %>% 
#     filter(!RCWid %in% year_subset$RCWid) %>% 
#     select(RCWid, Sex)
#   
#   if (nrow(candidate_breeders_subset) > 0) {
#     breeders_df_list[[index]] <- rbind(
#       year_subset,
#       data.frame(RCWid = candidate_breeders_subset$RCWid,
#                  Parent = ifelse(candidate_breeders_subset$Sex == 'M', 'MaleID', 'FemaleID'),
#                  Year = i,
#                  EggNum = 0,
#                  HatchNum = 0,
#                  FldgNum = 0)
#     )
#   } else {
#     breeders_df_list[[index]] <- year_subset
#   }
#   
#   index <- index + 1
# }
# 
# 
# breeders_df_addnonbreeders <- bind_rows(breeders_df_list)
# 
# 
# lifetime_breed_info <- breeders_df_addnonbreeders %>% 
#   group_by(RCWid) %>% 
#   summarize(Parent = first(Parent),
#             EggNum = sum(EggNum, na.rm = TRUE),
#             HatchNum = sum(HatchNum, na.rm = TRUE),
#             FldgNum = sum(FldgNum, na.rm = TRUE)) %>% 
#   left_join(., rcws[, c('RCWid', 'Translocation')], by = 'RCWid') %>% 
#   mutate(mig_status = case_when(Translocation == 'Inter-population' ~ 'migrant',
#                                 TRUE ~ 'else'))
# 
# lifetime_breed_info %>% 
#   filter(Parent == 'MaleID') %>% 
#   pull(FldgNum)
# 
# M_index(lifetime_breed_info %>% 
#           filter(Parent == 'MaleID') %>% 
#           pull(FldgNum),
#         rep(1, nrow(lifetime_breed_info %>% 
#                         filter(Parent == 'MaleID'))) )
# 
# M_index(lifetime_breed_info %>% 
#           filter(Parent == 'MaleID' & mig_status != 'migrant') %>% 
#           pull(FldgNum),
#         rep(1, nrow(lifetime_breed_info %>% 
#                       filter(Parent == 'MaleID' & mig_status != 'migrant'))) )
# 
# M_index(lifetime_breed_info %>% 
#           filter(Parent == 'FemaleID') %>% 
#           pull(HatchNum),
#         rep(1, nrow(lifetime_breed_info %>% 
#                       filter(Parent == 'FemaleID'))) )
# 
# 
# lifetime_breed_info %>% 
#   filter(Parent == 'FemaleID') %>% 
#   pull(FldgNum)
# 
# 
# lifetime_breed_info %>% 
#   ggplot() +
#   geom_histogram(aes(x = EggNum)) +
#   facet_wrap(~Parent)
# 
# 
# multinom_index_df <- data.frame(
#   Year = sort(unique(breeders_df_addnonbreeders$Year)),
#   male = NA,
#   female = NA
# )
# 
# 
# 
# for (i in multinom_index_df$Year) {
#   
#   male_subset <- breeders_df_addnonbreeders %>% 
#     filter(Parent == 'MaleID' & Year == i) %>% 
#     filter(!is.na(FldgNum))
#   
#   female_subset <- breeders_df_addnonbreeders %>% 
#     filter(Parent == 'FemaleID' & Year == i) %>% 
#     filter(!is.na(FldgNum))
#   
#   multinom_index_df[multinom_index_df$Year == i, 'male'] <- 
#     M_index(male_subset$FldgNum,
#             rep(1, length(male_subset$FldgNum)) )
#   
#   multinom_index_df[multinom_index_df$Year == i, 'female'] <- 
#     M_index(female_subset$FldgNum,
#             rep(1, length(female_subset$FldgNum)) )
#   message(i)
#   
# }
# 
# table(breeders_df_addnonbreeders$Parent,
#       breeders_df_addnonbreeders$Year)
# 
# multinom_index_df %>%
#   filter(Year != 2022) %>% 
#   ggplot() +
#   geom_point(aes(x = Year, y = male)) +
#   geom_hline(yintercept = 0, color = 'black') +
#   theme_bw()
# 
# multinom_index_df %>%
#   filter(Year != 2022) %>% 
#   ggplot() +
#   geom_point(aes(x = Year, y = female)) +
#   geom_hline(yintercept = 0, color = 'black') +
#   theme_bw()