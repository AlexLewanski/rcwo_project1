#####################
### SCRIPT SET-UP ###
#####################

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
                                "ne_f_df", "rcw_inbr_merge", 
                                "rcws_ancestry_info",
                                "rcws_ancestry_source_by_year",
                                "rcw_partial_founder_summed_processed",
                                "rcw_contr_info_df",
                                "fped_prop_by_group",
                                "inbr_founder_color")
                         )


results_list <- lapply(results_names, function(x) read.csv(here('results', paste0(x, '.csv'))))



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

rcw_inbr_merge_processed <- results_list$rcw_inbr_merge %>% 
  group_by(year) %>% 
  summarize(perc_10 = quantile(f_ped, prob = 0.05),
            mean = mean(f_ped),
            perc_90 = quantile(f_ped, prob = 0.95),
            .groups = 'drop') %>% 
  arrange(year) %>% 
  mutate(change_color = case_when(mean > lag(mean) ~ '#ff4c4c',
                                  mean < lag(mean) ~ '#329932',
                                  TRUE ~ '#4c4c4c'))

fped_plot <- results_list$rcw_inbr_merge %>% 
  ggplot() +
  geom_point(aes(x = year, y = f_ped),
             shape = 21, colour = '#b7b7b7', fill = '#e3e3e3', 
             size = 1, alpha = 0.8,
             position = position_jitter(w = 0.16, h = 0)) +
  geom_segment(data = rcw_inbr_merge_processed,
               aes(x = year, xend = year, 
                   y = perc_10, yend = perc_90),
               colour = '#4c4c4c', linewidth = 1.15) +
  geom_point(data = rcw_inbr_merge_processed,
             aes(x = year, y = mean),
             shape = 21, color = rcw_inbr_merge_processed$change_color, 
             fill = 'white', 
             size = 2, stroke = 1.4) +
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
        axis.title.y = element_text(margin = margin(0, 5.5, 0, 0))
  ) +
  scale_x_continuous(limits = c(1993.5, 2022.5), position = "top",
                     expand = c(0.01, 0)) +
  ylab(expression(italic(F)["p"])) +
  coord_cartesian(clip = 'off')


top_bar <- results_list$rcw_inbr_merge %>%
  mutate(count = 1) %>%
  group_by(year) %>%
  arrange(f_ped) %>%
  ggplot() +
  geom_bar(aes(x = year, y = count, fill = f_ped),
           position = "stack", stat = "identity", width = 0.98) +
  scale_fill_gradient(name = expression(italic(F)["p"]), 
                      low = '#f4eaf4', high = '#730073') +
  geom_point(shape = "l", 
             #shape = 108, 
             data = results_list$transloc_info_color,
             aes(x = year, y = 0, color = RCWid), size = 5.5,
             alpha = 1,
             position = position_jitter(width = 0.44, height = 0, seed = 195),
             show.legend = FALSE) +
  scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
                                       nm = results_list$transloc_info_color$RCWid)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 1))) +
  theme_bw() +
  theme(legend.title = element_text(size = 10),
        panel.grid.major.y = element_blank(),
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

bottom_bar <- results_list$rcw_inbr_merge %>%
  mutate(count = 1) %>%
  group_by(year) %>%
  arrange(anc_prop) %>%
  ggplot() +
  geom_bar(aes(x = year, y = count, fill = anc_prop),
           position = "stack", stat = "identity", width = 0.98) +
  scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = '#ededed', high = '#4c4c4c') +
  geom_point(shape = "l", data = results_list$transloc_info_color,
             aes(x = year, y = 0, color = RCWid), size = 5.5,
             alpha = 1,
             position = position_jitter(width = 0.44, height = 0, seed = 195),
             show.legend = FALSE) +
  scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
                                       nm = results_list$transloc_info_color$RCWid)) +
  theme_bw() +
  theme(legend.title = element_text(size = 10),
        panel.grid.major.y = element_blank(),
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


transloc_ancestry <- results_list$rcws_ancestry_source_by_year %>%
  mutate(group = factor(group, 
                levels = rev(c('ANF', 'CBJTC', 'FTB', 'FTS', 'ONF', 'WSF-CITRUS', 'non-transloc')))) %>% 
  ggplot() +
  geom_bar(aes(x = year, y = ancestry_prop, fill = group), 
           color = 'white', linewidth = 0.15, width = 0.72,
           position = "stack", stat = "identity") + 
  #scale_fill_manual(values = c('#ededed', '#4c4c4c')) +
  scale_fill_manual(values = group_color_vec) +
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
        legend.title = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0))
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0.0", "0.5", "1.0")) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
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
               size = 2, stroke = 1.4) +
    theme_bw() +
    theme(legend.key = element_blank(),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent", color = "transparent")) +
    scale_color_manual(substitute(paste('Change in mean(', italic(F)["p"], ")")),
                       values = c('#ff4c4c', '#4c4c4c', '#329932'))
  
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
            size = 2.8,
            hjust = 0) +
  theme_void() +
  theme(legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  xlim(-0.32, 4) +
  geom_text(aes(x = 2.76, y = 207, label = "Founder proportions"),
            size = 3.5)


### add in legends and RCW image ###
rcw_pop_overview <- ggpubr::annotate_figure(egg::ggarrange(fped_plot,
                                                           top_bar,
                                                           bottom_bar,
                                                           transloc_ancestry,
                                                           nrow = 4,
                                                           labels = c('A', 'B', '', 'C'),
                                                           label.args = list(gp = grid::gpar(font = 2, cex = 1.2)),
                                                           heights = c(30, 50, 50, 22)),
                                            left = grid::textGrob("Individuals", rot = 90, hjust = 0.58, vjust = 5.2, gp = grid::gpar(cex = 1))) +
  theme(plot.margin = unit(c(0.3, -0.8, 0, -3), "cm"))

rcw_pop_overview_rcw <- rcw_pop_overview +
  cowplot::draw_image(flying_rcw,
                      scale = 0.12,
                      x = 0.325, y = 0.432) +
  patchwork::inset_element(founder_prop_legend,
                           left = 0.849, bottom = -0.01, right = 0.98, top = 0.20
                           #left = 0.849, bottom = -0.01, right = 0.972, top = 0.20
                           ) +
  patchwork::inset_element(fped_plot_legend,
                           #left = 0, bottom = 0.9, right = 0.38, top = 0.95,
                           left = 0.891, bottom = 0.81, right = 0.94, top = 0.88)


if (isTRUE(output_figs)) {
  cowplot::ggsave2(filename = here('figures', 'main_paper', 'rcw_pop_overview_fig.png'),
                   plot = ggdraw(rcw_pop_overview_rcw) + theme(plot.margin = unit(c(-0.3, 0.3, 0.25, -0.6), "cm")),
                   width = 10*1, height = 8.1*1, bg = 'white')
  
  # cowplot::ggsave2(filename = here('figures', 'main_paper', 'rcw_pop_overview_fig.png'),
  #                  plot = rcw_pop_overview_rcw,
  #                  width = 10*1, height = 8.1*1, bg = 'white')
}



#################################################################
### FIG 2: GENETIC CONTRIBUTIONS AND MORE INFO. ON INBREEDING ###
#################################################################

fig2_list <- list()

for (i in c('sex', 'source')) {
  
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
          plot.margin = unit(c(0.32,1, 0, 1), "cm"),
          legend.position = 'none',
          axis.title.y = element_text(margin = margin(0, 12.5, 0, 0))
    ) +
    scale_x_continuous(expand = c(0.01, 0), 
                       limits = c(1993.5, 2022.5),
                       position = "top") +
    xlab('Year') +
    ylab("Exp. genetic contribution")
  
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
      plot.margin = unit(c(0.25,1,-0.05,1), "cm"),
      legend.position = 'none') +
    ylab(substitute(paste("Scaled ", italic(F)["p"]))) +
    scale_x_continuous(expand = c(0.01, 0), limits = c(1993.5, 2022.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(clip = 'off')
  
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
  
  fig2_list[[i]][['panel1']] <- egg::ggarrange(fig2_list[[i]][['contr_plot']],
                                               fig2_list[[i]][['fped_prop_by_indiv_founder_plot']],
                                               fig2_list[[i]][['fped_prop_by_group_plot']],
                                               nrow = 3,
                                               labels = c('A', 'B', 'C'),
                                               label.args = list(gp = grid::gpar(font = 2, cex = 1.2)),
                                               heights = c(65, 65, 15))
  
  fig2_list[[i]][['panel2']] <- cowplot::ggdraw(fig2_list[[i]][['panel1']]) +
    cowplot::draw_image(perched_rcw,
                        scale = 0.195,
                        x = -0.3695, y = -.13) +
    theme(plot.margin = unit(c(0.2, -0.7, 0, 0.1), "cm"))
  
  if (isTRUE(output_figs)) {
    cowplot::ggsave2(filename = here('figures', 'main_paper', paste0('indiv_fped_gencontr_panel_', i, '.png') ),
                     plot = fig2_list[[i]][['panel2']],
                     width = 10*1, height = 6.5*1, bg = 'white')
  }
}








#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# ancestry_plot_lewontin <- results_list$rcw_inbr_merge %>%
#   mutate(count = 1) %>%
#   group_by(year) %>%
#   arrange(anc_prop) %>%
#   ggplot() +
#   geom_bar(aes(x = year, y = count, fill = anc_prop),
#            position = "stack", stat = "identity", width = 0.98) +
#   scale_fill_gradient(name = 'Expected\ntranslocated\nancestry', low = '#ededed', high = '#4c4c4c') +
#   geom_point(shape = "l", data = results_list$transloc_info_color,
#              aes(x = year, y = 0), size = 5.5, color = 'black',
#              alpha = 1,
#              position = position_jitter(width = 0.37, height = 0, seed = 5345),
#              show.legend = FALSE) +
#   scale_color_manual(values = setNames(results_list$transloc_info_color$source_color,# $sex_color,
#                                        nm = results_list$transloc_info_color$RCWid)) +
#   theme_bw() +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_line(linetype = 'dashed', color = '#d8d8d8'),
#     panel.border = element_blank(),
#     #axis.line.x = element_blank(),
#     #axis.ticks.x = element_blank(),
#     axis.line = element_line(color = '#4c4c4c', linewidth = 1.5)
#   ) +
#   #scale_y_reverse(expand = expansion(mult = c(0, 0), add = c(1, 0))) +
#   theme(
#     plot.margin=unit(c(-0.05,1, 0,1), "cm"),
#     axis.title = element_text(size = 18),
#     #axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 15)#,
#   ) +
#   ylab('Individuals') +
#   xlab('Year') + 
#   coord_cartesian(clip = 'off')
# 
# 
# inbreeding_plot_lewontin <- results_list$rcw_inbr_merge %>% 
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
#              shape = 21, colour = rcw_inbr_merge_processed$change_color, 
#              fill = 'white', 
#              size = 2, stroke = 1.4) +
#   theme_bw() +
#   theme(axis.title.y = element_text(size = 18),
#         axis.text.x = element_text(size = 15),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(linetype = 'dashed', 
#                                           color = '#d8d8d8'),
#         panel.border = element_blank(),
#         axis.line.x = element_blank(),
#         axis.line.y = element_line(color = '#4c4c4c', linewidth = 1.5),
#         axis.ticks = element_line(color = '#808080', linewidth = 0.7),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.margin = unit(c(0.1, 1, 0.2, 1), "cm"),
#         
#   ) +
#   scale_x_continuous(limits = c(1993.5, 2022.5), position = "top",
#                      expand = c(0.01, 0)) +
#   ylab('Inbreeding') +
#   coord_cartesian(clip = 'off')
# 
# 
# lewontin_fig <- cowplot::plot_grid(inbreeding_plot_lewontin,
#                    ancestry_plot_lewontin,
#                    nrow = 2,
#                    align = "v",
#                    axis = 'tblr',
#                    labels = c('A', 'B'),
#                    label_x = 0.09,
#                    label_y = c(0.89,1),
#                    label_size = 23) +
#   cowplot::draw_image(flying_rcw,
#                       scale = 0.25,
#                       x = 0.32, y = 0.35) +
#   theme(plot.margin = unit(c(0, -1, 0, -0.7), "cm"))
# 
# lewontin_fig <- ggdraw(
#   egg::ggarrange(inbreeding_plot_lewontin,
#                  ancestry_plot_lewontin,
#                  labels = c('A', 'B'),
#                  label.args = list(gp = grid::gpar(font = 1, cex = 1.4)),
#           #hjust = c(0.2, 0.9),
#           nrow = 2,
#           heights = c(40, 60))
#           ) +
#   cowplot::draw_image(flying_rcw,
#                       scale = 0.25,
#                       x = 0.32, y = 0.35) +
#   theme(plot.margin = unit(c(0, -1, 0, -0.7), "cm"))
# 
# 
# cowplot::ggsave2(filename = '/Users/alexlewanski/Documents/michigan_state/funding/lewontin_2024/lewontin_fig.png',
#                  plot = lewontin_fig,
#                  width = 9.6*1, height = 6.5*1, bg = 'white')
