##############################################################################################
### SCRIPT NAME: cmr_model.R
### PURPOSE: Conduct capture-mark-recapture modeling and output results
### PRODUCT:
###     translocation_model_selection_table.txt: table of model fit info for translocation mods 
###     anytransloc_model_selection_table: table of model fit info for translocation ancestry mods 
###     cjs_translocation_top_mod.csv:  parameter estimates info for top translocation mod
###     cjs_anytransloc_top_mod.csv: parameter estimates info for top translocation ancestry mod
###     markrecap_anytransloc_top_mod_plot.png:
### INFO: Translocation mods examine predictors of survival between translocated individuals and
###      individuals without any translocation ancestry. Translocation ancestry mods examine 
###      predictors of survival between individuals with any expected translocation ancestry and 
###      individuals without translocation ancestry.
##############################################################################################


#####################
### SCRIPT SET-UP ###
#####################

### PACKAGES ###
library(here)
library(RMark)
library(dplyr)
library(ggplot2)
library(kableExtra)


### LOADING AND INITIAL PROCESSING OF DATA ###
rcw_dat <- read.csv(here('data', 'feb2024_databasemarch2023_processed', 'rcw_mcr.csv'),
                    colClasses=c("character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

rcw_dat1 <- rcw_dat %>%
  mutate(ch = as.character(ch),
         any_transloc = factor(if_else(transloc > 0, 1, 0)),
         translocated = factor(translocated)) %>%
  filter(sex_unknown != 1) %>%
  select(ch, male, any_transloc, translocated) %>%
  rename(sex = male) %>%
  mutate(sex = factor(sex))



###################################################################################
### ANALYSIS 1: SURVIVAL BETWEEN BIRDS WITH NO TRANSLOCATION ANCESTRY VS. BIRDS ###
###             WITH ANY TRANSLOCATION ANCESTRY                                 ###
###################################################################################

#Process encounter history dataframe for MARK analysis
processed_data <- process.data(rcw_dat1, 
                               model="CJS", 
                               begin.time=1, 
                               groups=c("sex", "any_transloc"))

#Create design dataframes for MARK model specification
rcw.ddl <- make.design.data(processed_data)


### DEFINITION COMPONENTS OF MODEL FORMULA ###
# Define range of models for Phi
Phidot <- list(formula=~1)
Phitime <- list(formula=~time)
Phisex <- list(formula=~sex)
Phitrans <- list(formula=~any_transloc)
Phisextime <- list(formula=~sex+time)
Phisextime_int <- list(formula=~sex*time)
Phisextrans <- list(formula=~sex+any_transloc)
Phisextrans_int <- list(formula=~sex*any_transloc)
Phitimetrans <- list(formula=~time+any_transloc)
Phitimetrans_int <- list(formula=~time*any_transloc)
Phialladd <- list(formula=~sex+any_transloc+time)
Phifull <- list(formula=~sex*any_transloc*time)

# Define range of models for p
pdot=list(formula=~1)
ptime=list(formula=~time)

### RUN MODELS ###
rcw.phidot.pdot          =mark(processed_data,rcw.ddl,
                               model.parameters=list(Phi=Phidot,p=pdot),delete=TRUE)
rcw.phidot.ptime          =mark(processed_data,rcw.ddl,
                                model.parameters=list(Phi=Phidot,p=ptime),delete=TRUE)

rcw.phitime.pdot          =mark(processed_data,rcw.ddl,
                                model.parameters=list(Phi=Phitime,p=pdot),delete=TRUE)
rcw.phitime.ptime          =mark(processed_data,rcw.ddl,
                                 model.parameters=list(Phi=Phitime,p=ptime),delete=TRUE)

rcw.phisex.pdot          =mark(processed_data,rcw.ddl,
                               model.parameters=list(Phi=Phisex,p=pdot),delete=TRUE)
rcw.phisex.ptime          =mark(processed_data,rcw.ddl,
                                model.parameters=list(Phi=Phisex,p=ptime),delete=TRUE)

rcw.phitrans.pdot          =mark(processed_data,rcw.ddl,
                                 model.parameters=list(Phi=Phitrans,p=pdot),delete=TRUE)
rcw.phitrans.ptime          =mark(processed_data,rcw.ddl,
                                  model.parameters=list(Phi=Phitrans,p=ptime),delete=TRUE)

rcw.phisextime.pdot          =mark(processed_data,rcw.ddl,
                                   model.parameters=list(Phi=Phisextime,p=pdot),delete=TRUE)
rcw.phisextime.ptime          =mark(processed_data,rcw.ddl,
                                    model.parameters=list(Phi=Phisextime,p=ptime),delete=TRUE)

rcw.phisextime_int.pdot          =mark(processed_data,rcw.ddl,
                                       model.parameters=list(Phi=Phisextime_int,p=pdot),delete=TRUE)
rcw.phisextime_int.ptime          =mark(processed_data,rcw.ddl,
                                        model.parameters=list(Phi=Phisextime_int,p=ptime),delete=TRUE)

rcw.phisextrans.pdot          =mark(processed_data,rcw.ddl,
                                    model.parameters=list(Phi=Phisextrans,p=pdot),delete=TRUE)
rcw.phisextrans.ptime          =mark(processed_data,rcw.ddl,
                                     model.parameters=list(Phi=Phisextrans,p=ptime),delete=TRUE)

rcw.phisextrans_int.pdot          =mark(processed_data,rcw.ddl,
                                        model.parameters=list(Phi=Phisextrans_int,p=pdot),delete=TRUE)
rcw.phisextrans_int.ptime          =mark(processed_data,rcw.ddl,
                                         model.parameters=list(Phi=Phisextrans_int,p=ptime),delete=TRUE)

rcw.phitimetrans.pdot          =mark(processed_data,rcw.ddl,
                                     model.parameters=list(Phi=Phitimetrans,p=pdot),delete=TRUE)
rcw.phitimetrans.ptime          =mark(processed_data,rcw.ddl,
                                      model.parameters=list(Phi=Phitimetrans,p=ptime),delete=TRUE)

rcw.phitimetrans_int.pdot          =mark(processed_data,rcw.ddl,
                                         model.parameters=list(Phi=Phitimetrans_int,p=pdot),delete=TRUE)
rcw.phitimetrans_int.ptime          =mark(processed_data,rcw.ddl,
                                          model.parameters=list(Phi=Phitimetrans_int,p=ptime),delete=TRUE)

rcw.phialladd.pdot          =mark(processed_data,rcw.ddl,
                                  model.parameters=list(Phi=Phialladd,p=pdot),delete=TRUE)
rcw.phialladd.ptime          =mark(processed_data,rcw.ddl,
                                   model.parameters=list(Phi=Phialladd,p=ptime),delete=TRUE)

rcw.phifull.pdot          =mark(processed_data,rcw.ddl,
                                model.parameters=list(Phi=Phifull,p=pdot),delete=TRUE)
rcw.phifull.ptime          =mark(processed_data,rcw.ddl,
                                 model.parameters=list(Phi=Phifull,p=ptime),delete=TRUE)


### PROCESS MODEL OUTPUT ###
#collect models into a list and construct table of model results
anytransloc_model_list <- collect.models(
  lx = c("rcw.phidot.pdot",
         "rcw.phidot.ptime",
         "rcw.phitime.pdot",
         "rcw.phitime.ptime",
         "rcw.phisex.pdot",
         "rcw.phisex.ptime",
         "rcw.phitrans.pdot",
         "rcw.phitrans.ptime",
         "rcw.phisextime.pdot",
         "rcw.phisextime.ptime",
         "rcw.phisextime_int.pdot",
         "rcw.phisextime_int.ptime",
         "rcw.phisextrans.pdot",
         "rcw.phisextrans.ptime",
         "rcw.phisextrans_int.pdot",
         "rcw.phisextrans_int.ptime",
         "rcw.phitimetrans.pdot",
         "rcw.phitimetrans.ptime",
         "rcw.phitimetrans_int.pdot",
         "rcw.phitimetrans_int.ptime",
         "rcw.phialladd.pdot",
         "rcw.phialladd.ptime",
         "rcw.phifull.pdot",
         "rcw.phifull.ptime"),
  type = NULL,
  table = TRUE,
  adjust = TRUE,
  external = FALSE
)

anytransloc_model_selection <- model.table(anytransloc_model_list)
#write.table(model_selection,"anyTransAICc.txt", sep="\t")

#output top model
#rcw.phisextrans.pdot          <- mark(processed_data,rcw.ddl,
#                                    model.parameters=list(Phi=Phisextrans,p=pdot),
#                                    delete=FALSE)

anytransloc_top_mod <- rcw.phisextrans.pdot$results$real



#########################################################################
### ANALYSIS 2: SURVIVAL BETWEEN TRANSLOCATED BIRDS AND BIRDS WITH NO ### 
###             TRANSLOCATION ANCESTRY                                ###
#########################################################################

processed_data <- process.data(rcw_dat1 %>%
                                 filter(translocated == '1' | any_transloc == '0'), 
                               model="CJS", 
                               begin.time=1, 
                               groups=c("sex", "translocated"))

#Create design dataframes for MARK model specification
rcw.ddl <- make.design.data(processed_data)


# Define range of models for Phi
Phidot <- list(formula=~1)
Phitime <- list(formula=~time)
Phisex <- list(formula=~sex)
Phitrans <- list(formula=~translocated)
Phisextime <- list(formula=~sex+time)
Phisextime_int <- list(formula=~sex*time)
Phisextrans <- list(formula=~sex+translocated)
Phisextrans_int <- list(formula=~sex*translocated)
Phitimetrans <- list(formula=~time+translocated)
Phitimetrans_int <- list(formula=~time*translocated)
Phialladd <- list(formula=~sex+translocated+time)
Phifull <- list(formula=~sex*translocated*time)


# Define range of models for p
pdot=list(formula=~1)
ptime=list(formula=~time)


### RUN MODELS ###
rcw.phidot.pdot          =mark(processed_data,rcw.ddl,
                               model.parameters=list(Phi=Phidot,p=pdot),delete=TRUE)
rcw.phidot.ptime          =mark(processed_data,rcw.ddl,
                                model.parameters=list(Phi=Phidot,p=ptime),delete=TRUE)

rcw.phitime.pdot          =mark(processed_data,rcw.ddl,
                                model.parameters=list(Phi=Phitime,p=pdot),delete=TRUE)
rcw.phitime.ptime          =mark(processed_data,rcw.ddl,
                                 model.parameters=list(Phi=Phitime,p=ptime),delete=TRUE)

rcw.phisex.pdot          =mark(processed_data,rcw.ddl,
                               model.parameters=list(Phi=Phisex,p=pdot),delete=TRUE)
rcw.phisex.ptime          =mark(processed_data,rcw.ddl,
                                model.parameters=list(Phi=Phisex,p=ptime),delete=TRUE)

rcw.phitrans.pdot          =mark(processed_data,rcw.ddl,
                                 model.parameters=list(Phi=Phitrans,p=pdot),delete=TRUE)
rcw.phitrans.ptime          =mark(processed_data,rcw.ddl,
                                  model.parameters=list(Phi=Phitrans,p=ptime),delete=TRUE)

rcw.phisextime.pdot          =mark(processed_data,rcw.ddl,
                                   model.parameters=list(Phi=Phisextime,p=pdot),delete=TRUE)
rcw.phisextime.ptime          =mark(processed_data,rcw.ddl,
                                    model.parameters=list(Phi=Phisextime,p=ptime),delete=TRUE)

rcw.phisextime_int.pdot          =mark(processed_data,rcw.ddl,
                                       model.parameters=list(Phi=Phisextime_int,p=pdot),delete=TRUE)
rcw.phisextime_int.ptime          =mark(processed_data,rcw.ddl,
                                        model.parameters=list(Phi=Phisextime_int,p=ptime),delete=TRUE)

rcw.phisextrans.pdot          =mark(processed_data,rcw.ddl,
                                    model.parameters=list(Phi=Phisextrans,p=pdot),delete=TRUE)
rcw.phisextrans.ptime          =mark(processed_data,rcw.ddl,
                                     model.parameters=list(Phi=Phisextrans,p=ptime),delete=TRUE)

rcw.phisextrans_int.pdot          =mark(processed_data,rcw.ddl,
                                        model.parameters=list(Phi=Phisextrans_int,p=pdot),delete=TRUE)
rcw.phisextrans_int.ptime          =mark(processed_data,rcw.ddl,
                                         model.parameters=list(Phi=Phisextrans_int,p=ptime),delete=TRUE)

rcw.phitimetrans.pdot          =mark(processed_data,rcw.ddl,
                                     model.parameters=list(Phi=Phitimetrans,p=pdot),delete=TRUE)
rcw.phitimetrans.ptime          =mark(processed_data,rcw.ddl,
                                      model.parameters=list(Phi=Phitimetrans,p=ptime),delete=TRUE)

rcw.phitimetrans_int.pdot          =mark(processed_data,rcw.ddl,
                                         model.parameters=list(Phi=Phitimetrans_int,p=pdot),delete=TRUE)
rcw.phitimetrans_int.ptime          =mark(processed_data,rcw.ddl,
                                          model.parameters=list(Phi=Phitimetrans_int,p=ptime),delete=TRUE)

rcw.phialladd.pdot          =mark(processed_data,rcw.ddl,
                                  model.parameters=list(Phi=Phialladd,p=pdot),delete=TRUE)
rcw.phialladd.ptime          =mark(processed_data,rcw.ddl,
                                   model.parameters=list(Phi=Phialladd,p=ptime),delete=TRUE)

rcw.phifull.pdot          =mark(processed_data,rcw.ddl,
                                model.parameters=list(Phi=Phifull,p=pdot),delete=TRUE)
rcw.phifull.ptime          =mark(processed_data,rcw.ddl,
                                 model.parameters=list(Phi=Phifull,p=ptime),delete=TRUE)


#collect models into a list and construct table of model results
translocation_model_list <- collect.models(
  lx = c("rcw.phidot.pdot",
         "rcw.phidot.ptime",
         "rcw.phitime.pdot",
         "rcw.phitime.ptime",
         "rcw.phisex.pdot",
         "rcw.phisex.ptime",
         "rcw.phitrans.pdot",
         "rcw.phitrans.ptime",
         "rcw.phisextime.pdot",
         "rcw.phisextime.ptime",
         "rcw.phisextime_int.pdot",
         "rcw.phisextime_int.ptime",
         "rcw.phisextrans.pdot",
         "rcw.phisextrans.ptime",
         "rcw.phisextrans_int.pdot",
         "rcw.phisextrans_int.ptime",
         "rcw.phitimetrans.pdot",
         "rcw.phitimetrans.ptime",
         "rcw.phitimetrans_int.pdot",
         "rcw.phitimetrans_int.ptime",
         "rcw.phialladd.pdot",
         "rcw.phialladd.ptime",
         "rcw.phifull.pdot",
         "rcw.phifull.ptime"),
  type = NULL,
  table = TRUE,
  adjust = TRUE,
  external = FALSE
)

translocation_top_mod <- rcw.phisextrans.pdot$results$real



##############################################
### PROCESSING AND OUTPUT OF MODEL RESULTS ###
##############################################

### ANALYSIS 1 RESULTS PROCESSING ###

anytransloc_top_mod$param_type <- gsub('([Pp]*)\\s.*', '\\1', rownames(anytransloc_top_mod))
anytransloc_top_mod$sex <- substr(gsub('.*g([01][01]).*', '\\1', rownames(anytransloc_top_mod)), 1, 1)
anytransloc_top_mod$transloc <- substr(gsub('.*g([01][01]).*', '\\1', rownames(anytransloc_top_mod)), 2, 2)

write.csv(anytransloc_top_mod,
          here('results', 'cjs_anytransloc_top_mod.csv'),
          row.names = FALSE)

anytransloc_model_selection_table <- anytransloc_model_selection %>%
  mutate(AICc = round(AICc, 2),
         DeltaAICc = round(DeltaAICc, 2),
         weight = signif(weight, digits = 3),
         Deviance = round(Deviance, 2)) %>%
  select(model, AICc, DeltaAICc, weight, npar, Deviance) %>% 
  kbl('latex', booktabs = TRUE, 
      align = "c",
      linesep = "",
      row.names = FALSE) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position"),
                position = "center",
                font_size = 7)

cat(anytransloc_model_selection_table, 
    file = here('tables', 'anytransloc_model_selection_table.txt'), 
    append = FALSE)


anytransloc_top_mod_param_table <- anytransloc_top_mod %>%
  mutate(group = paste(if_else(sex == '0', 'female', 'male'),
                       if_else(transloc == '0', '0% transloc. ancestry', '>0% transloc. ancestry')),
         estimate = signif(estimate, 2),
         se = signif(se, 2),
         lcl = signif(lcl, 2),
         ucl = signif(ucl, 2)) %>%
  select(param_type, group, estimate, se, lcl, ucl) %>%
  rename(Parameter = param_type,
         Group = group,
         Estimate = estimate,
         `Standard error` = se,
         `lower 95% CI` = lcl,
         `upper 95% CI` = ucl) %>% 
  mutate(Group = if_else(Parameter == 'p', 'all', Group)) %>%
  kbl('latex', booktabs = TRUE, 
      align = "c",
      linesep = "",
      row.names = FALSE) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position"),
                position = "center",
                font_size = 7)

cat(anytransloc_top_mod_param_table, 
    file = here('tables', 'anytransloc_top_mod_param_table.txt'), 
    append = FALSE)



### ANALYSIS 2 RESULTS PROCESSING ###
translocation_model_selection <- model.table(translocation_model_list)

translocation_model_selection_table <- translocation_model_selection %>%
  mutate(AICc = round(AICc, 2),
         DeltaAICc = round(DeltaAICc, 2),
         weight = signif(weight, digits = 3),
         Deviance = round(Deviance, 2)) %>%
  select(model, AICc, DeltaAICc, weight, npar, Deviance) %>% 
  kbl('latex', booktabs = TRUE, 
      align = "c",
      linesep = "",
      row.names = FALSE) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position"),
                position = "center",
                font_size = 7)

cat(translocation_model_selection_table, 
    file = here('tables', 'translocation_model_selection_table.txt'), 
    append = FALSE)

#output top model
#rcw.phisextrans.pdot          <- mark(processed_data,rcw.ddl,
#                                    model.parameters=list(Phi=Phisextrans,p=pdot),
#                                    delete=FALSE)


translocation_top_mod$param_type <- gsub('([Pp]*)\\s.*', '\\1', rownames(translocation_top_mod))
translocation_top_mod$sex <- substr(gsub('.*g([01][01]).*', '\\1', rownames(translocation_top_mod)), 1, 1)
translocation_top_mod$transloc <- substr(gsub('.*g([01][01]).*', '\\1', rownames(translocation_top_mod)), 2, 2)

write.csv(translocation_top_mod,
          here('results', 'cjs_translocation_top_mod.csv'),
          row.names = FALSE)

translocation_model_selection_table <- translocation_model_selection %>%
  mutate(AICc = round(AICc, 2),
         DeltaAICc = round(DeltaAICc, 2),
         weight = signif(weight, digits = 3),
         Deviance = round(Deviance, 2)) %>%
  select(model, AICc, DeltaAICc, weight, npar, Deviance) %>% 
  kbl('latex', booktabs = TRUE, 
      align = "c",
      linesep = "",
      row.names = FALSE) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position"),
                position = "center",
                font_size = 7)

cat(translocation_model_selection_table, 
    file = here('tables', 'translocation_model_selection_table.txt'), 
    append = FALSE)


translocation_top_mod_param_table <- translocation_top_mod %>%
  mutate(group = paste(if_else(sex == '0', 'female', 'male'),
                       if_else(transloc == '0', 'non-translocated (0% transloc. ancestry)', 'translocated')),
         estimate = signif(estimate, 2),
         se = signif(se, 2),
         lcl = signif(lcl, 2),
         ucl = signif(ucl, 2)) %>%
  select(param_type, group, estimate, se, lcl, ucl) %>%
  rename(Parameter = param_type,
         Group = group,
         Estimate = estimate,
         `Standard error` = se,
         `lower 95% CI` = lcl,
         `upper 95% CI` = ucl) %>% 
  mutate(Group = if_else(Parameter == 'p', 'all', Group)) %>%
  kbl('latex', booktabs = TRUE, 
      align = "c",
      linesep = "",
      row.names = FALSE) %>% 
  kable_styling(latex_options = c("scale_down", "hold_position"),
                position = "center",
                font_size = 7)

cat(translocation_top_mod_param_table, 
    file = here('tables', 'translocation_top_mod_param_table.txt'), 
    append = FALSE)


### Plot ##
#translocation model results are visualized in the transloc. viz figure in the main text
anytransloc_top_mod_plot <- anytransloc_top_mod %>%
  mutate(sex_update = if_else(sex == '0', 'female', 'male'),
         transloc_update = factor(if_else(transloc == '0', '0% translocation\nancestry', '>%0 Translocation\nancestry'),
                                  levels = c('>%0 Translocation\nancestry', '0% translocation\nancestry'))) %>%
  filter(param_type == 'Phi') %>%
  ggplot(aes(x = sex_update, y = estimate)) +
  # geom_segment(aes(x = level2, y = low95, yend = high95, color = level1), 
  #              linewidth = 6, lineend = 'round') +
  # geom_point(aes(shape = level2), 
  #            size = 3.5, color = 'white') +
  geom_segment(aes(x = sex_update, y = lcl, yend = ucl, color = transloc_update), 
               linewidth = 5.25, lineend = 'round', alpha = 0.65) +
  geom_point(aes(shape = sex_update, color = transloc_update, fill = transloc_update), 
             size = 9, stroke = 0) + 
  scale_shape_manual(name = 'Sex',
                     values = c(21, 22)) + #female: circle; male: square
  facet_wrap(~transloc_update, strip.position = "top", nrow = 1) +
  scale_color_manual(values = c('#818181', "#d8d8d8")) +
  scale_fill_manual(values = c('#4c4c4c', "#acacac")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(size = 12, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = 'transparent', color = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = 'transparent'), 
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dashed', color = '#d8d8d8'),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = 'black', linewidth = 1),
        #axis.ticks.y = element_blank(),
        legend.position = 'none',
        axis.text.y = element_text(size = 13),
        strip.background = element_blank()
  ) +
  ylab('Apparent survival') +
  coord_cartesian(clip = 'off')

ggsave(filename = here('figures', 'supplement', 'figures', 'markrecap_anytransloc_top_mod_plot.png' ),
                 plot = anytransloc_top_mod_plot,
                 width = 5, height = 4, bg = 'white')



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

#translocation_top_mod %>%
#  mutate(sex_update = if_else(sex == '0', 'female', 'male'),
#         transloc_update = if_else(transloc == '0', 'Non-translocated\n(0% transloc. ancestry)', 'Translocated')) %>%
#  filter(param_type == 'Phi') %>%
#  ggplot(aes(x = sex_update, y = estimate)) +
#  # geom_segment(aes(x = level2, y = low95, yend = high95, color = level1), 
#  #              linewidth = 6, lineend = 'round') +
#  # geom_point(aes(shape = level2), 
#  #            size = 3.5, color = 'white') +
#  geom_segment(aes(x = sex_update, y = lcl, yend = ucl, color = transloc_update), 
#               linewidth = 5.25, lineend = 'round', alpha = 0.65) +
#  geom_point(aes(shape = sex_update, color = transloc_update, fill = transloc_update), 
#             size = 9, stroke = 0) + 
#  scale_shape_manual(name = 'Sex',
#                     values = c(21, 22)) + #female: circle; male: square
#  facet_wrap(~transloc_update, strip.position = "top", nrow = 1) +
#  scale_color_manual(values = c("#d8d8d8", '#4c4c4c')) +
#  scale_fill_manual(values = c("#d8d8d8", '#4c4c4c')) +
#  theme_bw() +
#  theme(axis.title.x = element_blank(), 
#        strip.placement = "outside",
#        strip.text = element_text(size = 12, vjust = 1),
#        axis.ticks.x = element_blank(),
#        panel.background = element_rect(fill = 'transparent', color = 'transparent'),
#        plot.background = element_rect(fill = 'transparent', color = 'transparent'), 
#        axis.text.x = element_text(size = 10),
#        axis.title.y = element_text(size = 12),
#        panel.border = element_blank(),
#        panel.grid.major.x = element_blank(),
#        panel.grid.major.y = element_line(linetype = 'dashed', color = '#d8d8d8'),
#        panel.grid.minor = element_blank(),
#        axis.line.y = element_line(color = 'black', linewidth = 1),
#        #axis.ticks.y = element_blank(),
#        legend.position = 'none',
#        axis.text.y = element_text(size = 13),
#        strip.background = element_blank()
#  ) +
#  ylab('Apparent survival') +
#  coord_cartesian(clip = 'off')


# #plot real parameter estimates and 95% CI from top supported model
# # Create the data frame with the provided data
# phi_results <- data.frame(
#   group = c('female_residents', 'male_residents', 'female_translocated', 'male_translocated'),
#   estimate = c(0.706, 0.768, 0.765, 0.818),
#   se = c(0.016, 0.013, 0.015, 0.012),
#   low95 = c(0.675, 0.742, 0.735, 0.794),
#   high95 = c(0.735, 0.792, 0.793, 0.840)
# )
# 
# 
# ggplot(phi_results, aes(x=group, y=estimate)) +
#   geom_point(size=3) +
#   geom_segment(aes(ymin=low95, ymax=high95), width=0.2) +
#   #labs(title="survival estimates with 95% Confidence Intervals", x="", y="phi") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=45, hjust=1))
# 
# 
# # Plot the estimates with 95% confidence intervals
# ggplot(phi_results, aes(x=group, y=estimate)) +
#   geom_point(size=3) +
#   geom_errorbar(aes(ymin=low95, ymax=high95), width=0.2) +
#   labs(title="survival estimates with 95% Confidence Intervals", x="", y="phi") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=45, hjust=1))
# 
# ###################################
# 
# #ORIGIN dataset
# setwd("~/Desktop/rcw_mark/origin")
# 
# 
# data <- read.table("origin.txt", header=TRUE, stringsAsFactors=FALSE,
#                    colClasses=c("character", "factor", "factor"))
# 
# 
# processed_data <- process.data(data, model="CJS", begin.time=1, groups=c("sex", "origin"))
# rcw.ddl=make.design.data(processed_data)
# 
# # Define and run the models
# 
# #  Define range of models for Phi
# #
# Phidot=list(formula=~1)
# Phitime=list(formula=~time)
# Phisex=list(formula=~sex)
# Phitrans=list(formula=~origin)
# Phisextime=list(formula=~sex+time)
# Phisextime_int=list(formula=~sex*time)
# Phisextrans=list(formula=~sex+origin)
# Phisextrans_int=list(formula=~sex*origin)
# Phitimetrans=list(formula=~time+origin)
# Phitimetrans_int=list(formula=~time*origin)
# Phialladd=list(formula=~sex+origin+time)
# Phifull=list(formula=~sex*origin*time)
# 
# 
# #
# #  Define range of models for p
# #
# pdot=list(formula=~1)
# ptime=list(formula=~time)
# 
# # Run assortment of models
# #
# rcw.phidot.pdot          =mark(processed_data,rcw.ddl,
#                                model.parameters=list(Phi=Phidot,p=pdot),delete=TRUE)
# rcw.phidot.ptime          =mark(processed_data,rcw.ddl,
#                                 model.parameters=list(Phi=Phidot,p=ptime),delete=TRUE)
# 
# rcw.phitime.pdot          =mark(processed_data,rcw.ddl,
#                                 model.parameters=list(Phi=Phitime,p=pdot),delete=TRUE)
# rcw.phitime.ptime          =mark(processed_data,rcw.ddl,
#                                  model.parameters=list(Phi=Phitime,p=ptime),delete=TRUE)
# 
# rcw.phisex.pdot          =mark(processed_data,rcw.ddl,
#                                model.parameters=list(Phi=Phisex,p=pdot),delete=TRUE)
# rcw.phisex.ptime          =mark(processed_data,rcw.ddl,
#                                 model.parameters=list(Phi=Phisex,p=ptime),delete=TRUE)
# 
# rcw.phitrans.pdot          =mark(processed_data,rcw.ddl,
#                                  model.parameters=list(Phi=Phitrans,p=pdot),delete=TRUE)
# rcw.phitrans.ptime          =mark(processed_data,rcw.ddl,
#                                   model.parameters=list(Phi=Phitrans,p=ptime),delete=TRUE)
# 
# rcw.phisextime.pdot          =mark(processed_data,rcw.ddl,
#                                    model.parameters=list(Phi=Phisextime,p=pdot),delete=TRUE)
# rcw.phisextime.ptime          =mark(processed_data,rcw.ddl,
#                                     model.parameters=list(Phi=Phisextime,p=ptime),delete=TRUE)
# 
# rcw.phisextime_int.pdot          =mark(processed_data,rcw.ddl,
#                                        model.parameters=list(Phi=Phisextime_int,p=pdot),delete=TRUE)
# rcw.phisextime_int.ptime          =mark(processed_data,rcw.ddl,
#                                         model.parameters=list(Phi=Phisextime_int,p=ptime),delete=TRUE)
# 
# rcw.phisextrans.pdot          =mark(processed_data,rcw.ddl,
#                                     model.parameters=list(Phi=Phisextrans,p=pdot),delete=TRUE)
# rcw.phisextrans.ptime          =mark(processed_data,rcw.ddl,
#                                      model.parameters=list(Phi=Phisextrans,p=ptime),delete=TRUE)
# 
# rcw.phisextrans_int.pdot          =mark(processed_data,rcw.ddl,
#                                         model.parameters=list(Phi=Phisextrans_int,p=pdot),delete=TRUE)
# rcw.phisextrans_int.ptime          =mark(processed_data,rcw.ddl,
#                                          model.parameters=list(Phi=Phisextrans_int,p=ptime),delete=TRUE)
# 
# rcw.phitimetrans.pdot          =mark(processed_data,rcw.ddl,
#                                      model.parameters=list(Phi=Phitimetrans,p=pdot),delete=TRUE)
# rcw.phitimetrans.ptime          =mark(processed_data,rcw.ddl,
#                                       model.parameters=list(Phi=Phitimetrans,p=ptime),delete=TRUE)
# 
# rcw.phitimetrans_int.pdot          =mark(processed_data,rcw.ddl,
#                                          model.parameters=list(Phi=Phitimetrans_int,p=pdot),delete=TRUE)
# rcw.phitimetrans_int.ptime          =mark(processed_data,rcw.ddl,
#                                           model.parameters=list(Phi=Phitimetrans_int,p=ptime),delete=TRUE)
# 
# rcw.phialladd.pdot          =mark(processed_data,rcw.ddl,
#                                   model.parameters=list(Phi=Phialladd,p=pdot),delete=TRUE)
# rcw.phialladd.ptime          =mark(processed_data,rcw.ddl,
#                                    model.parameters=list(Phi=Phialladd,p=ptime),delete=TRUE)
# 
# rcw.phifull.pdot          =mark(processed_data,rcw.ddl,
#                                 model.parameters=list(Phi=Phifull,p=pdot),delete=TRUE)
# rcw.phifull.ptime          =mark(processed_data,rcw.ddl,
#                                  model.parameters=list(Phi=Phifull,p=ptime),delete=TRUE)
# 
# 
# #collect models into a list and construct table of model results
# model_list<-collect.models(
#   lx = NULL,
#   type = NULL,
#   table = TRUE,
#   adjust = TRUE,
#   external = FALSE
# )
# 
# model_selection <- model.table(model_list)
# write.table(model_selection,"originAICc.txt", sep="\t")
# 
# #output top model
# rcw.phisextrans.pdot          =mark(processed_data,rcw.ddl,
#                                     model.parameters=list(Phi=Phisextrans,p=pdot),
#                                     delete=FALSE,
#                                     prefix = 'test')
# 
# rcw.phisextrans.pdot$results$real
# 
# #plot real parameter estimates and 95% CI from top supported model
# # Create the data frame with the provided data
# phi_results <- data.frame(
#   group = c('female_resident', 'male_resident', 'female_translocated', 'male_translocated'),
#   estimate = c(0.709, 0.765, 0.799, 0.842),
#   se = c(0.017, 0.014, 0.032, 0.026),
#   low95 = c(0.674, 0.736, 0.730, 0.784),
#   high95 = c(0.741, 0.792, 0.855, 0.887)
# )
# 
# # Plot the estimates with 95% confidence intervals
# ggplot(phi_results, aes(x=group, y=estimate)) +
#   geom_point(size=3) +
#   geom_errorbar(aes(ymin=low95, ymax=high95), width=0.2) +
#   labs(title="survival estimates with 95% Confidence Intervals", x="", y="phi") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=45, hjust=1))
# 
# 
# # table(rcw_dat$sex_unknown)
# # 
# # rcw_dat1 %>%
# #   group_by(any_transloc) %>%
# #   summarize(count = n())
# # 
# # rcw_dat1 %>%
# #   filter(translocated == '1' | any_transloc == '0') %>%
# #   group_by(translocated) %>%
# #   summarize(count = n())
