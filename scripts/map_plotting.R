##############################################################################################
### SCRIPT NAME: map_plotting.R
### PURPOSE: create map of RCW distribution and study populations
### PRODUCT:
###     rcw_range_pop_map.png: study system map (included in the supplementary material)
##############################################################################################


#######################
### BACKGROUND INFO ###
#######################

### NAMES ###
#APAFR: Avon Park
#ONF: Osceola National Forest
#FTS: Fort Stewart
#WSF-CITRUS: Withlacoochee State Forest - Citrus
#FTB: Fort Benning
#ANF: Apalachicola National Forest
#CBJTC: Camp Blanding Joint Training Center

### COORDINATES FOR MAPPING ###
#APAFR: 27.644521, -81.343642
#ONF: 30.270959, -82.495517
#FTS: 31.875167, -81.615220
#WSF-CITRUS: 28.718578, -82.490166
#FTB: 32.359721, -84.949391
#ANF: 30.083168, -84.746376
#CBJTC: 29.985615, -81.993976



##########################
### SCRIPT PREPARATION ###
##########################

### libraries ###
#library(fields)
#remotes::install_github("ebird/ebirdst")
library(here)
library(tidyverse)
library(ebirdst)
library(rnaturalearth)
library(sf)
library(terra)
library(raster)


### DOWNLOAD DATA ###
set_ebirdst_access_key('olf2omjlol90')
ebirdst_download_status(species = "Red-cockaded Woodpecker")



######################################
### RANGE-WIDE MAP (IN SUPPLEMENT) ###
######################################

### LOADING AND PROCESSING DATA ###
crs_laea <- paste0("+proj=laea +lat_0=", 31,
                   " +lon_0=", -86)

rcw_resident <- load_raster("Red-cockaded Woodpecker", 
                            product = "abundance", 
                            period = "seasonal",
                            metric = "mean",
                            resolution = "3km")[['resident']]

us_mexico_map <- rnaturalearth::ne_states(country = c("United States of America", "Mexico"),
                       returnclass = "sf") %>% 
  st_transform(crs = crs_laea) %>% 
  st_geometry()

rcw_resident_laea <- terra::project(rcw_resident, crs_laea, method = 'near') |> 
  # remove areas of the raster containing no data
  terra::trim()

pop_order <- c('APAFR', 'ONF', 'FTS', 'WSF-CITRUS', 'FTB', 'ANF', 'CBJTC')
pop_df <- data.frame(pop = c("ANF", "CBJTC", "FTB", "FTS", "ONF", "WSF-CITRUS", "APAFR"),
                     lat = c(30.083168, 29.985615, 32.359721, 31.875167, 30.270959, 28.718578, 27.644521), 
                     lon = c(-84.746376, -81.993976, -84.949391, -81.615220, -82.495517, -82.490166, -81.343642),
                     size = c(21, 1, 4, 16, 6, 6, 126),
                     color = c("#f0597d", "#ffdb19", "#595959", "#06d6a0", "#118ab2", "#800080", 'black'))

rcw_pop_vec <- vect(cbind(pop_df$lon[pop_df$pop != 'APAFR'],
                         pop_df$lat[pop_df$pop != 'APAFR']), 
                   crs="+proj=longlat")
rcw_pop_vec_proj <- terra::project(rcw_pop_vec, 
                                 crs_laea)

apafr_vec <- vect(cbind(pop_df[pop_df$pop == 'APAFR',,drop=FALSE]$lon,
                          pop_df[pop_df$pop == 'APAFR',,drop=FALSE]$lat), 
                    crs="+proj=longlat")
apafr_vec_proj <- terra::project(apafr_vec, 
                                   crs_laea)

rcw_load_map_params <- load_fac_map_parameters("Red-cockaded Woodpecker", 
                                               path = ebirdst_data_dir())
color_generator <- colorRampPalette(c("#eac7b7","#814d36"))


### CREATING MAP ###
png(here('figures', 'supplement', 'figures', 'rcw_range_pop_map.png'), 
    width = 8, height = 6.4, units = "in", res = 1000)
par(mar = c(0.5, 0.5, 0.5, 2)) #bottom, left, top, right
plot(us_mexico_map, xlim  = c(-1.2e06, 1e06), ylim = c(-8.5e05, 8e05),
     col = '#f9f9f9', bg = 'white', grid = TRUE, border = "#d8d8d8", lwd = 1.9)
plot(rcw_resident_laea, grid = FALSE,
     #breaks = c(0, pars$seasonal_bins),
     #col = c(NA, rep(c("#dca388", "#90563c", "#673e2b"), each = 7)),
     breaks = c(0, rcw_load_map_params$seasonal_bins),
     #col = c(NA, rep(c("#dca388", "#90563c", "#673e2b"), each = 7)),
     col = c(NA, color_generator(length(rcw_load_map_params$seasonal_bins))),
     xlim  = c(-1.2e06, 1e06), ylim = c(-9e05, 8e05),
     legend = FALSE, axes = FALSE, maxpixels = ncell(rcw_resident_laea),
     add = TRUE)
points(rcw_pop_vec_proj, bg = pop_df$color[pop_df$pop != 'APAFR'], col = "#595959", 
       pch = 21,
       cex = 0.25 + sqrt(pop_df$size[pop_df$pop != 'APAFR']/pi))
points(apafr_vec_proj, bg = c('black'), col = 'black', 
       pch = 18,
       cex = 2)
sbar(500000, xy = c(6e5, -8.8e5), 
     cex = 0.8, below = "km",
     label = c(0, 250, 500), ticks = FALSE, 
     lwd = 2)
legend(x = 7.5e5, 
       y = 2e5,
        legend = pop_order,
        pt.bg = c(NA, pop_df$color[match(pop_order, pop_df$pop)][-1]),
       pt.cex = c(1.8, rep(1.25, 6)),
       col = 'black', cex = 0.8,
       pch = c(18, rep(21, 6)), 
       bty = "n")

text(x = 7.75e5, 
     y = -4.3e5,
     adj = c(0, 0),
     cex = 0.85,
     'Relative\nabundance')

legend(x = 7.45e5, 
       y = -4.5e5,
       legend = c('high', rep(NA, length(rcw_load_map_params$seasonal_bins) - 2 ), 'low'),
       fill = rev(color_generator(length(rcw_load_map_params$seasonal_bins))),
       border = NA,
       y.intersp = 0.14,
       cex = 0.8,
       text.font = 1,
       bty = "n",
       #title = 'Abundance',
       horiz = FALSE)
dev.off()


### RESOURCES ###
#https://stackoverflow.com/questions/75928646/cannot-use-project-function-to-reproject-spatvector-to-epsg3035
#https://ebird.github.io/ebird-best-practices/ebird.html
#https://ebird.github.io/ebirdst/articles/trends.html
#https://strimas.com/bioee-4751_ebirdst/
#https://stackoverflow.com/questions/2579995/control-the-size-of-points-in-an-r-scatterplot



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# states <- ne_states(iso_a2 = c("US", "CA"), returnclass = "sf")
# 
# map(rcw_resident)
# 
# plot(rcw_resident, axes = FALSE)
# 
# library(tidyterra)
# 
# 
# ggplot() +
#   geom_spatraster(data = rcw_resident) +
#   # You can use coord_sf
#   coord_sf(crs = 3857) +
#   scale_fill_hypso_c()
# library(maps)
# library(mapdata)
# 
# rcw_range_2022 <- read_sf('/Users/alexlewanski/Documents/michigan_state/research/rcwo_project1/figures/supplement/figures/recwoo_range_2022/recwoo_range_2022.gpkg')
# 
# class(rcw_range_2022)
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# 
# rcw_range_2022
# 
# 
# ggplot() +
#   geom_sf(data = world) +
#   geom_sf(data = st_transform(rcw_range_2022),
#           fill = '#4c4c4c') +
#   coord_sf(xlim = c(-100, -60),
#                   ylim = c(20, 40))
# 
# ggplot() + 
#   geom_polygon(data=world %>% 
#                  filter(region %in% c('USA', 'Canada', 'Mexico')) %>% 
#                  filter(long > -200 & long < -50), 
#                aes(x=long, y=lat, group=group),
#                fill='gray', color = 'white', linewidth = 0.5) +
#   #coord_fixed(1.3) +
#   # geom_sf(data = rcw_range_2022,
#   #         fill = '#4c4c4c') +
#   # theme(panel.background = element_blank(),
#   #       panel.grid.major = element_line(linetype = 'dashed', 
#   #                                       color = 'gray')) +
#   coord_cartesian(xlim = c(-100, -60),
#                   ylim = c(20, 40)) +
#   geom_sf(data = st_transform(rcw_range_2022),
#           fill = '#4c4c4c')
#   
#   #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
#   #      axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
# 
# class(rcw_range_2022)
#   
#   geom_sf(data = rcw_range_2022,
#           fill = '#4c4c4c') +
#   theme(panel.background = element_blank()) 
#           
# plot(st_geometry(rcw_range_2022))
