###########################
###########################
### STUDY SYSTEM MAP(S) ###
###########################
###########################

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
                            resolution = "9km")[['resident']]

us_mexico_map <- rnaturalearth::ne_states(country = c("United States of America", "Mexico"),
                       returnclass = "sf") %>% 
  st_transform(crs = crs_laea) %>% 
  st_geometry()

rcw_resident_laea <- terra::project(rcw_resident, crs_laea, method = 'near') |> 
  # remove areas of the raster containing no data
  terra::trim()

pop_order <- c('APAFR', 'ONF', 'FTS', 'WSF-CITRUS', 'FTB', 'ANF', 'CBJTC')
pop_df <- data.frame(pop = c("ANF", "CBJTC", "FTB", "FTS", "ONF", "WSF-CITRUS", "APAFR"),
                     lat = c(30.207865, 29.984685, 31.447247, 33.504599, 30.200350, 28.716566, 27.644521), 
                     lon = c(-84.663056, -81.989727, -84.345396, -83.389586, -82.442920, -82.493976, -81.343642),
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


### CREATING MAP ###
png(here('figures', 'supplement', 'figures', 'rcw_range_pop_map.png'), 
    width = 8, height = 6.4, units = "in", res = 500)
par(mar = c(0.5, 0.5, 0.5, 2)) #bottom, left, top, right
plot(us_mexico_map, xlim  = c(-1.2e06, 1e06), ylim = c(-8.5e05, 8e05),
     col = '#f2f2f2', bg = 'white', grid = TRUE, border = "#d8d8d8", lwd = 1.9)
plot(rcw_resident_laea, grid = FALSE,
     breaks = c(0, pars$seasonal_bins),
     col = c(NA, rep(c("#dca388", "#90563c", "#673e2b"), each = 7)),
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
legend(x = 7e5, 
       y = 2e5,
        legend = pop_order,
        pt.bg = c(NA, pop_df$color[match(pop_order, pop_df$pop)][-1]),
       pt.cex = c(1.8, rep(1.25, 6)),
       col = 'black', cex = 0.8,
       pch = c(18, rep(21, 6)), bty = "n")
dev.off()



### POPULATON INFO ###
#"ANF", "CBJTC", "FTB", "FTS", "ONF", "WSF-CITRUS", "APAFR"

#ANF (Apalachicola National Forest?): 30.207865, -84.663056
#CBJTC (Camp Blanding Joint Training Center?): 29.984685, -81.989727
#FTB (UNCLEAR): 31.447247, -84.345396
#FTS (UNCLEAR): 33.504599, -83.389586
#ONF (Osceola National Forest?): 30.200350, -82.442920
#WSF-CITRUS (Withlacoochee State Forest?): 28.716566, -82.493976
#APAFR: 27.644521, -81.343642

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









