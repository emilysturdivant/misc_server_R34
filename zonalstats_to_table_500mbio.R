library(raster)
library(dplyr)
library(ggplot2)
library(reshape2)
library(flextable)
library(officer)
library(rasterVis)
library(gridExtra)
library(tidyverse)

params <- list(
  shp = '/mnt/data3/esturdivant/data/shp/hih_sites/all_sites.shp',
  years = c(2003, 2018),
  units = 'sqkm',
  output_fp = '/mnt/data3/esturdivant/products/for_hih/carbon_report_2018.csv'
)


rasterOptions(tmpdir = '~/tmp')

# input shapefile
shp.file <- params$shp
sp <- shapefile(shp.file)

# time period
start.year <- params$years[1] # 2003
end.year <- params$years[2] # 2017 or 2018
years <- seq(from = start.year, to = end.year, by = 1)

# units for varying (english or metric)
sq.units <- params$units # 'sqkm' or 'sqmi'

# 2003-2018 (V1)
agb.mgha <- stack('/mnt/data3/biomass/final/bio_change/2003_2018/v1/annfit/c95/f03/global.vrt', 
                  bands = 1:length(years))

# reproject shapefile to match modis
sp.p <- spTransform(sp, crs(agb.mgha))

# get area of polygon
area.sqm <- raster::area(sp.p)
tot.area.ha <- sum(area.sqm)/10000

# crop/mask raster stack
agb.mgha.c <- crop(agb.mgha, sp.p)
agb.mgha.cm <- mask(agb.mgha.c, sp.p)

# compute pixel area in hectares
x.res.m <- xres(agb.mgha)
y.res.m <- yres(agb.mgha)
px.ha <- x.res.m * y.res.m * 1e-4

# extract first and last year in stack
agb.mgha.ly <- subset(agb.mgha, nlayers(agb.mgha))

# compute total metric tons of carbon in last year
tot.agb.mg <- raster::extract(agb.mgha.ly * px.ha, sp.p, 
                      fun = sum, na.rm=TRUE, df=TRUE)
tot.agb.mg %>% saveRDS('zonal_sums.rds')
poly.px.ct <- raster::extract(agb.mgha.ly, sp.p, 
                      fun = function(x, ...) length(x), 
                      na.rm = TRUE, df=TRUE)
poly.px.ct %>% saveRDS('zonal_cts.rds')

# Get DF with carbon measurements for each polygon
tot.c.df <- tot.agb.mg %>% 
  
  # get number of pixels in polygon
  left_join(poly.px.ct, by = 'ID') %>% 
  
  # compute average carbon density, total carbon, and total CO2
  rename(agb_mg = 2, ct = 3) %>% 
  mutate(c_mg = agb_mg / 2, 
         c_mgcha = c_mg / px.ha / ct,
         co2_mg = c_mg * 3.664,
         C_Mt = c_mg / 1e6,
         CO2_Mt = co2_mg / 1e6) %>% 
  
  # get site name
  cbind(name = sp.p$name,
        HIH_site = sp.p$HIH_site, 
        type = sp.p$type)

# Save DF to CSV
tot.c.df %>% 
  dplyr::select(HIH_site, 
                Name = name, 
                Type = type, 
                C_dens_MgCha = c_mgcha, 
                C_stock_Mg = c_mg, 
                CO2_Mg = co2_mg) %>% 
  write_csv(params$output_fp)

# Display table
tot.c.df %>% 
  dplyr::select(name, c_mgcha, c_mg, co2_mg) %>% 
  mutate(c_mgcha = round(c_mgcha, 1),
         site = stringr::str_to_title(site)) %>% 
  flextable() %>% 
  colformat_num(digits = 1) %>% 
  font(fontname = 'Times New Roman', part = 'all') %>% 
  fontsize(size = 12, part = 'all') %>% 
  set_header_labels(site = paste0('Division')) %>% 
  mk_par(j = 'c_mgcha', part = 'header', 
          value = as_paragraph('Carbon density (average tons C ha', 
                               as_sup('-1'), ')')) %>% 
  mk_par(j = 'c_mg', part = 'header', 
          value = as_paragraph('Total Carbon (tons C)')) %>% 
  mk_par(j = 'co2_mg', part = 'header', 
          value = as_paragraph('Carbon dioxide equivalent (tons CO', 
                               as_sub(as.integer(2)), ')')) %>% 
  hline_top(border = fp_border(width = 1), part = 'all') %>% 
  hline_bottom(border = fp_border(width = 1), part = 'all') %>% 
  set_table_properties(layout = "autofit", width = 0.5) %>% 
  set_caption(caption = paste0('Table 1: Total aboveground carbon as of ', 
                               end.year, ' for each internal division.'), 
              autonum = run_autonum('tab', '', ''))