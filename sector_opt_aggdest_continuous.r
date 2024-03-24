library(tidyverse)
library(ompr)
library(ompr.roi)
library(ompr.highs)
library(leaflet)
library(ggforce)
library(sf)


rm(list=ls())
gc(T,T,T)

options(scipen=9999)

array_2d_multiplication_fcn <- function(static_array, row_variable, column_variable){ 
  vapply(seq_along(static_array), function(k) static_array[
    row_variable[k], column_variable[k]], numeric(1L))  }

highs_solver_parameters <- list(mip_rel_gap=0.25,log_to_console=TRUE,time_limit = 180)
avg_sector_metric_pct_tolerance <- 0.10  #set tighter than the actual tolerance. Need some wiggle room
#for zips that are assigned to multiple sectors


sector_count <- 48 #48 prod recruiters + 12 SNCOICs in RS TWC


#import data from MCRC structure workbook
all_zips <- readxl::read_excel("C:/Users/r_j_a/Documents/USMC/MCRC Structure/FY20 - FY35 MCRC Structure Data for Distribution (Feb 2024).xlsx", 
                                  sheet = "MCRC Zip DATA")

rs_zips <- all_zips %>%
  janitor::clean_names() %>%
  dplyr::select(rs,rss,state,county,fips,zip,lat,long,metric=qep_2025) %>%
  dplyr::filter(rs == 'TWIN CITIES') %>%
  dplyr::mutate(zip4=substr(zip,1,4),
                  zip3=substr(zip,1,3)
               ) %>%
  dplyr::filter(!is.na(lat),!is.na(long),!is.nan(lat),!is.nan(long)) %>%
  dplyr::filter(metric>0) %>%
  dplyr::mutate(zip = as.character(zip))

rss_df <- cbind.data.frame(
  rss = c('BLOOMINGTON','BUFFALO','BURNSVILLE','COON RAPIDS','DULUTH','EAU CLAIRE','FARGO','MANKATO','ROSEVILLE','SIOUX CITY','SIOUX FALLS','ST CLOUD'),
  rss_abb = c('BL','BF','BU','CO','DU','EC','FA','MK','RV','SY','SF','SC'),
  lat = c(44.8114267,45.1782933,44.7733899,45.1740642,46.8049533,44.8071073,46.8640532,44.1728541,45.0060767,42.4503635,43.4994024,45.5537568),
  lng = c(-93.3291527,-93.8667234,-93.2766312,-93.3444056,-92.166503,-91.468971,-96.8528707,-93.9537317,-93.1566107,-96.352813,-96.7728009,-94.1913827),
  sectors = c(4,3,4,5,3,5,5,3,4,3,5,4),
  sncoic = c(1,1,1,1,1,1,1,1,1,1,1,1),
  rop = c(5,4,5,6,4,6,6,4,5,4,6,5)
)


#summarize destinations at zip4 level
dest_summarized <- rs_zips %>% 
  dplyr::group_by(zip4) %>%
  dplyr::summarise(lat=sum(lat*metric)/sum(metric),
                   long=sum(long*metric)/sum(metric)) %>%
  dplyr::rename(zip=zip4) 
  



zipcount <- nrow(rs_zips)
destcount <- nrow(dest_summarized)
rsscount <- nrow(rss_df)

bigM <- sum(rs_zips$metric) * 1.1
avg_sector_metric <- sum(rs_zips$metric)/sector_count


#zip to sector dest distance matrix, multiplied by volume for objective function
zip_dest_dist <- geosphere::distm(
  x=cbind(rs_zips$long,rs_zips$lat),
  y=cbind(dest_summarized$long,dest_summarized$lat)) %>%
  measurements::conv_unit('m','mi') * 0.0001 #scale down for numeric performance in solver

#sector dest to rss distance matrix, multiplied by volume for objective function
dest_rss_dist <- geosphere::distm(
  x=cbind(dest_summarized$long,dest_summarized$lat),
  y=cbind(rss_df$lng,rss_df$lat)) %>%
  measurements::conv_unit('m','mi') * 0.0001 #scale down for numeric performance in solver



mymodel <- ompr::MILPModel() %>%
  add_variable(zip_dest_assign[orig,dest],orig=1:zipcount,dest=1:destcount,
               type='continuous',lb=0) %>%
  add_variable(zipdest_popn[dest], dest=1:destcount,type='continuous',lb=0) %>%
  add_variable(zipdest_open[dest],dest=1:destcount,type='binary') %>%
  add_variable(dest_rss_assign[dest,rss],dest=1:destcount,rss=1:rsscount,type='binary') %>%
  
#constraint: for each orig zip, sum of metric equals metric
  add_constraint(
    sum_expr(zip_dest_assign[orig,dest],dest=1:destcount) == rs_zips$metric[orig],
    orig = 1:zipcount) %>%

#constraint: for each dest zip, relate zipdest_popn to sum of assigned popn
  add_constraint(
    sum_expr(zip_dest_assign[orig,dest],orig=1:zipcount) == zipdest_popn[dest],
    dest = 1:destcount) %>%
    
  #constraint: activate zipdest_open if zipdest_popn greater than zero
  add_constraint(zipdest_popn[dest] <= bigM*zipdest_open[dest], dest=1:destcount) %>%
  
  #constraint: must open exactly sector_count destinations
  add_constraint(sum_expr(zipdest_open[dest],dest=1:destcount) == sector_count) %>%

#constraint: if opening sector, constrain popn to within X% of average

add_constraint(zipdest_popn[dest] >= (1-avg_sector_metric_pct_tolerance) * 
                 avg_sector_metric * zipdest_open[dest], dest=1:destcount) %>%
add_constraint(zipdest_popn[dest] <= (1+avg_sector_metric_pct_tolerance) * 
                 avg_sector_metric * zipdest_open[dest], dest=1:destcount) %>%
  
  #constraint: each dest assigned to 0 RSSs if dest is not opened; 1 if opened
  add_constraint(sum_expr(dest_rss_assign[dest,rss],rss=1:rsscount)== zipdest_open[dest], 
                 dest=1:destcount) %>%

  #constraint: each RSS has same sector count as recruiters assigned
  add_constraint(sum_expr(dest_rss_assign[dest,rss],dest=1:destcount) == rss_df$sectors[rss],
                 rss=1:rsscount) %>%
  
#objective function:

  set_objective(
    #minimize zip_dest_dist * assigned volume
    sum_expr( ompr::colwise(
    array_2d_multiplication_fcn(static_array=zip_dest_dist,
row_variable=orig,column_variable=dest)) *
  zip_dest_assign[orig,dest], orig=1:zipcount,dest=1:destcount) +
  
  #minimize dest_rss_dist * dest to rss assignments
  
  sum_expr( ompr::colwise(
    array_2d_multiplication_fcn(static_array=dest_rss_dist,
                                row_variable=dest,column_variable=rss)) *
      dest_rss_assign[dest,rss],dest=1:destcount,rss=1:rsscount)
  
  
  ,sense='min')

solve_out <-  mymodel %>% solve_model(highs_optimizer(control=highs_solver_parameters))


 #retrieve solutions
zip_dest_assign_soln <- get_solution(solve_out , zip_dest_assign[orig,dest]) %>%
#assign each origin to only a single destination. winner takes all.
dplyr::filter(value>0) %>%
dplyr::group_by(orig) %>%
 dplyr::mutate(pct_popn = value/sum(value),
               count_dest = length(unique(dest)),
               total_value = sum(value) ) # %>%

# test not using winner takes all logic
  #winner takes all logic for zips assigned to multiple sectors
  #dplyr::mutate(max_pct = pct_popn == max(pct_popn)) %>%
  #dplyr::filter(max_pct == TRUE) %>%
  #dplyr::mutate(value = total_value) %>%
  #dplyr::select(-total_value)
  

#check by origin to ensure all assigned
origin_soln <- zip_dest_assign_soln %>% dplyr::group_by(orig )%>%
  dplyr::summarise(solver_metric=sum(value)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(metric = rs_zips$metric[orig]) %>%
  dplyr::mutate(absdiff = abs(solver_metric-metric))

summary(origin_soln$absdiff)

dest_soln  <- zip_dest_assign_soln %>% dplyr::group_by(dest )%>%
  dplyr::summarise(solver_metric=sum(value)) %>%
  dplyr::filter(solver_metric>0.01) %>%
  dplyr::mutate(pct = solver_metric / avg_sector_metric)

#check if any dest_soln are outside the allotted bounds. May happen when
#reassigning ZIPs to a single destination with winner-takes-all
summary(dest_soln$pct)


dest_open_soln <- get_solution(solve_out, zipdest_open[dest]) %>%
  dplyr::mutate(value=round(value)) %>%  dplyr::filter(value==1)

dest_rss_soln <- get_solution(solve_out, dest_rss_assign[dest,rss]) %>%
  dplyr::mutate(value=round(value)) %>%  dplyr::filter(value==1)

#mapping

#colors <- RColorBrewer::brewer.pal(n=sector_count,"Set1")
colors <- c('orange','blue','green','red','yellow','purple')
colors <- rep_len(colors,sector_count)
colordf <- cbind.data.frame(dest = unique(zip_dest_assign_soln$dest),colors,stringsAsFactors=FALSE)

zip_df <- zip_dest_assign_soln %>% 
  #since we are allowing a single zip to be aligned to multiple sectors for now,
  #limit to just a single row in zip_df and signify multi-alignment by black marker
  dplyr::group_by(orig) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::left_join(colordf) %>%
  dplyr::mutate(colors = dplyr::if_else(count_dest>1,"black",colors)) %>%
  dplyr::left_join(rs_zips %>% 
                     dplyr::ungroup() %>%
                     dplyr::mutate(orig = dplyr::row_number())) %>%
  dplyr::select(orig, zip,lat,long,metric,dest,count_dest,colors)






#choropleth map of all zips

#shapefile from:
#https://catalog.data.gov/dataset/tiger-line-shapefile-2019-2010-nation-u-s-2010-census-5-digit-zip-code-tabulation-area-zcta5-na
# Load the shapes

postal.codes=sf::st_read(
  "C:/Users/r_j_a/Documents/USMC/MCRC Structure/tl_2019_us_zcta510/tl_2019_us_zcta510.shp") %>%
  dplyr::filter(ZCTA5CE10 %in% rs_zips$zip)

postal.codes2 <- postal.codes %>% 
  dplyr::left_join(zip_df %>% 
  #dplyr::mutate(ZCTA5CE10 = as.character(zip)) %>%
    dplyr::ungroup() %>%
    dplyr::select(ZCTA5CE10 = zip, dest, colors)) 

leaflet_all_sectors <- leaflet(postal.codes2) %>% 
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(color=~colors,weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5)


#identify any sectors that consist of non-contiguous components

#use sf::st_union to check for disconnectedness
dest_soln$union_class <- NA

for(destindex in 1:sector_count){
  print(destindex)
  sector_postal_codes <- postal.codes %>% 
    filter(ZCTA5CE10 %in% zip_df$zip[zip_df$dest == dest_soln$dest[destindex]])
  
  sector_union <- sf::st_union(sector_postal_codes)
  
  dest_soln$union_class[destindex] <- class(sector_union)[1]
}  

#make leaflet maps for each sector
sector_maps <- list()

for(i in 1:sector_count){
  print(i)
  
  dest_ <- dest_soln$dest[i]
  union_class_ <- dest_soln$union_class[i]
  
  data_ <- postal.codes2 %>% dplyr::filter(dest == dest_)
  
  leaflet_sector <- leaflet(data_) %>% 
    addProviderTiles(providers$CartoDB.Positron) %>%
    addPolygons(color=~colors,weight = 1, smoothFactor = 0.5,
                opacity = 1.0, fillOpacity = 0.5) %>%
    addControl(paste(
      "Dest ID: ", dest_,"<br>",
      "Union Class: ",union_class_,"<br>",
      "ZCTA Count: ", nrow(data_),
      sep= " "),position = "topleft") %>%
    #minimap
    leaflet::addMiniMap()
  
  sector_maps[[i]] <- leaflet_sector
  
}

leafsync::latticeview(sector_maps)
